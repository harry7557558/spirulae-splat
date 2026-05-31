#pragma once

// EngineState — all persistent device-resident state shared across the
// engine_* entrypoints. State is a singleton struct accessed via `engine()`.
// Sub-structs group related fields (world, camera, fwd, gt, grad, optim,
// bilagrid_rgb/depth/normal, ppisp). The singleton lives in EngineState.cpp.

#include "Tensor.h"

#include "BackgroundSphericalHarmonics.cuh"
#include "BilagridUtils.cuh"
#include "Densify.cuh"
#include "FusedProjectionBwdOptim.cuh"
#include "FusedSSIM.cuh"
#include "IntersectTile.cuh"
#include "Optimizer.cuh"
#include "PerPixelLoss.cuh"
#include "PerSplatLoss.cuh"
#include "PixelWise.cuh"
#include "Primitive.cuh"
#include "ProjectionBwd.cuh"
#include "ProjectionFwd.cuh"
#include "ProjectionPackedFwd.cuh"
#include "RasterizationBwd.cuh"
#include "RasterizationEval3DBwd.cuh"
#include "RasterizationEval3DFwd.cuh"
#include "RasterizationFwd.cuh"
#include "RasterizationSortedEval3DBwd.cuh"
#include "RasterizationSortedEval3DFwd.cuh"
#include "SVHash.cuh"
#include "SplatTileIntersector.cuh"
#include "Visualizer.cuh"

#include "gsplat/Common.h"

#include <cstdint>
#include <string>
#include <vector>


// World splat parameters (allocated once at init; persistent on device).
struct WorldSplats {
    bool initialized = false;
    DeviceVector<float3>   means;
    DeviceVector<float4>   quats;
    DeviceVector<float3>   scales;
    DeviceVector<float>    opacities;
    DeviceVector<float3>   features_dc;
    DeviceTensor2D<float3> features_sh;   // [max_N, num_sh]
};

// Camera table (re-copied each frame into pool).
struct CameraTable {
    int32_t num    = 0;
    int32_t width  = 0;
    int32_t height = 0;
    ssplat::CameraModelType model = (ssplat::CameraModelType)-1;
    std::string model_str;
    DeviceTensor2D<float4> viewmats;      // [C, 4]
    DeviceVector<float4>   intrins;       // [C]
    DeviceTensor2D<float>  dist_coeffs;   // [C, 10]
};

// Forward-pass intermediates retained for backward.
struct ForwardCache {
    DeviceVector<int32_t>             camera_ids;
    DeviceVector<int32_t>             gaussian_ids;
    DeviceTensor2D<float4>            aabb;         // [nnz,1] packed or [C,N] non-packed
    std::vector<DeviceTensorFloatND>  splats_w;
    std::vector<DeviceTensorFloatND>  splats_s;
    DeviceTensor3D<int32_t>           tile_offsets;
    DeviceVector<int32_t>             flatten_ids;
    DeviceTensor3D<float>             render_Ts;
    DeviceTensor3D<int32_t>           last_ids;
    RenderOutput::TensorTuple         renders;
    DeviceVector<float>               accum_weight; // [max_num_splats] per-splat score from raster bwd
    // Screen-space gradient handed from rasterize_*_bwd. The default path
    // consumes this immediately inside projection_*_backward; the fused
    // projection-bwd+optim path stashes it here for the optim step.
    std::vector<DeviceTensorFloatND>  v_splats_s;
};

// Training ground truth (re-copied each batch).
struct GTData {
    DeviceTensor3D<float3> rgb;
    DeviceTensor3D<float>  depth;
    DeviceTensor3D<float3> normal;
    DeviceTensor3D<bool>   alpha;
    bool has_gt   = false;
    bool has_mask = false;
};

// Splat parameter gradients (pool-backed, zeroed each backward).
struct SplatGrad {
    DeviceVector<float3>   means;
    DeviceVector<float4>   quats;
    DeviceVector<float3>   scales;
    DeviceVector<float>    opacities;
    DeviceVector<float3>   features_dc;
    DeviceTensor2D<float3> features_sh;
};

// Adam moment + densify aux state (pool-backed, persistent across steps).
struct SplatOptim {
    bool initialized = false;
    bool quantize_sh = false;
    bool use_per_splat_bias_correction = false;

    DeviceVector<float3>   g1_means,         g2_means;
    DeviceVector<float4>   g1_quats,         g2_quats;
    DeviceVector<float3>   g1_scales,        g2_scales;
    DeviceVector<float>    g1_opacities,     g2_opacities;
    DeviceVector<float3>   g1_features_dc,   g2_features_dc;
    DeviceTensor2D<float3> g1_features_sh,   g2_features_sh;       // when !quantize
    QuantizedAdamState<8, 256> sh_quant_state;                      // when  quantize
    // FPBO uses a different per-bound granularity (one float4 per 256 splats,
    // covering all 3*K cells per splat), so the kernel cannot reuse
    // `sh_quant_state`. Allocated only when use_fused_proj_bwd_optim + quantize_sh.
    QuantizedAdamState<8, 256> sh_quant_state_fpbo;

    DeviceVector<float>    radii;                  // [max_N]
    DeviceVector<float2>   accum_buffer;           // [max_N]
    DeviceVector<int32_t>  bias_correction_steps;  // [max_N], or empty

    // Set per-step from cfg.optim.use_fused_proj_bwd_optim before forward/loss
    // so engine_compute_loss_backward knows to skip projection_*_backward and
    // stash v_splats_s for the subsequent fused projection-bwd-and-optim call.
    bool   use_fused_proj_bwd_optim = false;
    // Mirrors use_fused_proj_bwd_optim at the time _ensure_optim_state last
    // ran. Drives reallocation when the fused flag flips between steps (the
    // FPBO SH-quant layout differs from the regular Optimizer.cu layout).
    bool   fused_state_active       = false;
};

// One bilagrid channel (RGB / depth / normal).
struct BilagridRGB {
    DeviceTensor5D<float>      grids;
    DeviceTensor5D<float>      image_grad;
    DeviceTensor5D<float>      g1, g2;
    QuantizedAdamState<8, 256> quant_state;
    bool        quantize_optim     = false;
    bool        enabled            = false;
    bool        optim_initialized  = false;
    DeviceTensor3D<float3>     fwd_pre;            // pre-bilagrid render
    std::string type;
    int         C = 0;
};
struct BilagridDepth {
    DeviceTensor5D<float>      grids;
    DeviceTensor5D<float>      image_grad;
    DeviceTensor5D<float>      g1, g2;
    QuantizedAdamState<8, 256> quant_state;
    bool        quantize_optim     = false;
    bool        enabled            = false;
    bool        optim_initialized  = false;
    DeviceTensor3D<float>      fwd_pre;
    DeviceVector<float>        scalars;            // per-camera median quantile
};
struct BilagridNormal {
    DeviceTensor5D<float>      grids;
    DeviceTensor5D<float>      image_grad;
    DeviceTensor5D<float>      g1, g2;
    QuantizedAdamState<8, 256> quant_state;
    bool        quantize_optim     = false;
    bool        enabled            = false;
    bool        optim_initialized  = false;
    DeviceTensor3D<float3>     fwd_pre;
};

// Background blending. Applied BEFORE bilagrid/PPISP. Two modes:
//   - Noise: random per-pixel color (warmup-weighted). No persistent state.
//   - Sh:    skybox = SH(world ray dir) + DC color. DC + L1+ coeffs trained.
struct EngineBackground {
    enum class Mode { None = 0, Noise = 1, Sh = 2 };
    Mode mode    = Mode::None;
    bool enabled = false;

    // Common config (set at init time)
    bool splat_color_is_linear = false;  // noise mode: sRGB->linear conversion

    // SH mode config
    int  sh_degree       = 0;            // 0..4

    // SH parameters: index 0 is DC (the "background_color"), 1..(deg+1)^2-1 are
    // higher-order SH bands. Stored once on device, persistent across steps.
    DeviceVector<float3> sh_coeffs;

    // SH Adam moments (allocated lazily on first optim_step when train_sh_color).
    DeviceVector<float3> sh_g1, sh_g2;
    bool sh_optim_initialized = false;

    // Per-iter buffers (resized each forward).
    // - fwd_pre_blend_rgb: saved pre-blend rendered RGB (Sh mode; needed by
    //   blend backward to compute v_background -> v_sh).
    // - fwd_background:    skybox image (Sh mode).
    DeviceTensor3D<float3> fwd_pre_blend_rgb;
    DeviceTensor3D<float3> fwd_background;

    // Stashed at forward time so backward can reconstruct the random skybox
    // (Noise mode) and so the optim step can run after backward (both modes).
    uint32_t cur_seed             = 0;
    float    cur_randomize_weight = 0.0f;
};

// PPISP (RGB-only photometric correction; applied AFTER bilagrid).
struct PpispState {
    DeviceTensor2D<float>  params;
    DeviceTensor2D<float>  grads;
    DeviceTensor2D<float>  g1, g2;
    DeviceTensor3D<float3> fwd_pre;
    std::string param_type;
    int  num_params         = 0;
    bool enabled            = false;
    bool optim_initialized  = false;
};


struct EngineState {
    // Top-level scalars + flags shared across sections.
    std::string primitive;
    int64_t cur_num_splats = 0;
    int64_t max_num_splats = 0;
    int     num_sh         = 0;
    int     sh_degree      = 0;
    bool    packed         = false;

    WorldSplats    world;
    CameraTable    camera;
    ForwardCache   fwd;
    GTData         gt;
    SplatGrad      grad;
    SplatOptim     optim;

    BilagridRGB    bilagrid_rgb;
    BilagridDepth  bilagrid_depth;
    BilagridNormal bilagrid_normal;
    // Per-image-in-batch camera indices for the current step (shared by
    // background, bilagrid, and PPISP). Empty -> kernels fall back to identity.
    DeviceVector<int32_t> bilagrid_cur_cam_indices;

    EngineBackground background;
    PpispState       ppisp;
};


// Function-local-static singleton. Future evolution: a second instance for
// eval / viewer / multi-model can be added by swapping the accessor's body.
EngineState& engine();
