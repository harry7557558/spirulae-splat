#include "Engine.h"
#include "Tensor.h"
#include "Primitive.cuh"
#include "BilagridBindings.h"
#include "ProjectionFwd.cuh"
#include "ProjectionPackedFwd.cuh"
#include "ProjectionBwd.cuh"
#include "IntersectTile.cuh"
#include "RasterizationFwd.cuh"
#include "RasterizationBwd.cuh"
#include "RasterizationEval3DFwd.cuh"
#include "RasterizationEval3DBwd.cuh"

#include "npy.hpp"
#include <filesystem>
#include <fstream>
#include <cmath>

// Forward decls for the raw-pointer image conversion helpers in PixelWise.cu.
// Not in PixelWise.cuh because they bypass the DeviceTensor3D shape carrier.
void uint8_image_to_float_raw(const uint8_t* d_in, float* d_out, int B, int H, int W, int C);
void uint16_image_to_float_raw(const uint16_t* d_in, float* d_out, int B, int H, int W, int C);
void uint8_normal_to_float_raw(const uint8_t* d_in, float* d_out, int B, int H, int W, int C);
void uint16_depth_to_float_raw(const uint16_t* d_in, float* d_out, int B, int H, int W, int C);

// Fused image-grad + TV-loss + Adam optimizer step for bilagrid (Phase B).
// Lives in BilagridFusedAdam.cu.
void fused_bilagrid_tv_adam(
    float* grids,
    float*   g1_f,    float*   g2_f,
    uint8_t* packed,                    // AoS packed (u, sqrt_g2) when quantize
    float4*  quant_bounds,
    const float* image_grad,
    const int*   cam_indices,
    int N_grids, int C_batch, int C,
    int L, int H, int W,
    float lr,
    float tv_weight,
    int32_t adam_step,
    bool quantize,
    cudaStream_t stream
);


namespace Buffers {

std::string primitive;
int64_t cur_num_splats;
int64_t max_num_splats;
int sh_degree;

// World splats — allocated once at init, persistent on device
bool splats_initialized = false;
DeviceVector<float3> world_means;
DeviceVector<float4> world_quats;
DeviceVector<float3> world_scales;
DeviceVector<float> world_opacities;
DeviceVector<float3> world_features_dc;
DeviceTensor2D<float3> world_features_sh;  // [max_N, num_sh]

int32_t num_cameras;
int32_t image_width;
int32_t image_height;
ssplat::CameraModelType camera_model;
std::string camera_model_str;
// Camera params — re-copied each frame to pool
DeviceTensor2D<float4> camera_viewmats;     // [C, 4]  (4 rows of float4)
DeviceVector<float4> camera_intrins;        // [C]
DeviceTensor2D<float> camera_dist_coeffs;   // [C, 10]

bool packed;

// Forward intermediate results (kept for backward pass)
DeviceVector<int32_t> fwd_camera_ids;
DeviceVector<int32_t> fwd_gaussian_ids;
DeviceTensor2D<float4> fwd_aabb;   // [nnz,1] packed or [C,N] non-packed
std::vector<DeviceTensorFloatND> fwd_splats_w;
std::vector<DeviceTensorFloatND> fwd_splats_s;
DeviceTensor3D<int32_t> fwd_tile_offsets;
DeviceVector<int32_t> fwd_flatten_ids;
DeviceTensor3D<float> fwd_render_Ts;
DeviceTensor3D<int32_t> fwd_last_ids;
RenderOutput::TensorTuple fwd_renders;  // for 3dgut backward
DeviceVector<float> fwd_accum_weight;  // [max_num_splats] per-splat score from raster bwd

// Training ground truth (re-copied each batch to pool).
// The 4 mask buffers (rgb_mask, depth_mask, normal_mask, alpha_mask) were
// dropped: they are now derived inside per_pixel_losses.slang from gt_alpha
// (RGB mask), gt_depth (depth + alpha masks), and gt_normal (normal mask).
DeviceTensor3D<float3> gt_rgb;         // [C, H, W]
DeviceTensor3D<float> gt_depth;        // [C, H, W]
DeviceTensor3D<float3> gt_normal;      // [C, H, W]
DeviceTensor3D<bool> gt_alpha;         // [C, H, W] external mask; sets `has_mask` for the slang kernel
bool has_gt = false;
bool has_mask = false;                 // gt_alpha buffer present (controls per-pixel mask in slang kernel)

// Pool-backed training state
int num_sh = 0;
bool optim_initialized = false;
bool train_quantize_sh = false;

// Gradients (pool-backed, zeroed each backward)
DeviceVector<float3> grad_means;
DeviceVector<float4> grad_quats;
DeviceVector<float3> grad_scales;
DeviceVector<float> grad_opacities;
DeviceVector<float3> grad_features_dc;
DeviceTensor2D<float3> grad_features_sh;

// Optimizer state g1/g2 (pool-backed, persistent across steps)
DeviceVector<float3> g1_means, g2_means;
DeviceVector<float4> g1_quats, g2_quats;
DeviceVector<float3> g1_scales, g2_scales;
DeviceVector<float> g1_opacities, g2_opacities;
DeviceVector<float3> g1_features_dc, g2_features_dc;
DeviceTensor2D<float3> g1_features_sh, g2_features_sh;       // when !quantize
// When quantize: joint (u, sqrt(g2)) packed Adam state for SH features.
// One thread-block worth of SH cells (BLOCK_SIZE splats * K coeffs * 3 ch)
// shares a single float4 bounds slot; codec primitives live on the class.
QuantizedAdamState<8, 256> sh_quant_state;

// Training aux (pool-backed)
DeviceVector<float> radii;                     // [max_N]
DeviceVector<float2> accum_buffer;             // [max_N]
DeviceVector<int32_t> bias_correction_steps;   // [max_N], or empty
bool use_per_splat_bias_correction = false;

// ============================================================================
// Bilagrid state (RGB / depth / normal).
//
// Per type we hold: grid parameters, accumulated gradients (zeroed each step),
// Adam moments g1 / g2 (persistent), and a saved pre-bilagrid forward output
// used to drive the backward pass.
//
// RGB applies to the rendered prediction (fwd_renders.rgb). Depth and normal
// apply to the GT side (gt_depth / gt_normal), matching the Python flow.
// ============================================================================

// Bilagrid Adam state can be stored quantized (uint8 + per-block float4 bounds)
// like the SH path. quantize_optim is a per-type flag captured at init time.

// Phase B: bilagrid_*_image_grad is the SPARSE per-batch image-loss gradient,
// sized [C_batch, C, L, H, W] (re-sized each iter; pool high-water-marks).
// Replaces the dense [N_grids, C, L, H, W] grads table — saves ~195 MB on
// dome-scale (16k+ cameras) datasets. The TV-loss gradient is computed inline
// in the fused optimizer kernel rather than materialized to a buffer.

// RGB. Channel-last layout: 9 / 12 channels are innermost so the trilinear
// sampler reads C contiguous floats per corner instead of stride-L*H*W ones.
DeviceTensor5D<float> bilagrid_rgb_grids;   // [N_cam, L, H, W, C_rgb]
DeviceTensor5D<float> bilagrid_rgb_image_grad;  // [C_batch, L, H, W, C_rgb] (sparse over cams)
DeviceTensor5D<float> bilagrid_rgb_g1, bilagrid_rgb_g2;        // float (when !quantize)
QuantizedAdamState<8, 256> bilagrid_rgb_quant_state;           // when quantize
bool bilagrid_rgb_quantize_optim = false;
DeviceTensor3D<float3> fwd_rgb_pre_bilagrid; // [C, H, W]
std::string bilagrid_rgb_type;               // "affine" | "ppisp" | "loglinear"
int bilagrid_rgb_C = 0;                      // 12 for affine, 9 for ppisp/loglinear
bool bilagrid_rgb_enabled = false;
bool bilagrid_rgb_optim_initialized = false;

// Depth (gt-side)
DeviceTensor5D<float> bilagrid_depth_grids;  // [N_cam, L, H, W, 2]
DeviceTensor5D<float> bilagrid_depth_image_grad;
DeviceTensor5D<float> bilagrid_depth_g1, bilagrid_depth_g2;
QuantizedAdamState<8, 256> bilagrid_depth_quant_state;
bool bilagrid_depth_quantize_optim = false;
DeviceTensor3D<float> fwd_depth_pre_bilagrid;
DeviceVector<float> bilagrid_depth_scalars;  // [C], median-quantile of gt_depth
bool bilagrid_depth_enabled = false;
bool bilagrid_depth_optim_initialized = false;

// Normal (gt-side)
DeviceTensor5D<float> bilagrid_normal_grids; // [N_cam, L, H, W, 3]
DeviceTensor5D<float> bilagrid_normal_image_grad;
DeviceTensor5D<float> bilagrid_normal_g1, bilagrid_normal_g2;
QuantizedAdamState<8, 256> bilagrid_normal_quant_state;
bool bilagrid_normal_quantize_optim = false;
DeviceTensor3D<float3> fwd_normal_pre_bilagrid;
bool bilagrid_normal_enabled = false;
bool bilagrid_normal_optim_initialized = false;

// Per-image-in-batch camera indices for the current step (used by both
// bilagrid and PPISP). Length equals current batch size (Buffers::num_cameras).
// Empty when no cam_idx was supplied — kernels then fall back to identity
// (image `bid` -> grid slot `bid`).
DeviceVector<int32_t> bilagrid_cur_cam_indices;

// ============================================================================
// PPISP state (per-camera RGB-only photometric correction).
//
// Mirrors the bilagrid layout: parameter table, gradient accumulator, Adam
// moments, and a saved pre-PPISP RGB used to drive backward. PPISP is applied
// AFTER bilagrid in forward, so its backward runs BEFORE bilagrid's backward.
// ============================================================================
DeviceTensor2D<float> ppisp_params;          // [N_cam, P]
DeviceTensor2D<float> ppisp_grads;           // [N_cam, P]
DeviceTensor2D<float> ppisp_g1, ppisp_g2;    // Adam moments
DeviceTensor3D<float3> fwd_rgb_pre_ppisp;    // [C, H, W]
std::string ppisp_param_type;                // "original" | "rqs"
int ppisp_num_params = 0;                    // 36 (original) | 39 (rqs)
bool ppisp_enabled = false;
bool ppisp_optim_initialized = false;
} // namespace Buffers

// Forward declaration for affine identity initializer (BilagridInit.cu).
void bilagrid_affine_identity_init(
    float* grids, int N, int L, int H, int W, cudaStream_t stream);

// Forward declaration for PPISP "original" default initializer (PpispInit.cu).
void ppisp_original_default_init(
    float* params, int N, cudaStream_t stream);

// Forward declaration for elementwise grad accumulation (PpispInit.cu).
void ppisp_add_into_grad(
    const float* src, float* dst, size_t n, cudaStream_t stream);

// Forward declaration for scattering per-image scalars into a per-camera table.
void bilagrid_scatter_floats(
    const float* src, const int* indices, int n,
    float* dst, cudaStream_t stream);


// Create DeviceTensorFloatND from TorchTensorView by appending trailing 1 (for float-typed access)
static DeviceTensorFloatND tv_to_fnd(const TorchTensorView& tv) {
    auto shape = std::get<2>(tv);
    shape.push_back(1LL);
    return DeviceTensorFloatND(TorchTensorView{std::get<0>(tv), std::get<1>(tv), shape});
}

// Create DeviceTensor2D<float4> view from DeviceVector<float4> (for packed projection backward)
static DeviceTensor2D<float4> vec_to_2d_float4(const DeviceVector<float4>& vec) {
    TorchTensorView tv{(uint64_t)vec.data_ptr(), (uint32_t)sizeof(float), {vec.size(), 1LL, 4LL}};
    return DeviceTensor2D<float4>(tv);
}

// --- Conversion helpers: DeviceVector/DeviceTensor -> TorchTensorView (for external functions) ---

template<typename T>
static TorchTensorView _dv_tv(const DeviceVector<T>& dv) {
    static_assert(sizeof(T) % 4 == 0 || std::is_same<T, uint8_t>::value || std::is_same<T, bool>::value,
                  "Unsupported element type");
    if constexpr (sizeof(T) % 4 == 0) {
        return TorchTensorView((uint64_t)dv.data_ptr(), 4, {dv.size(), (int64_t)(sizeof(T) / 4)});
    } else {
        return TorchTensorView((uint64_t)dv.data_ptr(), 1, {dv.size(), (int64_t)sizeof(T)});
    }
}

template<typename T>
static TorchTensorView _dt2d_tv(const DeviceTensor2D<T>& dt) {
    // For scalar types (float, int, etc.) use 2D shape [n0, n1];
    // for vector types (float3, float4, ...) use 3D shape [n0, n1, channels].
    if constexpr (std::is_arithmetic_v<T>) {
        uint32_t es = (uint32_t)sizeof(T);
        return TorchTensorView((uint64_t)dt.data_ptr(), es,
            {dt.template size<0>(), dt.template size<1>()});
    } else if constexpr (sizeof(T) % 4 == 0) {
        return TorchTensorView((uint64_t)dt.data_ptr(), 4,
            {dt.template size<0>(), dt.template size<1>(), (int64_t)(sizeof(T) / 4)});
    } else {
        return TorchTensorView((uint64_t)dt.data_ptr(), 1,
            {dt.template size<0>(), dt.template size<1>(), (int64_t)sizeof(T)});
    }
}

template<typename T>
static TorchTensorView _dt3d_tv(const DeviceTensor3D<T>& dt) {
    if constexpr (sizeof(T) % 4 == 0) {
        return TorchTensorView((uint64_t)dt.data_ptr(), 4,
            {dt.template size<0>(), dt.template size<1>(), dt.template size<2>(),
             (int64_t)(sizeof(T) / 4)});
    } else {
        return TorchTensorView((uint64_t)dt.data_ptr(), 1,
            {dt.template size<0>(), dt.template size<1>(), dt.template size<2>(),
             (int64_t)sizeof(T)});
    }
}

// --- Source pointer location detection ---
// Returns true iff `ptr` refers to device-accessible memory (not pageable host).
static bool _is_device_ptr(const void* ptr) {
    if (ptr == nullptr) return false;
    cudaPointerAttributes attr{};
    cudaError_t err = cudaPointerGetAttributes(&attr, ptr);
    if (err != cudaSuccess) {
        // Pageable host pointer is unregistered -> reset error and treat as host.
        cudaGetLastError();
        return false;
    }
    return attr.type == cudaMemoryTypeDevice || attr.type == cudaMemoryTypeManaged;
}

// --- Pool allocation + H->D copy from a host TorchTensorView, OR zero-copy view of device source ---
// If source pointer is already device memory, returns a view (no pool allocation, no copy).
// If source is host, copies into pool with given key and returns view of pool buffer.

template<typename T>
static DeviceVector<T> _hv_to_dv(const std::string& key, const TorchTensorView& src_tv) {
    if (std::get<0>(src_tv) == 0) return DeviceVector<T>();
    int64_t n = std::get<2>(src_tv)[0];
    uint64_t src_ptr = std::get<0>(src_tv);
    T* ptr;
    if (_is_device_ptr((void*)src_ptr)) {
        ptr = (T*)src_ptr;  // zero-copy view
    } else {
        ptr = DevicePool::global().acquire<T>(key, (size_t)n);
        cudaMemcpy(ptr, (void*)src_ptr, n * sizeof(T), cudaMemcpyHostToDevice);
    }
    TorchTensorView dv_tv((uint64_t)ptr,
        (uint32_t)(sizeof(T) % 4 == 0 ? 4 : 1),
        {n, (int64_t)(sizeof(T) % 4 == 0 ? sizeof(T) / 4 : sizeof(T))});
    return DeviceVector<T>(dv_tv);
}

template<typename T>
static DeviceTensor2D<T> _hv_to_dt2d(const std::string& key, const TorchTensorView& src_tv) {
    if (std::get<0>(src_tv) == 0) return DeviceTensor2D<T>();
    auto& shape = std::get<2>(src_tv);
    int64_t n0 = shape[0], n1 = shape[1];
    uint64_t src_ptr = std::get<0>(src_tv);
    T* ptr;
    if (_is_device_ptr((void*)src_ptr)) {
        ptr = (T*)src_ptr;
    } else {
        ptr = DevicePool::global().acquire<T>(key, (size_t)(n0 * n1));
        cudaMemcpy(ptr, (void*)src_ptr, n0 * n1 * sizeof(T), cudaMemcpyHostToDevice);
    }
    TorchTensorView tv((uint64_t)ptr,
        (uint32_t)(sizeof(T) % 4 == 0 ? 4 : 1),
        {n0, n1, (int64_t)(sizeof(T) % 4 == 0 ? sizeof(T) / 4 : sizeof(T))});
    return DeviceTensor2D<T>(tv);
}

template<typename T>
static DeviceTensor3D<T> _hv_to_dt3d(const std::string& key, const TorchTensorView& src_tv) {
    if (std::get<0>(src_tv) == 0) return DeviceTensor3D<T>();
    auto& shape = std::get<2>(src_tv);
    int64_t n0 = shape[0], n1 = shape[1], n2 = shape[2];
    uint64_t src_ptr = std::get<0>(src_tv);
    T* ptr;
    if (_is_device_ptr((void*)src_ptr)) {
        ptr = (T*)src_ptr;
    } else {
        ptr = DevicePool::global().acquire<T>(key, (size_t)(n0 * n1 * n2));
        cudaMemcpy(ptr, (void*)src_ptr, n0 * n1 * n2 * sizeof(T), cudaMemcpyHostToDevice);
    }
    TorchTensorView tv((uint64_t)ptr,
        (uint32_t)(sizeof(T) % 4 == 0 ? 4 : 1),
        {n0, n1, n2, (int64_t)(sizeof(T) % 4 == 0 ? sizeof(T) / 4 : sizeof(T))});
    return DeviceTensor3D<T>(tv);
}

// Copy DeviceVector -> host TorchTensorView (D->H)
template<typename T>
static void _dv_to_host(const DeviceVector<T>& dv, const TorchTensorView& host_tv) {
    if (dv.data_ptr() == nullptr || std::get<0>(host_tv) == 0) return;
    int64_t n = dv.size();
    cudaMemcpy((void*)std::get<0>(host_tv), dv.data_ptr(), n * sizeof(T), cudaMemcpyDeviceToHost);
}

template<typename T>
static void _dt3d_to_host(const DeviceTensor3D<T>& dt, const TorchTensorView& host_tv) {
    if (dt.data_ptr() == nullptr || std::get<0>(host_tv) == 0) return;
    cudaMemcpy((void*)std::get<0>(host_tv), dt.data_ptr(), dt.numel() * sizeof(T), cudaMemcpyDeviceToHost);
}


void set_data_3dgs(
    int64_t num_splats,
    TorchTensorView means,
    TorchTensorView quats,
    TorchTensorView scales,
    TorchTensorView opacities,
    TorchTensorView features_dc,
    TorchTensorView features_sh
) {
    int64_t max_num_splats = std::get<2>(means)[0];
    if (std::get<2>(quats)[0] != max_num_splats ||
        std::get<2>(scales)[0] != max_num_splats ||
        std::get<2>(opacities)[0] != max_num_splats ||
        std::get<2>(features_dc)[0] != max_num_splats ||
        std::get<2>(features_sh)[0] != max_num_splats)
        throw std::runtime_error("setData3DGS: max_num_splats mismatch across splat tensors");
    if (max_num_splats < num_splats)
        throw std::runtime_error("setData3DGS: tensor size smaller than num_splats");

    Buffers::cur_num_splats = num_splats;
    Buffers::max_num_splats = max_num_splats;

    auto sh_shape = std::get<2>(features_sh);
    Buffers::num_sh = (sh_shape.size() >= 2) ? (int)sh_shape[1] : 0;

    // Splats are set once and persist on device; subsequent calls are no-ops.
    if (Buffers::splats_initialized) return;

    Buffers::world_means       = _hv_to_dv<float3>("world.means", means);
    Buffers::world_quats       = _hv_to_dv<float4>("world.quats", quats);
    Buffers::world_scales      = _hv_to_dv<float3>("world.scales", scales);
    Buffers::world_opacities   = _hv_to_dv<float>("world.opacities", opacities);
    Buffers::world_features_dc = _hv_to_dv<float3>("world.features_dc", features_dc);
    Buffers::world_features_sh = _hv_to_dt2d<float3>("world.features_sh", features_sh);

    Buffers::splats_initialized = true;
}


void set_camera_params(
    int width,
    int height,
    std::string camera_model,
    TorchTensorView viewmats,
    TorchTensorView intrins,
    TorchTensorView dist_coeffs
) {
    Buffers::image_width = width;
    Buffers::image_height = height;
    Buffers::camera_model = cmt(camera_model);
    Buffers::camera_model_str = camera_model;
    Buffers::num_cameras = std::get<2>(viewmats)[0];

    if (std::get<2>(intrins)[0] != Buffers::num_cameras ||
        std::get<2>(dist_coeffs)[0] != Buffers::num_cameras)
        throw std::runtime_error("setCameraParams: num_cameras mismatch");

    // viewmats: PyTorch shape [C, 4, 4] -> treat as DeviceTensor2D<float4> [C, 4]
    int64_t C = std::get<2>(viewmats)[0];
    {
        uint64_t src_ptr = std::get<0>(viewmats);
        float4* ptr;
        if (_is_device_ptr((void*)src_ptr)) {
            ptr = (float4*)src_ptr;
        } else {
            ptr = DevicePool::global().acquire<float4>("cam.viewmats", (size_t)(C * 4));
            cudaMemcpy(ptr, (void*)src_ptr, C * 4 * sizeof(float4), cudaMemcpyHostToDevice);
        }
        TorchTensorView dv_tv((uint64_t)ptr, 4, {C, 4LL, 4LL});
        Buffers::camera_viewmats = DeviceTensor2D<float4>(dv_tv);
    }
    Buffers::camera_intrins     = _hv_to_dv<float4>("cam.intrins", intrins);
    Buffers::camera_dist_coeffs = _hv_to_dt2d<float>("cam.dist_coeffs", dist_coeffs);
}


void forward_3dgs(
    std::string primitive,   // "3dgs", "mip", "3dgut"
    int sh_degree,
    bool packed
) {
    Buffers::primitive = primitive;
    Buffers::sh_degree = sh_degree;
    Buffers::packed = packed;

    // Build splats as DeviceTensorFloatND from typed device buffers.
    // For DeviceVector<float> (opacities), build a [N, 1] shape FND manually.
    DeviceTensorFloatND fnd_opac;
    {
        TorchTensorView opac_tv((uint64_t)Buffers::world_opacities.data_ptr(), 4,
            {Buffers::world_opacities.size(), 1LL, 1LL});
        fnd_opac = DeviceTensorFloatND(opac_tv);  // ndim=2, shape=[N, 1]
    }
    std::vector<DeviceTensorFloatND> in_splats = {
        DeviceTensorFloatND(Buffers::world_means),         // [N, 3]
        DeviceTensorFloatND(Buffers::world_quats),         // [N, 4]
        DeviceTensorFloatND(Buffers::world_scales),        // [N, 3]
        fnd_opac,                                          // [N, 1]
        DeviceTensorFloatND(Buffers::world_features_dc),   // [N, 3]
        DeviceTensorFloatND(Buffers::world_features_sh),   // [N, K, 3]
    };
    Buffers::fwd_splats_w = in_splats;

    DeviceTensorFloatND aabb_nd, depths_nd;

    // Allocate and zero radii buffer before projection (kernel uses atomicMax)
    Buffers::radii.resize("eng.radii", Buffers::max_num_splats);
    Buffers::radii.zero();

    // --- Projection ---
    if (packed) {
        DeviceVector<int32_t> cam_ids, gauss_ids;
        DeviceVector<float4> aabb_vec;
        DeviceVector<float> depths_vec;

        if (primitive == "3dgs") {
            auto [a, b, c, d, e] = projection_3dgs_packed_forward(
                Buffers::cur_num_splats, sh_degree, in_splats,
                _dt2d_tv(Buffers::camera_viewmats), _dv_tv(Buffers::camera_intrins),
                (uint32_t)Buffers::image_width, (uint32_t)Buffers::image_height,
                Buffers::camera_model_str, _dt2d_tv(Buffers::camera_dist_coeffs),
                Buffers::radii);
            cam_ids = a; gauss_ids = b; aabb_vec = c; depths_vec = d;
            Buffers::fwd_splats_s = e;
        } else if (primitive == "mip") {
            auto [a, b, c, d, e] = projection_mip_packed_forward(
                Buffers::cur_num_splats, sh_degree, in_splats,
                _dt2d_tv(Buffers::camera_viewmats), _dv_tv(Buffers::camera_intrins),
                (uint32_t)Buffers::image_width, (uint32_t)Buffers::image_height,
                Buffers::camera_model_str, _dt2d_tv(Buffers::camera_dist_coeffs),
                Buffers::radii);
            cam_ids = a; gauss_ids = b; aabb_vec = c; depths_vec = d;
            Buffers::fwd_splats_s = e;
        } else if (primitive == "3dgut") {
            auto [a, b, c, d, e] = projection_3dgut_packed_forward(
                Buffers::cur_num_splats, sh_degree, in_splats,
                _dt2d_tv(Buffers::camera_viewmats), _dv_tv(Buffers::camera_intrins),
                (uint32_t)Buffers::image_width, (uint32_t)Buffers::image_height,
                Buffers::camera_model_str, _dt2d_tv(Buffers::camera_dist_coeffs),
                Buffers::radii);
            cam_ids = a; gauss_ids = b; aabb_vec = c; depths_vec = d;
            Buffers::fwd_splats_s = e;
        } else {
            throw std::runtime_error("engine_forward: unknown primitive: " + primitive);
        }

        Buffers::fwd_camera_ids = cam_ids;
        Buffers::fwd_gaussian_ids = gauss_ids;
        Buffers::fwd_aabb = vec_to_2d_float4(aabb_vec);  // [nnz, 1] view for backward
        aabb_nd = DeviceTensorFloatND(aabb_vec);          // [nnz, 4] for intersect
        depths_nd = DeviceTensorFloatND(depths_vec);      // [nnz]   for intersect
    } else {
        DeviceTensor2D<float4> aabb_2d;
        DeviceTensor2D<float> depths_2d;  // [C, N] — sorting depths per camera

        if (primitive == "3dgs") {
            auto [a, b, c] = projection_3dgs_forward(
                Buffers::cur_num_splats, sh_degree, in_splats,
                _dt2d_tv(Buffers::camera_viewmats), _dv_tv(Buffers::camera_intrins),
                (uint32_t)Buffers::image_width, (uint32_t)Buffers::image_height,
                Buffers::camera_model_str, _dt2d_tv(Buffers::camera_dist_coeffs),
                Buffers::radii);
            aabb_2d = a; depths_2d = b;
            Buffers::fwd_splats_s = c;
        } else if (primitive == "mip") {
            auto [a, b, c] = projection_mip_forward(
                Buffers::cur_num_splats, sh_degree, in_splats,
                _dt2d_tv(Buffers::camera_viewmats), _dv_tv(Buffers::camera_intrins),
                (uint32_t)Buffers::image_width, (uint32_t)Buffers::image_height,
                Buffers::camera_model_str, _dt2d_tv(Buffers::camera_dist_coeffs),
                Buffers::radii);
            aabb_2d = a; depths_2d = b;
            Buffers::fwd_splats_s = c;
        } else if (primitive == "3dgut") {
            auto [a, b, c] = projection_3dgut_forward(
                Buffers::cur_num_splats, sh_degree, in_splats,
                _dt2d_tv(Buffers::camera_viewmats), _dv_tv(Buffers::camera_intrins),
                (uint32_t)Buffers::image_width, (uint32_t)Buffers::image_height,
                Buffers::camera_model_str, _dt2d_tv(Buffers::camera_dist_coeffs),
                Buffers::radii);
            aabb_2d = a; depths_2d = b;
            Buffers::fwd_splats_s = c;
        } else {
            throw std::runtime_error("engine_forward: unknown primitive: " + primitive);
        }

        Buffers::fwd_camera_ids = DeviceVector<int32_t>();
        Buffers::fwd_gaussian_ids = DeviceVector<int32_t>();
        Buffers::fwd_aabb = aabb_2d;                           // [C, N] for backward
        aabb_nd = DeviceTensorFloatND(aabb_2d);                // [C, N, 4] for intersect
        depths_nd = DeviceTensorFloatND(depths_2d);            // [C, N] -> numel=C*N for intersect
    }

    // --- Tile intersection (AABB mode) ---
    DeviceVector<int32_t>* image_ids_ptr = nullptr;
    if (packed && Buffers::fwd_camera_ids.data_ptr() != nullptr)
        image_ids_ptr = &Buffers::fwd_camera_ids;

    auto [isect_ids, flatten_ids, tile_offsets] = do_intersect_tile_generic(
        aabb_nd, depths_nd,
        nullptr, nullptr, nullptr,
        (uint32_t)Buffers::num_cameras,
        _dv_tv(Buffers::camera_intrins),
        (uint32_t)Buffers::image_width,
        (uint32_t)Buffers::image_height,
        image_ids_ptr
    );

    Buffers::fwd_tile_offsets = tile_offsets;
    Buffers::fwd_flatten_ids = flatten_ids;

    // --- Rasterization forward ---
    RenderOutput::TensorTuple renders;
    DeviceTensor3D<float> render_Ts;
    DeviceTensor3D<int32_t> last_ids;

    if (primitive == "3dgs") {
        auto [r, rTs, lids, r2, dist] = rasterize_to_pixels_3dgs_fwd(
            in_splats, Buffers::fwd_splats_s, Buffers::fwd_gaussian_ids,
            (uint32_t)Buffers::image_width, (uint32_t)Buffers::image_height,
            tile_offsets, flatten_ids, false);
        renders = r; render_Ts = rTs; last_ids = lids;
    } else if (primitive == "mip") {
        auto [r, rTs, lids, r2, dist] = rasterize_to_pixels_mip_fwd(
            in_splats, Buffers::fwd_splats_s, Buffers::fwd_gaussian_ids,
            (uint32_t)Buffers::image_width, (uint32_t)Buffers::image_height,
            tile_offsets, flatten_ids, false);
        renders = r; render_Ts = rTs; last_ids = lids;
    } else if (primitive == "3dgut") {
        auto [r, rTs, lids, r2, dist] = rasterize_to_pixels_3dgut_fwd(
            in_splats, Buffers::fwd_splats_s, Buffers::fwd_gaussian_ids,
            _dt2d_tv(Buffers::camera_viewmats), _dv_tv(Buffers::camera_intrins),
            Buffers::camera_model_str, _dt2d_tv(Buffers::camera_dist_coeffs),
            (uint32_t)Buffers::image_width, (uint32_t)Buffers::image_height,
            tile_offsets, flatten_ids, false);
        renders = r; render_Ts = rTs; last_ids = lids;
    }

    Buffers::fwd_renders = renders;
    Buffers::fwd_render_Ts = render_Ts;
    Buffers::fwd_last_ids = last_ids;

    // Results stay in pool — use engine_copy_render_to_host to fetch
}


void engine_debug_forward(
    TorchTensorView override_features_dc,
    int override_sh_degree,
    TorchTensorView out_rgb
) {
    if (std::get<0>(Buffers::fwd_renders).data_ptr() == nullptr)
        throw std::runtime_error("engine_debug_forward: forward_3dgs must be called first");

    // Swap in overrides (H->D copy for host tensor)
    DeviceVector<float3> saved_dc = Buffers::world_features_dc;
    int saved_sh = Buffers::sh_degree;
    if (std::get<0>(override_features_dc) != 0)
        Buffers::world_features_dc = _hv_to_dv<float3>("debug.features_dc", override_features_dc);
    if (override_sh_degree >= 0)
        Buffers::sh_degree = override_sh_degree;

    forward_3dgs(Buffers::primitive, Buffers::sh_degree, Buffers::packed);

    // Copy rgb result D->H
    auto& rgb = std::get<0>(Buffers::fwd_renders);
    if (rgb.data_ptr() && std::get<0>(out_rgb) != 0) {
        cudaMemcpy((void*)std::get<0>(out_rgb), rgb.data_ptr(),
                   rgb.numel() * sizeof(float3), cudaMemcpyDeviceToHost);
    }

    // Restore
    Buffers::world_features_dc = saved_dc;
    Buffers::sh_degree = saved_sh;
}


void engine_copy_accum_buffer(TorchTensorView dst) {
    if (Buffers::accum_buffer.data_ptr() == nullptr || std::get<0>(dst) == 0)
        return;
    int64_t dst_n = std::get<2>(dst)[0];
    int64_t src_n = Buffers::accum_buffer.size();
    int64_t n = std::min(dst_n, src_n);
    cudaMemcpy((void*)std::get<0>(dst),
               Buffers::accum_buffer.data_ptr(),
               n * sizeof(float2), cudaMemcpyDeviceToHost);
}


int64_t engine_get_cur_num_splats() {
    return Buffers::cur_num_splats;
}


void engine_copy_render_to_host(
    TorchTensorView out_rgb,
    TorchTensorView out_depth,
    TorchTensorView out_Ts
) {
    auto& renders = Buffers::fwd_renders;
    auto& rgb = std::get<0>(renders);
    auto& depth = std::get<1>(renders);
    if (rgb.data_ptr() && std::get<0>(out_rgb) != 0) {
        cudaMemcpy((void*)std::get<0>(out_rgb), rgb.data_ptr(),
                   rgb.numel() * sizeof(float3), cudaMemcpyDeviceToHost);
    }
    if (depth.data_ptr() && std::get<0>(out_depth) != 0) {
        cudaMemcpy((void*)std::get<0>(out_depth), depth.data_ptr(),
                   depth.numel() * sizeof(float), cudaMemcpyDeviceToHost);
    }
    if (Buffers::fwd_render_Ts.data_ptr() && std::get<0>(out_Ts) != 0) {
        cudaMemcpy((void*)std::get<0>(out_Ts), Buffers::fwd_render_Ts.data_ptr(),
                   Buffers::fwd_render_Ts.numel() * sizeof(float), cudaMemcpyDeviceToHost);
    }
}


void engine_copy_splats_to_host(
    TorchTensorView means,
    TorchTensorView quats,
    TorchTensorView scales,
    TorchTensorView opacities,
    TorchTensorView features_dc,
    TorchTensorView features_sh
) {
    _dv_to_host(Buffers::world_means, means);
    _dv_to_host(Buffers::world_quats, quats);
    _dv_to_host(Buffers::world_scales, scales);
    _dv_to_host(Buffers::world_opacities, opacities);
    _dv_to_host(Buffers::world_features_dc, features_dc);
    // features_sh: DeviceTensor2D<float3>
    if (Buffers::world_features_sh.data_ptr() && std::get<0>(features_sh) != 0) {
        cudaMemcpy((void*)std::get<0>(features_sh), Buffers::world_features_sh.data_ptr(),
                   Buffers::world_features_sh.numel() * sizeof(float3), cudaMemcpyDeviceToHost);
    }
}


std::vector<std::tuple<std::string, size_t, size_t>> engine_get_pool_breakdown() {
    return DevicePool::global().getBreakdown();
}

size_t engine_get_scratch_bytes() {
    return DeviceScratch::global().capBytes();
}


// Common ground-truth upload + GPU-side type conversion.
//
// kind picks the right per-element conversion:
//   "rgb"    : uint8 -> float in [0, 1] (divide by 255)
//   "normal" : uint8 -> float in [-1, 1] (x/127.5 - 1)
//   "depth"  : uint16 -> float (cast only, no scaling — matches `.float()`)
// For all kinds, `elem_size == sizeof(float)` is treated as already-float and
// goes through the regular _hv_to_dt3d<T>() zero-copy-or-H2D path.
//
// The shape carrier type T (float / float3) determines the DeviceTensor3D
// view returned. Conversion always produces a contiguous float buffer.
template<typename T>
static DeviceTensor3D<T> _hv_to_dt3d_gt(
    const TorchTensorView& src_tv,
    const std::string& key,
    const std::string& kind
) {
    if (std::get<0>(src_tv) == 0) return DeviceTensor3D<T>();
    uint32_t elem_size = std::get<1>(src_tv);
    if (elem_size == 4) {
        // Already float: existing zero-copy / pageable H2D path.
        return _hv_to_dt3d<T>(key, src_tv);
    }

    auto& shape = std::get<2>(src_tv);
    if (shape.size() != 4)
        throw std::runtime_error(key + ": expected [B, H, W, C] shape");
    int64_t B = shape[0], H = shape[1], W = shape[2], C = shape[3];
    int64_t numel = B * H * W * C;
    uint64_t src_ptr = std::get<0>(src_tv);
    bool src_is_device = _is_device_ptr((void*)src_ptr);

    // Float output buffer (consumed downstream by loss / bilagrid kernels).
    float* d_f = DevicePool::global().acquire<float>(key, (size_t)numel);

    // Shared H2D staging slot, reused across gt_rgb / gt_normal / gt_depth.
    // These calls happen sequentially in set_training_data on the default
    // stream, so each conversion kernel reads the staging buffer before the
    // next H2D copy overwrites it. One slot replaces the previous per-key
    // stagings (gt.rgb_u8, gt.normal_u8, gt.depth_u16) and saves ~32 MB per
    // step on dome-scale images.
    if (elem_size == 1) {
        const uint8_t* d_u8;
        if (src_is_device) {
            d_u8 = (const uint8_t*)src_ptr;
        } else {
            uint8_t* staging = DevicePool::global().acquire<uint8_t>(
                "gt.staging_u8", (size_t)numel);
            cudaMemcpy(staging, (void*)src_ptr, numel * sizeof(uint8_t), cudaMemcpyHostToDevice);
            d_u8 = staging;
        }
        if (kind == "rgb") {
            uint8_image_to_float_raw(d_u8, d_f, (int)B, (int)H, (int)W, (int)C);
        } else if (kind == "normal") {
            uint8_normal_to_float_raw(d_u8, d_f, (int)B, (int)H, (int)W, (int)C);
        } else {
            throw std::runtime_error(key + ": uint8 not supported for kind '" + kind + "'");
        }
    } else if (elem_size == 2) {
        const uint16_t* d_u16;
        if (src_is_device) {
            d_u16 = (const uint16_t*)src_ptr;
        } else {
            uint16_t* staging = DevicePool::global().acquire<uint16_t>(
                "gt.staging_u16", (size_t)numel);
            cudaMemcpy(staging, (void*)src_ptr, numel * sizeof(uint16_t), cudaMemcpyHostToDevice);
            d_u16 = staging;
        }
        if (kind == "rgb") {
            uint16_image_to_float_raw(d_u16, d_f, (int)B, (int)H, (int)W, (int)C);
        } else if (kind == "depth") {
            uint16_depth_to_float_raw(d_u16, d_f, (int)B, (int)H, (int)W, (int)C);
        } else {
            throw std::runtime_error(key + ": uint16 not supported for kind '" + kind + "'");
        }
    } else {
        throw std::runtime_error(key + ": unsupported element_size (expected 1, 2, or 4)");
    }

    TorchTensorView tv((uint64_t)d_f, 4, {B, H, W, C});
    return DeviceTensor3D<T>(tv);
}

void set_training_data(
    TorchTensorView gt_rgb,
    TorchTensorView gt_depth,
    TorchTensorView gt_normal,
    TorchTensorView gt_alpha
) {
    Buffers::gt_rgb    = _hv_to_dt3d_gt<float3>(gt_rgb,    "gt.rgb",    "rgb");
    Buffers::gt_depth  = _hv_to_dt3d_gt<float>(gt_depth,   "gt.depth",  "depth");
    Buffers::gt_normal = _hv_to_dt3d_gt<float3>(gt_normal, "gt.normal", "normal");
    // gt_alpha: small bool/uint8 buffer (the external mask). No conversion;
    // the slang kernel reads bool per pixel. Drives Buffers::has_mask.
    Buffers::gt_alpha  = _hv_to_dt3d<bool>("gt.alpha", gt_alpha);
    Buffers::has_gt    = (std::get<0>(gt_rgb) != 0);
    Buffers::has_mask  = (std::get<0>(gt_alpha) != 0);
}


static inline TorchTensorView _tv_null() {
    return TorchTensorView(0, 4, {0});
}

static inline bool _tv_valid(const TorchTensorView& tv) {
    return std::get<0>(tv) != 0;
}

static inline TorchTensorView _pool_tv(const std::string& key, int64_t B, int64_t H, int64_t W, int64_t C) {
    float* p = DevicePool::global().acquire<float>(key, (size_t)(B * H * W * C));
    return TorchTensorView((uint64_t)p, 4, {B, H, W, C});
}

static inline TorchTensorView _pool_tv_zero(const std::string& key, int64_t B, int64_t H, int64_t W, int64_t C) {
    float* p = DevicePool::global().acquire<float>(key, (size_t)(B * H * W * C));
    cudaMemset(p, 0, B * H * W * C * sizeof(float));
    return TorchTensorView((uint64_t)p, 4, {B, H, W, C});
}

static inline TorchTensorView _pool_tv_1d(const std::string& key, int64_t N) {
    float* p = DevicePool::global().acquire<float>(key, (size_t)N);
    return TorchTensorView((uint64_t)p, 4, {N});
}

static inline TorchTensorView _pool_tv_1d_zero(const std::string& key, int64_t N) {
    float* p = DevicePool::global().acquire<float>(key, (size_t)N);
    cudaMemset(p, 0, N * sizeof(float));
    return TorchTensorView((uint64_t)p, 4, {N});
}

static inline TorchTensorView _pool_tv_2d(const std::string& key, int64_t N, int64_t C) {
    float* p = DevicePool::global().acquire<float>(key, (size_t)(N * C));
    return TorchTensorView((uint64_t)p, 4, {N, C});
}

static inline TorchTensorView _pool_tv_3d_f(const std::string& key, int64_t N, int64_t K, int64_t C) {
    float* p = DevicePool::global().acquire<float>(key, (size_t)(N * K * C));
    return TorchTensorView((uint64_t)p, 4, {N, K, C});
}

static inline TorchTensorView _pool_tv_3d_u8(const std::string& key, int64_t N, int64_t K, int64_t C) {
    uint8_t* p = DevicePool::global().acquire<uint8_t>(key, (size_t)(N * K * C));
    return TorchTensorView((uint64_t)p, 1, {N, K, C});
}

static void _zero_tv(const TorchTensorView& tv) {
    if (std::get<0>(tv) == 0) return;
    int64_t bytes = std::get<1>(tv);
    for (auto s : std::get<2>(tv)) bytes *= s;
    cudaMemset((void*)std::get<0>(tv), 0, bytes);
}

static void _alloc_grad_buffers() {
    int64_t N = Buffers::max_num_splats;
    int64_t K = Buffers::num_sh;
    Buffers::grad_means.resize("eng.v_means", N);
    Buffers::grad_quats.resize("eng.v_quats", N);
    Buffers::grad_scales.resize("eng.v_scales", N);
    Buffers::grad_opacities.resize("eng.v_opacities", N);
    Buffers::grad_features_dc.resize("eng.v_features_dc", N);
    Buffers::grad_features_sh.resize("eng.v_features_sh", N, K);
    Buffers::grad_means.zero();
    Buffers::grad_quats.zero();
    Buffers::grad_scales.zero();
    Buffers::grad_opacities.zero();
    Buffers::grad_features_dc.zero();
    Buffers::grad_features_sh.zero();
}

static void _ensure_optim_state(bool quantize_sh, bool use_per_splat_bias_correction = false) {
    if (Buffers::optim_initialized && Buffers::train_quantize_sh == quantize_sh
        && Buffers::use_per_splat_bias_correction == use_per_splat_bias_correction)
        return;

    int64_t N = Buffers::max_num_splats;
    int64_t K = Buffers::num_sh;
    Buffers::train_quantize_sh = quantize_sh;

    // g1 (exp_avg)
    Buffers::g1_means.resize("eng.g1_means", N);
    Buffers::g1_quats.resize("eng.g1_quats", N);
    Buffers::g1_scales.resize("eng.g1_scales", N);
    Buffers::g1_opacities.resize("eng.g1_opacities", N);
    Buffers::g1_features_dc.resize("eng.g1_features_dc", N);
    if (!quantize_sh)
        Buffers::g1_features_sh.resize("eng.g1_features_sh", N, K);

    // g2 (exp_avg_sq)
    Buffers::g2_means.resize("eng.g2_means", N);
    Buffers::g2_quats.resize("eng.g2_quats", N);
    Buffers::g2_scales.resize("eng.g2_scales", N);
    Buffers::g2_opacities.resize("eng.g2_opacities", N);
    Buffers::g2_features_dc.resize("eng.g2_features_dc", N);
    if (!quantize_sh)
        Buffers::g2_features_sh.resize("eng.g2_features_sh", N, K);

    // radii [max_N]
    Buffers::radii.resize("eng.radii", N);

    // accum_buffer [max_N]
    Buffers::accum_buffer.resize("eng.accum_buffer", N);

    // Joint (u, sqrt(g2)) Adam state for SH features. The engine optim path
    // uses Optimizer.cu's fused_adam_with_steps_8bit_kernel which launches
    // one thread per cell (= per (splat, coef, channel)), grouped 256 threads
    // per CUDA block. So one float4 bounds slot covers 256 contiguous cells.
    if (quantize_sh) {
        constexpr int64_t BLOCK_SIZE = QuantizedAdamState<8, 256>::kBlockSize;
        int64_t sh_cells = (int64_t)N * K * 3;
        int64_t sh_bounds = (sh_cells + BLOCK_SIZE - 1) / BLOCK_SIZE;
        Buffers::sh_quant_state.resize("eng.sh_quant", sh_cells, sh_bounds);
    }

    // bias_correction_steps
    Buffers::use_per_splat_bias_correction = use_per_splat_bias_correction;
    if (use_per_splat_bias_correction) {
        Buffers::bias_correction_steps.resize("eng.bias_correction_steps", N);
    } else {
        Buffers::bias_correction_steps = DeviceVector<int32_t>();
    }

    // Zero everything on first init
    Buffers::g1_means.zero();       Buffers::g2_means.zero();
    Buffers::g1_quats.zero();       Buffers::g2_quats.zero();
    Buffers::g1_scales.zero();      Buffers::g2_scales.zero();
    Buffers::g1_opacities.zero();   Buffers::g2_opacities.zero();
    Buffers::g1_features_dc.zero(); Buffers::g2_features_dc.zero();
    if (quantize_sh) {
        Buffers::sh_quant_state.zero();
    } else {
        Buffers::g1_features_sh.zero();
        Buffers::g2_features_sh.zero();
    }
    Buffers::accum_buffer.zero();
    Buffers::bias_correction_steps.zero();

    Buffers::optim_initialized = true;
}


// Forward declarations for bilagrid helpers defined later in the file.
static void _ensure_bilagrid_optim_state();
static void _ensure_bilagrid_batch_grad(
    DeviceTensor5D<float>& grad,
    const DeviceTensor5D<float>& grids,
    int C_batch,
    const std::string& key);
static void _engine_bilagrid_backward_hook(
    TorchTensorView v_render_rgb,
    TorchTensorView v_ref_depth,
    TorchTensorView v_ref_normal);

// Forward declarations for PPISP helpers defined later in the file.
static void _ensure_ppisp_optim_state();
static void _engine_ppisp_backward_hook(TorchTensorView v_render_rgb);


std::map<std::string, float> engine_compute_loss_backward(
    int step,
    std::array<float, (int)LossWeightIndex::length> loss_weights,
    float w_ssim,
    int num_loss_scales,
    bool compute_loss_map
) {
    // Validate that forward was run
    if (std::get<0>(Buffers::fwd_renders) .data_ptr() == nullptr)
        throw std::runtime_error("engine_compute_loss_backward: forward_3dgs must be called first");
    if (!Buffers::has_gt)
        throw std::runtime_error("engine_compute_loss_backward: set_training_data must be called first");

    // Allocate and zero gradient buffers from pool
    _alloc_grad_buffers();

    int64_t C = Buffers::num_cameras;
    int64_t H = Buffers::image_height;
    int64_t W = Buffers::image_width;

    // Pool-allocate intermediates for loss computation
    TorchTensorView loss_map_buf = compute_loss_map ?
        _pool_tv_zero("eng.loss_map", C, H, W, 1) : _tv_null();

    // v_losses: constant vector [1, 0, 1, 1, ...] (1 for all, 0 for psnr slot)
    // Initialized once; pool never shrinks so pointer is stable
    static bool v_losses_initialized = false;
    TorchTensorView v_losses_buf = _pool_tv_1d("eng.v_losses", (int)LossIndex::length);
    if (!v_losses_initialized) {
        float h_v[(int)LossIndex::length];
        for (int i = 0; i < (int)LossIndex::length; i++) h_v[i] = 1.0f;
        h_v[(int)LossIndex::RgbPSNR] = 0.0f;
        cudaMemcpy((void*)std::get<0>(v_losses_buf), h_v, sizeof(h_v), cudaMemcpyHostToDevice);
        v_losses_initialized = true;
    }

    // Render outputs from forward pass (pool-backed, already populated)
    TorchTensorView render_rgb = TorchTensorView(
        (uint64_t)std::get<0>(Buffers::fwd_renders).data_ptr(), 4, {C, H, W, 3});
    TorchTensorView render_depth = TorchTensorView(
        (uint64_t)std::get<1>(Buffers::fwd_renders).data_ptr(), 4, {C, H, W, 1});
    TorchTensorView render_Ts = TorchTensorView(
        (uint64_t)Buffers::fwd_render_Ts.data_ptr(), 4, {C, H, W, 1});

    // Render normal: not yet available from forward (needs distortion output)
    // TODO: pass distortion outputs from forward when enabled
    TorchTensorView render_normal = _tv_null();
    TorchTensorView depth_normal = _tv_null();
    TorchTensorView rgb_dist = _tv_null();
    TorchTensorView depth_dist = _tv_null();
    TorchTensorView normal_dist = _tv_null();

    // Depth -> normal: derive depth_normal from rendered depth when gt_normal is provided
    // (matches training_losses.py logic: pred_normal is None, pred_depth exists, gt_normal exists).
    bool compute_depth_normal = (Buffers::gt_normal.data_ptr() != nullptr);
    bool is_ray_depth = (Buffers::primitive != "3dgs" && Buffers::primitive != "mip");
    if (compute_depth_normal) {
        depth_normal = _pool_tv("eng.depth_normal", C, H, W, 3);
        depth_to_normal_forward(
            Buffers::camera_model_str,
            _dv_tv(Buffers::camera_intrins),
            _dt2d_tv(Buffers::camera_dist_coeffs),
            is_ray_depth,
            DeviceTensor3D<float>(render_depth),
            DeviceTensor3D<float3>(depth_normal)
        );
    }

    std::vector<bool> needs_input_grad = {
        true,                                  // pred_rgb
        false,                                 // gt_rgb
        true,                                  // pred_depth
        Buffers::bilagrid_depth_enabled,       // gt_depth (true when bilagrid depth)
        true,                                  // pred_normal
        true,                                  // pred_depth_normal
        Buffers::bilagrid_normal_enabled,      // gt_normal (true when bilagrid normal)
        true,                                  // pred_transmittance
        true, true, true,                      // distortion
    };

    PerPixelGrads pixel_grads = {};

    // Allocate per-pixel gradient outputs
    pixel_grads.v_render_rgb  = _pool_tv("eng.v_rgb",   C, H, W, 3);
    pixel_grads.v_render_depth = _pool_tv("eng.v_depth", C, H, W, 1);
    pixel_grads.v_render_Ts   = _pool_tv("eng.v_Ts",    C, H, W, 1);
    if (compute_depth_normal) {
        pixel_grads.v_depth_normal = _pool_tv("eng.v_depth_normal", C, H, W, 3);
    }
    if (Buffers::bilagrid_depth_enabled) {
        pixel_grads.v_ref_depth = _pool_tv("eng.v_ref_depth", C, H, W, 1);
    }
    if (Buffers::bilagrid_normal_enabled) {
        pixel_grads.v_ref_normal = _pool_tv("eng.v_ref_normal", C, H, W, 3);
    }
    // TODO: allocate normal/distortion grads when those features are enabled

    // --- Compute per-pixel losses + SSIM, get gradients ---
    LossValues lv = compute_multi_scale_per_pixel_losses(
        num_loss_scales,
        render_rgb,
        _dt3d_tv(Buffers::gt_rgb),
        render_depth,
        _dt3d_tv(Buffers::gt_depth),
        render_normal,
        depth_normal,
        _dt3d_tv(Buffers::gt_normal),
        render_Ts,
        rgb_dist,
        depth_dist,
        normal_dist,
        _dt3d_tv(Buffers::gt_alpha),
        Buffers::has_mask,
        loss_weights,
        w_ssim,
        v_losses_buf,
        needs_input_grad,
        -1,  // num_train_images: -1 means use batch size
        _tv_null(),  // camera_indices: null means identity
        loss_map_buf,
        pixel_grads
    );

    // --- PPISP backward hook ---
    // Forward order is render -> bilagrid -> PPISP -> loss, so PPISP backward
    // runs FIRST: it rewrites v_render_rgb (post-PPISP -> pre-PPISP) and
    // accumulates the per-camera PPISP parameter gradient. The bilagrid hook
    // below then consumes the pre-PPISP v_render_rgb.
    if (Buffers::ppisp_enabled) {
        _ensure_ppisp_optim_state();
        _engine_ppisp_backward_hook(pixel_grads.v_render_rgb);
    }

    // --- Bilagrid backward hook ---
    // Transforms v_render_rgb (post-bilagrid -> pre-bilagrid) and accumulates
    // gradients into bilagrid_*_grads for any enabled bilagrid types. The
    // updated v_render_rgb then flows into rasterization backward as usual.
    // For depth / normal (gt-side), the loss-side gradient is consumed here
    // and only the bilagrid grid gradient is retained.
    if (Buffers::bilagrid_rgb_enabled || Buffers::bilagrid_depth_enabled ||
        Buffers::bilagrid_normal_enabled) {
        _ensure_bilagrid_optim_state();
        _engine_bilagrid_backward_hook(
            pixel_grads.v_render_rgb,
            pixel_grads.v_ref_depth,
            pixel_grads.v_ref_normal);
    }

    // TODO: color space conversion (rgb_to_srgb) forward/backward
    // TODO: background blending forward/backward

    // Depth -> normal backward: propagate v_depth_normal grads into v_render_depth (in-place add)
    if (compute_depth_normal) {
        depth_to_normal_backward(
            Buffers::camera_model_str,
            _dv_tv(Buffers::camera_intrins),
            _dt2d_tv(Buffers::camera_dist_coeffs),
            is_ray_depth,
            DeviceTensor3D<float>(render_depth),
            DeviceTensor3D<float3>(pixel_grads.v_depth_normal),
            DeviceTensor3D<float>(pixel_grads.v_render_depth)
        );
    }

    // --- Rasterization + projection backward ---
    // Build v_render_outputs from pixel_grads
    RenderOutput::TensorTuple v_render_outputs = std::make_tuple(
        DeviceTensor3D<float3>(pixel_grads.v_render_rgb),
        DeviceTensor3D<float>(pixel_grads.v_render_depth),
        DeviceTensor3D<float3>()  // no normal gradient yet
    );
    DeviceTensor3D<float> v_render_Ts(pixel_grads.v_render_Ts);

    DeviceTensorFloatND v_fnd_opac;
    {
        TorchTensorView opac_tv((uint64_t)Buffers::grad_opacities.data_ptr(), 4,
            {Buffers::grad_opacities.size(), 1LL, 1LL});
        v_fnd_opac = DeviceTensorFloatND(opac_tv);  // ndim=2, shape=[N, 1]
    }
    std::vector<DeviceTensorFloatND> v_splats_w = {
        DeviceTensorFloatND(Buffers::grad_means),         // [N, 3]
        DeviceTensorFloatND(Buffers::grad_quats),         // [N, 4]
        DeviceTensorFloatND(Buffers::grad_scales),        // [N, 3]
        v_fnd_opac,                                       // [N, 1]
        DeviceTensorFloatND(Buffers::grad_features_dc),   // [N, 3]
        DeviceTensorFloatND(Buffers::grad_features_sh),   // [N, K, 3]
    };

    // Build accum_weight_map from loss_map (pixel-space -> per-splat mapping in raster bwd)
    DeviceTensor3D<float> accum_weight_map;
    if (compute_loss_map && _tv_valid(loss_map_buf)) {
        accum_weight_map = DeviceTensor3D<float>(loss_map_buf);
    }

    std::vector<DeviceTensorFloatND> v_splats_w_out, v_splats_s_out;

    if (Buffers::primitive == "3dgs" || Buffers::primitive == "mip") {
        auto [vw, vs, accum_weight] = rasterize_to_pixels_3dgs_bwd(
            Buffers::cur_num_splats,
            Buffers::fwd_splats_w,
            Buffers::fwd_splats_s,
            Buffers::fwd_gaussian_ids,
            (uint32_t)Buffers::image_width,
            (uint32_t)Buffers::image_height,
            Buffers::fwd_tile_offsets,
            Buffers::fwd_flatten_ids,
            Buffers::fwd_render_Ts,
            Buffers::fwd_last_ids,
            Buffers::fwd_renders,
            accum_weight_map,
            v_render_outputs,
            v_render_Ts,
            std::make_optional(v_splats_w),
            std::nullopt
        );
        v_splats_w_out = vw;
        v_splats_s_out = vs;
        if (accum_weight.data_ptr()) Buffers::fwd_accum_weight = accum_weight;
    } else if (Buffers::primitive == "3dgut") {
        auto [vw, vs, vviewmats, accum_weight] = rasterize_to_pixels_3dgut_bwd(
            Buffers::cur_num_splats,
            Buffers::fwd_splats_w,
            Buffers::fwd_splats_s,
            Buffers::fwd_gaussian_ids,
            _dt2d_tv(Buffers::camera_viewmats),
            _dv_tv(Buffers::camera_intrins),
            Buffers::camera_model_str,
            _dt2d_tv(Buffers::camera_dist_coeffs),
            (uint32_t)Buffers::image_width,
            (uint32_t)Buffers::image_height,
            Buffers::fwd_tile_offsets,
            Buffers::fwd_flatten_ids,
            Buffers::fwd_render_Ts,
            Buffers::fwd_last_ids,
            Buffers::fwd_renders,
            std::nullopt,  // render2_outputs
            DeviceTensor3D<float>(),  // loss_map
            accum_weight_map,
            v_render_outputs,
            v_render_Ts,
            std::nullopt,  // v_distortion_outputs
            std::make_optional(v_splats_w),
            std::nullopt,
            false
        );
        v_splats_w_out = vw;
        v_splats_s_out = vs;
        if (accum_weight.data_ptr()) Buffers::fwd_accum_weight = accum_weight;
    } else {
        throw std::runtime_error("engine_compute_loss_backward: unknown primitive: " + Buffers::primitive);
    }

    // --- Projection backward ---
    if (Buffers::primitive == "3dgs") {
        projection_3dgs_backward(
            Buffers::cur_num_splats, Buffers::sh_degree, Buffers::fwd_splats_w,
            _dt2d_tv(Buffers::camera_viewmats), _dv_tv(Buffers::camera_intrins),
            (uint32_t)Buffers::image_width, (uint32_t)Buffers::image_height,
            Buffers::camera_model_str, _dt2d_tv(Buffers::camera_dist_coeffs),
            Buffers::fwd_camera_ids, Buffers::fwd_gaussian_ids,
            Buffers::fwd_aabb, v_splats_s_out, v_splats_w_out, nullptr);
    } else if (Buffers::primitive == "mip") {
        projection_mip_backward(
            Buffers::cur_num_splats, Buffers::sh_degree, Buffers::fwd_splats_w,
            _dt2d_tv(Buffers::camera_viewmats), _dv_tv(Buffers::camera_intrins),
            (uint32_t)Buffers::image_width, (uint32_t)Buffers::image_height,
            Buffers::camera_model_str, _dt2d_tv(Buffers::camera_dist_coeffs),
            Buffers::fwd_camera_ids, Buffers::fwd_gaussian_ids,
            Buffers::fwd_aabb, v_splats_s_out, v_splats_w_out, nullptr);
    } else if (Buffers::primitive == "3dgut") {
        projection_3dgut_backward(
            Buffers::cur_num_splats, Buffers::sh_degree, Buffers::fwd_splats_w,
            _dt2d_tv(Buffers::camera_viewmats), _dv_tv(Buffers::camera_intrins),
            (uint32_t)Buffers::image_width, (uint32_t)Buffers::image_height,
            Buffers::camera_model_str, _dt2d_tv(Buffers::camera_dist_coeffs),
            Buffers::fwd_camera_ids, Buffers::fwd_gaussian_ids,
            Buffers::fwd_aabb, v_splats_s_out, v_splats_w_out, nullptr);
    }

    // --- Build loss dict for display ---
    auto sdiv = [](float x, float y) -> float { return y != 0.0f ? x / y : 0.0f; };

    std::map<std::string, float> loss_dict;
    loss_dict["rgb_loss"] = lv.rgb_loss + w_ssim * (1.0f - lv.ssim);
    loss_dict["psnr"] = lv.rgb_psnr;
    loss_dict["ssim"] = lv.ssim;
    loss_dict["depth_loss"] = sdiv(lv.depth_sup,
        loss_weights[(int)LossWeightIndex::DepthSup]);
    loss_dict["normal_loss"] = sdiv(lv.normal_sup,
        loss_weights[(int)LossWeightIndex::NormalSup]);
    loss_dict["alpha_loss"] = sdiv(lv.alpha_sup,
        loss_weights[(int)LossWeightIndex::AlphaSup] + loss_weights[(int)LossWeightIndex::AlphaSupUnder]);
    loss_dict["normal_reg"] = sdiv(lv.normal_reg,
        loss_weights[(int)LossWeightIndex::NormalReg]);
    loss_dict["alpha_reg"] = sdiv(lv.alpha_reg,
        loss_weights[(int)LossWeightIndex::AlphaReg]);
    loss_dict["rgb_dist_reg"] = sdiv(lv.rgb_dist_reg,
        loss_weights[(int)LossWeightIndex::RgbDistReg]);
    loss_dict["depth_dist_reg"] = sdiv(lv.depth_dist_reg,
        loss_weights[(int)LossWeightIndex::DepthDistReg]);
    loss_dict["normal_dist_reg"] = sdiv(lv.normal_dist_reg,
        loss_weights[(int)LossWeightIndex::NormalDistReg]);

    return loss_dict;
}


// ============================================================================
// Bilagrid integration
// ============================================================================

// Pool keys are scoped under "eng.bg.*" so the pool reports them as engine state.
static constexpr cudaStream_t kBilagridStream = (cudaStream_t)0;

void engine_init_bilagrid_rgb(int n_grids, std::string type, int L, int H, int W,
                              bool quantize_optim) {
    int C;
    if (type == "affine") C = 12;
    else if (type == "ppisp" || type == "loglinear") C = 9;
    else throw std::runtime_error("engine_init_bilagrid_rgb: unknown type " + type);
    if (n_grids <= 0)
        throw std::runtime_error("engine_init_bilagrid_rgb: n_grids must be > 0");

    Buffers::bilagrid_rgb_type = type;
    Buffers::bilagrid_rgb_C = C;
    Buffers::bilagrid_rgb_quantize_optim = quantize_optim;
    Buffers::bilagrid_rgb_grids.resize("eng.bg.rgb.grids", n_grids, L, H, W, C);
    if (type == "affine") {
        Buffers::bilagrid_rgb_grids.zero();
        bilagrid_affine_identity_init(
            Buffers::bilagrid_rgb_grids.data_ptr(), n_grids, L, H, W, kBilagridStream);
    } else {
        Buffers::bilagrid_rgb_grids.zero();
    }
    Buffers::bilagrid_rgb_enabled = true;
    Buffers::bilagrid_rgb_optim_initialized = false;
}

void engine_init_bilagrid_depth(int n_grids, int L, int H, int W,
                                bool quantize_optim) {
    if (n_grids <= 0)
        throw std::runtime_error("engine_init_bilagrid_depth: n_grids must be > 0");
    Buffers::bilagrid_depth_quantize_optim = quantize_optim;
    Buffers::bilagrid_depth_grids.resize("eng.bg.depth.grids", n_grids, L, H, W, 2);
    Buffers::bilagrid_depth_grids.zero();
    Buffers::bilagrid_depth_scalars.resize("eng.bg.depth.scalars", n_grids);
    Buffers::bilagrid_depth_scalars.zero();
    Buffers::bilagrid_depth_enabled = true;
    Buffers::bilagrid_depth_optim_initialized = false;
}

void engine_init_bilagrid_normal(int n_grids, int L, int H, int W,
                                 bool quantize_optim) {
    if (n_grids <= 0)
        throw std::runtime_error("engine_init_bilagrid_normal: n_grids must be > 0");
    Buffers::bilagrid_normal_quantize_optim = quantize_optim;
    Buffers::bilagrid_normal_grids.resize("eng.bg.normal.grids", n_grids, L, H, W, 3);
    Buffers::bilagrid_normal_grids.zero();
    Buffers::bilagrid_normal_enabled = true;
    Buffers::bilagrid_normal_optim_initialized = false;
}

// Compute TV losses for each enabled bilagrid into a pool-backed 3-float
// device buffer at [rgb_tv, depth_tv, normal_tv]; unset entries stay 0.
static void _engine_bilagrid_tv_into(float* tv_buf3_device) {
    cudaMemsetAsync(tv_buf3_device, 0, 3 * sizeof(float), kBilagridStream);
    auto run = [&](DeviceTensor5D<float>& grids, int slot) {
        if (grids.data_ptr() == nullptr) return;
        int N = (int)grids.size<0>();
        int L = (int)grids.size<1>(), H = (int)grids.size<2>(), W = (int)grids.size<3>();
        int C = (int)grids.size<4>();
        tv_loss_forward(grids.data_ptr(), tv_buf3_device + slot,
                        N, C, L, H, W, kBilagridStream);
    };
    if (Buffers::bilagrid_rgb_enabled) run(Buffers::bilagrid_rgb_grids, 0);
    if (Buffers::bilagrid_depth_enabled) run(Buffers::bilagrid_depth_grids, 1);
    if (Buffers::bilagrid_normal_enabled) run(Buffers::bilagrid_normal_grids, 2);
}

// Populate Buffers::bilagrid_cur_cam_indices from a host or device int32
// tensor. Empty/null tensor -> leave it null (kernels fall back to identity).
// Called by engine_bilagrid_forward / engine_ppisp_forward, and directly by
// engine_train_step before the forwards.
static void _set_cur_cam_indices(TorchTensorView tv) {
    if (std::get<0>(tv) == 0 || std::get<2>(tv).empty()) {
        Buffers::bilagrid_cur_cam_indices = DeviceVector<int32_t>();
    } else {
        Buffers::bilagrid_cur_cam_indices = _hv_to_dv<int32_t>(
            "eng.bg.cam_indices", tv);
    }
}

// Apply bilagrid forward in-place for each enabled type for the current batch.
// Saves a pre-bilagrid copy used by the backward pass.
//
// The kernel uses indirect grid indexing: each image in the batch (index `ni`
// inside the kernel) reads from the grid slot `cam_indices[ni]`. This allows
// a multi-image batch with different per-image cam_idx values to be processed
// in a single kernel launch without any gather of grid params. When
// `cam_indices` is empty/null, kernels fall back to identity (ni -> ni).
void engine_bilagrid_forward(TorchTensorView cam_indices) {
    _set_cur_cam_indices(cam_indices);
    int H = Buffers::image_height;
    int W = Buffers::image_width;
    int C_batch = (int)Buffers::num_cameras;
    if (C_batch <= 0) return;

    const int* cam_idx_dev = Buffers::bilagrid_cur_cam_indices.data_ptr();

    // --- RGB: applies to Buffers::fwd_renders.rgb ---
    if (Buffers::bilagrid_rgb_enabled) {
        if (std::get<0>(Buffers::fwd_renders).data_ptr() == nullptr)
            throw std::runtime_error("engine_bilagrid_forward: forward_3dgs must run first");

        // Save pre-bilagrid render and apply in place.
        Buffers::fwd_rgb_pre_bilagrid.resize("eng.bg.rgb.pre", C_batch, H, W);
        float3* render_rgb = std::get<0>(Buffers::fwd_renders).data_ptr();
        cudaMemcpyAsync(
            Buffers::fwd_rgb_pre_bilagrid.data_ptr(),
            render_rgb,
            (size_t)C_batch * H * W * sizeof(float3),
            cudaMemcpyDeviceToDevice, kBilagridStream);

        int L = (int)Buffers::bilagrid_rgb_grids.size<1>();
        int gH = (int)Buffers::bilagrid_rgb_grids.size<2>();
        int gW = (int)Buffers::bilagrid_rgb_grids.size<3>();
        // Pass the FULL grid table; kernel does the per-image indirect lookup.
        float* grid_ptr = Buffers::bilagrid_rgb_grids.data_ptr();

        // In-place: rgb input and output are the same buffer (the kernel reads
        // and writes index-by-index per pixel, so this is safe).
        if (Buffers::bilagrid_rgb_type == "affine") {
            bilagrid_uniform_sample_forward(
                grid_ptr, (const float*)render_rgb, (float*)render_rgb,
                /*N=*/C_batch, L, gH, gW, /*m=*/1, H, W,
                kBilagridStream, cam_idx_dev);
        } else if (Buffers::bilagrid_rgb_type == "ppisp") {
            bilagrid_ppisp_uniform_sample_forward(
                grid_ptr, (const float*)render_rgb, (float*)render_rgb,
                /*N=*/C_batch, L, gH, gW, /*m=*/1, H, W,
                kBilagridStream, cam_idx_dev);
        } else if (Buffers::bilagrid_rgb_type == "loglinear") {
            bilagrid_loglinear_uniform_sample_forward(
                grid_ptr, (const float*)render_rgb, (float*)render_rgb,
                /*N=*/C_batch, L, gH, gW, /*m=*/1, H, W,
                kBilagridStream, cam_idx_dev);
        }
    }

    // --- Depth (gt side) ---
    if (Buffers::bilagrid_depth_enabled && Buffers::gt_depth.data_ptr() != nullptr) {
        Buffers::fwd_depth_pre_bilagrid.resize("eng.bg.depth.pre", C_batch, H, W);
        float* gt_depth_ptr = Buffers::gt_depth.data_ptr();
        cudaMemcpyAsync(
            Buffers::fwd_depth_pre_bilagrid.data_ptr(),
            gt_depth_ptr,
            (size_t)C_batch * H * W * sizeof(float),
            cudaMemcpyDeviceToDevice, kBilagridStream);

        // Compute per-image scalar (median quantile) over pre-bilagrid depth
        // and scatter into bilagrid_depth_scalars at the cam_indices slots.
        // compute_depth_scalars_tensor writes scalars[batch_i] = quantile_i —
        // we want to land them at scalars[cam_indices[batch_i]] in the full
        // table. Simplest: compute into a temp [C_batch] buffer then scatter.
        float* scalar_full = Buffers::bilagrid_depth_scalars.data_ptr();
        float* tmp_scalars = DevicePool::global().acquire<float>(
            "eng.bg.depth.tmp_scalars", (size_t)C_batch);
        TorchTensorView depth_tv((uint64_t)Buffers::fwd_depth_pre_bilagrid.data_ptr(),
            4, {C_batch, H, W, 1, 1});
        TorchTensorView scalar_tv((uint64_t)tmp_scalars, 4, {C_batch});
        compute_depth_scalars_tensor(depth_tv, /*patched=*/true, scalar_tv);
        // Scatter tmp_scalars -> scalar_full at cam_indices.
        if (cam_idx_dev != nullptr) {
            bilagrid_scatter_floats(tmp_scalars, cam_idx_dev, C_batch,
                                    scalar_full, kBilagridStream);
        } else {
            // Identity: tmp_scalars[i] -> scalar_full[i].
            cudaMemcpyAsync(scalar_full, tmp_scalars,
                C_batch * sizeof(float), cudaMemcpyDeviceToDevice, kBilagridStream);
        }

        int L = (int)Buffers::bilagrid_depth_grids.size<1>();
        int gH = (int)Buffers::bilagrid_depth_grids.size<2>();
        int gW = (int)Buffers::bilagrid_depth_grids.size<3>();
        float* grid_ptr = Buffers::bilagrid_depth_grids.data_ptr();
        bilagrid_depth_uniform_sample_forward(
            grid_ptr,
            (const float*)Buffers::fwd_depth_pre_bilagrid.data_ptr(),
            scalar_full,  // kernel reads scalars[cam_indices[ni]]
            gt_depth_ptr,
            /*N=*/C_batch, L, gH, gW, /*m=*/1, H, W,
            kBilagridStream, cam_idx_dev);
    }

    // --- Normal (gt side) ---
    if (Buffers::bilagrid_normal_enabled && Buffers::gt_normal.data_ptr() != nullptr) {
        Buffers::fwd_normal_pre_bilagrid.resize("eng.bg.normal.pre", C_batch, H, W);
        float3* gt_normal_ptr = Buffers::gt_normal.data_ptr();
        cudaMemcpyAsync(
            Buffers::fwd_normal_pre_bilagrid.data_ptr(),
            gt_normal_ptr,
            (size_t)C_batch * H * W * sizeof(float3),
            cudaMemcpyDeviceToDevice, kBilagridStream);

        int L = (int)Buffers::bilagrid_normal_grids.size<1>();
        int gH = (int)Buffers::bilagrid_normal_grids.size<2>();
        int gW = (int)Buffers::bilagrid_normal_grids.size<3>();
        float* grid_ptr = Buffers::bilagrid_normal_grids.data_ptr();
        bilagrid_normal_uniform_sample_forward(
            grid_ptr,
            (const float*)Buffers::fwd_normal_pre_bilagrid.data_ptr(),
            (float*)gt_normal_ptr,
            /*N=*/C_batch, L, gH, gW, /*m=*/1, H, W,
            kBilagridStream, cam_idx_dev);
    }
}

// Bilagrid backward hook (Phase B): image-loss gradient now writes a
// SPARSE per-batch buffer (`eng.bg.*.image_grad` sized [C_batch, C, L, H, W])
// instead of the full N_grids-wide dense table. Reads of the bilagrid grid
// parameters still go through cam_indices for the true camera lookup. The
// per-batch buffer is zeroed up front so the kernel's atomicAdds start clean.
// The TV-loss contribution is no longer materialized here — it folds into the
// fused optimizer step (BilagridFusedAdam.cu).
//
// For RGB the kernel additionally overwrites v_render_rgb (post-bilagrid ->
// pre-bilagrid) so it can flow into rasterization backward.
static void _engine_bilagrid_backward_hook(
    TorchTensorView v_render_rgb,   // [C, H, W, 3] — overwritten when RGB enabled
    TorchTensorView v_ref_depth,    // [C, H, W, 1] — input only
    TorchTensorView v_ref_normal    // [C, H, W, 3] — input only
) {
    int H = Buffers::image_height;
    int W = Buffers::image_width;
    int C_batch = (int)Buffers::num_cameras;
    const int* cam_idx_dev = Buffers::bilagrid_cur_cam_indices.data_ptr();

    // --- RGB backward ---
    if (Buffers::bilagrid_rgb_enabled && std::get<0>(v_render_rgb) != 0 &&
        Buffers::fwd_rgb_pre_bilagrid.data_ptr() != nullptr)
    {
        int L = (int)Buffers::bilagrid_rgb_grids.size<1>();
        int gH = (int)Buffers::bilagrid_rgb_grids.size<2>();
        int gW = (int)Buffers::bilagrid_rgb_grids.size<3>();
        float* grid_ptr = Buffers::bilagrid_rgb_grids.data_ptr();
        _ensure_bilagrid_batch_grad(Buffers::bilagrid_rgb_image_grad,
                                    Buffers::bilagrid_rgb_grids, C_batch,
                                    "eng.bg.rgb.image_grad");
        float* grad_grid_ptr = Buffers::bilagrid_rgb_image_grad.data_ptr();

        // v_render_rgb is the gradient w.r.t. POST-bilagrid rgb (what entered loss).
        // The backward overwrites it with the gradient w.r.t. PRE-bilagrid rgb.
        float* v_rgb_ptr = (float*)std::get<0>(v_render_rgb);
        const float* rgb_pre_ptr = (const float*)Buffers::fwd_rgb_pre_bilagrid.data_ptr();

        if (Buffers::bilagrid_rgb_type == "affine") {
            bilagrid_uniform_sample_backward_v1(
                grid_ptr, rgb_pre_ptr, v_rgb_ptr,
                grad_grid_ptr, v_rgb_ptr,
                /*N=*/C_batch, L, gH, gW, /*m=*/1, H, W,
                /*block_x=*/8, /*block_y=*/8, /*target_tile_size=*/5,
                kBilagridStream, cam_idx_dev);
        } else if (Buffers::bilagrid_rgb_type == "ppisp") {
            bilagrid_ppisp_uniform_sample_backward_v1(
                grid_ptr, rgb_pre_ptr, v_rgb_ptr,
                grad_grid_ptr, v_rgb_ptr,
                /*N=*/C_batch, L, gH, gW, /*m=*/1, H, W,
                8, 8, 5, kBilagridStream, cam_idx_dev);
        } else if (Buffers::bilagrid_rgb_type == "loglinear") {
            bilagrid_loglinear_uniform_sample_backward_v1(
                grid_ptr, rgb_pre_ptr, v_rgb_ptr,
                grad_grid_ptr, v_rgb_ptr,
                /*N=*/C_batch, L, gH, gW, /*m=*/1, H, W,
                8, 8, 5, kBilagridStream, cam_idx_dev);
        }
    }

    // --- Depth backward (accumulate into bilagrid grad; discard input grad) ---
    // We pass v_depth=nullptr so the bilagrid backward skips the kernel that
    // would compute the ∂L/∂pre-bilagrid (= ∂L/∂gt_depth) buffer. That kernel's
    // output is never consumed downstream (GT isn't a parameter), and skipping
    // it eliminates a full-resolution C*H*W*4 byte scratch + one kernel launch.
    if (Buffers::bilagrid_depth_enabled && std::get<0>(v_ref_depth) != 0 &&
        Buffers::fwd_depth_pre_bilagrid.data_ptr() != nullptr)
    {
        int L = (int)Buffers::bilagrid_depth_grids.size<1>();
        int gH = (int)Buffers::bilagrid_depth_grids.size<2>();
        int gW = (int)Buffers::bilagrid_depth_grids.size<3>();
        float* grid_ptr = Buffers::bilagrid_depth_grids.data_ptr();
        _ensure_bilagrid_batch_grad(Buffers::bilagrid_depth_image_grad,
                                    Buffers::bilagrid_depth_grids, C_batch,
                                    "eng.bg.depth.image_grad");
        float* grad_grid_ptr = Buffers::bilagrid_depth_image_grad.data_ptr();
        float* scalars = Buffers::bilagrid_depth_scalars.data_ptr();
        bilagrid_depth_uniform_sample_backward_v1(
            grid_ptr,
            (const float*)Buffers::fwd_depth_pre_bilagrid.data_ptr(),
            scalars,
            (const float*)std::get<0>(v_ref_depth),
            grad_grid_ptr, /*v_depth=*/nullptr,
            /*N=*/C_batch, L, gH, gW, /*m=*/1, H, W,
            8, 8, 5, kBilagridStream, cam_idx_dev);
    }

    // --- Normal backward (accumulate; discard input grad) ---
    // Same reasoning as depth: v_rgb=nullptr skips the unused gt-side grad
    // kernel, saving a C*H*W*12 byte scratch + one kernel launch.
    if (Buffers::bilagrid_normal_enabled && std::get<0>(v_ref_normal) != 0 &&
        Buffers::fwd_normal_pre_bilagrid.data_ptr() != nullptr)
    {
        int L = (int)Buffers::bilagrid_normal_grids.size<1>();
        int gH = (int)Buffers::bilagrid_normal_grids.size<2>();
        int gW = (int)Buffers::bilagrid_normal_grids.size<3>();
        float* grid_ptr = Buffers::bilagrid_normal_grids.data_ptr();
        _ensure_bilagrid_batch_grad(Buffers::bilagrid_normal_image_grad,
                                    Buffers::bilagrid_normal_grids, C_batch,
                                    "eng.bg.normal.image_grad");
        float* grad_grid_ptr = Buffers::bilagrid_normal_image_grad.data_ptr();
        bilagrid_normal_uniform_sample_backward_v1(
            grid_ptr,
            (const float*)Buffers::fwd_normal_pre_bilagrid.data_ptr(),
            (const float*)std::get<0>(v_ref_normal),
            grad_grid_ptr, /*v_rgb=*/nullptr,
            /*N=*/C_batch, L, gH, gW, /*m=*/1, H, W,
            8, 8, 5, kBilagridStream, cam_idx_dev);
    }
}


// Ensure bilagrid Adam moment buffers are allocated and zeroed for each
// enabled bilagrid type. When quantize_optim is true, g1/g2 storage is uint8
// (1 byte/cell) plus a float4 per 256-cell block holding the {g1 min, g1 max,
// g2 min, g2 max} ranges — same scheme as the SH optimizer state.
//
// Note: the per-batch image_grad buffer is sized + zeroed per-iter in
// _ensure_bilagrid_batch_grad (called from _engine_bilagrid_backward_hook),
// since its size depends on the current batch size.
static void _ensure_bilagrid_optim_state() {
    auto setup = [](DeviceTensor5D<float>& grids,
                    DeviceTensor5D<float>& g1,                 // used when !quantize
                    DeviceTensor5D<float>& g2,                 // used when !quantize
                    QuantizedAdamState<8, 256>& quant_state,   // used when quantize
                    bool quantize_optim,
                    const std::string& key_prefix,
                    bool& done) {
        if (done || grids.data_ptr() == nullptr) return;
        int64_t N = grids.size<0>();
        int64_t L = grids.size<1>(), H = grids.size<2>(), W = grids.size<3>();
        int64_t C = grids.size<4>();
        int64_t numel = N * L * H * W * C;
        if (quantize_optim) {
            constexpr int64_t BLOCK_SIZE = QuantizedAdamState<8, 256>::kBlockSize;
            int64_t n_blocks = (numel + BLOCK_SIZE - 1) / BLOCK_SIZE;
            quant_state.resize(key_prefix + ".bg_quant", numel, n_blocks);
            quant_state.zero();
        } else {
            g1.resize(key_prefix + ".g1", N, L, H, W, C);
            g1.zero();
            g2.resize(key_prefix + ".g2", N, L, H, W, C);
            g2.zero();
        }
        done = true;
    };
    if (Buffers::bilagrid_rgb_enabled) {
        setup(Buffers::bilagrid_rgb_grids,
              Buffers::bilagrid_rgb_g1, Buffers::bilagrid_rgb_g2,
              Buffers::bilagrid_rgb_quant_state,
              Buffers::bilagrid_rgb_quantize_optim,
              "eng.bg.rgb", Buffers::bilagrid_rgb_optim_initialized);
    }
    if (Buffers::bilagrid_depth_enabled) {
        setup(Buffers::bilagrid_depth_grids,
              Buffers::bilagrid_depth_g1, Buffers::bilagrid_depth_g2,
              Buffers::bilagrid_depth_quant_state,
              Buffers::bilagrid_depth_quantize_optim,
              "eng.bg.depth", Buffers::bilagrid_depth_optim_initialized);
    }
    if (Buffers::bilagrid_normal_enabled) {
        setup(Buffers::bilagrid_normal_grids,
              Buffers::bilagrid_normal_g1, Buffers::bilagrid_normal_g2,
              Buffers::bilagrid_normal_quant_state,
              Buffers::bilagrid_normal_quantize_optim,
              "eng.bg.normal", Buffers::bilagrid_normal_optim_initialized);
    }
}

// Allocate + zero a per-batch image_grad buffer for one bilagrid type. Called
// at the top of _engine_bilagrid_backward_hook so the image-loss backward
// kernel can atomicAdd into a clean buffer. Sized [C_batch, C, L, H, W].
static void _ensure_bilagrid_batch_grad(
    DeviceTensor5D<float>& grad,
    const DeviceTensor5D<float>& grids,
    int C_batch,
    const std::string& key
) {
    if (grids.data_ptr() == nullptr) return;
    int L = (int)grids.size<1>();
    int H = (int)grids.size<2>();
    int W = (int)grids.size<3>();
    int C = (int)grids.size<4>();
    grad.resize(key, C_batch, L, H, W, C);
    grad.zero();
}

// Helper: build a DeviceTensorFloatND from a DeviceTensor5D<float> by viewing
// it as a flat 1D float tensor of `numel` elements. This is what fused_adam_step
// expects (it only cares about the float count).
static DeviceTensorFloatND _bg_flat(const DeviceTensor5D<float>& t) {
    if (t.data_ptr() == nullptr) return DeviceTensorFloatND();
    int64_t numel = t.numel();
    TorchTensorView tv((uint64_t)t.data_ptr(), 4, {numel, 1LL});
    return DeviceTensorFloatND(tv);
}

void engine_bilagrid_optim_step(
    int step,
    float lr_rgb, float lr_depth, float lr_normal,
    float tv_weight_rgb, float tv_weight_depth, float tv_weight_normal
) {
    _ensure_bilagrid_optim_state();
    int32_t adam_step = step + 1;
    int C_batch = (int)Buffers::num_cameras;
    const int* cam_idx_dev = Buffers::bilagrid_cur_cam_indices.data_ptr();

    // Single fused pass per bilagrid type: reads the sparse per-batch
    // image_grad, computes the local TV gradient from `grids` neighbors, and
    // applies Adam (float or uint8-quantized) in one kernel. The previous
    // separate tv_loss_backward + zeroing of a dense grad table is gone —
    // TV is inlined, image_grad is the only persistent gradient state and
    // it gets re-zeroed at the start of the *next* iter's backward hook.
    auto run = [&](DeviceTensor5D<float>& grids,
                   DeviceTensor5D<float>& image_grad,
                   DeviceTensor5D<float>& g1,
                   DeviceTensor5D<float>& g2,
                   QuantizedAdamState<8, 256>& quant_state,
                   bool quantize_optim,
                   float lr, float tv_weight) {
        if (grids.data_ptr() == nullptr || lr <= 0.0f) return;
        if (image_grad.data_ptr() == nullptr) return;
        int N = (int)grids.size<0>();
        int L = (int)grids.size<1>(), H = (int)grids.size<2>(), W = (int)grids.size<3>();
        int C = (int)grids.size<4>();
        fused_bilagrid_tv_adam(
            grids.data_ptr(),
            quantize_optim ? nullptr : g1.data_ptr(),
            quantize_optim ? nullptr : g2.data_ptr(),
            quantize_optim ? quant_state.packed_ptr() : nullptr,
            quantize_optim ? quant_state.bounds_ptr() : nullptr,
            image_grad.data_ptr(),
            cam_idx_dev,
            N, C_batch, C, L, H, W,
            lr, tv_weight, adam_step,
            quantize_optim, kBilagridStream);
    };

    run(Buffers::bilagrid_rgb_grids, Buffers::bilagrid_rgb_image_grad,
        Buffers::bilagrid_rgb_g1, Buffers::bilagrid_rgb_g2,
        Buffers::bilagrid_rgb_quant_state, Buffers::bilagrid_rgb_quantize_optim,
        lr_rgb, tv_weight_rgb);
    run(Buffers::bilagrid_depth_grids, Buffers::bilagrid_depth_image_grad,
        Buffers::bilagrid_depth_g1, Buffers::bilagrid_depth_g2,
        Buffers::bilagrid_depth_quant_state, Buffers::bilagrid_depth_quantize_optim,
        lr_depth, tv_weight_depth);
    run(Buffers::bilagrid_normal_grids, Buffers::bilagrid_normal_image_grad,
        Buffers::bilagrid_normal_g1, Buffers::bilagrid_normal_g2,
        Buffers::bilagrid_normal_quant_state, Buffers::bilagrid_normal_quantize_optim,
        lr_normal, tv_weight_normal);
}


// ============================================================================
// PPISP integration
// ============================================================================

static constexpr cudaStream_t kPpispStream = (cudaStream_t)0;

// Build a TorchTensorView for the full [N_cam, P] PPISP parameter table.
static TorchTensorView _ppisp_params_full_tv(DeviceTensor2D<float>& table) {
    int64_t N = table.size<0>(), P = table.size<1>();
    return TorchTensorView(
        (uint64_t)table.data_ptr(), 4, {N, P});
}

// Build a TorchTensorView for the current batch's camera_intrins ([C_batch, 4]).
static TorchTensorView _ppisp_intrins_batch_tv() {
    int C_batch = (int)Buffers::num_cameras;
    return TorchTensorView(
        (uint64_t)Buffers::camera_intrins.data_ptr(), 4, {(int64_t)C_batch, 4LL});
}

// Build a TorchTensorView wrapper around the current batch's cam_indices
// DeviceVector, or null when no cam_indices were supplied (identity mode).
static TorchTensorView _ppisp_cam_indices_tv() {
    const int* ptr = Buffers::bilagrid_cur_cam_indices.data_ptr();
    if (ptr == nullptr) return TorchTensorView(0, 4, {});
    int64_t n = Buffers::bilagrid_cur_cam_indices.size();
    return TorchTensorView((uint64_t)ptr, 4, {n});
}

void engine_init_ppisp(int n_grids, std::string param_type) {
    if (n_grids <= 0)
        throw std::runtime_error("engine_init_ppisp: n_grids must be > 0");
    int P;
    if (param_type == "original" || param_type == "") P = 36;
    else if (param_type == "rqs") P = 39;
    else throw std::runtime_error(
        "engine_init_ppisp: unknown param_type \"" + param_type +
        "\", must be \"original\" or \"rqs\"");

    Buffers::ppisp_param_type = (param_type == "" ? std::string("original") : param_type);
    Buffers::ppisp_num_params = P;
    Buffers::ppisp_params.resize("eng.ppisp.params", n_grids, P);
    if (Buffers::ppisp_param_type == "original") {
        ppisp_original_default_init(
            Buffers::ppisp_params.data_ptr(), n_grids, kPpispStream);
    } else {
        Buffers::ppisp_params.zero();
    }
    Buffers::ppisp_enabled = true;
    Buffers::ppisp_optim_initialized = false;
}

// Apply PPISP forward in place on the current rendered RGB. Reads the
// per-image PPISP parameter slot via Buffers::bilagrid_cur_cam_indices (shared
// with bilagrid). When that index buffer is empty, fall back to identity.
void engine_ppisp_forward(TorchTensorView cam_indices) {
    if (!Buffers::ppisp_enabled) return;
    // Repopulating with the same tensor is cheap (DeviceVector<int32_t> takes
    // a view of device memory when the source is on device).
    _set_cur_cam_indices(cam_indices);

    int H = Buffers::image_height;
    int W = Buffers::image_width;
    int C_batch = (int)Buffers::num_cameras;
    if (C_batch <= 0) return;
    if (std::get<0>(Buffers::fwd_renders).data_ptr() == nullptr)
        throw std::runtime_error("engine_ppisp_forward: forward_3dgs must run first");

    // Save pre-PPISP render.
    Buffers::fwd_rgb_pre_ppisp.resize("eng.ppisp.rgb_pre", C_batch, H, W);
    float3* render_rgb = std::get<0>(Buffers::fwd_renders).data_ptr();
    cudaMemcpyAsync(
        Buffers::fwd_rgb_pre_ppisp.data_ptr(),
        render_rgb,
        (size_t)C_batch * H * W * sizeof(float3),
        cudaMemcpyDeviceToDevice, kPpispStream);

    // Build views: in_image and out_image are the same buffer (kernel is
    // pixel-local, so in-place is safe).
    TorchTensorView rgb_tv((uint64_t)render_rgb, 4, {(int64_t)C_batch, H, W, 3});
    DeviceTensor3D<float3> in_img(rgb_tv);
    DeviceTensor3D<float3> out_img(rgb_tv);

    ppisp_forward(
        in_img,
        _ppisp_params_full_tv(Buffers::ppisp_params),
        _ppisp_intrins_batch_tv(),
        (float)W, (float)H,
        Buffers::ppisp_param_type,
        _ppisp_cam_indices_tv(),
        out_img
    );
}

// Ensure PPISP Adam moment and gradient buffers are allocated + zeroed.
static void _ensure_ppisp_optim_state() {
    if (Buffers::ppisp_optim_initialized) return;
    if (Buffers::ppisp_params.data_ptr() == nullptr) return;
    int64_t N = Buffers::ppisp_params.size<0>();
    int64_t P = Buffers::ppisp_params.size<1>();
    Buffers::ppisp_grads.resize("eng.ppisp.grads", N, P);
    Buffers::ppisp_grads.zero();
    Buffers::ppisp_g1.resize("eng.ppisp.g1", N, P);
    Buffers::ppisp_g1.zero();
    Buffers::ppisp_g2.resize("eng.ppisp.g2", N, P);
    Buffers::ppisp_g2.zero();
    Buffers::ppisp_optim_initialized = true;
}

// Backward hook: rewrites v_render_rgb (post-PPISP -> pre-PPISP) in place and
// accumulates the per-camera PPISP parameter gradient into ppisp_grads. The
// kernel writes via atomicAdd, so duplicate cam_indices entries are accumulated
// correctly into the same per-camera slot.
static void _engine_ppisp_backward_hook(TorchTensorView v_render_rgb) {
    int H = Buffers::image_height;
    int W = Buffers::image_width;
    int C_batch = (int)Buffers::num_cameras;
    if (std::get<0>(v_render_rgb) == 0) return;
    if (Buffers::fwd_rgb_pre_ppisp.data_ptr() == nullptr) return;

    TorchTensorView pre_tv(
        (uint64_t)Buffers::fwd_rgb_pre_ppisp.data_ptr(), 4,
        {(int64_t)C_batch, H, W, 3});
    DeviceTensor3D<float3> in_img(pre_tv);
    DeviceTensor3D<float3> v_out_img(v_render_rgb);
    DeviceTensor3D<float3> v_in_img(v_render_rgb);  // overwrite in place

    ppisp_backward(
        in_img,
        _ppisp_params_full_tv(Buffers::ppisp_params),
        _ppisp_intrins_batch_tv(),
        (float)W, (float)H,
        v_out_img,
        Buffers::ppisp_param_type,
        _ppisp_cam_indices_tv(),
        v_in_img,
        _ppisp_params_full_tv(Buffers::ppisp_grads)
    );
}

// Compute the 6 PPISP regularization losses, scaled by their loss weights.
// Returns the [6]-float device buffer (pool-backed). When `compute_grad=true`,
// also accumulates the per-camera parameter gradient into Buffers::ppisp_grads
// using v_losses = ones (matching training_losses.py).
static float* _engine_ppisp_reg_loss_into(
    const std::array<float, (int)PPISPRegLossIndex::length>& loss_weights,
    bool compute_grad
) {
    int N = (int)Buffers::ppisp_params.size<0>();
    int kRaw = (Buffers::ppisp_param_type == "rqs")
        ? (int)RawPPISPRegLossIndexRQS::length
        : (int)RawPPISPRegLossIndex::length;
    int kLoss = (int)PPISPRegLossIndex::length;

    // Output losses (zeroed each call so the in-kernel write is a clean store).
    float* losses_buf = DevicePool::global().acquire<float>("eng.ppisp.reg_losses", kLoss);
    cudaMemsetAsync(losses_buf, 0, kLoss * sizeof(float), kPpispStream);

    // Raw losses scratch ([N+1, kRaw], pre-zeroed).
    float* raw_losses_buf = DevicePool::global().acquire<float>(
        "eng.ppisp.reg_raw_losses", (size_t)(N + 1) * kRaw);
    cudaMemsetAsync(raw_losses_buf, 0, (size_t)(N + 1) * kRaw * sizeof(float),
                    kPpispStream);

    TorchTensorView params_tv(
        (uint64_t)Buffers::ppisp_params.data_ptr(), 4,
        {(int64_t)N, (int64_t)Buffers::ppisp_num_params});
    TorchTensorView losses_tv((uint64_t)losses_buf, 4, {(int64_t)kLoss});
    TorchTensorView raw_tv((uint64_t)raw_losses_buf, 4,
        {(int64_t)(N + 1), (int64_t)kRaw});

    compute_ppsip_regularization_forward(
        params_tv, loss_weights, Buffers::ppisp_param_type,
        losses_tv, raw_tv);

    if (compute_grad) {
        // v_losses = ones[kLoss]: gradient flows back through reg-loss sum.
        float* v_losses = DevicePool::global().acquire<float>(
            "eng.ppisp.v_reg_losses", kLoss);
        std::vector<float> h_ones(kLoss, 1.0f);
        cudaMemcpyAsync(v_losses, h_ones.data(), kLoss * sizeof(float),
                        cudaMemcpyHostToDevice, kPpispStream);
        TorchTensorView v_losses_tv((uint64_t)v_losses, 4, {(int64_t)kLoss});

        // Accumulate into ppisp_grads (the regularization backward writes a
        // fresh tensor, so use a scratch buffer and add in afterwards).
        float* v_params_scratch = DevicePool::global().acquire<float>(
            "eng.ppisp.v_reg_params",
            (size_t)N * Buffers::ppisp_num_params);
        cudaMemsetAsync(v_params_scratch, 0,
            (size_t)N * Buffers::ppisp_num_params * sizeof(float), kPpispStream);
        TorchTensorView v_params_tv(
            (uint64_t)v_params_scratch, 4,
            {(int64_t)N, (int64_t)Buffers::ppisp_num_params});

        compute_ppsip_regularization_backward(
            params_tv, loss_weights, raw_tv, v_losses_tv,
            Buffers::ppisp_param_type, v_params_tv);

        // ppisp_grads += v_params_scratch (over all N * P floats).
        size_t total = (size_t)N * Buffers::ppisp_num_params;
        ppisp_add_into_grad(v_params_scratch,
            Buffers::ppisp_grads.data_ptr(), total, kPpispStream);
    }

    return losses_buf;
}

void engine_ppisp_optim_step(
    int step,
    float lr,
    float reg_exposure_mean,
    float reg_vig_center,
    float reg_vig_non_pos,
    float reg_vig_channel_var,
    float reg_color_mean,
    float reg_crf_channel_var
) {
    if (!Buffers::ppisp_enabled) return;
    _ensure_ppisp_optim_state();
    if (lr <= 0.0f) {
        // Still zero accumulated grads for next iteration to match other paths.
        Buffers::ppisp_grads.zero();
        return;
    }

    // Fold the 6 regularization losses' gradients into ppisp_grads.
    std::array<float, (int)PPISPRegLossIndex::length> loss_weights = {
        reg_exposure_mean,
        reg_vig_center,
        reg_vig_non_pos,
        reg_vig_channel_var,
        reg_color_mean,
        reg_crf_channel_var
    };
    bool any_reg = false;
    for (float w : loss_weights) if (w > 0.0f) { any_reg = true; break; }
    if (any_reg) {
        (void)_engine_ppisp_reg_loss_into(loss_weights, /*compute_grad=*/true);
    }

    int64_t N = Buffers::ppisp_params.size<0>();
    int64_t P = Buffers::ppisp_params.size<1>();
    int64_t numel = N * P;
    auto flat_view = [numel](float* p) {
        TorchTensorView tv((uint64_t)p, 4, {numel, 1LL});
        return DeviceTensorFloatND(tv);
    };
    DeviceVector<int32_t> no_per_splat_steps;
    fused_adam_step(
        numel,
        flat_view(Buffers::ppisp_params.data_ptr()),
        flat_view(Buffers::ppisp_grads.data_ptr()),
        flat_view(Buffers::ppisp_g1.data_ptr()),
        flat_view(Buffers::ppisp_g2.data_ptr()),
        lr, step + 1, no_per_splat_steps,
        /*l2_reg=*/0.0f, /*l2_reg_offset=*/0.0f);
    Buffers::ppisp_grads.zero();
}


void engine_optim_step(
    int step,
    float lr_means, float lr_quats, float lr_scales, float lr_opacities,
    float lr_features_dc, float lr_features_sh,
    float max_gauss_ratio, float scale_regularization_weight,
    float mcmc_opacity_reg_weight, float mcmc_scale_reg_weight,
    float erank_reg_weight, float erank_reg_weight_s3,
    float quat_norm_reg_weight, float sh_reg_weight,
    bool use_scale_agnostic_mean, bool quantize_sh,
    bool use_per_splat_bias_correction
) {
    _ensure_optim_state(quantize_sh, use_per_splat_bias_correction);

    int64_t N = Buffers::cur_num_splats;

    // Per-splat bias correction: increment all steps by 1 before use
    DeviceVector<int32_t> per_splat_steps;
    if (use_per_splat_bias_correction) {
        increment_int32_inplace(Buffers::bias_correction_steps, Buffers::max_num_splats);
        per_splat_steps = Buffers::bias_correction_steps;
    }

    fused_optim_3dgs_geometry(
        N,
        Buffers::world_means,     Buffers::grad_means,     Buffers::g1_means,     Buffers::g2_means,
        Buffers::world_quats,     Buffers::grad_quats,     Buffers::g1_quats,     Buffers::g2_quats,
        Buffers::world_scales,    Buffers::grad_scales,    Buffers::g1_scales,    Buffers::g2_scales,
        Buffers::world_opacities, Buffers::grad_opacities, Buffers::g1_opacities, Buffers::g2_opacities,
        Buffers::radii,
        lr_means, lr_quats, lr_scales, lr_opacities,
        max_gauss_ratio, scale_regularization_weight,
        mcmc_opacity_reg_weight, mcmc_scale_reg_weight,
        erank_reg_weight, erank_reg_weight_s3, quat_norm_reg_weight,
        use_scale_agnostic_mean,
        step + 1, per_splat_steps
    );

    fused_adam_step(N,
        DeviceTensorFloatND(Buffers::world_features_dc),
        DeviceTensorFloatND(Buffers::grad_features_dc),
        DeviceTensorFloatND(Buffers::g1_features_dc),
        DeviceTensorFloatND(Buffers::g2_features_dc),
        lr_features_dc, step + 1, per_splat_steps,
        sh_reg_weight, 0.5f / 0.28209479177387814f);

    if (quantize_sh && Buffers::sh_quant_state.initialized()) {
        fused_adam_step_8bit(N,
            DeviceTensorFloatND(Buffers::world_features_sh),
            DeviceTensorFloatND(Buffers::grad_features_sh),
            Buffers::sh_quant_state.packed_ptr(),
            Buffers::sh_quant_state.bounds_ptr(),
            lr_features_sh, step + 1, per_splat_steps,
            sh_reg_weight, 0.0f);
    } else {
        fused_adam_step(N,
            DeviceTensorFloatND(Buffers::world_features_sh),
            DeviceTensorFloatND(Buffers::grad_features_sh),
            DeviceTensorFloatND(Buffers::g1_features_sh),
            DeviceTensorFloatND(Buffers::g2_features_sh),
            lr_features_sh, step + 1, per_splat_steps,
            sh_reg_weight, 0.0f);
    }
}


// Helper: build DeviceVector<T> from TorchTensorView (non-owning view)
template<typename T>
static inline DeviceVector<T> _dv(const TorchTensorView& tv) {
    return DeviceVector<T>(tv);
}
// For null-able TorchTensorViews, return default (null) DeviceVector
template<typename T>
static inline DeviceVector<T> _dv_or_null(const TorchTensorView& tv) {
    if (std::get<0>(tv) == 0) return DeviceVector<T>();
    return DeviceVector<T>(tv);
}
// For [N, K, C] tensors like features_sh: flatten to [N*K, C] then build DeviceVector
template<typename T>
static inline DeviceVector<T> _dv_flat(const TorchTensorView& tv) {
    if (std::get<0>(tv) == 0) return DeviceVector<T>();
    const auto& shape = std::get<2>(tv);
    int64_t flat_count = 1;
    for (size_t i = 0; i + 1 < shape.size(); i++) flat_count *= shape[i];
    TorchTensorView flat_tv(std::get<0>(tv), std::get<1>(tv),
        {flat_count, (int64_t)(sizeof(T) / std::get<1>(tv))});
    return DeviceVector<T>(flat_tv);
}


int engine_densify_step(
    int step,
    int max_steps,
    int refine_start_iter,
    int refine_stop_num_iter,
    int refine_every,
    float growth_factor,
    float min_opacity,
    float max_screen_size,
    float max_screen_size_clip_hardness,
    float max_world_size,
    float noise_lr,
    float noise_lr_final,
    float relocate_heuristic_weight
) {
    // _ensure_optim_state(Buffers::train_quantize_sh);

    int64_t cur_num_splats = Buffers::cur_num_splats;
    int64_t max_num_splats = Buffers::max_num_splats;
    bool quantize_sh = Buffers::train_quantize_sh;
    int num_sh = Buffers::num_sh;

    bool use_revised = (relocate_heuristic_weight >= 1.0f);
    bool densify_ongoing = (step < max_steps - refine_stop_num_iter);
    bool do_densify = densify_ongoing && (step > refine_start_iter && step % refine_every == 0);
    float progress = ((float)step + 0.5f) / (float)max_steps;

    // Use pool-backed DeviceVector/DeviceTensor from Buffers directly
    auto& dv_means = Buffers::world_means;
    auto& dv_quats = Buffers::world_quats;
    auto& dv_scales = Buffers::world_scales;
    auto& dv_opacs = Buffers::world_opacities;
    auto& dv_features_dc = Buffers::world_features_dc;
    // features_sh: DeviceTensor2D<float3> [N, K] -> flatten to DeviceVector<float3> [N*K]
    DeviceVector<float3> dv_features_sh;
    {
        auto& t = Buffers::world_features_sh;
        TorchTensorView tv((uint64_t)t.data_ptr(), 4,
            {t.size<0>() * t.size<1>(), 3LL});
        dv_features_sh = DeviceVector<float3>(tv);
    }
    auto& dv_g1_means = Buffers::g1_means;
    auto& dv_g1_quats = Buffers::g1_quats;
    auto& dv_g1_scales = Buffers::g1_scales;
    auto& dv_g1_opacs = Buffers::g1_opacities;
    auto& dv_g1_features_dc = Buffers::g1_features_dc;
    auto& dv_g2_means = Buffers::g2_means;
    auto& dv_g2_quats = Buffers::g2_quats;
    auto& dv_g2_scales = Buffers::g2_scales;
    auto& dv_g2_opacs = Buffers::g2_opacities;
    auto& dv_g2_features_dc = Buffers::g2_features_dc;
    // g1/g2 features_sh as DeviceVector<float3> with numel=N*K. The densify kernel
    // reinterprets the data as uchar3 when quantize_sh, so the element_size carrier
    // type is just for the constructor's shape check.
    DeviceVector<float3> dv_g1_features_sh, dv_g2_features_sh;
    int64_t sh_flat = (int64_t)Buffers::max_num_splats * (int64_t)Buffers::num_sh;
    if (quantize_sh) {
        // Joint (u, sqrt(g2)) AoS layout: bytes for both halves of every cell
        // sit in a single packed buffer. Densify only zeros bytes on dst
        // splats, so we pass the same packed pointer for both the g1 and g2
        // wrappers; the kernel template is instantiated with `short3` (6 bytes
        // per SH coef = 3 channels x 2 bytes), so each write zeros the joint
        // (u, sqrt_g2) record for one (splat, coef). The second write
        // (g2_features_sh[i] = 0) targets the same memory and is a no-op.
        uint8_t* packed = Buffers::sh_quant_state.packed_ptr();
        TorchTensorView tv((uint64_t)packed, 1, {sh_flat, 12LL});
        dv_g1_features_sh = DeviceVector<float3>(tv);
        dv_g2_features_sh = DeviceVector<float3>(tv);
    } else {
        // Non-quantized: storage is float3 N*K elements (12 bytes each), shape [N*K, 3], element_size=4
        TorchTensorView tv1((uint64_t)Buffers::g1_features_sh.data_ptr(), 4, {sh_flat, 3LL});
        TorchTensorView tv2((uint64_t)Buffers::g2_features_sh.data_ptr(), 4, {sh_flat, 3LL});
        dv_g1_features_sh = DeviceVector<float3>(tv1);
        dv_g2_features_sh = DeviceVector<float3>(tv2);
    }
    auto& dv_radii = Buffers::radii;
    auto& dv_accum_buf = Buffers::accum_buffer;
    auto& dv_bias_steps = Buffers::bias_correction_steps;

    // Clip large splats
    if (std::isfinite(max_screen_size) || std::isfinite(max_world_size)) {
        densify_clip_scale_tensor(
            cur_num_splats,
            dv_radii, dv_scales,
            nullptr,  // don't clip opacities
            max_screen_size,
            max_screen_size_clip_hardness,
            max_world_size
        );
    }

    // Update densification score using accum_weight from rasterization backward
    if (densify_ongoing && use_revised && dv_accum_buf.data_ptr() != nullptr) {
        DeviceVector<float> score;
        if (Buffers::fwd_accum_weight.data_ptr() != nullptr) {
            score = Buffers::fwd_accum_weight;
        } else {
            score = Buffers::grad_opacities.data_ptr() ? Buffers::grad_opacities : dv_opacs;
        }

        densify_update_weight_tensor(
            cur_num_splats,
            dv_radii,
            nullptr,
            (float*)dv_opacs.data_ptr(),
            score,
            dv_accum_buf,
            false
        );
    }

    int num_added = 0;

    if (do_densify && use_revised) {
        // Revised relocation (long axis split)
        relocate_splats_with_long_axis_split_tensor(
            cur_num_splats, min_opacity,
            dv_means, dv_quats, dv_scales, dv_opacs, dv_features_dc, dv_features_sh,
            dv_g1_means, dv_g1_quats, dv_g1_scales, dv_g1_opacs, dv_g1_features_dc, dv_g1_features_sh,
            dv_g2_means, dv_g2_quats, dv_g2_scales, dv_g2_opacs, dv_g2_features_dc, dv_g2_features_sh,
            dv_accum_buf, dv_bias_steps,
            quantize_sh, num_sh, 2 * step + 0
        );

        // Add more splats
        int64_t n_target = std::min(max_num_splats, (int64_t)(growth_factor * cur_num_splats));
        num_added = (int)std::max((int64_t)0, n_target - cur_num_splats);
        if (num_added > 0) {
            add_splats_with_long_axis_split_tensor(
                cur_num_splats, num_added,
                dv_means, dv_quats, dv_scales, dv_opacs, dv_features_dc, dv_features_sh,
                dv_g1_means, dv_g1_quats, dv_g1_scales, dv_g1_opacs, dv_g1_features_dc, dv_g1_features_sh,
                dv_g2_means, dv_g2_quats, dv_g2_scales, dv_g2_opacs, dv_g2_features_dc, dv_g2_features_sh,
                dv_accum_buf, dv_bias_steps,
                quantize_sh, num_sh, 2 * step + 1
            );
        }
    } else if (do_densify) {
        // MCMC relocation
        relocate_splats_mcmc_tensor(
            cur_num_splats, min_opacity,
            dv_means, dv_quats, dv_scales, dv_opacs, dv_features_dc, dv_features_sh,
            dv_g1_means, dv_g1_quats, dv_g1_scales, dv_g1_opacs, dv_g1_features_dc, dv_g1_features_sh,
            dv_g2_means, dv_g2_quats, dv_g2_scales, dv_g2_opacs, dv_g2_features_dc, dv_g2_features_sh,
            dv_bias_steps,
            quantize_sh, num_sh, 2 * step + 0
        );

        // MCMC sample add
        int64_t n_target = std::min(max_num_splats, (int64_t)(growth_factor * cur_num_splats));
        num_added = (int)std::max((int64_t)0, n_target - cur_num_splats);
        if (num_added > 0) {
            add_splats_mcmc_tensor(
                cur_num_splats, num_added, min_opacity,
                dv_means, dv_quats, dv_scales, dv_opacs, dv_features_dc, dv_features_sh,
                dv_g1_means, dv_g1_quats, dv_g1_scales, dv_g1_opacs, dv_g1_features_dc, dv_g1_features_sh,
                dv_g2_means, dv_g2_quats, dv_g2_scales, dv_g2_opacs, dv_g2_features_dc, dv_g2_features_sh,
                dv_bias_steps,
                quantize_sh, num_sh, 2 * step + 1
            );
        }
    }

    // Reset accum buffer after densification step
    if (do_densify && dv_accum_buf.data_ptr() != nullptr) {
        cudaMemset(dv_accum_buf.data_ptr(), 0, dv_accum_buf.size() * sizeof(float2));
    }

    // Add MCMC noise
    if (noise_lr > 0.0f && noise_lr_final > 0.0f) {
        float noise_scalar = noise_lr * powf(noise_lr_final / noise_lr, progress);
        if (use_revised) {
            revised_add_noise_tensor(
                cur_num_splats + num_added, noise_scalar,
                dv_radii, dv_means, dv_scales, dv_quats, dv_opacs);
        } else {
            mcmc_add_noise_tensor(
                cur_num_splats + num_added, noise_scalar,
                dv_means, dv_scales, dv_quats, dv_opacs);
        }
    }

    Buffers::cur_num_splats = cur_num_splats + num_added;
    return num_added;
}


// ============================================================
// engine_train_step — single fused training step
// Combines set_camera_params + set_training_data + forward_3dgs +
//          compute_loss_backward + optim_step + densify_step
// Returns loss_dict (+ num_splats, num_added) for verbose.
// All tensor inputs are host (CPU); H->D happens inside; no D->H of large data.
// ============================================================

std::map<std::string, float> engine_train_step(
    int step, int max_steps,
    // Forward config
    std::string primitive,
    int sh_degree,
    bool packed,
    // Camera params (host)
    int width, int height, std::string camera_model,
    TorchTensorView viewmats,
    TorchTensorView intrins,
    TorchTensorView dist_coeffs,
    // GT data (host). Per-pixel masks are derived inside the slang kernel from
    // gt_depth (depth/alpha masks), gt_normal (normal mask sentinel), and
    // gt_alpha (the external mask -> RGB mask). apply_loss_for_mask is folded
    // into alpha_loss_weight / alpha_loss_weight_under on the Python side.
    TorchTensorView gt_rgb,
    TorchTensorView gt_depth,
    TorchTensorView gt_normal,
    TorchTensorView gt_alpha,
    // Loss config
    std::array<float, (int)LossWeightIndex::length> loss_weights,
    float w_ssim,
    int num_loss_scales,
    bool compute_loss_map,
    // Optimizer LRs
    float lr_means, float lr_quats, float lr_scales, float lr_opacities,
    float lr_features_dc, float lr_features_sh,
    // Optimizer regularization
    float max_gauss_ratio, float scale_regularization_weight,
    float mcmc_opacity_reg_weight, float mcmc_scale_reg_weight,
    float erank_reg_weight, float erank_reg_weight_s3,
    float quat_norm_reg_weight, float sh_reg_weight,
    bool use_scale_agnostic_mean,
    bool quantize_sh,
    bool use_per_splat_bias_correction,
    // Densify config
    int refine_start_iter, int refine_stop_num_iter, int refine_every,
    float growth_factor, float min_opacity,
    float max_screen_size, float max_screen_size_clip_hardness,
    float max_world_size,
    float noise_lr, float noise_lr_final,
    float relocate_heuristic_weight,
    // Bilagrid config
    TorchTensorView bilagrid_cam_indices,
    float bilagrid_lr_rgb, float bilagrid_lr_depth, float bilagrid_lr_normal,
    float bilagrid_tv_weight_rgb, float bilagrid_tv_weight_depth, float bilagrid_tv_weight_normal,
    // PPISP config (RGB only; cam_idx reuses bilagrid_cam_idx)
    float ppisp_lr,
    float ppisp_reg_exposure_mean,
    float ppisp_reg_vig_center,
    float ppisp_reg_vig_non_pos,
    float ppisp_reg_vig_channel_var,
    float ppisp_reg_color_mean,
    float ppisp_reg_crf_channel_var
) {
    // Camera + GT: H->D copy into pool
    set_camera_params(width, height, camera_model, viewmats, intrins, dist_coeffs);
    set_training_data(gt_rgb, gt_depth, gt_normal, gt_alpha);

    // Forward (writes to pool; no D->H)
    forward_3dgs(primitive, sh_degree, packed);

    // Bilagrid forward (between rendering and loss). No-op when disabled.
    // Both forwards populate the shared cam_indices buffer; the H->D copy of
    // a [C_batch] int32 tensor is in the microseconds.
    if (Buffers::bilagrid_rgb_enabled || Buffers::bilagrid_depth_enabled ||
        Buffers::bilagrid_normal_enabled) {
        engine_bilagrid_forward(bilagrid_cam_indices);
    }

    // PPISP forward (in place on rendered RGB, AFTER bilagrid). No-op when
    // disabled.
    if (Buffers::ppisp_enabled) {
        engine_ppisp_forward(bilagrid_cam_indices);
    }

    // Loss + backward (reads pool, writes pool grads; D->H only for scalar loss values)
    std::map<std::string, float> loss_dict = engine_compute_loss_backward(
        step, loss_weights, w_ssim, num_loss_scales, compute_loss_map);

    // Optimizer (in-place on pool buffers, no copies)
    engine_optim_step(
        step,
        lr_means, lr_quats, lr_scales, lr_opacities,
        lr_features_dc, lr_features_sh,
        max_gauss_ratio, scale_regularization_weight,
        mcmc_opacity_reg_weight, mcmc_scale_reg_weight,
        erank_reg_weight, erank_reg_weight_s3,
        quat_norm_reg_weight, sh_reg_weight,
        use_scale_agnostic_mean, quantize_sh,
        use_per_splat_bias_correction);

    // Bilagrid Adam step + TV regularization (after splat optim, before densify).
    if (Buffers::bilagrid_rgb_enabled || Buffers::bilagrid_depth_enabled ||
        Buffers::bilagrid_normal_enabled) {
        engine_bilagrid_optim_step(
            step,
            bilagrid_lr_rgb, bilagrid_lr_depth, bilagrid_lr_normal,
            bilagrid_tv_weight_rgb, bilagrid_tv_weight_depth, bilagrid_tv_weight_normal);

        // Read post-update TV losses for verbose display. Async pattern:
        // read previous iter's slot now (cheap event-sync), queue this iter's
        // D->H, never block the host on cudaMemcpy.
        float* tv_buf = DevicePool::global().acquire<float>("eng.bg.tv_readout", 3);
        _engine_bilagrid_tv_into(tv_buf);
        static AsyncReadout<float> tv_readout(3);
        const float* h_tv = tv_readout.read_previous();
        if (h_tv) {
            if (Buffers::bilagrid_rgb_enabled)    loss_dict["bilagrid_tv"]        = h_tv[0];
            if (Buffers::bilagrid_depth_enabled)  loss_dict["bilagrid_depth_tv"]  = h_tv[1];
            if (Buffers::bilagrid_normal_enabled) loss_dict["bilagrid_normal_tv"] = h_tv[2];
        }
        tv_readout.issue(tv_buf);
    }

    // PPISP Adam step + regularization (after splat optim, before densify).
    // Also read post-update regularization losses for verbose display.
    if (Buffers::ppisp_enabled) {
        engine_ppisp_optim_step(
            step, ppisp_lr,
            ppisp_reg_exposure_mean, ppisp_reg_vig_center,
            ppisp_reg_vig_non_pos, ppisp_reg_vig_channel_var,
            ppisp_reg_color_mean, ppisp_reg_crf_channel_var);

        std::array<float, (int)PPISPRegLossIndex::length> reg_w = {
            ppisp_reg_exposure_mean, ppisp_reg_vig_center,
            ppisp_reg_vig_non_pos, ppisp_reg_vig_channel_var,
            ppisp_reg_color_mean, ppisp_reg_crf_channel_var
        };
        float* losses_buf = _engine_ppisp_reg_loss_into(reg_w, /*compute_grad=*/false);
        // Async readout: previous iter's slot now, queue this iter's D->H.
        static AsyncReadout<float> ppisp_reg_readout((int)PPISPRegLossIndex::length);
        const float* h_losses = ppisp_reg_readout.read_previous();
        if (h_losses) {
            loss_dict["ppisp_reg_exposure_mean"]   = h_losses[(int)PPISPRegLossIndex::ExposureMean];
            loss_dict["ppisp_reg_vig_center"]      = h_losses[(int)PPISPRegLossIndex::VignettingCenter];
            loss_dict["ppisp_reg_vig_non_pos"]     = h_losses[(int)PPISPRegLossIndex::VignettingNonPositivity];
            loss_dict["ppisp_reg_vig_channel_var"] = h_losses[(int)PPISPRegLossIndex::VignettingChannelVariance];
            loss_dict["ppisp_reg_color_mean"]      = h_losses[(int)PPISPRegLossIndex::ColorMean];
            loss_dict["ppisp_reg_crf_channel_var"] = h_losses[(int)PPISPRegLossIndex::CRFChannelVariance];
        }
        ppisp_reg_readout.issue(losses_buf);
    }

    // Densify (in-place on pool buffers, no copies)
    int num_added = engine_densify_step(
        step, max_steps,
        refine_start_iter, refine_stop_num_iter, refine_every,
        growth_factor, min_opacity,
        max_screen_size, max_screen_size_clip_hardness, max_world_size,
        noise_lr, noise_lr_final, relocate_heuristic_weight);

    loss_dict["num_added"] = (float)num_added;
    loss_dict["cur_num_splats"] = (float)Buffers::cur_num_splats;
    loss_dict["max_num_splats"] = (float)Buffers::max_num_splats;
    return loss_dict;
}


// ============================================================================
// Checkpoint save (PLY + npy)
// ============================================================================

namespace {

namespace fs = std::filesystem;

static void _ckpt_mkdir(const fs::path& p) {
    std::error_code ec;
    fs::create_directories(p, ec);
    if (ec) throw std::runtime_error("create_directories failed: " + p.string() + " (" + ec.message() + ")");
}

// D->H copy into a freshly sized host buffer.
template<typename T>
static std::vector<T> _ckpt_d2h(const T* d_ptr, size_t n) {
    std::vector<T> host(n);
    if (n > 0 && d_ptr != nullptr) {
        cudaMemcpy(host.data(), d_ptr, n * sizeof(T), cudaMemcpyDeviceToHost);
    }
    return host;
}

// Save a device-side buffer of `numel` elements as a .npy with the given shape.
template<typename Scalar>
static void _ckpt_save_npy(const fs::path& path, const Scalar* d_ptr, size_t numel,
                           const npy::shape_t& shape) {
    if (d_ptr == nullptr || numel == 0) return;
    auto host = _ckpt_d2h<Scalar>(d_ptr, numel);
    npy::npy_data_ptr<Scalar> data{host.data(), shape, false};
    npy::write_npy(path.string(), data);
}

template<typename T>
static npy::shape_t _ckpt_shape_5d(const DeviceTensor5D<T>& t) {
    return {
        (unsigned long)t.template size<0>(),
        (unsigned long)t.template size<1>(),
        (unsigned long)t.template size<2>(),
        (unsigned long)t.template size<3>(),
        (unsigned long)t.template size<4>(),
    };
}

} // anon namespace


void engine_save_checkpoint(
    std::string output_dir,
    bool full_dump,
    int step
) {
    using namespace Buffers;

    fs::path out_root(output_dir);
    _ckpt_mkdir(out_root);

    cudaDeviceSynchronize();

    const int64_t N = cur_num_splats;
    const int K = num_sh;

    // --- D->H copy world splat tensors at cur_num_splats ---
    auto h_means        = _ckpt_d2h<float3>(world_means.data_ptr(),        (size_t)N);
    auto h_quats        = _ckpt_d2h<float4>(world_quats.data_ptr(),        (size_t)N);
    auto h_scales       = _ckpt_d2h<float3>(world_scales.data_ptr(),       (size_t)N);
    auto h_opacities    = _ckpt_d2h<float >(world_opacities.data_ptr(),    (size_t)N);
    auto h_features_dc  = _ckpt_d2h<float3>(world_features_dc.data_ptr(),  (size_t)N);
    // features_sh is stored [max_N, K] contiguous; copying first N rows = N*K float3.
    std::vector<float3> h_features_sh;
    if (K > 0) {
        h_features_sh = _ckpt_d2h<float3>(world_features_sh.data_ptr(), (size_t)(N * K));
    }

    // --- Filter mask: NaN/Inf + low-opacity (logit < logit(1/255) ~= -5.5373) ---
    const float OPA_MIN = -5.5373f;
    auto fin1 = [](float v) { return std::isfinite(v); };
    auto fin3 = [&](const float3& v) { return fin1(v.x) && fin1(v.y) && fin1(v.z); };
    auto fin4 = [&](const float4& v) { return fin1(v.x) && fin1(v.y) && fin1(v.z) && fin1(v.w); };

    std::vector<char> keep((size_t)N, 1);
    int64_t kept = 0, nan_dropped = 0, lowopa_dropped = 0;
    for (int64_t i = 0; i < N; i++) {
        bool ok = fin3(h_means[i]) && fin4(h_quats[i]) && fin3(h_scales[i])
                  && fin1(h_opacities[i]) && fin3(h_features_dc[i]);
        if (ok && K > 0) {
            for (int j = 0; j < K; j++) {
                if (!fin3(h_features_sh[i*K + j])) { ok = false; break; }
            }
        }
        if (!ok) { keep[i] = 0; nan_dropped++; continue; }
        if (h_opacities[i] < OPA_MIN) { keep[i] = 0; lowopa_dropped++; continue; }
        kept++;
    }

    // --- PLY (binary little-endian) ---
    const int vertex_floats = 3 + 3 + 3 + 3 * K + 1 + 3 + 4;  // pos + normal + dc + sh + opa + scale + rot
    fs::path ply_path = out_root / "splat.ply";
    {
        std::ofstream ply(ply_path, std::ios::binary);
        if (!ply) throw std::runtime_error("Failed to open PLY for write: " + ply_path.string());
        ply << "ply\n";
        ply << "format binary_little_endian 1.0\n";
        ply << "element vertex " << kept << "\n";
        ply << "property float x\nproperty float y\nproperty float z\n";
        ply << "property float nx\nproperty float ny\nproperty float nz\n";
        ply << "property float f_dc_0\nproperty float f_dc_1\nproperty float f_dc_2\n";
        for (int j = 0; j < 3 * K; j++) ply << "property float f_rest_" << j << "\n";
        ply << "property float opacity\n";
        ply << "property float scale_0\nproperty float scale_1\nproperty float scale_2\n";
        ply << "property float rot_0\nproperty float rot_1\nproperty float rot_2\nproperty float rot_3\n";
        ply << "end_header\n";

        // Buffer multiple rows per write for throughput.
        constexpr int ROWS_PER_FLUSH = 1024;
        std::vector<float> buf((size_t)vertex_floats * ROWS_PER_FLUSH);
        int rows_in_buf = 0;
        auto flush = [&]() {
            if (rows_in_buf == 0) return;
            ply.write(reinterpret_cast<const char*>(buf.data()),
                      (std::streamsize)rows_in_buf * vertex_floats * sizeof(float));
            rows_in_buf = 0;
        };
        for (int64_t i = 0; i < N; i++) {
            if (!keep[i]) continue;
            float* row = buf.data() + (size_t)rows_in_buf * vertex_floats;
            int p = 0;
            row[p++] = h_means[i].x; row[p++] = h_means[i].y; row[p++] = h_means[i].z;
            row[p++] = 0.f; row[p++] = 0.f; row[p++] = 0.f;  // nx, ny, nz
            row[p++] = h_features_dc[i].x; row[p++] = h_features_dc[i].y; row[p++] = h_features_dc[i].z;
            // SH rest in (3, K) order: r0..r_{K-1}, g0..g_{K-1}, b0..b_{K-1}
            for (int j = 0; j < K; j++) row[p++] = h_features_sh[i*K + j].x;
            for (int j = 0; j < K; j++) row[p++] = h_features_sh[i*K + j].y;
            for (int j = 0; j < K; j++) row[p++] = h_features_sh[i*K + j].z;
            row[p++] = h_opacities[i];
            row[p++] = h_scales[i].x; row[p++] = h_scales[i].y; row[p++] = h_scales[i].z;
            row[p++] = h_quats[i].x; row[p++] = h_quats[i].y; row[p++] = h_quats[i].z; row[p++] = h_quats[i].w;
            rows_in_buf++;
            if (rows_in_buf == ROWS_PER_FLUSH) flush();
        }
        flush();
    }

    // --- Top-level bilagrid / PPISP grids (param tables only, no optimizer state) ---
    if (bilagrid_rgb_enabled) {
        _ckpt_save_npy<float>(out_root / "bilagrid_rgb.npy",
                              bilagrid_rgb_grids.data_ptr(),
                              (size_t)bilagrid_rgb_grids.numel(),
                              _ckpt_shape_5d(bilagrid_rgb_grids));
    }
    if (bilagrid_depth_enabled) {
        _ckpt_save_npy<float>(out_root / "bilagrid_depth.npy",
                              bilagrid_depth_grids.data_ptr(),
                              (size_t)bilagrid_depth_grids.numel(),
                              _ckpt_shape_5d(bilagrid_depth_grids));
        if (bilagrid_depth_scalars.data_ptr() && bilagrid_depth_scalars.size() > 0) {
            _ckpt_save_npy<float>(out_root / "bilagrid_depth_scalars.npy",
                                  bilagrid_depth_scalars.data_ptr(),
                                  (size_t)bilagrid_depth_scalars.size(),
                                  {(unsigned long)bilagrid_depth_scalars.size()});
        }
    }
    if (bilagrid_normal_enabled) {
        _ckpt_save_npy<float>(out_root / "bilagrid_normal.npy",
                              bilagrid_normal_grids.data_ptr(),
                              (size_t)bilagrid_normal_grids.numel(),
                              _ckpt_shape_5d(bilagrid_normal_grids));
    }
    if (ppisp_enabled) {
        _ckpt_save_npy<float>(out_root / "ppisp.npy",
                              ppisp_params.data_ptr(),
                              (size_t)(ppisp_params.size<0>() * ppisp_params.size<1>()),
                              {(unsigned long)ppisp_params.size<0>(),
                               (unsigned long)ppisp_params.size<1>()});
    }

    // --- meta.txt ---
    {
        std::ofstream meta((out_root / "meta.txt").string());
        meta << "step=" << step << "\n";
        meta << "primitive=" << primitive << "\n";
        meta << "cur_num_splats=" << cur_num_splats << "\n";
        meta << "max_num_splats=" << max_num_splats << "\n";
        meta << "num_sh=" << num_sh << "\n";
        meta << "sh_degree=" << sh_degree << "\n";
        meta << "num_cameras=" << num_cameras << "\n";
        meta << "image_width=" << image_width << "\n";
        meta << "image_height=" << image_height << "\n";
        meta << "camera_model=" << camera_model_str << "\n";
        meta << "train_quantize_sh=" << (train_quantize_sh ? 1 : 0) << "\n";
        // QuantizedAdamState codec signature — tells offline tools which
        // byte-to-(g1, g2) inverse to use when reading the *_packed.npy
        // streams. Currently every quantized buffer (SH + bilagrid) uses the
        // same scheme: 8-bit AoS, 256-cell blocks, joint (u, log_s) primitives
        // with u linear and log_s = log1p(sqrt_g2 / eps), eps = 1e-15.
        meta << "quant_codec=joint_u_log_sqrt_g2_v1\n";
        meta << "quant_bits=8\n";
        meta << "quant_block_size=256\n";
        meta << "quant_eps=1e-15\n";
        meta << "use_per_splat_bias_correction=" << (use_per_splat_bias_correction ? 1 : 0) << "\n";
        meta << "ply_kept=" << kept << "\n";
        meta << "ply_nan_dropped=" << nan_dropped << "\n";
        meta << "ply_low_opacity_dropped=" << lowopa_dropped << "\n";
        if (bilagrid_rgb_enabled) {
            meta << "bilagrid_rgb_enabled=1\n";
            meta << "bilagrid_rgb_type=" << bilagrid_rgb_type << "\n";
            meta << "bilagrid_rgb_C=" << bilagrid_rgb_C << "\n";
            meta << "bilagrid_rgb_L=" << bilagrid_rgb_grids.size<1>() << "\n";
            meta << "bilagrid_rgb_H=" << bilagrid_rgb_grids.size<2>() << "\n";
            meta << "bilagrid_rgb_W=" << bilagrid_rgb_grids.size<3>() << "\n";
            meta << "bilagrid_rgb_quantize_optim=" << (bilagrid_rgb_quantize_optim ? 1 : 0) << "\n";
        }
        if (bilagrid_depth_enabled) {
            meta << "bilagrid_depth_enabled=1\n";
            meta << "bilagrid_depth_L=" << bilagrid_depth_grids.size<1>() << "\n";
            meta << "bilagrid_depth_H=" << bilagrid_depth_grids.size<2>() << "\n";
            meta << "bilagrid_depth_W=" << bilagrid_depth_grids.size<3>() << "\n";
            meta << "bilagrid_depth_quantize_optim=" << (bilagrid_depth_quantize_optim ? 1 : 0) << "\n";
        }
        if (bilagrid_normal_enabled) {
            meta << "bilagrid_normal_enabled=1\n";
            meta << "bilagrid_normal_L=" << bilagrid_normal_grids.size<1>() << "\n";
            meta << "bilagrid_normal_H=" << bilagrid_normal_grids.size<2>() << "\n";
            meta << "bilagrid_normal_W=" << bilagrid_normal_grids.size<3>() << "\n";
            meta << "bilagrid_normal_quantize_optim=" << (bilagrid_normal_quantize_optim ? 1 : 0) << "\n";
        }
        if (ppisp_enabled) {
            meta << "ppisp_enabled=1\n";
            meta << "ppisp_param_type=" << ppisp_param_type << "\n";
            meta << "ppisp_num_params=" << ppisp_num_params << "\n";
        }
    }

    if (!full_dump) return;

    // --- full/ subfolder: every persistent buffer, one .npy per buffer ---
    fs::path full_root = out_root / "full";
    _ckpt_mkdir(full_root);

    // World splat parameters at max_num_splats (full table).
    auto save_dv_f3 = [&](const char* name, const DeviceVector<float3>& v) {
        if (!v.data_ptr() || v.size() == 0) return;
        _ckpt_save_npy<float>(full_root / name,
                              reinterpret_cast<const float*>(v.data_ptr()),
                              (size_t)v.size() * 3,
                              {(unsigned long)v.size(), 3ul});
    };
    auto save_dv_f4 = [&](const char* name, const DeviceVector<float4>& v) {
        if (!v.data_ptr() || v.size() == 0) return;
        _ckpt_save_npy<float>(full_root / name,
                              reinterpret_cast<const float*>(v.data_ptr()),
                              (size_t)v.size() * 4,
                              {(unsigned long)v.size(), 4ul});
    };
    auto save_dv_f1 = [&](const char* name, const DeviceVector<float>& v) {
        if (!v.data_ptr() || v.size() == 0) return;
        _ckpt_save_npy<float>(full_root / name, v.data_ptr(),
                              (size_t)v.size(), {(unsigned long)v.size()});
    };
    auto save_dt2d_f3 = [&](const char* name, const DeviceTensor2D<float3>& t) {
        if (!t.data_ptr()) return;
        _ckpt_save_npy<float>(full_root / name,
                              reinterpret_cast<const float*>(t.data_ptr()),
                              (size_t)t.size<0>() * t.size<1>() * 3,
                              {(unsigned long)t.size<0>(), (unsigned long)t.size<1>(), 3ul});
    };
    auto save_dt2d_uc3 = [&](const char* name, const DeviceTensor2D<uchar3>& t) {
        if (!t.data_ptr()) return;
        _ckpt_save_npy<uint8_t>(full_root / name,
                                reinterpret_cast<const uint8_t*>(t.data_ptr()),
                                (size_t)t.size<0>() * t.size<1>() * 3,
                                {(unsigned long)t.size<0>(), (unsigned long)t.size<1>(), 3ul});
    };
    auto save_dt2d_f = [&](const char* name, const DeviceTensor2D<float>& t) {
        if (!t.data_ptr()) return;
        _ckpt_save_npy<float>(full_root / name, t.data_ptr(),
                              (size_t)t.size<0>() * t.size<1>(),
                              {(unsigned long)t.size<0>(), (unsigned long)t.size<1>()});
    };
    auto save_5d = [&](const char* name, const DeviceTensor5D<float>& t) {
        if (!t.data_ptr()) return;
        _ckpt_save_npy<float>(full_root / name, t.data_ptr(),
                              (size_t)t.numel(), _ckpt_shape_5d(t));
    };
    auto save_dv_u8 = [&](const char* name, const DeviceVector<uint8_t>& v) {
        if (!v.data_ptr() || v.size() == 0) return;
        _ckpt_save_npy<uint8_t>(full_root / name, v.data_ptr(),
                                (size_t)v.size(), {(unsigned long)v.size()});
    };

    // World params (max_num_splats sized).
    save_dv_f3("world_means.npy",        world_means);
    save_dv_f4("world_quats.npy",        world_quats);
    save_dv_f3("world_scales.npy",       world_scales);
    save_dv_f1("world_opacities.npy",    world_opacities);
    save_dv_f3("world_features_dc.npy",  world_features_dc);
    save_dt2d_f3("world_features_sh.npy", world_features_sh);

    // World Adam optimizer state.
    save_dv_f3("g1_means.npy",        g1_means);
    save_dv_f3("g2_means.npy",        g2_means);
    save_dv_f4("g1_quats.npy",        g1_quats);
    save_dv_f4("g2_quats.npy",        g2_quats);
    save_dv_f3("g1_scales.npy",       g1_scales);
    save_dv_f3("g2_scales.npy",       g2_scales);
    save_dv_f1("g1_opacities.npy",    g1_opacities);
    save_dv_f1("g2_opacities.npy",    g2_opacities);
    save_dv_f3("g1_features_dc.npy",  g1_features_dc);
    save_dv_f3("g2_features_dc.npy",  g2_features_dc);
    if (train_quantize_sh) {
        // Joint (u, sqrt_g2) AoS: single packed byte stream + per-block bounds.
        if (sh_quant_state.packed.data_ptr())
            save_dv_u8("sh_quant_packed.npy", sh_quant_state.packed);
        if (sh_quant_state.bounds.data_ptr())
            save_dv_f4("sh_quant_bounds.npy", sh_quant_state.bounds);
    } else {
        save_dt2d_f3("g1_features_sh.npy", g1_features_sh);
        save_dt2d_f3("g2_features_sh.npy", g2_features_sh);
    }

    // Per-splat aux.
    if (use_per_splat_bias_correction
        && bias_correction_steps.data_ptr() && bias_correction_steps.size() > 0) {
        _ckpt_save_npy<int32_t>(full_root / "bias_correction_steps.npy",
                                bias_correction_steps.data_ptr(),
                                (size_t)bias_correction_steps.size(),
                                {(unsigned long)bias_correction_steps.size()});
    }
    if (accum_buffer.data_ptr() && accum_buffer.size() > 0) {
        _ckpt_save_npy<float>(full_root / "accum_buffer.npy",
                              reinterpret_cast<const float*>(accum_buffer.data_ptr()),
                              (size_t)accum_buffer.size() * 2,
                              {(unsigned long)accum_buffer.size(), 2ul});
    }

    // Bilagrid grids + optimizer state, per type.
    auto save_bilagrid = [&](const char* prefix, bool enabled, bool quantize,
                             const DeviceTensor5D<float>& grids,
                             const DeviceTensor5D<float>& g1_f,
                             const DeviceTensor5D<float>& g2_f,
                             const QuantizedAdamState<8, 256>& quant_state) {
        if (!enabled) return;
        save_5d((std::string(prefix) + ".npy").c_str(), grids);
        if (quantize) {
            // AoS packed (u, sqrt_g2) byte stream + per-block bounds.
            save_dv_u8((std::string(prefix) + "_packed.npy").c_str(), quant_state.packed);
            save_dv_f4((std::string(prefix) + "_quant_bounds.npy").c_str(), quant_state.bounds);
        } else {
            save_5d((std::string(prefix) + "_g1.npy").c_str(), g1_f);
            save_5d((std::string(prefix) + "_g2.npy").c_str(), g2_f);
        }
    };
    save_bilagrid("bilagrid_rgb",   bilagrid_rgb_enabled,   bilagrid_rgb_quantize_optim,
                  bilagrid_rgb_grids, bilagrid_rgb_g1, bilagrid_rgb_g2,
                  bilagrid_rgb_quant_state);
    save_bilagrid("bilagrid_depth", bilagrid_depth_enabled, bilagrid_depth_quantize_optim,
                  bilagrid_depth_grids, bilagrid_depth_g1, bilagrid_depth_g2,
                  bilagrid_depth_quant_state);
    if (bilagrid_depth_enabled) {
        save_dv_f1("bilagrid_depth_scalars.npy", bilagrid_depth_scalars);
    }
    save_bilagrid("bilagrid_normal", bilagrid_normal_enabled, bilagrid_normal_quantize_optim,
                  bilagrid_normal_grids, bilagrid_normal_g1, bilagrid_normal_g2,
                  bilagrid_normal_quant_state);

    // PPISP table + Adam moments.
    if (ppisp_enabled) {
        save_dt2d_f("ppisp.npy",    ppisp_params);
        save_dt2d_f("ppisp_g1.npy", ppisp_g1);
        save_dt2d_f("ppisp_g2.npy", ppisp_g2);
    }
}
