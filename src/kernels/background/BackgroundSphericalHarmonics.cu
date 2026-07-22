// Background spherical-harmonics renderer - produces a per-pixel skybox color
// by evaluating an SH function along the world-space ray direction for each
// camera pixel. SH math lives in harmonics.slang (sh{N}_to_color_dir{,_vjp_*});
// ray generation + camera-to-world rotation come from projection_utils.slang
// (generate_ray + transform_ray_d), the same path used by the eval3d
// rasterizer. This keeps the convention consistent: viewmats are world->camera
// row-major (the engine's stored form), the Y/Z axis flip is already baked
// into them upstream, and distortion is handled per the camera model.
//
// Inputs and outputs are caller-allocated TorchTensorViews; no ATen
// allocations and no transient R_wc scratch buffer.

#include "kernels/background/BackgroundSphericalHarmonics.cuh"

#include "core/Common.cuh"

#include "generated/slang.cuh"
namespace SlangProjectionUtils {
#include "generated/set_namespace.cuh"
#include "generated/projection_utils.cuh"
}
namespace SlangHarmonics {
#include "generated/set_namespace.cuh"
#include "generated/harmonics.cuh"
}

#include <stdexcept>


// ----- Common ray-from-viewmat path used by fwd + bwd ------------------------

// Load viewmat[:3,:3] (row-major) as a slang-style float3x3, and intrinsics
// and (optional) distortion for one camera.
struct BgCamera {
    float3x3 R;            // world->camera (slang convention)
    float fx, fy, cx, cy;
    CameraDistortionCoeffs dist;
};

static __device__ __forceinline__ BgCamera _load_bg_camera(
    const float*  viewmats,    // [N_cam, 16] row-major
    const float4* intrins,     // [N_cam]
    const float*  dist_coeffs, // [N_cam, 10] or null
    int cam
) {
    const float* vm = viewmats + (int64_t)cam * 16;
    BgCamera c;
    c.R = float3x3{
        vm[0], vm[1], vm[2],   // row 0
        vm[4], vm[5], vm[6],   // row 1
        vm[8], vm[9], vm[10],  // row 2
    };
    float4 intr = intrins[cam];
    c.fx = intr.x; c.fy = intr.y; c.cx = intr.z; c.cy = intr.w;
    if (dist_coeffs == nullptr) {
        #pragma unroll
        for (int i = 0; i < 10; ++i) c.dist.m_data[i] = 0.0f;
    } else {
        const float* d = dist_coeffs + (int64_t)cam * 10;
        #pragma unroll
        for (int i = 0; i < 10; ++i) c.dist.m_data[i] = d[i];
    }
    return c;
}

// Compute world-space ray direction at pixel (px, py) for the given camera.
// Returns false if the pixel maps to an invalid ray (distortion model rejected).
static __device__ __forceinline__ bool _bg_pixel_world_ray(
    const BgCamera& c, int camera_model, int px, int py, float3& out_world_dir
) {
    const float2 uv = {
        ((float)px + 0.5f - c.cx) / c.fx,
        ((float)py + 0.5f - c.cy) / c.fy,
    };
    float3 raydir;
    if (!SlangProjectionUtils::generate_ray(uv, camera_model, c.dist, &raydir))
        return false;
    out_world_dir = SlangProjectionUtils::transform_ray_d(c.R, raydir);
    return true;
}


// ---- Forward kernel ---------------------------------------------------------

template<int SH_DEGREE>
__global__ void render_background_sh_forward_kernel(
    const dim3 img_size,
    int camera_model,
    const float*  __restrict__ viewmats,    // [B, 16] (per-batch image)
    const float4* __restrict__ intrins,     // [B]
    const float*  __restrict__ dist_coeffs, // [B, 10] or null
    const float3* __restrict__ sh_coeffs,   // [(SH_DEGREE+1)^2]
    float3*       __restrict__ out_img      // [B, H, W]
) {
    unsigned bi = blockIdx.z * blockDim.z + threadIdx.z;
    unsigned py = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned px = blockIdx.x * blockDim.x + threadIdx.x;
    if (py >= img_size.y || px >= img_size.x || bi >= img_size.z) return;
    int32_t pix_id = (bi * img_size.y + py) * img_size.x + px;

    // viewmats / intrins / dist_coeffs are sized to the current batch (the
    // caller has already gathered per-image data into a [B, ...] layout), so
    // we index directly by bi — NOT by some global camera id. Indexing by
    // `cam_indices[bi]` (the global-cam-id pattern used by bilagrid/PPISP) is
    // wrong here: that would read past the buffer's [B] worth of entries.
    BgCamera c = _load_bg_camera(viewmats, intrins, dist_coeffs, (int)bi);

    float3 world_dir;
    if (!_bg_pixel_world_ray(c, camera_model, (int)px, (int)py, world_dir)) {
        out_img[pix_id] = {0.0f, 0.0f, 0.0f};
        return;
    }

    const float3 dc = sh_coeffs[0];
    float3* coeffs_l1 = (float3*)(sh_coeffs + 1);

    float3 color;
    if constexpr (SH_DEGREE == 0)
        color = SlangHarmonics::sh0_to_color_dir(world_dir, dc, coeffs_l1);
    else if constexpr (SH_DEGREE == 1)
        color = SlangHarmonics::sh1_to_color_dir(world_dir, dc, coeffs_l1);
    else if constexpr (SH_DEGREE == 2)
        color = SlangHarmonics::sh2_to_color_dir(world_dir, dc, coeffs_l1);
    else if constexpr (SH_DEGREE == 3)
        color = SlangHarmonics::sh3_to_color_dir(world_dir, dc, coeffs_l1);
    else if constexpr (SH_DEGREE == 4)
        color = SlangHarmonics::sh4_to_color_dir(world_dir, dc, coeffs_l1);

    color.x = fmaxf(color.x + 0.5f, 0.0f);
    color.y = fmaxf(color.y + 0.5f, 0.0f);
    color.z = fmaxf(color.z + 0.5f, 0.0f);
    out_img[pix_id] = color;
}


// ---- Backward kernel --------------------------------------------------------

template<int SH_DEGREE>
__global__ void __launch_bounds__(512) render_background_sh_backward_kernel(
    const dim3 img_size,
    int camera_model,
    const float*  __restrict__ viewmats,    // [B, 16]
    const float4* __restrict__ intrins,     // [B]
    const float*  __restrict__ dist_coeffs, // [B, 10] or null
    const float3* __restrict__ sh_coeffs,
    const float3* __restrict__ out_color,
    const float3* __restrict__ v_out_color,
    float3*       __restrict__ v_sh_coeffs
) {
    // Flattened (pix, cam) layout: improves coalescing for the atomicAdd-heavy
    // backward (one block reduces into a small SH coefficient table).
    unsigned bi = blockIdx.y * blockDim.y + threadIdx.y;
    int32_t  pix_id = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned px = pix_id % img_size.x, py = pix_id / img_size.x;
    pix_id += bi * img_size.y * img_size.x;
    bool inside = (py < img_size.y && px < img_size.x && bi < img_size.z);

    // dL/d(rgb_post_clamp) -> dL/d(rgb_pre_clamp). Zero where the clamp kicked in.
    float3 v_color = {0.0f, 0.0f, 0.0f};
    if (inside) {
        float3 color = out_color[pix_id];
        v_color = v_out_color[pix_id];
        if (color.x == 0.0f || !isfinite(v_color.x)) v_color.x = 0.0f;
        if (color.y == 0.0f || !isfinite(v_color.y)) v_color.y = 0.0f;
        if (color.z == 0.0f || !isfinite(v_color.z)) v_color.z = 0.0f;
    }

    // See forward kernel: index directly by bi (per-batch buffers).
    BgCamera c = _load_bg_camera(viewmats, intrins, dist_coeffs, inside ? (int)bi : 0);

    float3 world_dir = {0.0f, 0.0f, 1.0f};
    if (inside) {
        if (!_bg_pixel_world_ray(c, camera_model, (int)px, (int)py, world_dir)) {
            inside = false;
            v_color = {0.0f, 0.0f, 0.0f};
        }
    }

    const float3 dc = sh_coeffs[0];
    float3* coeffs_l1   = (float3*)(sh_coeffs   + 1);
    float3* v_coeffs_l1 = (float3*)(v_sh_coeffs + 1);

    float3 v_dc  = {0.0f, 0.0f, 0.0f};
    float3 v_dir = {0.0f, 0.0f, 0.0f};

    // Block-reduced atomic mode: 512 threads cooperate to issue ONE atomicAdd
    // per channel per block (instead of one per pixel) for every higher-order
    // SH band. Requires BLOCK_SIZE == 512 and that every thread in the block
    // reach EVERY internal helper call (true here: the kernel has no early
    // returns; out-of-bounds threads just pass v_color=0 through).
    //
    // The trailing 0xFFFFFFFFu is the warp active mask: slang's CUDA codegen
    // surfaces it for any wave-intrinsic-using function, and it's safe to use
    // a full mask here since no warp lane is masked off above.
    const uint tid_in_block = threadIdx.x
        + threadIdx.y * blockDim.x
        + threadIdx.z * blockDim.x * blockDim.y;
    if constexpr (SH_DEGREE == 0)
        SlangHarmonics::sh0_to_color_dir_vjp_block_atomic(
            world_dir, dc, coeffs_l1, v_color, &v_dc, v_coeffs_l1, &v_dir, tid_in_block);
    else if constexpr (SH_DEGREE == 1)
        SlangHarmonics::sh1_to_color_dir_vjp_block_atomic(
            world_dir, dc, coeffs_l1, v_color, &v_dc, v_coeffs_l1, &v_dir, tid_in_block, 0xFFFFFFFFu);
    else if constexpr (SH_DEGREE == 2)
        SlangHarmonics::sh2_to_color_dir_vjp_block_atomic(
            world_dir, dc, coeffs_l1, v_color, &v_dc, v_coeffs_l1, &v_dir, tid_in_block, 0xFFFFFFFFu);
    else if constexpr (SH_DEGREE == 3)
        SlangHarmonics::sh3_to_color_dir_vjp_block_atomic(
            world_dir, dc, coeffs_l1, v_color, &v_dc, v_coeffs_l1, &v_dir, tid_in_block, 0xFFFFFFFFu);
    else if constexpr (SH_DEGREE == 4)
        SlangHarmonics::sh4_to_color_dir_vjp_block_atomic(
            world_dir, dc, coeffs_l1, v_color, &v_dc, v_coeffs_l1, &v_dir, tid_in_block, 0xFFFFFFFFu);

    // DC gradient: block-reduce one float3 atomic to the global table.
    atomicAddFVec<512>(&v_sh_coeffs[0], v_dc);

    // Camera gradient is not currently propagated (engine path holds cameras
    // fixed; v_dir is discarded). A v_viewmats output can be added later via
    // SlangProjectionUtils::transform_ray_d_vjp.
}


// ---- Host-side dispatchers --------------------------------------------------

static inline int64_t _batch_count(const TorchTensorView& t, int64_t per_item) {
    int64_t n = 1;
    for (auto s : std::get<2>(t)) n *= s;
    return n / per_item;
}

static inline int _camera_model_int(const std::string& camera_model) {
    auto cm = cmt(camera_model);
    if (cm == (CameraModelType)-1)
        throw std::runtime_error("Camera model " + camera_model + " is not supported for skybox");
    return (int)cm;
}


/*[AutoHeaderGeneratorExport]*/
void render_background_sh_forward(
    int w,
    int h,
    std::string camera_model,
    int sh_degree,                       // actual SH degree (0..4)
    TorchTensorView viewmats,            // [B, 4, 4] row-major world->camera (per-batch)
    TorchTensorView intrins,             // [B, 4]
    TorchTensorView dist_coeffs,         // [B, 10]; null/empty -> zeros
    TorchTensorView sh_coeffs,           // [(sh_degree+1)^2, 3]
    TorchTensorView out_color            // [B, H, W, 3]  pre-allocated
) {
    if (sh_degree < 0 || sh_degree > 4)
        throw std::runtime_error("render_background_sh_forward: sh_degree must be in [0, 4]");

    // Batch count of output = product of out_color leading dims / (H*W*3).
    int64_t b = _batch_count(out_color, (int64_t)h * w * 3);
    if (b * h * w == 0) return;
    if (_batch_count(viewmats, 16) != b || _batch_count(intrins, 4) != b)
        throw std::runtime_error("viewmats/intrins must have batch dim matching out_color");

    const dim3 img_size = {(uint32_t)w, (uint32_t)h, (uint32_t)b};
    const float*  p_vm        = (const float*)std::get<0>(viewmats);
    const float4* p_intrins   = (const float4*)std::get<0>(intrins);
    const float*  p_dist      = (const float*)std::get<0>(dist_coeffs);  // null -> zero default
    const float3* p_sh_coeffs = (const float3*)std::get<0>(sh_coeffs);
    float3*       p_out       = (float3*)std::get<0>(out_color);

    int cm = _camera_model_int(camera_model);

    #define LAUNCH(DEG) \
        render_background_sh_forward_kernel<DEG> \
            <<<_LAUNCH_ARGS_3D(w, h, b, 16, 16, 1)>>>( \
                img_size, cm, p_vm, p_intrins, p_dist, p_sh_coeffs, p_out)
    switch (sh_degree) {
        case 0: LAUNCH(0); break;
        case 1: LAUNCH(1); break;
        case 2: LAUNCH(2); break;
        case 3: LAUNCH(3); break;
        case 4: LAUNCH(4); break;
    }
    #undef LAUNCH
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


/*[AutoHeaderGeneratorExport]*/
void render_background_sh_backward(
    int w,
    int h,
    std::string camera_model,
    int sh_degree,
    TorchTensorView viewmats,
    TorchTensorView intrins,
    TorchTensorView dist_coeffs,
    TorchTensorView sh_coeffs,
    TorchTensorView out_color,
    TorchTensorView v_out_color,
    TorchTensorView v_sh_coeffs           // [(sh_degree+1)^2, 3] zeroed
) {
    if (sh_degree < 0 || sh_degree > 4)
        throw std::runtime_error("render_background_sh_backward: sh_degree must be in [0, 4]");

    const auto& oc = std::get<2>(out_color);
    if (oc.size() < 3 || oc[oc.size() - 1] != 3 ||
        oc[oc.size() - 2] != (int64_t)w || oc[oc.size() - 3] != (int64_t)h)
        throw std::runtime_error("out_color must be (..., H, W, 3)");
    const auto& vc = std::get<2>(v_out_color);
    if (vc.size() < 3 || vc[vc.size() - 1] != 3 ||
        vc[vc.size() - 2] != (int64_t)w || vc[vc.size() - 3] != (int64_t)h)
        throw std::runtime_error("v_out_color must be (..., H, W, 3)");

    int64_t b = _batch_count(out_color, (int64_t)h * w * 3);
    if (b * h * w == 0) return;

    const dim3 img_size = {(uint32_t)w, (uint32_t)h, (uint32_t)b};
    const float*  p_vm          = (const float*)std::get<0>(viewmats);
    const float4* p_intrins     = (const float4*)std::get<0>(intrins);
    const float*  p_dist        = (const float*)std::get<0>(dist_coeffs);
    const float3* p_sh_coeffs   = (const float3*)std::get<0>(sh_coeffs);
    const float3* p_out_color   = (const float3*)std::get<0>(out_color);
    const float3* p_v_out_color = (const float3*)std::get<0>(v_out_color);
    float3*       p_v_sh        = (float3*)std::get<0>(v_sh_coeffs);

    int cm = _camera_model_int(camera_model);

    #define LAUNCH(DEG) \
        render_background_sh_backward_kernel<DEG> \
            <<<_LAUNCH_ARGS_2D((uint32_t)(w * h), b, 512, 1)>>>( \
                img_size, cm, p_vm, p_intrins, p_dist, p_sh_coeffs, \
                p_out_color, p_v_out_color, p_v_sh)
    switch (sh_degree) {
        case 0: LAUNCH(0); break;
        case 1: LAUNCH(1); break;
        case 2: LAUNCH(2); break;
        case 3: LAUNCH(3); break;
        case 4: LAUNCH(4); break;
    }
    #undef LAUNCH
    CHECK_DEVICE_ERROR(cudaGetLastError());
}
