#include "FusedProjectionBwdOptim.cuh"

#include <Common.cuh>

#include <cooperative_groups.h>
namespace cg = cooperative_groups;

#include <cub/cub.cuh>


template<
    typename SplatPrimitive,
    CameraModelType camera_model,
    HessianDiagonalOutputMode hessian_diagonal_output_mode,
    const bool use_scale_agnostic_mean,
    const bool use_color_trust_region,
    const bool color_is_linear
>
void fused_projection_bwd_optimizer_3dgs_kernel_wrapper(
    cudaStream_t stream,
    // fwd inputs
    const uint32_t C,
    const uint32_t N,
    const uint32_t num_sh_buffer,
    typename SplatPrimitive::WorldBuffer splats_world,
    const float *__restrict__ viewmats, // [C, 4, 4]
    const float4 *__restrict__ intrins,  // [C, 4], fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const uint32_t image_width,
    const uint32_t image_height,
    // fwd outputs
    const int32_t *__restrict__ camera_id_bounds,   // [N+1]
    const int32_t *__restrict__ camera_ids,   // [nnz] -- ORIGINAL (unsorted) order
    const int32_t *__restrict__ perm,         // [nnz] -- sorted_pos -> original_pos
    const float4 *__restrict__ aabb,   // [C, N, 4] or [nnz, 4]
    // grad outputs from rasterization
    typename SplatPrimitive::WorldBuffer v_splats_world,
    typename SplatPrimitive::WorldBuffer vr_splats_world,
    typename SplatPrimitive::WorldBuffer h_splats_world,
    typename SplatPrimitive::ScreenBuffer v_splats_screen,
    typename SplatPrimitive::ScreenBuffer vr_splats_screen,
    typename SplatPrimitive::ScreenBuffer h_splats_screen,
    // optimizer states
    typename SplatPrimitive::WorldBuffer g1_splats_world,
    typename SplatPrimitive::WorldBuffer g2_splats_world,
    const uint8_t* __restrict__ sh_packed,      // AoS (u, sqrt_g2) packed SH state
    float4* __restrict__ sh_quant_bounds,
    const uint8_t* __restrict__ sh_value_packed,
    float2* __restrict__ sh_value_bounds,
    // float *__restrict__ v_viewmats // [C, 4, 4] optional
    // optimizer params
    const float* __restrict__ radii,
    const float lr_means,
    const float lr_quats,
    const float lr_scales,
    const float lr_opacs,
    const float lr_features_dc,
    const float lr_features_sh,
    const float max_gauss_ratio,
    const float scale_regularization_weight,
    const float mcmc_opacity_reg_weight,
    const float mcmc_scale_reg_weight,
    const float erank_reg_weight,
    const float erank_reg_weight_s3,
    const float quat_norm_reg_weight,
    const float sh_reg_weight,
    const float eps_tr,
    const int32_t scalar_step,
    const int32_t* __restrict__ steps,
    const int sh_quant_bits,
    const int sh_value_bits
);


// Fill buffer[i] = i for i in [0, n). Used to build the identity permutation
// that we sort alongside gaussian_ids so the FPBO kernel can map sorted
// positions back to the original nnz position of aabb / v_splats_screen.
__global__ static void fpbo_iota_kernel(int64_t n, int32_t* __restrict__ buf) {
    int64_t i = (int64_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) buf[i] = (int32_t)i;
}


__global__ void camera_id_bounds_kernel(
    int64_t nnz,
    int64_t N,
    const int32_t* __restrict__ gaussian_ids,  // [nnz]
    int32_t* __restrict__ camera_id_bounds  // [N+1]
) {
    int64_t gid = blockIdx.x * blockDim.x + threadIdx.x;
    if (gid > nnz)
        return;

    int32_t cid = (gid == nnz) ? N : gaussian_ids[gid];
    int32_t cid_0 = (gid == 0) ? 0 : gaussian_ids[gid-1]+1;
    for (int32_t i = cid_0; i <= cid; ++i) {
        camera_id_bounds[i] = gid;
    }
}


template<
    typename SplatPrimitive,
    HessianDiagonalOutputMode hessian_diagonal_output_mode,
    bool use_scale_agnostic_mean,
    bool use_color_trust_region,
    bool color_is_linear
>
inline void launch_fused_projection_bwd_optimizer_3dgs_kernel(
    // fwd inputs
    const int64_t N,
    const uint32_t num_sh_buffer,
    std::vector<DeviceTensorFloatND> splats_world,
    TorchTensorView viewmats,  // [..., C, 4, 4]
    TorchTensorView intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const CameraModelType camera_model,
    const TorchTensorView dist_coeffs,
    // fwd outputs
    DeviceVector<int32_t> camera_ids,
    DeviceVector<int32_t> gaussian_ids,
    DeviceTensorFloatND aabb,
    // grad outputs
    const std::vector<DeviceTensorFloatND> v_splats_world,
    const std::optional<std::vector<DeviceTensorFloatND>> vr_splats_world,
    const std::optional<std::vector<DeviceTensorFloatND>> h_splats_world,
    const std::vector<DeviceTensorFloatND> v_splats_screen,
    const std::optional<std::vector<DeviceTensorFloatND>> vr_splats_screen,
    const std::optional<std::vector<DeviceTensorFloatND>> h_splats_screen,
    // optimizer states
    const std::vector<DeviceTensorFloatND> g1_splats_world,
    const std::vector<DeviceTensorFloatND> g2_splats_world,
    const std::optional<TorchTensorView> sh_packed,         // AoS packed SH state
    const std::optional<TorchTensorView> sh_quant_bounds,
    const std::optional<TorchTensorView> sh_value_packed,
    const std::optional<TorchTensorView> sh_value_bounds,
    // optimizer params
    DeviceVector<float> radii,
    const float lr_means,
    const float lr_quats,
    const float lr_scales,
    const float lr_opacs,
    const float lr_features_dc,
    const float lr_features_sh,
    const float max_gauss_ratio,
    const float scale_regularization_weight,
    const float mcmc_opacity_reg_weight,
    const float mcmc_scale_reg_weight,
    const float erank_reg_weight,
    const float erank_reg_weight_s3,
    const float quat_norm_reg_weight,
    const float sh_reg_weight,
    const float eps_tr,
    const int32_t scalar_step,
    const std::optional<TorchTensorView> steps,
    const int sh_quant_bits,
    const int sh_value_bits
) {
    uint32_t C = (uint32_t)std::get<2>(viewmats)[0]; // number of cameras (first dim)
    // Note: viewmats shape is [C, 4, 4], so index 0 = C

    if (N == 0)
        return;

    bool packed = camera_ids.data_ptr() && gaussian_ids.data_ptr();

    DeviceVector<int32_t> camera_id_bounds;
    DeviceVector<int32_t> gaussian_ids_sorted_buf;
    DeviceVector<int32_t> perm_buf;
    DeviceVector<int32_t> perm_sorted_buf;
    const int32_t* sorted_gaussian_ids_ptr = nullptr;
    const int32_t* sorted_perm_ptr = nullptr;

    if (packed) {
        long nnz = camera_ids.size();
        // Sort by gaussian_id with an identity permutation as the *value*
        // (NOT camera_ids). The original camera_ids / aabb / v_splats_screen
        // arrays produced by the forward + raster_bwd are indexed by the
        // out_idx assigned during ProjectionPackedFwd (the intersection-
        // mask prefix scan), which is NOT the same as the sorted-by-gid
        // position. Sorting camera_ids alongside gaussian_ids would keep
        // (sorted_gid, sorted_cid) paired but `aabb` and `v_splats_screen`
        // would still live in original order, so indexing them with the
        // sorted position cid_t reads the wrong intersection. Sorting a
        // permutation instead lets the kernel recover the original
        // out_idx = perm[cid_t] for every aabb / screen-grad load.
        gaussian_ids_sorted_buf.resize("fused_proj_bwd.gauss_sorted", nnz);
        perm_buf.resize("fused_proj_bwd.perm", nnz);
        perm_sorted_buf.resize("fused_proj_bwd.perm_sorted", nnz);

        fpbo_iota_kernel<<<_LAUNCH_ARGS_1D(nnz, 256)>>>(
            nnz, perm_buf.data_ptr()
        );
        CHECK_DEVICE_ERROR(cudaGetLastError());

        cub::DoubleBuffer<int32_t> d_keys(
            gaussian_ids.data_ptr(), gaussian_ids_sorted_buf.data_ptr()
        );
        cub::DoubleBuffer<int32_t> d_values(
            perm_buf.data_ptr(), perm_sorted_buf.data_ptr()
        );
        int n_bits = 0;
        while ((1U << n_bits) <= N)
            ++n_bits;
        CUB_WRAPPER(cub::DeviceRadixSort::SortPairs, d_keys, d_values, nnz, 0, n_bits);
        CHECK_DEVICE_ERROR(cudaGetLastError());

        sorted_gaussian_ids_ptr = d_keys.selector ? gaussian_ids_sorted_buf.data_ptr() : gaussian_ids.data_ptr();
        sorted_perm_ptr = d_values.selector ? perm_sorted_buf.data_ptr() : perm_buf.data_ptr();

        camera_id_bounds.resize("fused_proj_bwd.cam_bounds", (int64_t)(N+1));
        camera_id_bounds_kernel<<<_LAUNCH_ARGS_1D(nnz+1, 256)>>>(
            nnz, N,
            sorted_gaussian_ids_ptr,
            camera_id_bounds.data_ptr()
        );
        CHECK_DEVICE_ERROR(cudaGetLastError());
    }

    const float* viewmats_ptr = (const float*)std::get<0>(viewmats);
    const float4* intrins_ptr = (const float4*)std::get<0>(intrins);
    const int32_t* steps_ptr = steps.has_value() ? (const int32_t*)std::get<0>(steps.value()) : nullptr;

    #define _LAUNCH_ARGS ( \
            (cudaStream_t)0, C, N, num_sh_buffer, \
            splats_world, viewmats_ptr, intrins_ptr, dist_coeffs, \
            image_width, image_height, \
            packed ? camera_id_bounds.data_ptr() : nullptr, \
            packed ? camera_ids.data_ptr() : nullptr, \
            packed ? sorted_perm_ptr : nullptr, \
            (float4*)aabb.data_ptr(), \
            v_splats_world, \
            hessian_diagonal_output_mode != HessianDiagonalOutputMode::None ? vr_splats_world.value() : typename SplatPrimitive::WorldBuffer{}, \
            hessian_diagonal_output_mode != HessianDiagonalOutputMode::None ? h_splats_world.value() : typename SplatPrimitive::WorldBuffer{}, \
            v_splats_screen, \
            hessian_diagonal_output_mode != HessianDiagonalOutputMode::None ? vr_splats_screen.value() : typename SplatPrimitive::ScreenBuffer{}, \
            hessian_diagonal_output_mode != HessianDiagonalOutputMode::None ? h_splats_screen.value() : typename SplatPrimitive::ScreenBuffer{}, \
            g1_splats_world, g2_splats_world, \
            sh_packed.has_value() ? (const uint8_t*)std::get<0>(sh_packed.value()) : nullptr, \
            sh_quant_bounds.has_value() ? (float4*)std::get<0>(sh_quant_bounds.value()) : nullptr, \
            sh_value_packed.has_value() ? (const uint8_t*)std::get<0>(sh_value_packed.value()) : nullptr, \
            sh_value_bounds.has_value() ? (float2*)std::get<0>(sh_value_bounds.value()) : nullptr, \
            /*v_viewmats.has_value() ? v_viewmats.value().data_ptr<float>() : nullptr */ \
            radii.data_ptr(), \
            lr_means, lr_quats, lr_scales, lr_opacs, lr_features_dc, lr_features_sh, \
            max_gauss_ratio, scale_regularization_weight, \
            mcmc_opacity_reg_weight, mcmc_scale_reg_weight, erank_reg_weight, erank_reg_weight_s3, quat_norm_reg_weight, sh_reg_weight, \
            eps_tr, \
            scalar_step, steps_ptr, \
            sh_quant_bits, sh_value_bits \
        )

    if (camera_model == CameraModelType::PINHOLE)
        fused_projection_bwd_optimizer_3dgs_kernel_wrapper<SplatPrimitive, CameraModelType::PINHOLE, hessian_diagonal_output_mode, use_scale_agnostic_mean, use_color_trust_region, color_is_linear> _LAUNCH_ARGS;
    else if (camera_model == CameraModelType::FISHEYE)
        fused_projection_bwd_optimizer_3dgs_kernel_wrapper<SplatPrimitive, CameraModelType::FISHEYE, hessian_diagonal_output_mode, use_scale_agnostic_mean, use_color_trust_region, color_is_linear> _LAUNCH_ARGS;
    else if (camera_model == CameraModelType::EQUISOLID)
        fused_projection_bwd_optimizer_3dgs_kernel_wrapper<SplatPrimitive, CameraModelType::EQUISOLID, hessian_diagonal_output_mode, use_scale_agnostic_mean, use_color_trust_region, color_is_linear> _LAUNCH_ARGS;
    else
        throw std::runtime_error("Unsupported camera model");
    CHECK_DEVICE_ERROR(cudaGetLastError());

    #undef _LAUNCH_ARGS

}



// ================
// Per-primitive host wrappers
// ================
//
// All three primitives (3DGS, MipSplatting, 3DGUT) share the same kernel and
// host-launcher body — they only differ in which `PrimT<sh_degree>` template
// alias is plugged into the templated `launch_*` call. Factor the dispatch
// into one helper templated on the primitive class template, then expose a
// thin wrapper per primitive.

template<template<int> class PrimT>
static inline void _fused_projection_bwd_optimizer_dispatch(
    // fwd inputs
    const int64_t num_splats,
    const int max_sh_degree,
    std::vector<DeviceTensorFloatND> splats_world,
    TorchTensorView viewmats,
    TorchTensorView intrins,
    const uint32_t image_width,
    const uint32_t image_height,
    const std::string camera_model,
    const TorchTensorView dist_coeffs,
    // fwd outputs
    const DeviceVector<int32_t> camera_ids,
    const DeviceVector<int32_t> gaussian_ids,
    DeviceTensorFloatND aabb,
    // grad outputs
    const std::vector<DeviceTensorFloatND> v_splats_world,
    const std::optional<std::vector<DeviceTensorFloatND>> vr_splats_world,
    const std::optional<std::vector<DeviceTensorFloatND>> h_splats_world,
    const std::vector<DeviceTensorFloatND> v_splats_screen,
    const std::optional<std::vector<DeviceTensorFloatND>> vr_splats_screen,
    const std::optional<std::vector<DeviceTensorFloatND>> h_splats_screen,
    // optimizer states
    const std::vector<DeviceTensorFloatND> g1_splats_world,
    const std::vector<DeviceTensorFloatND> g2_splats_world,
    const std::optional<TorchTensorView> sh_packed,
    const std::optional<TorchTensorView> sh_quant_bounds,
    const std::optional<TorchTensorView> sh_value_packed,
    const std::optional<TorchTensorView> sh_value_bounds,
    // optimizer params
    DeviceVector<float> radii,
    const float lr_means,
    const float lr_quats,
    const float lr_scales,
    const float lr_opacs,
    const float lr_features_dc,
    const float lr_features_sh,
    const float max_gauss_ratio,
    const float scale_regularization_weight,
    const float mcmc_opacity_reg_weight,
    const float mcmc_scale_reg_weight,
    const float erank_reg_weight,
    const float erank_reg_weight_s3,
    const float quat_norm_reg_weight,
    const float sh_reg_weight,
    bool use_scale_agnostic_mean,
    bool use_color_trust_region,
    bool color_is_linear,
    float eps_tr,
    std::variant<int32_t, TorchTensorView> step,
    int sh_quant_bits,
    int sh_value_bits
) {
    int32_t scalar_step = std::get_if<int32_t>(&step) ? std::get<int32_t>(step) : -1;
    std::optional<TorchTensorView> steps_view = std::get_if<TorchTensorView>(&step) ?
        (std::optional<TorchTensorView>)std::get<TorchTensorView>(step) : std::nullopt;

    // Buffer's actual SH coef count (3 * (buffer_sh_degree*(buffer_sh_degree+2))).
    // The launch dispatch below caps the kernel's *template* sh degree at the
    // runtime sh_degree_to_use; the packed SH-Adam buffer is sized for the
    // model's max degree (engine().num_sh) and must be indexed with that
    // stride regardless of warmup.
    const int _buf_sh_deg = typename PrimT<0>::WorldBuffer(splats_world).sh_degree();
    const uint32_t num_sh_buffer = (uint32_t)(_buf_sh_deg * (_buf_sh_deg + 2));
    #define _ARGS ( \
        num_splats, num_sh_buffer, \
        splats_world, \
        viewmats, \
        intrins, \
        image_width, \
        image_height, \
        cmt(camera_model), \
        dist_coeffs, \
        camera_ids, \
        gaussian_ids, \
        aabb, \
        v_splats_world, \
        vr_splats_world, \
        h_splats_world, \
        v_splats_screen, \
        vr_splats_screen, \
        h_splats_screen, \
        g1_splats_world, \
        g2_splats_world, \
        sh_packed, \
        sh_quant_bounds, \
        sh_value_packed, \
        sh_value_bounds, \
        radii, \
        lr_means, \
        lr_quats, \
        lr_scales, \
        lr_opacs, \
        lr_features_dc, \
        lr_features_sh, \
        max_gauss_ratio, \
        scale_regularization_weight, \
        mcmc_opacity_reg_weight, \
        mcmc_scale_reg_weight, \
        erank_reg_weight, \
        erank_reg_weight_s3, \
        quat_norm_reg_weight, \
        sh_reg_weight, \
        eps_tr, \
        scalar_step, \
        steps_view, \
        sh_quant_bits, \
        sh_value_bits \
    )
    // Match projection_*_backward: cap kernel-template SH degree at the
    // runtime sh_degree_to_use so forward/backward agree on which bands
    // contribute (otherwise SH backward computes a spurious gradient on
    // unused bands, which propagates through `v_viewdir` into `v_mean`).
    int sh_degree = typename PrimT<0>::WorldBuffer(splats_world).sh_degree();
    sh_degree = std::min(sh_degree, max_sh_degree);
    // Pack the 3 bools (use_scale_agnostic_mean, use_color_trust_region,
    // color_is_linear) into a single 3-bit key for the dispatch table.
    const int dispatch_key =
        (int)use_scale_agnostic_mean
        | ((int)use_color_trust_region  << 1)
        | ((int)color_is_linear         << 2);
    #define LAUNCH_BOOL3(n, sam, tr, cl) \
        if (sh_degree == (n) && dispatch_key == ((int)(sam) | ((int)(tr)<<1) | ((int)(cl)<<2))) \
            return (void)launch_fused_projection_bwd_optimizer_3dgs_kernel<PrimT<n>, HessianDiagonalOutputMode::None, sam, tr, cl> _ARGS;
    #define LAUNCH(n) \
        LAUNCH_BOOL3(n, false, false, false) \
        LAUNCH_BOOL3(n, true,  false, false) \
        LAUNCH_BOOL3(n, false, true,  false) \
        LAUNCH_BOOL3(n, true,  true,  false) \
        LAUNCH_BOOL3(n, false, true,  true) \
        LAUNCH_BOOL3(n, true,  true,  true) \
        LAUNCH_BOOL3(n, false, false, true) \
        LAUNCH_BOOL3(n, true,  false, true)
    LAUNCH(3) LAUNCH(2) LAUNCH(1) LAUNCH(0) LAUNCH(4)
    #undef LAUNCH
    #undef LAUNCH_BOOL3
    #undef _ARGS
}


// Per-primitive thin wrappers. Identical signatures — the auto-header
// generator extracts each declaration via regex (no preprocessor), so the
// argument list is spelled out verbatim rather than macro-folded.

/*[AutoHeaderGeneratorExport]*/
void fused_projection_bwd_optimizer_3dgs(
    const int64_t num_splats,
    const int max_sh_degree,
    std::vector<DeviceTensorFloatND> splats_world,
    TorchTensorView viewmats,
    TorchTensorView intrins,
    const uint32_t image_width,
    const uint32_t image_height,
    const std::string camera_model,
    const TorchTensorView dist_coeffs,
    const DeviceVector<int32_t> camera_ids,
    const DeviceVector<int32_t> gaussian_ids,
    DeviceTensorFloatND aabb,
    const std::vector<DeviceTensorFloatND> v_splats_world,
    const std::optional<std::vector<DeviceTensorFloatND>> vr_splats_world,
    const std::optional<std::vector<DeviceTensorFloatND>> h_splats_world,
    const std::vector<DeviceTensorFloatND> v_splats_screen,
    const std::optional<std::vector<DeviceTensorFloatND>> vr_splats_screen,
    const std::optional<std::vector<DeviceTensorFloatND>> h_splats_screen,
    const std::vector<DeviceTensorFloatND> g1_splats_world,
    const std::vector<DeviceTensorFloatND> g2_splats_world,
    const std::optional<TorchTensorView> sh_packed,
    const std::optional<TorchTensorView> sh_quant_bounds,
    const std::optional<TorchTensorView> sh_value_packed,
    const std::optional<TorchTensorView> sh_value_bounds,
    DeviceVector<float> radii,
    const float lr_means,
    const float lr_quats,
    const float lr_scales,
    const float lr_opacs,
    const float lr_features_dc,
    const float lr_features_sh,
    const float max_gauss_ratio,
    const float scale_regularization_weight,
    const float mcmc_opacity_reg_weight,
    const float mcmc_scale_reg_weight,
    const float erank_reg_weight,
    const float erank_reg_weight_s3,
    const float quat_norm_reg_weight,
    const float sh_reg_weight,
    bool use_scale_agnostic_mean,
    bool use_color_trust_region,
    bool color_is_linear,
    float eps_tr,
    std::variant<int32_t, TorchTensorView> step,
    int sh_quant_bits,
    int sh_value_bits
) {
    _fused_projection_bwd_optimizer_dispatch<Vanilla3DGS>(
        num_splats, max_sh_degree, splats_world, viewmats, intrins, image_width, image_height,
        camera_model, dist_coeffs, camera_ids, gaussian_ids, aabb,
        v_splats_world, vr_splats_world, h_splats_world,
        v_splats_screen, vr_splats_screen, h_splats_screen,
        g1_splats_world, g2_splats_world, sh_packed, sh_quant_bounds,
        sh_value_packed, sh_value_bounds,
        radii, lr_means, lr_quats, lr_scales, lr_opacs, lr_features_dc,
        lr_features_sh, max_gauss_ratio, scale_regularization_weight,
        mcmc_opacity_reg_weight, mcmc_scale_reg_weight,
        erank_reg_weight, erank_reg_weight_s3, quat_norm_reg_weight,
        sh_reg_weight, use_scale_agnostic_mean,
        use_color_trust_region, color_is_linear, eps_tr, step,
        sh_quant_bits, sh_value_bits);
}

/*[AutoHeaderGeneratorExport]*/
void fused_projection_bwd_optimizer_mip(
    const int64_t num_splats,
    const int max_sh_degree,
    std::vector<DeviceTensorFloatND> splats_world,
    TorchTensorView viewmats,
    TorchTensorView intrins,
    const uint32_t image_width,
    const uint32_t image_height,
    const std::string camera_model,
    const TorchTensorView dist_coeffs,
    const DeviceVector<int32_t> camera_ids,
    const DeviceVector<int32_t> gaussian_ids,
    DeviceTensorFloatND aabb,
    const std::vector<DeviceTensorFloatND> v_splats_world,
    const std::optional<std::vector<DeviceTensorFloatND>> vr_splats_world,
    const std::optional<std::vector<DeviceTensorFloatND>> h_splats_world,
    const std::vector<DeviceTensorFloatND> v_splats_screen,
    const std::optional<std::vector<DeviceTensorFloatND>> vr_splats_screen,
    const std::optional<std::vector<DeviceTensorFloatND>> h_splats_screen,
    const std::vector<DeviceTensorFloatND> g1_splats_world,
    const std::vector<DeviceTensorFloatND> g2_splats_world,
    const std::optional<TorchTensorView> sh_packed,
    const std::optional<TorchTensorView> sh_quant_bounds,
    const std::optional<TorchTensorView> sh_value_packed,
    const std::optional<TorchTensorView> sh_value_bounds,
    DeviceVector<float> radii,
    const float lr_means,
    const float lr_quats,
    const float lr_scales,
    const float lr_opacs,
    const float lr_features_dc,
    const float lr_features_sh,
    const float max_gauss_ratio,
    const float scale_regularization_weight,
    const float mcmc_opacity_reg_weight,
    const float mcmc_scale_reg_weight,
    const float erank_reg_weight,
    const float erank_reg_weight_s3,
    const float quat_norm_reg_weight,
    const float sh_reg_weight,
    bool use_scale_agnostic_mean,
    bool use_color_trust_region,
    bool color_is_linear,
    float eps_tr,
    std::variant<int32_t, TorchTensorView> step,
    int sh_quant_bits,
    int sh_value_bits
) {
    _fused_projection_bwd_optimizer_dispatch<MipSplatting>(
        num_splats, max_sh_degree, splats_world, viewmats, intrins, image_width, image_height,
        camera_model, dist_coeffs, camera_ids, gaussian_ids, aabb,
        v_splats_world, vr_splats_world, h_splats_world,
        v_splats_screen, vr_splats_screen, h_splats_screen,
        g1_splats_world, g2_splats_world, sh_packed, sh_quant_bounds,
        sh_value_packed, sh_value_bounds,
        radii, lr_means, lr_quats, lr_scales, lr_opacs, lr_features_dc,
        lr_features_sh, max_gauss_ratio, scale_regularization_weight,
        mcmc_opacity_reg_weight, mcmc_scale_reg_weight,
        erank_reg_weight, erank_reg_weight_s3, quat_norm_reg_weight,
        sh_reg_weight, use_scale_agnostic_mean,
        use_color_trust_region, color_is_linear, eps_tr, step,
        sh_quant_bits, sh_value_bits);
}

/*[AutoHeaderGeneratorExport]*/
void fused_projection_bwd_optimizer_3dgut(
    const int64_t num_splats,
    const int max_sh_degree,
    std::vector<DeviceTensorFloatND> splats_world,
    TorchTensorView viewmats,
    TorchTensorView intrins,
    const uint32_t image_width,
    const uint32_t image_height,
    const std::string camera_model,
    const TorchTensorView dist_coeffs,
    const DeviceVector<int32_t> camera_ids,
    const DeviceVector<int32_t> gaussian_ids,
    DeviceTensorFloatND aabb,
    const std::vector<DeviceTensorFloatND> v_splats_world,
    const std::optional<std::vector<DeviceTensorFloatND>> vr_splats_world,
    const std::optional<std::vector<DeviceTensorFloatND>> h_splats_world,
    const std::vector<DeviceTensorFloatND> v_splats_screen,
    const std::optional<std::vector<DeviceTensorFloatND>> vr_splats_screen,
    const std::optional<std::vector<DeviceTensorFloatND>> h_splats_screen,
    const std::vector<DeviceTensorFloatND> g1_splats_world,
    const std::vector<DeviceTensorFloatND> g2_splats_world,
    const std::optional<TorchTensorView> sh_packed,
    const std::optional<TorchTensorView> sh_quant_bounds,
    const std::optional<TorchTensorView> sh_value_packed,
    const std::optional<TorchTensorView> sh_value_bounds,
    DeviceVector<float> radii,
    const float lr_means,
    const float lr_quats,
    const float lr_scales,
    const float lr_opacs,
    const float lr_features_dc,
    const float lr_features_sh,
    const float max_gauss_ratio,
    const float scale_regularization_weight,
    const float mcmc_opacity_reg_weight,
    const float mcmc_scale_reg_weight,
    const float erank_reg_weight,
    const float erank_reg_weight_s3,
    const float quat_norm_reg_weight,
    const float sh_reg_weight,
    bool use_scale_agnostic_mean,
    bool use_color_trust_region,
    bool color_is_linear,
    float eps_tr,
    std::variant<int32_t, TorchTensorView> step,
    int sh_quant_bits,
    int sh_value_bits
) {
    _fused_projection_bwd_optimizer_dispatch<Vanilla3DGUT>(
        num_splats, max_sh_degree, splats_world, viewmats, intrins, image_width, image_height,
        camera_model, dist_coeffs, camera_ids, gaussian_ids, aabb,
        v_splats_world, vr_splats_world, h_splats_world,
        v_splats_screen, vr_splats_screen, h_splats_screen,
        g1_splats_world, g2_splats_world, sh_packed, sh_quant_bounds,
        sh_value_packed, sh_value_bounds,
        radii, lr_means, lr_quats, lr_scales, lr_opacs, lr_features_dc,
        lr_features_sh, max_gauss_ratio, scale_regularization_weight,
        mcmc_opacity_reg_weight, mcmc_scale_reg_weight,
        erank_reg_weight, erank_reg_weight_s3, quat_norm_reg_weight,
        sh_reg_weight, use_scale_agnostic_mean,
        use_color_trust_region, color_is_linear, eps_tr, step,
        sh_quant_bits, sh_value_bits);
}


