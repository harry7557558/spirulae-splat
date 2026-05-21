#include "ProjectionBwd.cuh"

#include <gsplat/Utils.cuh>

#include <c10/cuda/CUDAStream.h>
#include <cooperative_groups.h>
namespace cg = cooperative_groups;

#include <cub/cub.cuh>


template<
    typename SplatPrimitive,
    ssplat::CameraModelType camera_model,
    HessianDiagonalOutputMode hessian_diagonal_output_mode,
    const bool use_scale_agnostic_mean
>
void fused_projection_bwd_optimizer_3dgs_kernel_wrapper(
    cudaStream_t stream,
    // fwd inputs
    const uint32_t C,
    const uint32_t N,
    typename SplatPrimitive::WorldBuffer splats_world,
    const float *__restrict__ viewmats, // [C, 4, 4]
    const float4 *__restrict__ intrins,  // [C, 4], fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const uint32_t image_width,
    const uint32_t image_height,
    // fwd outputs
    const int32_t *__restrict__ camera_id_bounds,   // [N]
    const int32_t *__restrict__ camera_ids,   // [nnz]
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
    const uint8_t* __restrict__ g1_features_sh,
    const uint8_t* __restrict__ g2_features_sh,
    float4* __restrict__ sh_quant_bounds,
    // float *__restrict__ v_viewmats // [C, 4, 4] optional
    // optimizer params
    const float* __restrict__ radii,
    const float lr_means,
    const float lr_quats,
    const float lr_scales,
    const float lr_opacs,
    const float lr_features_dc,
    const float lr_features_sh,
    const float mcmc_noise_scalar,
    const float min_opacity,
    const float max_gauss_ratio,
    const float scale_regularization_weight,
    const float mcmc_opacity_reg_weight,
    const float mcmc_scale_reg_weight,
    const float erank_reg_weight,
    const float erank_reg_weight_s3,
    const float quat_norm_reg_weight,
    const float sh_reg_weight,
    const float mrnf_opacity_decay_factor,
    const float mrnf_scale_decay_factor,
    const int32_t scalar_step,
    const int32_t* __restrict__ steps
);


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
    bool use_scale_agnostic_mean
>
inline void launch_fused_projection_bwd_optimizer_3dgs_kernel(
    // fwd inputs
    const int64_t N,
    TensorList splats_world,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const ssplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    std::optional<at::Tensor> camera_ids,
    std::optional<at::Tensor> gaussian_ids,
    const at::Tensor aabb,
    // grad outputs
    const TensorList v_splats_world,
    const std::optional<TensorList> vr_splats_world,
    const std::optional<TensorList> h_splats_world,
    const TensorList v_splats_screen,
    const std::optional<TensorList> vr_splats_screen,
    const std::optional<TensorList> h_splats_screen,
    // optimizer states
    const TensorList g1_splats_world,
    const TensorList g2_splats_world,
    const std::optional<at::Tensor> g1_features_sh,
    const std::optional<at::Tensor> g2_features_sh,
    const std::optional<at::Tensor> sh_quant_bounds,
    // optimizer params
    const at::Tensor radii,
    const float lr_means,
    const float lr_quats,
    const float lr_scales,
    const float lr_opacs,
    const float lr_features_dc,
    const float lr_features_sh,
    const float mcmc_noise_scalar,
    const float min_opacity,
    const float max_gauss_ratio,
    const float scale_regularization_weight,
    const float mcmc_opacity_reg_weight,
    const float mcmc_scale_reg_weight,
    const float erank_reg_weight,
    const float erank_reg_weight_s3,
    const float quat_norm_reg_weight,
    const float sh_reg_weight,
    const float mrnf_opacity_decay_factor,
    const float mrnf_scale_decay_factor,
    const int32_t scalar_step,
    const std::optional<at::Tensor> steps
) {
    uint32_t C = viewmats.size(-3); // number of cameras

    if (N == 0)
        return;

    at::Tensor camera_id_bounds;

    bool packed = camera_ids.has_value() && gaussian_ids.has_value();
    if (packed) {
        at::Tensor camera_ids_sorted = at::empty_like(camera_ids.value());
        at::Tensor gaussian_ids_sorted = at::empty_like(gaussian_ids.value());
        cub::DoubleBuffer<int32_t> d_keys(
            gaussian_ids.value().data_ptr<int32_t>(), gaussian_ids_sorted.data_ptr<int32_t>()
        );
        cub::DoubleBuffer<int32_t> d_values(
            camera_ids.value().data_ptr<int32_t>(), camera_ids_sorted.data_ptr<int32_t>()
        );
        long nnz = camera_ids.value().numel();
        int n_bits = 0;
        while ((1U << n_bits) <= N)
            ++n_bits;
        CUB_WRAPPER(
            cub::DeviceRadixSort::SortPairs,
            d_keys,
            d_values,
            nnz,
            0,
            n_bits,
            at::cuda::getCurrentCUDAStream()
        );
        CHECK_DEVICE_ERROR(cudaGetLastError());
        switch (d_keys.selector) {
        case 0:
            break;
        case 1:
            gaussian_ids = gaussian_ids_sorted;
            break;
        }
        switch (d_values.selector) {
        case 0:
            break;
        case 1:
            camera_ids = camera_ids_sorted;
            break;
        }

        camera_id_bounds = at::empty({N+1}, kTensorOptionI32());
        camera_id_bounds_kernel<<<_LAUNCH_ARGS_1D(nnz+1, 256)>>>(
            nnz, N,
            gaussian_ids.value().data_ptr<int32_t>(),
            camera_id_bounds.data_ptr<int32_t>()
        );
        CHECK_DEVICE_ERROR(cudaGetLastError());
    }

    #define _LAUNCH_ARGS ( \
            (cudaStream_t)at::cuda::getCurrentCUDAStream(), C, N, \
            splats_world, viewmats.data_ptr<float>(), (float4*)intrins.data_ptr<float>(), dist_coeffs, \
            image_width, image_height, \
            packed ? camera_id_bounds.data_ptr<int32_t>() : nullptr, \
            packed ? camera_ids.value().data_ptr<int32_t>() : nullptr, \
            (float4*)aabb.data_ptr<float>(), \
            v_splats_world, \
            hessian_diagonal_output_mode != HessianDiagonalOutputMode::None ? vr_splats_world.value() : typename SplatPrimitive::WorldBuffer{}, \
            hessian_diagonal_output_mode != HessianDiagonalOutputMode::None ? h_splats_world.value() : typename SplatPrimitive::WorldBuffer{}, \
            v_splats_screen, \
            hessian_diagonal_output_mode != HessianDiagonalOutputMode::None ? vr_splats_screen.value() : typename SplatPrimitive::ScreenBuffer{}, \
            hessian_diagonal_output_mode != HessianDiagonalOutputMode::None ? h_splats_screen.value() : typename SplatPrimitive::ScreenBuffer{}, \
            g1_splats_world, g2_splats_world, \
            g1_features_sh.has_value() ? g1_features_sh.value().data_ptr<uint8_t>() : nullptr, \
            g2_features_sh.has_value() ? g2_features_sh.value().data_ptr<uint8_t>() : nullptr, \
            sh_quant_bounds.has_value() ? (float4*)sh_quant_bounds.value().data_ptr<float>() : nullptr, \
            /*v_viewmats.has_value() ? v_viewmats.value().data_ptr<float>() : nullptr */ \
            radii.data_ptr<float>(), \
            lr_means, lr_quats, lr_scales, lr_opacs, lr_features_dc, lr_features_sh, \
            mcmc_noise_scalar, min_opacity, max_gauss_ratio, scale_regularization_weight, \
            mcmc_opacity_reg_weight, mcmc_scale_reg_weight, erank_reg_weight, erank_reg_weight_s3, quat_norm_reg_weight, sh_reg_weight, \
            mrnf_opacity_decay_factor, mrnf_scale_decay_factor, \
            scalar_step, steps.has_value() ? steps.value().data_ptr<int32_t>() : nullptr \
        )

    if (camera_model == ssplat::CameraModelType::PINHOLE)
        fused_projection_bwd_optimizer_3dgs_kernel_wrapper<SplatPrimitive, ssplat::CameraModelType::PINHOLE, hessian_diagonal_output_mode, use_scale_agnostic_mean> _LAUNCH_ARGS;
    else if (camera_model == ssplat::CameraModelType::FISHEYE)
        fused_projection_bwd_optimizer_3dgs_kernel_wrapper<SplatPrimitive, ssplat::CameraModelType::FISHEYE, hessian_diagonal_output_mode, use_scale_agnostic_mean> _LAUNCH_ARGS;
    else if (camera_model == ssplat::CameraModelType::EQUISOLID)
        fused_projection_bwd_optimizer_3dgs_kernel_wrapper<SplatPrimitive, ssplat::CameraModelType::EQUISOLID, hessian_diagonal_output_mode, use_scale_agnostic_mean> _LAUNCH_ARGS;
    else
        throw std::runtime_error("Unsupported camera model");
    CHECK_DEVICE_ERROR(cudaGetLastError());

    #undef _LAUNCH_ARGS

}



// // ================
// // Vanilla3DGS
// // ================

// /*[AutoHeaderGeneratorExport]*/
// void projection_3dgs_backward_tensor(
//     // fwd inputs
//     const TensorList &splats_world,
//     const at::Tensor viewmats,  // [..., C, 4, 4]
//     const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
//     const uint32_t image_width,
//     const uint32_t image_height,
//     const std::string camera_model,
//     const CameraDistortionCoeffsTensor dist_coeffs,
//     // fwd outputs
//     const std::optional<at::Tensor> camera_ids,  // [nnz]
//     const std::optional<at::Tensor> gaussian_ids,  // [nnz]
//     const at::Tensor aabb,  // [..., C, N, 2]
//     // grad outputs
//     const TensorList &v_splats_screen,
//     // returns
//     const TensorList &v_splats_world,
//     const std::optional<at::Tensor> &v_viewmats
// ) {
//     launch_projection_projection_fused_bwd_kernel<Vanilla3DGS, HessianDiagonalOutputMode::None>(
//         splats_world, viewmats, intrins, image_width, image_height, cmt(camera_model), dist_coeffs,
//         camera_ids, gaussian_ids, aabb, v_splats_screen, {}, {},
//         v_splats_world, v_viewmats, std::nullopt, std::nullopt);
// }

// /*[AutoHeaderGeneratorExport]*/
// void projection_3dgs_backward_with_hessian_diagonal_tensor(
//     // fwd inputs
//     const TensorList &splats_world,
//     const at::Tensor viewmats,  // [..., C, 4, 4]
//     const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
//     const uint32_t image_width,
//     const uint32_t image_height,
//     const std::string camera_model,
//     const CameraDistortionCoeffsTensor dist_coeffs,
//     // fwd outputs
//     const std::optional<at::Tensor> camera_ids,  // [nnz]
//     const std::optional<at::Tensor> gaussian_ids,  // [nnz]
//     const at::Tensor aabb,                       // [..., C, N, 2]
//     // grad outputs
//     const TensorList &v_splats_screen,
//     const TensorList &vr_splats_screen,
//     const TensorList &h_splats_screen,
//     // returns
//     const TensorList &v_splats_world,
//     const std::optional<at::Tensor> &v_viewmats,
//     const TensorList &vr_splats_world,
//     const TensorList &h_splats_world
// ) {
//     launch_projection_projection_fused_bwd_kernel<Vanilla3DGS, HessianDiagonalOutputMode::AllReasonable>(
//         splats_world, viewmats, intrins, image_width, image_height, cmt(camera_model), dist_coeffs,
//         camera_ids, gaussian_ids, aabb, v_splats_screen, vr_splats_screen, h_splats_screen,
//         v_splats_world, v_viewmats, vr_splats_world, h_splats_world);
// }

// /*[AutoHeaderGeneratorExport]*/
// void projection_3dgs_backward_with_position_hessian_diagonal_tensor(
//     // fwd inputs
//     const TensorList &splats_world,
//     const at::Tensor viewmats,  // [..., C, 4, 4]
//     const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
//     const uint32_t image_width,
//     const uint32_t image_height,
//     const std::string camera_model,
//     const CameraDistortionCoeffsTensor dist_coeffs,
//     // fwd outputs
//     const std::optional<at::Tensor> camera_ids,  // [nnz]
//     const std::optional<at::Tensor> gaussian_ids,  // [nnz]
//     const at::Tensor aabb,                       // [..., C, N, 2]
//     // grad outputs
//     const TensorList &v_splats_screen,
//     const TensorList &vr_splats_screen,
//     const TensorList &h_splats_screen,
//     // returns
//     const TensorList &v_splats_world,
//     const std::optional<at::Tensor> &v_viewmats,
//     const at::Tensor &vr_splats_world,
//     const at::Tensor &h_splats_world
// ) {
//     launch_projection_projection_fused_bwd_kernel<Vanilla3DGS, HessianDiagonalOutputMode::Position>(
//         splats_world, viewmats, intrins, image_width, image_height, cmt(camera_model), dist_coeffs,
//         camera_ids, gaussian_ids, aabb, v_splats_screen, vr_splats_screen, h_splats_screen,
//         v_splats_world, v_viewmats, vr_splats_world, h_splats_world);
// }



// // ================
// // MipSplatting
// // ================

// /*[AutoHeaderGeneratorExport]*/
// void projection_mip_backward_tensor(
//     // fwd inputs
//     const TensorList &splats_world,
//     const at::Tensor viewmats,  // [..., C, 4, 4]
//     const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
//     const uint32_t image_width,
//     const uint32_t image_height,
//     const std::string camera_model,
//     const CameraDistortionCoeffsTensor dist_coeffs,
//     // fwd outputs
//     const std::optional<at::Tensor> camera_ids,  // [nnz]
//     const std::optional<at::Tensor> gaussian_ids,  // [nnz]
//     const at::Tensor aabb,                       // [..., C, N, 2]
//     // grad outputs
//     const TensorList &v_splats_screen,
//     // returns
//     const TensorList &v_splats_world,
//     const std::optional<at::Tensor> &v_viewmats
// ) {
//     launch_projection_projection_fused_bwd_kernel<MipSplatting, HessianDiagonalOutputMode::None>(
//         splats_world, viewmats, intrins, image_width, image_height, cmt(camera_model), dist_coeffs,
//         camera_ids, gaussian_ids, aabb, v_splats_screen, {}, {},
//         v_splats_world, v_viewmats, std::nullopt, std::nullopt);
// }

// /*[AutoHeaderGeneratorExport]*/
// void projection_mip_backward_with_hessian_diagonal_tensor(
//     // fwd inputs
//     const TensorList &splats_world,
//     const at::Tensor viewmats,  // [..., C, 4, 4]
//     const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
//     const uint32_t image_width,
//     const uint32_t image_height,
//     const std::string camera_model,
//     const CameraDistortionCoeffsTensor dist_coeffs,
//     // fwd outputs
//     const std::optional<at::Tensor> camera_ids,  // [nnz]
//     const std::optional<at::Tensor> gaussian_ids,  // [nnz]
//     const at::Tensor aabb,                       // [..., C, N, 2]
//     // grad outputs
//     const TensorList &v_splats_screen,
//     const TensorList &vr_splats_screen,
//     const TensorList &h_splats_screen,
//     // returns
//     const TensorList &v_splats_world,
//     const std::optional<at::Tensor> &v_viewmats,
//     const TensorList &vr_splats_world,
//     const TensorList &h_splats_world
// ) {
//     launch_projection_projection_fused_bwd_kernel<MipSplatting, HessianDiagonalOutputMode::AllReasonable>(
//         splats_world, viewmats, intrins, image_width, image_height, cmt(camera_model), dist_coeffs,
//         camera_ids, gaussian_ids, aabb, v_splats_screen, vr_splats_screen, h_splats_screen,
//         v_splats_world, v_viewmats, vr_splats_world, h_splats_world);
// }

// /*[AutoHeaderGeneratorExport]*/
// void projection_mip_backward_with_position_hessian_diagonal_tensor(
//     // fwd inputs
//     const TensorList &splats_world,
//     const at::Tensor viewmats,  // [..., C, 4, 4]
//     const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
//     const uint32_t image_width,
//     const uint32_t image_height,
//     const std::string camera_model,
//     const CameraDistortionCoeffsTensor dist_coeffs,
//     // fwd outputs
//     const std::optional<at::Tensor> camera_ids,  // [nnz]
//     const std::optional<at::Tensor> gaussian_ids,  // [nnz]
//     const at::Tensor aabb,                       // [..., C, N, 2]
//     // grad outputs
//     const TensorList &v_splats_screen,
//     const TensorList &vr_splats_screen,
//     const TensorList &h_splats_screen,
//     // returns
//     const TensorList &v_splats_world,
//     const std::optional<at::Tensor> &v_viewmats,
//     const at::Tensor &vr_splats_world,
//     const at::Tensor &h_splats_world
// ) {
//     launch_projection_projection_fused_bwd_kernel<MipSplatting, HessianDiagonalOutputMode::Position>(
//         splats_world, viewmats, intrins, image_width, image_height, cmt(camera_model), dist_coeffs,
//         camera_ids, gaussian_ids, aabb, v_splats_screen, vr_splats_screen, h_splats_screen,
//         v_splats_world, v_viewmats, vr_splats_world, h_splats_world);
// }


// ================
// Vanilla3DGUT
// ================

/*[AutoHeaderGeneratorExport]*/
void fused_projection_bwd_optimizer_3dgut_tensor(
    // fwd inputs
    const int64_t num_splats,
    TensorList splats_world,
    const at::Tensor viewmats,
    const at::Tensor intrins,
    const uint32_t image_width,
    const uint32_t image_height,
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const std::optional<at::Tensor> camera_ids,
    const std::optional<at::Tensor> gaussian_ids,
    const at::Tensor aabb,
    // grad outputs
    const TensorList v_splats_world,
    const std::optional<TensorList> vr_splats_world,
    const std::optional<TensorList> h_splats_world,
    const TensorList v_splats_screen,
    const std::optional<TensorList> vr_splats_screen,
    const std::optional<TensorList> h_splats_screen,
    // optimizer states
    const TensorList g1_splats_world,
    const TensorList g2_splats_world,
    const std::optional<at::Tensor> g1_features_sh,
    const std::optional<at::Tensor> g2_features_sh,
    const std::optional<at::Tensor> sh_quant_bounds,
    // optimizer params
    const at::Tensor radii,
    const float lr_means,
    const float lr_quats,
    const float lr_scales,
    const float lr_opacs,
    const float lr_features_dc,
    const float lr_features_sh,
    const float mcmc_noise_scalar,
    const float min_opacity,
    const float max_gauss_ratio,
    const float scale_regularization_weight,
    const float mcmc_opacity_reg_weight,
    const float mcmc_scale_reg_weight,
    const float erank_reg_weight,
    const float erank_reg_weight_s3,
    const float quat_norm_reg_weight,
    const float sh_reg_weight,
    const float mrnf_opacity_decay_factor,
    const float mrnf_scale_decay_factor,
    bool use_scale_agnostic_mean,
    std::variant<int32_t, at::Tensor> step
) {
    #define _ARGS ( \
        num_splats, \
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
        g1_features_sh, \
        g2_features_sh, \
        sh_quant_bounds, \
        radii, \
        lr_means, \
        lr_quats, \
        lr_scales, \
        lr_opacs, \
        lr_features_dc, \
        lr_features_sh, \
        mcmc_noise_scalar, \
        min_opacity, \
        max_gauss_ratio, \
        scale_regularization_weight, \
        mcmc_opacity_reg_weight, \
        mcmc_scale_reg_weight, \
        erank_reg_weight, \
        erank_reg_weight_s3, \
        quat_norm_reg_weight, \
        sh_reg_weight, \
        mrnf_opacity_decay_factor, \
        mrnf_scale_decay_factor, \
        std::get_if<int32_t>(&step) ? std::get<int32_t>(step) : -1, \
        std::get_if<at::Tensor>(&step) ? std::get<at::Tensor>(step) : (std::optional<at::Tensor>)std::nullopt \
    )
    int sh_degree = Vanilla3DGUT<0>::WorldBuffer(splats_world).sh_degree();
    #define LAUNCH(n) if (sh_degree == (n)) \
        (use_scale_agnostic_mean ? \
            launch_fused_projection_bwd_optimizer_3dgs_kernel<Vanilla3DGUT<n>, HessianDiagonalOutputMode::None, true> : \
            launch_fused_projection_bwd_optimizer_3dgs_kernel<Vanilla3DGUT<n>, HessianDiagonalOutputMode::None, false> \
        ) _ARGS;
    LAUNCH(3) LAUNCH(2) LAUNCH(1) LAUNCH(0) LAUNCH(4)
    #undef LAUNCH
    #undef _ARGS
}



// // ================
// // SphericalVoronoi3DGUT
// // ================

// /*[AutoHeaderGeneratorExport]*/
// void projection_3dgut_sv_backward_tensor(
//     // fwd inputs
//     const SphericalVoronoi3DGUT_Default::World::TensorTuple &splats_world,
//     const at::Tensor viewmats,  // [..., C, 4, 4]
//     const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
//     const uint32_t image_width,
//     const uint32_t image_height,
//     const std::string camera_model,
//     const CameraDistortionCoeffsTensor dist_coeffs,
//     // fwd outputs
//     const std::optional<at::Tensor> camera_ids,  // [nnz]
//     const std::optional<at::Tensor> gaussian_ids,  // [nnz]
//     const at::Tensor aabb,                       // [..., C, N, 2]
//     // grad outputs
//     const SphericalVoronoi3DGUT_Default::Screen::TensorTupleProj &v_splats_screen,
//     // returns
//     const SphericalVoronoi3DGUT_Default::World::TensorTuple &v_splats_world,
//     const std::optional<at::Tensor> &v_viewmats
// ) {
//     int num_sv = std::get<5>(splats_world).size(-2);
//     #define _CASE(n) \
//         if (num_sv == n) launch_projection_projection_fused_bwd_kernel<SphericalVoronoi3DGUT<n>, HessianDiagonalOutputMode::None>( \
//             splats_world, viewmats, intrins, image_width, image_height, cmt(camera_model), dist_coeffs, \
//             camera_ids, gaussian_ids, aabb, v_splats_screen, nullptr, nullptr, \
//             v_splats_world, v_viewmats, std::nullopt, std::nullopt);
//     _CASE(2) _CASE(3) _CASE(4) _CASE(5) _CASE(6) _CASE(7) _CASE(8)
//     #undef _CASE
//     throw std::invalid_argument("Unsupported num_sv");
// }



// // ================
// // OpaqueTriangle
// // ================


// /*[AutoHeaderGeneratorExport]*/
// void projection_opaque_triangle_backward_tensor(
//     // fwd inputs
//     const OpaqueTriangle::World::TensorTuple &splats_world,
//     const at::Tensor viewmats,  // [..., C, 4, 4]
//     const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
//     const uint32_t image_width,
//     const uint32_t image_height,
//     const std::string camera_model,
//     const CameraDistortionCoeffsTensor dist_coeffs,
//     // fwd outputs
//     const std::optional<at::Tensor> camera_ids,  // [nnz]
//     const std::optional<at::Tensor> gaussian_ids,  // [nnz]
//     const at::Tensor aabb,                       // [..., C, N, 2]
//     // grad outputs
//     const OpaqueTriangle::Screen::TensorTupleProj &v_splats_screen,
//     // returns
//     const OpaqueTriangle::World::TensorTuple &v_splats_world,
//     const std::optional<at::Tensor> &v_viewmats
// ) {
//     launch_projection_projection_fused_bwd_kernel<OpaqueTriangle, HessianDiagonalOutputMode::None>(
//         splats_world, viewmats, intrins, image_width, image_height, cmt(camera_model), dist_coeffs,
//         camera_ids, gaussian_ids, aabb, v_splats_screen, nullptr, nullptr,
//         v_splats_world, v_viewmats, std::nullopt, std::nullopt);
// }


// // ================
// // VoxelPrimitive
// // ================


// /*[AutoHeaderGeneratorExport]*/
// void projection_voxel_backward_tensor(
//     // fwd inputs
//     const VoxelPrimitive::World::TensorTuple &splats_world,
//     const at::Tensor viewmats,  // [..., C, 4, 4]
//     const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
//     const uint32_t image_width,
//     const uint32_t image_height,
//     const std::string camera_model,
//     const CameraDistortionCoeffsTensor dist_coeffs,
//     // fwd outputs
//     const std::optional<at::Tensor> camera_ids,  // [nnz]
//     const std::optional<at::Tensor> gaussian_ids,  // [nnz]
//     const at::Tensor aabb,                       // [..., C, N, 2]
//     // grad outputs
//     const VoxelPrimitive::Screen::TensorTupleProj &v_splats_screen,
//     // returns
//     const VoxelPrimitive::World::TensorTuple &v_splats_world,
//     const std::optional<at::Tensor> &v_viewmats
// ) {
//     launch_projection_projection_fused_bwd_kernel<VoxelPrimitive, HessianDiagonalOutputMode::None>(
//         splats_world, viewmats, intrins, image_width, image_height, cmt(camera_model), dist_coeffs,
//         camera_ids, gaussian_ids, aabb, v_splats_screen, nullptr, nullptr,
//         v_splats_world, v_viewmats, std::nullopt, std::nullopt);
// }

