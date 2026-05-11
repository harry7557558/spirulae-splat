#include "ProjectionBwd.cuh"

#include <gsplat/Utils.cuh>

#include <c10/cuda/CUDAStream.h>
#include <cooperative_groups.h>
namespace cg = cooperative_groups;


template<
    typename SplatPrimitive,
    ssplat::CameraModelType camera_model,
    HessianDiagonalOutputMode hessian_diagonal_output_mode
>
void projection_fused_bwd_kernel_wrapper(
    cudaStream_t stream,
    // fwd inputs
    const uint32_t B,
    const uint32_t C,
    const uint32_t N,
    const typename SplatPrimitive::WorldBuffer splats_world,
    const float * viewmats, // [B, C, 4, 4]
    const float4 * intrins,  // [B, C, 4], fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const uint32_t image_width,
    const uint32_t image_height,
    // fwd outputs
    const int32_t * camera_ids,          // [nnz, 4]
    const int32_t * gaussian_ids,          // [nnz, 4]
    const float4 * aabb,          // [B, C, N, 4]
    // grad outputs
    typename SplatPrimitive::ScreenBuffer v_splats_screen,
    typename SplatPrimitive::ScreenBuffer vr_splats_screen,
    typename SplatPrimitive::ScreenBuffer h_splats_screen,
    // grad inputs
    typename SplatPrimitive::WorldBuffer v_splats_world,
    float3* vr_world_pos_buffer,
    float3* h_world_pos_buffer,
    typename SplatPrimitive::WorldBuffer vr_splats_world,
    typename SplatPrimitive::WorldBuffer h_splats_world,
    float * v_viewmats // [B, C, 4, 4] optional
);


template<typename SplatPrimitive, HessianDiagonalOutputMode hessian_diagonal_output_mode>
inline void launch_projection_projection_fused_bwd_kernel(
    // fwd inputs
    const TensorList &splats_world_tuple,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const ssplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const std::optional<at::Tensor> camera_ids,          // [nnz, 4]
    const std::optional<at::Tensor> gaussian_ids,          // [nnz, 4]
    const at::Tensor aabb,                       // [..., C, N, 2]
    // grad outputs
    const TensorList v_splats_screen,
    const TensorList vr_splats_screen,
    const TensorList h_splats_screen,
    // returns
    TensorList v_splats_world,
    const std::optional<at::Tensor> v_viewmats,
    const std::optional<std::variant<at::Tensor, TensorList>> vr_splats_world_tuple,
    const std::optional<std::variant<at::Tensor, TensorList>> h_splats_world_tuple
) {
    typename SplatPrimitive::WorldBuffer splats_world(splats_world_tuple);
    uint32_t N = splats_world.size();    // number of gaussians
    uint32_t C = viewmats.size(-3); // number of cameras
    // uint32_t B = splats_world.batchSize();    // number of batches
    uint32_t B = 1;  // TODO

    at::Tensor vr_world_pos;
    at::Tensor h_world_pos;
    if (hessian_diagonal_output_mode == HessianDiagonalOutputMode::Position) {
        vr_world_pos = std::get<0>(vr_splats_world_tuple.value());
        h_world_pos = std::get<0>(h_splats_world_tuple.value());
    }
    typename SplatPrimitive::WorldBuffer vr_splats_world;
    typename SplatPrimitive::WorldBuffer h_splats_world;
    if (hessian_diagonal_output_mode == HessianDiagonalOutputMode::AllReasonable) {
        vr_splats_world = std::get<1>(vr_splats_world_tuple.value());
        h_splats_world = std::get<1>(h_splats_world_tuple.value());
    }

    if (camera_ids.has_value() && gaussian_ids.has_value()) {  // packed
        N = camera_ids.value().numel();
        B = C = 1;
    }

    if (B*C*N == 0)
        return;

    #define _LAUNCH_ARGS ( \
            (cudaStream_t)at::cuda::getCurrentCUDAStream(), B, C, N, \
            splats_world, viewmats.data_ptr<float>(), (float4*)intrins.data_ptr<float>(), dist_coeffs, \
            image_width, image_height, \
            camera_ids.has_value() ? camera_ids.value().data_ptr<int32_t>() : nullptr, \
            gaussian_ids.has_value() ? gaussian_ids.value().data_ptr<int32_t>() : nullptr, \
            (float4*)aabb.data_ptr<float>(), \
            v_splats_screen, \
            hessian_diagonal_output_mode != HessianDiagonalOutputMode::None ? vr_splats_screen : typename SplatPrimitive::ScreenBuffer{}, \
            hessian_diagonal_output_mode != HessianDiagonalOutputMode::None ? h_splats_screen : typename SplatPrimitive::ScreenBuffer{}, \
            v_splats_world, \
            hessian_diagonal_output_mode == HessianDiagonalOutputMode::Position ? (float3*)vr_world_pos.data_ptr<float>() : nullptr, \
            hessian_diagonal_output_mode == HessianDiagonalOutputMode::Position ? (float3*)h_world_pos.data_ptr<float>() : nullptr, \
            hessian_diagonal_output_mode == HessianDiagonalOutputMode::AllReasonable ? vr_splats_world : typename SplatPrimitive::WorldBuffer{}, \
            hessian_diagonal_output_mode == HessianDiagonalOutputMode::AllReasonable ? h_splats_world : typename SplatPrimitive::WorldBuffer{}, \
            v_viewmats.has_value() ? v_viewmats.value().data_ptr<float>() : nullptr \
        )

    if (camera_model == ssplat::CameraModelType::PINHOLE)
        projection_fused_bwd_kernel_wrapper<SplatPrimitive, ssplat::CameraModelType::PINHOLE, hessian_diagonal_output_mode> _LAUNCH_ARGS;
    else if (camera_model == ssplat::CameraModelType::FISHEYE)
        projection_fused_bwd_kernel_wrapper<SplatPrimitive, ssplat::CameraModelType::FISHEYE, hessian_diagonal_output_mode> _LAUNCH_ARGS;
    else
        throw std::runtime_error("Unsupported camera model");
    CHECK_DEVICE_ERROR(cudaGetLastError());

    #undef _LAUNCH_ARGS

}



// ================
// Vanilla3DGS
// ================


/*[AutoHeaderGeneratorExport]*/
void projection_3dgs_backward_tensor(
    // fwd inputs
    const TensorList &splats_world,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const std::optional<at::Tensor> camera_ids,  // [nnz]
    const std::optional<at::Tensor> gaussian_ids,  // [nnz]
    const at::Tensor aabb,  // [..., C, N, 2]
    // grad outputs
    const TensorList &v_splats_screen,
    // returns
    const TensorList &v_splats_world,
    const std::optional<at::Tensor> &v_viewmats
) {
    launch_projection_projection_fused_bwd_kernel<Vanilla3DGS, HessianDiagonalOutputMode::None>(
        splats_world, viewmats, intrins, image_width, image_height, cmt(camera_model), dist_coeffs,
        camera_ids, gaussian_ids, aabb, v_splats_screen, {}, {},
        v_splats_world, v_viewmats, std::nullopt, std::nullopt);
}

/*[AutoHeaderGeneratorExport]*/
void projection_3dgs_backward_with_hessian_diagonal_tensor(
    // fwd inputs
    const TensorList &splats_world,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const std::optional<at::Tensor> camera_ids,  // [nnz]
    const std::optional<at::Tensor> gaussian_ids,  // [nnz]
    const at::Tensor aabb,                       // [..., C, N, 2]
    // grad outputs
    const TensorList &v_splats_screen,
    const TensorList &vr_splats_screen,
    const TensorList &h_splats_screen,
    // returns
    const TensorList &v_splats_world,
    const std::optional<at::Tensor> &v_viewmats,
    const TensorList &vr_splats_world,
    const TensorList &h_splats_world
) {
    launch_projection_projection_fused_bwd_kernel<Vanilla3DGS, HessianDiagonalOutputMode::AllReasonable>(
        splats_world, viewmats, intrins, image_width, image_height, cmt(camera_model), dist_coeffs,
        camera_ids, gaussian_ids, aabb, v_splats_screen, vr_splats_screen, h_splats_screen,
        v_splats_world, v_viewmats, vr_splats_world, h_splats_world);
}

/*[AutoHeaderGeneratorExport]*/
void projection_3dgs_backward_with_position_hessian_diagonal_tensor(
    // fwd inputs
    const TensorList &splats_world,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const std::optional<at::Tensor> camera_ids,  // [nnz]
    const std::optional<at::Tensor> gaussian_ids,  // [nnz]
    const at::Tensor aabb,                       // [..., C, N, 2]
    // grad outputs
    const TensorList &v_splats_screen,
    const TensorList &vr_splats_screen,
    const TensorList &h_splats_screen,
    // returns
    const TensorList &v_splats_world,
    const std::optional<at::Tensor> &v_viewmats,
    const at::Tensor &vr_splats_world,
    const at::Tensor &h_splats_world
) {
    launch_projection_projection_fused_bwd_kernel<Vanilla3DGS, HessianDiagonalOutputMode::Position>(
        splats_world, viewmats, intrins, image_width, image_height, cmt(camera_model), dist_coeffs,
        camera_ids, gaussian_ids, aabb, v_splats_screen, vr_splats_screen, h_splats_screen,
        v_splats_world, v_viewmats, vr_splats_world, h_splats_world);
}



// ================
// MipSplatting
// ================

/*[AutoHeaderGeneratorExport]*/
void projection_mip_backward_tensor(
    // fwd inputs
    const TensorList &splats_world,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const std::optional<at::Tensor> camera_ids,  // [nnz]
    const std::optional<at::Tensor> gaussian_ids,  // [nnz]
    const at::Tensor aabb,                       // [..., C, N, 2]
    // grad outputs
    const TensorList &v_splats_screen,
    // returns
    const TensorList &v_splats_world,
    const std::optional<at::Tensor> &v_viewmats
) {
    launch_projection_projection_fused_bwd_kernel<MipSplatting, HessianDiagonalOutputMode::None>(
        splats_world, viewmats, intrins, image_width, image_height, cmt(camera_model), dist_coeffs,
        camera_ids, gaussian_ids, aabb, v_splats_screen, {}, {},
        v_splats_world, v_viewmats, std::nullopt, std::nullopt);
}

/*[AutoHeaderGeneratorExport]*/
void projection_mip_backward_with_hessian_diagonal_tensor(
    // fwd inputs
    const TensorList &splats_world,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const std::optional<at::Tensor> camera_ids,  // [nnz]
    const std::optional<at::Tensor> gaussian_ids,  // [nnz]
    const at::Tensor aabb,                       // [..., C, N, 2]
    // grad outputs
    const TensorList &v_splats_screen,
    const TensorList &vr_splats_screen,
    const TensorList &h_splats_screen,
    // returns
    const TensorList &v_splats_world,
    const std::optional<at::Tensor> &v_viewmats,
    const TensorList &vr_splats_world,
    const TensorList &h_splats_world
) {
    launch_projection_projection_fused_bwd_kernel<MipSplatting, HessianDiagonalOutputMode::AllReasonable>(
        splats_world, viewmats, intrins, image_width, image_height, cmt(camera_model), dist_coeffs,
        camera_ids, gaussian_ids, aabb, v_splats_screen, vr_splats_screen, h_splats_screen,
        v_splats_world, v_viewmats, vr_splats_world, h_splats_world);
}

/*[AutoHeaderGeneratorExport]*/
void projection_mip_backward_with_position_hessian_diagonal_tensor(
    // fwd inputs
    const TensorList &splats_world,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const std::optional<at::Tensor> camera_ids,  // [nnz]
    const std::optional<at::Tensor> gaussian_ids,  // [nnz]
    const at::Tensor aabb,                       // [..., C, N, 2]
    // grad outputs
    const TensorList &v_splats_screen,
    const TensorList &vr_splats_screen,
    const TensorList &h_splats_screen,
    // returns
    const TensorList &v_splats_world,
    const std::optional<at::Tensor> &v_viewmats,
    const at::Tensor &vr_splats_world,
    const at::Tensor &h_splats_world
) {
    launch_projection_projection_fused_bwd_kernel<MipSplatting, HessianDiagonalOutputMode::Position>(
        splats_world, viewmats, intrins, image_width, image_height, cmt(camera_model), dist_coeffs,
        camera_ids, gaussian_ids, aabb, v_splats_screen, vr_splats_screen, h_splats_screen,
        v_splats_world, v_viewmats, vr_splats_world, h_splats_world);
}



// ================
// Vanilla3DGUT
// ================

/*[AutoHeaderGeneratorExport]*/
void projection_3dgut_backward_tensor(
    // fwd inputs
    const TensorList &splats_world,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const std::optional<at::Tensor> camera_ids,  // [nnz]
    const std::optional<at::Tensor> gaussian_ids,  // [nnz]
    const at::Tensor aabb,                       // [..., C, N, 2]
    // grad outputs
    const TensorList &v_splats_screen,
    // returns
    const TensorList &v_splats_world,
    const std::optional<at::Tensor> &v_viewmats
) {
    launch_projection_projection_fused_bwd_kernel<Vanilla3DGUT, HessianDiagonalOutputMode::None>(
        splats_world, viewmats, intrins, image_width, image_height, cmt(camera_model), dist_coeffs,
        camera_ids, gaussian_ids, aabb, v_splats_screen, {}, {},
        v_splats_world, v_viewmats, std::nullopt, std::nullopt);
}

/*[AutoHeaderGeneratorExport]*/
void projection_3dgut_backward_with_hessian_diagonal_tensor(
    // fwd inputs
    const TensorList &splats_world,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const std::optional<at::Tensor> camera_ids,  // [nnz]
    const std::optional<at::Tensor> gaussian_ids,  // [nnz]
    const at::Tensor aabb,                       // [..., C, N, 2]
    // grad outputs
    const TensorList &v_splats_screen,
    const TensorList &vr_splats_screen,
    const TensorList &h_splats_screen,
    // returns
    const TensorList &v_splats_world,
    const std::optional<at::Tensor> &v_viewmats,
    const TensorList &vr_splats_world,
    const TensorList &h_splats_world
) {
    launch_projection_projection_fused_bwd_kernel<Vanilla3DGUT, HessianDiagonalOutputMode::AllReasonable>(
        splats_world, viewmats, intrins, image_width, image_height, cmt(camera_model), dist_coeffs,
        camera_ids, gaussian_ids, aabb, v_splats_screen, vr_splats_screen, h_splats_screen,
        v_splats_world, v_viewmats, vr_splats_world, h_splats_world);
}

/*[AutoHeaderGeneratorExport]*/
void projection_3dgut_backward_with_position_hessian_diagonal_tensor(
    // fwd inputs
    const TensorList &splats_world,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const std::optional<at::Tensor> camera_ids,  // [nnz]
    const std::optional<at::Tensor> gaussian_ids,  // [nnz]
    const at::Tensor aabb,                       // [..., C, N, 2]
    // grad outputs
    const TensorList &v_splats_screen,
    const TensorList &vr_splats_screen,
    const TensorList &h_splats_screen,
    // returns
    const TensorList &v_splats_world,
    const std::optional<at::Tensor> &v_viewmats,
    const at::Tensor &vr_splats_world,
    const at::Tensor &h_splats_world
) {
    launch_projection_projection_fused_bwd_kernel<Vanilla3DGUT, HessianDiagonalOutputMode::Position>(
        splats_world, viewmats, intrins, image_width, image_height, cmt(camera_model), dist_coeffs,
        camera_ids, gaussian_ids, aabb, v_splats_screen, vr_splats_screen, h_splats_screen,
        v_splats_world, v_viewmats, vr_splats_world, h_splats_world);
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

