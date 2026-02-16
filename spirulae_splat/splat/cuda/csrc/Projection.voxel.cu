#include "ProjectionFwd.cuh"
#include "ProjectionBwd.cuh"

#include "PrimitiveVoxel.cuh"

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    at::Tensor,  // aabb
    VoxelPrimitive::Screen::TensorTupleProj  // out splats
> projection_voxel_forward_tensor(
    // inputs
    const VoxelPrimitive::World::TensorTuple &in_splats,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,   // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const float near_plane,
    const float far_plane,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs
) {
    return launch_projection_fused_fwd_kernel<VoxelPrimitive>(
        in_splats, viewmats, intrins, image_width, image_height, near_plane, far_plane, camera_model, dist_coeffs);
}


/*[AutoHeaderGeneratorExport]*/
std::tuple<
    VoxelPrimitive::World::TensorTuple,  // v_splats
    at::Tensor  // v_viewmats
> projection_voxel_backward_tensor(
    // fwd inputs
    const VoxelPrimitive::World::TensorTuple &splats_world,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const at::Tensor aabb,                       // [..., C, N, 2]
    // grad outputs
    const VoxelPrimitive::Screen::TensorTupleProj &v_splats_screen,
    const bool viewmats_requires_grad
) {
    return launch_projection_projection_fused_bwd_kernel<VoxelPrimitive>(
        splats_world, viewmats, intrins, image_width, image_height, camera_model, dist_coeffs,
        aabb, v_splats_screen, viewmats_requires_grad);
}

