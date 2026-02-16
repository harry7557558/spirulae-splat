#include "ProjectionBwd.cuh"

#include "Primitive3DGS.cuh"

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    MipSplatting::World::TensorTuple,  // v_splats
    at::Tensor,  // v_viewmats
    std::variant<at::Tensor, MipSplatting::World::TensorTuple>,  // vr_world_pos or vr_splats
    std::variant<at::Tensor, MipSplatting::World::TensorTuple>  // h_world_pos or h_splats
> projection_mip_backward_with_position_hessian_diagonal_tensor(
    // fwd inputs
    const MipSplatting::World::TensorTuple &splats_world,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const at::Tensor aabb,                       // [..., C, N, 2]
    // grad outputs
    const MipSplatting::Screen::TensorTupleProj &v_splats_screen,
    const MipSplatting::Screen::TensorTupleProj &vr_splats_screen,
    const MipSplatting::Screen::TensorTupleProj &h_splats_screen,
    const bool viewmats_requires_grad
) {
    return _launch_projection_projection_fused_bwd_kernel<MipSplatting, HessianDiagonalOutputMode::Position>(
        splats_world, viewmats, intrins, image_width, image_height, camera_model, dist_coeffs,
        aabb, v_splats_screen, &vr_splats_screen, &h_splats_screen, viewmats_requires_grad);
}
