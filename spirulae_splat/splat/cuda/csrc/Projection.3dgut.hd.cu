#include "ProjectionFwd.cuh"
#include "ProjectionBwd.cuh"


/*[AutoHeaderGeneratorExport]*/
std::tuple<
    Vanilla3DGUT::World::TensorTuple,  // v_splats
    at::Tensor,  // v_viewmats
    at::Tensor  // h_world_pos
> projection_3dgut_backward_with_position_hessian_diagonal_tensor(
    // fwd inputs
    const Vanilla3DGUT::World::TensorTuple &splats_world,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const at::Tensor aabb,                       // [..., C, N, 2]
    // grad outputs
    const Vanilla3DGUT::Screen::TensorTupleProj &v_splats_screen,
    const Vanilla3DGUT::Screen::TensorTupleProj &h_splats_screen,
    const bool viewmats_requires_grad
) {
    return _launch_projection_projection_fused_bwd_kernel<Vanilla3DGUT, true>(
        splats_world, viewmats, intrins, image_width, image_height, camera_model, dist_coeffs,
        aabb, v_splats_screen, &h_splats_screen, viewmats_requires_grad);
}
