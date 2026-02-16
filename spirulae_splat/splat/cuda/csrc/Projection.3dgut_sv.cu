#include "ProjectionFwd.cuh"
#include "ProjectionBwd.cuh"

#include "Primitive3DGUT_SV.cuh"


/*[AutoHeaderGeneratorExport]*/
std::tuple<
    at::Tensor,  // aabb
    SphericalVoronoi3DGUT_Default::Screen::TensorTupleProj  // out splats
> projection_3dgut_sv_forward_tensor(
    // inputs
    const SphericalVoronoi3DGUT_Default::World::TensorTuple &in_splats,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const float near_plane,
    const float far_plane,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs
) {
    int num_sv = std::get<5>(in_splats).size(-2);
    #define _CASE(n) \
        if (num_sv == n) return launch_projection_fused_fwd_kernel<SphericalVoronoi3DGUT<n>>( \
            in_splats, viewmats, intrins, image_width, image_height, near_plane, far_plane, camera_model, dist_coeffs); \
    _CASE(2) _CASE(3) _CASE(4) _CASE(5) _CASE(6) _CASE(7) _CASE(8)
    #undef _CASE
    throw std::invalid_argument("Unsupported num_sv");
}

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    SphericalVoronoi3DGUT_Default::World::TensorTuple,  // v_splats
    at::Tensor  // v_viewmats
> projection_3dgut_sv_backward_tensor(
    // fwd inputs
    const SphericalVoronoi3DGUT_Default::World::TensorTuple &splats_world,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const at::Tensor aabb,                       // [..., C, N, 2]
    // grad outputs
    const SphericalVoronoi3DGUT_Default::Screen::TensorTupleProj &v_splats_screen,
    const bool viewmats_requires_grad
) {
    int num_sv = std::get<5>(splats_world).size(-2);
    #define _CASE(n) \
        if (num_sv == n) return launch_projection_projection_fused_bwd_kernel<SphericalVoronoi3DGUT<n>>( \
            splats_world, viewmats, intrins, image_width, image_height, camera_model, dist_coeffs, \
            aabb, v_splats_screen, viewmats_requires_grad);
    _CASE(2) _CASE(3) _CASE(4) _CASE(5) _CASE(6) _CASE(7) _CASE(8)
    #undef _CASE
    throw std::invalid_argument("Unsupported num_sv");
}

