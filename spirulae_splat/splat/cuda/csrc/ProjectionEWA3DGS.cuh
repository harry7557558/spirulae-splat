#include <cuda_runtime.h>
#include <cstdint>

#include <torch/types.h>

#include <gsplat/Common.h>

#include "Primitive3DGS.cuh"


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



std::tuple<
    at::Tensor,  // radii
    at::Tensor,  // depths
    Vanilla3DGS::Screen::TensorTuple  // out splats
> projection_ewa_3dgs_forward_tensor(
    const bool antialiased,
    // inputs
    const Vanilla3DGS::World::TensorTuple &in_splats,
    const at::Tensor viewmats,             // [..., C, 4, 4]
    const at::Tensor Ks,                   // [..., C, 3, 3]
    const uint32_t image_width,
    const uint32_t image_height,
    const float eps2d,
    const float near_plane,
    const float far_plane,
    const float radius_clip,
    const gsplat::CameraModelType camera_model
);


std::tuple<
    Vanilla3DGS::World::TensorTuple,  // v_splats
    at::Tensor  // v_viewmats
> projection_ewa_3dgs_backward_tensor(
    const bool antialiased,
    // fwd inputs
    const Vanilla3DGS::World::TensorTuple &splats_world_tuple,
    const at::Tensor viewmats,             // [..., C, 4, 4]
    const at::Tensor Ks,                   // [..., C, 3, 3]
    const uint32_t image_width,
    const uint32_t image_height,
    const float eps2d,
    const gsplat::CameraModelType camera_model,
    // fwd outputs
    const at::Tensor radii,                       // [..., C, N, 2]
    // grad outputs
    const at::Tensor v_depths,                      // [..., C, N]
    const Vanilla3DGS::Screen::TensorTuple &v_splats_screen_tuple,
    const bool viewmats_requires_grad
);
