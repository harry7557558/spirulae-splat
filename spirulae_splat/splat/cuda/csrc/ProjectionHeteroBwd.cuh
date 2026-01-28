#include <cuda_runtime.h>
#include <cstdint>

#include <ATen/Tensor.h>

#include <gsplat/Common.h>

#include "Primitive3DGS.cuh"
#include "Primitive3DGUT.cuh"
#include "PrimitiveOpaqueTriangle.cuh"


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



std::tuple<
    Vanilla3DGS::World::TensorTuple,  // v_splats
    at::Tensor  // v_viewmats
> projection_3dgs_hetero_backward_tensor(
    // fwd inputs
    const Vanilla3DGS::World::TensorTuple &splats_world_tuple,
    const at::Tensor viewmats, // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_width,
    const uint32_t tile_height,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const at::Tensor camera_ids, // [nnz]
    const at::Tensor gaussian_ids, // [nnz]
    const at::Tensor aabb,  // [nnz, 4]
    // grad outputs
    const Vanilla3DGS::Screen::TensorTuple &v_splats_proj_tuple,
    const bool viewmats_requires_grad,
    const bool sparse_grad
);


std::tuple<
    MipSplatting::World::TensorTuple,  // v_splats
    at::Tensor  // v_viewmats
> projection_mip_hetero_backward_tensor(
    // fwd inputs
    const MipSplatting::World::TensorTuple &splats_world_tuple,
    const at::Tensor viewmats, // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_width,
    const uint32_t tile_height,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const at::Tensor camera_ids, // [nnz]
    const at::Tensor gaussian_ids, // [nnz]
    const at::Tensor aabb,  // [nnz, 4]
    // grad outputs
    const MipSplatting::Screen::TensorTuple &v_splats_proj_tuple,
    const bool viewmats_requires_grad,
    const bool sparse_grad
);


std::tuple<
    Vanilla3DGUT::World::TensorTuple,  // v_splats
    at::Tensor  // v_viewmats
> projection_3dgut_hetero_backward_tensor(
    // fwd inputs
    const Vanilla3DGUT::World::TensorTuple &splats_world_tuple,
    const at::Tensor viewmats, // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_width,
    const uint32_t tile_height,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const at::Tensor camera_ids, // [nnz]
    const at::Tensor gaussian_ids, // [nnz]
    const at::Tensor aabb,  // [nnz, 4]
    // grad outputs
    const Vanilla3DGUT::Screen::TensorTuple &v_splats_proj_tuple,
    const bool viewmats_requires_grad,
    const bool sparse_grad
);


std::tuple<
    OpaqueTriangle::World::TensorTuple,  // v_splats
    at::Tensor  // v_viewmats
> projection_opaque_triangle_hetero_backward_tensor(
    // fwd inputs
    const OpaqueTriangle::World::TensorTuple &splats_world_tuple,
    const at::Tensor viewmats, // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_width,
    const uint32_t tile_height,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const at::Tensor camera_ids, // [nnz]
    const at::Tensor gaussian_ids, // [nnz]
    const at::Tensor aabb,  // [nnz, 4]
    // grad outputs
    const OpaqueTriangle::Screen::TensorTuple &v_splats_proj_tuple,
    const bool viewmats_requires_grad,
    const bool sparse_grad
);
