#include <cuda_runtime.h>
#include <cstdint>

#include <torch/types.h>

#include <gsplat/Common.h>

#include "Primitive3DGS.cuh"
#include "PrimitiveOpaqueTriangle.cuh"


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



std::tuple<
    at::Tensor,  // camera_ids
    at::Tensor,  // gaussian_ids
    at::Tensor,  // aabb
    Vanilla3DGS::Screen::TensorTupleProj  // out splats
> projection_3dgs_hetero_forward_tensor(
    // inputs
    const Vanilla3DGS::World::TensorTuple &in_splats_tensor,
    const at::Tensor viewmats,             // [..., C, 4, 4]
    const at::Tensor Ks,                   // [..., C, 3, 3]
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_width,
    const uint32_t tile_height,
    const float near_plane,
    const float far_plane,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    const at::Tensor intersection_count_map,  // [C+1]
    const at::Tensor intersection_splat_id  // [nnz]
);


std::tuple<
    at::Tensor,  // camera_ids
    at::Tensor,  // gaussian_ids
    at::Tensor,  // aabb
    OpaqueTriangle::Screen::TensorTupleProj  // out splats
> projection_opaque_triangle_hetero_forward_tensor(
    // inputs
    const OpaqueTriangle::World::TensorTuple &in_splats_tensor,
    const at::Tensor viewmats,             // [..., C, 4, 4]
    const at::Tensor Ks,                   // [..., C, 3, 3]
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_width,
    const uint32_t tile_height,
    const float near_plane,
    const float far_plane,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    const at::Tensor intersection_count_map,  // [C+1]
    const at::Tensor intersection_splat_id  // [nnz]
);


template<typename SplatPrimitive, gsplat::CameraModelType camera_model>
__global__ void projection_3dgs_hetero_backward_kernel(
    // fwd inputs
    const long C,
    const long N,
    const uint32_t nnz,
    const typename SplatPrimitive::World::Buffer splats_world,
    const float *__restrict__ viewmats, // [C, 4, 4]
    const float *__restrict__ Ks,       // [C, 3, 3]
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_width,
    const uint32_t tile_height,
    // fwd outputs
    const int64_t *__restrict__ camera_ids,     // [nnz]
    const int64_t *__restrict__ gaussian_ids,   // [nnz]
    const int4 *__restrict__ aabbs,          // [B, C, N, 4]
    // grad outputs
    typename SplatPrimitive::Screen::Buffer v_splats_proj,
    const bool sparse_grad, // whether the outputs are in COO format [nnz, ...]
    // grad inputs
    typename SplatPrimitive::World::Buffer v_splats_world,
    float *__restrict__ v_viewmats // [C, 4, 4]
);


std::tuple<
    Vanilla3DGS::World::TensorTuple,  // v_splats
    at::Tensor  // v_viewmats
> projection_3dgs_hetero_backward_tensor(
    // fwd inputs
    const Vanilla3DGS::World::TensorTuple &splats_world_tuple,
    const at::Tensor viewmats, // [..., C, 4, 4]
    const at::Tensor Ks, // [..., C, 3, 3]
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
    const Vanilla3DGS::Screen::TensorTupleProj &v_splats_proj_tuple,
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
    const at::Tensor Ks, // [..., C, 3, 3]
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
    const OpaqueTriangle::Screen::TensorTupleProj &v_splats_proj_tuple,
    const bool viewmats_requires_grad,
    const bool sparse_grad
);
