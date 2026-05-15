#include <cuda_runtime.h>
#include <cstdint>

#include <ATen/Tensor.h>

#include <gsplat/Common.h>

#include "types.cuh"

#include "Primitive3DGS.cuh"
#include "Primitive3DGUT.cuh"
#include "Primitive3DGUT_SV.cuh"
#include "PrimitiveOpaqueTriangle.cuh"
#include "PrimitiveVoxel.cuh"


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



std::tuple<
    at::Tensor,  // camera_ids
    at::Tensor,  // gaussian_ids
    at::Tensor,  // aabb
    at::Tensor,  // sorting depths
    TensorList  // out splats
> projection_3dgs_hetero_forward_tensor(
    // inputs
    const TensorList &in_splats_tensor,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_width,
    const uint32_t tile_height,
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    const at::Tensor intersection_count_map,  // [C+1]
    const at::Tensor intersection_splat_id  // [nnz]
);


std::tuple<
    TensorList,  // v_splats
    at::Tensor  // v_viewmats
> projection_3dgs_hetero_backward_tensor(
    // fwd inputs
    const TensorList &splats_world_tuple,
    const at::Tensor viewmats, // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_width,
    const uint32_t tile_height,
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const at::Tensor camera_ids, // [nnz]
    const at::Tensor gaussian_ids, // [nnz]
    const at::Tensor aabb,  // [nnz, 4]
    // grad outputs
    const TensorList &v_splats_proj_tuple,
    const bool viewmats_requires_grad
);


std::tuple<
    at::Tensor,  // camera_ids
    at::Tensor,  // gaussian_ids
    at::Tensor,  // aabb
    at::Tensor,  // sorting_depths
    TensorList  // out splats
> projection_mip_hetero_forward_tensor(
    // inputs
    const TensorList &in_splats_tensor,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_width,
    const uint32_t tile_height,
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    const at::Tensor intersection_count_map,  // [C+1]
    const at::Tensor intersection_splat_id  // [nnz]
);


std::tuple<
    TensorList,  // v_splats
    at::Tensor  // v_viewmats
> projection_mip_hetero_backward_tensor(
    // fwd inputs
    const TensorList &splats_world_tuple,
    const at::Tensor viewmats, // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_width,
    const uint32_t tile_height,
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const at::Tensor camera_ids, // [nnz]
    const at::Tensor gaussian_ids, // [nnz]
    const at::Tensor aabb,  // [nnz, 4]
    // grad outputs
    const TensorList &v_splats_proj_tuple,
    const bool viewmats_requires_grad
);


std::tuple<
    at::Tensor,  // camera_ids
    at::Tensor,  // gaussian_ids
    at::Tensor,  // aabb
    at::Tensor,  // sorting depths
    TensorList  // out splats
> projection_3dgut_hetero_forward_tensor(
    // inputs
    const TensorList &in_splats_tensor,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_width,
    const uint32_t tile_height,
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    const at::Tensor intersection_count_map,  // [C+1]
    const at::Tensor intersection_splat_id  // [nnz]
);


std::tuple<
    TensorList,  // v_splats
    at::Tensor  // v_viewmats
> projection_3dgut_hetero_backward_tensor(
    // fwd inputs
    const TensorList &splats_world_tuple,
    const at::Tensor viewmats, // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_width,
    const uint32_t tile_height,
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const at::Tensor camera_ids, // [nnz]
    const at::Tensor gaussian_ids, // [nnz]
    const at::Tensor aabb,  // [nnz, 4]
    // grad outputs
    const TensorList &v_splats_proj_tuple,
    const bool viewmats_requires_grad
);


// std::tuple<
//     at::Tensor,  // camera_ids
//     at::Tensor,  // gaussian_ids
//     at::Tensor,  // aabb
//     OpaqueTriangle::Screen::TensorTuple  // out splats
// > projection_opaque_triangle_hetero_forward_tensor(
//     // inputs
//     const OpaqueTriangle::World::TensorTuple &in_splats_tensor,
//     const at::Tensor viewmats,  // [..., C, 4, 4]
//     const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
//     const uint32_t image_width,
//     const uint32_t image_height,
//     const uint32_t tile_width,
//     const uint32_t tile_height,
//     const std::string camera_model,
//     const CameraDistortionCoeffsTensor dist_coeffs,
//     const at::Tensor intersection_count_map,  // [C+1]
//     const at::Tensor intersection_splat_id  // [nnz]
// );


// std::tuple<
//     OpaqueTriangle::World::TensorTuple,  // v_splats
//     at::Tensor  // v_viewmats
// > projection_opaque_triangle_hetero_backward_tensor(
//     // fwd inputs
//     const OpaqueTriangle::World::TensorTuple &splats_world_tuple,
//     const at::Tensor viewmats, // [..., C, 4, 4]
//     const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
//     const uint32_t image_width,
//     const uint32_t image_height,
//     const uint32_t tile_width,
//     const uint32_t tile_height,
//     const std::string camera_model,
//     const CameraDistortionCoeffsTensor dist_coeffs,
//     // fwd outputs
//     const at::Tensor camera_ids, // [nnz]
//     const at::Tensor gaussian_ids, // [nnz]
//     const at::Tensor aabb,  // [nnz, 4]
//     // grad outputs
//     const OpaqueTriangle::Screen::TensorTuple &v_splats_proj_tuple,
//     const bool viewmats_requires_grad
// );
