#pragma once

#include <cuda_runtime.h>
#include <cstdint>

#include <ATen/Tensor.h>

#include <gsplat/Common.h>

#include "Primitive3DGS.cuh"
#include "Primitive3DGUT.cuh"
#include "Primitive3DGUT_SV.cuh"
#include "PrimitiveOpaqueTriangle.cuh"
#include "PrimitiveVoxel.cuh"

#include "common.cuh"
#include "types.cuh"


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



std::tuple<
    at::Tensor,  // aabb
    at::Tensor,  // sorting_depths
    at::Tensor,  // radii
    TensorList  // out splats
> projection_3dgs_forward_tensor(
    // inputs
    const TensorList &in_splats,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs
);


std::tuple<
    at::Tensor,  // aabb
    at::Tensor,  // sorting_depths
    at::Tensor,  // radii
    TensorList  // out splats
> projection_mip_forward_tensor(
    // inputs
    const TensorList &in_splats,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs
);


std::tuple<
    at::Tensor,  // aabb
    at::Tensor,  // sorting_depths
    at::Tensor,  // radii
    TensorList  // out splats
> projection_3dgut_forward_tensor(
    // inputs
    const TensorList &in_splats,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs
);


// std::tuple<
//     at::Tensor,  // aabb
//     TensorList  // out splats
// > projection_3dgut_sv_forward_tensor(
//     // inputs
//     const TensorList &in_splats,
//     const at::Tensor viewmats,  // [..., C, 4, 4]
//     const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
//     const uint32_t image_width,
//     const uint32_t image_height,
//     const std::string camera_model,
//     const CameraDistortionCoeffsTensor dist_coeffs
// );


// std::tuple<
//     at::Tensor,  // aabb
//     OpaqueTriangle::Screen::TensorTupleProj  // out splats
// > projection_opaque_triangle_forward_tensor(
//     // inputs
//     const OpaqueTriangle::World::TensorTuple &in_splats,
//     const at::Tensor viewmats,  // [..., C, 4, 4]
//     const at::Tensor intrins,   // [..., C, 4], fx, fy, cx, cy
//     const uint32_t image_width,
//     const uint32_t image_height,
//     const std::string camera_model,
//     const CameraDistortionCoeffsTensor dist_coeffs
// );


// std::tuple<
//     at::Tensor,  // aabb
//     VoxelPrimitive::Screen::TensorTupleProj  // out splats
// > projection_voxel_forward_tensor(
//     // inputs
//     const VoxelPrimitive::World::TensorTuple &in_splats,
//     const at::Tensor viewmats,  // [..., C, 4, 4]
//     const at::Tensor intrins,   // [..., C, 4], fx, fy, cx, cy
//     const uint32_t image_width,
//     const uint32_t image_height,
//     const std::string camera_model,
//     const CameraDistortionCoeffsTensor dist_coeffs
// );
