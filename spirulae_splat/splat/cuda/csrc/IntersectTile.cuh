#pragma once

#include <cstdint>

#include <ATen/Tensor.h>
#include <c10/core/TensorOptions.h>
#include <ATen/Device.h>

#include "Primitive3DGS.cuh"
#include "Primitive3DGUT.cuh"
#include "Primitive3DGUT_SV.cuh"
#include "PrimitiveOpaqueTriangle.cuh"
#include "PrimitiveVoxel.cuh"


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



std::tuple<
    at::Tensor,  // isect_ids, [n_isects], int64
    at::Tensor,  // flatten_ids, [n_isects], int32
    at::Tensor,  // offsets, [I * n_tiles], int32
    at::Tensor  // radii, [N], float32
> do_intersect_tile_generic(
    at::Tensor aabb,  // [..., N, 4], float32, xyxy in pixels
    at::Tensor depths,  // [..., N], float32
    const uint32_t I,
    const uint32_t image_width,
    const uint32_t image_height,
    std::optional<at::Tensor> image_ids
);


std::tuple<
    at::Tensor,  // isect_ids, [n_isects], int64
    at::Tensor,  // flatten_ids, [n_isects], int32
    at::Tensor,  // offsets, [I * n_tiles], int32
    at::Tensor  // radii, [N], float32
> intersect_tile_3dgs_tensor(
    at::Tensor aabb,  // [..., N, 4], float, xyxy in pixels
    at::Tensor depths,  // [..., N], float32
    const uint32_t I,
    const uint32_t image_width,
    const uint32_t image_height,
    typename Vanilla3DGS::Screen::TensorTuple splats,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs
);


std::tuple<
    at::Tensor,  // isect_ids, [n_isects], int64
    at::Tensor,  // flatten_ids, [n_isects], int32
    at::Tensor,  // offsets, [I * n_tiles], int32
    at::Tensor  // radii, [N], float32
> intersect_tile_mip_tensor(
    at::Tensor aabb,  // [..., N, 4], float, xyxy in pixels
    at::Tensor depths,  // [..., N], float32
    const uint32_t I,
    const uint32_t image_width,
    const uint32_t image_height,
    typename MipSplatting::Screen::TensorTuple splats,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs
);


std::tuple<
    at::Tensor,  // isect_ids, [n_isects], int64
    at::Tensor,  // flatten_ids, [n_isects], int32
    at::Tensor,  // offsets, [I * n_tiles], int32
    at::Tensor  // radii, [N], float32
> intersect_tile_3dgut_tensor(
    at::Tensor aabb,  // [..., N, 4], float, xyxy in pixels
    at::Tensor depths,  // [..., N], float32
    const uint32_t I,
    const uint32_t image_width,
    const uint32_t image_height,
    typename Vanilla3DGUT::Screen::TensorTuple splats,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs
);


std::tuple<
    at::Tensor,  // isect_ids, [n_isects], int64
    at::Tensor,  // flatten_ids, [n_isects], int32
    at::Tensor,  // offsets, [I * n_tiles], int32
    at::Tensor  // radii, [N], float32
> intersect_tile_3dgut_sv_tensor(
    at::Tensor aabb,  // [..., N, 4], float, xyxy in pixels
    at::Tensor depths,  // [..., N], float32
    const uint32_t I,
    const uint32_t image_width,
    const uint32_t image_height,
    typename SphericalVoronoi3DGUT_Default::Screen::TensorTuple splats,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs
);


std::tuple<
    at::Tensor,  // isect_ids, [n_isects], int64
    at::Tensor,  // flatten_ids, [n_isects], int32
    at::Tensor,  // offsets, [I * n_tiles], int32
    at::Tensor  // radii, [N], float32
> intersect_tile_opaque_triangle_tensor(
    at::Tensor aabb,  // [..., N, 4], float, xyxy in pixels
    at::Tensor depths,  // [..., N], float32
    const uint32_t I,
    const uint32_t image_width,
    const uint32_t image_height,
    typename OpaqueTriangle::Screen::TensorTuple splats,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs
);


std::tuple<
    at::Tensor,  // isect_ids, [n_isects], int64
    at::Tensor,  // flatten_ids, [n_isects], int32
    at::Tensor,  // offsets, [I * n_tiles], int32
    at::Tensor  // radii, [N], float32
> intersect_tile_voxel_tensor(
    at::Tensor aabb,  // [..., N, 4], float, xyxy in pixels
    at::Tensor depths,  // [..., N], float32
    const uint32_t I,
    const uint32_t image_width,
    const uint32_t image_height,
    typename VoxelPrimitive::Screen::TensorTuple splats,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs
);
