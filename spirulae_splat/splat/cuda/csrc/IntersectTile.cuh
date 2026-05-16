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
    at::Tensor  // offsets, [I * n_tiles], int32
> do_intersect_tile_generic(
    at::Tensor aabb,  // [..., N, 4], float32, xyxy in pixels
    at::Tensor depths,  // [..., N], float32
    std::optional<std::tuple<at::optional<at::Tensor>, at::Tensor, at::Tensor>> splats_proj,
    const uint32_t I,
    at::Tensor intrins,
    const uint32_t image_width,
    const uint32_t image_height,
    std::optional<at::Tensor> image_ids
);
