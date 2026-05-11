#pragma once

#include <cuda.h>
#include <cuda_runtime.h>
#include <cstdint>

#include <ATen/Tensor.h>

#include "Primitive3DGS.cuh"
// #include "Primitive3DGUT.cuh"
// #include "Primitive3DGUT_SV.cuh"
// #include "PrimitiveOpaqueTriangle.cuh"
// #include "PrimitiveVoxel.cuh"

#include "common.cuh"


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



std::tuple<
    RenderOutput::TensorTuple,
    at::Tensor,
    at::Tensor,
    std::optional<RenderOutput::TensorTuple>,
    std::optional<RenderOutput::TensorTuple>
> rasterize_to_pixels_3dgs_fwd(
    // Gaussian parameters
    TensorList splats_w,
    TensorList splats_s,
    std::optional<at::Tensor> gaussian_ids,
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    // intersections
    const at::Tensor tile_offsets, // [..., tile_height, tile_width]
    const at::Tensor flatten_ids,   // [n_isects]
    bool output_distortion
);


std::tuple<
    RenderOutput::TensorTuple,
    at::Tensor,
    at::Tensor,
    std::optional<RenderOutput::TensorTuple>,
    std::optional<RenderOutput::TensorTuple>
> rasterize_to_pixels_mip_fwd(
    // Gaussian parameters
    TensorList splats_w,
    TensorList splats_s,
    std::optional<at::Tensor> gaussian_ids,
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    // intersections
    const at::Tensor tile_offsets, // [..., tile_height, tile_width]
    const at::Tensor flatten_ids,   // [n_isects]
    bool output_distortion
);
