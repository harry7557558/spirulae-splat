#pragma once

#include <cuda.h>
#include <cuda_runtime.h>
#include <cstdint>

#include <ATen/Tensor.h>

// #include "Primitive3DGS.cuh"
#include "Primitive3DGUT.cuh"
#include "Primitive3DGUT_SV.cuh"
// #include "PrimitiveOpaqueTriangle.cuh"
#include "PrimitiveVoxel.cuh"

#include "types.cuh"
#include "common.cuh"


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



std::tuple<
    Vanilla3DGUT::RenderOutput::TensorTuple,
    at::Tensor,
    at::Tensor,
    std::optional<Vanilla3DGUT::RenderOutput::TensorTuple>,
    std::optional<Vanilla3DGUT::RenderOutput::TensorTuple>,
    std::optional<at::Tensor>
> rasterize_to_pixels_3dgut_fwd(
    // Gaussian parameters
    Vanilla3DGUT::Screen::TensorTuple splats_tuple,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    const std::optional<at::Tensor> backgrounds, // [..., channels]
    const std::optional<at::Tensor> masks,       // [..., tile_height, tile_width]
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    // intersections
    const at::Tensor tile_offsets, // [..., tile_height, tile_width]
    const at::Tensor flatten_ids,   // [n_isects]
    bool output_distortion
);


std::tuple<
    SphericalVoronoi3DGUT_Default::RenderOutput::TensorTuple,
    at::Tensor,
    at::Tensor,
    std::optional<SphericalVoronoi3DGUT_Default::RenderOutput::TensorTuple>,
    std::optional<SphericalVoronoi3DGUT_Default::RenderOutput::TensorTuple>,
    std::optional<at::Tensor>
> rasterize_to_pixels_3dgut_sv_fwd(
    // Gaussian parameters
    SphericalVoronoi3DGUT_Default::Screen::TensorTuple splats_tuple,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    const std::optional<at::Tensor> backgrounds, // [..., channels]
    const std::optional<at::Tensor> masks,       // [..., tile_height, tile_width]
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    // intersections
    const at::Tensor tile_offsets, // [..., tile_height, tile_width]
    const at::Tensor flatten_ids,   // [n_isects]
    bool output_distortion
);


std::tuple<
    VoxelPrimitive::RenderOutput::TensorTuple,
    at::Tensor,
    at::Tensor,
    std::optional<VoxelPrimitive::RenderOutput::TensorTuple>,
    std::optional<VoxelPrimitive::RenderOutput::TensorTuple>,
    std::optional<at::Tensor>
> rasterize_to_pixels_voxel_eval3d_fwd(
    // Gaussian parameters
    VoxelPrimitive::Screen::TensorTuple splats_tuple,
    const at::Tensor viewmats,      // [..., C, 4, 4]
    const at::Tensor intrins,       // [..., C, 4], fx, fy, cx, cy
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    const std::optional<at::Tensor> backgrounds, // [..., channels]
    const std::optional<at::Tensor> masks,       // [..., tile_height, tile_width]
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    // intersections
    const at::Tensor tile_offsets, // [..., tile_height, tile_width]
    const at::Tensor flatten_ids   // [n_isects]
);
