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



Vanilla3DGS::Screen::TensorTuple rasterize_to_pixels_3dgs_bwd(
    // Gaussian parameters
    Vanilla3DGS::Screen::TensorTuple splats_tuple,
    const std::optional<at::Tensor> backgrounds, // [..., channels]
    const std::optional<at::Tensor> masks,       // [..., tile_height, tile_width]
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    // intersections
    const at::Tensor tile_offsets, // [..., tile_height, tile_width]
    const at::Tensor flatten_ids,  // [n_isects]
    // forward outputs
    const at::Tensor render_Ts, // [..., image_height, image_width, 1]
    const at::Tensor last_ids,      // [..., image_height, image_width]
    // gradients of outputs
    Vanilla3DGS::RenderOutput::TensorTuple v_render_outputs,
    const at::Tensor v_render_Ts // [..., image_height, image_width, 1]
);


std::tuple<
    Vanilla3DGS::Screen::TensorTuple,
    std::optional<Vanilla3DGS::Screen::TensorTuple>,  // jacobian residual product
    std::optional<Vanilla3DGS::Screen::TensorTuple>  // hessian diagonal
> rasterize_to_pixels_3dgs_bwd_with_hessian_diagonal(
    // Gaussian parameters
    Vanilla3DGS::Screen::TensorTuple splats_tuple,
    const std::optional<at::Tensor> backgrounds, // [..., channels]
    const std::optional<at::Tensor> masks,       // [..., tile_height, tile_width]
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    // intersections
    const at::Tensor tile_offsets, // [..., tile_height, tile_width]
    const at::Tensor flatten_ids,  // [n_isects]
    // forward outputs
    const at::Tensor render_Ts, // [..., image_height, image_width, 1]
    const at::Tensor last_ids,      // [..., image_height, image_width]
    std::optional<typename Vanilla3DGS::RenderOutput::TensorTuple> render_outputs,
    std::optional<typename Vanilla3DGS::RenderOutput::TensorTuple> render2_outputs,
    std::optional<at::Tensor> loss_map,  // [..., image_height, image_width, 1]
    // gradients of outputs
    Vanilla3DGS::RenderOutput::TensorTuple v_render_outputs,
    const at::Tensor v_render_Ts, // [..., image_height, image_width, 1]
    std::optional<typename Vanilla3DGS::RenderOutput::TensorTuple> v_distortion_outputs
);


MipSplatting::Screen::TensorTuple rasterize_to_pixels_mip_bwd(
    // Gaussian parameters
    MipSplatting::Screen::TensorTuple splats_tuple,
    const std::optional<at::Tensor> backgrounds, // [..., channels]
    const std::optional<at::Tensor> masks,       // [..., tile_height, tile_width]
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    // intersections
    const at::Tensor tile_offsets, // [..., tile_height, tile_width]
    const at::Tensor flatten_ids,  // [n_isects]
    // forward outputs
    const at::Tensor render_Ts, // [..., image_height, image_width, 1]
    const at::Tensor last_ids,      // [..., image_height, image_width]
    // gradients of outputs
    MipSplatting::RenderOutput::TensorTuple v_render_outputs,
    const at::Tensor v_render_Ts // [..., image_height, image_width, 1]
);


std::tuple<
    MipSplatting::Screen::TensorTuple,
    std::optional<MipSplatting::Screen::TensorTuple>,  // jacobian residual product
    std::optional<MipSplatting::Screen::TensorTuple>  // hessian diagonal
> rasterize_to_pixels_mip_bwd_with_hessian_diagonal(
    // Gaussian parameters
    MipSplatting::Screen::TensorTuple splats_tuple,
    const std::optional<at::Tensor> backgrounds, // [..., channels]
    const std::optional<at::Tensor> masks,       // [..., tile_height, tile_width]
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    // intersections
    const at::Tensor tile_offsets, // [..., tile_height, tile_width]
    const at::Tensor flatten_ids,  // [n_isects]
    // forward outputs
    const at::Tensor render_Ts, // [..., image_height, image_width, 1]
    const at::Tensor last_ids,      // [..., image_height, image_width]
    std::optional<typename MipSplatting::RenderOutput::TensorTuple> render_outputs,
    std::optional<typename MipSplatting::RenderOutput::TensorTuple> render2_outputs,
    std::optional<at::Tensor> loss_map,  // [..., image_height, image_width, 1]
    // gradients of outputs
    MipSplatting::RenderOutput::TensorTuple v_render_outputs,
    const at::Tensor v_render_Ts, // [..., image_height, image_width, 1]
    std::optional<typename MipSplatting::RenderOutput::TensorTuple> v_distortion_outputs
);
