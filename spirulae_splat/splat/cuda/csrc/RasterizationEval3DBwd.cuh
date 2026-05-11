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
    TensorList, TensorList,  // gradient
    std::optional<at::Tensor>,  // v_viewmats
    std::optional<at::Tensor>  // accum_weight
> rasterize_to_pixels_3dgut_bwd(
    // Gaussian parameters
    int64_t num_splats,
    TensorList splats_w,
    TensorList splats_s,
    std::optional<at::Tensor> gaussian_ids,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    // intersections
    const at::Tensor tile_offsets, // [..., tile_height, tile_width]
    const at::Tensor flatten_ids,  // [n_isects]
    // forward outputs
    const at::Tensor render_Ts, // [..., image_height, image_width, 1]
    const at::Tensor last_ids,      // [..., image_height, image_width]
    std::optional<RenderOutput::TensorTuple> render_outputs,
    std::optional<RenderOutput::TensorTuple> render2_outputs,
    std::optional<at::Tensor> loss_map,  // [..., image_height, image_width, 1]
    std::optional<at::Tensor> accum_weight_map,  // [..., image_height, image_width, 1]
    // gradients of outputs
    RenderOutput::TensorTuple v_render_outputs,
    const at::Tensor v_render_Ts, // [..., image_height, image_width, 1]
    std::optional<RenderOutput::TensorTuple> v_distortion_outputs,
    bool need_viewmat_grad
);


std::tuple<
    TensorList, TensorList,  // gradient
    std::optional<at::Tensor>,  // v_viewmats
    std::optional<TensorList>, std::optional<TensorList>,  // jacobian residual product
    std::optional<TensorList>, std::optional<TensorList>,  // hessian diagonal
    std::optional<at::Tensor>  // accum_weight
> rasterize_to_pixels_3dgut_bwd_with_hessian_diagonal(
    // Gaussian parameters
    int64_t num_splats,
    TensorList splats_w,
    TensorList splats_s,
    std::optional<at::Tensor> gaussian_ids,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    // intersections
    const at::Tensor tile_offsets, // [..., tile_height, tile_width]
    const at::Tensor flatten_ids,  // [n_isects]
    // forward outputs
    const at::Tensor render_Ts, // [..., image_height, image_width, 1]
    const at::Tensor last_ids,      // [..., image_height, image_width]
    std::optional<RenderOutput::TensorTuple> render_outputs,
    std::optional<RenderOutput::TensorTuple> render2_outputs,
    std::optional<at::Tensor> loss_map,  // [..., image_height, image_width, 1]
    std::optional<at::Tensor> accum_weight_map,  // [..., image_height, image_width, 1]
    // gradients of outputs
    RenderOutput::TensorTuple v_render_outputs,
    const at::Tensor v_render_Ts, // [..., image_height, image_width, 1]
    std::optional<RenderOutput::TensorTuple> v_distortion_outputs,
    bool need_viewmat_grad
);


// std::tuple<
//     VoxelPrimitive::Screen::TensorTuple,
//     std::optional<at::Tensor>,  // accum_weight
//     std::optional<at::Tensor>  // v_viewmats
// > rasterize_to_pixels_voxel_eval3d_bwd(
//     // Gaussian parameters
//     VoxelPrimitive::Screen::TensorTuple splats_tuple,
//     std::optional<at::Tensor> gaussian_ids,
//     const at::Tensor viewmats,  // [..., C, 4, 4]
//     const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
//     const std::string camera_model,
//     const CameraDistortionCoeffsTensor dist_coeffs,
//     // image size
//     const uint32_t image_width,
//     const uint32_t image_height,
//     // intersections
//     const at::Tensor tile_offsets, // [..., tile_height, tile_width]
//     const at::Tensor flatten_ids,  // [n_isects]
//     // forward outputs
//     const at::Tensor render_Ts, // [..., image_height, image_width, 1]
//     const at::Tensor last_ids,      // [..., image_height, image_width]
//     std::optional<RenderOutput::TensorTuple> render_outputs,
//     std::optional<RenderOutput::TensorTuple> render2_outputs,
//     std::optional<at::Tensor> loss_map,  // [..., image_height, image_width, 1]
//     std::optional<at::Tensor> accum_weight_map,  // [..., image_height, image_width, 1]
//     // gradients of outputs
//     RenderOutput::TensorTuple v_render_outputs,
//     const at::Tensor v_render_Ts, // [..., image_height, image_width, 1]
//     std::optional<RenderOutput::TensorTuple> v_distortion_outputs,
//     bool need_viewmat_grad
// );
