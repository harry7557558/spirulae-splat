#include <cuda.h>
#include <cuda_runtime.h>
#include <cstdint>

#include <torch/types.h>

#include "Primitive3DGUT.cuh"
#include "PrimitiveOpaqueTriangle.cuh"
#include "PrimitiveVoxel.cuh"

#include "types.cuh"


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



std::tuple<
    Vanilla3DGUT::Screen::TensorTuple,
    std::optional<at::Tensor>  // v_viewmats
> rasterize_to_pixels_3dgut_bwd(
    // Gaussian parameters
    Vanilla3DGUT::Screen::TensorTuple splats_tuple,
    const at::Tensor viewmats,             // [..., C, 4, 4]
    const at::Tensor Ks,                   // [..., C, 3, 3]
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    const std::optional<at::Tensor> backgrounds, // [..., channels]
    const std::optional<at::Tensor> masks,       // [..., tile_height, tile_width]
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_size,
    // intersections
    const at::Tensor tile_offsets, // [..., tile_height, tile_width]
    const at::Tensor flatten_ids,  // [n_isects]
    // forward outputs
    const at::Tensor render_Ts, // [..., image_height, image_width, 1]
    const at::Tensor last_ids,      // [..., image_height, image_width]
    std::optional<typename Vanilla3DGUT::RenderOutput::TensorTuple> render_outputs,
    std::optional<typename Vanilla3DGUT::RenderOutput::TensorTuple> render2_outputs,
    // gradients of outputs
    Vanilla3DGUT::RenderOutput::TensorTuple v_render_outputs,
    const at::Tensor v_render_alphas, // [..., image_height, image_width, 1]
    std::optional<typename Vanilla3DGUT::RenderOutput::TensorTuple> v_distortion_outputs
);


std::tuple<
    OpaqueTriangle::Screen::TensorTuple,
    std::optional<at::Tensor>  // v_viewmats
> rasterize_to_pixels_opaque_triangle_eval3d_bwd(
    // Gaussian parameters
    OpaqueTriangle::Screen::TensorTuple splats_tuple,
    const at::Tensor viewmats,             // [..., C, 4, 4]
    const at::Tensor Ks,                   // [..., C, 3, 3]
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    const std::optional<at::Tensor> backgrounds, // [..., channels]
    const std::optional<at::Tensor> masks,       // [..., tile_height, tile_width]
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_size,
    // intersections
    const at::Tensor tile_offsets, // [..., tile_height, tile_width]
    const at::Tensor flatten_ids,  // [n_isects]
    // forward outputs
    const at::Tensor render_Ts, // [..., image_height, image_width, 1]
    const at::Tensor last_ids,      // [..., image_height, image_width]
    std::optional<typename OpaqueTriangle::RenderOutput::TensorTuple> render_outputs,
    std::optional<typename OpaqueTriangle::RenderOutput::TensorTuple> render2_outputs,
    // gradients of outputs
    OpaqueTriangle::RenderOutput::TensorTuple v_render_outputs,
    const at::Tensor v_render_alphas, // [..., image_height, image_width, 1]
    std::optional<typename OpaqueTriangle::RenderOutput::TensorTuple> v_distortion_outputs
);


std::tuple<
    VoxelPrimitive::Screen::TensorTuple,
    std::optional<at::Tensor>  // v_viewmats
> rasterize_to_pixels_voxel_eval3d_bwd(
    // Gaussian parameters
    VoxelPrimitive::Screen::TensorTuple splats_tuple,
    const at::Tensor viewmats,             // [..., C, 4, 4]
    const at::Tensor Ks,                   // [..., C, 3, 3]
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    const std::optional<at::Tensor> backgrounds, // [..., channels]
    const std::optional<at::Tensor> masks,       // [..., tile_height, tile_width]
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_size,
    // intersections
    const at::Tensor tile_offsets, // [..., tile_height, tile_width]
    const at::Tensor flatten_ids,  // [n_isects]
    // forward outputs
    const at::Tensor render_Ts, // [..., image_height, image_width, 1]
    const at::Tensor last_ids,      // [..., image_height, image_width]
    std::optional<typename VoxelPrimitive::RenderOutput::TensorTuple> render_outputs,
    std::optional<typename VoxelPrimitive::RenderOutput::TensorTuple> render2_outputs,
    // gradients of outputs
    VoxelPrimitive::RenderOutput::TensorTuple v_render_outputs,
    const at::Tensor v_render_alphas, // [..., image_height, image_width, 1]
    std::optional<typename VoxelPrimitive::RenderOutput::TensorTuple> v_distortion_outputs
);
