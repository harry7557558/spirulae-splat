#include <cuda.h>
#include <cuda_runtime.h>
#include <cstdint>

#include <torch/types.h>

#include "Primitive3DGS.cuh"
#include "PrimitiveOpaqueTriangle.cuh"


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



template <typename SplatPrimitive>
inline void launch_rasterize_to_pixels_bwd_kernel(
    // Gaussian parameters
    typename SplatPrimitive::Screen::Tensor splats,
    const std::optional<at::Tensor> backgrounds, // [..., 3]
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
    typename SplatPrimitive::RenderOutput::Tensor v_render_outputs,
    const at::Tensor v_render_alphas, // [..., image_height, image_width, 1]
    // outputs
    typename SplatPrimitive::Screen::Tensor v_splats
);


template<typename SplatPrimitive>
inline std::tuple<
    typename SplatPrimitive::Screen::TensorTuple,
    std::optional<at::Tensor>  // absgrad
> rasterize_to_pixels_bwd_tensor(
    // Gaussian parameters
    typename SplatPrimitive::Screen::TensorTuple splats_tuple,
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
    // gradients of outputs
    typename SplatPrimitive::RenderOutput::TensorTuple v_render_outputs,
    const at::Tensor v_render_alphas, // [..., image_height, image_width, 1]
    // options
    bool absgrad
);


std::tuple<
    Vanilla3DGS::Screen::TensorTuple,
    std::optional<at::Tensor>  // absgrad
> rasterize_to_pixels_3dgs_bwd(
    // Gaussian parameters
    Vanilla3DGS::Screen::TensorTuple splats_tuple,
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
    // gradients of outputs
    Vanilla3DGS::RenderOutput::TensorTuple v_render_outputs,
    const at::Tensor v_render_alphas, // [..., image_height, image_width, 1]
    // options
    bool absgrad
);


std::tuple<
    OpaqueTriangle::Screen::TensorTuple,
    at::Tensor,  // v_colors
    std::optional<at::Tensor>  // absgrad
> rasterize_to_pixels_opaque_triangle_bwd(
    // Gaussian parameters
    OpaqueTriangle::Screen::TensorTuple splats_tuple,
    const at::Tensor colors,                    // [..., N, channels] or [nnz, channels]
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
    // gradients of outputs
    const at::Tensor v_render_colors, // [..., image_height, image_width, channels]
    const at::Tensor v_render_alphas, // [..., image_height, image_width, 1]
    // options
    bool absgrad
);
