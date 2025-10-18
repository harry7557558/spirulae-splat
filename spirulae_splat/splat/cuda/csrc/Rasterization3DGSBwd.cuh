#include <cuda.h>
#include <cuda_runtime.h>
#include <cstdint>

#include <torch/types.h>

#include "Primitive3DGS.cuh"


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



template <uint32_t CDIM>
void launch_rasterize_to_pixels_3dgs_bwd_kernel(
    // Gaussian parameters
    Vanilla3DGS::Screen::Tensor splats,
    const at::Tensor colors,                    // [..., N, 3] or [nnz, 3]
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
    const at::Tensor v_render_colors, // [..., image_height, image_width, 3]
    const at::Tensor v_render_alphas, // [..., image_height, image_width, 1]
    // outputs
    Vanilla3DGS::Screen::Tensor v_splats,
    at::Tensor v_colors                    // [..., N, 3] or [nnz, 3]
);


std::tuple<
    Vanilla3DGS::Screen::TensorTuple,
    at::Tensor,  // v_colors
    std::optional<at::Tensor>  // absgrad
> rasterize_to_pixels_3dgs_bwd(
    // Gaussian parameters
    Vanilla3DGS::Screen::TensorTuple splats_tuple,
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
