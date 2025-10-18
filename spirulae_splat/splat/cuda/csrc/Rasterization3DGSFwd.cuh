#include <cuda.h>
#include <cuda_runtime.h>
#include <cstdint>

#include <torch/types.h>

#include "Primitive3DGS.cuh"


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



template <uint32_t CDIM>
void launch_rasterize_to_pixels_3dgs_fwd_kernel(
    // Gaussian parameters
    Vanilla3DGS::Screen::Tensor splats,
    const at::Tensor colors,    // [..., N, channels] or [nnz, channels]
    const std::optional<at::Tensor> backgrounds, // [..., channels]
    const std::optional<at::Tensor> masks,       // [..., tile_height, tile_width]
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_size,
    // intersections
    const at::Tensor tile_offsets, // [..., tile_height, tile_width]
    const at::Tensor flatten_ids,  // [n_isects]
    // outputs
    at::Tensor renders, // [..., image_height, image_width, channels]
    at::Tensor transmittances,  // [..., image_height, image_width]
    at::Tensor last_ids // [..., image_height, image_width]
);


std::tuple<at::Tensor, at::Tensor, at::Tensor> rasterize_to_pixels_3dgs_fwd(
    // Gaussian parameters
    Vanilla3DGS::Screen::TensorTuple splats_tuple,
    const at::Tensor colors,    // [..., N, channels] or [nnz, channels]
    const std::optional<at::Tensor> backgrounds, // [..., channels]
    const std::optional<at::Tensor> masks,       // [..., tile_height, tile_width]
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_size,
    // intersections
    const at::Tensor tile_offsets, // [..., tile_height, tile_width]
    const at::Tensor flatten_ids   // [n_isects]
);
