#include "RasterizationFwd.cuh"
#include "RasterizationBwd.cuh"


/*[AutoHeaderGeneratorExport]*/
std::tuple<Vanilla3DGS::RenderOutput::TensorTuple, at::Tensor, at::Tensor>
rasterize_to_pixels_3dgs_fwd(
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
    const at::Tensor flatten_ids   // [n_isects]
) {
    return rasterize_to_pixels_fwd_tensor<Vanilla3DGS>(
        splats_tuple, backgrounds, masks,
        image_width, image_height, tile_size,
        tile_offsets, flatten_ids
    );
}

/*[AutoHeaderGeneratorExport]*/
Vanilla3DGS::Screen::TensorTuple rasterize_to_pixels_3dgs_bwd(
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
    const at::Tensor v_render_alphas // [..., image_height, image_width, 1]
) {
    return rasterize_to_pixels_bwd_tensor<Vanilla3DGS>(
        splats_tuple, backgrounds, masks,
        image_width, image_height, tile_size, tile_offsets, flatten_ids,
        render_Ts, last_ids, v_render_outputs, v_render_alphas
    );
}
