#include "RasterizationFwd.cuh"
#include "RasterizationBwd.cuh"

#include "Primitive3DGS.cuh"


/*[AutoHeaderGeneratorExport]*/
std::tuple<MipSplatting::RenderOutput::TensorTuple, at::Tensor, at::Tensor>
rasterize_to_pixels_mip_fwd(
    // Gaussian parameters
    MipSplatting::Screen::TensorTuple splats_tuple,
    const std::optional<at::Tensor> backgrounds, // [..., channels]
    const std::optional<at::Tensor> masks,       // [..., tile_height, tile_width]
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    // intersections
    const at::Tensor tile_offsets, // [..., tile_height, tile_width]
    const at::Tensor flatten_ids   // [n_isects]
) {
    return rasterize_to_pixels_fwd_tensor<MipSplatting>(
        splats_tuple, backgrounds, masks,
        image_width, image_height,
        tile_offsets, flatten_ids
    );
}

/*[AutoHeaderGeneratorExport]*/
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
    const at::Tensor v_render_alphas // [..., image_height, image_width, 1]
) {
    return rasterize_to_pixels_bwd_tensor<MipSplatting>(
        splats_tuple, backgrounds, masks,
        image_width, image_height, tile_offsets, flatten_ids,
        render_Ts, last_ids, v_render_outputs, v_render_alphas
    );
}
