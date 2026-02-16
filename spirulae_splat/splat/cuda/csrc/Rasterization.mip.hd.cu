#include "RasterizationBwd.cuh"

#include "Primitive3DGS.cuh"


/*[AutoHeaderGeneratorExport]*/
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
    const at::Tensor v_render_alphas, // [..., image_height, image_width, 1]
    std::optional<typename MipSplatting::RenderOutput::TensorTuple> v_distortion_outputs
) {
    if (v_distortion_outputs.has_value())
        return _rasterize_to_pixels_bwd_tensor<MipSplatting, true, true>(
            splats_tuple, backgrounds, masks,
            image_width, image_height, tile_offsets, flatten_ids,
            render_Ts, last_ids, render_outputs, render2_outputs, loss_map,
            v_render_outputs, v_render_alphas, v_distortion_outputs
        );
    return _rasterize_to_pixels_bwd_tensor<MipSplatting, false, true>(
        splats_tuple, backgrounds, masks,
        image_width, image_height, tile_offsets, flatten_ids,
        render_Ts, last_ids, render_outputs, render2_outputs, loss_map,
        v_render_outputs, v_render_alphas, v_distortion_outputs
    );
}
