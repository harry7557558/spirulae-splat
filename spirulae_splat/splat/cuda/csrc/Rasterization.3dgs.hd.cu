#include "RasterizationFwd.cuh"
#include "RasterizationBwd.cuh"


/*[AutoHeaderGeneratorExport]*/
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
    const at::Tensor v_render_alphas, // [..., image_height, image_width, 1]
    std::optional<typename Vanilla3DGS::RenderOutput::TensorTuple> v_distortion_outputs
) {
    if (v_distortion_outputs.has_value())
        return _rasterize_to_pixels_bwd_tensor<Vanilla3DGS, true, true>(
            splats_tuple, backgrounds, masks,
            image_width, image_height, tile_offsets, flatten_ids,
            render_Ts, last_ids, render_outputs, render2_outputs, loss_map,
            v_render_outputs, v_render_alphas, v_distortion_outputs
        );
    return _rasterize_to_pixels_bwd_tensor<Vanilla3DGS, false, true>(
        splats_tuple, backgrounds, masks,
        image_width, image_height, tile_offsets, flatten_ids,
        render_Ts, last_ids, render_outputs, render2_outputs, loss_map,
        v_render_outputs, v_render_alphas, v_distortion_outputs
    );
}
