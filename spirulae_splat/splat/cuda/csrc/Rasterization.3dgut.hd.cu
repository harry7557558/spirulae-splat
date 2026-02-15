#include "RasterizationEval3DFwd.cuh"
#include "RasterizationEval3DBwd.cuh"


/*[AutoHeaderGeneratorExport]*/
std::tuple<
    Vanilla3DGUT::Screen::TensorTuple,
    std::optional<at::Tensor>,  // v_viewmats
    std::optional<Vanilla3DGUT::Screen::TensorTuple>,  // jacobian residual product
    std::optional<Vanilla3DGUT::Screen::TensorTuple>  // hessian diagonal
> rasterize_to_pixels_3dgut_bwd_with_hessian_diagonal(
    // Gaussian parameters
    Vanilla3DGUT::Screen::TensorTuple splats_tuple,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
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
    std::optional<at::Tensor> loss_map,  // [..., image_height, image_width, 1]
    // gradients of outputs
    Vanilla3DGUT::RenderOutput::TensorTuple v_render_outputs,
    const at::Tensor v_render_alphas, // [..., image_height, image_width, 1]
    std::optional<typename Vanilla3DGUT::RenderOutput::TensorTuple> v_distortion_outputs,
    bool need_viewmat_grad
) {
    if (v_distortion_outputs.has_value())
        return _rasterize_to_pixels_eval3d_bwd_tensor<Vanilla3DGUT, true, true>(
            splats_tuple,
            viewmats, intrins, camera_model, dist_coeffs,
            backgrounds, masks,
            image_width, image_height, tile_size, tile_offsets, flatten_ids,
            render_Ts, last_ids, render_outputs, render2_outputs, loss_map,
            v_render_outputs, v_render_alphas, v_distortion_outputs,
            need_viewmat_grad
        );
    return _rasterize_to_pixels_eval3d_bwd_tensor<Vanilla3DGUT, false, true>(
        splats_tuple,
        viewmats, intrins, camera_model, dist_coeffs,
        backgrounds, masks,
        image_width, image_height, tile_size, tile_offsets, flatten_ids,
        render_Ts, last_ids, render_outputs, render2_outputs, loss_map,
        v_render_outputs, v_render_alphas, v_distortion_outputs,
        need_viewmat_grad
    );
}
