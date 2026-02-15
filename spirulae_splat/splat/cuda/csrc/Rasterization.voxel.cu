#include "RasterizationEval3DFwd.cuh"
#include "RasterizationEval3DBwd.cuh"


/*[AutoHeaderGeneratorExport]*/
std::tuple<
    VoxelPrimitive::RenderOutput::TensorTuple,
    at::Tensor,
    at::Tensor,
    std::optional<VoxelPrimitive::RenderOutput::TensorTuple>,
    std::optional<VoxelPrimitive::RenderOutput::TensorTuple>,
    std::optional<at::Tensor>
> rasterize_to_pixels_voxel_eval3d_fwd(
    // Gaussian parameters
    VoxelPrimitive::Screen::TensorTuple splats_tuple,
    const at::Tensor viewmats,      // [..., C, 4, 4]
    const at::Tensor intrins,       // [..., C, 4], fx, fy, cx, cy
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
    const at::Tensor flatten_ids   // [n_isects]
) {
    return rasterize_to_pixels_eval3d_fwd_tensor<VoxelPrimitive, true, true>(
        splats_tuple,
        viewmats, intrins, camera_model, dist_coeffs,
        backgrounds, masks,
        image_width, image_height, tile_size,
        tile_offsets, flatten_ids
    );
}

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    VoxelPrimitive::Screen::TensorTuple,
    std::optional<at::Tensor>  // v_viewmats
> rasterize_to_pixels_voxel_eval3d_bwd(
    // Gaussian parameters
    VoxelPrimitive::Screen::TensorTuple splats_tuple,
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
    std::optional<typename VoxelPrimitive::RenderOutput::TensorTuple> render_outputs,
    std::optional<typename VoxelPrimitive::RenderOutput::TensorTuple> render2_outputs,
    std::optional<at::Tensor> loss_map,  // [..., image_height, image_width, 1]
    // gradients of outputs
    VoxelPrimitive::RenderOutput::TensorTuple v_render_outputs,
    const at::Tensor v_render_alphas, // [..., image_height, image_width, 1]
    std::optional<typename VoxelPrimitive::RenderOutput::TensorTuple> v_distortion_outputs,
    bool need_viewmat_grad
) {
    return rasterize_to_pixels_eval3d_bwd_tensor<VoxelPrimitive, false>(
        splats_tuple,
        viewmats, intrins, camera_model, dist_coeffs,
        backgrounds, masks,
        image_width, image_height, tile_size, tile_offsets, flatten_ids,
        render_Ts, last_ids, render_outputs, render2_outputs, loss_map,
        v_render_outputs, v_render_alphas, v_distortion_outputs,
        need_viewmat_grad
    );
}
