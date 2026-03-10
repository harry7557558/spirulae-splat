#include "RasterizationBwd.cuh"

#include <gsplat/Utils.cuh>

#include "types.cuh"
#include "common.cuh"

#include <c10/cuda/CUDAStream.h>


template <
    typename SplatPrimitive,
    bool output_distortion,
    bool output_hessian_diagonal
>
void rasterize_to_pixels_bwd_kernel_wrapper(
    cudaStream_t stream,
    const uint32_t I,
    const uint32_t n_isects,
    // fwd inputs
    typename SplatPrimitive::Screen::Buffer splat_buffer,
    const float *__restrict__ backgrounds, // [..., CDIM] or [nnz, CDIM]
    const bool *__restrict__ masks,           // [..., tile_height, tile_width]
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_width,
    const uint32_t tile_height,
    const int32_t *__restrict__ tile_offsets, // [..., tile_height, tile_width]
    const int32_t *__restrict__ flatten_ids,  // [n_isects]
    // fwd outputs
    const float
        *__restrict__ render_Ts,      // [..., image_height, image_width, 1]
    const int32_t *__restrict__ last_ids, // [..., image_height, image_width]
    typename SplatPrimitive::RenderOutput::Buffer render_output_buffer,
    typename SplatPrimitive::RenderOutput::Buffer render2_output_buffer,
    const float *__restrict__ loss_map_buffer,           // [..., image_height, image_width, 1]
    // grad outputs
    typename SplatPrimitive::RenderOutput::Buffer v_render_output_buffer,
    const float *__restrict__ v_render_alphas, // [..., image_height, image_width, 1]
    typename SplatPrimitive::RenderOutput::Buffer v_distortions_output_buffer,
    // grad inputs
    typename SplatPrimitive::Screen::Buffer v_splat_buffer,
    typename SplatPrimitive::Screen::Buffer vr_splat_buffer,
    typename SplatPrimitive::Screen::Buffer h_splat_buffer
);

template <typename SplatPrimitive, bool output_distortion, bool output_hessian_diagonal>
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
    typename SplatPrimitive::RenderOutput::Tensor *render_outputs,
    typename SplatPrimitive::RenderOutput::Tensor *render2_outputs,
    const at::Tensor *loss_map,           // [..., image_height, image_width, 1]
    // gradients of outputs
    typename SplatPrimitive::RenderOutput::Tensor v_render_outputs,
    const at::Tensor v_render_alphas, // [..., image_height, image_width, 1]
    typename SplatPrimitive::RenderOutput::Tensor *v_distortion_outputs,
    // outputs
    typename SplatPrimitive::Screen::Tensor v_splats,
    typename SplatPrimitive::Screen::Tensor *vr_splats,
    typename SplatPrimitive::Screen::Tensor *h_splats
) {
    bool packed = splats.isPacked();
    uint32_t N = packed ? 0 : splats.size(); // number of gaussians
    uint32_t I = render_Ts.numel() / (image_height * image_width); // number of images
    uint32_t tile_height = tile_offsets.size(-2);
    uint32_t tile_width = tile_offsets.size(-1);
    uint32_t n_isects = flatten_ids.size(0);

    if (n_isects == 0) {
        // skip the kernel launch if there are no elements
        return;
    }

    typename SplatPrimitive::Screen::Buffer vr_splats_buffer;
    typename SplatPrimitive::Screen::Buffer h_splats_buffer;
    if (output_hessian_diagonal) {
        vr_splats_buffer = vr_splats->buffer();
        h_splats_buffer = h_splats->buffer();
    }

    rasterize_to_pixels_bwd_kernel_wrapper<SplatPrimitive, output_distortion, output_hessian_diagonal>(
        (cudaStream_t)at::cuda::getCurrentCUDAStream(), I, n_isects,
        splats.buffer(),
        backgrounds.has_value() ? backgrounds.value().data_ptr<float>() : nullptr,
        masks.has_value() ? masks.value().data_ptr<bool>() : nullptr,
        image_width, image_height, tile_width, tile_height,
        tile_offsets.data_ptr<int32_t>(), flatten_ids.data_ptr<int32_t>(),
        render_Ts.data_ptr<float>(), last_ids.data_ptr<int32_t>(),
        output_distortion ? render_outputs->buffer() : typename SplatPrimitive::RenderOutput::Buffer(),
        output_distortion ? render2_outputs->buffer() : typename SplatPrimitive::RenderOutput::Buffer(),
        output_hessian_diagonal ? loss_map->data_ptr<float>() : nullptr,
        v_render_outputs.buffer(), v_render_alphas.data_ptr<float>(),
        output_distortion ? v_distortion_outputs->buffer() : typename SplatPrimitive::RenderOutput::Buffer(),
        v_splats.buffer(), vr_splats_buffer, h_splats_buffer
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


template<typename SplatPrimitive, bool output_distortion, bool output_hessian_diagonal>
inline std::tuple<
    typename SplatPrimitive::Screen::TensorTuple,
    std::optional<typename SplatPrimitive::Screen::TensorTuple>,  // jacobian residual product
    std::optional<typename SplatPrimitive::Screen::TensorTuple>  // hessian diagonal
> _rasterize_to_pixels_bwd_tensor(
    // Gaussian parameters
    typename SplatPrimitive::Screen::TensorTuple splats_tuple,
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
    std::optional<typename SplatPrimitive::RenderOutput::TensorTuple> render_outputs_tuple,
    std::optional<typename SplatPrimitive::RenderOutput::TensorTuple> render2_outputs_tuple,
    std::optional<at::Tensor> loss_map,  // [..., image_height, image_width, 1]
    // gradients of outputs
    typename SplatPrimitive::RenderOutput::TensorTuple v_render_outputs,
    const at::Tensor v_render_alphas, // [..., image_height, image_width, 1]
    std::optional<typename SplatPrimitive::RenderOutput::TensorTuple> v_distortion_outputs_tuple
) {
    DEVICE_GUARD(tile_offsets);
    CHECK_INPUT(tile_offsets);
    CHECK_INPUT(flatten_ids);
    CHECK_INPUT(render_Ts);
    CHECK_INPUT(last_ids);
    CHECK_INPUT(v_render_alphas);
    if (backgrounds.has_value())
        CHECK_INPUT(backgrounds.value());
    if (masks.has_value())
        CHECK_INPUT(masks.value());
    if (loss_map.has_value())
        CHECK_INPUT(loss_map.value());

    typename SplatPrimitive::Screen::Tensor splats(splats_tuple);
    typename SplatPrimitive::Screen::Tensor v_splats = splats.allocRasterBwd();

    std::optional<typename SplatPrimitive::RenderOutput::Tensor> render_outputs = std::nullopt;
    std::optional<typename SplatPrimitive::RenderOutput::Tensor> render2_outputs = std::nullopt;
    std::optional<typename SplatPrimitive::RenderOutput::Tensor> v_distortion_outputs = std::nullopt;
    if (output_distortion) {
        render_outputs = render_outputs_tuple;
        render2_outputs = render2_outputs_tuple;
        v_distortion_outputs = v_distortion_outputs_tuple;
    }
    std::optional<typename SplatPrimitive::Screen::Tensor> vr_splats = std::nullopt;
    std::optional<typename SplatPrimitive::Screen::Tensor> h_splats = std::nullopt;
    if (output_hessian_diagonal) {
        vr_splats = splats.allocRasterBwd();
        h_splats = splats.allocRasterBwd();
    }

    launch_rasterize_to_pixels_bwd_kernel
    <SplatPrimitive, output_distortion, output_hessian_diagonal>(
        splats,
        backgrounds, masks,
        image_width, image_height, tile_offsets, flatten_ids,
        render_Ts, last_ids,
        output_distortion ? &render_outputs.value() : nullptr,
        output_distortion ? &render2_outputs.value() : nullptr,
        output_hessian_diagonal ? &loss_map.value() : nullptr,
        v_render_outputs, v_render_alphas,
        output_distortion ? &v_distortion_outputs.value() : nullptr,
        v_splats,
        output_hessian_diagonal ? &vr_splats.value() : nullptr,
        output_hessian_diagonal ? &h_splats.value() : nullptr
    );

    if (output_hessian_diagonal)
        return std::make_tuple(v_splats.tupleRasterBwd(),
            (std::optional<typename SplatPrimitive::Screen::TensorTuple>)vr_splats.value().tupleRasterBwd(),
            (std::optional<typename SplatPrimitive::Screen::TensorTuple>)h_splats.value().tupleRasterBwd());
    return std::make_tuple(v_splats.tupleRasterBwd(),
        (std::optional<typename SplatPrimitive::Screen::TensorTuple>)std::nullopt,
        (std::optional<typename SplatPrimitive::Screen::TensorTuple>)std::nullopt);
}

template<typename SplatPrimitive>
typename SplatPrimitive::Screen::TensorTuple
inline rasterize_to_pixels_bwd_tensor(
    // Gaussian parameters
    typename SplatPrimitive::Screen::TensorTuple splats_tuple,
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
    typename SplatPrimitive::RenderOutput::TensorTuple v_render_outputs,
    const at::Tensor v_render_alphas // [..., image_height, image_width, 1]
) {
    // TODO: add interface for output_distortion
    auto [v_splats, vr_splats, h_splats] =
        _rasterize_to_pixels_bwd_tensor<SplatPrimitive, false, false>
    (
        splats_tuple,
        backgrounds, masks,
        image_width, image_height, tile_offsets, flatten_ids,
        render_Ts, last_ids, std::nullopt, std::nullopt, std::nullopt,
        v_render_outputs, v_render_alphas, std::nullopt
    );
    return v_splats;
}



// ================
// Vanilla3DGS
// ================

/*[AutoHeaderGeneratorExport]*/
Vanilla3DGS::Screen::TensorTuple rasterize_to_pixels_3dgs_bwd(
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
    // gradients of outputs
    Vanilla3DGS::RenderOutput::TensorTuple v_render_outputs,
    const at::Tensor v_render_alphas // [..., image_height, image_width, 1]
) {
    return rasterize_to_pixels_bwd_tensor<Vanilla3DGS>(
        splats_tuple, backgrounds, masks,
        image_width, image_height, tile_offsets, flatten_ids,
        render_Ts, last_ids, v_render_outputs, v_render_alphas
    );
}


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



// ================
// MipSplatting
// ================

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
    // return rasterize_to_pixels_bwd_tensor<MipSplatting>(
    return rasterize_to_pixels_bwd_tensor<Vanilla3DGS>(
        splats_tuple, backgrounds, masks,
        image_width, image_height, tile_offsets, flatten_ids,
        render_Ts, last_ids, v_render_outputs, v_render_alphas
    );
}

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
        // return _rasterize_to_pixels_bwd_tensor<MipSplatting, true, true>(
        return _rasterize_to_pixels_bwd_tensor<Vanilla3DGS, true, true>(
            splats_tuple, backgrounds, masks,
            image_width, image_height, tile_offsets, flatten_ids,
            render_Ts, last_ids, render_outputs, render2_outputs, loss_map,
            v_render_outputs, v_render_alphas, v_distortion_outputs
        );
    // return _rasterize_to_pixels_bwd_tensor<MipSplatting, false, true>(
    return _rasterize_to_pixels_bwd_tensor<Vanilla3DGS, false, true>(
        splats_tuple, backgrounds, masks,
        image_width, image_height, tile_offsets, flatten_ids,
        render_Ts, last_ids, render_outputs, render2_outputs, loss_map,
        v_render_outputs, v_render_alphas, v_distortion_outputs
    );
}
