#include "RasterizationBwd.cuh"

#include <gsplat/Utils.cuh>

#include "types.cuh"
#include "common.cuh"

#include <c10/cuda/CUDAStream.h>


template <
    typename SplatPrimitive,
    bool output_distortion,
    bool output_hessian_diagonal,
    bool output_accum_weight
>
void rasterize_to_pixels_bwd_kernel_wrapper(
    cudaStream_t stream,
    const uint32_t I,
    const uint32_t N,   // zero if packed
    const uint32_t n_isects,
    // fwd inputs
    const uint32_t *__restrict__ gaussian_ids,  // [nnz] optional, for packed mode
    const typename SplatPrimitive::WorldBuffer splat_wbuffer,
    const typename SplatPrimitive::ScreenBuffer splat_sbuffer,
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_width,
    const uint32_t tile_height,
    const int32_t *__restrict__ tile_offsets, // [..., tile_height, tile_width]
    const int32_t *__restrict__ flatten_ids,  // [n_isects]
    // fwd outputs
    const float *__restrict__ render_Ts,      // [..., image_height, image_width, 1]
    const int32_t *__restrict__ last_ids, // [..., image_height, image_width]
    RenderOutput::Buffer render_output_buffer,
    RenderOutput::Buffer render2_output_buffer,
    const float *__restrict__ loss_map_buffer,           // [..., image_height, image_width, 1]
    const float *__restrict__ accum_weight_map_buffer,           // [..., image_height, image_width, 1]
    // grad outputs
    RenderOutput::Buffer v_render_output_buffer,
    const float *__restrict__ v_render_Ts, // [..., image_height, image_width, 1]
    RenderOutput::Buffer v_distortions_output_buffer,
    // grad inputs
    typename SplatPrimitive::WorldBuffer v_splat_wbuffer,
    typename SplatPrimitive::ScreenBuffer v_splat_sbuffer,
    typename SplatPrimitive::WorldBuffer vr_splat_wbuffer,
    typename SplatPrimitive::ScreenBuffer vr_splat_sbuffer,
    typename SplatPrimitive::WorldBuffer h_splat_wbuffer,
    typename SplatPrimitive::ScreenBuffer h_splat_sbuffer,
    float *__restrict__ o_accum_weight
);

template <typename SplatPrimitive, bool output_distortion, bool output_hessian_diagonal, bool output_accum_weight>
inline void launch_rasterize_to_pixels_bwd_kernel(
    // Gaussian parameters
    const typename SplatPrimitive::WorldBuffer splat_wbuffer,
    const typename SplatPrimitive::ScreenBuffer splat_sbuffer,
    std::optional<at::Tensor> gaussian_ids,
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    // intersections
    const at::Tensor tile_offsets, // [..., tile_height, tile_width]
    const at::Tensor flatten_ids,  // [n_isects]
    // forward outputs
    const at::Tensor render_Ts, // [..., image_height, image_width, 1]
    const at::Tensor last_ids,      // [..., image_height, image_width]
    RenderOutput::Tensor *render_outputs,
    RenderOutput::Tensor *render2_outputs,
    const at::Tensor *loss_map,           // [..., image_height, image_width, 1]
    const at::Tensor *accum_weight_map,           // [..., image_height, image_width, 1]
    // gradients of outputs
    RenderOutput::Tensor v_render_outputs,
    const at::Tensor v_render_Ts, // [..., image_height, image_width, 1]
    RenderOutput::Tensor *v_distortion_outputs,
    // outputs
    typename SplatPrimitive::WorldBuffer v_splat_wbuffer,
    typename SplatPrimitive::ScreenBuffer v_splat_sbuffer,
    typename SplatPrimitive::WorldBuffer vr_splat_wbuffer,
    typename SplatPrimitive::ScreenBuffer vr_splat_sbuffer,
    typename SplatPrimitive::WorldBuffer h_splat_wbuffer,
    typename SplatPrimitive::ScreenBuffer h_splat_sbuffer,
    std::optional<at::Tensor> o_accum_weight
) {
    uint32_t N = gaussian_ids.has_value() ? 0 : splat_wbuffer.size();
    uint32_t I = render_Ts.numel() / (image_height * image_width); // number of images
    uint32_t tile_height = tile_offsets.size(-2);
    uint32_t tile_width = tile_offsets.size(-1);
    uint32_t n_isects = flatten_ids.size(0);
    // printf("N=%u I=%u th=%u tw=%u ni=%u\n", N, I, tile_height, tile_width, n_isects);

    if (n_isects == 0) {
        // skip the kernel launch if there are no elements
        return;
    }

    rasterize_to_pixels_bwd_kernel_wrapper<SplatPrimitive, output_distortion, output_hessian_diagonal, output_accum_weight>(
        (cudaStream_t)at::cuda::getCurrentCUDAStream(), I, N, n_isects,
        gaussian_ids.has_value() ? (uint32_t*)gaussian_ids.value().data_ptr<int32_t>() : nullptr,
        splat_wbuffer, splat_sbuffer,
        image_width, image_height, tile_width, tile_height,
        tile_offsets.data_ptr<int32_t>(), flatten_ids.data_ptr<int32_t>(),
        render_Ts.data_ptr<float>(), last_ids.data_ptr<int32_t>(),
        output_distortion ? *render_outputs : RenderOutput::Buffer(),
        output_distortion ? *render2_outputs : RenderOutput::Buffer(),
        output_hessian_diagonal ? loss_map->data_ptr<float>() : nullptr,
        output_accum_weight ? accum_weight_map->data_ptr<float>() : nullptr,
        v_render_outputs, v_render_Ts.data_ptr<float>(),
        output_distortion ? *v_distortion_outputs : RenderOutput::Buffer(),
        v_splat_wbuffer, v_splat_sbuffer, vr_splat_wbuffer, vr_splat_sbuffer, h_splat_wbuffer, h_splat_sbuffer,
        output_accum_weight ? o_accum_weight.value().data_ptr<float>() : nullptr
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


template<typename SplatPrimitive, bool output_distortion, bool output_hessian_diagonal, bool output_accum_weight>
inline std::tuple<
    TensorList, TensorList,  // gradient
    std::optional<TensorList>, std::optional<TensorList>,  // jacobian residual product
    std::optional<TensorList>, std::optional<TensorList>,  // hessian diagonal
    std::optional<at::Tensor>  // accum_weight
> _rasterize_to_pixels_bwd_tensor(
    // Gaussian parameters
    int64_t num_splats,
    TensorList splats_w,
    TensorList splats_s,
    std::optional<at::Tensor> gaussian_ids,
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    // intersections
    const at::Tensor tile_offsets, // [..., tile_height, tile_width]
    const at::Tensor flatten_ids,  // [n_isects]
    // forward outputs
    const at::Tensor render_Ts, // [..., image_height, image_width, 1]
    const at::Tensor last_ids,      // [..., image_height, image_width]
    std::optional<RenderOutput::TensorTuple> render_outputs_tuple,
    std::optional<RenderOutput::TensorTuple> render2_outputs_tuple,
    std::optional<at::Tensor> loss_map,  // [..., image_height, image_width, 1]
    std::optional<at::Tensor> accum_weight_map,  // [..., image_height, image_width, 1]
    // gradients of outputs
    RenderOutput::TensorTuple v_render_outputs,
    const at::Tensor v_render_Ts, // [..., image_height, image_width, 1]
    std::optional<RenderOutput::TensorTuple> v_distortion_outputs_tuple
) {
    DEVICE_GUARD(tile_offsets);
    CHECK_INPUT(tile_offsets);
    CHECK_INPUT(flatten_ids);
    CHECK_INPUT(render_Ts);
    CHECK_INPUT(last_ids);
    CHECK_INPUT(v_render_Ts);
    if (loss_map.has_value())
        CHECK_INPUT(loss_map.value());

    TensorList v_splats_w = SplatPrimitive::WorldBuffer::zeros_like(splats_w);
    TensorList v_splats_s = SplatPrimitive::ScreenBuffer::zeros_like(splats_s);

    std::optional<RenderOutput::Tensor> render_outputs = std::nullopt;
    std::optional<RenderOutput::Tensor> render2_outputs = std::nullopt;
    std::optional<RenderOutput::Tensor> v_distortion_outputs = std::nullopt;
    if (output_distortion) {
        render_outputs = render_outputs_tuple;
        render2_outputs = render2_outputs_tuple;
        v_distortion_outputs = v_distortion_outputs_tuple;
    }
    TensorList vr_splats_w;
    TensorList vr_splats_s;
    TensorList h_splats_w;
    TensorList h_splats_s;
    if (output_hessian_diagonal) {
        vr_splats_w = SplatPrimitive::WorldBuffer::zeros_like(splats_w);
        vr_splats_s = SplatPrimitive::ScreenBuffer::zeros_like(splats_s);
        h_splats_w = SplatPrimitive::WorldBuffer::zeros_like(splats_w);
        h_splats_s = SplatPrimitive::ScreenBuffer::zeros_like(splats_s);
    }
    std::optional<at::Tensor> o_accum_weight = std::nullopt;
    if (output_accum_weight) {
        o_accum_weight = at::empty({num_splats}, accum_weight_map.value().options());
        set_zero_tensor(o_accum_weight.value());
    }

    launch_rasterize_to_pixels_bwd_kernel
    <SplatPrimitive, output_distortion, output_hessian_diagonal, output_accum_weight>(
        splats_w, splats_s, gaussian_ids,
        image_width, image_height, tile_offsets, flatten_ids,
        render_Ts, last_ids,
        output_distortion ? &render_outputs.value() : nullptr,
        output_distortion ? &render2_outputs.value() : nullptr,
        output_hessian_diagonal ? &loss_map.value() : nullptr,
        output_accum_weight ? &accum_weight_map.value() : nullptr,
        v_render_outputs, v_render_Ts,
        output_distortion ? &v_distortion_outputs.value() : nullptr,
        v_splats_w, v_splats_s,
        vr_splats_w, vr_splats_s, h_splats_w, h_splats_s,
        o_accum_weight
    );

    if (output_hessian_diagonal)
        return std::make_tuple(v_splats_w, v_splats_s,
            (std::optional<TensorList>)vr_splats_w,
            (std::optional<TensorList>)vr_splats_s,
            (std::optional<TensorList>)h_splats_w,
            (std::optional<TensorList>)h_splats_s,
            o_accum_weight);
    return std::make_tuple(v_splats_w, v_splats_s,
        (std::optional<TensorList>)std::nullopt,
        (std::optional<TensorList>)std::nullopt,
        (std::optional<TensorList>)std::nullopt,
        (std::optional<TensorList>)std::nullopt,
        o_accum_weight);
}

template<typename SplatPrimitive>
inline std::tuple<
    TensorList, TensorList,
    std::optional<at::Tensor>  // accum_weight
> rasterize_to_pixels_bwd_tensor(
    // Gaussian parameters
    int64_t num_splats,
    TensorList splats_w,
    TensorList splats_s,
    std::optional<at::Tensor> gaussian_ids,
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    // intersections
    const at::Tensor tile_offsets, // [..., tile_height, tile_width]
    const at::Tensor flatten_ids,  // [n_isects]
    // forward outputs
    const at::Tensor render_Ts, // [..., image_height, image_width, 1]
    const at::Tensor last_ids,      // [..., image_height, image_width]
    std::optional<at::Tensor> accum_weight_map,  // [..., image_height, image_width, 1]
    // gradients of outputs
    RenderOutput::TensorTuple v_render_outputs,
    const at::Tensor v_render_Ts // [..., image_height, image_width, 1]
) {
    // TODO: add interface for output_distortion
    auto [v_splats_w, v_splats_s, vr_splats_w, vr_splats_s, h_splats_w, h_splats_s, accum_weight] = (
        accum_weight_map.has_value() ?
        _rasterize_to_pixels_bwd_tensor<SplatPrimitive, false, false, true> :
        _rasterize_to_pixels_bwd_tensor<SplatPrimitive, false, false, false>
    )(
        num_splats, splats_w, splats_s, gaussian_ids,
        image_width, image_height, tile_offsets, flatten_ids,
        render_Ts, last_ids, std::nullopt, std::nullopt, std::nullopt, accum_weight_map,
        v_render_outputs, v_render_Ts, std::nullopt
    );
    return std::make_tuple(v_splats_w, v_splats_s, accum_weight);
}



// ================
// Vanilla3DGS and Mip-Splatting
// ================

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    TensorList, TensorList,  // gradient
    std::optional<at::Tensor>  // accum_weight
> rasterize_to_pixels_3dgs_bwd(
    // Gaussian parameters
    int64_t num_splats,
    TensorList splats_w,
    TensorList splats_s,
    std::optional<at::Tensor> gaussian_ids,
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    // intersections
    const at::Tensor tile_offsets, // [..., tile_height, tile_width]
    const at::Tensor flatten_ids,  // [n_isects]
    // forward outputs
    const at::Tensor render_Ts, // [..., image_height, image_width, 1]
    const at::Tensor last_ids,      // [..., image_height, image_width]
    std::optional<at::Tensor> accum_weight_map,  // [..., image_height, image_width, 1]
    // gradients of outputs
    RenderOutput::TensorTuple v_render_outputs,
    const at::Tensor v_render_Ts // [..., image_height, image_width, 1]
) {
    return rasterize_to_pixels_bwd_tensor<Vanilla3DGS>(
        num_splats, splats_w, splats_s, gaussian_ids,
        image_width, image_height, tile_offsets, flatten_ids,
        render_Ts, last_ids, accum_weight_map, v_render_outputs, v_render_Ts
    );
}


/*[AutoHeaderGeneratorExport]*/
std::tuple<
    TensorList, TensorList,  // gradient
    std::optional<TensorList>, std::optional<TensorList>,  // jacobian residual product
    std::optional<TensorList>, std::optional<TensorList>,  // hessian diagonal
    std::optional<at::Tensor>  // accum_weight
> rasterize_to_pixels_3dgs_bwd_with_hessian_diagonal(
    // Gaussian parameters
    int64_t num_splats,
    TensorList splats_w,
    TensorList splats_s,
    std::optional<at::Tensor> gaussian_ids,
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    // intersections
    const at::Tensor tile_offsets, // [..., tile_height, tile_width]
    const at::Tensor flatten_ids,  // [n_isects]
    // forward outputs
    const at::Tensor render_Ts, // [..., image_height, image_width, 1]
    const at::Tensor last_ids,      // [..., image_height, image_width]
    std::optional<RenderOutput::TensorTuple> render_outputs,
    std::optional<RenderOutput::TensorTuple> render2_outputs,
    std::optional<at::Tensor> loss_map,  // [..., image_height, image_width, 1]
    std::optional<at::Tensor> accum_weight_map,  // [..., image_height, image_width, 1]
    // gradients of outputs
    RenderOutput::TensorTuple v_render_outputs,
    const at::Tensor v_render_Ts, // [..., image_height, image_width, 1]
    std::optional<RenderOutput::TensorTuple> v_distortion_outputs
) {
    using Fn = decltype(&_rasterize_to_pixels_bwd_tensor<Vanilla3DGS, false, true, false>);
    static constexpr Fn funcs[2][2] = { {
        _rasterize_to_pixels_bwd_tensor<Vanilla3DGS, false, true, false>,
        _rasterize_to_pixels_bwd_tensor<Vanilla3DGS, false, true, true>,
    }, {
        _rasterize_to_pixels_bwd_tensor<Vanilla3DGS, true, true, false>,
        _rasterize_to_pixels_bwd_tensor<Vanilla3DGS, true, true, true>,
    } };
    return funcs[v_distortion_outputs.has_value()][accum_weight_map.has_value()](
        num_splats, splats_w, splats_s, gaussian_ids,
        image_width, image_height, tile_offsets, flatten_ids,
        render_Ts, last_ids, render_outputs, render2_outputs, loss_map, accum_weight_map,
        v_render_outputs, v_render_Ts, v_distortion_outputs
    );
}
