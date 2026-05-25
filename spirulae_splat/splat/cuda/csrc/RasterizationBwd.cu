#include "RasterizationBwd.cuh"

#include <gsplat/Utils.cuh>

#include "types.cuh"
#include "common.cuh"



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
    std::optional<DeviceVector<int32_t>> gaussian_ids,
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    // intersections
    const DeviceTensor3D<int32_t> tile_offsets, // [I, tile_height, tile_width]
    const DeviceVector<int32_t> flatten_ids,    // [n_isects]
    // forward outputs
    const DeviceTensor3D<float> render_Ts,   // [I, image_height, image_width]
    const DeviceTensor3D<int32_t> last_ids,  // [I, image_height, image_width]
    RenderOutput::Tensor render_outputs,
    RenderOutput::Tensor *render2_outputs,
    const DeviceTensor3D<float> *loss_map,   // [..., image_height, image_width]
    const DeviceTensor3D<float> *accum_weight_map,  // [I, H, W]
    // gradients of outputs
    RenderOutput::Tensor v_render_outputs,
    const DeviceTensor3D<float> v_render_Ts, // [..., image_height, image_width]
    RenderOutput::Tensor *v_distortion_outputs,
    // outputs
    typename SplatPrimitive::WorldBuffer v_splat_wbuffer,
    typename SplatPrimitive::ScreenBuffer v_splat_sbuffer,
    typename SplatPrimitive::WorldBuffer vr_splat_wbuffer,
    typename SplatPrimitive::ScreenBuffer vr_splat_sbuffer,
    typename SplatPrimitive::WorldBuffer h_splat_wbuffer,
    typename SplatPrimitive::ScreenBuffer h_splat_sbuffer,
    std::optional<DeviceTensor3D<float>> *o_accum_weight
) {
    uint32_t N = gaussian_ids.has_value() ? 0 : splat_wbuffer.size();
    uint32_t I = render_Ts.size<0>(); // number of images
    uint32_t tile_height = tile_offsets.size<1>();
    uint32_t tile_width = tile_offsets.size<2>();
    uint32_t n_isects = flatten_ids.size();

    if (n_isects == 0) {
        // skip the kernel launch if there are no elements
        return;
    }

    rasterize_to_pixels_bwd_kernel_wrapper<SplatPrimitive, output_distortion, output_hessian_diagonal, output_accum_weight>(
        (cudaStream_t)0, I, N, n_isects,
        gaussian_ids.has_value() ? (uint32_t*)gaussian_ids.value().data_ptr() : nullptr,
        splat_wbuffer, splat_sbuffer,
        image_width, image_height, tile_width, tile_height,
        tile_offsets.data_ptr(), flatten_ids.data_ptr(),
        render_Ts.data_ptr(), last_ids.data_ptr(), render_outputs,
        output_distortion ? *render2_outputs : RenderOutput::Buffer(),
        output_hessian_diagonal ? loss_map->data_ptr() : nullptr,
        output_accum_weight ? accum_weight_map->data_ptr() : nullptr,
        v_render_outputs, v_render_Ts.data_ptr(),
        output_distortion ? *v_distortion_outputs : RenderOutput::Buffer(),
        v_splat_wbuffer, v_splat_sbuffer, vr_splat_wbuffer, vr_splat_sbuffer, h_splat_wbuffer, h_splat_sbuffer,
        output_accum_weight ? o_accum_weight->value().data_ptr() : nullptr
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


template<typename SplatPrimitive, bool output_distortion, bool output_hessian_diagonal, bool output_accum_weight>
inline std::tuple<
    std::vector<DeviceTensorFloatND>, std::vector<DeviceTensorFloatND>,  // gradient
    std::optional<std::vector<DeviceTensorFloatND>>, std::optional<std::vector<DeviceTensorFloatND>>,  // jacobian residual product
    std::optional<std::vector<DeviceTensorFloatND>>, std::optional<std::vector<DeviceTensorFloatND>>,  // hessian diagonal
    std::optional<DeviceTensor3D<float>>  // accum_weight
> _rasterize_to_pixels_bwd_tensor(
    // Gaussian parameters
    int64_t num_splats,
    std::vector<DeviceTensorFloatND> splats_w,
    std::vector<DeviceTensorFloatND> splats_s,
    std::optional<DeviceVector<int32_t>> gaussian_ids,
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    // intersections
    const DeviceTensor3D<int32_t> tile_offsets, // [I, tile_height, tile_width]
    const DeviceVector<int32_t> flatten_ids,    // [n_isects]
    // forward outputs
    const DeviceTensor3D<float> render_Ts,  // [I, image_height, image_width]
    const DeviceTensor3D<int32_t> last_ids, // [I, image_height, image_width]
    RenderOutput::TensorTuple render_outputs_tuple,
    std::optional<RenderOutput::TensorTuple> render2_outputs_tuple,
    std::optional<DeviceTensor3D<float>> loss_map,
    std::optional<DeviceTensor3D<float>> accum_weight_map,  // [I, H, W]
    // gradients of outputs
    RenderOutput::TensorTuple v_render_outputs,
    const DeviceTensor3D<float> v_render_Ts,
    std::optional<RenderOutput::TensorTuple> v_distortion_outputs_tuple,
    std::optional<std::vector<DeviceTensorFloatND>> v_splats_w,
    std::optional<std::vector<DeviceTensorFloatND>> v_splats_s
) {
    if (!v_splats_w.has_value())
        v_splats_w = SplatPrimitive::WorldBuffer::zeros_pool(splats_w, "raster_bwd.v_world");
    if (!v_splats_s.has_value())
        v_splats_s = SplatPrimitive::ScreenBuffer::zeros_pool(splats_s, "raster_bwd.v_screen");

    RenderOutput::Tensor render_outputs = render_outputs_tuple;
    std::optional<RenderOutput::Tensor> render2_outputs = std::nullopt;
    std::optional<RenderOutput::Tensor> v_distortion_outputs = std::nullopt;
    if (output_distortion) {
        render2_outputs = render2_outputs_tuple;
        v_distortion_outputs = v_distortion_outputs_tuple;
    }
    std::vector<DeviceTensorFloatND> vr_splats_w;
    std::vector<DeviceTensorFloatND> vr_splats_s;
    std::vector<DeviceTensorFloatND> h_splats_w;
    std::vector<DeviceTensorFloatND> h_splats_s;
    if (output_hessian_diagonal) {
        vr_splats_w = SplatPrimitive::WorldBuffer::zeros_pool(splats_w, "raster_bwd.vr_world");
        vr_splats_s = SplatPrimitive::ScreenBuffer::zeros_pool(splats_s, "raster_bwd.vr_screen");
        h_splats_w = SplatPrimitive::WorldBuffer::zeros_pool(splats_w, "raster_bwd.h_world");
        h_splats_s = SplatPrimitive::ScreenBuffer::zeros_pool(splats_s, "raster_bwd.h_screen");
    }
    std::optional<DeviceTensor3D<float>> o_accum_weight = std::nullopt;
    if (output_accum_weight) {
        DeviceTensor3D<float> _aw;
        _aw.resize("raster_bwd.accum_weight",
            accum_weight_map.value().size<0>(),
            accum_weight_map.value().size<1>(),
            accum_weight_map.value().size<2>());
        _aw.zero();
        o_accum_weight = _aw;
    }

    launch_rasterize_to_pixels_bwd_kernel
    <SplatPrimitive, output_distortion, output_hessian_diagonal, output_accum_weight>(
        splats_w, splats_s, gaussian_ids,
        image_width, image_height, tile_offsets, flatten_ids,
        render_Ts, last_ids, render_outputs,
        output_distortion ? &render2_outputs.value() : nullptr,
        output_hessian_diagonal ? &loss_map.value() : nullptr,
        output_accum_weight ? &accum_weight_map.value() : nullptr,
        v_render_outputs, v_render_Ts,
        output_distortion ? &v_distortion_outputs.value() : nullptr,
        v_splats_w.value(), v_splats_s.value(),
        vr_splats_w, vr_splats_s, h_splats_w, h_splats_s,
        &o_accum_weight
    );

    if (output_hessian_diagonal)
        return std::make_tuple(
            v_splats_w.value(), v_splats_s.value(),
            (std::optional<std::vector<DeviceTensorFloatND>>)vr_splats_w,
            (std::optional<std::vector<DeviceTensorFloatND>>)vr_splats_s,
            (std::optional<std::vector<DeviceTensorFloatND>>)h_splats_w,
            (std::optional<std::vector<DeviceTensorFloatND>>)h_splats_s,
            o_accum_weight);
    return std::make_tuple(
        v_splats_w.value(), v_splats_s.value(),
        (std::optional<std::vector<DeviceTensorFloatND>>)std::nullopt,
        (std::optional<std::vector<DeviceTensorFloatND>>)std::nullopt,
        (std::optional<std::vector<DeviceTensorFloatND>>)std::nullopt,
        (std::optional<std::vector<DeviceTensorFloatND>>)std::nullopt,
        o_accum_weight);
}

template<typename SplatPrimitive>
inline std::tuple<
    std::vector<DeviceTensorFloatND>, std::vector<DeviceTensorFloatND>,
    std::optional<DeviceTensor3D<float>>  // accum_weight
> rasterize_to_pixels_bwd_tensor(
    // Gaussian parameters
    int64_t num_splats,
    std::vector<DeviceTensorFloatND> splats_w,
    std::vector<DeviceTensorFloatND> splats_s,
    std::optional<DeviceVector<int32_t>> gaussian_ids,
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    // intersections
    const DeviceTensor3D<int32_t> tile_offsets, // [I, tile_height, tile_width]
    const DeviceVector<int32_t> flatten_ids,    // [n_isects]
    // forward outputs
    const DeviceTensor3D<float> render_Ts,  // [I, image_height, image_width]
    const DeviceTensor3D<int32_t> last_ids, // [I, image_height, image_width]
    RenderOutput::TensorTuple render_outputs_tuple,
    std::optional<DeviceTensor3D<float>> accum_weight_map,  // [I, H, W]
    // gradients of outputs
    RenderOutput::TensorTuple v_render_outputs,
    const DeviceTensor3D<float> v_render_Ts,
    std::optional<std::vector<DeviceTensorFloatND>> v_splats_w,
    std::optional<std::vector<DeviceTensorFloatND>> v_splats_s
) {
    // TODO: add interface for output_distortion
    auto [v_splats_w_1, v_splats_s_1, vr_splats_w, vr_splats_s, h_splats_w, h_splats_s, accum_weight] = (
        accum_weight_map.has_value() ?
        _rasterize_to_pixels_bwd_tensor<SplatPrimitive, false, false, true> :
        _rasterize_to_pixels_bwd_tensor<SplatPrimitive, false, false, false>
    )(
        num_splats, splats_w, splats_s, gaussian_ids,
        image_width, image_height, tile_offsets, flatten_ids,
        render_Ts, last_ids, render_outputs_tuple, std::nullopt, std::nullopt, accum_weight_map,
        v_render_outputs, v_render_Ts, std::nullopt, v_splats_w, v_splats_s
    );
    return std::make_tuple(v_splats_w_1, v_splats_s_1, accum_weight);
}



// ================
// Vanilla3DGS and Mip-Splatting
// ================

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    std::vector<DeviceTensorFloatND>, std::vector<DeviceTensorFloatND>,  // gradient
    std::optional<DeviceTensor3D<float>>  // accum_weight
> rasterize_to_pixels_3dgs_bwd(
    // Gaussian parameters
    int64_t num_splats,
    std::vector<DeviceTensorFloatND> splats_w,
    std::vector<DeviceTensorFloatND> splats_s,
    std::optional<DeviceVector<int32_t>> gaussian_ids,
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    // intersections
    const DeviceTensor3D<int32_t> tile_offsets, // [I, tile_height, tile_width]
    const DeviceVector<int32_t> flatten_ids,    // [n_isects]
    // forward outputs
    const DeviceTensor3D<float> render_Ts,  // [I, image_height, image_width]
    const DeviceTensor3D<int32_t> last_ids, // [I, image_height, image_width]
    RenderOutput::TensorTuple render_outputs_tuple,
    std::optional<DeviceTensor3D<float>> accum_weight_map,  // [I, H, W]
    // gradients of outputs
    RenderOutput::TensorTuple v_render_outputs,
    const DeviceTensor3D<float> v_render_Ts,
    std::optional<std::vector<DeviceTensorFloatND>> v_splats_w,
    std::optional<std::vector<DeviceTensorFloatND>> v_splats_s
) {
    return rasterize_to_pixels_bwd_tensor<Vanilla3DGS<0>>(
        num_splats, splats_w, splats_s, gaussian_ids,
        image_width, image_height, tile_offsets, flatten_ids,
        render_Ts, last_ids, render_outputs_tuple, accum_weight_map, v_render_outputs, v_render_Ts,
        v_splats_w, v_splats_s
    );
}


/*[AutoHeaderGeneratorExport]*/
std::tuple<
    std::vector<DeviceTensorFloatND>, std::vector<DeviceTensorFloatND>,  // gradient
    std::optional<std::vector<DeviceTensorFloatND>>, std::optional<std::vector<DeviceTensorFloatND>>,  // jacobian residual product
    std::optional<std::vector<DeviceTensorFloatND>>, std::optional<std::vector<DeviceTensorFloatND>>,  // hessian diagonal
    std::optional<DeviceTensor3D<float>>  // accum_weight
> rasterize_to_pixels_3dgs_bwd_with_hessian_diagonal(
    // Gaussian parameters
    int64_t num_splats,
    std::vector<DeviceTensorFloatND> splats_w,
    std::vector<DeviceTensorFloatND> splats_s,
    std::optional<DeviceVector<int32_t>> gaussian_ids,
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    // intersections
    const DeviceTensor3D<int32_t> tile_offsets, // [I, tile_height, tile_width]
    const DeviceVector<int32_t> flatten_ids,    // [n_isects]
    // forward outputs
    const DeviceTensor3D<float> render_Ts,  // [I, image_height, image_width]
    const DeviceTensor3D<int32_t> last_ids, // [I, image_height, image_width]
    RenderOutput::TensorTuple render_outputs,
    std::optional<RenderOutput::TensorTuple> render2_outputs,
    std::optional<DeviceTensor3D<float>> loss_map,
    std::optional<DeviceTensor3D<float>> accum_weight_map,  // [I, H, W]
    // gradients of outputs
    RenderOutput::TensorTuple v_render_outputs,
    const DeviceTensor3D<float> v_render_Ts,
    std::optional<RenderOutput::TensorTuple> v_distortion_outputs,
    std::optional<std::vector<DeviceTensorFloatND>> v_splats_w,
    std::optional<std::vector<DeviceTensorFloatND>> v_splats_s
) {
    using Fn = decltype(&_rasterize_to_pixels_bwd_tensor<Vanilla3DGS<0>, false, true, false>);
    static constexpr Fn funcs[2][2] = { {
        _rasterize_to_pixels_bwd_tensor<Vanilla3DGS<0>, false, true, false>,
        _rasterize_to_pixels_bwd_tensor<Vanilla3DGS<0>, false, true, true>,
    }, {
        _rasterize_to_pixels_bwd_tensor<Vanilla3DGS<0>, true, true, false>,
        _rasterize_to_pixels_bwd_tensor<Vanilla3DGS<0>, true, true, true>,
    } };
    return funcs[v_distortion_outputs.has_value()][accum_weight_map.has_value()](
        num_splats, splats_w, splats_s, gaussian_ids,
        image_width, image_height, tile_offsets, flatten_ids,
        render_Ts, last_ids, render_outputs, render2_outputs, loss_map, accum_weight_map,
        v_render_outputs, v_render_Ts, v_distortion_outputs, v_splats_w, v_splats_s
    );
}
