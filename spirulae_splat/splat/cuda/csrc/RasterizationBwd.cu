#include "RasterizationBwd.cuh"

#include <Common.cuh>



template <
    typename SplatPrimitive,
    bool output_distortion,
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
    float *__restrict__ o_accum_weight
);

template <typename SplatPrimitive, bool output_distortion, bool output_accum_weight>
inline void launch_rasterize_to_pixels_bwd_kernel(
    // Gaussian parameters
    int64_t num_splats,  // = cur_num_splats; non-packed projection layout stride
    const typename SplatPrimitive::WorldBuffer splat_wbuffer,
    const typename SplatPrimitive::ScreenBuffer splat_sbuffer,
    DeviceVector<int32_t> gaussian_ids,
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
    RenderOutput::Tensor render2_outputs,
    const DeviceTensor3D<float> loss_map,
    const DeviceTensor3D<float> accum_weight_map,  // [I, H, W]
    // gradients of outputs
    RenderOutput::Tensor v_render_outputs,
    const DeviceTensor3D<float> v_render_Ts, // [..., image_height, image_width]
    RenderOutput::Tensor v_distortion_outputs,
    // outputs
    typename SplatPrimitive::WorldBuffer v_splat_wbuffer,
    typename SplatPrimitive::ScreenBuffer v_splat_sbuffer,
    DeviceTensor3D<float> o_accum_weight
) {
    // Use the runtime cur_num_splats for the non-packed stride. The forward
    // splat buffers (splat_wbuffer) are pre-allocated to max_num_splats, but
    // the projection kernel only writes the first cur_num_splats per camera;
    // splat_sid = cid * cur_num_splats + gid in flatten_ids, so the modulo
    // base must be cur_num_splats to recover gid (and keep o_accum_weight
    // indices < cur_num_splats).
    uint32_t N = gaussian_ids.data_ptr() ? 0 : (uint32_t)num_splats;
    uint32_t I = render_Ts.size<0>(); // number of images
    uint32_t tile_height = tile_offsets.size<1>();
    uint32_t tile_width = tile_offsets.size<2>();
    uint32_t n_isects = flatten_ids.size();

    if (n_isects == 0) return;

    rasterize_to_pixels_bwd_kernel_wrapper<SplatPrimitive, output_distortion, output_accum_weight>(
        (cudaStream_t)0, I, N, n_isects,
        (uint32_t*)gaussian_ids.data_ptr(),
        splat_wbuffer, splat_sbuffer,
        image_width, image_height, tile_width, tile_height,
        tile_offsets.data_ptr(), flatten_ids.data_ptr(),
        render_Ts.data_ptr(), last_ids.data_ptr(), render_outputs,
        render2_outputs.has_value() ? render2_outputs : RenderOutput::Buffer(),
        loss_map.data_ptr(),
        accum_weight_map.data_ptr(),
        v_render_outputs, v_render_Ts.data_ptr(),
        v_distortion_outputs.has_value() ? v_distortion_outputs : RenderOutput::Buffer(),
        v_splat_wbuffer, v_splat_sbuffer,
        o_accum_weight.data_ptr()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


template<typename SplatPrimitive, bool output_distortion, bool output_accum_weight>
inline std::tuple<
    std::vector<DeviceTensorFloatND>, std::vector<DeviceTensorFloatND>,  // gradient
    DeviceVector<float>  // accum_weight
> _rasterize_to_pixels_bwd_tensor(
    // Gaussian parameters
    int64_t num_splats,
    std::vector<DeviceTensorFloatND> splats_w,
    std::vector<DeviceTensorFloatND> splats_s,
    DeviceVector<int32_t> gaussian_ids,
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
    DeviceTensor3D<float> loss_map,
    DeviceTensor3D<float> accum_weight_map,  // [I, H, W]
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
    RenderOutput::Tensor render2_outputs;
    RenderOutput::Tensor v_distortion_outputs;
    if (output_distortion) {
        render2_outputs = render2_outputs_tuple.value();
        v_distortion_outputs = v_distortion_outputs_tuple.value();
    }
    DeviceVector<float> o_accum_weight;
    if (output_accum_weight) {
        o_accum_weight.resize("raster_bwd.accum_weight", num_splats);
        o_accum_weight.zero();
    }

    // launch kernel needs DeviceTensor3D for o_accum_weight to pass data_ptr
    // but it's really a flat per-splat buffer
    DeviceTensor3D<float> o_accum_weight_3d;
    if (output_accum_weight) {
        // Wrap the flat buffer as a dummy 3D tensor just for the launch interface
        TorchTensorView tv((uint64_t)o_accum_weight.data_ptr(), 4, {num_splats, 1, 1, 1});
        o_accum_weight_3d = DeviceTensor3D<float>(tv);
    }

    launch_rasterize_to_pixels_bwd_kernel
    <SplatPrimitive, output_distortion, output_accum_weight>(
        num_splats,
        splats_w, splats_s, gaussian_ids,
        image_width, image_height, tile_offsets, flatten_ids,
        render_Ts, last_ids, render_outputs,
        render2_outputs, loss_map, accum_weight_map,
        v_render_outputs, v_render_Ts,
        v_distortion_outputs,
        v_splats_w.value(), v_splats_s.value(),
        o_accum_weight_3d
    );

    return std::make_tuple(
        v_splats_w.value(), v_splats_s.value(),
        o_accum_weight);
}

template<typename SplatPrimitive>
inline std::tuple<
    std::vector<DeviceTensorFloatND>, std::vector<DeviceTensorFloatND>,
    DeviceVector<float>  // accum_weight
> rasterize_to_pixels_bwd_tensor(
    // Gaussian parameters
    int64_t num_splats,
    std::vector<DeviceTensorFloatND> splats_w,
    std::vector<DeviceTensorFloatND> splats_s,
    DeviceVector<int32_t> gaussian_ids,
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
    DeviceTensor3D<float> accum_weight_map,  // [I, H, W]
    // gradients of outputs
    RenderOutput::TensorTuple v_render_outputs,
    const DeviceTensor3D<float> v_render_Ts,
    std::optional<std::vector<DeviceTensorFloatND>> v_splats_w,
    std::optional<std::vector<DeviceTensorFloatND>> v_splats_s
) {
    // TODO: add interface for output_distortion
    auto [v_splats_w_1, v_splats_s_1, accum_weight] = (
        accum_weight_map.data_ptr() != nullptr ?
        _rasterize_to_pixels_bwd_tensor<SplatPrimitive, false, true> :
        _rasterize_to_pixels_bwd_tensor<SplatPrimitive, false, false>
    )(
        num_splats, splats_w, splats_s, gaussian_ids,
        image_width, image_height, tile_offsets, flatten_ids,
        render_Ts, last_ids, render_outputs_tuple, std::nullopt,
        DeviceTensor3D<float>(), accum_weight_map,
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
    DeviceVector<float>  // accum_weight
> rasterize_to_pixels_3dgs_bwd(
    // Gaussian parameters
    int64_t num_splats,
    std::vector<DeviceTensorFloatND> splats_w,
    std::vector<DeviceTensorFloatND> splats_s,
    DeviceVector<int32_t> gaussian_ids,
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
    DeviceTensor3D<float> accum_weight_map,  // [I, H, W]
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

