#include "RasterizationEval3DBwd.cuh"

#include <gsplat/Utils.cuh>

#include "types.cuh"
#include "common.cuh"



template <
    typename SplatPrimitive,
    ssplat::CameraModelType camera_model,
    bool output_distortion,
    bool output_viewmat_grad,
    bool output_hessian_diagonal,
    bool output_accum_weight
>
void rasterize_to_pixels_eval3d_bwd_kernel_wrapper(
    cudaStream_t stream,
    const uint32_t I,
    const uint32_t N,   // zero if packed
    const uint32_t n_isects,
    // fwd inputs
    const uint32_t *__restrict__ gaussian_ids,  // [nnz] optional, for packed mode
    const typename SplatPrimitive::WorldBuffer splat_wbuffer,
    const typename SplatPrimitive::ScreenBuffer splat_sbuffer,
    const float *__restrict__ viewmats, // [B, C, 4, 4]
    const float4 *__restrict__ intrins,  // [B, C, 4], fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
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
    float *__restrict__ o_accum_weight,
    float *__restrict__ v_viewmats // [B, C, 4, 4]
);


template <typename SplatPrimitive, bool output_distortion, bool output_hessian_diagonal, bool output_accum_weight>
inline void launch_rasterize_to_pixels_eval3d_bwd_kernel(
    // Gaussian parameters
    int64_t num_splats,  // = cur_num_splats; non-packed projection layout stride
    const typename SplatPrimitive::WorldBuffer splat_wbuffer,
    const typename SplatPrimitive::ScreenBuffer splat_sbuffer,
    DeviceVector<int32_t> gaussian_ids,
    TorchTensorView viewmats,              // [..., C, 4, 4]
    TorchTensorView intrins,  // [..., C, 4], fx, fy, cx, cy
    const ssplat::CameraModelType camera_model,
    const TorchTensorView dist_coeffs,
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    // intersections
    const DeviceTensor3D<int32_t> tile_offsets, // [I, tile_height, tile_width]
    const DeviceVector<int32_t> flatten_ids,    // [n_isects]
    // forward outputs
    const DeviceTensor3D<float> render_Ts,  // [I, image_height, image_width]
    const DeviceTensor3D<int32_t> last_ids, // [I, image_height, image_width]
    RenderOutput::Tensor render_outputs,
    RenderOutput::Tensor render2_outputs,
    const DeviceTensor3D<float> loss_map,
    const DeviceTensor3D<float> accum_weight_map,  // [I, H, W]
    // gradients of outputs
    RenderOutput::Tensor v_render_outputs,
    const DeviceTensor3D<float> v_render_Ts, // [..., image_height, image_width, 1]
    RenderOutput::Tensor v_distortion_outputs,
    // outputs
    typename SplatPrimitive::WorldBuffer v_splat_wbuffer,
    typename SplatPrimitive::ScreenBuffer v_splat_sbuffer,
    typename SplatPrimitive::WorldBuffer vr_splat_wbuffer,
    typename SplatPrimitive::ScreenBuffer vr_splat_sbuffer,
    typename SplatPrimitive::WorldBuffer h_splat_wbuffer,
    typename SplatPrimitive::ScreenBuffer h_splat_sbuffer,
    DeviceTensor3D<float> o_accum_weight,
    DeviceTensor2D<float> v_viewmats
) {
    // See RasterizationBwd.cu: splat_wbuffer.size() returns max_num_splats
    // (pre-allocated), but the projection layout uses cur_num_splats per
    // camera stride. Use the explicit num_splats so non-packed splat_wid
    // recovery (splat_sid % N) returns valid gids.
    uint32_t N = gaussian_ids.data_ptr() ? 0 : (uint32_t)num_splats;
    uint32_t I = render_Ts.size<0>(); // number of images
    uint32_t tile_height = tile_offsets.size<1>();
    uint32_t tile_width = tile_offsets.size<2>();
    uint32_t n_isects = flatten_ids.size();

    if (n_isects == 0) {
        // skip the kernel launch if there are no elements
        return;
    }

    #define _LAUNCH_ARGS ( \
            (cudaStream_t)0, I, N, n_isects, \
            (uint32_t*)gaussian_ids.data_ptr(), \
            splat_wbuffer, splat_sbuffer, \
            (const float*)std::get<0>(viewmats), (const float4*)std::get<0>(intrins), dist_coeffs, \
            image_width, image_height, tile_width, tile_height, \
            tile_offsets.data_ptr(), flatten_ids.data_ptr(), \
            render_Ts.data_ptr(), last_ids.data_ptr(), \
            render_outputs, \
            render2_outputs.has_value() ? render2_outputs : RenderOutput::Buffer(), \
            loss_map.data_ptr(), \
            accum_weight_map.data_ptr(), \
            v_render_outputs, v_render_Ts.data_ptr(), \
            v_distortion_outputs.has_value() ? v_distortion_outputs : RenderOutput::Buffer(), \
            v_splat_wbuffer, v_splat_sbuffer, vr_splat_wbuffer, vr_splat_sbuffer, h_splat_wbuffer, h_splat_sbuffer, \
            o_accum_weight.data_ptr(), \
            v_viewmats.data_ptr() \
        )

    if (camera_model == ssplat::CameraModelType::PINHOLE) {
        if (v_viewmats.data_ptr() != nullptr)
            rasterize_to_pixels_eval3d_bwd_kernel_wrapper<SplatPrimitive,
                ssplat::CameraModelType::PINHOLE, output_distortion, true, output_hessian_diagonal, output_accum_weight> _LAUNCH_ARGS;
        else
            rasterize_to_pixels_eval3d_bwd_kernel_wrapper<SplatPrimitive,
                ssplat::CameraModelType::PINHOLE, output_distortion, false, output_hessian_diagonal, output_accum_weight> _LAUNCH_ARGS;
    }
    else if (camera_model == ssplat::CameraModelType::FISHEYE) {
        if (v_viewmats.data_ptr() != nullptr)
            rasterize_to_pixels_eval3d_bwd_kernel_wrapper<SplatPrimitive,
                ssplat::CameraModelType::FISHEYE, output_distortion, true, output_hessian_diagonal, output_accum_weight> _LAUNCH_ARGS;
        else
            rasterize_to_pixels_eval3d_bwd_kernel_wrapper<SplatPrimitive,
                ssplat::CameraModelType::FISHEYE, output_distortion, false, output_hessian_diagonal, output_accum_weight> _LAUNCH_ARGS;
    }
    else if (camera_model == ssplat::CameraModelType::EQUISOLID) {
        if (v_viewmats.data_ptr() != nullptr)
            rasterize_to_pixels_eval3d_bwd_kernel_wrapper<SplatPrimitive,
                ssplat::CameraModelType::EQUISOLID, output_distortion, true, output_hessian_diagonal, output_accum_weight> _LAUNCH_ARGS;
        else
            rasterize_to_pixels_eval3d_bwd_kernel_wrapper<SplatPrimitive,
                ssplat::CameraModelType::EQUISOLID, output_distortion, false, output_hessian_diagonal, output_accum_weight> _LAUNCH_ARGS;
    }
    else
        throw std::runtime_error("Unsupported camera model");
    CHECK_DEVICE_ERROR(cudaGetLastError());

    #undef _LAUNCH_ARGS
}


template<typename SplatPrimitive, bool output_distortion, bool output_hessian_diagonal, bool output_accum_weight>
inline std::tuple<
    std::vector<DeviceTensorFloatND>, std::vector<DeviceTensorFloatND>,  // gradient
    DeviceTensor2D<float>,  // v_viewmats
    std::optional<std::vector<DeviceTensorFloatND>>, std::optional<std::vector<DeviceTensorFloatND>>,  // jacobian residual product
    std::optional<std::vector<DeviceTensorFloatND>>, std::optional<std::vector<DeviceTensorFloatND>>,  // hessian diagonal
    DeviceVector<float>  // accum_weight
> _rasterize_to_pixels_eval3d_bwd_tensor(
    // Gaussian parameters
    int64_t num_splats,
    std::vector<DeviceTensorFloatND> splats_w,
    std::vector<DeviceTensorFloatND> splats_s,
    DeviceVector<int32_t> gaussian_ids,
    TorchTensorView viewmats,  // [..., C, 4, 4]
    TorchTensorView intrins,  // [..., C, 4], fx, fy, cx, cy
    const ssplat::CameraModelType camera_model,
    const TorchTensorView dist_coeffs,
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
    DeviceTensor3D<float> loss_map,  // [..., image_height, image_width, 1]
    DeviceTensor3D<float> accum_weight_map,  // [I, H, W]
    // gradients of outputs
    RenderOutput::TensorTuple v_render_outputs,
    const DeviceTensor3D<float> v_render_Ts, // [..., image_height, image_width, 1]
    std::optional<RenderOutput::TensorTuple> v_distortion_outputs_tuple,
    std::optional<std::vector<DeviceTensorFloatND>> v_splats_w,
    std::optional<std::vector<DeviceTensorFloatND>> v_splats_s,
    bool need_viewmat_grad
) {
    if (!v_splats_w.has_value())
        v_splats_w = SplatPrimitive::WorldBuffer::zeros_pool(splats_w, "raster_bwd.v_world");
    if (!v_splats_s.has_value())
        v_splats_s = SplatPrimitive::ScreenBuffer::zeros_pool(splats_s, "raster_bwd.v_screen");

    DeviceTensor2D<float> v_viewmats_buf;
    if (need_viewmat_grad) {
        auto& shape = std::get<2>(viewmats);
        int64_t total = 1;
        for (auto s : shape) total *= s;
        v_viewmats_buf.resize("raster_bwd.v_viewmats", total, 1);
        v_viewmats_buf.zero();
    }

    RenderOutput::Tensor render_outputs = render_outputs_tuple;
    RenderOutput::Tensor render2_outputs;
    RenderOutput::Tensor v_distortion_outputs;
    if (output_distortion) {
        render2_outputs = render2_outputs_tuple.value();
        v_distortion_outputs = v_distortion_outputs_tuple.value();
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
    DeviceVector<float> o_accum_weight;
    if (output_accum_weight) {
        o_accum_weight.resize("raster_bwd.accum_weight", num_splats);
        o_accum_weight.zero();
    }
    DeviceTensor3D<float> o_accum_weight_3d;
    if (output_accum_weight) {
        TorchTensorView tv((uint64_t)o_accum_weight.data_ptr(), 4, {num_splats, 1, 1, 1});
        o_accum_weight_3d = DeviceTensor3D<float>(tv);
    }

    launch_rasterize_to_pixels_eval3d_bwd_kernel<SplatPrimitive, output_distortion, output_hessian_diagonal, output_accum_weight>(
        num_splats,
        splats_w, splats_s, gaussian_ids,
        viewmats, intrins, camera_model, dist_coeffs,
        image_width, image_height, tile_offsets, flatten_ids,
        render_Ts, last_ids, render_outputs,
        render2_outputs, loss_map, accum_weight_map,
        v_render_outputs, v_render_Ts,
        v_distortion_outputs,
        v_splats_w.value(), v_splats_s.value(),
        vr_splats_w, vr_splats_s, h_splats_w, h_splats_s,
        o_accum_weight_3d,
        v_viewmats_buf
    );

    if (output_hessian_diagonal)
        return std::make_tuple(
            v_splats_w.value(), v_splats_s.value(), v_viewmats_buf,
            (std::optional<std::vector<DeviceTensorFloatND>>)vr_splats_w,
            (std::optional<std::vector<DeviceTensorFloatND>>)vr_splats_s,
            (std::optional<std::vector<DeviceTensorFloatND>>)h_splats_w,
            (std::optional<std::vector<DeviceTensorFloatND>>)h_splats_s,
            o_accum_weight);
    return std::make_tuple(
        v_splats_w.value(), v_splats_s.value(), v_viewmats_buf,
        (std::optional<std::vector<DeviceTensorFloatND>>)std::nullopt,
        (std::optional<std::vector<DeviceTensorFloatND>>)std::nullopt,
        (std::optional<std::vector<DeviceTensorFloatND>>)std::nullopt,
        (std::optional<std::vector<DeviceTensorFloatND>>)std::nullopt,
        o_accum_weight);
}


template<typename SplatPrimitive, bool output_distortion, bool output_accum_weight>
inline std::tuple<
    std::vector<DeviceTensorFloatND>, std::vector<DeviceTensorFloatND>,
    DeviceTensor2D<float>,  // v_viewmats
    DeviceVector<float>  // accum_weight
> rasterize_to_pixels_eval3d_bwd_tensor(
    // Gaussian parameters
    int64_t num_splats,
    std::vector<DeviceTensorFloatND> splats_w,
    std::vector<DeviceTensorFloatND> splats_s,
    DeviceVector<int32_t> gaussian_ids,
    TorchTensorView viewmats,  // [..., C, 4, 4]
    TorchTensorView intrins,  // [..., C, 4], fx, fy, cx, cy
    const ssplat::CameraModelType camera_model,
    const TorchTensorView dist_coeffs,
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
    DeviceTensor3D<float> loss_map,  // [..., image_height, image_width, 1]
    DeviceTensor3D<float> accum_weight_map,  // [I, H, W]
    // gradients of outputs
    RenderOutput::TensorTuple v_render_outputs,
    const DeviceTensor3D<float> v_render_Ts, // [..., image_height, image_width, 1]
    std::optional<RenderOutput::TensorTuple> v_distortion_outputs,
    std::optional<std::vector<DeviceTensorFloatND>> v_splats_w,
    std::optional<std::vector<DeviceTensorFloatND>> v_splats_s,
    bool need_viewmat_grad
) {
    auto [v_splats_w_1, v_splats_s_1, v_viewmats, vr_splats_w, vr_splats_s, h_splats_w, h_splats_s, accum_weight] =
        _rasterize_to_pixels_eval3d_bwd_tensor<SplatPrimitive, output_distortion, false, output_accum_weight>
    (
        num_splats, splats_w, splats_s, gaussian_ids,
        viewmats, intrins, camera_model, dist_coeffs,
        image_width, image_height, tile_offsets, flatten_ids,
        render_Ts, last_ids, render_outputs, render2_outputs, loss_map, accum_weight_map,
        v_render_outputs, v_render_Ts, v_distortion_outputs, v_splats_w, v_splats_s,
        need_viewmat_grad
    );
    return std::make_tuple(v_splats_w_1, v_splats_s_1, v_viewmats, accum_weight);
}


// ================
// Vanilla3DGUT
// ================

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    std::vector<DeviceTensorFloatND>, std::vector<DeviceTensorFloatND>,  // gradient
    DeviceTensor2D<float>,  // v_viewmats
    DeviceVector<float>  // accum_weight
> rasterize_to_pixels_3dgut_bwd(
    // Gaussian parameters
    int64_t num_splats,
    std::vector<DeviceTensorFloatND> splats_w,
    std::vector<DeviceTensorFloatND> splats_s,
    DeviceVector<int32_t> gaussian_ids,
    TorchTensorView viewmats,  // [..., C, 4, 4]
    TorchTensorView intrins,  // [..., C, 4], fx, fy, cx, cy
    const std::string camera_model,
    const TorchTensorView dist_coeffs,
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
    DeviceTensor3D<float> loss_map,  // [..., image_height, image_width, 1]
    DeviceTensor3D<float> accum_weight_map,  // [I, H, W]
    // gradients of outputs
    RenderOutput::TensorTuple v_render_outputs,
    const DeviceTensor3D<float> v_render_Ts, // [..., image_height, image_width, 1]
    std::optional<RenderOutput::TensorTuple> v_distortion_outputs,
    std::optional<std::vector<DeviceTensorFloatND>> v_splats_w,
    std::optional<std::vector<DeviceTensorFloatND>> v_splats_s,
    bool need_viewmat_grad
) {
    using Fn = decltype(&rasterize_to_pixels_eval3d_bwd_tensor<Vanilla3DGUT<0>, false, false>);
    static constexpr Fn funcs[2][2] = { {
        rasterize_to_pixels_eval3d_bwd_tensor<Vanilla3DGUT<0>, false, false>,
        rasterize_to_pixels_eval3d_bwd_tensor<Vanilla3DGUT<0>, false, true>,
    }, {
        rasterize_to_pixels_eval3d_bwd_tensor<Vanilla3DGUT<0>, true, false>,
        rasterize_to_pixels_eval3d_bwd_tensor<Vanilla3DGUT<0>, true, true>,
    } };
    return funcs[v_distortion_outputs.has_value()][accum_weight_map.data_ptr() != nullptr](
        num_splats, splats_w, splats_s, gaussian_ids,
        viewmats, intrins, cmt(camera_model), dist_coeffs,
        image_width, image_height, tile_offsets, flatten_ids,
        render_Ts, last_ids, render_outputs, render2_outputs, loss_map, accum_weight_map,
        v_render_outputs, v_render_Ts, v_distortion_outputs, v_splats_w, v_splats_s,
        need_viewmat_grad
    );
}

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    std::vector<DeviceTensorFloatND>, std::vector<DeviceTensorFloatND>,  // gradient
    DeviceTensor2D<float>,  // v_viewmats
    std::optional<std::vector<DeviceTensorFloatND>>, std::optional<std::vector<DeviceTensorFloatND>>,  // jacobian residual product
    std::optional<std::vector<DeviceTensorFloatND>>, std::optional<std::vector<DeviceTensorFloatND>>,  // hessian diagonal
    DeviceVector<float>  // accum_weight
> rasterize_to_pixels_3dgut_bwd_with_hessian_diagonal(
    // Gaussian parameters
    int64_t num_splats,
    std::vector<DeviceTensorFloatND> splats_w,
    std::vector<DeviceTensorFloatND> splats_s,
    DeviceVector<int32_t> gaussian_ids,
    TorchTensorView viewmats,  // [..., C, 4, 4]
    TorchTensorView intrins,  // [..., C, 4], fx, fy, cx, cy
    const std::string camera_model,
    const TorchTensorView dist_coeffs,
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
    DeviceTensor3D<float> loss_map,  // [..., image_height, image_width, 1]
    DeviceTensor3D<float> accum_weight_map,  // [I, H, W]
    // gradients of outputs
    RenderOutput::TensorTuple v_render_outputs,
    const DeviceTensor3D<float> v_render_Ts, // [..., image_height, image_width, 1]
    std::optional<RenderOutput::TensorTuple> v_distortion_outputs,
    std::optional<std::vector<DeviceTensorFloatND>> v_splats_w,
    std::optional<std::vector<DeviceTensorFloatND>> v_splats_s,
    bool need_viewmat_grad
) {
    using Fn = decltype(&_rasterize_to_pixels_eval3d_bwd_tensor<Vanilla3DGUT<0>, false, true, false>);
    static constexpr Fn funcs[2][2] = { {
        _rasterize_to_pixels_eval3d_bwd_tensor<Vanilla3DGUT<0>, false, true, false>,
        _rasterize_to_pixels_eval3d_bwd_tensor<Vanilla3DGUT<0>, false, true, true>,
    }, {
        _rasterize_to_pixels_eval3d_bwd_tensor<Vanilla3DGUT<0>, true, true, false>,
        _rasterize_to_pixels_eval3d_bwd_tensor<Vanilla3DGUT<0>, true, true, true>,
    } };
    return funcs[v_distortion_outputs.has_value()][accum_weight_map.data_ptr() != nullptr](
        num_splats, splats_w, splats_s, gaussian_ids,
        viewmats, intrins, cmt(camera_model), dist_coeffs,
        image_width, image_height, tile_offsets, flatten_ids,
        render_Ts, last_ids, render_outputs, render2_outputs, loss_map, accum_weight_map,
        v_render_outputs, v_render_Ts, v_distortion_outputs, v_splats_w, v_splats_s,
        need_viewmat_grad
    );
}


