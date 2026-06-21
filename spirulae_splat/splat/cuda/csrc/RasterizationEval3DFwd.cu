#include "RasterizationEval3DFwd.cuh"

#include <Common.cuh>


#include <Tensor.h>


template<
    typename SplatPrimitive,
    CameraModelType camera_model,
    bool output_distortion,
    bool output_median
>
void rasterize_to_pixels_eval3d_fwd_kernel_wrapper(
    cudaStream_t stream,
    const uint32_t I,
    const uint32_t N,
    const uint32_t n_isects,
    const uint32_t *__restrict__ gaussian_ids,  // [nnz] optional, for packed mode
    const typename SplatPrimitive::WorldBuffer splat_wbuffer,
    const typename SplatPrimitive::ScreenBuffer splat_sbuffer,
    const float *__restrict__ viewmats, // [B, C, 4, 4]
    const float4 *__restrict__ intrins,  // [B, C, 4], fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const float4 *__restrict__ aabb,  // [..., N] projected 2D AABB
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_width,
    const uint32_t tile_height,
    const int32_t *__restrict__ tile_offsets, // [I, tile_height, tile_width]
    const int32_t *__restrict__ flatten_ids,  // [n_isects]
    RenderOutput::Buffer render_colors, // [I, image_height, image_width, ...]
    float *__restrict__ render_Ts, // [I, image_height, image_width, 1]
    int32_t *__restrict__ last_ids, // [I, image_height, image_width]
    RenderOutput::Buffer render_colors2, // [I, image_height, image_width, ...]
    RenderOutput::Buffer render_distortions, // [I, image_height, image_width, ...]
    float *__restrict__ render_median, // [I, image_height, image_width, 1], optional
    float2 *__restrict__ render_median_anchor // [I, image_height, image_width, 2], optional
);


template <typename SplatPrimitive, bool output_distortion, bool output_median>
inline void launch_rasterize_to_pixels_eval3d_fwd_kernel(
    // Gaussian parameters
    int64_t num_splats,  // = cur_num_splats; non-packed projection layout stride
    typename SplatPrimitive::WorldBuffer splats_w,
    typename SplatPrimitive::ScreenBuffer splats_s,
    DeviceVector<int32_t> gaussian_ids,
    TorchTensorView viewmats,  // [..., C, 4, 4]
    TorchTensorView intrins,  // [..., C, 4], fx, fy, cx, cy
    const CameraModelType camera_model,
    const TorchTensorView dist_coeffs,
    DeviceTensor2D<float4> aabb,  // [..., N] projected 2D AABB, for sub-tile culling
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    // intersections
    const DeviceTensor3D<int32_t> tile_offsets, // [I, tile_height, tile_width]
    const DeviceVector<int32_t> flatten_ids,    // [n_isects]
    // outputs
    RenderOutput::Tensor renders,
    DeviceTensor3D<float> transmittances,
    DeviceTensor3D<int32_t> last_ids,
    RenderOutput::Tensor renders2,
    RenderOutput::Tensor distortions,
    DeviceTensor3D<float> render_median,
    DeviceTensor3D<float2> render_median_anchor
) {
    // splats_w.size() returns max_num_splats (pre-allocated); the projection
    // layout uses cur_num_splats per camera stride. See RasterizationBwd.cu.
    uint32_t N = gaussian_ids.data_ptr() ? 0 : (uint32_t)num_splats;
    uint32_t I = transmittances.size<0>();  // number of images
    uint32_t tile_height = tile_offsets.size<1>();
    uint32_t tile_width = tile_offsets.size<2>();
    uint32_t n_isects = flatten_ids.size();

    const float* viewmats_ptr = (const float*)std::get<0>(viewmats);
    const float4* intrins_ptr = (const float4*)std::get<0>(intrins);
    const float4* aabb_ptr = (const float4*)aabb.data_ptr();

    #define _LAUNCH_ARGS ( \
            (cudaStream_t)0, I, N, n_isects, \
            (uint32_t*)gaussian_ids.data_ptr(), \
            splats_w, splats_s, \
            viewmats_ptr, intrins_ptr, dist_coeffs, aabb_ptr, \
            image_width, image_height, tile_width, tile_height, \
            tile_offsets.data_ptr(), flatten_ids.data_ptr(), \
            renders, transmittances.data_ptr(), last_ids.data_ptr(), \
            output_distortion ? renders2.buffer() : RenderOutput::Buffer(), \
            output_distortion ? distortions.buffer() : RenderOutput::Buffer(), \
            output_median ? render_median.data_ptr() : nullptr, \
            output_median ? render_median_anchor.data_ptr() : nullptr \
        )

    if (camera_model == CameraModelType::PINHOLE)
        rasterize_to_pixels_eval3d_fwd_kernel_wrapper<SplatPrimitive,
            CameraModelType::PINHOLE, output_distortion, output_median> _LAUNCH_ARGS;
    else if (camera_model == CameraModelType::FISHEYE)
        rasterize_to_pixels_eval3d_fwd_kernel_wrapper<SplatPrimitive,
            CameraModelType::FISHEYE, output_distortion, output_median> _LAUNCH_ARGS;
    else if (camera_model == CameraModelType::EQUISOLID)
        rasterize_to_pixels_eval3d_fwd_kernel_wrapper<SplatPrimitive,
            CameraModelType::EQUISOLID, output_distortion, output_median> _LAUNCH_ARGS;
    else
        throw std::runtime_error("Unsupported camera model");
    CHECK_DEVICE_ERROR(cudaGetLastError());

    #undef _LAUNCH_ARGS
}


template <typename SplatPrimitive, bool output_distortion, bool output_median>
inline std::tuple<
    RenderOutput::TensorTuple,  // renders
    DeviceTensor3D<float>,  // transmittances
    DeviceTensor3D<int32_t>,  // last_ids
    RenderOutput::TensorTuple,  // renders2, optional
    RenderOutput::TensorTuple,  // distortions, optional
    DeviceTensor3D<float>,  // median depth, optional
    DeviceTensor3D<float2>  // median near-anchor (z1, v1), optional
> rasterize_to_pixels_eval3d_fwd_tensor(
    // Gaussian parameters
    int64_t num_splats,  // = cur_num_splats
    std::vector<DeviceTensorFloatND> splats_w,
    std::vector<DeviceTensorFloatND> splats_s,
    DeviceVector<int32_t> gaussian_ids,
    TorchTensorView viewmats,  // [..., C, 4, 4]
    TorchTensorView intrins,  // [..., C, 4], fx, fy, cx, cy
    const CameraModelType camera_model,
    const TorchTensorView dist_coeffs,
    DeviceTensor2D<float4> aabb,  // [..., N] projected 2D AABB, for sub-tile culling
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    // intersections
    const DeviceTensor3D<int32_t> tile_offsets, // [I, tile_height, tile_width]
    const DeviceVector<int32_t> flatten_ids     // [n_isects]
) {
    int64_t batch = tile_offsets.size<0>();

    RenderOutput::TensorTuple renders, renders2, distortions;
    RenderOutput::resize<SplatPrimitive::pixelType>(
        renders, batch, image_height, image_width, "renders");
    if (output_distortion) {
        RenderOutput::resize<SplatPrimitive::pixelType>(
            renders2, batch, image_height, image_width, "renders2");
        RenderOutput::resize<SplatPrimitive::pixelType>(
            distortions, batch, image_height, image_width, "distortions");
    }

    DeviceTensor3D<float> render_Ts;
    DeviceTensor3D<int32_t> render_last_ids;
    render_Ts.resize("render.Ts", batch, image_height, image_width);
    render_last_ids.resize("render.last_ids", batch, image_height, image_width);

    DeviceTensor3D<float> render_median;
    DeviceTensor3D<float2> render_median_anchor;
    if (output_median) {
        render_median.resize("render.median", batch, image_height, image_width);
        render_median_anchor.resize("render.median_anchor", batch, image_height, image_width);
    }

    launch_rasterize_to_pixels_eval3d_fwd_kernel<SplatPrimitive, output_distortion, output_median>(
        num_splats,
        splats_w, splats_s, gaussian_ids,
        viewmats, intrins, camera_model, dist_coeffs, aabb,
        image_width, image_height, tile_offsets, flatten_ids,
        renders, render_Ts, render_last_ids, renders2, distortions,
        render_median, render_median_anchor
    );

    return std::make_tuple(
        renders, render_Ts, render_last_ids,
        renders2, distortions, render_median, render_median_anchor
    );
}



// ================
// Vanilla3DGUT
// ================

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    RenderOutput::TensorTuple,  // renders
    DeviceTensor3D<float>,  // transmittances
    DeviceTensor3D<int32_t>,  // last_ids
    RenderOutput::TensorTuple,  // renders2, optional
    RenderOutput::TensorTuple,  // distortions, optional
    DeviceTensor3D<float>,  // median depth, optional
    DeviceTensor3D<float2>  // median near-anchor (z1, v1), optional
> rasterize_to_pixels_3dgut_fwd(
    // Gaussian parameters
    int64_t num_splats,
    std::vector<DeviceTensorFloatND> splats_w,
    std::vector<DeviceTensorFloatND> splats_s,
    DeviceVector<int32_t> gaussian_ids,
    TorchTensorView viewmats,  // [..., C, 4, 4]
    TorchTensorView intrins,  // [..., C, 4], fx, fy, cx, cy
    const std::string camera_model,
    const TorchTensorView dist_coeffs,
    DeviceTensor2D<float4> aabb,  // [..., N] projected 2D AABB, for sub-tile culling
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    // intersections
    const DeviceTensor3D<int32_t> tile_offsets, // [I, tile_height, tile_width]
    const DeviceVector<int32_t> flatten_ids,    // [n_isects]
    bool output_distortion,
    bool output_median
) {
    auto dispatch = output_distortion ?
        (output_median ?
            rasterize_to_pixels_eval3d_fwd_tensor<Vanilla3DGUT<0>, true, true> :
            rasterize_to_pixels_eval3d_fwd_tensor<Vanilla3DGUT<0>, true, false>) :
        (output_median ?
            rasterize_to_pixels_eval3d_fwd_tensor<Vanilla3DGUT<0>, false, true> :
            rasterize_to_pixels_eval3d_fwd_tensor<Vanilla3DGUT<0>, false, false>);
    return dispatch(
        num_splats,
        splats_w, splats_s, gaussian_ids,
        viewmats, intrins, cmt(camera_model), dist_coeffs, aabb,
        image_width, image_height,
        tile_offsets, flatten_ids
    );
}


