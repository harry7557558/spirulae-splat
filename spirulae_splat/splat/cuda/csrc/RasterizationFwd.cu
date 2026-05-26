#include "RasterizationFwd.cuh"

#include <gsplat/Utils.cuh>

#include "types.cuh"
#include "common.cuh"


#include <Tensor.h>


template <typename SplatPrimitive, bool output_distortion>
void rasterize_to_pixels_fwd_kernel_wrapper(
    cudaStream_t stream,
    const uint32_t I,
    const uint32_t N,
    const uint32_t n_isects,
    const uint32_t *__restrict__ gaussian_ids,  // [nnz] optional, for packed mode
    const typename SplatPrimitive::WorldBuffer splat_wbuffer,
    const typename SplatPrimitive::ScreenBuffer splat_sbuffer,
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_width,
    const uint32_t tile_height,
    const int32_t *__restrict__ tile_offsets, // [I, tile_height, tile_width]
    const int32_t *__restrict__ flatten_ids,  // [n_isects]
    RenderOutput::Buffer render_colors, // [I, image_height, image_width, 3]
    float *__restrict__ render_Ts, // [I, image_height, image_width, 1]
    int32_t *__restrict__ last_ids,        // [I, image_height, image_width]
    RenderOutput::Buffer render_colors2, // [I, image_height, image_width, ...]
    RenderOutput::Buffer render_distortions // [I, image_height, image_width, ...]
);

template <typename SplatPrimitive, bool output_distortion>
inline void launch_rasterize_to_pixels_fwd_kernel(
    // Gaussian parameters
    typename SplatPrimitive::WorldBuffer splats_w,
    typename SplatPrimitive::ScreenBuffer splats_s,
    DeviceVector<int32_t> gaussian_ids,
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
    RenderOutput::Tensor distortions
) {
    uint32_t N = gaussian_ids.data_ptr() ? 0 : splats_w.size(); // number of gaussians
    uint32_t I = transmittances.size<0>();  // number of images
    uint32_t tile_height = tile_offsets.size<1>();
    uint32_t tile_width = tile_offsets.size<2>();
    uint32_t n_isects = flatten_ids.size();

    rasterize_to_pixels_fwd_kernel_wrapper<SplatPrimitive, output_distortion>(
        (cudaStream_t)0,
        I, N, n_isects,
        (uint32_t*)gaussian_ids.data_ptr(),
        splats_w, splats_s,
        image_width,
        image_height,
        tile_width,
        tile_height,
        tile_offsets.data_ptr(),
        flatten_ids.data_ptr(),
        renders,
        transmittances.data_ptr(),
        last_ids.data_ptr(),
        output_distortion ? renders2.buffer() : RenderOutput::Buffer(),
        output_distortion ? distortions.buffer() : RenderOutput::Buffer()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


template <typename SplatPrimitive, bool output_distortion>
inline std::tuple<
    RenderOutput::TensorTuple,  // renders
    DeviceTensor3D<float>,  // transmittances
    DeviceTensor3D<int32_t>,  // last_ids
    RenderOutput::TensorTuple,  // renders2, optional
    RenderOutput::TensorTuple  // distortions, optional
> rasterize_to_pixels_fwd_tensor(
    // Gaussian parameters
    std::vector<DeviceTensorFloatND> splats_w,
    std::vector<DeviceTensorFloatND> splats_s,
    DeviceVector<int32_t> gaussian_ids,
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

    launch_rasterize_to_pixels_fwd_kernel<SplatPrimitive, output_distortion>(
        splats_w, splats_s, gaussian_ids,
        image_width,
        image_height,
        tile_offsets,
        flatten_ids,
        renders,
        render_Ts,
        render_last_ids,
        renders2,
        distortions
    );

    return std::make_tuple(
        renders, render_Ts, render_last_ids,
        renders2, distortions
    );
}




// ================
// Vanilla3DGS
// ================

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    RenderOutput::TensorTuple,  // renders
    DeviceTensor3D<float>,  // transmittances
    DeviceTensor3D<int32_t>,  // last_ids
    RenderOutput::TensorTuple,  // renders2, optional
    RenderOutput::TensorTuple  // distortions, optional
> rasterize_to_pixels_3dgs_fwd(
    // Gaussian parameters
    std::vector<DeviceTensorFloatND> splats_w,
    std::vector<DeviceTensorFloatND> splats_s,
    DeviceVector<int32_t> gaussian_ids,
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    // intersections
    const DeviceTensor3D<int32_t> tile_offsets,
    const DeviceVector<int32_t> flatten_ids,
    bool output_distortion
) {
    return (output_distortion ?
        rasterize_to_pixels_fwd_tensor<Vanilla3DGS<0>, true> :
        rasterize_to_pixels_fwd_tensor<Vanilla3DGS<0>, false>
    )(
        splats_w, splats_s, gaussian_ids,
        image_width, image_height,
        tile_offsets, flatten_ids
    );
}



// ================
// MipSplatting
// ================

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    RenderOutput::TensorTuple,  // renders
    DeviceTensor3D<float>,  // transmittances
    DeviceTensor3D<int32_t>,  // last_ids
    RenderOutput::TensorTuple,  // renders2, optional
    RenderOutput::TensorTuple  // distortions, optional
> rasterize_to_pixels_mip_fwd(
    // Gaussian parameters
    std::vector<DeviceTensorFloatND> splats_w,
    std::vector<DeviceTensorFloatND> splats_s,
    DeviceVector<int32_t> gaussian_ids,
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    // intersections
    const DeviceTensor3D<int32_t> tile_offsets,
    const DeviceVector<int32_t> flatten_ids,
    bool output_distortion
) {
    return (output_distortion ?
        rasterize_to_pixels_fwd_tensor<Vanilla3DGS<0>, true> :
        rasterize_to_pixels_fwd_tensor<Vanilla3DGS<0>, false>
    )(
        splats_w, splats_s, gaussian_ids,
        image_width, image_height,
        tile_offsets, flatten_ids
    );
}
