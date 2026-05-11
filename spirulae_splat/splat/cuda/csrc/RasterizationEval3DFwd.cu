#include "RasterizationEval3DFwd.cuh"

#include <gsplat/Utils.cuh>

#include "types.cuh"
#include "common.cuh"

#include <c10/cuda/CUDAStream.h>


template<
    typename SplatPrimitive,
    ssplat::CameraModelType camera_model,
    bool output_distortion
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
    RenderOutput::Buffer render_distortions // [I, image_height, image_width, ...]
);


template <typename SplatPrimitive, bool output_distortion>
inline void launch_rasterize_to_pixels_eval3d_fwd_kernel(
    // Gaussian parameters
    typename SplatPrimitive::WorldBuffer splats_w,
    typename SplatPrimitive::ScreenBuffer splats_s,
    std::optional<at::Tensor> gaussian_ids,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const ssplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    // intersections
    const at::Tensor tile_offsets, // [..., tile_height, tile_width]
    const at::Tensor flatten_ids,  // [n_isects]
    // outputs
    RenderOutput::Tensor renders,
    at::Tensor transmittances,  // [..., image_height, image_width]
    at::Tensor last_ids, // [..., image_height, image_width]
    RenderOutput::Tensor *renders2,
    RenderOutput::Tensor *distortions
) {
    uint32_t N = gaussian_ids.has_value() ? 0 : splats_w.size();
    uint32_t I = transmittances.numel() / (image_height * image_width); // number of images
    uint32_t tile_height = tile_offsets.size(-2);
    uint32_t tile_width = tile_offsets.size(-1);
    uint32_t n_isects = flatten_ids.size(0);

    #define _LAUNCH_ARGS ( \
            (cudaStream_t)at::cuda::getCurrentCUDAStream(), I, N, n_isects, \
            gaussian_ids.has_value() ? (uint32_t*)gaussian_ids.value().data_ptr<int32_t>() : nullptr, \
            splats_w, splats_s, \
            viewmats.data_ptr<float>(), (float4*)intrins.data_ptr<float>(), dist_coeffs, \
            image_width, image_height, tile_width, tile_height, \
            tile_offsets.data_ptr<int32_t>(), flatten_ids.data_ptr<int32_t>(), \
            renders, transmittances.data_ptr<float>(), last_ids.data_ptr<int32_t>(), \
            output_distortion ? renders2->buffer() : RenderOutput::Buffer(), \
            output_distortion ? distortions->buffer() : RenderOutput::Buffer() \
        )

    if (camera_model == ssplat::CameraModelType::PINHOLE)
        rasterize_to_pixels_eval3d_fwd_kernel_wrapper<SplatPrimitive,
            ssplat::CameraModelType::PINHOLE, output_distortion> _LAUNCH_ARGS;
    else if (camera_model == ssplat::CameraModelType::FISHEYE)
        rasterize_to_pixels_eval3d_fwd_kernel_wrapper<SplatPrimitive,
            ssplat::CameraModelType::FISHEYE, output_distortion> _LAUNCH_ARGS;
    else
        throw std::runtime_error("Unsupported camera model");
    CHECK_DEVICE_ERROR(cudaGetLastError());

    #undef _LAUNCH_ARGS
}


template <typename SplatPrimitive, bool output_distortion>
inline std::tuple<
    RenderOutput::TensorTuple,
    at::Tensor,
    at::Tensor,
    std::optional<RenderOutput::TensorTuple>,
    std::optional<RenderOutput::TensorTuple>
> rasterize_to_pixels_eval3d_fwd_tensor(
    // Gaussian parameters
    TensorList splats_w,
    TensorList splats_s,
    std::optional<at::Tensor> gaussian_ids,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const ssplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    // intersections
    const at::Tensor tile_offsets, // [..., tile_height, tile_width]
    const at::Tensor flatten_ids   // [n_isects]
) {
    DEVICE_GUARD(tile_offsets);
    CHECK_INPUT(tile_offsets);
    CHECK_INPUT(flatten_ids);
    CHECK_INPUT(viewmats);
    CHECK_INPUT(intrins);
    
    at::DimVector image_dims(tile_offsets.sizes().slice(0, tile_offsets.dim() - 2));

    at::DimVector renders_dims(image_dims);
    renders_dims.append({image_height, image_width});
    RenderOutput::Tensor renders =
        RenderOutput::Tensor::empty<SplatPrimitive::pixelType>(renders_dims);

    std::optional<RenderOutput::Tensor> renders2 = std::nullopt;
    std::optional<RenderOutput::Tensor> distortions = std::nullopt;
    if (output_distortion) {
        renders2 = RenderOutput::Tensor::empty<SplatPrimitive::pixelType>(renders_dims);
        distortions = RenderOutput::Tensor::empty<SplatPrimitive::pixelType>(renders_dims);
    }

    at::DimVector transmittance_dims(image_dims);
    transmittance_dims.append({image_height, image_width, 1});
    at::Tensor transmittances = at::empty(transmittance_dims, kTensorOptionF32());

    at::DimVector last_ids_dims(image_dims);
    last_ids_dims.append({image_height, image_width});
    at::Tensor last_ids = at::empty(last_ids_dims, kTensorOptionI32());

    launch_rasterize_to_pixels_eval3d_fwd_kernel<SplatPrimitive, output_distortion>(
        splats_w, splats_s, gaussian_ids,
        viewmats, intrins, camera_model, dist_coeffs,
        image_width, image_height, tile_offsets, flatten_ids,
        renders, transmittances, last_ids,
        output_distortion ? &renders2.value() : nullptr,
        output_distortion ? &distortions.value() : nullptr
    );

    if (output_distortion)
        return std::make_tuple(
            renders.tuple(), transmittances, last_ids,
            (std::optional<RenderOutput::TensorTuple>)renders2.value().tuple(),
            (std::optional<RenderOutput::TensorTuple>)distortions.value().tuple()
        );
    return std::make_tuple(
        renders.tuple(), transmittances, last_ids,
        (std::optional<RenderOutput::TensorTuple>)std::nullopt,
        (std::optional<RenderOutput::TensorTuple>)std::nullopt
    );
}



// ================
// Vanilla3DGUT
// ================

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    RenderOutput::TensorTuple,
    at::Tensor,
    at::Tensor,
    std::optional<RenderOutput::TensorTuple>,
    std::optional<RenderOutput::TensorTuple>
> rasterize_to_pixels_3dgut_fwd(
    // Gaussian parameters
    TensorList splats_w,
    TensorList splats_s,
    std::optional<at::Tensor> gaussian_ids,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    // intersections
    const at::Tensor tile_offsets, // [..., tile_height, tile_width]
    const at::Tensor flatten_ids,   // [n_isects]
    bool output_distortion
) {
    return (output_distortion ?
        rasterize_to_pixels_eval3d_fwd_tensor<Vanilla3DGUT, true> :
        rasterize_to_pixels_eval3d_fwd_tensor<Vanilla3DGUT, false>
    )(
        splats_w, splats_s, gaussian_ids,
        viewmats, intrins, cmt(camera_model), dist_coeffs,
        image_width, image_height,
        tile_offsets, flatten_ids
    );
}



// // ================
// // SphericalVoronoi3DGUT
// // ================

// /*[AutoHeaderGeneratorExport]*/
// std::tuple<
//     RenderOutput::TensorTuple,
//     at::Tensor,
//     at::Tensor,
//     std::optional<RenderOutput::TensorTuple>,
//     std::optional<RenderOutput::TensorTuple>,
//     std::optional<at::Tensor>
// > rasterize_to_pixels_3dgut_sv_fwd(
//     // Gaussian parameters
//     SphericalVoronoi3DGUT_Default::Screen::TensorTuple splats_tuple,
//     std::optional<at::Tensor> gaussian_ids,
//     const at::Tensor viewmats,  // [..., C, 4, 4]
//     const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
//     const std::string camera_model,
//     const CameraDistortionCoeffsTensor dist_coeffs,
//     // image size
//     const uint32_t image_width,
//     const uint32_t image_height,
//     // intersections
//     const at::Tensor tile_offsets, // [..., tile_height, tile_width]
//     const at::Tensor flatten_ids,   // [n_isects]
//     bool output_distortion
// ) {
//     if (output_distortion)
//         // return rasterize_to_pixels_eval3d_fwd_tensor<SphericalVoronoi3DGUT_Default, true, false>(
//         return rasterize_to_pixels_eval3d_fwd_tensor<Vanilla3DGUT, true, false>(
//             splats_tuple, gaussian_ids,
//             viewmats, intrins, cmt(camera_model), dist_coeffs,
//             image_width, image_height,
//             tile_offsets, flatten_ids
//         );
//     // return rasterize_to_pixels_eval3d_fwd_tensor<SphericalVoronoi3DGUT_Default, false, false>(
//     return rasterize_to_pixels_eval3d_fwd_tensor<Vanilla3DGUT, false, false>(
//         splats_tuple, gaussian_ids,
//         viewmats, intrins, cmt(camera_model), dist_coeffs,
//         image_width, image_height,
//         tile_offsets, flatten_ids
//     );
// }


// // ================
// // VoxelPrimitive
// // ================

// /*[AutoHeaderGeneratorExport]*/
// std::tuple<
//     RenderOutput::TensorTuple,
//     at::Tensor,
//     at::Tensor,
//     std::optional<RenderOutput::TensorTuple>,
//     std::optional<RenderOutput::TensorTuple>,
//     std::optional<at::Tensor>
// > rasterize_to_pixels_voxel_eval3d_fwd(
//     // Gaussian parameters
//     VoxelPrimitive::Screen::TensorTuple splats_tuple,
//     std::optional<at::Tensor> gaussian_ids,
//     const at::Tensor viewmats,      // [..., C, 4, 4]
//     const at::Tensor intrins,       // [..., C, 4], fx, fy, cx, cy
//     const std::string camera_model,
//     const CameraDistortionCoeffsTensor dist_coeffs,
//     // image size
//     const uint32_t image_width,
//     const uint32_t image_height,
//     // intersections
//     const at::Tensor tile_offsets, // [..., tile_height, tile_width]
//     const at::Tensor flatten_ids   // [n_isects]
// ) {
//     return rasterize_to_pixels_eval3d_fwd_tensor<VoxelPrimitive, true, true>(
//         splats_tuple, gaussian_ids,
//         viewmats, intrins, cmt(camera_model), dist_coeffs,
//         image_width, image_height,
//         tile_offsets, flatten_ids
//     );
// }
