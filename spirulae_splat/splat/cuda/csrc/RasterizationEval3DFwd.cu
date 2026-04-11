#include "RasterizationEval3DFwd.cuh"

#include <gsplat/Utils.cuh>

#include "types.cuh"
#include "common.cuh"

#include <c10/cuda/CUDAStream.h>


template<
    typename SplatPrimitive,
    ssplat::CameraModelType camera_model,
    bool output_distortion,
    bool output_accum_weight,
    bool output_max_blending
>
void rasterize_to_pixels_eval3d_fwd_kernel_wrapper(
    cudaStream_t stream,
    const uint32_t I,
    const uint32_t N,
    const uint32_t n_isects,
    const uint32_t *__restrict__ gaussian_ids,  // [nnz] optional, for packed mode
    const typename SplatPrimitive::Screen::Buffer splat_buffer,
    const float *__restrict__ viewmats, // [B, C, 4, 4]
    const float4 *__restrict__ intrins,  // [B, C, 4], fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const float3 *__restrict__ backgrounds, // [I, 3]
    const float *__restrict__ accum_weight_map,  // [B, C, image_width, image_height]
    const bool *__restrict__ max_blending_masks,  // [B, C, image_width, image_height]
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_width,
    const uint32_t tile_height,
    const int32_t *__restrict__ tile_offsets, // [I, tile_height, tile_width]
    const int32_t *__restrict__ flatten_ids,  // [n_isects]
    typename SplatPrimitive::RenderOutput::Buffer render_colors, // [I, image_height, image_width, ...]
    float *__restrict__ render_Ts, // [I, image_height, image_width, 1]
    int32_t *__restrict__ last_ids, // [I, image_height, image_width]
    typename SplatPrimitive::RenderOutput::Buffer render_colors2, // [I, image_height, image_width, ...]
    typename SplatPrimitive::RenderOutput::Buffer render_distortions, // [I, image_height, image_width, ...]
    float* __restrict__ out_accum_weight,
    float* __restrict__ out_max_blending
);


template <typename SplatPrimitive, bool output_distortion, bool output_accum_weight, bool output_max_blending>
inline void launch_rasterize_to_pixels_eval3d_fwd_kernel(
    // Gaussian parameters
    typename SplatPrimitive::Screen::Tensor splats,
    std::optional<at::Tensor> gaussian_ids,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const ssplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    const std::optional<at::Tensor> backgrounds, // [..., channels]
    const std::optional<at::Tensor> accum_weight_map,  // [..., C, image_width, image_height]
    const std::optional<at::Tensor> max_blending_masks,  // [..., C, image_width, image_height]
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    // intersections
    const at::Tensor tile_offsets, // [..., tile_height, tile_width]
    const at::Tensor flatten_ids,  // [n_isects]
    // outputs
    typename SplatPrimitive::RenderOutput::Tensor renders,
    at::Tensor transmittances,  // [..., image_height, image_width]
    at::Tensor last_ids, // [..., image_height, image_width]
    typename SplatPrimitive::RenderOutput::Tensor *renders2,
    typename SplatPrimitive::RenderOutput::Tensor *distortions,
    std::optional<at::Tensor>& out_accum_weight,
    std::optional<at::Tensor>& out_max_blending
) {
    bool packed = splats.isPacked();
    uint32_t N = packed ? 0 : splats.size(); // number of gaussians
    uint32_t I = transmittances.numel() / (image_height * image_width); // number of images
    uint32_t tile_height = tile_offsets.size(-2);
    uint32_t tile_width = tile_offsets.size(-1);
    uint32_t n_isects = flatten_ids.size(0);

    #define _LAUNCH_ARGS ( \
            (cudaStream_t)at::cuda::getCurrentCUDAStream(), I, N, n_isects, \
            gaussian_ids.has_value() ? (uint32_t*)gaussian_ids.value().data_ptr<int32_t>() : nullptr, splats.buffer(), \
            viewmats.data_ptr<float>(), (float4*)intrins.data_ptr<float>(), dist_coeffs, \
            backgrounds.has_value() ? (float3*)backgrounds.value().data_ptr<float>() : nullptr, \
            (output_accum_weight && accum_weight_map.has_value()) ? accum_weight_map.value().data_ptr<float>() : nullptr, \
            (output_max_blending && max_blending_masks.has_value()) ? max_blending_masks.value().data_ptr<bool>() : nullptr, \
            image_width, image_height, tile_width, tile_height, \
            tile_offsets.data_ptr<int32_t>(), flatten_ids.data_ptr<int32_t>(), \
            renders, transmittances.data_ptr<float>(), last_ids.data_ptr<int32_t>(), \
            output_distortion ? renders2->buffer() : typename SplatPrimitive::RenderOutput::Buffer(), \
            output_distortion ? distortions->buffer() : typename SplatPrimitive::RenderOutput::Buffer(), \
            out_accum_weight.has_value() ? out_accum_weight.value().data_ptr<float>() : nullptr, \
            out_max_blending.has_value() ? out_max_blending.value().data_ptr<float>() : nullptr \
        )

    if (camera_model == ssplat::CameraModelType::PINHOLE)
        rasterize_to_pixels_eval3d_fwd_kernel_wrapper<SplatPrimitive,
            ssplat::CameraModelType::PINHOLE, output_distortion, output_accum_weight, output_max_blending> _LAUNCH_ARGS;
    else if (camera_model == ssplat::CameraModelType::FISHEYE)
        rasterize_to_pixels_eval3d_fwd_kernel_wrapper<SplatPrimitive,
            ssplat::CameraModelType::FISHEYE, output_distortion, output_accum_weight, output_max_blending> _LAUNCH_ARGS;
    else
        throw std::runtime_error("Unsupported camera model");
    CHECK_DEVICE_ERROR(cudaGetLastError());

    #undef _LAUNCH_ARGS
}


template <typename SplatPrimitive, bool output_distortion, bool output_max_blending>
inline std::tuple<
    typename SplatPrimitive::RenderOutput::TensorTuple,
    at::Tensor,
    at::Tensor,
    std::optional<typename SplatPrimitive::RenderOutput::TensorTuple>,
    std::optional<typename SplatPrimitive::RenderOutput::TensorTuple>,
    std::optional<at::Tensor>,
    std::optional<at::Tensor>
> rasterize_to_pixels_eval3d_fwd_tensor(
    // Gaussian parameters
    typename SplatPrimitive::Screen::TensorTuple splats_tuple,
    std::optional<at::Tensor> gaussian_ids,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const ssplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    const std::optional<at::Tensor> backgrounds, // [..., channels]
    const std::optional<at::Tensor> accum_weight_map,  // [..., C, image_width, image_height]
    const std::optional<at::Tensor> max_blending_masks,       // [..., tile_height, tile_width]
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
    if (backgrounds.has_value())
        CHECK_INPUT(backgrounds.value());
    if (accum_weight_map.has_value())
        CHECK_INPUT(accum_weight_map.value());
    if (max_blending_masks.has_value())
        CHECK_INPUT(max_blending_masks.value());
    
    typename SplatPrimitive::Screen::Tensor splats(splats_tuple);

    auto opt = splats.options();
    at::DimVector image_dims(tile_offsets.sizes().slice(0, tile_offsets.dim() - 2));

    at::DimVector renders_dims(image_dims);
    renders_dims.append({image_height, image_width});
    typename SplatPrimitive::RenderOutput::Tensor renders =
        SplatPrimitive::RenderOutput::Tensor::empty(renders_dims, opt);

    std::optional<typename SplatPrimitive::RenderOutput::Tensor> renders2 = std::nullopt;
    std::optional<typename SplatPrimitive::RenderOutput::Tensor> distortions = std::nullopt;
    if (output_distortion) {
        renders2 = SplatPrimitive::RenderOutput::Tensor::empty(renders_dims, opt);
        distortions = SplatPrimitive::RenderOutput::Tensor::empty(renders_dims, opt);
    }

    at::DimVector transmittance_dims(image_dims);
    transmittance_dims.append({image_height, image_width, 1});
    at::Tensor transmittances = at::empty(transmittance_dims, opt);

    at::DimVector last_ids_dims(image_dims);
    last_ids_dims.append({image_height, image_width});
    at::Tensor last_ids = at::empty(last_ids_dims, opt.dtype(at::kInt));

    std::optional<at::Tensor> out_max_blending, out_accum_weight;
    if (output_max_blending || accum_weight_map.has_value()) {
        out_max_blending = at::empty({splats.size()}, opt);
        set_zero<float>(out_max_blending.value());
    }
    if (accum_weight_map.has_value()) {
        out_accum_weight = at::empty({splats.size()}, opt);
        set_zero<float>(out_accum_weight.value());
    }

    if (accum_weight_map.has_value())
        launch_rasterize_to_pixels_eval3d_fwd_kernel<SplatPrimitive, output_distortion, true, output_max_blending>(
            splats, gaussian_ids,
            viewmats, intrins, camera_model, dist_coeffs,
            backgrounds, accum_weight_map, max_blending_masks,
            image_width, image_height, tile_offsets, flatten_ids,
            renders, transmittances, last_ids,
            output_distortion ? &renders2.value() : nullptr,
            output_distortion ? &distortions.value() : nullptr,
            out_accum_weight, out_max_blending
        );
    else
        launch_rasterize_to_pixels_eval3d_fwd_kernel<SplatPrimitive, output_distortion, false, output_max_blending>(
            splats, gaussian_ids,
            viewmats, intrins, camera_model, dist_coeffs,
            backgrounds, accum_weight_map, max_blending_masks,
            image_width, image_height, tile_offsets, flatten_ids,
            renders, transmittances, last_ids,
            output_distortion ? &renders2.value() : nullptr,
            output_distortion ? &distortions.value() : nullptr,
            out_accum_weight, out_max_blending
        );

    if (output_distortion)
        return std::make_tuple(renders.tuple(), transmittances, last_ids,
            renders2.value().tuple(), distortions.value().tuple(), out_accum_weight, out_max_blending);
    return std::make_tuple(renders.tuple(), transmittances, last_ids,
        std::nullopt, std::nullopt, out_accum_weight, out_max_blending);
}



// ================
// Vanilla3DGUT
// ================

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    Vanilla3DGUT::RenderOutput::TensorTuple,
    at::Tensor,
    at::Tensor,
    std::optional<Vanilla3DGUT::RenderOutput::TensorTuple>,
    std::optional<Vanilla3DGUT::RenderOutput::TensorTuple>,
    std::optional<at::Tensor>,
    std::optional<at::Tensor>
> rasterize_to_pixels_3dgut_fwd(
    // Gaussian parameters
    Vanilla3DGUT::Screen::TensorTuple splats_tuple,
    std::optional<at::Tensor> gaussian_ids,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    const std::optional<at::Tensor> backgrounds, // [..., channels]
    const std::optional<at::Tensor> accum_weight_map,       // [..., tile_height, tile_width]
    const std::optional<at::Tensor> max_blending_masks,       // [..., tile_height, tile_width]
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    // intersections
    const at::Tensor tile_offsets, // [..., tile_height, tile_width]
    const at::Tensor flatten_ids,   // [n_isects]
    bool output_distortion
) {
    if (output_distortion)
        return rasterize_to_pixels_eval3d_fwd_tensor<Vanilla3DGUT, true, false>(
            splats_tuple, gaussian_ids,
            viewmats, intrins, cmt(camera_model), dist_coeffs,
            backgrounds, accum_weight_map, max_blending_masks,
            image_width, image_height,
            tile_offsets, flatten_ids
        );
    return rasterize_to_pixels_eval3d_fwd_tensor<Vanilla3DGUT, false, false>(
        splats_tuple, gaussian_ids,
        viewmats, intrins, cmt(camera_model), dist_coeffs,
        backgrounds, accum_weight_map, max_blending_masks,
        image_width, image_height,
        tile_offsets, flatten_ids
    );
}



// ================
// SphericalVoronoi3DGUT
// ================

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    SphericalVoronoi3DGUT_Default::RenderOutput::TensorTuple,
    at::Tensor,
    at::Tensor,
    std::optional<SphericalVoronoi3DGUT_Default::RenderOutput::TensorTuple>,
    std::optional<SphericalVoronoi3DGUT_Default::RenderOutput::TensorTuple>,
    std::optional<at::Tensor>,
    std::optional<at::Tensor>
> rasterize_to_pixels_3dgut_sv_fwd(
    // Gaussian parameters
    SphericalVoronoi3DGUT_Default::Screen::TensorTuple splats_tuple,
    std::optional<at::Tensor> gaussian_ids,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    const std::optional<at::Tensor> backgrounds, // [..., channels]
    const std::optional<at::Tensor> accum_weight_map,       // [..., tile_height, tile_width]
    const std::optional<at::Tensor> max_blending_masks,       // [..., tile_height, tile_width]
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    // intersections
    const at::Tensor tile_offsets, // [..., tile_height, tile_width]
    const at::Tensor flatten_ids,   // [n_isects]
    bool output_distortion
) {
    if (output_distortion)
        // return rasterize_to_pixels_eval3d_fwd_tensor<SphericalVoronoi3DGUT_Default, true, false>(
        return rasterize_to_pixels_eval3d_fwd_tensor<Vanilla3DGUT, true, false>(
            splats_tuple, gaussian_ids,
            viewmats, intrins, cmt(camera_model), dist_coeffs,
            backgrounds, accum_weight_map, max_blending_masks,
            image_width, image_height,
            tile_offsets, flatten_ids
        );
    // return rasterize_to_pixels_eval3d_fwd_tensor<SphericalVoronoi3DGUT_Default, false, false>(
    return rasterize_to_pixels_eval3d_fwd_tensor<Vanilla3DGUT, false, false>(
        splats_tuple, gaussian_ids,
        viewmats, intrins, cmt(camera_model), dist_coeffs,
        backgrounds, accum_weight_map, max_blending_masks,
        image_width, image_height,
        tile_offsets, flatten_ids
    );
}


// ================
// VoxelPrimitive
// ================

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    VoxelPrimitive::RenderOutput::TensorTuple,
    at::Tensor,
    at::Tensor,
    std::optional<VoxelPrimitive::RenderOutput::TensorTuple>,
    std::optional<VoxelPrimitive::RenderOutput::TensorTuple>,
    std::optional<at::Tensor>,
    std::optional<at::Tensor>
> rasterize_to_pixels_voxel_eval3d_fwd(
    // Gaussian parameters
    VoxelPrimitive::Screen::TensorTuple splats_tuple,
    std::optional<at::Tensor> gaussian_ids,
    const at::Tensor viewmats,      // [..., C, 4, 4]
    const at::Tensor intrins,       // [..., C, 4], fx, fy, cx, cy
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    const std::optional<at::Tensor> backgrounds, // [..., channels]
    const std::optional<at::Tensor> accum_weight_map,       // [..., tile_height, tile_width]
    const std::optional<at::Tensor> max_blending_masks,       // [..., tile_height, tile_width]
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    // intersections
    const at::Tensor tile_offsets, // [..., tile_height, tile_width]
    const at::Tensor flatten_ids   // [n_isects]
) {
    return rasterize_to_pixels_eval3d_fwd_tensor<VoxelPrimitive, true, true>(
        splats_tuple, gaussian_ids,
        viewmats, intrins, cmt(camera_model), dist_coeffs,
        backgrounds, accum_weight_map, max_blending_masks,
        image_width, image_height,
        tile_offsets, flatten_ids
    );
}
