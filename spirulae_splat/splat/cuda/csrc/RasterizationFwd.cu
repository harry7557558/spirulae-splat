#include "RasterizationFwd.cuh"

#include <gsplat/Utils.cuh>

#include "types.cuh"
#include "common.cuh"

#include <c10/cuda/CUDAStream.h>

template <typename SplatPrimitive>
void rasterize_to_pixels_fwd_kernel_wrapper(
    cudaStream_t stream,
    const uint32_t I,
    const uint32_t N,
    const uint32_t n_isects,
    const bool packed,
    typename SplatPrimitive::Screen::Buffer splat_buffer,
    const float3 *__restrict__ backgrounds, // [I, 3]
    const bool *__restrict__ masks,           // [I, tile_height, tile_width]
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_width,
    const uint32_t tile_height,
    const int32_t *__restrict__ tile_offsets, // [I, tile_height, tile_width]
    const int32_t *__restrict__ flatten_ids,  // [n_isects]
    typename SplatPrimitive::RenderOutput::Buffer render_colors, // [I, image_height, image_width, 3]
    float *__restrict__ render_Ts, // [I, image_height, image_width, 1]
    int32_t *__restrict__ last_ids        // [I, image_height, image_width]
);

template <typename SplatPrimitive>
inline void launch_rasterize_to_pixels_fwd_kernel(
    // Gaussian parameters
    typename SplatPrimitive::Screen::Tensor splats,
    const std::optional<at::Tensor> backgrounds, // [..., channels]
    const std::optional<at::Tensor> masks,       // [..., tile_height, tile_width]
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    // intersections
    const at::Tensor tile_offsets, // [..., tile_height, tile_width]
    const at::Tensor flatten_ids,  // [n_isects]
    // outputs
    typename SplatPrimitive::RenderOutput::Tensor renders,
    at::Tensor transmittances,  // [..., image_height, image_width]
    at::Tensor last_ids // [..., image_height, image_width]
) {
    bool packed = splats.isPacked();
    uint32_t N = packed ? 0 : splats.size(); // number of gaussians
    uint32_t I = transmittances.numel() / (image_height * image_width); // number of images
    uint32_t tile_height = tile_offsets.size(-2);
    uint32_t tile_width = tile_offsets.size(-1);
    uint32_t n_isects = flatten_ids.size(0);

    rasterize_to_pixels_fwd_kernel_wrapper<SplatPrimitive>(
        (cudaStream_t)at::cuda::getCurrentCUDAStream(),
        I, N, n_isects,
        packed,
        splats.buffer(),
        backgrounds.has_value() ? (float3*)backgrounds.value().data_ptr<float>()
                                : nullptr,
        masks.has_value() ? masks.value().data_ptr<bool>() : nullptr,
        image_width,
        image_height,
        tile_width,
        tile_height,
        tile_offsets.data_ptr<int32_t>(),
        flatten_ids.data_ptr<int32_t>(),
        renders,
        transmittances.data_ptr<float>(),
        last_ids.data_ptr<int32_t>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


template <typename SplatPrimitive>
inline std::tuple<typename SplatPrimitive::RenderOutput::TensorTuple, at::Tensor, at::Tensor>
rasterize_to_pixels_fwd_tensor(
    // Gaussian parameters
    typename SplatPrimitive::Screen::TensorTuple splats_tuple,
    const std::optional<at::Tensor> backgrounds, // [..., channels]
    const std::optional<at::Tensor> masks,       // [..., tile_height, tile_width]
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
    if (backgrounds.has_value())
        CHECK_INPUT(backgrounds.value());
    if (masks.has_value())
        CHECK_INPUT(masks.value());
    
    typename SplatPrimitive::Screen::Tensor splats(splats_tuple);

    auto opt = splats.options();
    at::DimVector image_dims(tile_offsets.sizes().slice(0, tile_offsets.dim() - 2));

    at::DimVector renders_dims(image_dims);
    renders_dims.append({image_height, image_width});
    typename SplatPrimitive::RenderOutput::Tensor renders =
        SplatPrimitive::RenderOutput::Tensor::empty(renders_dims, opt);

    at::DimVector transmittance_dims(image_dims);
    transmittance_dims.append({image_height, image_width, 1});
    at::Tensor transmittances = at::empty(transmittance_dims, opt);

    at::DimVector last_ids_dims(image_dims);
    last_ids_dims.append({image_height, image_width});
    at::Tensor last_ids = at::empty(last_ids_dims, opt.dtype(at::kInt));

    launch_rasterize_to_pixels_fwd_kernel<SplatPrimitive>(
        splats,
        backgrounds,
        masks,
        image_width,
        image_height,
        tile_offsets,
        flatten_ids,
        renders,
        transmittances,
        last_ids
    );

    return std::make_tuple(renders.tuple(), transmittances, last_ids);
}




// ================
// Vanilla3DGS
// ================

/*[AutoHeaderGeneratorExport]*/
std::tuple<Vanilla3DGS::RenderOutput::TensorTuple, at::Tensor, at::Tensor>
rasterize_to_pixels_3dgs_fwd(
    // Gaussian parameters
    Vanilla3DGS::Screen::TensorTuple splats_tuple,
    const std::optional<at::Tensor> backgrounds, // [..., channels]
    const std::optional<at::Tensor> masks,       // [..., tile_height, tile_width]
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    // intersections
    const at::Tensor tile_offsets, // [..., tile_height, tile_width]
    const at::Tensor flatten_ids   // [n_isects]
) {
    return rasterize_to_pixels_fwd_tensor<Vanilla3DGS>(
        splats_tuple, backgrounds, masks,
        image_width, image_height,
        tile_offsets, flatten_ids
    );
}



// ================
// MipSplatting
// ================

/*[AutoHeaderGeneratorExport]*/
std::tuple<MipSplatting::RenderOutput::TensorTuple, at::Tensor, at::Tensor>
rasterize_to_pixels_mip_fwd(
    // Gaussian parameters
    MipSplatting::Screen::TensorTuple splats_tuple,
    const std::optional<at::Tensor> backgrounds, // [..., channels]
    const std::optional<at::Tensor> masks,       // [..., tile_height, tile_width]
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    // intersections
    const at::Tensor tile_offsets, // [..., tile_height, tile_width]
    const at::Tensor flatten_ids   // [n_isects]
) {
    // return rasterize_to_pixels_fwd_tensor<MipSplatting>(
    return rasterize_to_pixels_fwd_tensor<Vanilla3DGS>(
        splats_tuple, backgrounds, masks,
        image_width, image_height,
        tile_offsets, flatten_ids
    );
}
