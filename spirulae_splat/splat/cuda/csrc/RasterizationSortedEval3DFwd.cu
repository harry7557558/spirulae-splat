// Modified from https://github.com/nerfstudio-project/gsplat/blob/main/gsplat/cuda/csrc/RasterizeToPixels3DGSFwd.cu

#include "RasterizationEval3DFwd.cuh"

#include "common.cuh"

#include <gsplat/Common.h>
#include <gsplat/Utils.cuh>


template<uint MAX_SIZE>
struct MinPriorityQueue {
    uint2 arr[MAX_SIZE];   // x is key and y is sorting value
    uint size;
    
    __forceinline__ __device__ void init() { size = 0; }
    __forceinline__ __device__ bool empty() const { return size == 0; }
    __forceinline__ __device__ void clear() { size = 0; }
    __forceinline__ __device__ bool full() const { return size >= MAX_SIZE;  }
    
    inline __device__ void push(uint key, float value) {
        uint idx = size;
        arr[idx] = make_uint2(key, value <= 0.0f ? 0u : __float_as_uint(value));
        size++;
        while (idx > 0) {
            uint parent = (idx - 1) / 2;
            if (arr[idx].y >= arr[parent].y)
                break;
            uint2 temp = arr[idx];
            arr[idx] = arr[parent];
            arr[parent] = temp;
            idx = parent;
        }
    }

    inline __device__ uint2 pop() {
        uint2 result = arr[0];
        size--;
        if (size > 0) {
            arr[0] = arr[size];
            uint idx = 0;
            while (true) {
                uint left = 2 * idx + 1;
                uint right = 2 * idx + 2;
                uint smallest = idx;
                if (left < size && arr[left].y < arr[smallest].y)
                    smallest = left;
                if (right < size && arr[right].y < arr[smallest].y)
                    smallest = right;
                if (smallest == idx)
                    break;
                uint2 temp = arr[idx];
                arr[idx] = arr[smallest];
                arr[smallest] = temp;
                idx = smallest;
            }
        }
        return result;
    }

};


template <typename SplatPrimitive, gsplat::CameraModelType camera_model, bool output_distortion>
__global__ void rasterize_to_pixels_sorted_eval3d_fwd_kernel(
    const uint32_t I,
    const uint32_t N,
    const uint32_t n_isects,
    const bool packed,
    const typename SplatPrimitive::WorldEval3D::Buffer splat_buffer,
    const float *__restrict__ viewmats, // [B, C, 4, 4]
    const float *__restrict__ Ks,       // [B, C, 3, 3]
    const CameraDistortionCoeffsBuffer dist_coeffs,
    const float3 *__restrict__ backgrounds, // [I, 3]
    const bool *__restrict__ masks,           // [I, tile_height, tile_width]
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
    typename SplatPrimitive::RenderOutput::Buffer render_distortions // [I, image_height, image_width, ...]
) {
    // each thread draws one pixel, but also timeshares caching gaussians in a
    // shared tile

    auto block = cg::this_thread_block();
    int32_t image_id = block.group_index().x;
    int32_t tile_id =
        block.group_index().y * tile_width + block.group_index().z;

    tile_offsets += image_id * tile_height * tile_width;
    render_Ts += image_id * image_height * image_width;
    last_ids += image_id * image_height * image_width;
    if (backgrounds != nullptr) {
        backgrounds += image_id;
    }
    if (masks != nullptr) {
        masks += image_id * tile_height * tile_width;
    }

    // arrange 16x16 tile into 8x4 subtiles (one per warp)
    static_assert(TILE_SIZE == 16);
    static_assert(WARP_SIZE == 32);
    uint tile_idx = block.thread_index().y * TILE_SIZE + block.thread_index().x;
    uint warp_idx = tile_idx / WARP_SIZE,
         lane_idx = tile_idx % WARP_SIZE;
    uint32_t i = block.group_index().y * TILE_SIZE + (warp_idx / 2) * 4 + lane_idx / 8;
    uint32_t j = block.group_index().z * TILE_SIZE + (warp_idx % 2) * 8 + lane_idx % 8;

    float px = (float)j + 0.5f;
    float py = (float)i + 0.5f;
    int32_t pix_id = i * image_width + j;

    // Load camera
    viewmats += image_id * 16;
    Ks += image_id * 9;
    float3x3 R = {
        viewmats[0], viewmats[1], viewmats[2],  // 1st row
        viewmats[4], viewmats[5], viewmats[6],  // 2nd row
        viewmats[8], viewmats[9], viewmats[10],  // 3rd row
    };
    float3 t = { viewmats[3], viewmats[7], viewmats[11] };
    float fx = Ks[0], fy = Ks[4], cx = Ks[2], cy = Ks[5];
    float4 radial_coeffs = dist_coeffs.radial_coeffs != nullptr ?
        dist_coeffs.radial_coeffs[image_id] : make_float4(0.0f);
    float2 tangential_coeffs = dist_coeffs.tangential_coeffs != nullptr ?
        dist_coeffs.tangential_coeffs[image_id] : make_float2(0.0f);
    float2 thin_prism_coeffs = dist_coeffs.thin_prism_coeffs != nullptr ?
        dist_coeffs.thin_prism_coeffs[image_id] : make_float2(0.0f);

    float3 ray_o; float3 ray_d;
    generate_ray(
        R, t, {(px-cx)/fx, (py-cy)/fy}, camera_model == gsplat::CameraModelType::FISHEYE,
        radial_coeffs, tangential_coeffs, thin_prism_coeffs,
        &ray_o, &ray_d
    );

    // return if out of bounds
    // keep not rasterizing threads around for reading data
    bool inside = (i < image_height && j < image_width);
    bool done = !inside;

    // when the mask is provided, render the background color and return
    // if this tile is labeled as False
    if (masks != nullptr && inside && !masks[tile_id]) {
        // TODO
        // render_colors[pix_id] = backgrounds == nullptr ?
        //     SplatPrimitive::RenderOutput(make_float3(0.f)) :
        //     SplatPrimitive::RenderOutput(*backgrounds);
        return;
    }

    // have all threads in tile process the same gaussians in batches
    // first collect gaussians between range.x and range.y in batches
    // which gaussians to look through in this tile
    int32_t range_start = tile_offsets[tile_id];
    int32_t range_end =
        (image_id == I - 1) && (tile_id == tile_width * tile_height - 1)
            ? n_isects
            : tile_offsets[tile_id + 1];

    // current visibility left to render
    // transmittance is gonna be used in the backward pass which requires a high
    // numerical precision so we use double for it. However double make bwd 1.5x
    // slower so we stick with float for now.
    float T = 1.0f;
    // index of most recent gaussian to write to this thread's pixel
    uint32_t cur_idx = 0;

    typename SplatPrimitive::RenderOutput pix_out = SplatPrimitive::RenderOutput::zero();
    typename SplatPrimitive::RenderOutput pix2_out = SplatPrimitive::RenderOutput::zero();
    typename SplatPrimitive::RenderOutput distortion_out = SplatPrimitive::RenderOutput::zero();

    static constexpr uint MAX_PQUEUE_SIZE = 32;
    MinPriorityQueue<MAX_PQUEUE_SIZE> pqueue;
    pqueue.init();

    for (uint32_t t = range_start; t < range_end + MAX_PQUEUE_SIZE + 1; ++t) {

        done |= (t >= range_end && pqueue.empty());
        if (__ballot_sync(~0u, !done) == 0)
            break;

        bool hasSplat = false;
        float depth;
        if (!done && t < range_end) {
            uint32_t splat_idx = flatten_ids[t];
            typename SplatPrimitive::WorldEval3D splat =
                SplatPrimitive::WorldEval3D::loadWithPrecompute(splat_buffer, splat_idx);
            float alpha = splat.evaluate_alpha(ray_o, ray_d);
            hasSplat |= (alpha >= ALPHA_THRESHOLD);
            if (hasSplat)
                depth = splat.evaluate_sorting_depth(ray_o, ray_d);
        }

        if (pqueue.full() || (t >= range_end && !pqueue.empty())) {
            uint32_t t = pqueue.pop().x;
            uint32_t splat_idx = flatten_ids[t];
            typename SplatPrimitive::WorldEval3D splat =
                SplatPrimitive::WorldEval3D::loadWithPrecompute(splat_buffer, splat_idx);
            float alpha = splat.evaluate_alpha(ray_o, ray_d);

            const float next_T = T * (1.0f - alpha);
            if (next_T <= 1e-4f && false) { // this pixel is done: exclusive
                done = true;
                pqueue.clear();
            }
            else {
                const float vis = alpha * T;
                const typename SplatPrimitive::RenderOutput color = splat.evaluate_color(ray_o, ray_d);

                if (output_distortion) {
                    distortion_out += (
                        color * color * (1.0f - T)
                        + color * pix_out * -2.0f
                        + pix2_out
                    ) * vis;
                    pix2_out += color * color * vis;
                }
                pix_out += color * vis;
                cur_idx = t;

                T = next_T;
            }
        }

        if (__ballot_sync(~0u, !done && hasSplat) == 0)
            continue;
        if (hasSplat) {
            pqueue.push(t, depth);
        }
    }

    if (inside) {
        render_Ts[pix_id] = T;
        int pix_id_global = image_id * image_height * image_width + pix_id;
        // TODO: blend background
        pix_out.saveParamsToBuffer(render_colors, pix_id_global);
        // index in bin of last gaussian in this pixel
        last_ids[pix_id] = static_cast<int32_t>(cur_idx);
        // distortion
        if (output_distortion) {
            pix2_out.saveParamsToBuffer(render_colors2, pix_id_global);
            distortion_out.saveParamsToBuffer(render_distortions, pix_id_global);
        }
    }
}

template <typename SplatPrimitive, bool output_distortion>
inline void launch_rasterize_to_pixels_sorted_eval3d_fwd_kernel(
    // Gaussian parameters
    typename SplatPrimitive::WorldEval3D::Tensor splats,
    const at::Tensor viewmats,             // [..., C, 4, 4]
    const at::Tensor Ks,                   // [..., C, 3, 3]
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
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
    at::Tensor last_ids, // [..., image_height, image_width]
    typename SplatPrimitive::RenderOutput::Tensor *renders2,
    typename SplatPrimitive::RenderOutput::Tensor *distortions
) {
    bool packed = splats.isPacked();
    uint32_t N = packed ? 0 : splats.size(); // number of gaussians
    uint32_t I = transmittances.numel() / (image_height * image_width); // number of images
    uint32_t tile_height = tile_offsets.size(-2);
    uint32_t tile_width = tile_offsets.size(-1);
    uint32_t n_isects = flatten_ids.size(0);

    // Each block covers a tile on the image. In total there are
    // I * tile_height * tile_width blocks.
    dim3 threads = {TILE_SIZE, TILE_SIZE, 1};
    dim3 grid = {I, tile_height, tile_width};

    #define _LAUNCH_ARGS <<<grid, threads>>>( \
            I, N, n_isects, packed, \
            splats.buffer(), \
            viewmats.data_ptr<float>(), Ks.data_ptr<float>(), dist_coeffs, \
            backgrounds.has_value() ? (float3*)backgrounds.value().data_ptr<float>() : nullptr, \
            masks.has_value() ? masks.value().data_ptr<bool>() : nullptr, \
            image_width, image_height, tile_width, tile_height, \
            tile_offsets.data_ptr<int32_t>(), flatten_ids.data_ptr<int32_t>(), \
            renders, transmittances.data_ptr<float>(), last_ids.data_ptr<int32_t>(), \
            output_distortion ? renders2->buffer() : typename SplatPrimitive::RenderOutput::Buffer(), \
            output_distortion ? distortions->buffer() : typename SplatPrimitive::RenderOutput::Buffer() \
        )

    if (camera_model == gsplat::CameraModelType::PINHOLE)
        rasterize_to_pixels_sorted_eval3d_fwd_kernel<SplatPrimitive,
            gsplat::CameraModelType::PINHOLE, output_distortion> _LAUNCH_ARGS;
    else if (camera_model == gsplat::CameraModelType::FISHEYE)
        rasterize_to_pixels_sorted_eval3d_fwd_kernel<SplatPrimitive,
            gsplat::CameraModelType::FISHEYE, output_distortion> _LAUNCH_ARGS;
    else
        throw std::runtime_error("Unsupported camera model");

    #undef _LAUNCH_ARGS
}


template <typename SplatPrimitive, bool output_distortion>
inline std::tuple<
    typename SplatPrimitive::RenderOutput::TensorTuple,
    at::Tensor,
    at::Tensor,
    std::optional<typename SplatPrimitive::RenderOutput::TensorTuple>,
    std::optional<typename SplatPrimitive::RenderOutput::TensorTuple>
> rasterize_to_pixels_sorted_eval3d_fwd_tensor(
    // Gaussian parameters
    typename SplatPrimitive::WorldEval3D::TensorTuple splats_tuple,
    const at::Tensor viewmats,             // [..., C, 4, 4]
    const at::Tensor Ks,                   // [..., C, 3, 3]
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
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
    CHECK_INPUT(viewmats);
    CHECK_INPUT(Ks);
    if (backgrounds.has_value())
        CHECK_INPUT(backgrounds.value());
    if (masks.has_value())
        CHECK_INPUT(masks.value());
    
    typename SplatPrimitive::WorldEval3D::Tensor splats(splats_tuple);

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

    launch_rasterize_to_pixels_sorted_eval3d_fwd_kernel<SplatPrimitive, output_distortion>(
        splats,
        viewmats, Ks, camera_model, dist_coeffs,
        backgrounds, masks,
        image_width, image_height, tile_offsets, flatten_ids,
        renders, transmittances, last_ids,
        output_distortion ? &renders2.value() : nullptr,
        output_distortion ? &distortions.value() : nullptr
    );

    if (output_distortion)
        return std::make_tuple(renders.tuple(), transmittances, last_ids,
            renders2.value().tuple(), distortions.value().tuple());
    return std::make_tuple(renders.tuple(), transmittances, last_ids,
        std::nullopt, std::nullopt);
}


std::tuple<
    OpaqueTriangle::RenderOutput::TensorTuple,
    at::Tensor,
    at::Tensor,
    std::optional<OpaqueTriangle::RenderOutput::TensorTuple>,
    std::optional<OpaqueTriangle::RenderOutput::TensorTuple>
> rasterize_to_pixels_opaque_triangle_sorted_eval3d_fwd(
    // Gaussian parameters
    OpaqueTriangle::WorldEval3D::TensorTuple splats_tuple,
    const at::Tensor viewmats,             // [..., C, 4, 4]
    const at::Tensor Ks,                   // [..., C, 3, 3]
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    const std::optional<at::Tensor> backgrounds, // [..., channels]
    const std::optional<at::Tensor> masks,       // [..., tile_height, tile_width]
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_size,
    // intersections
    const at::Tensor tile_offsets, // [..., tile_height, tile_width]
    const at::Tensor flatten_ids   // [n_isects]
) {
    if (tile_size != TILE_SIZE)
        AT_ERROR("Tile size must be " + std::to_string(TILE_SIZE));
    return rasterize_to_pixels_sorted_eval3d_fwd_tensor<OpaqueTriangle, true>(
        splats_tuple,
        viewmats, Ks, camera_model, dist_coeffs,
        backgrounds, masks,
        image_width, image_height,
        tile_offsets, flatten_ids
    );
}
