// Modified from https://github.com/nerfstudio-project/gsplat/blob/main/gsplat/cuda/csrc/RasterizeToPixels3DGSBwd.cu

#include "RasterizationSortedEval3DBwd.cuh"

#include "common.cuh"

#include <gsplat/Common.h>
#include <gsplat/Utils.cuh>

#include <cub/cub.cuh>


constexpr uint BLOCK_SIZE = TILE_SIZE * TILE_SIZE;


template<uint MAX_SIZE>
struct MaxPriorityQueue {
    uint2 arr[MAX_SIZE];   // x is key and y is sorting value
    // uint size;  // commented to prevent spill to local memory

    inline __device__ void push(uint key, uint& size, float value) {
        uint idx = size;
        arr[idx] = make_uint2(key, value <= 0.0f ? 0u : __float_as_uint(value));
        size++;
        while (idx > 0) {
            uint parent = (idx - 1) / 2;
            if (arr[idx].y <= arr[parent].y)
                break;
            uint2 temp = arr[idx];
            arr[idx] = arr[parent];
            arr[parent] = temp;
            idx = parent;
        }
    }

    inline __device__ uint2 pop(uint& size) {
        uint2 result = arr[0];
        size--;
        if (size > 0) {
            uint idx = 0;
            uint2 mval = arr[size];
            arr[0] = mval;
            while (true) {
                uint left = 2 * idx + 1;
                uint right = 2 * idx + 2;
                uint midx = idx;
                if (left < size && arr[left].y > mval.y)
                    midx = left, mval = arr[left];
                if (right < size && arr[right].y > mval.y)
                    midx = right, mval = arr[right];
                if (midx == idx)
                    break;
                uint2 temp = arr[midx];
                mval = arr[midx] = arr[idx];
                arr[idx] = temp;
                idx = midx;
            }
        }
        return result;
    }

};


template <typename SplatPrimitive, gsplat::CameraModelType camera_model, bool output_distortion>
__global__ void rasterize_to_pixels_sorted_eval3d_bwd_kernel(
    const uint32_t I,
    const uint32_t N,
    const uint32_t n_isects,
    const bool packed,
    // fwd inputs
    typename SplatPrimitive::Screen::Buffer splat_buffer,
    const float *__restrict__ viewmats, // [B, C, 4, 4]
    const float4 *__restrict__ intrins,  // [B, C, 4], fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const float *__restrict__ backgrounds, // [..., CDIM] or [nnz, CDIM]
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_width,
    const uint32_t tile_height,
    const int32_t *__restrict__ tile_offsets, // [..., tile_height, tile_width]
    const int32_t *__restrict__ flatten_ids,  // [n_isects]
    // fwd outputs
    const float *__restrict__ render_Ts,      // [..., image_height, image_width, 1]
    const int32_t *__restrict__ last_ids, // [..., image_height, image_width]
    typename SplatPrimitive::RenderOutput::Buffer render_output_buffer,
    typename SplatPrimitive::RenderOutput::Buffer render2_output_buffer,
    // grad outputs
    typename SplatPrimitive::RenderOutput::Buffer v_render_output_buffer,
    const float *__restrict__ v_render_alphas, // [..., image_height, image_width, 1]
    typename SplatPrimitive::RenderOutput::Buffer v_distortions_output_buffer,
    // grad inputs
    typename SplatPrimitive::Screen::Buffer v_splat_buffer
) {
    auto block = cg::this_thread_block();
    cg::thread_block_tile<WARP_SIZE> warp = cg::tiled_partition<WARP_SIZE>(block);
    uint32_t image_id = block.group_index().x;
    uint32_t tile_id = block.group_index().y * tile_width + block.group_index().z;
    uint32_t thread_id = block.thread_rank();

    tile_offsets += image_id * tile_height * tile_width;
    render_Ts += image_id * image_height * image_width;
    last_ids += image_id * image_height * image_width;
    v_render_alphas += image_id * image_height * image_width;
    if (backgrounds != nullptr) {
        backgrounds += image_id;
    }

    // Load camera
    viewmats += image_id * 16;
    float4 intrin = intrins[image_id];
    float3x3 R = {
        viewmats[0], viewmats[1], viewmats[2],  // 1st row
        viewmats[4], viewmats[5], viewmats[6],  // 2nd row
        viewmats[8], viewmats[9], viewmats[10],  // 3rd row
    };
    float3 t = { viewmats[3], viewmats[7], viewmats[11] };
    float fx = intrin.x, fy = intrin.y, cx = intrin.z, cy = intrin.w;
    CameraDistortionCoeffs dist_coeffs = dist_coeffs_buffer.load(image_id);

    // arrange 16x16 tile into 8x4 subtiles (one per warp)
    static_assert(TILE_SIZE == 16);
    static_assert(WARP_SIZE == 32);
    uint tile_idx = block.thread_index().y * TILE_SIZE + block.thread_index().x;
    uint warp_idx = tile_idx / WARP_SIZE,
         lane_idx = tile_idx % WARP_SIZE;
    uint32_t pix_y = block.group_index().y * TILE_SIZE + (warp_idx / 2) * 4 + lane_idx / 8;
    uint32_t pix_x = block.group_index().z * TILE_SIZE + (warp_idx % 2) * 8 + lane_idx % 8;

    int32_t pix_id_global = pix_y * image_width + pix_x;
    uint pix_id_image_global = image_id * image_height * image_width + pix_id_global;
    bool inside = (pix_x < image_width && pix_y < image_height);

    const float px = (float)pix_x + 0.5f;
    const float py = (float)pix_y + 0.5f;
    float3 raydir;
    inside &= generate_ray(
        {(px-cx)/fx, (py-cy)/fy},
        camera_model == gsplat::CameraModelType::FISHEYE, dist_coeffs,
        &raydir
    );
    float3 ray_o = transform_ray_o(R, t);
    float3 ray_d = transform_ray_d(R, raydir);

    float2 pix_Ts_with_grad = {
        (inside ? render_Ts[pix_id_global] : 0.0f),
        (inside ? -v_render_alphas[pix_id_global] : 0.0f)
    };
    typename SplatPrimitive::RenderOutput v_pix_colors = (inside ?
            SplatPrimitive::RenderOutput::load(v_render_output_buffer, pix_id_image_global)
             : SplatPrimitive::RenderOutput::zero());
    // float pix_background[CDIM];  // TODO

    typename SplatPrimitive::RenderOutput pix_colors;
    typename SplatPrimitive::RenderOutput pix2_colors;
    typename SplatPrimitive::RenderOutput v_distortion_out;
    if (output_distortion) {
        pix_colors = (inside ?
            SplatPrimitive::RenderOutput::load(render_output_buffer, pix_id_image_global)
            : SplatPrimitive::RenderOutput::zero());
        pix2_colors = (inside ?
            SplatPrimitive::RenderOutput::load(render2_output_buffer, pix_id_image_global)
            : SplatPrimitive::RenderOutput::zero());
        v_distortion_out = (inside ?
            SplatPrimitive::RenderOutput::load(v_distortions_output_buffer, pix_id_image_global)
            : SplatPrimitive::RenderOutput::zero());
    }

    int32_t range_start = tile_offsets[tile_id];
    int32_t range_end =
        (image_id == I - 1) && (tile_id == tile_width * tile_height - 1)
            ? n_isects
            : tile_offsets[tile_id + 1];
    int32_t bin_final = (inside ? last_ids[pix_id_global] : 0x7fffffff);

    static constexpr uint MAX_PQUEUE_SIZE = 32;
    MaxPriorityQueue<MAX_PQUEUE_SIZE> pqueue;
    uint pqueue_size = 0;

    __shared__ typename SplatPrimitive::Screen splat_batch[BLOCK_SIZE];

    for (int32_t t = range_end-1; t >= (int)range_start - (int)MAX_PQUEUE_SIZE-1; t--) {

        // load splats into shared memory
        if (((range_end-1) - t) % BLOCK_SIZE == 0) {
            block.sync();
            int32_t t1 = t - (int)block.thread_rank();
            if (t1 >= range_start) {
                uint32_t splat_idx = flatten_ids[t1];
                typename SplatPrimitive::Screen splat =
                    SplatPrimitive::Screen::loadWithPrecompute(splat_buffer, splat_idx);
                splat_batch[block.thread_rank()] = splat;
            }
            block.sync();
        }

        // early skip if done
        bool active = inside && (t <= bin_final);
        active &= (t >= range_start || pqueue_size != 0);
        if (__ballot_sync(~0u, active) == 0)
            continue;

        bool hasSplat = false;
        float depth;
        if (active && t >= range_start) {
            // uint32_t splat_idx = flatten_ids[t];
            // typename SplatPrimitive::Screen splat =
            //     SplatPrimitive::Screen::loadWithPrecompute(splat_buffer, splat_idx);
            typename SplatPrimitive::Screen splat =
                splat_batch[((range_end-1) - t) % BLOCK_SIZE];
            float alpha = splat.evaluate_alpha(ray_o, ray_d);
            hasSplat |= (alpha >= ALPHA_THRESHOLD && dot(ray_d, ray_d) > 0.0f);
            if (hasSplat)
                depth = splat.evaluate_sorting_depth(ray_o, ray_d);
        }

        if (pqueue_size >= MAX_PQUEUE_SIZE || (t < range_start && pqueue_size != 0)) {
            uint32_t t = pqueue.pop(pqueue_size).x;
            uint32_t splat_idx = flatten_ids[t];
            typename SplatPrimitive::Screen splat =
                SplatPrimitive::Screen::loadWithPrecompute(splat_buffer, splat_idx);
            float alpha = splat.evaluate_alpha(ray_o, ray_d);

            // forward:
            // \left(c_{1},T_{1}\right)=\left(c_{0}+\alpha_{i}T_{0}c_{i},\ T_{0}\left(1-\alpha_{i}\right)\right)
            float T1 = pix_Ts_with_grad.x;
            float v_T1 = pix_Ts_with_grad.y;

            // undo pixel:
            // T_{0}=\frac{T_{1}}{1-\alpha_{i}}
            float ra = 1.0f / (1.0f - alpha);
            float T0 = T1 * ra;

            typename SplatPrimitive::RenderOutput color = splat.evaluate_color(ray_o, ray_d);
            typename SplatPrimitive::RenderOutput v_c = v_pix_colors;

            // gradient to alpha:
            // \frac{dL}{d\alpha_{i}}
            // = \frac{dL}{dc_{1}}\frac{dc_{1}}{d\alpha_{i}}+\frac{dL}{dT_{1}}\frac{dT_{1}}{d\alpha_{i}}
            // = T_{0}\frac{dL}{dc_{1}}c_{i}-\frac{dL}{dT_{1}}T_{0}
            float v_alpha = T0 * color.dot(v_c) -v_T1 * T0;

            // gradient to color:
            // \frac{dL}{dc_{i}}
            // = \frac{dL}{dc_{1}}\frac{dc_{1}}{dc_{i}}
            // = \alpha_{i}T_{0}\frac{dL}{dc_{1}}
            typename SplatPrimitive::RenderOutput v_color = v_c * (alpha * T0);

            // update pixel gradient:
            // \frac{dL}{dT_{0}}
            // = \frac{dL}{dc_{1}}\frac{dc_{1}}{dT_{0}}+\frac{dL}{dT_{1}}\frac{dT_{1}}{dT_{0}}
            // = \alpha_{i}\frac{dL}{dc_{1}}c_{i}+\frac{dL}{dT_{1}}\left(1-\alpha_{i}\right)
            float v_T0 = alpha * color.dot(v_c) + v_T1 * (1.0f - alpha);

            // distortion
            if (output_distortion) {
                // \left(d_{1},s_{1}\right)=\left(d_{0}+\alpha_{i}T_{0}\left(c_{i}^{2}\left(1-T_{0}\right)-2c_{i}c_{0}+s_{0}\right),\ s_{0}+\alpha_{i}T_{0}c_{i}^{2}\right)
                // \frac{dL}{ds}=0
                typename SplatPrimitive::RenderOutput v_dist = v_distortion_out;
                typename SplatPrimitive::RenderOutput c0 =
                    pix_colors + color * -alpha * T0;
                typename SplatPrimitive::RenderOutput s0 =
                    pix2_colors + color * color * -alpha * T0;
                // \frac{dL}{d\alpha_{i}}=\frac{dL}{dd_{1}}\frac{dd_{1}}{d\alpha_{i}}=T_{0}\left(c_{i}^{2}\left(1-T_{0}\right)-2c_{i}c_{0}+s_{0}\right)\frac{dL}{dd_{1}}
                v_alpha += T0 * (
                    color * color * (1.0f-T0) +
                    color * c0 * -2.0f + s0
                ).dot(v_dist);
                // \frac{dL}{dc_{i}}=\frac{dL}{dd_{1}}\frac{dd_{1}}{dc_{i}}=2\alpha_{i}T_{0}\left(c_{i}\left(1-T_{0}\right)-c_{0}\right)\frac{dL}{dd_{1}}
                v_color += (
                    color * (1.0f-T0) +
                    c0 * -1.0f
                ) * v_dist * (2.0f * alpha * T0);
                // \alpha_{i}\left(c_{i}^{2}\left(1-2T_{0}\right)-2c_{i}c_{0}+s_{0}\right)\frac{dL}{dd_{1}}
                v_T0 += alpha * (
                    color * color * (1.0f-2.0f*T0) +
                    color * c0 * -2.0f + s0
                ).dot(v_dist);
                // undo pixel state
                pix_colors = c0;
                pix2_colors = s0;
            }

            // backward diff splat
            float3 v_ray_o, v_ray_d;  // TODO
            typename SplatPrimitive::Screen v_splat = SplatPrimitive::Screen::zero();
            v_splat.addGradient(splat.evaluate_alpha_vjp(ray_o, ray_d, v_alpha, v_ray_o, v_ray_d));
            v_splat.addGradient(splat.evaluate_color_vjp(ray_o, ray_d, v_color, v_ray_o, v_ray_d));

            // update pixel states
            pix_Ts_with_grad = { T0, v_T0 };
            // v_pix_colors remains the same

            // accumulate gradient
            splat.atomicAddGradientToBuffer(v_splat, v_splat_buffer, splat_idx);
        }

        if (hasSplat) {
            pqueue.push(t, pqueue_size, depth);
        }
    }
    // TODO: gradient to viewmat
}


template <typename SplatPrimitive, bool output_distortion>
inline void launch_rasterize_to_pixels_sorted_eval3d_bwd_kernel(
    // Gaussian parameters
    typename SplatPrimitive::Screen::Tensor splats,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    const std::optional<at::Tensor> backgrounds, // [..., 3]
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    // intersections
    const at::Tensor tile_offsets, // [..., tile_height, tile_width]
    const at::Tensor flatten_ids,  // [n_isects]
    // forward outputs
    const at::Tensor render_Ts, // [..., image_height, image_width, 1]
    const at::Tensor last_ids,      // [..., image_height, image_width]
    typename SplatPrimitive::RenderOutput::Tensor *render_outputs,
    typename SplatPrimitive::RenderOutput::Tensor *render2_outputs,
    // gradients of outputs
    typename SplatPrimitive::RenderOutput::Tensor v_render_outputs,
    const at::Tensor v_render_alphas, // [..., image_height, image_width, 1]
    typename SplatPrimitive::RenderOutput::Tensor *v_distortion_outputs,
    // outputs
    typename SplatPrimitive::Screen::Tensor v_splats
) {
    bool packed = splats.isPacked();
    uint32_t N = packed ? 0 : splats.size(); // number of gaussians
    uint32_t I = render_Ts.numel() / (image_height * image_width); // number of images
    uint32_t tile_height = tile_offsets.size(-2);
    uint32_t tile_width = tile_offsets.size(-1);
    uint32_t n_isects = flatten_ids.size(0);

    // Each block covers a tile on the image. In total there are
    // I * tile_height * tile_width blocks.
    dim3 threads = {TILE_SIZE, TILE_SIZE, 1};
    dim3 grid = {I, tile_height, tile_width};

    if (n_isects == 0) {
        // skip the kernel launch if there are no elements
        return;
    }

    #define _LAUNCH_ARGS <<<grid, threads, 0, at::cuda::getCurrentCUDAStream()>>>( \
            I, N, n_isects, packed, \
            splats.buffer(), \
            viewmats.data_ptr<float>(), (float4*)intrins.data_ptr<float>(), dist_coeffs, \
            backgrounds.has_value() ? backgrounds.value().data_ptr<float>() : nullptr, \
            image_width, image_height, tile_width, tile_height, \
            tile_offsets.data_ptr<int32_t>(), flatten_ids.data_ptr<int32_t>(), \
            render_Ts.data_ptr<float>(), last_ids.data_ptr<int32_t>(), \
            output_distortion ? render_outputs->buffer() : typename SplatPrimitive::RenderOutput::Buffer(), \
            output_distortion ? render2_outputs->buffer() : typename SplatPrimitive::RenderOutput::Buffer(), \
            v_render_outputs.buffer(), v_render_alphas.data_ptr<float>(), \
            output_distortion ? v_distortion_outputs->buffer() : typename SplatPrimitive::RenderOutput::Buffer(), \
            v_splats.buffer() \
        )

    if (camera_model == gsplat::CameraModelType::PINHOLE)
        rasterize_to_pixels_sorted_eval3d_bwd_kernel<SplatPrimitive,
            gsplat::CameraModelType::PINHOLE, output_distortion> _LAUNCH_ARGS;
    else if (camera_model == gsplat::CameraModelType::FISHEYE)
        rasterize_to_pixels_sorted_eval3d_bwd_kernel<SplatPrimitive,
            gsplat::CameraModelType::FISHEYE, output_distortion> _LAUNCH_ARGS;
    else
        throw std::runtime_error("Unsupported camera model");
    CHECK_DEVICE_ERROR(cudaGetLastError());

    #undef _LAUNCH_ARGS
}


template<typename SplatPrimitive, bool output_distortion>
inline std::tuple<
    typename SplatPrimitive::Screen::TensorTuple,
    std::optional<at::Tensor>  // v_viewmats
> rasterize_to_pixels_sorted_eval3d_bwd_tensor(
    // Gaussian parameters
    typename SplatPrimitive::Screen::TensorTuple splats_tuple,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    const std::optional<at::Tensor> backgrounds, // [..., channels]
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_size,
    // intersections
    const at::Tensor tile_offsets, // [..., tile_height, tile_width]
    const at::Tensor flatten_ids,  // [n_isects]
    // forward outputs
    const at::Tensor render_Ts, // [..., image_height, image_width, 1]
    const at::Tensor last_ids,      // [..., image_height, image_width]
    std::optional<typename SplatPrimitive::RenderOutput::TensorTuple> render_outputs_tuple,
    std::optional<typename SplatPrimitive::RenderOutput::TensorTuple> render2_outputs_tuple,
    // gradients of outputs
    typename SplatPrimitive::RenderOutput::TensorTuple v_render_outputs,
    const at::Tensor v_render_alphas, // [..., image_height, image_width, 1]
    std::optional<typename SplatPrimitive::RenderOutput::TensorTuple> v_distortion_outputs_tuple
) {
    DEVICE_GUARD(tile_offsets);
    CHECK_INPUT(tile_offsets);
    CHECK_INPUT(flatten_ids);
    CHECK_INPUT(viewmats);
    CHECK_INPUT(intrins);
    CHECK_INPUT(render_Ts);
    CHECK_INPUT(last_ids);
    CHECK_INPUT(v_render_alphas);
    if (backgrounds.has_value())
        CHECK_INPUT(backgrounds.value());

    if (tile_size != TILE_SIZE)
        AT_ERROR("Unsupported tile size");

    typename SplatPrimitive::Screen::Tensor splats(splats_tuple);
    typename SplatPrimitive::Screen::Tensor v_splats = splats.allocRasterBwd();

    at::Tensor v_viewmats = at::zeros_like(viewmats);  // TODO

    std::optional<typename SplatPrimitive::RenderOutput::Tensor> render_outputs = std::nullopt;
    std::optional<typename SplatPrimitive::RenderOutput::Tensor> render2_outputs = std::nullopt;
    std::optional<typename SplatPrimitive::RenderOutput::Tensor> v_distortion_outputs = std::nullopt;
    if (output_distortion) {
        render_outputs = render_outputs_tuple;
        render2_outputs = render2_outputs_tuple;
        v_distortion_outputs = v_distortion_outputs_tuple;
    }

    launch_rasterize_to_pixels_sorted_eval3d_bwd_kernel<SplatPrimitive, output_distortion>(
        splats,
        viewmats, intrins, camera_model, dist_coeffs,
        backgrounds,
        image_width, image_height, tile_offsets, flatten_ids,
        render_Ts, last_ids,
        output_distortion ? &render_outputs.value() : nullptr,
        output_distortion ? &render2_outputs.value() : nullptr,
        v_render_outputs, v_render_alphas,
        output_distortion ? &v_distortion_outputs.value() : nullptr,
        v_splats
    );

    return std::make_tuple(v_splats.tupleRasterBwd(), v_viewmats);
}

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    OpaqueTriangle::Screen::TensorTuple,
    std::optional<at::Tensor>  // v_viewmats
> rasterize_to_pixels_opaque_triangle_sorted_bwd(
    // Gaussian parameters
    OpaqueTriangle::Screen::TensorTuple splats_tuple,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    const std::optional<at::Tensor> backgrounds, // [..., channels]
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_size,
    // intersections
    const at::Tensor tile_offsets, // [..., tile_height, tile_width]
    const at::Tensor flatten_ids,  // [n_isects]
    // forward outputs
    const at::Tensor render_Ts, // [..., image_height, image_width, 1]
    const at::Tensor last_ids,      // [..., image_height, image_width]
    std::optional<typename OpaqueTriangle::RenderOutput::TensorTuple> render_outputs,
    std::optional<typename OpaqueTriangle::RenderOutput::TensorTuple> render2_outputs,
    // gradients of outputs
    OpaqueTriangle::RenderOutput::TensorTuple v_render_outputs,
    const at::Tensor v_render_alphas, // [..., image_height, image_width, 1]
    std::optional<typename OpaqueTriangle::RenderOutput::TensorTuple> v_distortion_outputs
) {
    return rasterize_to_pixels_sorted_eval3d_bwd_tensor<OpaqueTriangle, true>(
        splats_tuple,
        viewmats, intrins, camera_model, dist_coeffs,
        backgrounds,
        image_width, image_height, tile_size, tile_offsets, flatten_ids,
        render_Ts, last_ids, render_outputs, render2_outputs,
        v_render_outputs, v_render_alphas, v_distortion_outputs
    );
}
