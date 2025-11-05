// Modified from https://github.com/nerfstudio-project/gsplat/blob/main/gsplat/cuda/csrc/RasterizeToPixels3DGSBwd.cu

#include "RasterizationEval3DBwd.cuh"

#include "common.cuh"

#include <gsplat/Common.h>
#include <gsplat/Utils.cuh>

#include <cub/cub.cuh>


<<<<<<< HEAD
template<uint MAX_SIZE>
struct MaxPriorityQueue {
    uint2 arr[MAX_SIZE];   // x is key and y is sorting value
    uint size;
    
    __forceinline__ __device__ void init() { size = 0; }
    __forceinline__ __device__ bool empty() const { return size == 0; }
    __forceinline__ __device__ void clear() { size = 0; }
    __forceinline__ __device__ bool full() const { return size >= MAX_SIZE;  }

    inline __device__ void push(uint key, float value) {
=======
constexpr uint BLOCK_SIZE = TILE_SIZE * TILE_SIZE;


template<uint MAX_SIZE>
struct MaxPriorityQueue {
    uint2 arr[MAX_SIZE];   // x is key and y is sorting value
    // uint size;  // commented to prevent spill to local memory

    inline __device__ void push(uint key, uint& size, float value) {
>>>>>>> 251104-triangle-splatting
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

<<<<<<< HEAD
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
                if (left < size && arr[left].y > arr[smallest].y)
                    smallest = left;
                if (right < size && arr[right].y > arr[smallest].y)
                    smallest = right;
                if (smallest == idx)
                    break;
                uint2 temp = arr[idx];
                arr[idx] = arr[smallest];
                arr[smallest] = temp;
                idx = smallest;
=======
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
>>>>>>> 251104-triangle-splatting
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
    typename SplatPrimitive::WorldEval3D::Buffer splat_buffer,
    const float *__restrict__ viewmats, // [B, C, 4, 4]
    const float *__restrict__ Ks,       // [B, C, 3, 3]
    const CameraDistortionCoeffsBuffer dist_coeffs,
    const float *__restrict__ backgrounds, // [..., CDIM] or [nnz, CDIM]
    const bool *__restrict__ masks,           // [..., tile_height, tile_width]
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
    typename SplatPrimitive::WorldEval3D::Buffer v_splat_buffer
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
    if (masks != nullptr) {
        masks += image_id * tile_height * tile_width;
    }

    // when the mask is provided, do nothing and return if
    // this tile is labeled as False
    if (masks != nullptr && !masks[tile_id]) {
        return;
    }

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

    float2 pix_Ts_with_grad = {
        (inside ? render_Ts[pix_id_global] : 0.0f),
        (inside ? -v_render_alphas[pix_id_global] : 0.0f)
    };
    typename SplatPrimitive::RenderOutput v_pix_colors = (inside ?
            SplatPrimitive::RenderOutput::load(v_render_output_buffer, pix_id_image_global)
             : SplatPrimitive::RenderOutput::zero());
    // float pix_background[CDIM];  // TODO

    const float px = (float)pix_x + 0.5f;
    const float py = (float)pix_y + 0.5f;
    float3 ray_o, ray_d;
    generate_ray(
        R, t, {(px-cx)/fx, (py-cy)/fy}, camera_model == gsplat::CameraModelType::FISHEYE,
        radial_coeffs, tangential_coeffs, thin_prism_coeffs,
        &ray_o, &ray_d
    );

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
<<<<<<< HEAD
    pqueue.init();

    for (int32_t t = range_end-1; t >= (int)range_start - (int)MAX_PQUEUE_SIZE-1; t--) {
        bool active = inside && (t <= bin_final);
        active &= (t >= range_start || !pqueue.empty());
=======
    uint pqueue_size = 0;

    __shared__ typename SplatPrimitive::WorldEval3D splat_batch[BLOCK_SIZE];

    for (int32_t t = range_end-1; t >= (int)range_start - (int)MAX_PQUEUE_SIZE-1; t--) {

        // load splats into shared memory
        if (((range_end-1) - t) % BLOCK_SIZE == 0) {
            block.sync();
            int32_t t1 = t - (int)block.thread_rank();
            if (t1 >= range_start) {
                uint32_t splat_idx = flatten_ids[t1];
                typename SplatPrimitive::WorldEval3D splat =
                    SplatPrimitive::WorldEval3D::loadWithPrecompute(splat_buffer, splat_idx);
                splat_batch[block.thread_rank()] = splat;
            }
            block.sync();
        }

        // early skip if done
        bool active = inside && (t <= bin_final);
        active &= (t >= range_start || pqueue_size != 0);
>>>>>>> 251104-triangle-splatting
        if (__ballot_sync(~0u, active) == 0)
            continue;

        bool hasSplat = false;
        float depth;
        if (active && t >= range_start) {
<<<<<<< HEAD
            uint32_t splat_idx = flatten_ids[t];
            typename SplatPrimitive::WorldEval3D splat =
                SplatPrimitive::WorldEval3D::loadWithPrecompute(splat_buffer, splat_idx);
=======
            // uint32_t splat_idx = flatten_ids[t];
            // typename SplatPrimitive::WorldEval3D splat =
            //     SplatPrimitive::WorldEval3D::loadWithPrecompute(splat_buffer, splat_idx);
            typename SplatPrimitive::WorldEval3D splat =
                splat_batch[((range_end-1) - t) % BLOCK_SIZE];
>>>>>>> 251104-triangle-splatting
            float alpha = splat.evaluate_alpha(ray_o, ray_d);
            hasSplat |= (alpha >= ALPHA_THRESHOLD);
            if (hasSplat)
                depth = splat.evaluate_sorting_depth(ray_o, ray_d);
        }

<<<<<<< HEAD
        if (pqueue.full() || (t < range_start && !pqueue.empty())) {
            uint32_t t = pqueue.pop().x;
=======
        if (pqueue_size >= MAX_PQUEUE_SIZE || (t < range_start && pqueue_size != 0)) {
            uint32_t t = pqueue.pop(pqueue_size).x;
>>>>>>> 251104-triangle-splatting
            uint32_t splat_idx = flatten_ids[t];
            typename SplatPrimitive::WorldEval3D splat =
                SplatPrimitive::WorldEval3D::loadWithPrecompute(splat_buffer, splat_idx);
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
            typename SplatPrimitive::WorldEval3D v_splat = SplatPrimitive::WorldEval3D::zero();
            v_splat += splat.evaluate_alpha_vjp(ray_o, ray_d, v_alpha, v_ray_o, v_ray_d);
            v_splat += splat.evaluate_color_vjp(ray_o, ray_d, v_color, v_ray_o, v_ray_d);

            // update pixel states
            pix_Ts_with_grad = { T0, v_T0 };
            // v_pix_colors remains the same

            // accumulate gradient
            splat.atomicAddGradientToBuffer(v_splat, v_splat_buffer, splat_idx);
        }

<<<<<<< HEAD
        if (__ballot_sync(~0u, active && hasSplat) == 0)
            continue;
        if (hasSplat) {
            pqueue.push(t, depth);
=======
        if (hasSplat) {
            pqueue.push(t, pqueue_size, depth);
>>>>>>> 251104-triangle-splatting
        }
    }
    // TODO: gradient to viewmat
}


template <typename SplatPrimitive, bool output_distortion>
inline void launch_rasterize_to_pixels_sorted_eval3d_bwd_kernel(
    // Gaussian parameters
    typename SplatPrimitive::WorldEval3D::Tensor splats,
    const at::Tensor viewmats,             // [..., C, 4, 4]
    const at::Tensor Ks,                   // [..., C, 3, 3]
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    const std::optional<at::Tensor> backgrounds, // [..., 3]
    const std::optional<at::Tensor> masks,       // [..., tile_height, tile_width]
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
    typename SplatPrimitive::WorldEval3D::Tensor v_splats
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

    #define _LAUNCH_ARGS <<<grid, threads>>>( \
            I, N, n_isects, packed, \
            splats.buffer(), \
            viewmats.data_ptr<float>(), Ks.data_ptr<float>(), dist_coeffs, \
            backgrounds.has_value() ? backgrounds.value().data_ptr<float>() : nullptr, \
            masks.has_value() ? masks.value().data_ptr<bool>() : nullptr, \
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
    typename SplatPrimitive::WorldEval3D::TensorTuple,
    std::optional<at::Tensor>  // v_viewmats
> rasterize_to_pixels_sorted_eval3d_bwd_tensor(
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
    CHECK_INPUT(Ks);
    CHECK_INPUT(render_Ts);
    CHECK_INPUT(last_ids);
    CHECK_INPUT(v_render_alphas);
    if (backgrounds.has_value())
        CHECK_INPUT(backgrounds.value());
    if (masks.has_value())
        CHECK_INPUT(masks.value());

    if (tile_size != TILE_SIZE)
        AT_ERROR("Unsupported tile size");

    typename SplatPrimitive::WorldEval3D::Tensor splats(splats_tuple);
    typename SplatPrimitive::WorldEval3D::Tensor v_splats = splats.zeros_like();

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
        viewmats, Ks, camera_model, dist_coeffs,
        backgrounds, masks,
        image_width, image_height, tile_offsets, flatten_ids,
        render_Ts, last_ids,
        output_distortion ? &render_outputs.value() : nullptr,
        output_distortion ? &render2_outputs.value() : nullptr,
        v_render_outputs, v_render_alphas,
        output_distortion ? &v_distortion_outputs.value() : nullptr,
        v_splats
    );

    return std::make_tuple(v_splats.tupleAll(), v_viewmats);
}

std::tuple<
    OpaqueTriangle::WorldEval3D::TensorTuple,
    std::optional<at::Tensor>  // v_viewmats
> rasterize_to_pixels_opaque_triangle_sorted_eval3d_bwd(
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
        viewmats, Ks, camera_model, dist_coeffs,
        backgrounds, masks,
        image_width, image_height, tile_size, tile_offsets, flatten_ids,
        render_Ts, last_ids, render_outputs, render2_outputs,
        v_render_outputs, v_render_alphas, v_distortion_outputs
    );
}
