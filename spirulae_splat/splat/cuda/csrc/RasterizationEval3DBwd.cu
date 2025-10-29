// Modified from https://github.com/nerfstudio-project/gsplat/blob/main/gsplat/cuda/csrc/RasterizeToPixels3DGSBwd.cu

#include "RasterizationEval3DBwd.cuh"

#include "common.cuh"

#include <gsplat/Common.h>
#include <gsplat/Utils.cuh>

#include <cub/cub.cuh>


constexpr uint BLOCK_SIZE = TILE_SIZE * TILE_SIZE;
constexpr uint SPLAT_BATCH_SIZE_NO_DISTORTION = WARP_SIZE;
constexpr uint SPLAT_BATCH_SIZE_WITH_DISTORTION = 128;


template <typename SplatPrimitive, gsplat::CameraModelType camera_model, bool output_distortion>
__global__ void rasterize_to_pixels_eval3d_bwd_kernel(
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

    // load pixels
    __shared__ float4 shared_ray_d_pix_bin_final[BLOCK_SIZE];
    __shared__ float2 pix_Ts_with_grad[BLOCK_SIZE];
    __shared__ typename SplatPrimitive::RenderOutput v_pix_colors[BLOCK_SIZE];
    // __shared__ float pix_background[CDIM];  // TODO

    __shared__ typename SplatPrimitive::RenderOutput pix_colors[output_distortion ? BLOCK_SIZE : 1];
    __shared__ typename SplatPrimitive::RenderOutput pix2_colors[output_distortion ? BLOCK_SIZE : 1];
    __shared__ typename SplatPrimitive::RenderOutput v_distortion_out[output_distortion ? BLOCK_SIZE : 1];

    float3 ray_o;

    constexpr uint SPLAT_BATCH_SIZE = output_distortion ?
        SPLAT_BATCH_SIZE_WITH_DISTORTION : SPLAT_BATCH_SIZE_NO_DISTORTION;

    #pragma unroll
    for (uint pix_id0 = 0; pix_id0 < BLOCK_SIZE; pix_id0 += SPLAT_BATCH_SIZE) {
        static_assert(BLOCK_SIZE % SPLAT_BATCH_SIZE == 0);
        uint pix_id_local = pix_id0 + thread_id;
        int pix_x = block.group_index().z * TILE_SIZE + pix_id_local % TILE_SIZE;
        int pix_y = block.group_index().y * TILE_SIZE + pix_id_local / TILE_SIZE;
        uint pix_id_global = pix_y * image_width + pix_x;
        uint pix_id_image_global = image_id * image_height * image_width + pix_id_global;
        bool inside = (pix_x < image_width && pix_y < image_height);
        
        int32_t bin_final = (inside ? last_ids[pix_id_global] : 0);
        pix_Ts_with_grad[pix_id_local] = {
            (inside ? render_Ts[pix_id_global] : 0.0f),
            (inside ? -v_render_alphas[pix_id_global] : 0.0f)
        };
        v_pix_colors[pix_id_local] = (inside ?
            SplatPrimitive::RenderOutput::load(v_render_output_buffer, pix_id_image_global)
             : SplatPrimitive::RenderOutput::zero());

        const float px = (float)pix_x + 0.5f;
        const float py = (float)pix_y + 0.5f;
        float3 ray_d;
        generate_ray(
            R, t, {(px-cx)/fx, (py-cy)/fy}, camera_model == gsplat::CameraModelType::FISHEYE,
            radial_coeffs, tangential_coeffs, thin_prism_coeffs,
            &ray_o, &ray_d
        );
        shared_ray_d_pix_bin_final[pix_id_local] =
            {ray_d.x, ray_d.y, ray_d.z, __int_as_float(bin_final)};

        if (output_distortion) {
            pix_colors[pix_id_local] = (inside ?
                SplatPrimitive::RenderOutput::load(render_output_buffer, pix_id_image_global)
                : SplatPrimitive::RenderOutput::zero());
            pix2_colors[pix_id_local] = (inside ?
                SplatPrimitive::RenderOutput::load(render2_output_buffer, pix_id_image_global)
                : SplatPrimitive::RenderOutput::zero());
            v_distortion_out[pix_id_local] = (inside ?
                SplatPrimitive::RenderOutput::load(v_distortions_output_buffer, pix_id_image_global)
                : SplatPrimitive::RenderOutput::zero());
        }
    }
    block.sync();

    // threads fist load splats, then swept through pixels
    // do this in batches

    int32_t range_start = tile_offsets[tile_id];
    int32_t range_end =
        (image_id == I - 1) && (tile_id == tile_width * tile_height - 1)
            ? n_isects
            : tile_offsets[tile_id + 1];
    const uint32_t num_splat_batches =
        _CEIL_DIV(range_end - range_start, SPLAT_BATCH_SIZE);

    // if (warp.thread_rank() == 0)
    //     printf("range_start=%d range_end=%d num_splat_batches=%u\n", range_start, range_end, num_splat_batches);
    for (uint32_t splat_b = 0; splat_b < num_splat_batches; ++splat_b) {
        const int32_t splat_batch_end = range_end - 1 - SPLAT_BATCH_SIZE * splat_b;
        const int32_t splat_batch_size = min(SPLAT_BATCH_SIZE, splat_batch_end + 1 - range_start);
        const int32_t splat_idx = splat_batch_end - thread_id;

        // load splats
        typename SplatPrimitive::WorldEval3D splat;
        uint32_t splat_gid;
        if (splat_idx >= range_start) {
            splat_gid = flatten_ids[splat_idx]; // flatten index in [I * N] or [nnz]
            splat = SplatPrimitive::WorldEval3D::loadWithPrecompute(splat_buffer, splat_gid);
        }

        // accumulate gradient
        typename SplatPrimitive::WorldEval3D v_splat = SplatPrimitive::WorldEval3D::zero();

        // thread 0 takes last splat, 1 takes second last, etc.
        // at t=0, thread 0 (splat -1) undo pixel 0
        // at t=1, thread 0 (splat -1) undo pixel 1, thread 1 (splat -2) undo pixel 0
        // ......

        // process gaussians in the current batch for this pixel
        // 0 index is the furthest back gaussian in the batch
        for (int t = 0; t < splat_batch_size + BLOCK_SIZE - 1; ++t, __syncwarp()) {
            int pix_id = t - thread_id;
            if (pix_id < 0 || pix_id >= BLOCK_SIZE || splat_idx < range_start)
                continue;
            float4 ray_d_pix_bin_final = shared_ray_d_pix_bin_final[pix_id];
            if (splat_idx > __float_as_int(ray_d_pix_bin_final.w))
                continue;

            // evaluate alpha and early skip
            float3 ray_d = {ray_d_pix_bin_final.x, ray_d_pix_bin_final.y, ray_d_pix_bin_final.z};
            float alpha = splat.evaluate_alpha(ray_o, ray_d);
            if (alpha <= ALPHA_THRESHOLD)
                continue;

            // printf("t=%d, thread %u, splat %d (%u), pix_id %d, pix %d %d\n", t, thread_id, splat_idx-range_start, splat_gid, pix_id, pix_global_x, pix_global_y);

            // forward:
            // \left(c_{1},T_{1}\right)=\left(c_{0}+\alpha_{i}T_{0}c_{i},\ T_{0}\left(1-\alpha_{i}\right)\right)
            float T1 = pix_Ts_with_grad[pix_id].x;
            float v_T1 = pix_Ts_with_grad[pix_id].y;

            // undo pixel:
            // T_{0}=\frac{T_{1}}{1-\alpha_{i}}
            float ra = 1.0f / (1.0f - alpha);
            float T0 = T1 * ra;

            typename SplatPrimitive::RenderOutput color = splat.evaluate_color(ray_o, ray_d);
            typename SplatPrimitive::RenderOutput v_c = v_pix_colors[pix_id];

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
                typename SplatPrimitive::RenderOutput v_dist = v_distortion_out[pix_id];
                typename SplatPrimitive::RenderOutput c0 =
                    pix_colors[pix_id] + color * -alpha * T0;
                typename SplatPrimitive::RenderOutput s0 =
                    pix2_colors[pix_id] + color * color * -alpha * T0;
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
                pix_colors[pix_id] = c0;
                pix2_colors[pix_id] = s0;
            }

            // backward diff splat
            float3 v_ray_o, v_ray_d;  // TODO
            v_splat += splat.evaluate_alpha_vjp(ray_o, ray_d, v_alpha, v_ray_o, v_ray_d);
            v_splat += splat.evaluate_color_vjp(ray_o, ray_d, v_color, v_ray_o, v_ray_d);

            // update pixel states
            pix_Ts_with_grad[pix_id] = { T0, v_T0 };
            // v_pix_colors remains the same

        }

        // accumulate gradient
        if (splat_idx >= range_start)
            splat.atomicAddGradientToBuffer(v_splat, v_splat_buffer, splat_gid);
    }
    // TODO: gradient to viewmat
}


template <typename SplatPrimitive, bool output_distortion>
inline void launch_rasterize_to_pixels_eval3d_bwd_kernel(
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
    dim3 threads = {output_distortion ?
        SPLAT_BATCH_SIZE_WITH_DISTORTION : SPLAT_BATCH_SIZE_NO_DISTORTION,
        1, 1};
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
        rasterize_to_pixels_eval3d_bwd_kernel<SplatPrimitive,
            gsplat::CameraModelType::PINHOLE, output_distortion> _LAUNCH_ARGS;
    else if (camera_model == gsplat::CameraModelType::FISHEYE)
        rasterize_to_pixels_eval3d_bwd_kernel<SplatPrimitive,
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
> rasterize_to_pixels_eval3d_bwd_tensor(
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

    launch_rasterize_to_pixels_eval3d_bwd_kernel<SplatPrimitive, output_distortion>(
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
    Vanilla3DGS::WorldEval3D::TensorTuple,
    std::optional<at::Tensor>  // v_viewmats
> rasterize_to_pixels_3dgs_eval3d_bwd(
    // Gaussian parameters
    Vanilla3DGS::WorldEval3D::TensorTuple splats_tuple,
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
    std::optional<typename Vanilla3DGS::RenderOutput::TensorTuple> render_outputs,
    std::optional<typename Vanilla3DGS::RenderOutput::TensorTuple> render2_outputs,
    // gradients of outputs
    Vanilla3DGS::RenderOutput::TensorTuple v_render_outputs,
    const at::Tensor v_render_alphas, // [..., image_height, image_width, 1]
    std::optional<typename Vanilla3DGS::RenderOutput::TensorTuple> v_distortion_outputs
) {
    return rasterize_to_pixels_eval3d_bwd_tensor<Vanilla3DGS, false>(
        splats_tuple,
        viewmats, Ks, camera_model, dist_coeffs,
        backgrounds, masks,
        image_width, image_height, tile_size, tile_offsets, flatten_ids,
        render_Ts, last_ids, render_outputs, render2_outputs,
        v_render_outputs, v_render_alphas, v_distortion_outputs
    );
}

std::tuple<
    OpaqueTriangle::WorldEval3D::TensorTuple,
    std::optional<at::Tensor>  // v_viewmats
> rasterize_to_pixels_opaque_triangle_eval3d_bwd(
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
    return rasterize_to_pixels_eval3d_bwd_tensor<OpaqueTriangle, true>(
        splats_tuple,
        viewmats, Ks, camera_model, dist_coeffs,
        backgrounds, masks,
        image_width, image_height, tile_size, tile_offsets, flatten_ids,
        render_Ts, last_ids, render_outputs, render2_outputs,
        v_render_outputs, v_render_alphas, v_distortion_outputs
    );
}
