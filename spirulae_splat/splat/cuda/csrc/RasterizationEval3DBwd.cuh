#pragma once

// Modified from https://github.com/nerfstudio-project/gsplat/blob/main/gsplat/cuda/csrc/RasterizeToPixels3DGSBwd.cu

#include <cuda.h>
#include <cuda_runtime.h>
#include <cstdint>

#include <ATen/Tensor.h>

#include "Primitive3DGUT.cuh"
#include "PrimitiveVoxel.cuh"

#include "types.cuh"

#include "common.cuh"

#include <gsplat/Common.h>
#include <gsplat/Utils.cuh>

#include <cub/cub.cuh>



template <
    typename SplatPrimitive,
    gsplat::CameraModelType camera_model,
    bool output_distortion,
    bool output_viewmat_grad,
    bool output_hessian_diagonal
>
__global__ void rasterize_to_pixels_eval3d_bwd_kernel(
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
    const float *__restrict__ loss_map_buffer,           // [..., image_height, image_width, 1]
    // grad outputs
    typename SplatPrimitive::RenderOutput::Buffer v_render_output_buffer,
    const float *__restrict__ v_render_alphas, // [..., image_height, image_width, 1]
    typename SplatPrimitive::RenderOutput::Buffer v_distortions_output_buffer,
    // grad inputs
    typename SplatPrimitive::Screen::Buffer v_splat_buffer,
    typename SplatPrimitive::Screen::Buffer vr_splat_buffer,
    typename SplatPrimitive::Screen::Buffer h_splat_buffer,
    float *__restrict__ v_viewmats // [B, C, 4, 4]
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
    viewmats += image_id * 16;  // world to camera
    float4 intrin = intrins[image_id];
    float3x3 R = {  // row major
        viewmats[0], viewmats[1], viewmats[2],  // 1st row
        viewmats[4], viewmats[5], viewmats[6],  // 2nd row
        viewmats[8], viewmats[9], viewmats[10],  // 3rd row
    };
    float3 t = { viewmats[3], viewmats[7], viewmats[11] };
    float fx = intrin.x, fy = intrin.y, cx = intrin.z, cy = intrin.w;
    CameraDistortionCoeffs dist_coeffs = dist_coeffs_buffer.load(image_id);

    // load pixels
    __shared__ float4 shared_ray_d_pix_bin_final[BLOCK_SIZE];
    __shared__ float2 pix_Ts_with_grad[BLOCK_SIZE];
    __shared__ typename SplatPrimitive::RenderOutput v_pix_colors[BLOCK_SIZE];
    // __shared__ float pix_background[CDIM];  // TODO

    __shared__ typename SplatPrimitive::RenderOutput pix_colors[output_distortion ? BLOCK_SIZE : 1];
    __shared__ typename SplatPrimitive::RenderOutput pix2_colors[output_distortion ? BLOCK_SIZE : 1];
    __shared__ typename SplatPrimitive::RenderOutput v_distortion_out[output_distortion ? BLOCK_SIZE : 1];

    __shared__ float residual_map[output_hessian_diagonal ? BLOCK_SIZE : 1];

    float3 ray_o = transform_ray_o(R, t);
    float3 total_v_ray_o = make_float3(0.0f, 0.0f, 0.0f);

    __shared__ float3 shared_v_ray_d[output_viewmat_grad ? BLOCK_SIZE : 1];

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
        float3 raydir;
        inside &= generate_ray(
            {(px-cx)/fx, (py-cy)/fy},
            camera_model == gsplat::CameraModelType::FISHEYE, dist_coeffs,
            &raydir
        );
        float3 ray_d = transform_ray_d(R, raydir);  // mul(raydir, R);
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

        if (output_viewmat_grad) {
            shared_v_ray_d[pix_id_local] = make_float3(0.0f, 0.0f, 0.0f);
        }

        if (output_hessian_diagonal) {
            residual_map[pix_id_local] = (loss_map_buffer != nullptr && inside) ?
                sqrtf(fmaxf(loss_map_buffer[pix_id_image_global], 0.0f)) : 1.0f;
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

    const float total_num_pixels = (float)(I * image_width * image_height);

    // if (warp.thread_rank() == 0)
    //     printf("range_start=%d range_end=%d num_splat_batches=%u\n", range_start, range_end, num_splat_batches);
    for (uint32_t splat_b = 0; splat_b < num_splat_batches; ++splat_b) {
        const int32_t splat_batch_end = range_end - 1 - SPLAT_BATCH_SIZE * splat_b;
        const int32_t splat_batch_size = min(SPLAT_BATCH_SIZE, splat_batch_end + 1 - range_start);
        const int32_t splat_idx = splat_batch_end - thread_id;

        // load splats
        typename SplatPrimitive::Screen splat;
        uint32_t splat_gid;
        if (splat_idx >= range_start) {
            splat_gid = flatten_ids[splat_idx]; // flatten index in [I * N] or [nnz]
            splat = SplatPrimitive::Screen::loadWithPrecompute(splat_buffer, splat_gid);
        }

        // accumulate gradient
        typename SplatPrimitive::Screen v_splat = SplatPrimitive::Screen::zero();
        typename SplatPrimitive::Screen vr_splat = SplatPrimitive::Screen::zero();

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
            if (alpha <= ALPHA_THRESHOLD || dot(ray_d, ray_d) == 0.0f)
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
            float3 v_ray_o_alpha, v_ray_d_alpha;
            float3 v_ray_o_color, v_ray_d_color;
            if (output_hessian_diagonal) {
                float weight = residual_map[pix_id] / sqrtf(fmaxf(v_c.dot(v_c) / 3.0f, 1e-30f));
                auto v_splat_i = splat.evaluate_alpha_vjp(ray_o, ray_d, v_alpha, v_ray_o_alpha, v_ray_d_alpha);
                v_splat.addGradient(v_splat_i);
                vr_splat.addGradient(v_splat_i, weight);
                v_splat_i = splat.evaluate_color_vjp(ray_o, ray_d, v_color, v_ray_o_color, v_ray_d_color);
                v_splat.addGradient(v_splat_i);
                vr_splat.addGradient(v_splat_i, weight);
            } else {
                v_splat.addGradient(splat.evaluate_alpha_vjp(ray_o, ray_d, v_alpha, v_ray_o_alpha, v_ray_d_alpha));
                v_splat.addGradient(splat.evaluate_color_vjp(ray_o, ray_d, v_color, v_ray_o_color, v_ray_d_color));
            }

            if (output_viewmat_grad) {
                total_v_ray_o += v_ray_o_alpha + v_ray_o_color;
                shared_v_ray_d[pix_id] += v_ray_d_alpha + v_ray_d_color;
            }

            // update pixel states
            pix_Ts_with_grad[pix_id] = { T0, v_T0 };
            // v_pix_colors remains the same

        }

        // accumulate gradient
        if (splat_idx >= range_start) {
            splat.atomicAddGradientToBuffer(v_splat, v_splat_buffer, splat_gid);
            if (output_hessian_diagonal) {
                splat.atomicAddGradientToBuffer(vr_splat, vr_splat_buffer, splat_gid);
                splat.atomicAddGaussNewtonHessianDiagonalToBuffer(v_splat, h_splat_buffer, splat_gid, total_num_pixels);
            }
        }
    }
    if (output_viewmat_grad) {
        // accumulate to viewmat gradient
        float3x3 v_R;
        float3 v_t;
        // gradient from ray_o (will fill v_R and v_t)
        transform_ray_o_vjp(R, t, total_v_ray_o, &v_R, &v_t);
        // gradient from ray_d
        #pragma unroll
        for (uint pix_id0 = 0; pix_id0 < BLOCK_SIZE; pix_id0 += SPLAT_BATCH_SIZE) {
            uint pix_id_local = pix_id0 + thread_id;
            float4 ray_d_pix_bin_final = shared_ray_d_pix_bin_final[pix_id_local];
            float3 raydir = undo_transform_ray_d(R,
                make_float3(
                    ray_d_pix_bin_final.x,
                    ray_d_pix_bin_final.y,
                    ray_d_pix_bin_final.z
                )
            );
            float3 v_ray_d = shared_v_ray_d[pix_id_local];
            float3x3 v_R_delta;
            float3 temp;
            transform_ray_d_vjp(R, raydir, v_ray_d, &v_R_delta, &temp);
            v_R = v_R + v_R_delta;
        }
        // atomic add to global viewmat gradient
        if (v_viewmats != nullptr) {
            float *v_viewmat = v_viewmats + image_id * 16;
            float temp;
            #define _ATOMIC_ADD(ptr, offset, value) do { \
                temp = isfinite(value) ? value : 0.0f; \
                warpSum(temp, warp); \
                if (warp.thread_rank() == 0 && temp != 0.0f) \
                    atomicAdd((ptr) + (offset), (temp)); \
            } while(0)
            _ATOMIC_ADD(v_viewmat, 0, v_R[0].x);
            _ATOMIC_ADD(v_viewmat, 1, v_R[0].y);
            _ATOMIC_ADD(v_viewmat, 2, v_R[0].z);
            _ATOMIC_ADD(v_viewmat, 3, v_t.x);
            _ATOMIC_ADD(v_viewmat, 4, v_R[1].x);
            _ATOMIC_ADD(v_viewmat, 5, v_R[1].y);
            _ATOMIC_ADD(v_viewmat, 6, v_R[1].z);
            _ATOMIC_ADD(v_viewmat, 7, v_t.y);
            _ATOMIC_ADD(v_viewmat, 8, v_R[2].x);
            _ATOMIC_ADD(v_viewmat, 9, v_R[2].y);
            _ATOMIC_ADD(v_viewmat, 10, v_R[2].z);
            _ATOMIC_ADD(v_viewmat, 11, v_t.z);
            #undef _ATOMIC_ADD
        }
    }
}


template <typename SplatPrimitive, bool output_distortion, bool output_hessian_diagonal>
inline void launch_rasterize_to_pixels_eval3d_bwd_kernel(
    // Gaussian parameters
    typename SplatPrimitive::Screen::Tensor splats,
    const at::Tensor viewmats,             // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
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
    const at::Tensor *loss_map,           // [..., image_height, image_width, 1]
    // gradients of outputs
    typename SplatPrimitive::RenderOutput::Tensor v_render_outputs,
    const at::Tensor v_render_alphas, // [..., image_height, image_width, 1]
    typename SplatPrimitive::RenderOutput::Tensor *v_distortion_outputs,
    // outputs
    typename SplatPrimitive::Screen::Tensor v_splats,
    typename SplatPrimitive::Screen::Tensor *vr_splats,
    typename SplatPrimitive::Screen::Tensor *h_splats,
    std::optional<at::Tensor> v_viewmats
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

    typename SplatPrimitive::Screen::Buffer vr_splats_buffer;
    typename SplatPrimitive::Screen::Buffer h_splats_buffer;
    if (output_hessian_diagonal) {
        vr_splats_buffer = vr_splats->buffer();
        h_splats_buffer = h_splats->buffer();
    }

    #define _LAUNCH_ARGS <<<grid, threads, 0, at::cuda::getCurrentCUDAStream()>>>( \
            I, N, n_isects, packed, \
            splats.buffer(), \
            viewmats.data_ptr<float>(), (float4*)intrins.data_ptr<float>(), dist_coeffs, \
            backgrounds.has_value() ? backgrounds.value().data_ptr<float>() : nullptr, \
            masks.has_value() ? masks.value().data_ptr<bool>() : nullptr, \
            image_width, image_height, tile_width, tile_height, \
            tile_offsets.data_ptr<int32_t>(), flatten_ids.data_ptr<int32_t>(), \
            render_Ts.data_ptr<float>(), last_ids.data_ptr<int32_t>(), \
            output_distortion ? render_outputs->buffer() : typename SplatPrimitive::RenderOutput::Buffer(), \
            output_distortion ? render2_outputs->buffer() : typename SplatPrimitive::RenderOutput::Buffer(), \
            output_hessian_diagonal ? loss_map->data_ptr<float>() : nullptr, \
            v_render_outputs.buffer(), v_render_alphas.data_ptr<float>(), \
            output_distortion ? v_distortion_outputs->buffer() : typename SplatPrimitive::RenderOutput::Buffer(), \
            v_splats.buffer(), vr_splats_buffer, h_splats_buffer, \
            v_viewmats.has_value() ? v_viewmats.value().data_ptr<float>() : nullptr \
        )

    if (camera_model == gsplat::CameraModelType::PINHOLE) {
        if (v_viewmats.has_value())
            rasterize_to_pixels_eval3d_bwd_kernel<SplatPrimitive,
                gsplat::CameraModelType::PINHOLE, output_distortion, true, output_hessian_diagonal> _LAUNCH_ARGS;
        else
            rasterize_to_pixels_eval3d_bwd_kernel<SplatPrimitive,
                gsplat::CameraModelType::PINHOLE, output_distortion, false, output_hessian_diagonal> _LAUNCH_ARGS;
    }
    else if (camera_model == gsplat::CameraModelType::FISHEYE) {
        if (v_viewmats.has_value())
            rasterize_to_pixels_eval3d_bwd_kernel<SplatPrimitive,
                gsplat::CameraModelType::FISHEYE, output_distortion, true, output_hessian_diagonal> _LAUNCH_ARGS;
        else
            rasterize_to_pixels_eval3d_bwd_kernel<SplatPrimitive,
                gsplat::CameraModelType::FISHEYE, output_distortion, false, output_hessian_diagonal> _LAUNCH_ARGS;
    }
    else
        throw std::runtime_error("Unsupported camera model");
    CHECK_DEVICE_ERROR(cudaGetLastError());

    #undef _LAUNCH_ARGS
}


template<typename SplatPrimitive, bool output_distortion, bool output_hessian_diagonal>
inline std::tuple<
    typename SplatPrimitive::Screen::TensorTuple,
    std::optional<at::Tensor>,  // v_viewmats
    std::optional<typename SplatPrimitive::Screen::TensorTuple>,  // jacobian residual product
    std::optional<typename SplatPrimitive::Screen::TensorTuple>  // hessian diagonal
> _rasterize_to_pixels_eval3d_bwd_tensor(
    // Gaussian parameters
    typename SplatPrimitive::Screen::TensorTuple splats_tuple,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
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
    std::optional<at::Tensor> loss_map,  // [..., image_height, image_width, 1]
    // gradients of outputs
    typename SplatPrimitive::RenderOutput::TensorTuple v_render_outputs,
    const at::Tensor v_render_alphas, // [..., image_height, image_width, 1]
    std::optional<typename SplatPrimitive::RenderOutput::TensorTuple> v_distortion_outputs_tuple,
    bool need_viewmat_grad
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
    if (masks.has_value())
        CHECK_INPUT(masks.value());

    if (tile_size != TILE_SIZE)
        AT_ERROR("Unsupported tile size");

    typename SplatPrimitive::Screen::Tensor splats(splats_tuple);
    typename SplatPrimitive::Screen::Tensor v_splats = splats.allocRasterBwd();

    std::optional<at::Tensor> v_viewmats = need_viewmat_grad ?
        (std::optional<at::Tensor>)at::zeros_like(viewmats) : (std::optional<at::Tensor>)std::nullopt;

    std::optional<typename SplatPrimitive::RenderOutput::Tensor> render_outputs = std::nullopt;
    std::optional<typename SplatPrimitive::RenderOutput::Tensor> render2_outputs = std::nullopt;
    std::optional<typename SplatPrimitive::RenderOutput::Tensor> v_distortion_outputs = std::nullopt;
    if (output_distortion) {
        render_outputs = render_outputs_tuple;
        render2_outputs = render2_outputs_tuple;
        v_distortion_outputs = v_distortion_outputs_tuple;
    }
    std::optional<typename SplatPrimitive::Screen::Tensor> vr_splats = std::nullopt;
    std::optional<typename SplatPrimitive::Screen::Tensor> h_splats = std::nullopt;
    if (output_hessian_diagonal) {
        vr_splats = splats.allocRasterBwd();
        h_splats = splats.allocRasterBwd();
    }

    launch_rasterize_to_pixels_eval3d_bwd_kernel<SplatPrimitive, output_distortion, output_hessian_diagonal>(
        splats,
        viewmats, intrins, camera_model, dist_coeffs,
        backgrounds, masks,
        image_width, image_height, tile_offsets, flatten_ids,
        render_Ts, last_ids,
        output_distortion ? &render_outputs.value() : nullptr,
        output_distortion ? &render2_outputs.value() : nullptr,
        output_hessian_diagonal ? &loss_map.value() : nullptr,
        v_render_outputs, v_render_alphas,
        output_distortion ? &v_distortion_outputs.value() : nullptr,
        v_splats,
        output_hessian_diagonal ? &vr_splats.value() : nullptr,
        output_hessian_diagonal ? &h_splats.value() : nullptr,
        v_viewmats
    );

    if (output_hessian_diagonal)
        return std::make_tuple(v_splats.tupleRasterBwd(), v_viewmats,
            (std::optional<typename SplatPrimitive::Screen::TensorTuple>)vr_splats.value().tupleRasterBwd(),
            (std::optional<typename SplatPrimitive::Screen::TensorTuple>)h_splats.value().tupleRasterBwd());
    return std::make_tuple(v_splats.tupleRasterBwd(), v_viewmats,
        (std::optional<typename SplatPrimitive::Screen::TensorTuple>)std::nullopt,
        (std::optional<typename SplatPrimitive::Screen::TensorTuple>)std::nullopt);
}


template<typename SplatPrimitive, bool output_distortion>
inline std::tuple<
    typename SplatPrimitive::Screen::TensorTuple,
    std::optional<at::Tensor>  // v_viewmats
> rasterize_to_pixels_eval3d_bwd_tensor(
    // Gaussian parameters
    typename SplatPrimitive::Screen::TensorTuple &splats_tuple,
    const at::Tensor &viewmats,  // [..., C, 4, 4]
    const at::Tensor &intrins,  // [..., C, 4], fx, fy, cx, cy
    const gsplat::CameraModelType &camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    const std::optional<at::Tensor> &backgrounds, // [..., channels]
    const std::optional<at::Tensor> &masks,       // [..., tile_height, tile_width]
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_size,
    // intersections
    const at::Tensor &tile_offsets, // [..., tile_height, tile_width]
    const at::Tensor &flatten_ids,  // [n_isects]
    // forward outputs
    const at::Tensor &render_Ts, // [..., image_height, image_width, 1]
    const at::Tensor &last_ids,      // [..., image_height, image_width]
    std::optional<typename SplatPrimitive::RenderOutput::TensorTuple> &render_outputs,
    std::optional<typename SplatPrimitive::RenderOutput::TensorTuple> &render2_outputs,
    std::optional<at::Tensor> loss_map,  // [..., image_height, image_width, 1]
    // gradients of outputs
    typename SplatPrimitive::RenderOutput::TensorTuple &v_render_outputs,
    const at::Tensor &v_render_alphas, // [..., image_height, image_width, 1]
    std::optional<typename SplatPrimitive::RenderOutput::TensorTuple> &v_distortion_outputs,
    bool need_viewmat_grad
) {
    auto [v_splats, v_viewmat, vr_splats, h_splats] =
        _rasterize_to_pixels_eval3d_bwd_tensor<SplatPrimitive, output_distortion, false>
    (
        splats_tuple,
        viewmats, intrins, camera_model, dist_coeffs,
        backgrounds, masks,
        image_width, image_height, tile_size, tile_offsets, flatten_ids,
        render_Ts, last_ids, render_outputs, render2_outputs, loss_map,
        v_render_outputs, v_render_alphas, v_distortion_outputs,
        need_viewmat_grad
    );
    return std::make_tuple(v_splats, v_viewmat);
}
