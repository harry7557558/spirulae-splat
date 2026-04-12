#pragma once

// Modified from https://github.com/nerfstudio-project/gsplat/blob/main/gsplat/cuda/csrc/RasterizeToPixels3DGSBwd.cu

#include <cuda_runtime.h>
#include <cstdint>

#include <cooperative_groups.h>
namespace cg = cooperative_groups;

#ifdef __CUDACC__
#include "generated/slang.cuh"
namespace SlangProjectionUtils {
#include "generated/set_namespace.cuh"
#include "generated/projection_utils.cuh"
}
#endif

#include <gsplat/Common.h>
#include <gsplat/Utils.cuh>

#ifndef NO_TORCH
#define NO_TORCH
#endif

#include "types.cuh"
#include "common.cuh"


constexpr uint SPLAT_BATCH_SIZE_NO_DISTORTION = WARP_SIZE;
constexpr uint SPLAT_BATCH_SIZE_WITH_DISTORTION = WARP_SIZE;


template <
    typename SplatPrimitive,
    ssplat::CameraModelType camera_model,
    bool output_distortion,
    bool output_viewmat_grad,
    bool output_hessian_diagonal,
    bool output_accum_weight
>
__global__ void rasterize_to_pixels_eval3d_bwd_kernel(
    const uint32_t I,
    const uint32_t n_isects,
    // fwd inputs
    const uint32_t *__restrict__ gaussian_ids,  // [nnz] optional, for packed mode
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
    const float *__restrict__ accum_weight_map_buffer,           // [..., image_height, image_width, 1]
    // grad outputs
    typename SplatPrimitive::RenderOutput::Buffer v_render_output_buffer,
    const float *__restrict__ v_render_Ts, // [..., image_height, image_width, 1]
    typename SplatPrimitive::RenderOutput::Buffer v_distortions_output_buffer,
    // grad inputs
    typename SplatPrimitive::Screen::Buffer v_splat_buffer,
    typename SplatPrimitive::Screen::Buffer vr_splat_buffer,
    typename SplatPrimitive::Screen::Buffer h_splat_buffer,
    float *__restrict__ o_accum_weight,
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
    v_render_Ts += image_id * image_height * image_width;
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

    constexpr uint BLOCK_SIZE = TILE_SIZE * TILE_SIZE;

    // load pixels
    __shared__ float4 shared_ray_d_pix_bin_final[BLOCK_SIZE];
    __shared__ float2 pix_Ts_with_grad[BLOCK_SIZE];
    __shared__ typename SplatPrimitive::RenderOutput v_pix_colors[BLOCK_SIZE];
    // __shared__ float pix_background[CDIM];  // TODO

    __shared__ typename SplatPrimitive::RenderOutput pix_colors[output_distortion ? BLOCK_SIZE : 1];
    __shared__ typename SplatPrimitive::RenderOutput pix2_colors[output_distortion ? BLOCK_SIZE : 1];
    __shared__ typename SplatPrimitive::RenderOutput v_distortion_out[output_distortion ? BLOCK_SIZE : 1];

    __shared__ float hess_weight_map[output_hessian_diagonal ? BLOCK_SIZE : 1];
    __shared__ float accum_weight_map[output_accum_weight ? BLOCK_SIZE : 1];

    float3 ray_o = SlangProjectionUtils::transform_ray_o(R, t);
    float3 total_v_ray_o = make_float3(0.0f, 0.0f, 0.0f);

    __shared__ float3 shared_v_ray_d[output_viewmat_grad ? BLOCK_SIZE : 1];

    constexpr uint SPLAT_BATCH_SIZE_CONST = output_distortion ?
        SPLAT_BATCH_SIZE_WITH_DISTORTION : SPLAT_BATCH_SIZE_NO_DISTORTION;

    #pragma unroll
    for (uint pix_id0 = 0; pix_id0 < BLOCK_SIZE; pix_id0 += SPLAT_BATCH_SIZE_CONST) {
        static_assert(BLOCK_SIZE % SPLAT_BATCH_SIZE_CONST == 0);
        uint pix_id_local = pix_id0 + thread_id;
        int pix_x = block.group_index().z * TILE_SIZE + pix_id_local % TILE_SIZE;
        int pix_y = block.group_index().y * TILE_SIZE + pix_id_local / TILE_SIZE;
        uint pix_id_global = pix_y * image_width + pix_x;
        uint pix_id_image_global = image_id * image_height * image_width + pix_id_global;
        bool inside = (pix_x < image_width && pix_y < image_height);
        
        int32_t bin_final = (inside ? last_ids[pix_id_global] : 0);
        pix_Ts_with_grad[pix_id_local] = {
            (inside ? render_Ts[pix_id_global] : 0.0f),
            // (inside ? -v_render_alphas[pix_id_global] : 0.0f)
            (inside ? v_render_Ts[pix_id_global] : 0.0f)
        };
        v_pix_colors[pix_id_local] = (inside ?
            SplatPrimitive::RenderOutput::load(v_render_output_buffer, pix_id_image_global)
             : SplatPrimitive::RenderOutput::zero());

        const float px = (float)pix_x + 0.5f;
        const float py = (float)pix_y + 0.5f;
        float3 raydir;
        inside &= SlangProjectionUtils::generate_ray(
            {(px-cx)/fx, (py-cy)/fy},
            camera_model == ssplat::CameraModelType::FISHEYE, dist_coeffs,
            &raydir
        );
        float3 ray_d = SlangProjectionUtils::transform_ray_d(R, raydir);  // mul(raydir, R);
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

        if (output_accum_weight) {
            accum_weight_map[pix_id_local] = (accum_weight_map_buffer != nullptr && inside) ?
                accum_weight_map_buffer[pix_id_image_global] : 0.0f;
        }
        if (output_hessian_diagonal) {
            // https://www.desmos.com/calculator/ld9wg7cuxz
            hess_weight_map[pix_id_local] = (loss_map_buffer != nullptr && inside) ?
                0.5f / fmaxf(loss_map_buffer[pix_id_image_global], 1e-6f) : 0.0f;
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
    uint SPLAT_BATCH_SIZE = SPLAT_BATCH_SIZE_CONST;
    if (SPLAT_BATCH_SIZE_CONST > WARP_SIZE) {
        SPLAT_BATCH_SIZE = (uint)sqrtf((float)(range_end - range_start) * (float)BLOCK_SIZE);
        // SPLAT_BATCH_SIZE = min(SPLAT_BATCH_SIZE_CONST, (SPLAT_BATCH_SIZE + WARP_SIZE) & ~(WARP_SIZE-1));
        SPLAT_BATCH_SIZE = min(SPLAT_BATCH_SIZE_CONST, max(SPLAT_BATCH_SIZE, 1u));
    }
    const uint32_t num_splat_batches =
        _CEIL_DIV(range_end - range_start, SPLAT_BATCH_SIZE);

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
            splat = SplatPrimitive::Screen::loadWithPrecompute(splat_buffer, splat_gid, gaussian_ids);
        }

        // accumulate gradient
        typename SplatPrimitive::Screen v_splat = SplatPrimitive::Screen::zero();
        typename SplatPrimitive::Screen vr_splat = SplatPrimitive::Screen::zero();
        typename SplatPrimitive::Screen h_splat = SplatPrimitive::Screen::zero();
        float accum_weight = 0.0f;

        // thread 0 takes last splat, 1 takes second last, etc.
        // at t=0, thread 0 (splat -1) undo pixel 0
        // at t=1, thread 0 (splat -1) undo pixel 1, thread 1 (splat -2) undo pixel 0
        // ......

        // process gaussians in the current batch for this pixel
        // 0 index is the furthest back gaussian in the batch
        for (int t = 0; t < splat_batch_size + BLOCK_SIZE - 1; ++t,
                (SPLAT_BATCH_SIZE_CONST <= WARP_SIZE ? __syncwarp() : __syncthreads())
        ) {
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
                typename SplatPrimitive::Screen v_splat_temp = SplatPrimitive::Screen::zero();
                v_splat_temp.addGradient(splat.evaluate_alpha_vjp(ray_o, ray_d, v_alpha, v_ray_o_alpha, v_ray_d_alpha));
                v_splat_temp.addGradient(splat.evaluate_color_vjp(ray_o, ray_d, v_color, v_ray_o_color, v_ray_d_color));
                v_splat.addGradient(v_splat_temp);
                // https://www.desmos.com/calculator/ld9wg7cuxz
                float weight = 1.0f / sqrtf(fmaxf(v_c.dot(v_c) / 3.0f, 1e-30f));
                vr_splat.addGradient(v_splat_temp, weight);
                weight = hess_weight_map[pix_id] * weight * weight;
                splat.precomputeBackward(v_splat_temp);
                h_splat.addGaussNewtonHessianDiagonal(v_splat_temp, weight);
            } else {
                v_splat.addGradient(splat.evaluate_alpha_vjp(ray_o, ray_d, v_alpha, v_ray_o_alpha, v_ray_d_alpha));
                v_splat.addGradient(splat.evaluate_color_vjp(ray_o, ray_d, v_color, v_ray_o_color, v_ray_d_color));
            }

            if (output_accum_weight) {
                accum_weight += accum_weight_map[pix_id] * alpha * T0;
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
            splat.precomputeBackward(v_splat);
            v_splat.atomicAddToBuffer(v_splat_buffer, splat_gid, gaussian_ids);
            if (output_hessian_diagonal) {
                splat.precomputeBackward(vr_splat);
                vr_splat.atomicAddToBuffer(vr_splat_buffer, splat_gid, gaussian_ids);
                h_splat.atomicAddToBuffer(h_splat_buffer, splat_gid, gaussian_ids);
            }
            if (output_accum_weight) {
                uint32_t idx = gaussian_ids ? gaussian_ids[splat_gid] : splat_gid % v_splat_buffer.size;
                atomicAddFVec(o_accum_weight + idx, accum_weight);
            }
        }
    }
    if (output_viewmat_grad) {
        // accumulate to viewmat gradient
        float3x3 v_R;
        float3 v_t;
        // gradient from ray_o (will fill v_R and v_t)
        SlangProjectionUtils::transform_ray_o_vjp(R, t, total_v_ray_o, &v_R, &v_t);
        // gradient from ray_d
        #pragma unroll
        for (uint pix_id0 = 0; pix_id0 < BLOCK_SIZE; pix_id0 += SPLAT_BATCH_SIZE_CONST) {
            uint pix_id_local = pix_id0 + thread_id;
            float4 ray_d_pix_bin_final = shared_ray_d_pix_bin_final[pix_id_local];
            float3 raydir = SlangProjectionUtils::undo_transform_ray_d(R,
                make_float3(
                    ray_d_pix_bin_final.x,
                    ray_d_pix_bin_final.y,
                    ray_d_pix_bin_final.z
                )
            );
            float3 v_ray_d = shared_v_ray_d[pix_id_local];
            float3x3 v_R_delta;
            float3 temp;
            SlangProjectionUtils::transform_ray_d_vjp(R, raydir, v_ray_d, &v_R_delta, &temp);
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


template <
    typename SplatPrimitive,
    ssplat::CameraModelType camera_model,
    bool output_distortion,
    bool output_viewmat_grad,
    bool output_hessian_diagonal,
    bool output_accum_weight
>
void rasterize_to_pixels_eval3d_bwd_kernel_wrapper(
    cudaStream_t stream,
    const uint32_t I,
    const uint32_t n_isects,
    // fwd inputs
    const uint32_t *__restrict__ gaussian_ids,  // [nnz] optional, for packed mode
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
    const float *__restrict__ accum_weight_map_buffer,           // [..., image_height, image_width, 1]
    // grad outputs
    typename SplatPrimitive::RenderOutput::Buffer v_render_output_buffer,
    const float *__restrict__ v_render_Ts, // [..., image_height, image_width, 1]
    typename SplatPrimitive::RenderOutput::Buffer v_distortions_output_buffer,
    // grad inputs
    typename SplatPrimitive::Screen::Buffer v_splat_buffer,
    typename SplatPrimitive::Screen::Buffer vr_splat_buffer,
    typename SplatPrimitive::Screen::Buffer h_splat_buffer,
    float *__restrict__ o_accum_weight,
    float *__restrict__ v_viewmats // [B, C, 4, 4]
) {
    dim3 threads = {output_distortion ?
        SPLAT_BATCH_SIZE_WITH_DISTORTION : SPLAT_BATCH_SIZE_NO_DISTORTION,
        1, 1};
    dim3 grid = {I, tile_height, tile_width};

    rasterize_to_pixels_eval3d_bwd_kernel<
        SplatPrimitive, camera_model, output_distortion, output_viewmat_grad, output_hessian_diagonal, output_accum_weight
    ><<<grid, threads, 0, stream>>>(
        I, n_isects,
        gaussian_ids, splat_buffer,
        viewmats, intrins, dist_coeffs_buffer, backgrounds, masks,
        image_width, image_height, tile_width, tile_height,
        tile_offsets, flatten_ids,
        render_Ts, last_ids,
        render_output_buffer, render2_output_buffer, loss_map_buffer, accum_weight_map_buffer,
        v_render_output_buffer, v_render_Ts,
        v_distortions_output_buffer, v_splat_buffer, vr_splat_buffer, h_splat_buffer,
        o_accum_weight, v_viewmats
    );

}
