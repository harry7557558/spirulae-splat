#pragma once

// Modified from https://github.com/nerfstudio-project/gsplat/blob/main/gsplat/cuda/csrc/RasterizeToPixels3DGSFwd.cu

#include <cuda.h>
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

#include "types.cuh"
#include "common.cuh"

#include <gsplat/Common.h>
#include <gsplat/Utils.cuh>


template<
    typename SplatPrimitive,
    ssplat::CameraModelType camera_model,
    bool output_distortion,
    bool output_max_blending
>
__global__ void rasterize_to_pixels_eval3d_fwd_kernel(
    const uint32_t I,
    const uint32_t N,
    const uint32_t n_isects,
    const uint32_t *__restrict__ gaussian_ids,  // [nnz] optional, for packed mode
    const typename SplatPrimitive::Screen::Buffer splat_buffer,
    const float *__restrict__ viewmats, // [B, C, 4, 4]
    const float4 *__restrict__ intrins,  // [B, C, 4], fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const float3 *__restrict__ backgrounds, // [I, 3]
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
    float* __restrict__ out_max_blending
) {
    // each thread draws one pixel, but also timeshares caching gaussians in a
    // shared tile

    auto block = cg::this_thread_block();
    int32_t image_id = block.group_index().x;
    int32_t tile_id =
        block.group_index().y * tile_width + block.group_index().z;
    uint32_t i = block.group_index().y * TILE_SIZE + block.thread_index().y;
    uint32_t j = block.group_index().z * TILE_SIZE + block.thread_index().x;

    tile_offsets += image_id * tile_height * tile_width;
    render_Ts += image_id * image_height * image_width;
    last_ids += image_id * image_height * image_width;
    if (backgrounds != nullptr) {
        backgrounds += image_id;
    }
    if (max_blending_masks != nullptr)
        max_blending_masks += image_id * image_height * image_width;

    float px = (float)j + 0.5f;
    float py = (float)i + 0.5f;
    int32_t pix_id = i * image_width + j;

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

    bool inside = (i < image_height && j < image_width);

    float3 raydir;
    inside &= SlangProjectionUtils::generate_ray(
        {(px-cx)/fx, (py-cy)/fy},
        camera_model == ssplat::CameraModelType::FISHEYE, dist_coeffs,
        &raydir
    );
    float3 ray_o = SlangProjectionUtils::transform_ray_o(R, t);
    float3 ray_d = SlangProjectionUtils::transform_ray_d(R, raydir);

    bool done = !inside;

    // have all threads in tile process the same gaussians in batches
    // first collect gaussians between range.x and range.y in batches
    // which gaussians to look through in this tile
    int32_t range_start = tile_offsets[tile_id];
    int32_t range_end =
        (image_id == I - 1) && (tile_id == tile_width * tile_height - 1)
            ? n_isects
            : tile_offsets[tile_id + 1];
    constexpr uint BLOCK_SIZE = TILE_SIZE * TILE_SIZE;
    uint32_t num_batches =
        (range_end - range_start + BLOCK_SIZE - 1) / BLOCK_SIZE;

    __shared__ typename SplatPrimitive::Screen splat_batch[BLOCK_SIZE];
    __shared__ uint32_t splat_idx_batch[output_max_blending ? BLOCK_SIZE : 1];

    // current visibility left to render
    // transmittance is gonna be used in the backward pass which requires a high
    // numerical precision so we use double for it. However double make bwd 1.5x
    // slower so we stick with float for now.
    float T = 1.0f;
    // index of most recent gaussian to write to this thread's pixel
    uint32_t cur_idx = 0;

    // collect and process batches of gaussians
    // each thread loads one gaussian at a time before rasterizing its
    // designated pixel
    uint32_t tr = block.thread_rank();

    bool max_blending_mask = (output_max_blending && max_blending_masks) ?
        (inside ? max_blending_masks[pix_id] : false) : inside;

    typename SplatPrimitive::RenderOutput pix_out = SplatPrimitive::RenderOutput::zero();
    typename SplatPrimitive::RenderOutput pix2_out = SplatPrimitive::RenderOutput::zero();
    typename SplatPrimitive::RenderOutput distortion_out = SplatPrimitive::RenderOutput::zero();

    for (uint32_t b = 0; b < num_batches; ++b) {
        // resync all threads before beginning next batch
        // end early if entire tile is done
        if (__syncthreads_count(done) >= BLOCK_SIZE) {
            break;
        }

        // each thread fetch 1 gaussian from front to back
        // index of gaussian to load
        uint32_t batch_start = range_start + BLOCK_SIZE * b;
        uint32_t idx = batch_start + tr;
        if (idx < range_end) {
            int32_t g = flatten_ids[idx]; // flatten index in [I * N] or [nnz]
            splat_batch[tr] = SplatPrimitive::Screen::loadWithPrecompute(splat_buffer, g, gaussian_ids);
            if (output_max_blending)
                splat_idx_batch[tr] = g;
        }

        // wait for other threads to collect the gaussians in batch
        block.sync();

        // process gaussians in the current batch for this pixel
        uint32_t batch_size = min(BLOCK_SIZE, range_end - batch_start);
        for (uint32_t t = 0; t < batch_size; ++t) {
            float vis = 0.0f;
            if (!done) {
                typename SplatPrimitive::Screen splat = splat_batch[t];
                float alpha = splat.evaluate_alpha(ray_o, ray_d);
                if (alpha > ALPHA_THRESHOLD) {
                    const float next_T = T * (1.0f - alpha);
                    if (next_T > 1e-4f) {
                        vis = alpha * T;
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
                        cur_idx = batch_start + t;
                        T = next_T;
                    }
                    else done = true;
                }
            }

            if (!max_blending_mask)
                vis = 0.0f;
            if (output_max_blending && out_max_blending != nullptr) {
                uint32_t splat_idx = splat_idx_batch[t];
                atomicMax(out_max_blending + (N == 0 ? splat_idx : splat_idx % N), vis);  // TODO: reduce
            }
        }
    }

    if (i < image_height && j < image_width) {
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

template<
    typename SplatPrimitive,
    ssplat::CameraModelType camera_model,
    bool output_distortion,
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
    float* __restrict__ out_max_blending
) {
    // Each block covers a tile on the image. In total there are
    // I * tile_height * tile_width blocks.
    dim3 threads = {TILE_SIZE, TILE_SIZE, 1};
    dim3 grid = {I, tile_height, tile_width};

    rasterize_to_pixels_eval3d_fwd_kernel<
        SplatPrimitive, camera_model, output_distortion, output_max_blending
    ><<<grid, threads, 0, stream>>>(
        I, N, n_isects,
        gaussian_ids, splat_buffer, viewmats, intrins, dist_coeffs_buffer, backgrounds, max_blending_masks,
        image_width, image_height, tile_width, tile_height, tile_offsets, flatten_ids,
        render_colors, render_Ts, last_ids,
        render_colors2, render_distortions, out_max_blending
    );
}
