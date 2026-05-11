#define IS_EVAL3D 0
#include "RasterizationEval3DFwd_kernel.cuh"

#if 0
#pragma once

// Modified from https://github.com/nerfstudio-project/gsplat/blob/main/gsplat/cuda/csrc/RasterizeToPixels3DGSFwd.cu

#include <cuda_runtime.h>
#include <cstdint>

#include <cooperative_groups.h>
namespace cg = cooperative_groups;

#include <gsplat/Common.h>
#include <gsplat/Utils.cuh>

#ifndef NO_TORCH
#define NO_TORCH
#endif

#include "types.cuh"
#include "common.cuh"

#include "Primitive.cuh"


template <typename SplatPrimitive>
__global__ void rasterize_to_pixels_fwd_kernel(
    const uint32_t I,
    const uint32_t N,
    const uint32_t n_isects,
    const bool packed,
    typename SplatPrimitive::WorldBuffer splat_wbuffer,
    typename SplatPrimitive::ScreenBuffer splat_sbuffer,
    const float3 *__restrict__ backgrounds, // [I, 3]
    const bool *__restrict__ masks,           // [I, tile_height, tile_width]
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_width,
    const uint32_t tile_height,
    const int32_t *__restrict__ tile_offsets, // [I, tile_height, tile_width]
    const int32_t *__restrict__ flatten_ids,  // [n_isects]
    RenderOutput::Buffer render_colors, // [I, image_height, image_width, 3]
    float *__restrict__ render_Ts, // [I, image_height, image_width, 1]
    int32_t *__restrict__ last_ids        // [I, image_height, image_width]
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
    if (masks != nullptr) {
        masks += image_id * tile_height * tile_width;
    }

    float px = (float)j + 0.5f;
    float py = (float)i + 0.5f;
    int32_t pix_id = i * image_width + j;

    // return if out of bounds
    // keep not rasterizing threads around for reading data
    bool inside = (i < image_height && j < image_width);
    bool done = !inside;

    // when the mask is provided, render the background color and return
    // if this tile is labeled as False
    if (masks != nullptr && inside && !masks[tile_id]) {
        // TODO
        // render_colors[pix_id] = backgrounds == nullptr ?
        //     RenderOutput(make_float3(0.f)) :
        //     RenderOutput(*backgrounds);
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
    const uint32_t BLOCK_SIZE = TILE_SIZE * TILE_SIZE;
    uint32_t num_batches =
        (range_end - range_start + BLOCK_SIZE - 1) / BLOCK_SIZE;

    __shared__ typename SplatPrimitive::Screen splat_batch[BLOCK_SIZE];

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

    RenderOutput pix_out = RenderOutput::zero();
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
            splat_batch[tr].load(splat_buffer, g);
        }

        // wait for other threads to collect the gaussians in batch
        block.sync();

        // process gaussians in the current batch for this pixel
        uint32_t batch_size = min(BLOCK_SIZE, range_end - batch_start);
        for (uint32_t t = 0; (t < batch_size) && !done; ++t) {
            typename SplatPrimitive::Screen splat = splat_batch[t];
            float alpha = splat.evaluate_alpha(px, py);
            if (alpha < ALPHA_THRESHOLD) {
                continue;
            }

            const float next_T = T * (1.0f - alpha);
            if (next_T <= 1e-4f) { // this pixel is done: exclusive
                done = true;
                break;
            }

            const float vis = alpha * T;
            const RenderOutput color = splat.evaluate_color(px, py);
            pix_out += color * vis;
            cur_idx = batch_start + t;

            T = next_T;
        }
    }

    if (inside) {
        render_Ts[pix_id] = T;
        // TODO: blend background
        pix_out.saveParamsToBuffer<SplatPrimitive::pixelType>(render_colors, image_id * image_height * image_width + pix_id);
        // index in bin of last gaussian in this pixel
        last_ids[pix_id] = static_cast<int32_t>(cur_idx);
    }
}


template <typename SplatPrimitive>
void rasterize_to_pixels_fwd_kernel_wrapper(
    cudaStream_t stream,
    const uint32_t I,
    const uint32_t N,
    const uint32_t n_isects,
    const bool packed,
    typename SplatPrimitive::ScreenBuffer splat_buffer,
    const float3 *__restrict__ backgrounds, // [I, 3]
    const bool *__restrict__ masks,           // [I, tile_height, tile_width]
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_width,
    const uint32_t tile_height,
    const int32_t *__restrict__ tile_offsets, // [I, tile_height, tile_width]
    const int32_t *__restrict__ flatten_ids,  // [n_isects]
    RenderOutput::Buffer render_colors, // [I, image_height, image_width, 3]
    float *__restrict__ render_Ts, // [I, image_height, image_width, 1]
    int32_t *__restrict__ last_ids        // [I, image_height, image_width]
) {
    // Each block covers a tile on the image. In total there are
    // I * tile_height * tile_width blocks.
    dim3 threads = {TILE_SIZE, TILE_SIZE, 1};
    dim3 grid = {I, tile_height, tile_width};

    rasterize_to_pixels_fwd_kernel<SplatPrimitive><<<grid, threads, 0, stream>>>(
        I, N, n_isects, packed,
        splat_buffer, backgrounds, masks,
        image_width, image_height, tile_width, tile_height, tile_offsets, flatten_ids,
        render_colors, render_Ts, last_ids
    );
}
#endif