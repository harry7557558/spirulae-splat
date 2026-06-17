#ifdef _KERNEL_CUH_INC
#error "Rasterization forward kernel must be included no more than once"
#endif

#define _KERNEL_CUH_INC

#ifndef NO_TORCH
#define NO_TORCH
#endif


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

#include <Common.cuh>

#include "Primitive.cuh"


#ifndef IS_EVAL3D
#define IS_EVAL3D 1
#endif

inline constexpr uint32_t NUM_WARPS = 10;
inline constexpr uint32_t NUM_THREADS = NUM_WARPS * WARP_SIZE;

// Test whether the unit ellipse (x^T inv_cov x = 1, centered at origin) overlaps
// an axis-aligned box [x0,x1]x[y0,y1]. Used to cull (gaussian, sub-tile) pairs.
static __device__ __forceinline__ bool ellipse_box_overlap_test(
    float3 inv_cov, float x0, float x1, float y0, float y1
) {
    float a = inv_cov.x, b = inv_cov.y, c = inv_cov.z;
    // parabola vertex on the 4 edges (uses the un-doubled off-diagonal)
    float wx = -b / c, wy = -b / a;
    float u0 = fminf(fmaxf(x0 * wx, y0), y1);
    float u1 = fminf(fmaxf(x1 * wx, y0), y1);
    float v0 = fminf(fmaxf(y0 * wy, x0), x1);
    float v1 = fminf(fmaxf(y1 * wy, x0), x1);
    b *= 2.0f;
    float mx = fminf(
        a*x0*x0 + b*x0*u0 + c*u0*u0,
        a*x1*x1 + b*x1*u1 + c*u1*u1
    );
    float my = fminf(
        a*v0*v0 + b*v0*y0 + c*y0*y0,
        a*v1*v1 + b*v1*y1 + c*y1*y1
    );
    // box contains the ellipse center, or the box boundary meets the ellipse
    float mc = fmaxf(fmaxf(x0, -x1), fmaxf(y0, -y1));
    return fminf(mc, fminf(mx, my) - 1.0f) < 0.f;
}

template<
    typename SplatPrimitive,
#if IS_EVAL3D
    CameraModelType camera_model,
#endif
    bool output_distortion
>
#if IS_EVAL3D
__global__ void rasterize_to_pixels_eval3d_fwd_kernel(
#else
__global__ void rasterize_to_pixels_fwd_kernel(
#endif
    const uint32_t I,
    const uint32_t N,   // zero if packed
    const uint32_t n_isects,
    const uint32_t *__restrict__ gaussian_ids,  // [nnz], nullptr if not packed
    const typename SplatPrimitive::WorldBuffer splat_wbuffer,
    const typename SplatPrimitive::ScreenBuffer splat_sbuffer,
#if IS_EVAL3D
    const float *__restrict__ viewmats, // [B, C, 4, 4]
    const float4 *__restrict__ intrins,  // [B, C, 4], fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const float4 *__restrict__ aabb,  // [..., N] projected 2D AABB (xmin,ymin,xmax,ymax)
#endif
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_width,
    const uint32_t tile_height,
    const int32_t *__restrict__ tile_offsets, // [I, tile_height, tile_width]
    const int32_t *__restrict__ flatten_ids,  // [n_isects]
    RenderOutput::Buffer render_colors, // [I, image_height, image_width, ...]
    float * render_Ts, // [I, image_height, image_width, 1]
    int32_t *__restrict__ last_ids, // [I, image_height, image_width]
    RenderOutput::Buffer render_colors2, // [I, image_height, image_width, ...]
    RenderOutput::Buffer render_distortions // [I, image_height, image_width, ...]
) {
    uint32_t tid = threadIdx.x;
    uint32_t wid = tid / WARP_SIZE;
    uint32_t lid = tid % WARP_SIZE;
    int32_t image_id = blockIdx.x;
    int32_t tile_id = blockIdx.y * tile_width + blockIdx.z;

    tile_offsets += image_id * tile_height * tile_width;
    render_Ts += image_id * image_height * image_width;
    last_ids += image_id * image_height * image_width;

#if IS_EVAL3D
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

    float3 ray_o = SlangProjectionUtils::transform_ray_o(R, t);
#endif

    // number of warp-batches preloaded per depth batch, sized so the whole
    // batch of gaussians fits in static shared memory (smaller for fatter
    // fragments, e.g. 3DGUT, larger for lighter primitives)
    constexpr uint32_t DEPTH_BATCH_SIZE = bool(IS_EVAL3D) ? 12 : 16;
    constexpr uint32_t DEPTH_BATCH_SPLATS = DEPTH_BATCH_SIZE * WARP_SIZE;

    // gaussians overlapping this macro tile, shared by all of its sub-tiles
    int32_t range_start = tile_offsets[tile_id];
    int32_t range_end =
        (image_id == I - 1) && (tile_id == tile_width * tile_height - 1)
            ? n_isects
            : tile_offsets[tile_id + 1];
    uint32_t num_batches =
        (range_end - range_start + WARP_SIZE - 1) / WARP_SIZE;
    uint32_t num_depth_batches = (num_batches + DEPTH_BATCH_SIZE - 1) / DEPTH_BATCH_SIZE;
    num_depth_batches = max(num_depth_batches, 1);
    // when there is more than one depth batch we round-trip per-pixel state
    // (transmittance / colors / last index) through global memory between
    // batches and have to remember which pixels already saturated.
    const bool needs_resume = num_depth_batches > 1;

    constexpr uint32_t MACRO_NUM_TILES = MACRO_TILE_SIZE_X * MACRO_TILE_SIZE_Y;
    constexpr uint32_t WARPS_PER_TILE = (TILE_SIZE_X * TILE_SIZE_Y) / WARP_SIZE;

    // a sub-tile is fully done once every one of its pixels is outside the image
    // or has saturated; such sub-tiles are dropped from later depth batches and
    // the remaining ones are compacted so warps never pick an idle sub-tile.
    __shared__ bool tile_done[MACRO_NUM_TILES];
    __shared__ int active_tiles[MACRO_NUM_TILES];
    __shared__ int numActiveTiles;
    __shared__ uint32_t numTilesProcessed;
    // a whole depth batch of gaussians is loaded once and shared by every
    // sub-tile of the macro tile (they all sweep the same gaussian range)
    __shared__ typename SplatPrimitive::FragmentFwd splat_batch[DEPTH_BATCH_SPLATS];
    // per-gaussian sub-tile intersection mask: bit s set => gaussian overlaps
    // sub-tile s. MACRO_NUM_TILES == WARP_SIZE, so one uint32 holds all sub-tiles.
    static_assert(MACRO_NUM_TILES <= 32);
    __shared__ uint32_t isect_mask[DEPTH_BATCH_SPLATS];

    if (tid < MACRO_NUM_TILES)
        tile_done[tid] = false;

  for (uint32_t depth_batch = 0; depth_batch < num_depth_batches; ++depth_batch) {

    __syncthreads();

    // compact the still-active sub-tiles into a dense list
    if (tid == 0) {
        int k = 0;
        for (int i = 0; i < (int)MACRO_NUM_TILES; ++i)
            if (!tile_done[i])
                active_tiles[k++] = i;
        numActiveTiles = k;
        numTilesProcessed = 0;
    }

    // preload this depth batch's gaussians into shared memory once
    uint32_t db_base = (uint32_t)range_start + depth_batch * DEPTH_BATCH_SPLATS;
    uint32_t db_count = min((uint32_t)range_end - db_base, DEPTH_BATCH_SPLATS);
    uint32_t num_batches_inner = (db_count + WARP_SIZE - 1) / WARP_SIZE;
    for (uint32_t k = tid; k < db_count; k += NUM_THREADS) {
        int32_t g = flatten_ids[db_base + k]; // flatten index in [I * N] or [nnz]
        splat_batch[k].load(splat_wbuffer, splat_sbuffer,
            gaussian_ids ? gaussian_ids[g] : g % N, g);
    }
#if !IS_EVAL3D
    __syncthreads();  // non-eval3d derives the mask from the just-loaded fragments
#endif

    // precompute, for each gaussian in the batch, the set of sub-tiles its
    // footprint overlaps. One warp per gaussian, lane == sub-tile index, so a
    // single __ballot packs all MACRO_NUM_TILES(<=32) sub-tiles into one word.
    for (uint32_t g = wid; g < db_count; g += NUM_WARPS) {
        bool hit = false;
        if (lid < MACRO_NUM_TILES) {
            uint32_t sx = lid % MACRO_TILE_SIZE_X;
            uint32_t sy = lid / MACRO_TILE_SIZE_X;
            float tx0 = (float)((blockIdx.z * MACRO_TILE_SIZE_X + sx) * TILE_SIZE_X);
            float ty0 = (float)((blockIdx.y * MACRO_TILE_SIZE_Y + sy) * TILE_SIZE_Y);
        #if IS_EVAL3D
            // 3dgut keeps the projected 2D conic in the screen "scale" channel
            // and the ellipse center is the AABB center; same alpha-threshold
            // contour test as 2D splatting, matching the tile intersector.
            int32_t sid = flatten_ids[db_base + g];
            float opac = splat_sbuffer.opacities(sid);
            if (opac > ALPHA_THRESHOLD) {
                float3 conic = splat_sbuffer.scales(sid);
                float4 bb = aabb[sid];
                float ecx = 0.5f * (bb.x + bb.z);
                float ecy = 0.5f * (bb.y + bb.w);
                float kk = 0.5f / __logf(opac / ALPHA_THRESHOLD);
                float3 inv_cov = { conic.x * kk, conic.y * kk, conic.z * kk };
                hit = ellipse_box_overlap_test(
                    inv_cov,
                    tx0 - ecx, tx0 + TILE_SIZE_X - ecx,
                    ty0 - ecy, ty0 + TILE_SIZE_Y - ecy
                );
            }
        #else
            // 2D splatting: exact ellipse vs sub-tile box at the alpha threshold
            typename SplatPrimitive::FragmentFwd sp = splat_batch[g];
            if (sp.opac > ALPHA_THRESHOLD) {
                float kk = 0.5f / __logf(sp.opac / ALPHA_THRESHOLD);
                float3 inv_cov = { sp.conic.x * kk, sp.conic.y * kk, sp.conic.z * kk };
                hit = ellipse_box_overlap_test(
                    inv_cov,
                    tx0 - sp.xy.x, tx0 + TILE_SIZE_X - sp.xy.x,
                    ty0 - sp.xy.y, ty0 + TILE_SIZE_Y - sp.xy.y
                );
            }
        #endif
        }
        uint32_t m = __ballot_sync(~0u, hit);
        if (lid == 0)
            isect_mask[g] = m;
    }
    __syncthreads();

  for (;;) {
    // each warp processes one sub-tile, queued so no warp idles
    uint32_t slot = 0;
    if (lid == 0)
        slot = atomicAdd(&numTilesProcessed, 1);
    slot = __shfl_sync(~0u, slot, 0);
    if (slot >= (uint32_t)numActiveTiles)
        break;
    uint32_t tileIdx = (uint32_t)active_tiles[slot];

    uint32_t tileOffsetY = blockIdx.y * MACRO_TILE_SIZE_Y + (tileIdx / MACRO_TILE_SIZE_X);
    uint32_t tileOffsetX = blockIdx.z * MACRO_TILE_SIZE_X + (tileIdx % MACRO_TILE_SIZE_X);

    // whether every pixel of this sub-tile is done after this depth batch
    bool sub_tile_done = true;

  #pragma unroll
  for (int tr_batch = 0; tr_batch < WARPS_PER_TILE; ++tr_batch) {
    int local_pid = tr_batch * WARP_SIZE + lid;
    uint32_t i = tileOffsetY * TILE_SIZE_Y + local_pid / TILE_SIZE_X;
    uint32_t j = tileOffsetX * TILE_SIZE_X + local_pid % TILE_SIZE_X;

    float px = (float)j + 0.5f;
    float py = (float)i + 0.5f;

    bool inside = (i < image_height && j < image_width);

#if IS_EVAL3D
    float3 raydir;
    inside &= SlangProjectionUtils::generate_ray(
        {(px-cx)/fx, (py-cy)/fy},
        (int)camera_model, dist_coeffs,
        &raydir
    );
    float3 ray_d = SlangProjectionUtils::transform_ray_d(R, raydir);
#endif

    bool done = !inside;
    // a pixel that saturated in an earlier depth batch must not resume
    // accumulating in this one (the gaussians here are further back)
    bool saturated = false;

    // current visibility left to render
    // transmittance is gonna be used in the backward pass which requires a high
    // numerical precision so we use double for it. However double make bwd 1.5x
    // slower so we stick with float for now.
    int32_t pix_id = i * image_width + j;
    int32_t pix_id_global = image_id * image_height * image_width + pix_id;
    float T = 1.0f;
    // index of most recent gaussian to write to this thread's pixel
    uint32_t cur_idx = 0;

    RenderOutput pix_out = RenderOutput::zero();
    RenderOutput pix2_out = RenderOutput::zero();
    RenderOutput distortion_out = RenderOutput::zero();
    if (depth_batch > 0 && inside) {
        // reload state persisted at the end of the previous depth batch
        T = render_Ts[pix_id];
        int32_t s = last_ids[pix_id];
        saturated = (s < 0);
        cur_idx = saturated ? (uint32_t)(-s - 1) : (uint32_t)s;
        done = done || saturated;
        pix_out = render_colors.load<SplatPrimitive::pixelType>(pix_id_global);
        if constexpr (RenderOutput::has_depth(SplatPrimitive::pixelType))
            pix_out.depth *= (1.0f - T);
        if constexpr (output_distortion) {
            pix2_out = render_colors2.load<SplatPrimitive::pixelType>(pix_id_global);
            distortion_out = render_distortions.load<SplatPrimitive::pixelType>(pix_id_global);
        }
    }

    for (uint32_t inner_batch = 0; inner_batch < num_batches_inner; ++inner_batch) {
        // end early if every pixel in this warp is done
        if (__popc(__ballot_sync(~0u, done)) >= WARP_SIZE)
            break;

        // gaussians for this batch are already in shared memory; read directly
        uint32_t local_base = inner_batch * WARP_SIZE;
        uint32_t batch_start = db_base + local_base;  // global flatten index base
        uint32_t batch_size = min((uint32_t)WARP_SIZE, db_count - local_base);
        // skip gaussians whose footprint never reaches this sub-tile; the
        // survivors are compacted via ballot so we only evaluate the hits
        bool lane_hit = (lid < batch_size) &&
            ((isect_mask[local_base + lid] >> tileIdx) & 1u);
        uint32_t surv = __ballot_sync(~0u, lane_hit);
        while (surv) {
            uint32_t t = __ffs(surv) - 1;
            surv &= surv - 1;
            if (done)
                continue;
            typename SplatPrimitive::FragmentFwd splat = splat_batch[local_base + t];
        #if IS_EVAL3D
            float alpha = splat.evaluate_alpha(ray_o, ray_d);
        #else
            float alpha = splat.evaluate_alpha(px, py);
        #endif
            if (alpha <= ALPHA_THRESHOLD)
                continue;

        #if IS_EVAL3D
            const RenderOutput color = splat.evaluate_color(ray_o, ray_d);
            if (color.depth <= 0.0f)
                continue;
        #else
            const RenderOutput color = splat.evaluate_color(px, py);
        #endif

            const float next_T = T * (1.0f - alpha);
            if (next_T > 1e-4f) {
                const float vis = alpha * T;
                if constexpr (output_distortion) {
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
            } else { done = true; saturated = true; }

        }  // while (surv)
    }

    // this sub-tile is only finished once every lane is done
    sub_tile_done &= (__popc(__ballot_sync(~0u, done)) >= WARP_SIZE);

    if (i < image_height && j < image_width) {
        render_Ts[pix_id] = T;
        if constexpr (RenderOutput::has_depth(SplatPrimitive::pixelType))
            pix_out.depth /= fmaxf(1.0f - T, 1e-10f);
        pix_out.saveParamsToBuffer<SplatPrimitive::pixelType>(render_colors, pix_id_global);
        // index in bin of last gaussian in this pixel; saturated pixels are
        // flagged with a negative encoding so a later depth batch knows to stop
        last_ids[pix_id] = (needs_resume && saturated)
            ? -static_cast<int32_t>(cur_idx) - 1
            : static_cast<int32_t>(cur_idx);
        // distortion
        if constexpr (output_distortion) {
            pix2_out.saveParamsToBuffer<SplatPrimitive::pixelType>(render_colors2, pix_id_global);
            distortion_out.saveParamsToBuffer<SplatPrimitive::pixelType>(render_distortions, pix_id_global);
        }
    }

  }  // for (int tr_batch = 0; tr_batch < WARPS_PER_TILE; ++tr_batch)

    if (lid == 0)
        tile_done[tileIdx] = sub_tile_done;

  }  // for (;;)

  }  // for (uint32_t depth_batch = 0; depth_batch < num_depth_batches; ++depth_batch)

    // undo the saturation encoding so last_ids holds plain gaussian indices
    if (needs_resume) {
        __syncthreads();
        constexpr uint32_t MACRO_PIX_H = MACRO_TILE_SIZE_Y * TILE_SIZE_Y;
        constexpr uint32_t MACRO_PIX_W = MACRO_TILE_SIZE_X * TILE_SIZE_X;
        uint32_t base_i = blockIdx.y * MACRO_PIX_H;
        uint32_t base_j = blockIdx.z * MACRO_PIX_W;
        for (uint32_t p = tid; p < MACRO_PIX_H * MACRO_PIX_W; p += NUM_THREADS) {
            uint32_t i = base_i + p / MACRO_PIX_W;
            uint32_t j = base_j + p % MACRO_PIX_W;
            if (i < image_height && j < image_width) {
                int32_t pix_id = i * image_width + j;
                int32_t s = last_ids[pix_id];
                if (s < 0)
                    last_ids[pix_id] = -s - 1;
            }
        }
    }

}

template<
    typename SplatPrimitive,
#if IS_EVAL3D
    CameraModelType camera_model,
#endif
    bool output_distortion
>
#if IS_EVAL3D
void rasterize_to_pixels_eval3d_fwd_kernel_wrapper(
#else
void rasterize_to_pixels_fwd_kernel_wrapper(
#endif
    cudaStream_t stream,
    const uint32_t I,
    const uint32_t N,
    const uint32_t n_isects,
    const uint32_t *__restrict__ gaussian_ids,  // [nnz] optional, for packed mode
    const typename SplatPrimitive::WorldBuffer splat_wbuffer,
    const typename SplatPrimitive::ScreenBuffer splat_sbuffer,
#if IS_EVAL3D
    const float *__restrict__ viewmats, // [B, C, 4, 4]
    const float4 *__restrict__ intrins,  // [B, C, 4], fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const float4 *__restrict__ aabb,  // [..., N] projected 2D AABB (xmin,ymin,xmax,ymax)
#endif
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_width,
    const uint32_t tile_height,
    const int32_t *__restrict__ tile_offsets, // [I, tile_height, tile_width]
    const int32_t *__restrict__ flatten_ids,  // [n_isects]
    RenderOutput::Buffer render_colors, // [I, image_height, image_width, ...]
    float *__restrict__ render_Ts, // [I, image_height, image_width, 1]
    int32_t *__restrict__ last_ids, // [I, image_height, image_width]
    RenderOutput::Buffer render_colors2, // [I, image_height, image_width, ...]
    RenderOutput::Buffer render_distortions // [I, image_height, image_width, ...]
) {
    // Each block covers a tile on the image. In total there are
    // I * tile_height * tile_width blocks.
    dim3 threads = {NUM_THREADS, 1, 1};
    dim3 grid = {I, tile_height, tile_width};

#if IS_EVAL3D
    rasterize_to_pixels_eval3d_fwd_kernel<
#else
    rasterize_to_pixels_fwd_kernel<
#endif
        SplatPrimitive,
    #if IS_EVAL3D
        camera_model,
    #endif
        output_distortion
    ><<<grid, threads, 0, stream>>>(
        I, N, n_isects,
        gaussian_ids, splat_wbuffer, splat_sbuffer,
    #if IS_EVAL3D
        viewmats, intrins, dist_coeffs_buffer, aabb,
    #endif
        image_width, image_height, tile_width, tile_height, tile_offsets, flatten_ids,
        render_colors, render_Ts, last_ids,
        render_colors2, render_distortions
    );
}

#undef IS_EVAL3D
