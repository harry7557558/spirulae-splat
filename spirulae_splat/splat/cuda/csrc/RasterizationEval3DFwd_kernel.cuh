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

// One CUDA block per micro-tile (TILE_SIZE_X x TILE_SIZE_Y pixels, one thread
// per pixel). Binning is still done at the coarser macro-tile granularity, so a
// block reads its macro tile's gaussian range and culls the ones that miss its
// micro-tile. Per-micro-tile blocks let the scheduler spread a dense macro
// tile's work across many SMs (a single huge macro tile no longer stalls the
// whole launch on one SM) and keep shared-memory tiny for high occupancy.
inline constexpr uint32_t TILE_AREA = TILE_SIZE_X * TILE_SIZE_Y;

template<
    typename SplatPrimitive,
#if IS_EVAL3D
    CameraModelType camera_model,
#endif
    bool output_distortion,
    bool output_median
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
    RenderOutput::Buffer render_distortions, // [I, image_height, image_width, ...]
    float * render_median, // [I, image_height, image_width, 1], optional
    float2 * render_median_anchor // [I, image_height, image_width, 2], optional
                                  // near-anchor (z1, v1) for the moment-line median backward
) {
    auto block = cg::this_thread_block();
    uint32_t tid = threadIdx.x;  // 0 .. TILE_AREA-1, one per pixel
    int32_t image_id = blockIdx.x;

    // one block per micro-tile; recover the macro tile it belongs to for binning
    uint32_t mt_y = blockIdx.y;  // micro-tile row
    uint32_t mt_x = blockIdx.z;  // micro-tile col
    int32_t tile_id = (mt_y / MACRO_TILE_SIZE_Y) * tile_width + (mt_x / MACRO_TILE_SIZE_X);

    tile_offsets += image_id * tile_height * tile_width;
    render_Ts += image_id * image_height * image_width;
    last_ids += image_id * image_height * image_width;
    if constexpr (output_median) {
        render_median += image_id * image_height * image_width;
        render_median_anchor += image_id * image_height * image_width;
    }

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

    // this thread's pixel within the micro-tile
    uint32_t i = mt_y * TILE_SIZE_Y + tid / TILE_SIZE_X;
    uint32_t j = mt_x * TILE_SIZE_X + tid % TILE_SIZE_X;
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

    // current visibility left to render
    int32_t pix_id = i * image_width + j;
    int32_t pix_id_global = image_id * image_height * image_width + pix_id;
    float T = 1.0f;
    // index of most recent gaussian to write to this thread's pixel
    uint32_t cur_idx = 0;

    // median depth: 1/2-isosurface of a linear occupancy model occ(z) = 1 - T(z),
    // fit through the two stable endpoints -- the first and last *contributing*
    // splats (occupancy v = 1 - T after that splat). z_med solves the line at
    // occ = 1/2: z_med = z1 + (1/2 - v1)/c1, c1 = (v2 - v1)/(z2 - z1). This is a
    // smooth (rational) function of the contributing opacities and the two anchor
    // depths -- far fewer discontinuities than a per-splat 1/2-crossing, which
    // hops between adjacent splats. Stays 0 if occupancy never reaches 1/2.
    // Only the near anchor (z1, v1) is stored for backward (2 floats/pixel); the
    // far anchor is recovered there (= last_ids splat, occupancy local at its step).
    float med_z1 = 0.0f, med_v1 = 0.0f;  // near anchor (first contributing)
    float med_z2 = 0.0f, med_v2 = 0.0f;  // far  anchor (last  contributing)
    bool  med_has1 = false;

    RenderOutput pix_out = RenderOutput::zero();
    RenderOutput pix2_out = RenderOutput::zero();
    RenderOutput distortion_out = RenderOutput::zero();

    // gaussians overlapping this micro-tile's macro tile
    int32_t range_start = tile_offsets[tile_id];
    int32_t range_end =
        (image_id == I - 1) && (tile_id == tile_width * tile_height - 1)
            ? n_isects
            : tile_offsets[tile_id + 1];
    uint32_t num_batches =
        (range_end - range_start + TILE_AREA - 1) / TILE_AREA;

    // micro-tile pixel box, for the footprint cull
    const float mt_bx0 = (float)(mt_x * TILE_SIZE_X);
    const float mt_by0 = (float)(mt_y * TILE_SIZE_Y);

    __shared__ typename SplatPrimitive::FragmentFwd splat_batch[TILE_AREA];
    // whether the loaded gaussian's footprint reaches this micro-tile (uniform
    // across the whole block, so skipping it is divergence-free)
    __shared__ bool splat_hit[TILE_AREA];

    // collect and process batches of gaussians, one thread loads one gaussian
    for (uint32_t b = 0; b < num_batches; ++b) {
        // resync all threads before the next batch overwrites shared, and stop
        // once every pixel in the block is done
        if (__syncthreads_count(done) >= (int)TILE_AREA)
            break;

        uint32_t batch_start = range_start + TILE_AREA * b;
        uint32_t idx = batch_start + tid;
        bool hit = false;
        if (idx < (uint32_t)range_end) {
            int32_t g = flatten_ids[idx]; // flatten index in [I * N] or [nnz]
            uint32_t wi = gaussian_ids ? gaussian_ids[g] : g % N;
        #if IS_EVAL3D
            // 3dgut keeps the projected 2D conic in the screen "scale" channel
            // and the ellipse center is the AABB center. Cull here, before the
            // (expensive) 3D fragment load, matching the tile intersector.
            float opac = splat_sbuffer.opacities(g);
            if (opac > ALPHA_THRESHOLD) {
                float3 conic = splat_sbuffer.scales(g);
                float4 bb = aabb[g];
                float ecx = 0.5f * (bb.x + bb.z);
                float ecy = 0.5f * (bb.y + bb.w);
                float kk = 0.5f / __logf(opac / ALPHA_THRESHOLD);
                float3 inv_cov = { conic.x * kk, conic.y * kk, conic.z * kk };
                hit = ellipse_box_overlap_test(inv_cov,
                    mt_bx0 - ecx, mt_bx0 + TILE_SIZE_X - ecx,
                    mt_by0 - ecy, mt_by0 + TILE_SIZE_Y - ecy);
            }
            if (hit)
                splat_batch[tid].load(splat_wbuffer, splat_sbuffer, wi, g);
        #else
            // 2D splatting: the fragment is just the screen data, so load it
            // (cheap) and cull from its own conic / opacity.
            splat_batch[tid].load(splat_wbuffer, splat_sbuffer, wi, g);
            const typename SplatPrimitive::FragmentFwd& sp = splat_batch[tid];
            if (sp.opac > ALPHA_THRESHOLD) {
                float kk = 0.5f / __logf(sp.opac / ALPHA_THRESHOLD);
                float3 inv_cov = { sp.conic.x * kk, sp.conic.y * kk, sp.conic.z * kk };
                hit = ellipse_box_overlap_test(inv_cov,
                    mt_bx0 - sp.xy.x, mt_bx0 + TILE_SIZE_X - sp.xy.x,
                    mt_by0 - sp.xy.y, mt_by0 + TILE_SIZE_Y - sp.xy.y);
            }
        #endif
        }
        splat_hit[tid] = hit;

        // wait for other threads to collect the gaussians in the batch
        block.sync();

        // process gaussians in the current batch for this pixel
        uint32_t batch_size = min((uint32_t)TILE_AREA, (uint32_t)(range_end - batch_start));
        for (uint32_t t = 0; t < batch_size; ++t) {
            if (!splat_hit[t])  // gaussian misses this micro-tile (uniform skip)
                continue;
            if (done)
                continue;
            typename SplatPrimitive::FragmentFwd splat = splat_batch[t];
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

                // moment-line median: record the first and last contributing
                // splats as the line anchors (same set as cur_idx / last_ids, so
                // the far anchor is recoverable in backward via last_ids). The
                // occupancy after this splat is v = 1 - next_T.
                if constexpr (output_median) {
                    const float v = 1.0f - next_T;
                    if (!med_has1) { med_z1 = color.depth; med_v1 = v; med_has1 = true; }
                    med_z2 = color.depth; med_v2 = v;
                }
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
            } else done = true;

        }  // for (uint32_t t = 0; t < batch_size; ++t)
    }

    if (i < image_height && j < image_width) {
        render_Ts[pix_id] = T;
        if constexpr (output_median) {
            // 1/2-isosurface of the occupancy line through (z1,v1)-(z2,v2).
            // Valid only if occupancy reaches 1/2 (v2 >= 1/2) and the anchors are
            // distinct in depth and occupancy; else median stays 0. Clamp the
            // crossing to [z1, z2] so a faint near floater can't extrapolate it
            // wildly (the clamp is mirrored in the backward).
            float median_depth = 0.0f;
            float dz = med_z2 - med_z1;
            float dv = med_v2 - med_v1;
            if (med_has1 && med_v2 >= 0.5f && fabsf(dz) > 1e-9f && fabsf(dv) > 1e-9f) {
                float c1 = dv / dz;
                float z = med_z1 + (0.5f - med_v1) / c1;
                median_depth = fminf(fmaxf(z, med_z1), med_z2);
            } else if (med_has1 && med_v2 >= 0.5f) {
                // degenerate line (single contributor / coincident anchors)
                median_depth = med_z1;
            }
            render_median[pix_id] = median_depth;
            render_median_anchor[pix_id] = make_float2(med_z1, med_v1);
        }
        if constexpr (RenderOutput::has_depth(SplatPrimitive::pixelType))
            pix_out.depth /= fmaxf(1.0f - T, 1e-10f);
        pix_out.saveParamsToBuffer<SplatPrimitive::pixelType>(render_colors, pix_id_global);
        // index in bin of last gaussian in this pixel
        last_ids[pix_id] = static_cast<int32_t>(cur_idx);
        // distortion
        if constexpr (output_distortion) {
            pix2_out.saveParamsToBuffer<SplatPrimitive::pixelType>(render_colors2, pix_id_global);
            distortion_out.saveParamsToBuffer<SplatPrimitive::pixelType>(render_distortions, pix_id_global);
        }
    }

}

template<
    typename SplatPrimitive,
#if IS_EVAL3D
    CameraModelType camera_model,
#endif
    bool output_distortion,
    bool output_median
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
    RenderOutput::Buffer render_distortions, // [I, image_height, image_width, ...]
    float *__restrict__ render_median, // [I, image_height, image_width, 1], optional
    float2 *__restrict__ render_median_anchor // [I, image_height, image_width, 2], optional
) {
    // One block per micro-tile. The macro tile spans MACRO_TILE_SIZE_{X,Y}
    // micro-tiles, so the grid is the macro-tile grid scaled up accordingly.
    dim3 threads = {TILE_AREA, 1, 1};
    dim3 grid = {I, tile_height * MACRO_TILE_SIZE_Y, tile_width * MACRO_TILE_SIZE_X};

#if IS_EVAL3D
    rasterize_to_pixels_eval3d_fwd_kernel<
#else
    rasterize_to_pixels_fwd_kernel<
#endif
        SplatPrimitive,
    #if IS_EVAL3D
        camera_model,
    #endif
        output_distortion,
        output_median
    ><<<grid, threads, 0, stream>>>(
        I, N, n_isects,
        gaussian_ids, splat_wbuffer, splat_sbuffer,
    #if IS_EVAL3D
        viewmats, intrins, dist_coeffs_buffer, aabb,
    #endif
        image_width, image_height, tile_width, tile_height, tile_offsets, flatten_ids,
        render_colors, render_Ts, last_ids,
        render_colors2, render_distortions,
        render_median, render_median_anchor
    );
}

#undef IS_EVAL3D
