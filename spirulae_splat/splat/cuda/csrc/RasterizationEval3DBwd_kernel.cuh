#ifdef _KERNEL_CUH_INC
#error "Rasterization backward kernel must be included no more than once"
#endif

#define _KERNEL_CUH_INC

#ifndef NO_TORCH
#define NO_TORCH
#endif


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

#include <Common.cuh>


#ifndef IS_EVAL3D
#define IS_EVAL3D 1
#endif


constexpr uint SPLAT_BATCH_SIZE_NO_DISTORTION = WARP_SIZE;
constexpr uint SPLAT_BATCH_SIZE_WITH_DISTORTION = WARP_SIZE;

constexpr uint TILE_SIZE_DX = 8;
static_assert(TILE_SIZE_DX > 0 && TILE_SIZE_DX <= TILE_SIZE_X && TILE_SIZE_X % TILE_SIZE_DX == 0);

constexpr uint TILE_SIZE_DY = 8;
static_assert(TILE_SIZE_DY > 0 && TILE_SIZE_DY <= TILE_SIZE_Y && TILE_SIZE_Y % TILE_SIZE_DY == 0);


template <
    typename SplatPrimitive,
#if IS_EVAL3D
    CameraModelType camera_model,
#endif
    bool output_distortion,
#if IS_EVAL3D
    bool output_viewmat_grad,
#endif
    bool output_accum_weight,
    bool output_median
>
#if IS_EVAL3D
__global__ void rasterize_to_pixels_eval3d_bwd_kernel(
#else
__global__ void rasterize_to_pixels_bwd_kernel(
#endif
    const uint32_t I,
    const uint32_t N,   // zero if packed
    const uint32_t n_isects,
    // fwd inputs
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
    const int32_t *__restrict__ tile_offsets, // [..., tile_height, tile_width]
    const int32_t *__restrict__ flatten_ids,  // [n_isects]
    // fwd outputs
    const float *__restrict__ render_Ts,      // [..., image_height, image_width, 1]
    const int32_t *__restrict__ last_ids, // [..., image_height, image_width]
    const float2 *__restrict__ render_median_anchor, // [..., image_height, image_width, 2], optional
    RenderOutput::Buffer render_output_buffer,
    RenderOutput::Buffer render2_output_buffer,
    const float *__restrict__ loss_map_buffer,           // [..., image_height, image_width, 1]
    const float *__restrict__ accum_weight_map_buffer,           // [..., image_height, image_width, 1]
    // grad outputs
    RenderOutput::Buffer v_render_output_buffer,
    const float *__restrict__ v_render_Ts, // [..., image_height, image_width, 1]
    const float *__restrict__ v_median, // [..., image_height, image_width, 1], optional
    RenderOutput::Buffer v_distortions_output_buffer,
    // grad inputs
    typename SplatPrimitive::WorldBuffer v_splat_wbuffer,
    typename SplatPrimitive::ScreenBuffer v_splat_sbuffer,
    float *__restrict__ o_accum_weight
#if IS_EVAL3D
    ,
    float *__restrict__ v_viewmats // [B, C, 4, 4]
#endif
) {
    auto block = cg::this_thread_block();
    cg::thread_block_tile<WARP_SIZE> warp = cg::tiled_partition<WARP_SIZE>(block);
    uint32_t image_id = block.group_index().x;
    uint32_t tile_id = (block.group_index().y * TILE_SIZE_DY / (TILE_SIZE_Y * MACRO_TILE_SIZE_Y)) * tile_width +
        (block.group_index().z * TILE_SIZE_DX / (TILE_SIZE_X * MACRO_TILE_SIZE_X));
    uint32_t thread_id = block.thread_rank();

    tile_offsets += image_id * tile_height * tile_width;
    render_Ts += image_id * image_height * image_width;
    last_ids += image_id * image_height * image_width;
    v_render_Ts += image_id * image_height * image_width;

#if IS_EVAL3D
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
#endif

    constexpr uint BLOCK_SIZE = TILE_SIZE_DX * TILE_SIZE_DY;

    // load pixels
#if IS_EVAL3D
    __shared__ float4 shared_ray_d_pix_bin_final[BLOCK_SIZE];
#else
    __shared__ int32_t pix_bin_final[BLOCK_SIZE];
#endif
    __shared__ float2 pix_Ts_with_grad[BLOCK_SIZE];
    __shared__ RenderOutput v_pix_colors[BLOCK_SIZE];

    __shared__ RenderOutput pix_colors[output_distortion ? BLOCK_SIZE : 1];
    __shared__ RenderOutput pix2_colors[output_distortion ? BLOCK_SIZE : 1];
    __shared__ RenderOutput v_distortion_out[output_distortion ? BLOCK_SIZE : 1];

    __shared__ float accum_weight_map[output_accum_weight ? BLOCK_SIZE : 1];

    // median depth backward state (per pixel). v_median = dL/d(median depth).
    // The median is the 1/2-isosurface of the occupancy line through the near
    // anchor (z1,v1) [first contributing splat, loaded from render_median_anchor]
    // and the far anchor (z2,v2) [last contributing splat = last_ids, recovered
    // locally at its step]. The full gradient is exact for a fixed anchor pair:
    // depth grads to both anchors, occupancy grads via v_T injection at both
    // steps. Because the bwd sweep is back-to-front, the far anchor is reached
    // first; it computes both anchors' grads and stashes the near anchor's into
    // median_pend_* (consumed when the near-anchor splat is later reached, matched
    // by bit-identical depth == z1). median_pend_set != 0 means "pending armed".
    __shared__ float median_v[output_median ? BLOCK_SIZE : 1];
    __shared__ float median_z1[output_median ? BLOCK_SIZE : 1];
    __shared__ float median_v1[output_median ? BLOCK_SIZE : 1];
    __shared__ float median_pend_zgrad[output_median ? BLOCK_SIZE : 1];
    __shared__ float median_pend_Tgrad[output_median ? BLOCK_SIZE : 1];
    __shared__ int   median_pend_set[output_median ? BLOCK_SIZE : 1];

#if IS_EVAL3D
    float3 ray_o = SlangProjectionUtils::transform_ray_o(R, t);
    float3 total_v_ray_o = make_float3(0.0f, 0.0f, 0.0f);

    __shared__ float3 shared_v_ray_d[output_viewmat_grad ? BLOCK_SIZE : 1];
#endif

    constexpr uint SPLAT_BATCH_SIZE_CONST = output_distortion ?
        SPLAT_BATCH_SIZE_WITH_DISTORTION : SPLAT_BATCH_SIZE_NO_DISTORTION;

    #pragma unroll
    for (uint pix_id0 = 0; pix_id0 < BLOCK_SIZE; pix_id0 += SPLAT_BATCH_SIZE_CONST) {
        static_assert(BLOCK_SIZE % SPLAT_BATCH_SIZE_CONST == 0);
        uint pix_id_local = pix_id0 + thread_id;
        int pix_x = blockIdx.z * TILE_SIZE_DX + pix_id_local % TILE_SIZE_DX;
        int pix_y = blockIdx.y * TILE_SIZE_DY + pix_id_local / TILE_SIZE_DX;
        uint pix_id_global = pix_y * image_width + pix_x;
        bool inside = (pix_x < image_width && pix_y < image_height);
        
        int32_t bin_final = (inside ? last_ids[pix_id_global] : -1);
    #if IS_EVAL3D
        const float px = (float)pix_x + 0.5f;
        const float py = (float)pix_y + 0.5f;
        float3 raydir;
        inside &= SlangProjectionUtils::generate_ray(
            {(px-cx)/fx, (py-cy)/fy},
            (int)camera_model, dist_coeffs,
            &raydir
        );
        float3 ray_d = SlangProjectionUtils::transform_ray_d(R, raydir);  // mul(raydir, R);
        shared_ray_d_pix_bin_final[pix_id_local] =
            {ray_d.x, ray_d.y, ray_d.z, __int_as_float(bin_final)};
    #else
        pix_bin_final[pix_id_local] = bin_final;
    #endif

        uint pix_id_image_global = image_id * image_height * image_width + pix_id_global;
        float render_Ts_local = (inside ? render_Ts[pix_id_global] : 0.0f);
        float v_render_Ts_local = (inside ? v_render_Ts[pix_id_global] : 0.0f);
        RenderOutput v_render_output_local = (inside ?
            v_render_output_buffer.load<SplatPrimitive::pixelType>(pix_id_image_global)
             : RenderOutput::zero());
        if constexpr (RenderOutput::has_depth(SplatPrimitive::pixelType)) {
            float inv_alpha = 1.0f / fmaxf(1.0f - render_Ts_local, 1e-10f);
            float exp_depth = (inside ? render_output_buffer.depths[pix_id_image_global] : 0.0f);
            v_render_Ts_local += exp_depth * v_render_output_local.depth * inv_alpha;
            v_render_output_local.depth *= inv_alpha;
        }
        pix_Ts_with_grad[pix_id_local] = {render_Ts_local, v_render_Ts_local};
        v_pix_colors[pix_id_local] = v_render_output_local;

        if (output_distortion) {
            pix_colors[pix_id_local] = (inside ?
                render_output_buffer.load<SplatPrimitive::pixelType>(pix_id_image_global)
                : RenderOutput::zero());
            if constexpr (RenderOutput::has_depth(SplatPrimitive::pixelType)) {
                float alpha = fmaxf(1.0f - render_Ts_local, 1e-10f);
                pix_colors[pix_id_local].depth *= alpha;
            }
            pix2_colors[pix_id_local] = (inside ?
                render2_output_buffer.load<SplatPrimitive::pixelType>(pix_id_image_global)
                : RenderOutput::zero());
            v_distortion_out[pix_id_local] = (inside ?
                v_distortions_output_buffer.load<SplatPrimitive::pixelType>(pix_id_image_global)
                : RenderOutput::zero());
        }

    #if IS_EVAL3D
        if (output_viewmat_grad) {
            shared_v_ray_d[pix_id_local] = make_float3(0.0f, 0.0f, 0.0f);
        }
    #endif

        if (output_accum_weight) {
            accum_weight_map[pix_id_local] = (accum_weight_map_buffer != nullptr && inside) ?
                accum_weight_map_buffer[pix_id_image_global] : 0.0f;
        }

        if constexpr (output_median) {
            median_v[pix_id_local] = (v_median != nullptr && inside) ?
                v_median[pix_id_image_global] : 0.0f;
            float2 anchor = (render_median_anchor != nullptr && inside) ?
                render_median_anchor[pix_id_image_global] : make_float2(0.0f, 0.0f);
            median_z1[pix_id_local] = anchor.x;
            median_v1[pix_id_local] = anchor.y;
            median_pend_zgrad[pix_id_local] = 0.0f;
            median_pend_Tgrad[pix_id_local] = 0.0f;
            median_pend_set[pix_id_local] = 0;
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
    // ---- microtile survivor compaction ------------------------------------
    // The block owns one TILE_SIZE_DX x TILE_SIZE_DY sub-tile, but the macro-tile
    // range it reads is >80% gaussians that never touch this sub-tile. The
    // diagonal sweep length is fixed by the number of gaussians processed, so
    // per-thread culling alone gives no speedup (culled lanes just idle while the
    // warp still iterates). Instead we scan the range back-to-front, compact the
    // survivors (same ellipse-vs-box test as the forward sub-tile mask) into
    // shared, and only sweep those. surv[] is kept strictly back-to-front so the
    // per-pixel undo stays a continuous back-to-front walk across batches.
    constexpr int SURV_CAP = 8 * (int)WARP_SIZE;
    __shared__ int32_t surv[SURV_CAP];
    const float cull_bx0 = (float)(blockIdx.z * TILE_SIZE_DX);
    const float cull_by0 = (float)(blockIdx.y * TILE_SIZE_DY);

  for (int32_t scan_end = range_end - 1; scan_end >= range_start; ) {

    // fill surv[] with up to SURV_CAP survivors (back-to-front). The loop bounds
    // (s, count) are warp-uniform, so __ballot_sync stays convergent.
    int count = 0;
    int32_t s = scan_end;
    while (s >= range_start && count <= SURV_CAP - (int)WARP_SIZE) {
        int32_t gi = s - (int32_t)thread_id;
        bool ok = false;
        if (gi >= range_start) {
            int32_t sid = flatten_ids[gi];
        #if IS_EVAL3D
            // 3dgut keeps the 2D conic in the screen "scale" channel; center is
            // the AABB center (matches the tile intersector / forward mask).
            float c_opac = splat_sbuffer.opacities(sid);
            float3 c_conic = splat_sbuffer.scales(sid);
            float4 c_bb = aabb[sid];
            float c_ex = 0.5f * (c_bb.x + c_bb.z);
            float c_ey = 0.5f * (c_bb.y + c_bb.w);
        #else
            float c_opac = splat_sbuffer.opac(sid);
            float3 c_conic = splat_sbuffer.conic(sid);
            float2 c_xy = splat_sbuffer.xy(sid);
            float c_ex = c_xy.x, c_ey = c_xy.y;
        #endif
            if (c_opac > ALPHA_THRESHOLD) {
                float kk = 0.5f / __logf(c_opac / ALPHA_THRESHOLD);
                float3 inv_cov = { c_conic.x*kk, c_conic.y*kk, c_conic.z*kk };
                ok = ellipse_box_overlap_test(inv_cov,
                    cull_bx0 - c_ex, cull_bx0 + TILE_SIZE_DX - c_ex,
                    cull_by0 - c_ey, cull_by0 + TILE_SIZE_DY - c_ey);
            }
        }
        uint32_t bal = __ballot_sync(~0u, ok);
        int pos = count + __popc(bal & (((uint32_t)1 << thread_id) - 1));
        if (ok) surv[pos] = gi;
        count += __popc(bal);
        s -= (int32_t)WARP_SIZE;
    }
    scan_end = s;  // next segment resumes here
    __syncwarp();

    // sweep the compacted survivors in batches of SPLAT_BATCH_SIZE
    const int num_seg_batches = (count + (int)SPLAT_BATCH_SIZE - 1) / (int)SPLAT_BATCH_SIZE;
    for (int splat_b = 0; splat_b < num_seg_batches; ++splat_b) {
        const int batch_base = splat_b * (int)SPLAT_BATCH_SIZE;
        const int splat_batch_size = min((int)SPLAT_BATCH_SIZE, count - batch_base);
        const int surv_pos = batch_base + (int)thread_id;
        const bool active = (surv_pos < count);
        // thread 0 owns the back-most survivor of the batch (surv is back-to-front)
        const int32_t splat_idx = active ? surv[surv_pos] : (range_start - 1);

        // load splats
        typename SplatPrimitive::FragmentBwd splat;
        uint32_t splat_wid, splat_sid;
        if (active) {
            splat_sid = flatten_ids[splat_idx]; // flatten index in [I * N] or [nnz]
            splat_wid = gaussian_ids ? gaussian_ids[splat_sid] : splat_sid % N;
            splat.load(splat_wbuffer, splat_sbuffer, splat_wid, splat_sid);
        }

        // accumulate gradient
        typename SplatPrimitive::FragmentBwd v_splat = SplatPrimitive::FragmentBwd::zero(splat);
        float accum_weight = 0.0f;

        // at t=0, thread 0 (back-most survivor) undoes pixel 0; at t=1 it undoes
        // pixel 1 while thread 1 undoes pixel 0; etc. -> each pixel sees survivors
        // strictly back-to-front, continuous across batches and segments.
        for (int t = 0; t < splat_batch_size + BLOCK_SIZE - 1; ++t,
                (SPLAT_BATCH_SIZE_CONST <= WARP_SIZE ? __syncwarp() : __syncthreads())
        ) {
            int pix_id = t - thread_id;
            if (pix_id < 0 || pix_id >= BLOCK_SIZE || !active)
                continue;
        #if IS_EVAL3D
            float4 ray_d_pix_bin_final = shared_ray_d_pix_bin_final[pix_id];
            if (splat_idx > __float_as_int(ray_d_pix_bin_final.w))
                continue;
        #else
            int pix_global_x = blockIdx.z * TILE_SIZE_DX + pix_id % TILE_SIZE_DX;
            int pix_global_y = blockIdx.y * TILE_SIZE_DY + pix_id / TILE_SIZE_DX;
            const float px = (float)pix_global_x + 0.5f;
            const float py = (float)pix_global_y + 0.5f;
            if (splat_idx > pix_bin_final[pix_id])
                continue;
        #endif

            // evaluate alpha and early skip
        #if IS_EVAL3D
            float3 ray_d = {ray_d_pix_bin_final.x, ray_d_pix_bin_final.y, ray_d_pix_bin_final.z};
            float alpha = splat.evaluate_alpha(ray_o, ray_d);
            if (alpha <= ALPHA_THRESHOLD || dot(ray_d, ray_d) == 0.0f)
                continue;
        #else
            float alpha = splat.evaluate_alpha(px, py);
            if (alpha <= ALPHA_THRESHOLD)
                continue;
        #endif

        #if IS_EVAL3D
            RenderOutput color = splat.evaluate_color(ray_o, ray_d);
            if (color.depth <= 0.0f)
                continue;
        #else
            RenderOutput color = splat.evaluate_color(px, py);
        #endif

            // printf("t=%d, thread %u, splat %d (%u), pix_id %d, pix %d %d\n", t, thread_id, splat_idx-range_start, splat_gid, pix_id, pix_global_x, pix_global_y);

            // forward:
            // \left(c_{1},T_{1}\right)=\left(c_{0}+\alpha_{i}T_{0}c_{i},\ T_{0}\left(1-\alpha_{i}\right)\right)
            float T1 = pix_Ts_with_grad[pix_id].x;
            float v_T1 = pix_Ts_with_grad[pix_id].y;

            // undo pixel:
            // T_{0}=\frac{T_{1}}{1-\alpha_{i}}
            float ra = 1.0f / (1.0f - alpha);
            float T0 = T1 * ra;

            RenderOutput v_c = v_pix_colors[pix_id];

            // ---- median depth backward (moment-line) --------------------
            // z_med = clamp(z1 + (1/2 - v1)/c1, z1, z2), c1 = (v2-v1)/(z2-z1),
            // with anchors (z1,v1)=first and (z2,v2)=last contributing splat,
            // occupancy v = 1 - T(post). Exact gradient for the fixed anchor
            // pair: depth grads to z1,z2; occupancy grads as v_T injections
            // (dL/dT = -v_v since v = 1 - T) at each anchor's step.
            float v_depth_median = 0.0f;
            if constexpr (output_median) {
            #if IS_EVAL3D
                const int32_t bin_final = __float_as_int(ray_d_pix_bin_final.w);
            #else
                const int32_t bin_final = pix_bin_final[pix_id];
            #endif
                const float v_med = median_v[pix_id];
                // Far anchor: this is the last contributing splat (== last_ids).
                if (splat_idx == bin_final && v_med != 0.0f) {
                    const float z1 = median_z1[pix_id];
                    const float v1 = median_v1[pix_id];
                    const float z2 = color.depth;
                    const float v2 = 1.0f - T1;     // occupancy after far splat
                    const float dz = z2 - z1;
                    const float dv = v2 - v1;
                    if (z2 == z1) {
                        // degenerate (single contributor): median = z1 = this splat
                        v_depth_median += v_med;
                    } else if (v2 >= 0.5f && fabsf(dz) > 1e-9f && fabsf(dv) > 1e-9f) {
                        const float c1 = dv / dz;
                        const float z_raw = z1 + (0.5f - v1) / c1;
                        if (z_raw <= z1) {
                            // clamped low -> median = z1: all grad to near depth
                            median_pend_zgrad[pix_id] = v_med;
                            median_pend_Tgrad[pix_id] = 0.0f;
                            median_pend_set[pix_id] = 1;
                        } else if (z_raw >= z2) {
                            // clamped high -> median = z2: grad to far depth here
                            v_depth_median += v_med;
                        } else {
                            // interior: exact moment-line gradient
                            const float g = (0.5f - v1) / dv;     // z_med = z1 + g*dz
                            const float inv_dv2 = 1.0f / (dv * dv);
                            const float v_z1 = v_med * (1.0f - g);
                            const float v_z2 = v_med * g;
                            const float v_v1 = v_med * dz * (0.5f - v2) * inv_dv2;
                            const float v_v2 = v_med * (-dz) * (0.5f - v1) * inv_dv2;
                            // far anchor (this step): depth + occupancy (-> dL/dT)
                            v_depth_median += v_z2;
                            v_T1 += -v_v2;
                            // near anchor: stash, consumed at its (nearer) step
                            median_pend_zgrad[pix_id] = v_z1;
                            median_pend_Tgrad[pix_id] = -v_v1;
                            median_pend_set[pix_id] = 1;
                        }
                    }
                }
                // Near anchor: matched by bit-identical depth (z1 came from this
                // splat's color.depth in the forward). Reached after the far step.
                if (median_pend_set[pix_id] && color.depth == median_z1[pix_id]) {
                    v_depth_median += median_pend_zgrad[pix_id];
                    v_T1 += median_pend_Tgrad[pix_id];
                    median_pend_set[pix_id] = 0;
                }
            }

            // gradient to alpha:
            // \frac{dL}{d\alpha_{i}}
            // = \frac{dL}{dc_{1}}\frac{dc_{1}}{d\alpha_{i}}+\frac{dL}{dT_{1}}\frac{dT_{1}}{d\alpha_{i}}
            // = T_{0}\frac{dL}{dc_{1}}c_{i}-\frac{dL}{dT_{1}}T_{0}
            float v_alpha = T0 * color.dot(v_c) -v_T1 * T0;

            // gradient to color:
            // \frac{dL}{dc_{i}}
            // = \frac{dL}{dc_{1}}\frac{dc_{1}}{dc_{i}}
            // = \alpha_{i}T_{0}\frac{dL}{dc_{1}}
            RenderOutput v_color = v_c * (alpha * T0);
            v_color.depth += v_depth_median;  // median's direct depth grad to this splat

            // update pixel gradient:
            // \frac{dL}{dT_{0}}
            // = \frac{dL}{dc_{1}}\frac{dc_{1}}{dT_{0}}+\frac{dL}{dT_{1}}\frac{dT_{1}}{dT_{0}}
            // = \alpha_{i}\frac{dL}{dc_{1}}c_{i}+\frac{dL}{dT_{1}}\left(1-\alpha_{i}\right)
            float v_T0 = alpha * color.dot(v_c) + v_T1 * (1.0f - alpha);

            // distortion
            if (output_distortion) {
                // \left(d_{1},s_{1}\right)=\left(d_{0}+\alpha_{i}T_{0}\left(c_{i}^{2}\left(1-T_{0}\right)-2c_{i}c_{0}+s_{0}\right),\ s_{0}+\alpha_{i}T_{0}c_{i}^{2}\right)
                // \frac{dL}{ds}=0
                RenderOutput v_dist = v_distortion_out[pix_id];
                RenderOutput c0 =
                    pix_colors[pix_id] + color * -alpha * T0;
                RenderOutput s0 =
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
        #if IS_EVAL3D
            float3 v_ray_o_alpha = make_float3(0), v_ray_d_alpha = make_float3(0);
            float3 v_ray_o_color = make_float3(0), v_ray_d_color = make_float3(0);
        #endif
            {
            #if IS_EVAL3D
                splat.evaluate_alpha_vjp(ray_o, ray_d, v_alpha, v_splat, v_ray_o_alpha, v_ray_d_alpha);
                splat.evaluate_color_vjp(ray_o, ray_d, v_color, v_splat, v_ray_o_color, v_ray_d_color);
            #else
                float old_opac = v_splat.opac;
                splat.evaluate_alpha_vjp(px, py, v_alpha, v_splat);
                splat.evaluate_color_vjp(px, py, v_color, v_splat);
            #endif
            }

            if constexpr (output_accum_weight) {
                // accum_weight += accum_weight_map[pix_id] * alpha * T0;
                accum_weight = fmaxf(accum_weight, accum_weight_map[pix_id] * alpha * T0);
            }

        #if IS_EVAL3D
            if (output_viewmat_grad) {
                total_v_ray_o += v_ray_o_alpha + v_ray_o_color;
                shared_v_ray_d[pix_id] += v_ray_d_alpha + v_ray_d_color;
            }
        #endif

            // update pixel states
            pix_Ts_with_grad[pix_id] = { T0, v_T0 };
            // v_pix_colors remains the same

        }

        // accumulate gradient (only survivors reach here)
        if (active) {
            v_splat.atomicStore(v_splat_wbuffer, v_splat_sbuffer, splat_wid, splat_sid);
            if constexpr (output_accum_weight) {
                // atomicAddFVec(o_accum_weight + splat_wid, accum_weight);
                atomicMax(o_accum_weight + splat_wid, accum_weight);
            }
        }
    }  // for splat_b (within segment)

    __syncwarp();  // done reading surv[] before the next segment overwrites it
  }  // segment loop
#if IS_EVAL3D
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
#endif
}


template <
    typename SplatPrimitive,
#if IS_EVAL3D
    CameraModelType camera_model,
#endif
    bool output_distortion,
#if IS_EVAL3D
    bool output_viewmat_grad,
#endif
    bool output_accum_weight,
    bool output_median
>
#if IS_EVAL3D
void rasterize_to_pixels_eval3d_bwd_kernel_wrapper(
#else
void rasterize_to_pixels_bwd_kernel_wrapper(
#endif
    cudaStream_t stream,
    const uint32_t I,
    const uint32_t N,   // zero if packed
    const uint32_t n_isects,
    // fwd inputs
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
    const int32_t *__restrict__ tile_offsets, // [..., tile_height, tile_width]
    const int32_t *__restrict__ flatten_ids,  // [n_isects]
    // fwd outputs
    const float *__restrict__ render_Ts,      // [..., image_height, image_width, 1]
    const int32_t *__restrict__ last_ids, // [..., image_height, image_width]
    const float2 *__restrict__ render_median_anchor, // [..., image_height, image_width, 2], optional
    RenderOutput::Buffer render_output_buffer,
    RenderOutput::Buffer render2_output_buffer,
    const float *__restrict__ loss_map_buffer,           // [..., image_height, image_width, 1]
    const float *__restrict__ accum_weight_map_buffer,           // [..., image_height, image_width, 1]
    // grad outputs
    RenderOutput::Buffer v_render_output_buffer,
    const float *__restrict__ v_render_Ts, // [..., image_height, image_width, 1]
    const float *__restrict__ v_median, // [..., image_height, image_width, 1], optional
    RenderOutput::Buffer v_distortions_output_buffer,
    // grad inputs
    typename SplatPrimitive::WorldBuffer v_splat_wbuffer,
    typename SplatPrimitive::ScreenBuffer v_splat_sbuffer,
    float *__restrict__ o_accum_weight
#if IS_EVAL3D
    ,
    float *__restrict__ v_viewmats // [B, C, 4, 4]
#endif
) {
    dim3 threads = {output_distortion ?
        SPLAT_BATCH_SIZE_WITH_DISTORTION : SPLAT_BATCH_SIZE_NO_DISTORTION,
        1, 1};
    dim3 grid = {
        I,
        tile_height * (TILE_SIZE_Y * MACRO_TILE_SIZE_Y / TILE_SIZE_DY),
        tile_width * (TILE_SIZE_X * MACRO_TILE_SIZE_X / TILE_SIZE_DX)
    };

#if IS_EVAL3D
    rasterize_to_pixels_eval3d_bwd_kernel<
#else
    rasterize_to_pixels_bwd_kernel<
#endif
        SplatPrimitive,
    #if IS_EVAL3D
        camera_model,
    #endif
        output_distortion,
    #if IS_EVAL3D
        output_viewmat_grad,
    #endif
        output_accum_weight,
        output_median
    ><<<grid, threads, 0, stream>>>(
        I, N, n_isects,
        gaussian_ids, splat_wbuffer, splat_sbuffer,
    #if IS_EVAL3D
        viewmats, intrins, dist_coeffs_buffer, aabb,
    #endif
        image_width, image_height, tile_width, tile_height,
        tile_offsets, flatten_ids,
        render_Ts, last_ids, render_median_anchor,
        render_output_buffer, render2_output_buffer, loss_map_buffer, accum_weight_map_buffer,
        v_render_output_buffer, v_render_Ts,
        v_median,
        v_distortions_output_buffer,
        v_splat_wbuffer, v_splat_sbuffer,
        o_accum_weight
    #if IS_EVAL3D
        , v_viewmats
    #endif
    );

}
