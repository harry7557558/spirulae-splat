// Grad-quant projection backward (Vulkan backend), mirroring
// csrc/ProjectionBwdQuantGrad_kernel.cuh: SPLAT-parallel (one thread per
// splat, 256/block); each thread loops its cameras in registers (no
// atomics), then folds the register gradient onto the running QUANTIZED
// grad accumulator via decode -> add -> block-reduce bounds -> encode.
// Codecs from optim_quant.slang (signed symmetric gradq; And+Or masked
// sub-word stores since per-splat cell runs need not be word-aligned).
//
// Spec constants 0..4 from proj_bwd_common.slang, plus:
//   kEval3d       = 1: 3DGUT screen-grad layout + VJP
//   kGqMask       = active-attribute bitmask (1 means, 2 quats, 4 scales,
//                   8 opac, 16 dc, 32 sh) -- folded so inactive attributes'
//                   loads never exist (llvmpipe speculated-load rule)
//   kWorldGradAdd = 1: add mean/quat/scale register grads onto the fp32
//                   world-grad buffer (3DGUT; one owner thread -> plain add)
//
// One thread per splat also for !inside tail threads: the camera loop is
// fully skipped for them (CUDA reads out-of-bounds aabb entries harmlessly
// from pool slack; BDA reads must stay in bounds), and every thread reaches
// every block reduce.

#define CudaDeviceExport ForceInline
#include "primitive_3dgs.slang"
#include "harmonics.slang"
#undef CudaDeviceExport
#include "sh_quant.slang"
#include "sh_vjp.slang"
#include "atomic_float.slang"
#include "proj_bwd_common.slang"
#include "optim_quant.slang"

[SpecializationConstant]
const int kEval3d = 0;
[SpecializationConstant]
const int kGqMask = 0;
[SpecializationConstant]
const int kWorldGradAdd = 0;

struct ProjectionQgradParams {
    // splat world source values
    float* means;
    float* quats;
    float* scales;
    float* opacities;
    float* features_dc;
    float* features_sh;
    // cameras
    float* viewmats;
    float4* intrins;
    float* dist_coeffs;
    // packed intersection ranges (kPacked; fallback otherwise)
    int* camera_id_bounds;  // [N+1]
    int* camera_ids;        // [nnz] ORIGINAL order
    int* perm;              // [nnz] sorted_pos -> original_pos
    float4* aabb;           // [C,N] or [nnz] (original order)
    // screen-space gradients: 3DGS xy/depth/conic/opac/rgb; 3DGUT
    // scale/opac/rgb in vs0..vs2 (vs3/vs4 fallback)
    float* vs0;
    float* vs1;
    float* vs2;
    float* vs3;
    float* vs4;
    // fp32 world grads (kWorldGradAdd: mean/quat/scale plain add)
    float* vw_means;
    float* vw_quats;
    float* vw_scales;
    // quantized grad accumulators (16-bit non-SH, 8-bit SH)
    uint* gq_means_packed;
    float2* gq_means_bounds;
    uint* gq_quats_packed;
    float2* gq_quats_bounds;
    uint* gq_scales_packed;
    float2* gq_scales_bounds;
    uint* gq_opac_packed;
    float2* gq_opac_bounds;
    uint* gq_dc_packed;
    float2* gq_dc_bounds;
    uint* gq_sh_packed;
    float2* gq_sh_bounds;
    // SH VALUE quant source (kShValueBits 8/16)
    uint* sh_value_packed;
    float2* sh_value_bounds;
    int64_t sh_value_bounds_stride;
    uint32_t sh_stride_src;
    uint32_t C;
    uint32_t N;  // num_splats (splat-parallel even when packed)
    uint32_t width;
    uint32_t height;
    uint32_t wgs_per_row;
    uint32_t num_sh_buffer;
    uint32_t _pad0;
};

// One non-SH attribute: decode running grad + add this sub-batch's
// register contribution, block-reduce fresh bounds, re-encode. The block
// reduce ALWAYS runs (uniform barriers); `active` is spec-folded.
[ForceInline]
void _qg_nonsh<let PRIMS : int>(bool active, uint* packed, float2* bounds,
                                uint block_id, uint gid, uint tid,
                                bool inside, bool blk_ok,
                                float vals_in[PRIMS]) {
    float vals[PRIMS];
    float2 mm_t = float2(1e30f, -1e30f);
    if (active && inside) {
        float2 old = bounds[block_id];
        [ForceUnroll]
        for (int i = 0; i < PRIMS; i++) {
            vals[i] = gradq_decode1<16>(packed, uint(PRIMS) * gid + uint(i),
                                        old) +
                      vals_in[i];
            mm_t.x = min(mm_t.x, vals[i]);
            mm_t.y = max(mm_t.y, vals[i]);
        }
    } else {
        [ForceUnroll]
        for (int i = 0; i < PRIMS; i++)
            vals[i] = 0.0f;
    }
    float2 mm = oq_block_reduce_mm2(mm_t, tid);
    if (active && tid == 0 && blk_ok)
        bounds[block_id] = mm;
    if (active && inside) {
        [ForceUnroll]
        for (int i = 0; i < PRIMS; i++)
            vkq_store_u16(packed, uint(PRIMS) * gid + uint(i),
                          gradq_encode1<16>(vals[i], mm));
    }
}

[shader("compute")]
[numthreads(256, 1, 1)]
void projection_qgrad(uint3 wg: SV_GroupID, uint tid: SV_GroupThreadID,
                      uniform ProjectionQgradParams* pp) {
    ProjectionQgradParams p = *pp;
    uint block_id = wg.y * p.wgs_per_row + wg.x;
    uint gid = block_id * 256 + tid;
    bool inside = gid < p.N;
    // Grid-fold tail: block_id can pad past the last quant block, so bounds
    // writes are blk_ok-gated (optimizer.slang convention).
    bool blk_ok = block_id < (p.N + 255u) / 256u;

    int cid_0 = 0, cid_1 = int(p.C);
    if (kPacked != 0) {
        cid_0 = inside ? p.camera_id_bounds[gid] : 0;
        cid_1 = inside ? p.camera_id_bounds[gid + 1] : 0;
    }

    PbGrads g = pb_grads_zero();
    if (inside) {
        PbSplat s = pb_load_splat(p.means, p.quats, p.scales, p.opacities,
                                  p.features_dc, gid);
        int64_t qbase =
            (int64_t)3 * (int64_t)p.num_sh_buffer * (int64_t)gid;
        for (int cid_t = cid_0; cid_t < cid_1; ++cid_t) {
            int idx = (kPacked != 0) ? p.perm[cid_t]
                                     : cid_t * int(p.N) + int(gid);
            int cid = (kPacked != 0) ? p.camera_ids[idx] : cid_t;
            float4 bb = p.aabb[idx];
            if (bb.z <= bb.x || bb.w <= bb.y)
                continue;
            PbCamera cam = pb_load_camera(p.viewmats, p.intrins,
                                          p.dist_coeffs, uint(cid));
            if (kEval3d != 0) {
                float3 v_ss = float3(p.vs0[3 * idx], p.vs0[3 * idx + 1],
                                     p.vs0[3 * idx + 2]);
                float v_so = p.vs1[idx];
                float3 v_rgb = float3(p.vs2[3 * idx], p.vs2[3 * idx + 1],
                                      p.vs2[3 * idx + 2]);
                pb_vjp_3dgut(s, cam, p.width, p.height, v_ss, v_so, g);
                pb_sh_vjp(s, cam, p.features_sh + gid * p.sh_stride_src,
                          p.sh_value_packed, p.sh_value_bounds, qbase,
                          p.sh_value_bounds_stride, v_rgb, g);
            } else {
                float2 v_xy = float2(p.vs0[2 * idx], p.vs0[2 * idx + 1]);
                float v_depth = p.vs1[idx];
                float3 v_conic = float3(p.vs2[3 * idx], p.vs2[3 * idx + 1],
                                        p.vs2[3 * idx + 2]);
                float v_opac = p.vs3[idx];
                float3 v_rgb = float3(p.vs4[3 * idx], p.vs4[3 * idx + 1],
                                      p.vs4[3 * idx + 2]);
                pb_vjp_2d(s, cam, p.width, p.height, v_xy, v_depth, v_conic,
                          v_opac, g);
                pb_sh_vjp(s, cam, p.features_sh + gid * p.sh_stride_src,
                          p.sh_value_packed, p.sh_value_bounds, qbase,
                          p.sh_value_bounds_stride, v_rgb, g);
            }
        }
    }

    // fp32 geometry writeback (3DGUT: add proj contribution onto the
    // raster-bwd world grads; one owner thread per gid -> race-free adds).
    if (kWorldGradAdd != 0 && inside) {
        [ForceUnroll]
        for (int k = 0; k < 3; k++)
            p.vw_means[3 * gid + k] += g.v_mean[k];
        [ForceUnroll]
        for (int k = 0; k < 4; k++)
            p.vw_quats[4 * gid + k] += g.v_quat[k];
        [ForceUnroll]
        for (int k = 0; k < 3; k++)
            p.vw_scales[3 * gid + k] += g.v_scale[k];
    }

    // Quantized attribute writebacks (each block-reduces unconditionally).
    {
        float vm[3] = {g.v_mean.x, g.v_mean.y, g.v_mean.z};
        _qg_nonsh<3>((kGqMask & 1) != 0, p.gq_means_packed,
                     p.gq_means_bounds, block_id, gid, tid, inside, blk_ok,
                     vm);
        float vq[4] = {g.v_quat.x, g.v_quat.y, g.v_quat.z, g.v_quat.w};
        _qg_nonsh<4>((kGqMask & 2) != 0, p.gq_quats_packed,
                     p.gq_quats_bounds, block_id, gid, tid, inside, blk_ok,
                     vq);
        float vs[3] = {g.v_scale.x, g.v_scale.y, g.v_scale.z};
        _qg_nonsh<3>((kGqMask & 4) != 0, p.gq_scales_packed,
                     p.gq_scales_bounds, block_id, gid, tid, inside, blk_ok,
                     vs);
        float vo[1] = {g.v_opacity};
        _qg_nonsh<1>((kGqMask & 8) != 0, p.gq_opac_packed, p.gq_opac_bounds,
                     block_id, gid, tid, inside, blk_ok, vo);
        float vd[3] = {g.v_dc.x, g.v_dc.y, g.v_dc.z};
        _qg_nonsh<3>((kGqMask & 16) != 0, p.gq_dc_packed, p.gq_dc_bounds,
                     block_id, gid, tid, inside, blk_ok, vd);
    }

    // SH grad (8-bit): two passes over the full 3*K cells so inactive
    // warmup bands (contribution 0) stay re-encoded against the fresh bound.
    if ((kGqMask & 32) != 0) {
        const int num_sh = kShDegree * (kShDegree + 2);
        const uint sh_base = 3u * p.num_sh_buffer * gid;
        const int total = int(3u * p.num_sh_buffer);
        const float2 sh_old =
            inside ? p.gq_sh_bounds[block_id] : float2(0.0f, 0.0f);
        float2 sh_mm = float2(1e30f, -1e30f);
        if (inside) {
            for (int c = 0; c < total; ++c) {
                float old = gradq_decode1<8>(p.gq_sh_packed,
                                             sh_base + uint(c), sh_old);
                float contrib =
                    (c < 3 * num_sh) ? g.v_sh[c / 3][c % 3] : 0.0f;
                float nv = old + contrib;
                sh_mm.x = min(sh_mm.x, nv);
                sh_mm.y = max(sh_mm.y, nv);
            }
        }
        sh_mm = oq_block_reduce_mm2(sh_mm, tid);
        if (tid == 0 && blk_ok)
            p.gq_sh_bounds[block_id] = sh_mm;
        if (inside) {
            for (int c = 0; c < total; ++c) {
                float old = gradq_decode1<8>(p.gq_sh_packed,
                                             sh_base + uint(c), sh_old);
                float contrib =
                    (c < 3 * num_sh) ? g.v_sh[c / 3][c % 3] : 0.0f;
                vkq_store_u8(p.gq_sh_packed, sh_base + uint(c),
                             gradq_encode1<8>(old + contrib, sh_mm));
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Packed-mode preprocessing helpers (mirror the FPBO/qgrad iota +
// camera_id_bounds kernels; the sort itself runs via backend::sort_pairs).
// ---------------------------------------------------------------------------

struct QgIotaParams {  // 16 bytes: pushed directly
    int* buf;
    uint32_t n;
    uint32_t wgs_per_row;
};

[shader("compute")]
[numthreads(256, 1, 1)]
void qgrad_iota(uint3 wg: SV_GroupID, uint tid: SV_GroupThreadID,
                uniform QgIotaParams p) {
    uint i = (wg.y * p.wgs_per_row + wg.x) * 256 + tid;
    if (i < p.n)
        p.buf[i] = int(i);
}

struct QgCamBoundsParams {  // 24 bytes: pushed directly
    int* gaussian_ids_sorted;  // [nnz], sorted by gid
    int* camera_id_bounds;     // [N+1]
    uint32_t nnz;
    uint32_t N;
    uint32_t wgs_per_row;
    uint32_t _pad0;
};

[shader("compute")]
[numthreads(256, 1, 1)]
void qgrad_camera_id_bounds(uint3 wg: SV_GroupID, uint tid: SV_GroupThreadID,
                            uniform QgCamBoundsParams p) {
    uint i = (wg.y * p.wgs_per_row + wg.x) * 256 + tid;
    if (i > p.nnz)
        return;
    int hi = (i == p.nnz) ? int(p.N) : p.gaussian_ids_sorted[i];
    int lo = (i == 0) ? 0 : p.gaussian_ids_sorted[i - 1] + 1;
    for (int k = lo; k <= hi; ++k)
        p.camera_id_bounds[k] = int(i);
}
