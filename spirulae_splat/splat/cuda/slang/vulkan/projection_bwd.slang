// Fused projection backward (Vulkan backend), mirroring
// csrc/ProjectionBwd_kernel.cuh: one thread per (camera, gaussian)
// intersection; per-intersection VJP machinery in proj_bwd_common.slang,
// gradients scattered with CAS float atomics (multiple cameras hit the same
// splat).
//
// Variants (specialization constants, matching the CUDA template axes;
// 0..4 declared in proj_bwd_common.slang):
//   kCameraModel  = CameraModelType (0..3)
//   kShDegree     = 0..4  (already capped at max_sh_degree by the host)
//   kAntialiased  = 0 (Vanilla3DGS) | 1 (MipSplatting); ignored by 3DGUT
//   kShValueBits  = 0 (fp32 SH source) | 8 | 16 (word-based codec source)
//   kPacked       = 1: camera_ids/gaussian_ids index nnz intersections
//   kViewmatGrad  = 1: accumulate v_R / v_t into v_viewmats (per-thread CAS
//                   adds -- CUDA's labeled-partition warp reduce is a perf
//                   shortcut with identical, order-nondeterministic sums)

#define CudaDeviceExport ForceInline
#include "primitive_3dgs.slang"
#include "harmonics.slang"
#undef CudaDeviceExport
#include "sh_quant.slang"
#include "sh_vjp.slang"
#include "atomic_float.slang"
#include "proj_bwd_common.slang"

[SpecializationConstant]
const int kViewmatGrad = 0;

// ---------------------------------------------------------------------------
// 3DGS / MipSplatting: ScreenBuffer grads {xy, depth, conic, opac, rgb}.
// ---------------------------------------------------------------------------

struct ProjectionBwd3dgsParams {
    // splat world source values (layouts [N,3] [N,4] [N,3] [N,1] [N,3] [N,K3])
    float* means;
    float* quats;
    float* scales;
    float* opacities;
    float* features_dc;
    float* features_sh;
    // cameras
    float* viewmats;     // [C,16] row-major
    float4* intrins;     // [C]
    float* dist_coeffs;  // [C,10] or fallback zeros
    // packed intersection ids (kPacked; fallback otherwise)
    int* camera_ids;    // [nnz]
    int* gaussian_ids;  // [nnz]
    float4* aabb;       // [C,N] (non-packed cull; unused when kPacked)
    // screen-space gradients (inputs)
    float* vs_xy;
    float* vs_depth;
    float* vs_conic;
    float* vs_opac;
    float* vs_rgb;
    // world-space gradients (outputs, CAS accumulate)
    float* vw_means;
    float* vw_quats;
    float* vw_scales;
    float* vw_opacities;
    float* vw_dc;
    float* vw_sh;
    float* v_viewmats;  // [C,16] (kViewmatGrad; fallback otherwise)
    // SH value-quant source (kShValueBits 8/16)
    uint* sh_packed;
    float2* sh_bounds;
    int64_t sh_bounds_stride;
    uint32_t sh_stride_src;  // floats per splat in features_sh
    uint32_t sh_stride_v;    // floats per splat in vw_sh
    uint32_t C;
    uint32_t N;  // num_splats, or nnz when kPacked
    uint32_t width;
    uint32_t height;
    uint32_t wgs_per_row;
    uint32_t num_sh_buffer;
};

void _pb_run_3dgs(ProjectionBwd3dgsParams p, uint idx) {
    uint cid, gid;
    if (kPacked != 0) {
        cid = uint(p.camera_ids[idx]);
        gid = uint(p.gaussian_ids[idx]);
    } else {
        float4 bb = p.aabb[idx];
        if (bb.z <= bb.x || bb.w <= bb.y)
            return;
        cid = (idx / p.N) % p.C;
        gid = idx % p.N;
    }

    PbCamera cam = pb_load_camera(p.viewmats, p.intrins, p.dist_coeffs, cid);
    PbSplat s = pb_load_splat(p.means, p.quats, p.scales, p.opacities,
                              p.features_dc, gid);

    float2 v_xy = float2(p.vs_xy[2 * idx], p.vs_xy[2 * idx + 1]);
    float v_depth = p.vs_depth[idx];
    float3 v_conic = float3(p.vs_conic[3 * idx], p.vs_conic[3 * idx + 1],
                            p.vs_conic[3 * idx + 2]);
    float v_opac = p.vs_opac[idx];
    float3 v_rgb = float3(p.vs_rgb[3 * idx], p.vs_rgb[3 * idx + 1],
                          p.vs_rgb[3 * idx + 2]);

    PbGrads g = pb_grads_zero();
    pb_vjp_2d(s, cam, p.width, p.height, v_xy, v_depth, v_conic, v_opac, g);
    int64_t qbase = (int64_t)3 * (int64_t)p.num_sh_buffer * (int64_t)gid;
    pb_sh_vjp(s, cam, p.features_sh + gid * p.sh_stride_src, p.sh_packed,
              p.sh_bounds, qbase, p.sh_bounds_stride, v_rgb, g);

    pb_add3(p.vw_means + 3 * gid, g.v_mean);
    pb_add4(p.vw_quats + 4 * gid, g.v_quat);
    pb_add3(p.vw_scales + 3 * gid, g.v_scale);
    pb_add1(p.vw_opacities + gid, g.v_opacity);
    pb_add3(p.vw_dc + 3 * gid, g.v_dc);
    pb_scatter_sh(p.vw_sh + gid * p.sh_stride_v, g.v_sh, kShDegree);

    if (kViewmatGrad != 0)
        pb_scatter_viewmat(p.v_viewmats, cid, g.v_R, g.v_t);
}

[shader("compute")]
[numthreads(128, 1, 1)]
void projection_bwd_3dgs(uint3 wg: SV_GroupID, uint tid: SV_GroupThreadID,
                         uniform ProjectionBwd3dgsParams* pp) {
    ProjectionBwd3dgsParams p = *pp;
    uint idx = (wg.y * p.wgs_per_row + wg.x) * 128 + tid;
    uint total = (kPacked != 0) ? p.N : p.C * p.N;
    if (idx >= total) return;
    _pb_run_3dgs(p, idx);
}

// ---------------------------------------------------------------------------
// 3DGUT: ScreenBuffer grads {scale, opacity, rgb}.
// ---------------------------------------------------------------------------

struct ProjectionBwd3dgutParams {
    float* means;
    float* quats;
    float* scales;
    float* opacities;
    float* features_dc;
    float* features_sh;
    float* viewmats;
    float4* intrins;
    float* dist_coeffs;
    int* camera_ids;
    int* gaussian_ids;
    float4* aabb;
    // screen grads {scale, opac, rgb}
    float* vs_scale;
    float* vs_opac;
    float* vs_rgb;
    // world grads (outputs)
    float* vw_means;
    float* vw_quats;
    float* vw_scales;
    float* vw_opacities;
    float* vw_dc;
    float* vw_sh;
    float* v_viewmats;
    uint* sh_packed;
    float2* sh_bounds;
    int64_t sh_bounds_stride;
    uint32_t sh_stride_src;
    uint32_t sh_stride_v;
    uint32_t C;
    uint32_t N;
    uint32_t width;
    uint32_t height;
    uint32_t wgs_per_row;
    uint32_t num_sh_buffer;
};

void _pb_run_3dgut(ProjectionBwd3dgutParams p, uint idx) {
    uint cid, gid;
    if (kPacked != 0) {
        cid = uint(p.camera_ids[idx]);
        gid = uint(p.gaussian_ids[idx]);
    } else {
        float4 bb = p.aabb[idx];
        if (bb.z <= bb.x || bb.w <= bb.y)
            return;
        cid = (idx / p.N) % p.C;
        gid = idx % p.N;
    }

    PbCamera cam = pb_load_camera(p.viewmats, p.intrins, p.dist_coeffs, cid);
    PbSplat s = pb_load_splat(p.means, p.quats, p.scales, p.opacities,
                              p.features_dc, gid);

    float3 v_sscale = float3(p.vs_scale[3 * idx], p.vs_scale[3 * idx + 1],
                             p.vs_scale[3 * idx + 2]);
    float v_sopac = p.vs_opac[idx];
    float3 v_rgb = float3(p.vs_rgb[3 * idx], p.vs_rgb[3 * idx + 1],
                          p.vs_rgb[3 * idx + 2]);

    PbGrads g = pb_grads_zero();
    pb_vjp_3dgut(s, cam, p.width, p.height, v_sscale, v_sopac, g);
    int64_t qbase = (int64_t)3 * (int64_t)p.num_sh_buffer * (int64_t)gid;
    pb_sh_vjp(s, cam, p.features_sh + gid * p.sh_stride_src, p.sh_packed,
              p.sh_bounds, qbase, p.sh_bounds_stride, v_rgb, g);

    pb_add3(p.vw_means + 3 * gid, g.v_mean);
    pb_add4(p.vw_quats + 4 * gid, g.v_quat);
    pb_add3(p.vw_scales + 3 * gid, g.v_scale);
    pb_add1(p.vw_opacities + gid, g.v_opacity);
    pb_add3(p.vw_dc + 3 * gid, g.v_dc);
    pb_scatter_sh(p.vw_sh + gid * p.sh_stride_v, g.v_sh, kShDegree);

    if (kViewmatGrad != 0)
        pb_scatter_viewmat(p.v_viewmats, cid, g.v_R, g.v_t);
}

[shader("compute")]
[numthreads(128, 1, 1)]
void projection_bwd_3dgut(uint3 wg: SV_GroupID, uint tid: SV_GroupThreadID,
                          uniform ProjectionBwd3dgutParams* pp) {
    ProjectionBwd3dgutParams p = *pp;
    uint idx = (wg.y * p.wgs_per_row + wg.x) * 128 + tid;
    uint total = (kPacked != 0) ? p.N : p.C * p.N;
    if (idx >= total) return;
    _pb_run_3dgut(p, idx);
}
