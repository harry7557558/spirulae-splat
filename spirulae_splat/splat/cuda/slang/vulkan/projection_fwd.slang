// Fused projection forward (Vulkan backend), mirroring
// csrc/ProjectionFwd_kernel.cuh: one thread per (camera, gaussian), all math
// via the shared slang modules (primitive_3dgs.slang / harmonics.slang).
//
// Variants (specialization constants, matching the CUDA template axes):
//   kCameraModel  = CameraModelType (0..3)
//   kShDegree     = 0..4
//   kAntialiased  = 0 (Vanilla3DGS) | 1 (MipSplatting)
//   kShValueBits  = 0 (fp32) | 8 | 16 (the CUDA VALUE_BITS template axis;
//                   decode via sh_quant.slang's word-based q8/q16 loads)
//
// Params exceed the 128-byte push floor -> delivered via the params ring
// (push constant = one device address).

#define CudaDeviceExport ForceInline
#include "primitive_3dgs.slang"
#include "harmonics.slang"
#undef CudaDeviceExport
#include "sh_quant.slang"

[SpecializationConstant]
const int kCameraModel = 0;
[SpecializationConstant]
const int kShDegree = 0;
[SpecializationConstant]
const int kAntialiased = 0;
[SpecializationConstant]
const int kShValueBits = 0;

// SH -> rgb for one splat, gated exactly like _Base3DGS::project /
// Primitive3DGUT::project: fp32 features_sh when kShValueBits == 0, else the
// q8/q16 codec path (base = gid * 3 * num_sh_buffer, matching the encode
// side; bounds_stride is resolved by the host launcher).
[ForceInline]
float3 eval_splat_color(float3 mean, float3x3 R, float3 t, float3 dc,
                        float* features_sh, uint sh_stride, uint gid,
                        uint* sh_packed, float2* sh_bounds,
                        uint num_sh_buffer, int64_t sh_bounds_stride) {
    if (kShValueBits == 8 || kShValueBits == 16) {
        int64_t base = (int64_t)3 * (int64_t)num_sh_buffer * (int64_t)gid;
        if (kShValueBits == 8)
            return sh_to_color_qw<8>(kShDegree, mean, R, t, dc, sh_packed,
                                     sh_bounds, base, sh_bounds_stride);
        return sh_to_color_qw<16>(kShDegree, mean, R, t, dc, sh_packed,
                                  sh_bounds, base, sh_bounds_stride);
    }
    float3* sh = (float3*)(features_sh + gid * sh_stride);
    switch (kShDegree) {
    case 0: return sh0_to_color(mean, R, t, dc, sh);
    case 1: return sh1_to_color(mean, R, t, dc, sh);
    case 2: return sh2_to_color(mean, R, t, dc, sh);
    case 3: return sh3_to_color(mean, R, t, dc, sh);
    case 4: return sh4_to_color(mean, R, t, dc, sh);
    }
    return float3(0.0f);
}

struct ProjectionFwdParams {
    // WorldBuffer (TensorArray<6>) flattened; layouts [N,3] [N,4] [N,3]
    // [N,1] [N,3] [N, sh_stride].
    float* means;
    float* quats;
    float* scales;
    float* opacities;
    float* features_dc;
    float* features_sh;
    // cameras
    float* viewmats;      // [C,16] row-major
    float4* intrins;      // [C]
    float* dist_coeffs;   // [C,10] or null
    // outputs
    float4* out_aabb;     // [C,N]
    float* out_depths;    // [C,N]
    uint32_t* out_radii;  // [N] float bits; atomic max (radius >= 0)
    // ScreenBuffer (TensorArray<5>) flattened; null members are skipped.
    float* s_xy;
    float* s_depth;
    float* s_conic;
    float* s_opac;
    float* s_rgb;
    // SH value-quant (kShValueBits 8/16); null / 0 in the fp32 path.
    uint* sh_packed;
    float2* sh_bounds;
    int64_t sh_bounds_stride;  // cells per bound, resolved host-side
    uint32_t sh_stride;   // floats per splat in features_sh
    uint32_t C;
    uint32_t N;
    uint32_t width;
    uint32_t height;
    uint32_t wgs_per_row; // grid fold (flat id = wg.y * wgs_per_row*128 + x)
    uint32_t num_sh_buffer;
    uint32_t _pad0;
};

void run_projection<let cam: CameraModelType>(ProjectionFwdParams p,
                                              uint idx) {
    uint cid = (idx / p.N) % p.C;
    uint gid = idx % p.N;

    // camera
    float* vm = p.viewmats + cid * 16;
    float3x3 R;
    float3 t;
    [ForceUnroll]
    for (int r = 0; r < 3; r++) {
        R[r] = float3(vm[4 * r + 0], vm[4 * r + 1], vm[4 * r + 2]);
        t[r] = vm[4 * r + 3];
    }
    float4 intrin = p.intrins[cid];
    float dist[10];
    [ForceUnroll]
    for (int i = 0; i < 10; i++)
        dist[i] = p.dist_coeffs == nullptr ? 0.0f
                                           : p.dist_coeffs[10 * cid + i];

    // splat
    float3 mean = float3(p.means[3 * gid], p.means[3 * gid + 1],
                         p.means[3 * gid + 2]);
    float4 quat = float4(p.quats[4 * gid], p.quats[4 * gid + 1],
                         p.quats[4 * gid + 2], p.quats[4 * gid + 3]);
    float3 scale = float3(p.scales[3 * gid], p.scales[3 * gid + 1],
                          p.scales[3 * gid + 2]);
    float in_opacity = p.opacities[gid];
    float3 dc = float3(p.features_dc[3 * gid], p.features_dc[3 * gid + 1],
                       p.features_dc[3 * gid + 2]);

    float4 aabb;
    float sorting_depth, radius;
    float2 s_xy;
    float s_depth;
    float3 s_conic;
    float s_opac;
    projection_3dgs<cam, false, false>(
        kAntialiased != 0,
        mean, quat, scale, in_opacity, R, t,
        intrin.x, intrin.y, intrin.z, intrin.w, dist,
        p.width, p.height,
        aabb, sorting_depth, radius, s_xy, s_depth, s_conic, s_opac);

    // SH -> rgb, gated exactly like _Base3DGS::project.
    float3 rgb = float3(0.0f);
    if (aabb.z > aabb.x && aabb.w > aabb.y)
        rgb = eval_splat_color(mean, R, t, dc, p.features_sh, p.sh_stride,
                               gid, p.sh_packed, p.sh_bounds,
                               p.num_sh_buffer, p.sh_bounds_stride);

    // store (mirrors the kernel epilogue: clamp, degenerate cull, atomics)
    aabb.x = clamp(aabb.x, 0.0f, float(p.width) - 1.0f);
    aabb.y = clamp(aabb.y, 0.0f, float(p.height) - 1.0f);
    aabb.z = clamp(aabb.z, 0.0f, float(p.width) - 1.0f);
    aabb.w = clamp(aabb.w, 0.0f, float(p.height) - 1.0f);
    if (aabb.z - aabb.x > 1e-3f && aabb.w - aabb.y > 1e-3f) {
        if (p.s_xy != nullptr) {
            p.s_xy[2 * idx] = s_xy.x;
            p.s_xy[2 * idx + 1] = s_xy.y;
        }
        if (p.s_depth != nullptr) p.s_depth[idx] = s_depth;
        if (p.s_conic != nullptr) {
            p.s_conic[3 * idx] = s_conic.x;
            p.s_conic[3 * idx + 1] = s_conic.y;
            p.s_conic[3 * idx + 2] = s_conic.z;
        }
        if (p.s_opac != nullptr) p.s_opac[idx] = s_opac;
        if (p.s_rgb != nullptr) {
            p.s_rgb[3 * idx] = rgb.x;
            p.s_rgb[3 * idx + 1] = rgb.y;
            p.s_rgb[3 * idx + 2] = rgb.z;
        }
        p.out_aabb[idx] = aabb;
        p.out_depths[idx] = sorting_depth;
        // float atomicMax on non-negative floats == uint max on bits.
        InterlockedMax(p.out_radii[gid], asuint(radius));
    } else {
        p.out_aabb[idx] = float4(0.0f);
        p.out_depths[idx] = 0.0f;
    }
}

[shader("compute")]
[numthreads(128, 1, 1)]
void projection_fwd_3dgs(uint3 wg: SV_GroupID, uint tid: SV_GroupThreadID,
                         uniform ProjectionFwdParams* pp) {
    ProjectionFwdParams p = *pp;
    uint idx = (wg.y * p.wgs_per_row + wg.x) * 128 + tid;
    if (idx >= p.C * p.N) return;
    switch (kCameraModel) {
    case 0: run_projection<CameraModelType::PINHOLE>(p, idx); break;
    case 1: run_projection<CameraModelType::FISHEYE>(p, idx); break;
    case 2: run_projection<CameraModelType::EQUISOLID>(p, idx); break;
    case 3: run_projection<CameraModelType::EQUIRECTANGULAR>(p, idx); break;
    }
}

// ---------------------------------------------------------------------------
// 3DGUT (eval3d): projection_3dgs<cam, false, true>, antialiased always
// false; the ScreenBuffer is TensorArray<3> {scale, opacity, rgb} where the
// "conic" output lands in the scale slot (the eval3d rasterizer culls with
// the projected 2D conic). Mirrors Primitive3DGUT.cuh project().
// ---------------------------------------------------------------------------

struct Projection3dgutParams {
    float* means;
    float* quats;
    float* scales;
    float* opacities;
    float* features_dc;
    float* features_sh;
    float* viewmats;
    float4* intrins;
    float* dist_coeffs;
    float4* out_aabb;
    float* out_depths;
    uint32_t* out_radii;
    float* s_scale;  // [C*N,3]
    float* s_opac;   // [C*N,1]
    float* s_rgb;    // [C*N,3]
    uint* sh_packed;
    float2* sh_bounds;
    int64_t sh_bounds_stride;
    uint32_t sh_stride;
    uint32_t C;
    uint32_t N;
    uint32_t width;
    uint32_t height;
    uint32_t wgs_per_row;
    uint32_t num_sh_buffer;
    uint32_t _pad0;
};

void run_projection_3dgut<let cam: CameraModelType>(Projection3dgutParams p,
                                                    uint idx) {
    uint cid = (idx / p.N) % p.C;
    uint gid = idx % p.N;

    float* vm = p.viewmats + cid * 16;
    float3x3 R;
    float3 t;
    [ForceUnroll]
    for (int r = 0; r < 3; r++) {
        R[r] = float3(vm[4 * r + 0], vm[4 * r + 1], vm[4 * r + 2]);
        t[r] = vm[4 * r + 3];
    }
    float4 intrin = p.intrins[cid];
    float dist[10];
    [ForceUnroll]
    for (int i = 0; i < 10; i++)
        dist[i] = p.dist_coeffs == nullptr ? 0.0f
                                           : p.dist_coeffs[10 * cid + i];

    float3 mean = float3(p.means[3 * gid], p.means[3 * gid + 1],
                         p.means[3 * gid + 2]);
    float4 quat = float4(p.quats[4 * gid], p.quats[4 * gid + 1],
                         p.quats[4 * gid + 2], p.quats[4 * gid + 3]);
    float3 scale = float3(p.scales[3 * gid], p.scales[3 * gid + 1],
                          p.scales[3 * gid + 2]);
    float in_opacity = p.opacities[gid];
    float3 dc = float3(p.features_dc[3 * gid], p.features_dc[3 * gid + 1],
                       p.features_dc[3 * gid + 2]);

    float4 aabb;
    float sorting_depth, radius;
    float2 xy;
    float depth;
    float3 s_scale;
    float s_opac;
    projection_3dgs<cam, false, true>(
        false,  // antialiased: always false for 3DGUT
        mean, quat, scale, in_opacity, R, t,
        intrin.x, intrin.y, intrin.z, intrin.w, dist,
        p.width, p.height,
        aabb, sorting_depth, radius, xy, depth, s_scale, s_opac);

    float3 rgb = float3(0.0f);
    if (aabb.z > aabb.x && aabb.w > aabb.y)
        rgb = eval_splat_color(mean, R, t, dc, p.features_sh, p.sh_stride,
                               gid, p.sh_packed, p.sh_bounds,
                               p.num_sh_buffer, p.sh_bounds_stride);

    aabb.x = clamp(aabb.x, 0.0f, float(p.width) - 1.0f);
    aabb.y = clamp(aabb.y, 0.0f, float(p.height) - 1.0f);
    aabb.z = clamp(aabb.z, 0.0f, float(p.width) - 1.0f);
    aabb.w = clamp(aabb.w, 0.0f, float(p.height) - 1.0f);
    if (aabb.z - aabb.x > 1e-3f && aabb.w - aabb.y > 1e-3f) {
        if (p.s_scale != nullptr) {
            p.s_scale[3 * idx] = s_scale.x;
            p.s_scale[3 * idx + 1] = s_scale.y;
            p.s_scale[3 * idx + 2] = s_scale.z;
        }
        if (p.s_opac != nullptr) p.s_opac[idx] = s_opac;
        if (p.s_rgb != nullptr) {
            p.s_rgb[3 * idx] = rgb.x;
            p.s_rgb[3 * idx + 1] = rgb.y;
            p.s_rgb[3 * idx + 2] = rgb.z;
        }
        p.out_aabb[idx] = aabb;
        p.out_depths[idx] = sorting_depth;
        InterlockedMax(p.out_radii[gid], asuint(radius));
    } else {
        p.out_aabb[idx] = float4(0.0f);
        p.out_depths[idx] = 0.0f;
    }
}

[shader("compute")]
[numthreads(128, 1, 1)]
void projection_fwd_3dgut(uint3 wg: SV_GroupID, uint tid: SV_GroupThreadID,
                          uniform Projection3dgutParams* pp) {
    Projection3dgutParams p = *pp;
    uint idx = (wg.y * p.wgs_per_row + wg.x) * 128 + tid;
    if (idx >= p.C * p.N) return;
    switch (kCameraModel) {
    case 0: run_projection_3dgut<CameraModelType::PINHOLE>(p, idx); break;
    case 1: run_projection_3dgut<CameraModelType::FISHEYE>(p, idx); break;
    case 2: run_projection_3dgut<CameraModelType::EQUISOLID>(p, idx); break;
    case 3: run_projection_3dgut<CameraModelType::EQUIRECTANGULAR>(p, idx);
        break;
    }
}

// ---------------------------------------------------------------------------
// Packed projection (two-pass compaction, ProjectionPackedFwd_kernel.cuh).
// Pass 1 writes an int32 0/1 mask (the CUDA bool mask would need 8-bit
// stores) which the host scans with backend::inclusive_sum<int32>; pass 2
// recomputes the projection and scatters into nnz-compacted outputs.
// Spec constant 3 (kEval3d) selects the 3DGUT projection in the mask pass.
// ---------------------------------------------------------------------------

[SpecializationConstant]
const int kEval3d = 0;

struct PackedMaskParams {  // fits in the 128-byte push floor: pushed direct
    float* means;
    float* quats;
    float* scales;
    float* opacities;
    float* viewmats;
    float4* intrins;
    float* dist_coeffs;
    int* out_mask;  // [C*N] 0/1
    uint32_t C;
    uint32_t N;
    uint32_t width;
    uint32_t height;
    uint32_t wgs_per_row;
};

void run_packed_mask<let cam: CameraModelType, let eval3d: bool>(
    PackedMaskParams p, uint idx) {
    uint cid = (idx / p.N) % p.C;
    uint gid = idx % p.N;

    float* vm = p.viewmats + cid * 16;
    float3x3 R;
    float3 t;
    [ForceUnroll]
    for (int r = 0; r < 3; r++) {
        R[r] = float3(vm[4 * r + 0], vm[4 * r + 1], vm[4 * r + 2]);
        t[r] = vm[4 * r + 3];
    }
    float4 intrin = p.intrins[cid];
    float dist[10];
    [ForceUnroll]
    for (int i = 0; i < 10; i++)
        dist[i] = p.dist_coeffs == nullptr ? 0.0f
                                           : p.dist_coeffs[10 * cid + i];

    float3 mean = float3(p.means[3 * gid], p.means[3 * gid + 1],
                         p.means[3 * gid + 2]);
    float4 quat = float4(p.quats[4 * gid], p.quats[4 * gid + 1],
                         p.quats[4 * gid + 2], p.quats[4 * gid + 3]);
    float3 scale = float3(p.scales[3 * gid], p.scales[3 * gid + 1],
                          p.scales[3 * gid + 2]);
    float in_opacity = p.opacities[gid];

    float4 aabb;
    float sorting_depth, radius;
    float2 xy;
    float depth;
    float3 conic;
    float opac;
    // The mask only consumes the AABB; SH evaluation is skipped (the CUDA
    // mask kernel computes and discards it).
    projection_3dgs<cam, false, eval3d>(
        !eval3d && kAntialiased != 0,
        mean, quat, scale, in_opacity, R, t,
        intrin.x, intrin.y, intrin.z, intrin.w, dist,
        p.width, p.height,
        aabb, sorting_depth, radius, xy, depth, conic, opac);

    aabb.x = clamp(aabb.x, 0.0f, float(p.width) - 1.0f);
    aabb.y = clamp(aabb.y, 0.0f, float(p.height) - 1.0f);
    aabb.z = clamp(aabb.z, 0.0f, float(p.width) - 1.0f);
    aabb.w = clamp(aabb.w, 0.0f, float(p.height) - 1.0f);
    p.out_mask[idx] =
        (aabb.z - aabb.x > 1e-3f && aabb.w - aabb.y > 1e-3f) ? 1 : 0;
}

[shader("compute")]
[numthreads(128, 1, 1)]
void projection_packed_mask(uint3 wg: SV_GroupID, uint tid: SV_GroupThreadID,
                            uniform PackedMaskParams p) {
    uint idx = (wg.y * p.wgs_per_row + wg.x) * 128 + tid;
    if (idx >= p.C * p.N) return;
    if (kEval3d != 0) {
        switch (kCameraModel) {
        case 0: run_packed_mask<CameraModelType::PINHOLE, true>(p, idx); break;
        case 1: run_packed_mask<CameraModelType::FISHEYE, true>(p, idx); break;
        case 2: run_packed_mask<CameraModelType::EQUISOLID, true>(p, idx); break;
        case 3: run_packed_mask<CameraModelType::EQUIRECTANGULAR, true>(p, idx); break;
        }
    } else {
        switch (kCameraModel) {
        case 0: run_packed_mask<CameraModelType::PINHOLE, false>(p, idx); break;
        case 1: run_packed_mask<CameraModelType::FISHEYE, false>(p, idx); break;
        case 2: run_packed_mask<CameraModelType::EQUISOLID, false>(p, idx); break;
        case 3: run_packed_mask<CameraModelType::EQUIRECTANGULAR, false>(p, idx); break;
        }
    }
}

struct PackedFwdParams {
    float* means;
    float* quats;
    float* scales;
    float* opacities;
    float* features_dc;
    float* features_sh;
    float* viewmats;
    float4* intrins;
    float* dist_coeffs;
    int* mask_scan;         // [C*N] inclusive scan of the mask
    int* out_camera_ids;    // [nnz]
    int* out_gaussian_ids;  // [nnz]
    float4* out_aabb;       // [nnz]
    float* out_depths;      // [nnz]
    uint32_t* out_radii;    // [N] float bits, atomic max
    float* s0;              // screen tensors: 3DGS xy/depth/conic/opac/rgb,
    float* s1;              // 3DGUT scale/opac/rgb (s3, s4 null)
    float* s2;
    float* s3;
    float* s4;
    uint* sh_packed;
    float2* sh_bounds;
    int64_t sh_bounds_stride;
    uint32_t sh_stride;
    uint32_t C;
    uint32_t N;
    uint32_t width;
    uint32_t height;
    uint32_t wgs_per_row;
    uint32_t num_sh_buffer;
    uint32_t _pad0;
};

void run_packed_fwd<let cam: CameraModelType, let eval3d: bool>(
    PackedFwdParams p, uint idx) {
    int out_idx = idx == 0 ? 0 : p.mask_scan[idx - 1];
    if (p.mask_scan[idx] == out_idx) return;

    uint cid = (idx / p.N) % p.C;
    uint gid = idx % p.N;

    float* vm = p.viewmats + cid * 16;
    float3x3 R;
    float3 t;
    [ForceUnroll]
    for (int r = 0; r < 3; r++) {
        R[r] = float3(vm[4 * r + 0], vm[4 * r + 1], vm[4 * r + 2]);
        t[r] = vm[4 * r + 3];
    }
    float4 intrin = p.intrins[cid];
    float dist[10];
    [ForceUnroll]
    for (int i = 0; i < 10; i++)
        dist[i] = p.dist_coeffs == nullptr ? 0.0f
                                           : p.dist_coeffs[10 * cid + i];

    float3 mean = float3(p.means[3 * gid], p.means[3 * gid + 1],
                         p.means[3 * gid + 2]);
    float4 quat = float4(p.quats[4 * gid], p.quats[4 * gid + 1],
                         p.quats[4 * gid + 2], p.quats[4 * gid + 3]);
    float3 scale = float3(p.scales[3 * gid], p.scales[3 * gid + 1],
                          p.scales[3 * gid + 2]);
    float in_opacity = p.opacities[gid];
    float3 dc = float3(p.features_dc[3 * gid], p.features_dc[3 * gid + 1],
                       p.features_dc[3 * gid + 2]);

    float4 aabb;
    float sorting_depth, radius;
    float2 xy;
    float depth;
    float3 conic;  // eval3d: carries the culling scale (see 3DGUT notes)
    float opac;
    projection_3dgs<cam, false, eval3d>(
        !eval3d && kAntialiased != 0,
        mean, quat, scale, in_opacity, R, t,
        intrin.x, intrin.y, intrin.z, intrin.w, dist,
        p.width, p.height,
        aabb, sorting_depth, radius, xy, depth, conic, opac);

    float3 rgb = float3(0.0f);
    if (aabb.z > aabb.x && aabb.w > aabb.y)
        rgb = eval_splat_color(mean, R, t, dc, p.features_sh, p.sh_stride,
                               gid, p.sh_packed, p.sh_bounds,
                               p.num_sh_buffer, p.sh_bounds_stride);

    p.out_camera_ids[out_idx] = int(cid);
    p.out_gaussian_ids[out_idx] = int(gid);
    aabb.x = clamp(aabb.x, 0.0f, float(p.width) - 1.0f);
    aabb.y = clamp(aabb.y, 0.0f, float(p.height) - 1.0f);
    aabb.z = clamp(aabb.z, 0.0f, float(p.width) - 1.0f);
    aabb.w = clamp(aabb.w, 0.0f, float(p.height) - 1.0f);
    p.out_aabb[out_idx] = aabb;
    p.out_depths[out_idx] = sorting_depth;
    if (eval3d) {  // Screen = {scale, opacity, rgb}
        if (p.s0 != nullptr) {
            p.s0[3 * out_idx] = conic.x;
            p.s0[3 * out_idx + 1] = conic.y;
            p.s0[3 * out_idx + 2] = conic.z;
        }
        if (p.s1 != nullptr) p.s1[out_idx] = opac;
        if (p.s2 != nullptr) {
            p.s2[3 * out_idx] = rgb.x;
            p.s2[3 * out_idx + 1] = rgb.y;
            p.s2[3 * out_idx + 2] = rgb.z;
        }
    } else {  // Screen = {xy, depth, conic, opac, rgb}
        if (p.s0 != nullptr) {
            p.s0[2 * out_idx] = xy.x;
            p.s0[2 * out_idx + 1] = xy.y;
        }
        if (p.s1 != nullptr) p.s1[out_idx] = depth;
        if (p.s2 != nullptr) {
            p.s2[3 * out_idx] = conic.x;
            p.s2[3 * out_idx + 1] = conic.y;
            p.s2[3 * out_idx + 2] = conic.z;
        }
        if (p.s3 != nullptr) p.s3[out_idx] = opac;
        if (p.s4 != nullptr) {
            p.s4[3 * out_idx] = rgb.x;
            p.s4[3 * out_idx + 1] = rgb.y;
            p.s4[3 * out_idx + 2] = rgb.z;
        }
    }
    InterlockedMax(p.out_radii[gid], asuint(radius));
}

[shader("compute")]
[numthreads(128, 1, 1)]
void projection_packed_fwd(uint3 wg: SV_GroupID, uint tid: SV_GroupThreadID,
                           uniform PackedFwdParams* pp) {
    PackedFwdParams p = *pp;
    uint idx = (wg.y * p.wgs_per_row + wg.x) * 128 + tid;
    if (idx >= p.C * p.N) return;
    if (kEval3d != 0) {
        switch (kCameraModel) {
        case 0: run_packed_fwd<CameraModelType::PINHOLE, true>(p, idx); break;
        case 1: run_packed_fwd<CameraModelType::FISHEYE, true>(p, idx); break;
        case 2: run_packed_fwd<CameraModelType::EQUISOLID, true>(p, idx); break;
        case 3: run_packed_fwd<CameraModelType::EQUIRECTANGULAR, true>(p, idx); break;
        }
    } else {
        switch (kCameraModel) {
        case 0: run_packed_fwd<CameraModelType::PINHOLE, false>(p, idx); break;
        case 1: run_packed_fwd<CameraModelType::FISHEYE, false>(p, idx); break;
        case 2: run_packed_fwd<CameraModelType::EQUISOLID, false>(p, idx); break;
        case 3: run_packed_fwd<CameraModelType::EQUIRECTANGULAR, false>(p, idx); break;
        }
    }
}
