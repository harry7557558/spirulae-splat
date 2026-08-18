// Normal (3-channel axis-angle rotation) bilateral-grid samplers (Vulkan
// backend), mirroring BilagridNormalSampleFwd_kernel.cuh and
// BilagridNormalSampleBwdV1_kernel.cuh (uniform + PATCHED). The axis-angle
// rotate forward/backward are hand-ports of BilagridAxisAngle.cuh — the CUDA
// backward is a hand-written Rodrigues-form gradient (NOT autodiff of the
// forward), so it is mirrored term by term.
//
// Spec constants (IDs by declaration order — dispatches must pass ALL):
//   0: kValueQuant  1: kPatched

#include "bilagrid_common.slang"

[SpecializationConstant]
const int kValueQuant = 0;
[SpecializationConstant]
const int kPatched = 0;

// ------------------------- axis-angle rotate (fwd) -------------------------
[ForceInline]
float3 aa_rotate_fwd(float3 w, float3 v) {
    const float eps = 1e-5f;
    float theta2 = dot(w, w);
    float theta = sqrt(theta2 + 1e-20f);

    float c, s, t;
    float3 u;
    if (theta < eps) {
        float inv6 = 1.0f / 6.0f;
        c = 1.0f - 0.5f * theta2;
        s = theta - inv6 * theta * theta2;
        t = 0.5f * theta2;
        float inv_theta = 1.0f / (theta + 1e-20f);
        u = w * inv_theta;
    } else {
        c = cos(theta);
        s = sin(theta);
        t = 1.0f - c;
        float inv_theta = 1.0f / theta;
        u = w * inv_theta;
    }

    float3 uv = cross(u, v);
    float uvd = dot(u, v);
    return v * c + uv * s + u * (uvd * t);
}

// ------------------------- axis-angle rotate (bwd) -------------------------
[ForceInline]
void aa_rotate_bwd(float3 w, float3 v, float3 grad_out, out float3 grad_w,
                   out float3 grad_v) {
    const float eps = 1e-5f;
    float theta2 = dot(w, w);
    float theta = sqrt(theta2 + 1e-20f);

    float A, B, dA_dtheta, dB_dtheta;
    if (theta < eps) {
        float t2 = theta2;
        float t4 = t2 * t2;
        A = 1.0f - t2 * (1.0f / 6.0f) + t4 * (1.0f / 120.0f);
        B = 0.5f - t2 * (1.0f / 24.0f) + t4 * (1.0f / 720.0f);
        dA_dtheta = -theta * (1.0f / 3.0f) + theta * t2 * (1.0f / 30.0f);
        dB_dtheta = -theta * (1.0f / 12.0f) + theta * t2 * (1.0f / 180.0f);
    } else {
        float s = sin(theta);
        float c = cos(theta);
        A = s / theta;
        B = (1.0f - c) / theta2;
        dA_dtheta = (theta * c - s) / theta2;
        dB_dtheta = (theta * s - 2.0f * (1.0f - c)) / (theta2 * theta);
    }

    float3 c1 = cross(w, v);
    float3 c2 = cross(w, c1);

    float g_c1 = dot(grad_out, c1);
    float g_c2 = dot(grad_out, c2);
    float inv_theta = 1.0f / (theta + 1e-20f);
    float alpha = (dA_dtheta * g_c1 + dB_dtheta * g_c2) * inv_theta;

    float w_dot_v = dot(w, v);
    float g_dot_v = dot(grad_out, v);

    float3 term_theta = w * alpha;
    float3 term_a = cross(v, grad_out) * A;
    float3 term_b = cross(c1, grad_out) * B;
    float3 term_c = (grad_out * w_dot_v + w * (-g_dot_v)) * B;
    grad_w = term_theta + (term_a + (term_b + term_c));

    grad_v = aa_rotate_fwd(w * (-1.0f), grad_out);
}

// ---------------------------------------------------------------------------
// uniform / patched forward
// ---------------------------------------------------------------------------

struct BnUniformParams {
    float* fp32;
    uint* q16;
    float2* vbounds;
    float* normal_in;   // [N,(m,)h,w,3]
    float* output;      // fwd out / bwd v_output
    int* offsets;
    int* grid_indices;
    int N, L, H, W, m, h, w, h0, w0;
    int has_grid_indices;
    uint total;
    uint wgs_per_row;
};

[ForceInline]
void _bn_pix(int idx, BnUniformParams p, out int ni, out int g_id,
             out int g_off3, out float gx, out float gy) {
    int tmp = idx;
    int wi = tmp % p.w; tmp /= p.w;
    int hi = tmp % p.h; tmp /= p.h;
    int mi = 0;
    if (kPatched != 0) { mi = tmp % p.m; tmp /= p.m; }
    ni = tmp;
    if (kPatched != 0) {
        g_id = ni;
        int o = (ni * p.m + mi) * 2;
        gx = float(wi + p.offsets[o + 0]) / float(p.w0 - 1);
        gy = float(hi + p.offsets[o + 1]) / float(p.h0 - 1);
        g_off3 = (((ni * p.m + mi) * p.h + hi) * p.w + wi) * 3;
    } else {
        g_id = p.has_grid_indices != 0 ? p.grid_indices[ni] : ni;
        gx = float(wi) / float(p.w - 1);
        gy = float(hi) / float(p.h - 1);
        g_off3 = ((ni * p.h + hi) * p.w + wi) * 3;
    }
}

[shader("compute")]
[numthreads(256, 1, 1)]
void bilagrid_normal_fwd(uint3 wg: SV_GroupID, uint tid: SV_GroupThreadID,
                         uniform BnUniformParams p) {
    uint idx = (wg.y * p.wgs_per_row + wg.x) * 256 + tid;
    if (idx >= p.total) return;
    int ni, g_id, g_off3;
    float gx, gy;
    _bn_pix(int(idx), p, ni, g_id, g_off3, gx, gy);
    BgReader r = {p.fp32, p.q16, p.vbounds};

    float sr = p.normal_in[g_off3 + 0];
    float sg = p.normal_in[g_off3 + 1];
    float sb = p.normal_in[g_off3 + 2];

    // "No normal here" must survive the grid; see the CUDA kernel.
    if (!(sr + sg + sb > -2.366f && sr * sr + sg * sg + sb * sb > 0.25f)) {
        p.output[g_off3 + 0] = sr;
        p.output[g_off3 + 1] = sg;
        p.output[g_off3 + 2] = sb;
        return;
    }

    float inv_norm = 1.0f / sqrt(sr * sr + sg * sg + sb * sb + 1e-20f);

    float gz = 0.5f + 0.5f * sb * inv_norm;
    gz = min(max(gz, 0.0f), 1.0f);
    float x = gx * float(p.W - 1);
    float y = gy * float(p.H - 1);
    float z = gz * float(p.L - 1);

    int base = g_id * p.L * p.H * p.W * 3;
    BgCells c = bg_cells_uniform(base, p.L, p.H, p.W, 3, x, y, z);

    float3 axis_angle;
    axis_angle.x = bg_trilerp(kValueQuant != 0, r, c, 0);
    axis_angle.y = bg_trilerp(kValueQuant != 0, r, c, 1);
    axis_angle.z = bg_trilerp(kValueQuant != 0, r, c, 2);

    float3 normal = float3(sr * inv_norm, sg * inv_norm, sb * inv_norm);
    normal = aa_rotate_fwd(axis_angle, normal);

    p.output[g_off3 + 0] = bg_finite(normal.x) ? normal.x : 0.5f;
    p.output[g_off3 + 1] = bg_finite(normal.y) ? normal.y : 0.5f;
    p.output[g_off3 + 2] = bg_finite(normal.z) ? normal.z : 0.5f;
}

// ---------------------------------------------------------------------------
// backward v1: grid-grad kernel (8x8x1, 3-channel shmem preload)
// ---------------------------------------------------------------------------

static const int kBnFootX = kBgV1BlockX + 2;
static const int kBnFootY = kBgV1BlockY + 2;
static const int kBnFootZ = 3;
static const int kBnChans = 3;
static const int kBnCells = kBnFootX * kBnFootY * kBnFootZ;

groupshared float _bn_shmem[kBnFootZ * kBnFootY * kBnFootX * kBnChans];

[ForceInline]
int _bn_sh_idx(int lz, int ly, int lx, int ci) {
    return ((lz * kBnFootY + ly) * kBnFootX + lx) * kBnChans + ci;
}

struct BnV1GridParams {
    float* fp32;
    uint* q16;
    float2* vbounds;
    float* normal_in;
    float* v_output;
    float* v_bilagrid;
    int* offsets;
    int* grid_indices;
    int N, L, H, W, m, h, w, h0, w0;
    int mult_x, mult_y, m_batch_stride;
    int has_grid_indices;
    int _pad0;
};

[shader("compute")]
// Flat 64-thread workgroup (X must be a multiple of the pinned subgroup
// size); logical 8x8 decode below matches the CUDA block geometry.
[numthreads(64, 1, 1)]
void bilagrid_normal_bwd_v1_grid(uint3 bid: SV_GroupID, uint ltid: SV_GroupThreadID,
                                 uniform BnV1GridParams p) {
    uint3 t = uint3(ltid % 8u, ltid / 8u, 0u);
    int block_z_total =
        kPatched != 0 ? p.N * p.m_batch_stride * p.L : p.N * p.L;
    if (int(bid.z) >= block_z_total) return;

    int idx_block = int(bid.z);
    int zi = idx_block % p.L; idx_block /= p.L;
    int m_batch_i = 0;
    if (kPatched != 0) {
        m_batch_i = idx_block % p.m_batch_stride;
        idx_block /= p.m_batch_stride;
    }
    int ni = idx_block;
    int g_id = kPatched != 0
                   ? ni
                   : (p.has_grid_indices != 0 ? p.grid_indices[ni] : ni);
    int corner_base = g_id * p.L * p.H * p.W * kBnChans;

    int base_xi = (int(bid.x) * kBgV1BlockX) / p.mult_x - 1;
    int base_yi = (int(bid.y) * kBgV1BlockY) / p.mult_y - 1;
    int base_zi = zi - 1;
    BgReader r = {p.fp32, p.q16, p.vbounds};

    for (int local = int(ltid); local < kBnCells; local += 64) {
        int lx = local % kBnFootX;
        int ly = (local / kBnFootX) % kBnFootY;
        int lz = local / (kBnFootX * kBnFootY);
        int gxc = min(max(base_xi + lx, 0), p.W - 1);
        int gyc = min(max(base_yi + ly, 0), p.H - 1);
        int gzc = min(max(base_zi + lz, 0), p.L - 1);
        int g_off = corner_base + ((gzc * p.H + gyc) * p.W + gxc) * kBnChans;
        for (int ci = 0; ci < kBnChans; ci++)
            _bn_shmem[_bn_sh_idx(lz, ly, lx, ci)] =
                bg_read(kValueQuant != 0, r, g_off + ci);
    }
    GroupMemoryBarrierWithGroupSync();

    int x_idx = int(bid.x) * kBgV1BlockX + int(t.x);
    int y_idx = int(bid.y) * kBgV1BlockY + int(t.y);
    int xi = x_idx / p.mult_x, xf = x_idx % p.mult_x;
    int yi = y_idx / p.mult_y, yf = y_idx % p.mult_y;
    bool inside = (xi >= 0 && xi < p.W && yi >= 0 && yi < p.H);

    bool slow_path = p.mult_x * p.mult_y == 1 ||
                     (p.mult_x % kBgV1BlockX != 0 ||
                      p.mult_y % kBgV1BlockY != 0);
    if (!inside && slow_path) return;

    float accum[kBnChans] = {0.0f, 0.0f, 0.0f};

    if (inside) {
        int mi_begin = kPatched != 0 ? m_batch_i : 0;
        int mi_end = kPatched != 0 ? p.m : 1;
        int mi_step = kPatched != 0 ? p.m_batch_stride : 1;
        for (int mi = mi_begin; mi < mi_end; mi += mi_step) {
            BgV1Range rng;
            int off_x = 0, off_y = 0;
            if (kPatched != 0) {
                int o = (ni * p.m + mi) * 2;
                off_x = p.offsets[o + 0];
                off_y = p.offsets[o + 1];
                rng = bg_v1_range(xi, yi, xf, yf, p.mult_x, p.mult_y, p.W, p.H,
                                  p.w0, p.h0, off_x, off_y,
                                  min(off_x + p.w, p.w0) - 1,
                                  min(off_y + p.h, p.h0) - 1);
                if (!rng.nonempty && slow_path) continue;
            } else {
                rng = bg_v1_range(xi, yi, xf, yf, p.mult_x, p.mult_y, p.W, p.H,
                                  p.w, p.h, 0, 0, p.w - 1, p.h - 1);
            }

            for (int hi = rng.hi0; hi < rng.hi1; hi++)
            for (int wi = rng.wi0; wi < rng.wi1; wi++) {
                int g_off;
                float x, y;
                if (kPatched != 0) {
                    g_off = (((ni * p.m + mi) * p.h + (hi - off_y)) * p.w +
                             (wi - off_x)) * 3;
                    x = float(wi) / float(p.w0 - 1) * float(p.W - 1);
                    y = float(hi) / float(p.h0 - 1) * float(p.H - 1);
                } else {
                    g_off = ((ni * p.h + hi) * p.w + wi) * 3;
                    x = float(wi) / float(p.w - 1) * float(p.W - 1);
                    y = float(hi) / float(p.h - 1) * float(p.H - 1);
                }
                float sr = p.normal_in[g_off + 0];
                float sg = p.normal_in[g_off + 1];
                float sb = p.normal_in[g_off + 2];
                float inv_norm =
                    1.0f / sqrt(sr * sr + sg * sg + sb * sb + 1e-20f);
                float z = 0.5f + 0.5f * sb * inv_norm;
                z = min(max(z, 0.0f), 1.0f) * float(p.L - 1);

                BgCells c = bg_cells_uniform(0, p.L, p.H, p.W, kBnChans, x, y,
                                             z);
                if (zi != c.z0 && zi != c.z1) continue;

                int sx0 = max(c.x0 - base_xi, 0);
                int sx1 = min(c.x1 - base_xi, kBnFootX - 1);
                int sy0 = max(c.y0 - base_yi, 0);
                int sy1 = min(c.y1 - base_yi, kBnFootY - 1);
                int sz0 = max(c.z0 - base_zi, 0);
                int sz1 = min(c.z1 - base_zi, kBnFootZ - 1);

                float3 axis_angle;
                [ForceUnroll]
                for (int ci = 0; ci < kBnChans; ci++) {
                    float v000 = _bn_shmem[_bn_sh_idx(sz0, sy0, sx0, ci)];
                    float v001 = _bn_shmem[_bn_sh_idx(sz0, sy0, sx1, ci)];
                    float v010 = _bn_shmem[_bn_sh_idx(sz0, sy1, sx0, ci)];
                    float v011 = _bn_shmem[_bn_sh_idx(sz0, sy1, sx1, ci)];
                    float v100 = _bn_shmem[_bn_sh_idx(sz1, sy0, sx0, ci)];
                    float v101 = _bn_shmem[_bn_sh_idx(sz1, sy0, sx1, ci)];
                    float v110 = _bn_shmem[_bn_sh_idx(sz1, sy1, sx0, ci)];
                    float v111 = _bn_shmem[_bn_sh_idx(sz1, sy1, sx1, ci)];
                    float c00 = v000 * (1.0f - c.fx) + v001 * c.fx;
                    float c01 = v010 * (1.0f - c.fx) + v011 * c.fx;
                    float c10 = v100 * (1.0f - c.fx) + v101 * c.fx;
                    float c11 = v110 * (1.0f - c.fx) + v111 * c.fx;
                    float c0 = c00 * (1.0f - c.fy) + c01 * c.fy;
                    float c1 = c10 * (1.0f - c.fy) + c11 * c.fy;
                    float val = c0 * (1.0f - c.fz) + c1 * c.fz;
                    if (ci == 0) axis_angle.x = val;
                    else if (ci == 1) axis_angle.y = val;
                    else axis_angle.z = val;
                }

                float3 grad_out = float3(p.v_output[g_off + 0],
                                         p.v_output[g_off + 1],
                                         p.v_output[g_off + 2]);
                float3 normal =
                    float3(sr * inv_norm, sg * inv_norm, sb * inv_norm);
                float3 grad_axis_angle, grad_normal;
                aa_rotate_bwd(axis_angle, normal, grad_out, grad_axis_angle,
                              grad_normal);

                float accum_t = bg_v1_accum_t(xi, yi, zi, c);
                accum[0] += accum_t * grad_axis_angle.x;
                accum[1] += accum_t * grad_axis_angle.y;
                accum[2] += accum_t * grad_axis_angle.z;
            }
        }
    }

    int out_idx_start = (((ni * p.L + zi) * p.H + yi) * p.W + xi) * kBnChans;

    if (slow_path) {
        [ForceUnroll]
        for (int ci = 0; ci < kBnChans; ci++) {
            if (bg_finite(accum[ci]) && accum[ci] != 0.0f)
                atomic_add_f32(&p.v_bilagrid[out_idx_start + ci], accum[ci]);
        }
        return;
    }

    [ForceUnroll]
    for (int ci = 0; ci < kBnChans; ci++)
        bg_block_add_64(&p.v_bilagrid[out_idx_start + ci], accum[ci], ltid);
}

// ---------------------------------------------------------------------------
// backward v1: normal-grad kernel (flat 256)
// ---------------------------------------------------------------------------

struct BnV1NormalParams {
    float* fp32;
    uint* q16;
    float2* vbounds;
    float* normal_in;
    float* v_output;
    float* v_normal_in;
    int* offsets;
    int* grid_indices;
    int N, L, H, W, m, h, w, h0, w0;
    int has_grid_indices;
    uint total;
    uint wgs_per_row;
};

[shader("compute")]
[numthreads(256, 1, 1)]
void bilagrid_normal_bwd_v1_normal(uint3 wg: SV_GroupID,
                                   uint tid: SV_GroupThreadID,
                                   uniform BnV1NormalParams p) {
    uint idx = (wg.y * p.wgs_per_row + wg.x) * 256 + tid;
    if (idx >= p.total) return;

    BnUniformParams u;
    u.fp32 = p.fp32; u.q16 = p.q16; u.vbounds = p.vbounds;
    u.normal_in = p.normal_in; u.output = p.v_output;
    u.offsets = p.offsets; u.grid_indices = p.grid_indices;
    u.N = p.N; u.L = p.L; u.H = p.H; u.W = p.W;
    u.m = p.m; u.h = p.h; u.w = p.w; u.h0 = p.h0; u.w0 = p.w0;
    u.has_grid_indices = p.has_grid_indices;
    u.total = p.total; u.wgs_per_row = p.wgs_per_row;
    int ni, g_id, g_off3;
    float gx, gy;
    _bn_pix(int(idx), u, ni, g_id, g_off3, gx, gy);
    BgReader r = {p.fp32, p.q16, p.vbounds};

    float sr = p.normal_in[g_off3 + 0];
    float sg = p.normal_in[g_off3 + 1];
    float sb = p.normal_in[g_off3 + 2];
    float inv_norm = 1.0f / sqrt(sr * sr + sg * sg + sb * sb + 1e-20f);
    float3 unit_normal = float3(sr * inv_norm, sg * inv_norm, sb * inv_norm);
    float3 grad_out = float3(p.v_output[g_off3 + 0], p.v_output[g_off3 + 1],
                             p.v_output[g_off3 + 2]);

    float x = gx * float(p.W - 1);
    float y = gy * float(p.H - 1);
    float gz_raw = 0.5f + 0.5f * sb * inv_norm;
    float gz = min(max(gz_raw, 0.0f), 1.0f);
    bool gz_in_range = (gz_raw >= 0.0f && gz_raw <= 1.0f);
    float z = gz * float(p.L - 1);

    int base = g_id * p.L * p.H * p.W * 3;
    BgCells c = bg_cells_uniform(base, p.L, p.H, p.W, 3, x, y, z);

    float3 axis_angle;
    axis_angle.x = bg_trilerp(kValueQuant != 0, r, c, 0);
    axis_angle.y = bg_trilerp(kValueQuant != 0, r, c, 1);
    axis_angle.z = bg_trilerp(kValueQuant != 0, r, c, 2);

    float3 grad_axis_angle, grad_unit_normal;
    aa_rotate_bwd(axis_angle, unit_normal, grad_out, grad_axis_angle,
                  grad_unit_normal);

    float d8[8];
    bg_dwdz(c, d8);
    float gz_grad = 0.0f;
    [ForceUnroll]
    for (int corner = 0; corner < 8; corner++) {
        float trilerp = 0.0f;
        [ForceUnroll]
        for (int ci = 0; ci < 3; ci++) {
            float gv = bg_read(kValueQuant != 0, r, c.o[corner] + ci);
            trilerp += gv * (ci == 0 ? grad_axis_angle.x
                                     : ci == 1 ? grad_axis_angle.y
                                               : grad_axis_angle.z);
        }
        gz_grad += d8[corner] * float(p.L - 1) * trilerp;
    }
    if (!gz_in_range) gz_grad = 0.0f;
    grad_unit_normal.z += 0.5f * gz_grad;

    // vjp of raw -> unit normalization.
    float3 grad_normal =
        (grad_unit_normal +
         unit_normal * (-dot(unit_normal, grad_unit_normal))) * inv_norm;
    p.v_normal_in[g_off3 + 0] = bg_finite(grad_normal.x) ? grad_normal.x : 0.0f;
    p.v_normal_in[g_off3 + 1] = bg_finite(grad_normal.y) ? grad_normal.y : 0.0f;
    p.v_normal_in[g_off3 + 2] = bg_finite(grad_normal.z) ? grad_normal.z : 0.0f;
}

// ---------------------------------------------------------------------------
// backward v2: thread-per-pixel SCATTER (grid-grad only; normal is GT-side, so
// needs_image_grad=false -> no v_normal_in, no gz chain). Read grid values from
// the g_id slot, scatter grad into the per-batch ni buffer. Transform-backward
// is the hand-written Rodrigues form (aa_rotate_bwd). See
// BilagridBackwardSelection.md and BilagridNormalSampleBwdV2_kernel.cuh.
// ---------------------------------------------------------------------------

struct BnV2Params {
    float* fp32;
    uint* q16;
    float2* vbounds;
    float* normal_in;
    float* v_output;
    float* v_bilagrid;
    int* grid_indices;
    int N, L, H, W, h, w;
    int has_grid_indices;
    uint total;
    uint wgs_per_row;
    uint _pad0;
};

[shader("compute")]
[numthreads(256, 1, 1)]
void bilagrid_normal_bwd_v2(uint3 wg: SV_GroupID, uint tid: SV_GroupThreadID,
                            uniform BnV2Params p) {
    uint idx = (wg.y * p.wgs_per_row + wg.x) * 256 + tid;
    if (idx >= p.total) return;
    int tmp = int(idx);
    int wi = tmp % p.w; tmp /= p.w;
    int hi = tmp % p.h; tmp /= p.h;
    int ni = tmp;
    int g_id = p.has_grid_indices != 0 ? p.grid_indices[ni] : ni;
    int g_off = ((ni * p.h + hi) * p.w + wi) * 3;
    BgReader r = {p.fp32, p.q16, p.vbounds};

    float sr = p.normal_in[g_off + 0];
    float sg = p.normal_in[g_off + 1];
    float sb = p.normal_in[g_off + 2];
    float inv_norm = 1.0f / sqrt(sr * sr + sg * sg + sb * sb + 1e-20f);
    float x = float(wi) / float(p.w - 1) * float(p.W - 1);
    float y = float(hi) / float(p.h - 1) * float(p.H - 1);
    float z = 0.5f + 0.5f * sb * inv_norm;
    z = min(max(z, 0.0f), 1.0f) * float(p.L - 1);

    int base = g_id * p.L * p.H * p.W * kBnChans;
    int out_base = ni * p.L * p.H * p.W * kBnChans;
    BgCells c = bg_cells_uniform(base, p.L, p.H, p.W, kBnChans, x, y, z);
    float w8[8];
    bg_corner_weights(c, w8);

    float3 axis_angle;
    axis_angle.x = bg_trilerp(kValueQuant != 0, r, c, 0);
    axis_angle.y = bg_trilerp(kValueQuant != 0, r, c, 1);
    axis_angle.z = bg_trilerp(kValueQuant != 0, r, c, 2);

    float3 grad_out = float3(p.v_output[g_off + 0], p.v_output[g_off + 1],
                             p.v_output[g_off + 2]);
    float3 normal = float3(sr * inv_norm, sg * inv_norm, sb * inv_norm);
    float3 grad_axis_angle, grad_normal;
    aa_rotate_bwd(axis_angle, normal, grad_out, grad_axis_angle, grad_normal);
    float grad_weight[kBnChans];
    grad_weight[0] = grad_axis_angle.x;
    grad_weight[1] = grad_axis_angle.y;
    grad_weight[2] = grad_axis_angle.z;

    [ForceUnroll]
    for (int corner = 0; corner < 8; corner++) {
        int wbase = c.o[corner] - base + out_base;
        [ForceUnroll]
        for (int ci = 0; ci < kBnChans; ci++) {
            float aa = w8[corner] * grad_weight[ci];
            if (aa != 0.0f && bg_finite(aa))
                atomic_add_f32(&p.v_bilagrid[wbase + ci], aa);
        }
    }
}
