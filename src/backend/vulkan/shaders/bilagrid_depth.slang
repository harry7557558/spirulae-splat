// Depth (2-channel log-linear) bilateral-grid samplers (Vulkan backend),
// mirroring BilagridDepthSampleFwd_kernel.cuh and
// BilagridDepthSampleBwdV1_kernel.cuh (uniform + PATCHED).
//
// Spec constants (IDs by declaration order — dispatches must pass ALL):
//   0: kValueQuant  1: kPatched

#include "bilagrid_common.slang"

[SpecializationConstant]
const int kValueQuant = 0;
[SpecializationConstant]
const int kPatched = 0;

struct BdUniformParams {
    float* fp32;
    uint* q16;
    float2* vbounds;
    float* depth;       // [N,(m,)h,w,1]
    float* scalars;     // [N] (uniform: scalars[g_id]; patched: scalars[0])
    float* output;      // fwd out / bwd v_output
    int* offsets;
    int* grid_indices;
    int N, L, H, W, m, h, w, h0, w0;
    int has_grid_indices;
    uint total;
    uint wgs_per_row;
};

[ForceInline]
void _bd_pix(int idx, BdUniformParams p, out int ni, out int g_id,
             out int g_off, out float gx, out float gy) {
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
        g_off = ((ni * p.m + mi) * p.h + hi) * p.w + wi;
    } else {
        g_id = p.has_grid_indices != 0 ? p.grid_indices[ni] : ni;
        gx = float(wi) / float(p.w - 1);
        gy = float(hi) / float(p.h - 1);
        g_off = (ni * p.h + hi) * p.w + wi;
    }
}

[shader("compute")]
[numthreads(256, 1, 1)]
void bilagrid_depth_fwd(uint3 wg: SV_GroupID, uint tid: SV_GroupThreadID,
                        uniform BdUniformParams p) {
    uint idx = (wg.y * p.wgs_per_row + wg.x) * 256 + tid;
    if (idx >= p.total) return;
    int ni, g_id, g_off;
    float gx, gy;
    _bd_pix(int(idx), p, ni, g_id, g_off, gx, gy);
    BgReader r = {p.fp32, p.q16, p.vbounds};

    float scalar = kPatched != 0 ? p.scalars[0] : p.scalars[g_id];
    float sr = p.depth[g_off] * scalar;

    float gz = sr / (sr + 1.0f);
    gz = min(max(gz, 0.0f), 1.0f);
    float x = gx * float(p.W - 1);
    float y = gy * float(p.H - 1);
    float z = gz * float(p.L - 1);

    int base = g_id * p.L * p.H * p.W * 2;
    BgCells c = bg_cells_uniform(base, p.L, p.H, p.W, 2, x, y, z);

    float dr = 0.0f;
    [ForceUnroll]
    for (int ci = 0; ci < 2; ci++) {
        float val = bg_trilerp(kValueQuant != 0, r, c, ci);
        dr += (ci == 0 ? exp(val) * log(sr) : val);
    }
    dr = exp(dr);
    p.output[g_off] = bg_finite(dr) ? dr : 0.5f;
}

// ---------------------------------------------------------------------------
// backward v1: grid-grad kernel (8x8x1, 2-channel shmem preload)
// ---------------------------------------------------------------------------

static const int kBdFootX = kBgV1BlockX + 2;
static const int kBdFootY = kBgV1BlockY + 2;
static const int kBdFootZ = 3;
static const int kBdChans = 2;
static const int kBdCells = kBdFootX * kBdFootY * kBdFootZ;

groupshared float _bd_shmem[kBdFootZ * kBdFootY * kBdFootX * kBdChans];

[ForceInline]
int _bd_sh_idx(int lz, int ly, int lx, int ci) {
    return ((lz * kBdFootY + ly) * kBdFootX + lx) * kBdChans + ci;
}

struct BdV1GridParams {
    float* fp32;
    uint* q16;
    float2* vbounds;
    float* depth;
    float* scalars;
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
void bilagrid_depth_bwd_v1_grid(uint3 bid: SV_GroupID, uint ltid: SV_GroupThreadID,
                                uniform BdV1GridParams p) {
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
    int g_id;
    float scalar;
    if (kPatched != 0) {
        g_id = ni;
        scalar = p.scalars[0];
    } else {
        g_id = p.has_grid_indices != 0 ? p.grid_indices[ni] : ni;
        scalar = p.scalars[g_id];
    }
    int corner_base = g_id * p.L * p.H * p.W * kBdChans;

    int base_xi = (int(bid.x) * kBgV1BlockX) / p.mult_x - 1;
    int base_yi = (int(bid.y) * kBgV1BlockY) / p.mult_y - 1;
    int base_zi = zi - 1;
    BgReader r = {p.fp32, p.q16, p.vbounds};

    for (int local = int(ltid); local < kBdCells; local += 64) {
        int lx = local % kBdFootX;
        int ly = (local / kBdFootX) % kBdFootY;
        int lz = local / (kBdFootX * kBdFootY);
        int gxc = min(max(base_xi + lx, 0), p.W - 1);
        int gyc = min(max(base_yi + ly, 0), p.H - 1);
        int gzc = min(max(base_zi + lz, 0), p.L - 1);
        int g_off = corner_base + ((gzc * p.H + gyc) * p.W + gxc) * kBdChans;
        for (int ci = 0; ci < kBdChans; ci++)
            _bd_shmem[_bd_sh_idx(lz, ly, lx, ci)] =
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

    float accum[kBdChans] = {0.0f, 0.0f};

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
                    g_off = ((ni * p.m + mi) * p.h + (hi - off_y)) * p.w +
                            (wi - off_x);
                    x = float(wi) / float(p.w0 - 1) * float(p.W - 1);
                    y = float(hi) / float(p.h0 - 1) * float(p.H - 1);
                } else {
                    g_off = (ni * p.h + hi) * p.w + wi;
                    x = float(wi) / float(p.w - 1) * float(p.W - 1);
                    y = float(hi) / float(p.h - 1) * float(p.H - 1);
                }
                float sr = p.depth[g_off] * scalar;
                float z = sr / (sr + 1.0f);
                z = min(max(z, 0.0f), 1.0f) * float(p.L - 1);

                BgCells c = bg_cells_uniform(0, p.L, p.H, p.W, kBdChans, x, y,
                                             z);
                if (zi != c.z0 && zi != c.z1) continue;

                float accum_t = bg_v1_accum_t(xi, yi, zi, c);

                int sx0 = max(c.x0 - base_xi, 0);
                int sx1 = min(c.x1 - base_xi, kBdFootX - 1);
                int sy0 = max(c.y0 - base_yi, 0);
                int sy1 = min(c.y1 - base_yi, kBdFootY - 1);
                int sz0 = max(c.z0 - base_zi, 0);
                int sz1 = min(c.z1 - base_zi, kBdFootZ - 1);

                float affine[kBdChans];
                [ForceUnroll]
                for (int ci = 0; ci < kBdChans; ci++) {
                    float v000 = _bd_shmem[_bd_sh_idx(sz0, sy0, sx0, ci)];
                    float v001 = _bd_shmem[_bd_sh_idx(sz0, sy0, sx1, ci)];
                    float v010 = _bd_shmem[_bd_sh_idx(sz0, sy1, sx0, ci)];
                    float v011 = _bd_shmem[_bd_sh_idx(sz0, sy1, sx1, ci)];
                    float v100 = _bd_shmem[_bd_sh_idx(sz1, sy0, sx0, ci)];
                    float v101 = _bd_shmem[_bd_sh_idx(sz1, sy0, sx1, ci)];
                    float v110 = _bd_shmem[_bd_sh_idx(sz1, sy1, sx0, ci)];
                    float v111 = _bd_shmem[_bd_sh_idx(sz1, sy1, sx1, ci)];
                    float c00 = v000 * (1.0f - c.fx) + v001 * c.fx;
                    float c01 = v010 * (1.0f - c.fx) + v011 * c.fx;
                    float c10 = v100 * (1.0f - c.fx) + v101 * c.fx;
                    float c11 = v110 * (1.0f - c.fx) + v111 * c.fx;
                    float c0 = c00 * (1.0f - c.fy) + c01 * c.fy;
                    float c1 = c10 * (1.0f - c.fy) + c11 * c.fy;
                    float val = c0 * (1.0f - c.fz) + c1 * c.fz;
                    affine[ci] = (ci == 0 ? exp(val) : val);
                }
                float log_sr = log(sr);
                float output = exp(affine[0] * log_sr + affine[1]);
                float dr = p.v_output[g_off] * output;

                [ForceUnroll]
                for (int ci = 0; ci < kBdChans; ci++) {
                    float r_coeff = (ci == 0 ? log_sr : 1.0f);
                    float grad_weight = r_coeff * dr;
                    if (ci == 0) grad_weight *= affine[ci];
                    accum[ci] += accum_t * grad_weight;
                }
            }
        }
    }

    int out_idx_start = (((ni * p.L + zi) * p.H + yi) * p.W + xi) * kBdChans;

    if (slow_path) {
        [ForceUnroll]
        for (int ci = 0; ci < kBdChans; ci++) {
            if (bg_finite(accum[ci]) && accum[ci] != 0.0f)
                atomic_add_f32(&p.v_bilagrid[out_idx_start + ci], accum[ci]);
        }
        return;
    }

    [ForceUnroll]
    for (int ci = 0; ci < kBdChans; ci++)
        bg_block_add_64(&p.v_bilagrid[out_idx_start + ci], accum[ci], ltid);
}

// ---------------------------------------------------------------------------
// backward v1: depth-grad kernel (flat 256)
// ---------------------------------------------------------------------------

struct BdV1DepthParams {
    float* fp32;
    uint* q16;
    float2* vbounds;
    float* depth;
    float* scalars;
    float* v_output;
    float* v_depth;
    int* offsets;
    int* grid_indices;
    int N, L, H, W, m, h, w, h0, w0;
    int has_grid_indices;
    uint total;
    uint wgs_per_row;
    uint _pad0;
};

[shader("compute")]
[numthreads(256, 1, 1)]
void bilagrid_depth_bwd_v1_depth(uint3 wg: SV_GroupID,
                                 uint tid: SV_GroupThreadID,
                                 uniform BdV1DepthParams p) {
    uint idx = (wg.y * p.wgs_per_row + wg.x) * 256 + tid;
    if (idx >= p.total) return;

    BdUniformParams u;
    u.fp32 = p.fp32; u.q16 = p.q16; u.vbounds = p.vbounds;
    u.depth = p.depth; u.scalars = p.scalars; u.output = p.v_output;
    u.offsets = p.offsets; u.grid_indices = p.grid_indices;
    u.N = p.N; u.L = p.L; u.H = p.H; u.W = p.W;
    u.m = p.m; u.h = p.h; u.w = p.w; u.h0 = p.h0; u.w0 = p.w0;
    u.has_grid_indices = p.has_grid_indices;
    u.total = p.total; u.wgs_per_row = p.wgs_per_row;
    int ni, g_id, g_off;
    float gx, gy;
    _bd_pix(int(idx), u, ni, g_id, g_off, gx, gy);
    BgReader r = {p.fp32, p.q16, p.vbounds};

    float scalar = kPatched != 0 ? p.scalars[0] : p.scalars[g_id];
    float sr = p.depth[g_off] * scalar;
    float dr = p.v_output[g_off];
    float vr = 0.0f;

    float x = gx * float(p.W - 1);
    float y = gy * float(p.H - 1);
    float gz_raw = sr / (sr + 1.0f);
    float gz = min(max(gz_raw, 0.0f), 1.0f);
    bool gz_in_range = (gz_raw >= 0.0f && gz_raw <= 1.0f);
    float z = gz * float(p.L - 1);

    int base = g_id * p.L * p.H * p.W * 2;
    BgCells c = bg_cells_uniform(base, p.L, p.H, p.W, 2, x, y, z);
    float w8[8];
    bg_corner_weights(c, w8);

    float affine[2];
    [ForceUnroll]
    for (int ci = 0; ci < 2; ci++) {
        float val = 0.0f;
        [ForceUnroll]
        for (int corner = 0; corner < 8; corner++)
            val += bg_read(kValueQuant != 0, r, c.o[corner] + ci) * w8[corner];
        if (ci == 0) val = exp(val);
        affine[ci] = val;
    }
    float log_sr = log(sr);
    dr *= exp(affine[0] * log_sr + affine[1]);
    vr += affine[0] * dr / sr;

    float d8[8];
    bg_dwdz(c, d8);
    float gz_grad = 0.0f;
    [ForceUnroll]
    for (int corner = 0; corner < 8; corner++) {
        float trilerp = 0.0f;
        [ForceUnroll]
        for (int ci = 0; ci < 2; ci++) {
            float gv = bg_read(kValueQuant != 0, r, c.o[corner] + ci);
            float r_coeff = (ci == 0 ? log_sr * affine[ci] : 1.0f);
            trilerp += gv * r_coeff * dr;
        }
        gz_grad += d8[corner] * float(p.L - 1) * trilerp;
    }
    if (!gz_in_range) gz_grad = 0.0f;
    vr += gz_grad / ((sr + 1.0f) * (sr + 1.0f));
    vr *= scalar;
    p.v_depth[g_off] = bg_finite(vr) ? vr : 0.0f;
}

// ---------------------------------------------------------------------------
// backward v2: thread-per-pixel SCATTER (grid-grad only; depth is GT-side, so
// needs_image_grad=false -> no v_depth, no gz chain). Read grid values from the
// g_id slot, scatter grad into the per-batch ni buffer. See
// BilagridBackwardSelection.md and BilagridDepthSampleBwdV2_kernel.cuh.
// ---------------------------------------------------------------------------

struct BdV2Params {
    float* fp32;
    uint* q16;
    float2* vbounds;
    float* depth;
    float* scalars;
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
void bilagrid_depth_bwd_v2(uint3 wg: SV_GroupID, uint tid: SV_GroupThreadID,
                           uniform BdV2Params p) {
    uint idx = (wg.y * p.wgs_per_row + wg.x) * 256 + tid;
    if (idx >= p.total) return;
    int tmp = int(idx);
    int wi = tmp % p.w; tmp /= p.w;
    int hi = tmp % p.h; tmp /= p.h;
    int ni = tmp;
    int g_id = p.has_grid_indices != 0 ? p.grid_indices[ni] : ni;
    float scalar = p.scalars[g_id];
    int g_off = (ni * p.h + hi) * p.w + wi;
    BgReader r = {p.fp32, p.q16, p.vbounds};

    float sr = p.depth[g_off] * scalar;
    float dr = p.v_output[g_off];
    float x = float(wi) / float(p.w - 1) * float(p.W - 1);
    float y = float(hi) / float(p.h - 1) * float(p.H - 1);
    float z = sr / (sr + 1.0f);
    z = min(max(z, 0.0f), 1.0f) * float(p.L - 1);

    int base = g_id * p.L * p.H * p.W * kBdChans;
    int out_base = ni * p.L * p.H * p.W * kBdChans;
    BgCells c = bg_cells_uniform(base, p.L, p.H, p.W, kBdChans, x, y, z);
    float w8[8];
    bg_corner_weights(c, w8);

    float affine[kBdChans];
    [ForceUnroll]
    for (int ci = 0; ci < kBdChans; ci++) {
        float v = bg_trilerp(kValueQuant != 0, r, c, ci);
        affine[ci] = (ci == 0) ? exp(v) : v;
    }
    float log_sr = log(sr);
    dr *= exp(affine[0] * log_sr + affine[1]);
    float grad_weight[kBdChans];
    grad_weight[0] = log_sr * dr * affine[0];
    grad_weight[1] = dr;

    [ForceUnroll]
    for (int corner = 0; corner < 8; corner++) {
        int wbase = c.o[corner] - base + out_base;
        [ForceUnroll]
        for (int ci = 0; ci < kBdChans; ci++) {
            float aa = w8[corner] * grad_weight[ci];
            if (aa != 0.0f && bg_finite(aa))
                atomic_add_f32(&p.v_bilagrid[wbase + ci], aa);
        }
    }
}
