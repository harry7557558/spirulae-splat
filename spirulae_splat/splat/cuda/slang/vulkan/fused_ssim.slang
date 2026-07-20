// Fused SSIM forward+backward (Vulkan backend), mirroring the
// memory_efficient_ssim_backward_kernel<inplace=true> path of
// csrc/FusedSSIM.cu -- the only variant the multi-scale loss host uses.
//
// CUDA runs 24x24 blocks; the resulting shared arrays (~68 KB) exceed the
// 32 KB llvmpipe groupshared floor, so this port runs 16x16 logical tiles as
// flat 256-wide workgroups (VUID-02757: workgroup X must be a multiple of
// the pinned 32-wide subgroup). Every per-pixel output (SSIM map, dL/dimg1)
// is a pure function of the zero-padded input image, so the tile size does
// not change results -- only the atomic ssim scalar's accumulation grouping
// differs (loose parity channel). Per-thread inner-loop summation order (5
// symmetric taps + center) matches CUDA exactly.
//
// Caveat mirrored from CUDA ("TODO: mask this"): a masked pixel OVERWRITES
// its thread's ssim_val accumulator with 3.0, so with a mask present the
// scalar depends on the thread->accumulated-entries mapping and is NOT
// cross-backend comparable. dL/dimg1 and the loss map stay deterministic.

#include "int8_compat.slang"
#include "atomic_float.slang"

static const int BX = 16;
static const int BY = 16;
static const int HALO = 5;
static const int SX = BX + 2 * HALO;    // 26
static const int SY = BY + 2 * HALO;    // 26
static const int CX = BX;               // 16
static const int CY = SY;               // 26
static const int LX = BX + 4 * HALO;    // 36
static const int LY = BY + 4 * HALO;    // 36

static const float kGauss[HALO + 1] = {
    0.001028380123898387f, 0.0075987582094967365f, 0.036000773310661316f,
    0.10936068743467331f,  0.21300552785396576f,   0.26601171493530273f,
};

static const float kC1 = 0.01f * 0.01f;
static const float kC2 = 0.03f * 0.03f;

// flags bits
static const uint kSsimHasMask = 1u << 0;
static const uint kSsimHasVal = 1u << 1;
static const uint kSsimWriteGrad = 1u << 2;
static const uint kSsimHasLossMap = 1u << 3;

struct SsimParams {
    float* img1;      // [B,H,W,3]
    float* img2;      // [B,H,W,3]
    uint* masks;      // [Bm,Hm,Wm] bool bytes
    float* dL_dimg1;  // [B,H,W,3], accumulated (+=)
    float* ssim_val;  // scalar accumulator
    float* loss_map;  // [B,H,W,1], accumulated (+=)
    int B, H, W;
    int Bm, Hm, Wm;
    float dL_dmap;
    float map_weight;
    int mode;         // SsimLossMapMode: 0 skip, 1 LCS, 2 CS, 3 structure
    uint flags;
};

[ForceInline]
float _fs_pix(float* img, int b, int y, int x, int ci, int H, int W) {
    if (uint(x) >= uint(W) || uint(y) >= uint(H)) return 0.0f;
    return img[(((b * H) + y) * W + x) * 3 + ci];
}

[ForceInline]
bool _fs_mask(SsimParams p, int b, int y, int x) {
    if ((p.flags & kSsimHasMask) == 0) return true;
    if (uint(b) >= uint(p.Bm) || uint(x) >= uint(p.Wm) ||
        uint(y) >= uint(p.Hm))
        return true;
    uint idx = (uint(b) * uint(p.Hm) + uint(y)) * uint(p.Wm) + uint(x);
    return u8_load(p.masks, idx) != 0u;
}

// Shared buffers reused across stages exactly like the CUDA kernel's
// _shared_buffer1/2 overlays.
groupshared float _fs_sb1[LY * LX * 2];             // sTile | sData (2592/2028)
groupshared float _fs_sb2[LY * SX * 5];             // xconv | sScratch (4680/1248)
groupshared float _fs_ssim[BY * BX];                // per-pixel map, then reduce

#define sTile(i, j, k) _fs_sb1[((i) * LX + (j)) * 2 + (k)]
#define xconv(i, j, k) _fs_sb2[((i) * SX + (j)) * 5 + (k)]
#define sData(k, i, j) _fs_sb1[((k) * SY + (i)) * SX + (j)]
#define sScratch(i, j, k) _fs_sb2[((i) * CX + (j)) * 3 + (k)]

[shader("compute")]
[numthreads(256, 1, 1)]
void ssim_bwd_inplace(uint3 wg: SV_GroupID, uint tid: SV_GroupThreadID,
                      uniform SsimParams p) {
    int tx = int(tid % 16u);
    int ty = int(tid / 16u);
    const int bIdx = int(wg.z);
    const int startY = int(wg.y) * BY;
    const int startX = int(wg.x) * BX;
    const int pix_y = startY + ty;
    const int pix_x = startX + tx;
    const int pix_id = pix_y * p.W + pix_x;
    const int num_pix = p.H * p.W;

    const bool lm_active = p.mode != 0 && (p.flags & kSsimHasLossMap) != 0;
    const bool has_val = (p.flags & kSsimHasVal) != 0;

    const float grad =
        p.dL_dmap / (3.0f * float(p.B) * float(p.W) * float(p.H));

    float ssim_val = 0.0f;

    for (int ci = 0; ci < 3; ci++) {
        if (ci != 0) GroupMemoryBarrierWithGroupSync();

        // (1) Load img1/img2 tile + double halo.
        for (int s = 0; s < (LX * LY + 255) / 256; s++) {
            int t = s * 256 + int(tid);
            if (t < LX * LY) {
                int ly = t / LX;
                int lx = t % LX;
                int gy = startY + ly - 2 * HALO;
                int gx = startX + lx - 2 * HALO;
                float X = _fs_pix(p.img1, bIdx, gy, gx, ci, p.H, p.W);
                float Y = _fs_pix(p.img2, bIdx, gy, gx, ci, p.H, p.W);
                if (!_fs_mask(p, bIdx, gy, gx)) X = Y;
                sTile(ly, lx, 0) = X;
                sTile(ly, lx, 1) = Y;
            }
        }
        GroupMemoryBarrierWithGroupSync();

        // (2) Horizontal 11-tap convolution of the 5 moment sums.
        for (int l = int(tid); l < SX * LY; l += 256) {
            int lx = l % SX + HALO;
            int ly = l / SX;
            float sumX = 0.0f, sumX2 = 0.0f, sumY = 0.0f, sumY2 = 0.0f,
                  sumXY = 0.0f;
            [ForceUnroll]
            for (int d = 1; d <= HALO; d++) {
                float w = kGauss[HALO - d];
                float Xl = sTile(ly, lx - d, 0);
                float Yl = sTile(ly, lx - d, 1);
                float Xr = sTile(ly, lx + d, 0);
                float Yr = sTile(ly, lx + d, 1);
                sumX += (Xl + Xr) * w;
                sumX2 += ((Xl * Xl) + (Xr * Xr)) * w;
                sumY += (Yl + Yr) * w;
                sumY2 += ((Yl * Yl) + (Yr * Yr)) * w;
                sumXY += ((Xl * Yl) + (Xr * Yr)) * w;
            }
            {
                float cX = sTile(ly, lx, 0);
                float cY = sTile(ly, lx, 1);
                float wc = kGauss[HALO];
                sumX += cX * wc;
                sumX2 += (cX * cX) * wc;
                sumY += cY * wc;
                sumY2 += (cY * cY) * wc;
                sumXY += (cX * cY) * wc;
            }
            xconv(ly, lx - HALO, 0) = sumX;
            xconv(ly, lx - HALO, 1) = sumX2;
            xconv(ly, lx - HALO, 2) = sumY;
            xconv(ly, lx - HALO, 3) = sumY2;
            xconv(ly, lx - HALO, 4) = sumXY;
        }
        GroupMemoryBarrierWithGroupSync();

        // (3) Vertical pass -> SSIM value + partial derivatives.
        for (int l = int(tid); l < SX * SY; l += 256) {
            int lx = l % SX;
            int ly = l / SX + HALO;
            float out0 = 0.0f, out1 = 0.0f, out2 = 0.0f, out3 = 0.0f,
                  out4 = 0.0f;
            [ForceUnroll]
            for (int d = 1; d <= HALO; d++) {
                float w = kGauss[HALO - d];
                out0 += (xconv(ly - d, lx, 0) + xconv(ly + d, lx, 0)) * w;
                out1 += (xconv(ly - d, lx, 1) + xconv(ly + d, lx, 1)) * w;
                out2 += (xconv(ly - d, lx, 2) + xconv(ly + d, lx, 2)) * w;
                out3 += (xconv(ly - d, lx, 3) + xconv(ly + d, lx, 3)) * w;
                out4 += (xconv(ly - d, lx, 4) + xconv(ly + d, lx, 4)) * w;
            }
            {
                float wC = kGauss[HALO];
                out0 += xconv(ly, lx, 0) * wC;
                out1 += xconv(ly, lx, 1) * wC;
                out2 += xconv(ly, lx, 2) * wC;
                out3 += xconv(ly, lx, 3) * wC;
                out4 += xconv(ly, lx, 4) * wC;
            }

            float mu1 = out0;
            float mu2 = out2;
            float mu1_sq = mu1 * mu1;
            float mu2_sq = mu2 * mu2;
            float sigma1_sq = out1 - mu1_sq;
            float sigma2_sq = out3 - mu2_sq;
            float sigma12 = out4 - mu1 * mu2;

            float A = mu1_sq + mu2_sq + kC1;
            float B_ = sigma1_sq + sigma2_sq + kC2;
            float C_ = 2.0f * mu1 * mu2 + kC1;
            float D_ = 2.0f * sigma12 + kC2;

            if ((has_val || lm_active) && lx < BX && (ly - HALO) < BY) {
                float ssim_v = (C_ * D_) / (A * B_);
                if (lm_active) {
                    float lm_v;
                    if (p.mode == 2) {  // contrast*structure
                        lm_v = D_ / B_;
                    } else if (p.mode == 3) {  // structure only
                        float s1s2 = sqrt(max(sigma1_sq, 0.0f) *
                                          max(sigma2_sq, 0.0f));
                        lm_v = D_ / (2.0f * s1s2 + kC2);
                    } else {  // full SSIM
                        lm_v = ssim_v;
                    }
                    if (ci == 0)
                        _fs_ssim[(ly - HALO) * BX + lx] = lm_v;
                    else
                        _fs_ssim[(ly - HALO) * BX + lx] += lm_v;
                }
                ssim_val += ssim_v;
            }

            float d_m_dmu1 = ((mu2 * 2.0f * D_) / (A * B_) -
                              (mu2 * 2.0f * C_) / (A * B_) -
                              (mu1 * 2.0f * C_ * D_) / (A * A * B_) +
                              (mu1 * 2.0f * C_ * D_) / (A * B_ * B_));
            float d_m_dsigma1_sq = (-C_ * D_) / (A * B_ * B_);
            float d_m_dsigma12 = (2.0f * C_) / (A * B_);

            int dpix_y = startY + ly - HALO;
            int dpix_x = startX + lx;
            float masked_grad =
                grad * float(dpix_x >= HALO && dpix_y >= HALO &&
                             dpix_x < p.W + HALO && dpix_y < p.H + HALO);

            sData(0, ly - HALO, lx) = d_m_dmu1 * masked_grad;
            sData(1, ly - HALO, lx) = d_m_dsigma1_sq * masked_grad;
            sData(2, ly - HALO, lx) = d_m_dsigma12 * masked_grad;
        }
        GroupMemoryBarrierWithGroupSync();

        // (4) Horizontal pass over the derivative maps.
        {
            int ly = ty;
            int lx = tx + HALO;
            for (int pass = 0; pass < 2; pass++) {
                int yy = ly + pass * BY;
                if (yy < CY) {
                    float a0 = 0.0f, a1 = 0.0f, a2 = 0.0f;
                    [ForceUnroll]
                    for (int d = 1; d <= HALO; d++) {
                        float w = kGauss[HALO - d];
                        a0 += (sData(0, yy, lx - d) + sData(0, yy, lx + d)) * w;
                        a1 += (sData(1, yy, lx - d) + sData(1, yy, lx + d)) * w;
                        a2 += (sData(2, yy, lx - d) + sData(2, yy, lx + d)) * w;
                    }
                    {
                        float wc = kGauss[HALO];
                        a0 += sData(0, yy, lx) * wc;
                        a1 += sData(1, yy, lx) * wc;
                        a2 += sData(2, yy, lx) * wc;
                    }
                    sScratch(yy, tx, 0) = a0;
                    sScratch(yy, tx, 1) = a1;
                    sScratch(yy, tx, 2) = a2;
                }
            }
        }
        GroupMemoryBarrierWithGroupSync();

        // (5) Vertical pass -> dL/d(img1).
        if (pix_x < p.W && pix_y < p.H) {
            int ly = ty + HALO;
            int lx = tx;
            float sum0 = 0.0f, sum1 = 0.0f, sum2 = 0.0f;
            [ForceUnroll]
            for (int d = 1; d <= HALO; d++) {
                float w = kGauss[HALO - d];
                sum0 += (sScratch(ly - d, lx, 0) + sScratch(ly + d, lx, 0)) * w;
                sum1 += (sScratch(ly - d, lx, 1) + sScratch(ly + d, lx, 1)) * w;
                sum2 += (sScratch(ly - d, lx, 2) + sScratch(ly + d, lx, 2)) * w;
            }
            {
                float wc = kGauss[HALO];
                sum0 += sScratch(ly, lx, 0) * wc;
                sum1 += sScratch(ly, lx, 1) * wc;
                sum2 += sScratch(ly, lx, 2) * wc;
            }

            float p1 = _fs_pix(p.img1, bIdx, pix_y, pix_x, ci, p.H, p.W);
            float p2 = _fs_pix(p.img2, bIdx, pix_y, pix_x, ci, p.H, p.W);
            bool mask = _fs_mask(p, bIdx, pix_y, pix_x);
            float dL_dpix = sum0 + (2.0f * p1) * sum1 + p2 * sum2;
            if (!mask) {
                dL_dpix = 0.0f;
                ssim_val = 3.0f;  // mirrors CUDA "TODO: mask this"
            }

            if ((p.flags & kSsimWriteGrad) != 0) {
                int out_idx = bIdx * num_pix + pix_id;
                p.dL_dimg1[out_idx * 3 + ci] += dL_dpix;
            }
        }
    }  // ci

    if (lm_active && pix_x < p.W && pix_y < p.H) {
        float ssim_v = _fs_ssim[ty * BX + tx] / 3.0f;
        int out_idx = bIdx * num_pix + pix_id;
        p.loss_map[out_idx] += p.map_weight * (1.0f - ssim_v);
    }

    if (has_val) {
        // Reuse _fs_ssim as the block-reduce scratch (map reads done above).
        GroupMemoryBarrierWithGroupSync();
        ssim_val /= 3.0f * float(p.B) * float(p.H) * float(p.W);
        _fs_ssim[tid] = ssim_val;
        GroupMemoryBarrierWithGroupSync();
        for (uint s = 128; s > 0; s >>= 1) {
            if (tid < s) _fs_ssim[tid] += _fs_ssim[tid + s];
            GroupMemoryBarrierWithGroupSync();
        }
        if (tid == 0 && _fs_ssim[0] != 0.0f)
            atomic_add_f32(p.ssim_val, _fs_ssim[0]);
    }
}

#undef sTile
#undef xconv
#undef sData
#undef sScratch
