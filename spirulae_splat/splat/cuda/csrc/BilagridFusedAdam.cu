// Fused bilagrid optimizer: image-loss gradient (sparse, per-batch) + TV loss
// gradient (computed inline from grid neighborhood) + Adam update (float or
// uint8 quantized). Replaces the previous chain of
//     dense bilagrid_*_grads (atomicAdd from image-loss backward)
//     -> tv_loss_backward (writes more into the same dense grad)
//     -> fused_adam_step{,_8bit} (reads dense grad, updates grids+moments)
//
// The image gradient is now stored sparse over the camera axis as
// [C_batch, C, L, H, W] -- usually 12 KB per batch image instead of the
// full [N_grids, C, L, H, W] table (often >100 MB on photogrammetry-scale
// scenes). The TV gradient is local per cell so we can recompute it on the
// fly during the optimizer pass; nothing needs to be materialized for it.
//
// One template instantiation per channel count and quantize-mode covers all
// bilagrid types (RGB affine=12, RGB ppisp/loglinear=9, normal=3, depth=2).
//
// Trade-off: an Adam update is now ONLY applied to cells whose camera appears
// in the current batch. Cells of untouched cameras get neither image-loss nor
// TV gradient that step. In practice training cycles through every camera, so
// untouched cameras catch up across iterations -- this is the "lossless" mode
// the user signed off on.

#include "BilagridConfig.cuh"
#include "Tensor.h"
#include <cuda_runtime.h>
#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>

namespace cg = cooperative_groups;

#ifndef WARP_SIZE
#define WARP_SIZE 32
#endif

template<int BLOCK_SIZE, int C, bool QUANT, int QUANT_BITS = 8>
__global__ void fused_bilagrid_tv_adam_kernel(
    // Parameter table -- kept on device permanently. Read AND written here.
    // Channel-last layout: the C channels of a single (l,h,w) cell are
    // contiguous, then (cam, l, h, w) sweep outer.
    float* __restrict__ grids,                  // [N_grids, L, H, W, C]
    // Adam moments. float path or QuantizedAdamState path (mutually exclusive).
    float*   __restrict__ g1_f,                 // when !QUANT
    float*   __restrict__ g2_f,                 // when !QUANT
    uint8_t* __restrict__ packed,               // when QUANT (AoS packed cells)
    float4*  __restrict__ quant_bounds,         // [n_blocks], when QUANT
    // Sparse image gradient. ni is the batch slot; the *camera* the gradient
    // belongs to is cam_indices[ni].
    const float* __restrict__ image_grad,       // [C_batch, L, H, W, C]
    const int*   __restrict__ cam_indices,      // [C_batch]
    int N_grids, int C_batch,
    int L, int H, int W,
    float lr,
    float tv_weight,
    int32_t adam_step
) {
    static constexpr float eps   = 1e-15f;
    static constexpr float beta1 = 0.9f;
    static constexpr float beta2 = 0.999f;

    const int64_t total_cells =
        (int64_t)N_grids * (int64_t)C * (int64_t)L * (int64_t)H * (int64_t)W;
    const int64_t idx =
        (int64_t)blockIdx.x * BLOCK_SIZE + (int64_t)threadIdx.x;
    const bool inside = idx < total_cells;

    // Decode (cam, l, h, w, ci) from the flat index. Layout: C innermost.
    int cam = 0, ci = 0, l = 0, h = 0, w = 0;
    if (inside) {
        int64_t t = idx;
        ci  = (int)(t % C); t /= C;
        w   = (int)(t % W); t /= W;
        h   = (int)(t % H); t /= H;
        l   = (int)(t % L); t /= L;
        cam = (int)t;
    }

    // ---- Image gradient ---------------------------------------------------
    // Linear scan over batch slots looking for any whose camera matches this
    // cell's `cam`. For the typical 1- to ~10-image batch this is trivial;
    // also covers duplicate cam_indices (multiple batch images on same cam
    // sum, just like the previous dense atomicAdd path would have done).
    //
    // When cam_indices is null we're in identity mode (batch slot b == cam b),
    // so the only batch slot contributing to cell `cam` is b == cam (if it's
    // in range).
    float image_g = 0.0f;
    if (inside) {
        const int64_t per_slot = (int64_t)L * H * W * C;
        const int64_t cell_off =
            ((((int64_t)l * H + h) * W + w) * C) + ci;
        if (cam_indices == nullptr) {
            if (cam < C_batch)
                image_g = image_grad[(int64_t)cam * per_slot + cell_off];
        } else {
            for (int b = 0; b < C_batch; ++b) {
                if (cam_indices[b] == cam)
                    image_g += image_grad[(int64_t)b * per_slot + cell_off];
            }
        }
    }

    // ---- TV gradient (inline, same formula as tv_loss_backward_kernel) ----
    // Spatial neighbors at the same channel are stride C apart in channel-last
    // (W neighbor), W*C (H neighbor), H*W*C (L neighbor).
    float tv_g = 0.0f;
    if (inside && tv_weight > 0.0f) {
        // Constants from tv_loss_backward.
        const float s  = tv_weight / (6.0f * (float)N_grids);
        const float sx = s / (float)(L * H * (W - 1));
        const float sy = s / (float)(L * (H - 1) * W);
        const float sz = s / (float)((L - 1) * H * W);

        const int64_t cam_base = (int64_t)cam * L * H * W * C;
        const int64_t self_off = (((int64_t)l * H + h) * W + w) * C + ci;
        const float val = grids[cam_base + self_off];

        const int64_t sx_step = (int64_t)C;
        const int64_t sy_step = (int64_t)W * C;
        const int64_t sz_step = (int64_t)H * W * C;

        if (w > 0)     tv_g += (val - grids[cam_base + self_off - sx_step]) * sx;
        if (w < W - 1) tv_g += (val - grids[cam_base + self_off + sx_step]) * sx;
        if (h > 0)     tv_g += (val - grids[cam_base + self_off - sy_step]) * sy;
        if (h < H - 1) tv_g += (val - grids[cam_base + self_off + sy_step]) * sy;
        if (l > 0)     tv_g += (val - grids[cam_base + self_off - sz_step]) * sz;
        if (l < L - 1) tv_g += (val - grids[cam_base + self_off + sz_step]) * sz;
    }

    float v = image_g + tv_g;
    if (!isfinite(v)) v = 0.0f;

    // ---- Adam update ------------------------------------------------------
    const float step_f   = (float)adam_step;
    const float inv_bc1  = 1.0f / (1.0f - powf(beta1, step_f));
    const float inv_bc2  = 1.0f / (1.0f - powf(beta2, step_f));

    float x = inside ? grids[idx] : 0.0f;
    float g1, g2;

    using QState = QuantizedAdamState<QUANT_BITS, BLOCK_SIZE>;
    float4 mm;
    if constexpr (QUANT) {
        // Joint (u, sqrt(g2)) AoS quantization via QuantizedAdamState codec.
        mm = quant_bounds[blockIdx.x];
        if (inside) {
            float2 g1g2 = QState::decode_g1g2(packed, idx, mm);
            g1 = g1g2.x;
            g2 = g1g2.y;
        } else {
            g1 = 0.0f;
            g2 = 0.0f;
        }
    } else {
        g1 = inside ? g1_f[idx] : 0.0f;
        g2 = inside ? g2_f[idx] : 0.0f;
    }

    g1 = beta1 * g1 + (1.0f - beta1) * v;
    g2 = beta2 * g2 + (1.0f - beta2) * v * v;

    x -= lr * inv_bc1 * g1 / (sqrtf(g2 * inv_bc2) + eps);
    if (inside) grids[idx] = x;

    // ---- Writeback (Adam moments) ----------------------------------------
    if constexpr (QUANT) {
        // Re-encode in the (u, sqrt(g2)) basis via QuantizedAdamState codec.
        float2 us = QState::g1g2_to_us(g1, g2);
        float u_new       = us.x;
        float sqrt_g2_new = us.y;

        // Block-reduce {u, sqrt(g2)} min/max for the new per-block bounds.
        cg::thread_block block = cg::this_thread_block();
        cg::thread_block_tile<WARP_SIZE> warp = cg::tiled_partition<WARP_SIZE>(block);
        mm = inside
            ? float4{u_new, u_new, sqrt_g2_new, sqrt_g2_new}
            : float4{1e30f, -1e30f, 1e30f, -1e30f};
        mm.x = cg::reduce(warp, mm.x, cg::less<float>());
        mm.y = cg::reduce(warp, mm.y, cg::greater<float>());
        mm.z = cg::reduce(warp, mm.z, cg::less<float>());
        mm.w = cg::reduce(warp, mm.w, cg::greater<float>());
        __shared__ float4 shared_reduce[BLOCK_SIZE / WARP_SIZE];
        if (threadIdx.x % WARP_SIZE == 0)
            shared_reduce[threadIdx.x / WARP_SIZE] = mm;
        __syncthreads();
        mm = (threadIdx.x < BLOCK_SIZE / WARP_SIZE)
            ? shared_reduce[threadIdx.x]
            : float4{1e30f, -1e30f, 1e30f, -1e30f};
        mm.x = cg::reduce(warp, mm.x, cg::less<float>());
        mm.y = cg::reduce(warp, mm.y, cg::greater<float>());
        mm.z = cg::reduce(warp, mm.z, cg::less<float>());
        mm.w = cg::reduce(warp, mm.w, cg::greater<float>());
        __syncthreads();
        if (threadIdx.x < BLOCK_SIZE / WARP_SIZE)
            shared_reduce[threadIdx.x] = mm;
        __syncthreads();
        mm = shared_reduce[threadIdx.x / WARP_SIZE];

        if (inside) {
            QState::encode_us(packed, idx, u_new, sqrt_g2_new, mm);
        }
        if (threadIdx.x == 0) quant_bounds[blockIdx.x] = mm;
    } else {
        if (inside) {
            g1_f[idx] = g1;
            g2_f[idx] = g2;
        }
    }
}


// Host-side dispatcher: picks the right template instantiation by channel
// count and quantize flag, computes the launch grid, and calls the kernel.
void fused_bilagrid_tv_adam(
    float* grids,
    float*   g1_f,    float*   g2_f,           // null when quantize
    uint8_t* packed,                           // null when !quantize (AoS packed)
    float4*  quant_bounds,                     // null when !quantize
    const float* image_grad,
    const int*   cam_indices,
    int N_grids, int C_batch, int C,
    int L, int H, int W,
    float lr,
    float tv_weight,
    int32_t adam_step,
    bool quantize,
    cudaStream_t stream
) {
    constexpr int BLOCK_SIZE = 256;
    const int64_t total_cells =
        (int64_t)N_grids * (int64_t)C * (int64_t)L * (int64_t)H * (int64_t)W;
    int blocks = (int)((total_cells + BLOCK_SIZE - 1) / BLOCK_SIZE);
    if (blocks == 0) return;

    #define LAUNCH(CC, QQ) \
        fused_bilagrid_tv_adam_kernel<BLOCK_SIZE, CC, QQ, 8> \
            <<<blocks, BLOCK_SIZE, 0, stream>>>( \
                grids, g1_f, g2_f, packed, quant_bounds, \
                image_grad, cam_indices, \
                N_grids, C_batch, L, H, W, \
                lr, tv_weight, adam_step)

    if (quantize) {
        if      (C == 12) { LAUNCH(12, true); }
        else if (C == 9)  { LAUNCH(9,  true); }
        else if (C == 3)  { LAUNCH(3,  true); }
        else if (C == 2)  { LAUNCH(2,  true); }
    } else {
        if      (C == 12) { LAUNCH(12, false); }
        else if (C == 9)  { LAUNCH(9,  false); }
        else if (C == 3)  { LAUNCH(3,  false); }
        else if (C == 2)  { LAUNCH(2,  false); }
    }
    #undef LAUNCH
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


// ============================================================================
// AdaGrad variant -- same image_grad + TV-from-neighbors pipeline, but with
// AdaGrad update:
//     accum += g*g
//     x    -= lr * g / (sqrt(accum) + eps)
// Hyperparameter mapping requested:
//   lr_decay = 0           (no per-step lr decay; caller does scheduling)
//   weight_decay = 0       (no L2)
//   initial_accumulator_value = 0  (matches all-zero state on first step)
//   eps = 1e-15            (matches QuantizedTensorLog::kEps)
//
// Quantization path stores `accum` (NOT sqrt(accum)) through
// QuantizedTensorLog<QUANT_BITS, BLOCK_SIZE> -- block-wise log encoding gives
// uniform relative precision across the orders of magnitude AdaGrad's
// accumulator spans during training. bounds is float2 per block (log-domain
// min/max).
// ============================================================================

template<int BLOCK_SIZE, int C, bool QUANT, int QUANT_BITS = 8>
__global__ void fused_bilagrid_tv_adagrad_kernel(
    float* __restrict__ grids,                  // [N_grids, L, H, W, C]
    float*   __restrict__ accum_f,              // when !QUANT
    uint8_t* __restrict__ packed,               // when QUANT
    float2*  __restrict__ quant_bounds,         // [n_blocks], when QUANT
    const float* __restrict__ image_grad,       // [C_batch, L, H, W, C]
    const int*   __restrict__ cam_indices,      // [C_batch]
    int N_grids, int C_batch,
    int L, int H, int W,
    float lr,
    float tv_weight
) {
    static constexpr float eps = 1e-15f;

    const int64_t total_cells =
        (int64_t)N_grids * (int64_t)C * (int64_t)L * (int64_t)H * (int64_t)W;
    const int64_t idx =
        (int64_t)blockIdx.x * BLOCK_SIZE + (int64_t)threadIdx.x;
    const bool inside = idx < total_cells;

    // Decode (cam, l, h, w, ci) -- same layout as Adam kernel.
    int cam = 0, ci = 0, l = 0, h = 0, w = 0;
    if (inside) {
        int64_t t = idx;
        ci  = (int)(t % C); t /= C;
        w   = (int)(t % W); t /= W;
        h   = (int)(t % H); t /= H;
        l   = (int)(t % L); t /= L;
        cam = (int)t;
    }

    // ---- Image gradient (sparse over batch slots) -----------------------
    float image_g = 0.0f;
    if (inside) {
        const int64_t per_slot = (int64_t)L * H * W * C;
        const int64_t cell_off =
            ((((int64_t)l * H + h) * W + w) * C) + ci;
        if (cam_indices == nullptr) {
            if (cam < C_batch)
                image_g = image_grad[(int64_t)cam * per_slot + cell_off];
        } else {
            for (int b = 0; b < C_batch; ++b) {
                if (cam_indices[b] == cam)
                    image_g += image_grad[(int64_t)b * per_slot + cell_off];
            }
        }
    }

    // ---- TV gradient (inline, same formula as Adam kernel) --------------
    float tv_g = 0.0f;
    if (inside && tv_weight > 0.0f) {
        const float s  = tv_weight / (6.0f * (float)N_grids);
        const float sx = s / (float)(L * H * (W - 1));
        const float sy = s / (float)(L * (H - 1) * W);
        const float sz = s / (float)((L - 1) * H * W);

        const int64_t cam_base = (int64_t)cam * L * H * W * C;
        const int64_t self_off = (((int64_t)l * H + h) * W + w) * C + ci;
        const float val = grids[cam_base + self_off];

        const int64_t sx_step = (int64_t)C;
        const int64_t sy_step = (int64_t)W * C;
        const int64_t sz_step = (int64_t)H * W * C;

        if (w > 0)     tv_g += (val - grids[cam_base + self_off - sx_step]) * sx;
        if (w < W - 1) tv_g += (val - grids[cam_base + self_off + sx_step]) * sx;
        if (h > 0)     tv_g += (val - grids[cam_base + self_off - sy_step]) * sy;
        if (h < H - 1) tv_g += (val - grids[cam_base + self_off + sy_step]) * sy;
        if (l > 0)     tv_g += (val - grids[cam_base + self_off - sz_step]) * sz;
        if (l < L - 1) tv_g += (val - grids[cam_base + self_off + sz_step]) * sz;
    }

    float v = image_g + tv_g;
    if (!isfinite(v)) v = 0.0f;

    // ---- AdaGrad update --------------------------------------------------
    float x = inside ? grids[idx] : 0.0f;
    float accum;

    using QTL = QuantizedTensorLog<QUANT_BITS, BLOCK_SIZE>;
    float2 mm;
    if constexpr (QUANT) {
        mm = quant_bounds[blockIdx.x];
        accum = inside ? QTL::decode_v(packed, idx, mm) : 0.0f;
    } else {
        accum = inside ? accum_f[idx] : 0.0f;
    }

    accum += v * v;
    x -= lr * v / (sqrtf(accum) + eps);
    if (inside) grids[idx] = x;

    // ---- Writeback -------------------------------------------------------
    if constexpr (QUANT) {
        float log_v = inside ? QTL::forward_v(accum) : 0.0f;

        // Block-reduce log-domain min/max for the new per-block bounds.
        cg::thread_block block = cg::this_thread_block();
        cg::thread_block_tile<WARP_SIZE> warp = cg::tiled_partition<WARP_SIZE>(block);
        float2 lm = inside ? float2{log_v, log_v}
                           : float2{1e30f, -1e30f};
        lm.x = cg::reduce(warp, lm.x, cg::less<float>());
        lm.y = cg::reduce(warp, lm.y, cg::greater<float>());
        __shared__ float2 shared_reduce[BLOCK_SIZE / WARP_SIZE];
        if (threadIdx.x % WARP_SIZE == 0)
            shared_reduce[threadIdx.x / WARP_SIZE] = lm;
        __syncthreads();
        lm = (threadIdx.x < BLOCK_SIZE / WARP_SIZE)
            ? shared_reduce[threadIdx.x]
            : float2{1e30f, -1e30f};
        lm.x = cg::reduce(warp, lm.x, cg::less<float>());
        lm.y = cg::reduce(warp, lm.y, cg::greater<float>());
        __syncthreads();
        if (threadIdx.x < BLOCK_SIZE / WARP_SIZE)
            shared_reduce[threadIdx.x] = lm;
        __syncthreads();
        mm = shared_reduce[threadIdx.x / WARP_SIZE];

        // 8-bit: 1 byte per cell; 4-bit packs two cells per byte and needs
        // serialization between odd/even partners. Only QUANT_BITS=8 wired
        // here for now (matches the BilagridFusedAdam 8-bit decision).
        if (inside) QTL::encode_log(packed, idx, log_v, mm);
        if (threadIdx.x == 0) quant_bounds[blockIdx.x] = mm;
    } else {
        if (inside) accum_f[idx] = accum;
    }
}


void fused_bilagrid_tv_adagrad(
    float* grids,
    float*   accum_f,                          // null when quantize
    uint8_t* packed,                           // null when !quantize
    float2*  quant_bounds,                     // null when !quantize
    const float* image_grad,
    const int*   cam_indices,
    int N_grids, int C_batch, int C,
    int L, int H, int W,
    float lr,
    float tv_weight,
    bool quantize,
    cudaStream_t stream
) {
    constexpr int BLOCK_SIZE = 256;
    const int64_t total_cells =
        (int64_t)N_grids * (int64_t)C * (int64_t)L * (int64_t)H * (int64_t)W;
    int blocks = (int)((total_cells + BLOCK_SIZE - 1) / BLOCK_SIZE);
    if (blocks == 0) return;

    #define LAUNCH(CC, QQ) \
        fused_bilagrid_tv_adagrad_kernel<BLOCK_SIZE, CC, QQ, 8> \
            <<<blocks, BLOCK_SIZE, 0, stream>>>( \
                grids, accum_f, packed, quant_bounds, \
                image_grad, cam_indices, \
                N_grids, C_batch, L, H, W, \
                lr, tv_weight)

    if (quantize) {
        if      (C == 12) { LAUNCH(12, true); }
        else if (C == 9)  { LAUNCH(9,  true); }
        else if (C == 3)  { LAUNCH(3,  true); }
        else if (C == 2)  { LAUNCH(2,  true); }
    } else {
        if      (C == 12) { LAUNCH(12, false); }
        else if (C == 9)  { LAUNCH(9,  false); }
        else if (C == 3)  { LAUNCH(3,  false); }
        else if (C == 2)  { LAUNCH(2,  false); }
    }
    #undef LAUNCH
    CHECK_DEVICE_ERROR(cudaGetLastError());
}
