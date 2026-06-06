// Bilagrid RGB color-shift regularizer (design 1).
//
// Goal: penalize the dataset-wide mean direction of the bilagrid correction
//
//     R = w * || E_data[ mean_per_pixel( sign(post - pre) ) ] ||^2
//
// where post = bilagrid_forward(pre, grid). Injects the surrogate gradient
//
//     v_render_rgb_post[p, c] += (2 * w / (N * (1 - beta^t))) * ema[c] * sign(c_p[c])
//
// on the POST-bilagrid loss-side gradient buffer, BEFORE bilagrid_*_backward
// runs. The bilagrid vjp then routes this contribution to (a) the grid
// parameter gradient (what we actually want regularized) and (b) a small leak
// into the pre-bilagrid (splat-side) gradient via the bilagrid's local
// Jacobian. The leak is ~0 when the bilagrid is near-identity and grows as the
// bilagrid acquires non-trivial gain - typically harmless given the splats are
// being driven by the main RGB loss simultaneously.
//
// EMA state is per-channel (3 floats on device). Updated each step as
//
//     ema = beta * ema + (1 - beta) * batch_mean
//
// with bias-correction (1 - beta^t) computed on host -- no D->H readbacks of
// the EMA itself, ever. The injection + per-channel batch-sum reduction is
// fused into one kernel; a tiny 3-thread kernel applies the EMA update.

#include "BilagridConfig.cuh"
#include "Tensor.h"
#include <cuda_runtime.h>
#include <cmath>


template<int BLOCK_SIZE>
__global__ void bilagrid_rgb_shift_reg_inject_kernel(
    float3* __restrict__ v_render_rgb,        // [N] gradient on POST render
    const float3* __restrict__ post_rgb,       // [N]
    const float3* __restrict__ pre_rgb,        // [N]
    const float* __restrict__ ema,              // [3] running EMA (NOT bias-corrected)
    float* __restrict__ batch_sum_out,          // [3] sum of sign(c) over batch
    int N,
    float reg_coef                              // 2 * w / (N * (1 - beta^t))
) {
    int idx = (int)blockIdx.x * BLOCK_SIZE + (int)threadIdx.x;

    // One shared cache of EMA[0..2] and a per-block partial sum, written to
    // global with one atomicAdd per channel at the end of the block (avoids
    // N_blocks * N_threads global atomics).
    __shared__ float ema_s[3];
    __shared__ float block_sum[3];
    if ((int)threadIdx.x < 3) {
        ema_s[threadIdx.x] = ema[threadIdx.x];
        block_sum[threadIdx.x] = 0.0f;
    }
    __syncthreads();

    float sx = 0.0f, sy = 0.0f, sz = 0.0f;

    if (idx < N) {
        float3 post = post_rgb[idx];
        float3 pre  = pre_rgb[idx];
        // sign(x) -> {-1, 0, +1}; (x>0)-(x<0) is branchless.
        sx = (float)((post.x > pre.x) - (post.x < pre.x));
        sy = (float)((post.y > pre.y) - (post.y < pre.y));
        sz = (float)((post.z > pre.z) - (post.z < pre.z));

        // Inject regularizer gradient on POST: + reg_coef * ema_hat[c] * sign(c[c]).
        float3 v = v_render_rgb[idx];
        v.x += reg_coef * ema_s[0] * sx;
        v.y += reg_coef * ema_s[1] * sy;
        v.z += reg_coef * ema_s[2] * sz;
        v_render_rgb[idx] = v;

        // Shared-memory atomicAdd is fast on modern GPUs; the 3-channel
        // contention is tiny and beats a full warp-then-block reduction here.
        atomicAdd(&block_sum[0], sx);
        atomicAdd(&block_sum[1], sy);
        atomicAdd(&block_sum[2], sz);
    }
    __syncthreads();

    // One thread per channel pushes the block sum to the global per-channel
    // accumulator.
    if ((int)threadIdx.x < 3) {
        atomicAdd(&batch_sum_out[threadIdx.x], block_sum[threadIdx.x]);
    }
}


__global__ void bilagrid_rgb_shift_reg_update_kernel(
    float* __restrict__ ema,           // [3] in-place
    float* __restrict__ batch_sum,      // [3] in-place (reset to 0)
    float beta,
    float inv_n_pixels                  // 1.0f / N
) {
    int c = (int)threadIdx.x;
    if (c < 3) {
        float batch_mean = batch_sum[c] * inv_n_pixels;
        ema[c] = beta * ema[c] + (1.0f - beta) * batch_mean;
        batch_sum[c] = 0.0f;
    }
}


// Host dispatcher. Caller (the bilagrid bwd hook) is responsible for
// (1) allocating + zero-initing ema[3] and batch_sum[3] on first use, and
// (2) tracking `step` (number of EMA updates performed so far, pre-this-call).
//
// step == 0 on the first call -> ema is zero -> injected gradient is zero
// (and the bias-correction divisor is also zero in the limit; we clamp it to
// avoid div-by-zero, but it doesn't matter because ema is zero anyway).
//
// All work runs on `stream`; no D->H copies, no host sync.
void bilagrid_rgb_shift_reg_step(
    float* v_render_rgb,
    const float* post_rgb,
    const float* pre_rgb,
    float* shift_reg_ema,
    float* shift_reg_batch_sum,
    int N_pixels,
    float weight,
    float beta,
    int64_t step,
    cudaStream_t stream
) {
    if (weight <= 0.0f || N_pixels <= 0) return;
    if (!shift_reg_ema || !shift_reg_batch_sum) return;

    // Standard Adam-style bias correction. Use (step + 1) since after this
    // call's update we'll have done that many EMA updates total.
    float bc = 1.0f - powf(beta, (float)(step + 1));
    if (bc < 1e-30f) bc = 1.0f;   // safety; ema is zero in this regime anyway

    float reg_coef = 2.0f * weight / ((float)N_pixels * bc);

    constexpr int BLOCK_SIZE = 256;
    int blocks = (N_pixels + BLOCK_SIZE - 1) / BLOCK_SIZE;

    bilagrid_rgb_shift_reg_inject_kernel<BLOCK_SIZE>
        <<<blocks, BLOCK_SIZE, 0, stream>>>(
            (float3*)v_render_rgb,
            (const float3*)post_rgb,
            (const float3*)pre_rgb,
            shift_reg_ema,
            shift_reg_batch_sum,
            N_pixels,
            reg_coef);

    bilagrid_rgb_shift_reg_update_kernel<<<1, 4, 0, stream>>>(
        shift_reg_ema, shift_reg_batch_sum, beta, 1.0f / (float)N_pixels);

    CHECK_DEVICE_ERROR(cudaGetLastError());
}
