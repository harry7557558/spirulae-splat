#include "BilagridConfig.cuh"

// Depth bilagrid backward v2 = thread-per-pixel SCATTER. Depth grids are
// GT-side (needs_image_grad = false: the engine passes v_depth = nullptr and
// discards the input gradient), so v2 only scatters the grid gradient -- no
// v_depth write, no gz chain. Reads grid VALUES from the g_id slot; scatters
// grad into the per-batch buffer indexed by ni. See BilagridBackwardSelection.md
// and BilagridDepthSampleBwdV1_kernel.cuh (the transform math is mirrored).
__global__ void bilagrid_depth_uniform_sample_backward_v2_kernel(
    BilagridReader bilagrid,              // [N_grids,L,H,W,2]
    const float* __restrict__ depth,      // [N,h,w]
    const float* __restrict__ scalars,    // [N_grids]
    const float* __restrict__ v_output,   // [N,h,w]
    float* __restrict__ v_bilagrid,       // [N,L,H,W,2] (per-batch, accumulated)
    int N, int L, int H, int W, int h, int w,
    const int* __restrict__ grid_indices  // [N], or nullptr -> identity
) {
    int idx = blockIdx.x * kBilagridBwdV1RgbThreads + threadIdx.x;
    int total = N * h * w;
    if (idx >= total) return;
    int tmp = idx;
    int wi = tmp % w; tmp /= w;
    int hi = tmp % h; tmp /= h;
    int ni = tmp;
    int g_id = grid_indices ? grid_indices[ni] : ni;
    const float scalar = scalars[g_id];
    int g_off = (ni * h + hi) * w + wi;
    float sr = depth[g_off] * scalar;
    float dr = v_output[g_off];

    float x = (float)wi / (float)(w - 1) * (float)(W - 1);
    float y = (float)hi / (float)(h - 1) * (float)(H - 1);
    float gz_raw = sr / (sr + 1.0f);
    float gz = fminf(fmaxf(gz_raw, 0.0f), 1.0f);
    float z = gz * (L - 1);
    int x0 = floorf(x), y0 = floorf(y), z0 = floorf(z);
    int x1 = min(x0 + 1, W - 1);
    int y1 = min(y0 + 1, H - 1);
    int z1 = min(z0 + 1, L - 1);
    float fx = x - x0, fy = y - y0, fz = z - z0;

    const float ww[8] = {
        (1 - fx) * (1 - fy) * (1 - fz), fx * (1 - fy) * (1 - fz),
        (1 - fx) * fy * (1 - fz),       fx * fy * (1 - fz),
        (1 - fx) * (1 - fy) * fz,       fx * (1 - fy) * fz,
        (1 - fx) * fy * fz,             fx * fy * fz,
    };
    const int cc[8] = {
        ((z0 * H + y0) * W + x0) * 2, ((z0 * H + y0) * W + x1) * 2,
        ((z0 * H + y1) * W + x0) * 2, ((z0 * H + y1) * W + x1) * 2,
        ((z1 * H + y0) * W + x0) * 2, ((z1 * H + y0) * W + x1) * 2,
        ((z1 * H + y1) * W + x0) * 2, ((z1 * H + y1) * W + x1) * 2,
    };
    int grid_base = g_id * L * H * W * 2;
    int out_base  = ni  * L * H * W * 2;

    // interpolate the 2 params (read grid values from the g_id slot)
    float affine[2];
    #pragma unroll
    for (int ci = 0; ci < 2; ci++) {
        float val = 0.0f;
        #pragma unroll
        for (int corner = 0; corner < 8; corner++)
            val += bilagrid[grid_base + cc[corner] + ci] * ww[corner];
        affine[ci] = (ci == 0) ? __expf(val) : val;
    }
    float log_sr = __logf(sr);
    dr *= __expf(affine[0] * log_sr + affine[1]);
    const float grad_weight[2] = { log_sr * dr * affine[0], dr };

    // scatter param grad to the 8 corners (per-batch ni buffer)
    #pragma unroll
    for (int corner = 0; corner < 8; corner++) {
        int wbase = out_base + cc[corner];
        #pragma unroll
        for (int ci = 0; ci < 2; ci++) {
            float aa = ww[corner] * grad_weight[ci];
            if (aa != 0.0f && isfinite(aa)) atomicAdd(v_bilagrid + wbase + ci, aa);
        }
    }
}
