#include "kernels/bilagrid/BilagridConfig.cuh"
#include "kernels/bilagrid/BilagridAxisAngle.cuh"

// Normal bilagrid backward v2 = thread-per-pixel SCATTER. Normal grids are
// GT-side (needs_image_grad = false: engine passes v_rgb = nullptr, input grad
// discarded), so v2 only scatters the grid gradient -- no v_rgb write, no gz
// chain. Reads grid VALUES from the g_id slot; scatters grad into the per-batch
// buffer indexed by ni. Transform-backward is the hand-written Rodrigues form
// (axis_angle_rotate_bwd), mirroring BilagridNormalSampleBwdV1_kernel.cuh.
__global__ void bilagrid_normal_uniform_sample_backward_v2_kernel(
    BilagridReader bilagrid,              // [N_grids,L,H,W,3]
    const float* __restrict__ normal_in,  // [N,h,w,3]
    const float* __restrict__ v_output,   // [N,h,w,3]
    float* __restrict__ v_bilagrid,       // [N,L,H,W,3] (per-batch, accumulated)
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
    int g_off = ((ni * h + hi) * w + wi) * 3;

    float sr = normal_in[g_off + 0];
    float sg = normal_in[g_off + 1];
    float sb = normal_in[g_off + 2];
    float inv_norm = rsqrtf(sr * sr + sg * sg + sb * sb + 1e-20f);
    float3 unit_normal = {sr * inv_norm, sg * inv_norm, sb * inv_norm};
    float dr = v_output[g_off + 0];
    float dg = v_output[g_off + 1];
    float db = v_output[g_off + 2];

    float x = (float)wi / (float)(w - 1) * (float)(W - 1);
    float y = (float)hi / (float)(h - 1) * (float)(H - 1);
    float gz_raw = 0.5f + 0.5f * sb * inv_norm;
    float gz = fminf(fmaxf(gz_raw, 0.0f), 1.0f);
    float z = gz * (L - 1);
    int x0 = floorf(x), y0 = floorf(y), z0 = floorf(z);
    int x1 = min(x0 + 1, W - 1);
    int y1 = min(y0 + 1, H - 1);
    int z1 = min(z0 + 1, L - 1);
    float fx = x - x0, fy = y - y0, fz = z - z0;

    int grid_base = g_id * L * H * W * 3;
    int out_base  = ni  * L * H * W * 3;
    const int cc[8] = {
        ((z0 * H + y0) * W + x0) * 3, ((z0 * H + y0) * W + x1) * 3,
        ((z0 * H + y1) * W + x0) * 3, ((z0 * H + y1) * W + x1) * 3,
        ((z1 * H + y0) * W + x0) * 3, ((z1 * H + y0) * W + x1) * 3,
        ((z1 * H + y1) * W + x0) * 3, ((z1 * H + y1) * W + x1) * 3,
    };

    // interpolate the 3 axis-angle params (nested lerp, from the g_id slot)
    float3 axis_angle;
    #pragma unroll
    for (int ci = 0; ci < 3; ci++) {
        float v000 = bilagrid[grid_base + cc[0] + ci];
        float v001 = bilagrid[grid_base + cc[1] + ci];
        float v010 = bilagrid[grid_base + cc[2] + ci];
        float v011 = bilagrid[grid_base + cc[3] + ci];
        float v100 = bilagrid[grid_base + cc[4] + ci];
        float v101 = bilagrid[grid_base + cc[5] + ci];
        float v110 = bilagrid[grid_base + cc[6] + ci];
        float v111 = bilagrid[grid_base + cc[7] + ci];
        float c00 = v000 * (1.0f - fx) + v001 * fx;
        float c01 = v010 * (1.0f - fx) + v011 * fx;
        float c10 = v100 * (1.0f - fx) + v101 * fx;
        float c11 = v110 * (1.0f - fx) + v111 * fx;
        float c0 = c00 * (1.0f - fy) + c01 * fy;
        float c1 = c10 * (1.0f - fy) + c11 * fy;
        float val = c0 * (1.0f - fz) + c1 * fz;
        (ci == 0 ? axis_angle.x : ci == 1 ? axis_angle.y : axis_angle.z) = val;
    }

    float3 grad_axis_angle = {0.0f, 0.0f, 0.0f};
    float3 grad_unit_normal = {0.0f, 0.0f, 0.0f};
    axis_angle_rotate_bwd(axis_angle, unit_normal, {dr, dg, db},
                          grad_axis_angle, grad_unit_normal);
    const float grad_weight[3] = {grad_axis_angle.x, grad_axis_angle.y,
                                  grad_axis_angle.z};

    const float wcx[2] = {1.0f - fx, fx};
    const float wcy[2] = {1.0f - fy, fy};
    const float wcz[2] = {1.0f - fz, fz};
    #pragma unroll
    for (int corner = 0; corner < 8; corner++) {
        float wc = wcx[corner & 1] * wcy[(corner >> 1) & 1] * wcz[(corner >> 2) & 1];
        int wbase = out_base + cc[corner];
        #pragma unroll
        for (int ci = 0; ci < 3; ci++) {
            float aa = wc * grad_weight[ci];
            if (aa != 0.0f && isfinite(aa)) atomicAdd(v_bilagrid + wbase + ci, aa);
        }
    }
}
