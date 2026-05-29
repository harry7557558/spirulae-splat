#include "BilagridConfig.cuh"


#ifdef PATCHED
__global__ void bilagrid_patched_sample_forward_kernel(
#else
__global__ void bilagrid_uniform_sample_forward_kernel(
#endif
    const float* __restrict__ bilagrid, // [N_grids,L,H,W,12] or [N,L,H,W,12]
    const float* __restrict__ rgb,  // [N,m,h,w,3]
    float* __restrict__ output,  // [N,m,h,w,3]
    int N, int L, int H, int W,
    int m, int h, int w
#ifdef PATCHED
    , int h0, int w0,
    const int* __restrict__ offsets  // [N,m,2]
#endif
#ifndef PATCHED
    , const int* __restrict__ grid_indices  // [N], or nullptr -> identity
#endif
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int total = N * m * h * w;
    if (idx >= total) return;

    int tmp = idx;
    int wi = tmp % w; tmp /= w;
    int hi = tmp % h; tmp /= h;
    int mi = tmp % m; tmp /= m;
    int ni = tmp;
#ifndef PATCHED
    int g_id = grid_indices ? grid_indices[ni] : ni;
#else
    int g_id = ni;
#endif

    // read coords
    int g_offset = (((ni * m + mi) * h + hi) * w + wi) * 3;

    // input and output colors
    float sr = rgb[g_offset+0];
    float sg = rgb[g_offset+1];
    float sb = rgb[g_offset+2];
    float dr = 0.0, dg = 0.0, db = 0.0;

    // grid coords
#ifdef PATCHED
    offsets += (ni * m + mi) * 2;
    float gx = (float)(wi + offsets[0]) / (float)(w0-1);
    float gy = (float)(hi + offsets[1]) / (float)(h0-1);
#else
    float gx = (float)wi / (float)(w-1);
    float gy = (float)hi / (float)(h-1);
#endif
    float gz = kC2G_r * sr + kC2G_g * sg + kC2G_b * sb;
    float x = gx * (W - 1);
    float y = gy * (H - 1);
    float z = gz * (L - 1);

    // find corner indices
    int x0 = (int)floorf(x);
    int y0 = (int)floorf(y);
    int z0 = (int)floorf(z);
    int x1 = min(x0+1, W-1);
    int y1 = min(y0+1, H-1);
    int z1 = z0 + 1;
    z0 = min(max(z0, 0), L-1);
    z1 = min(max(z1, 0), L-1);

    // interpolation parameters
    float fx = x - (float)x0;
    float fy = y - (float)y0;
    float fz = z - (float)z0;

    // Channel-last: the 12 channels of one corner are contiguous in memory,
    // so each ci iter loads from 8 nearby cache lines instead of 8 different
    // (L*H*W)-strided planes.
    int corner_base = g_id * L * H * W * 12;
    int off000 = corner_base + ((z0*H + y0)*W + x0) * 12;
    int off001 = corner_base + ((z0*H + y0)*W + x1) * 12;
    int off010 = corner_base + ((z0*H + y1)*W + x0) * 12;
    int off011 = corner_base + ((z0*H + y1)*W + x1) * 12;
    int off100 = corner_base + ((z1*H + y0)*W + x0) * 12;
    int off101 = corner_base + ((z1*H + y0)*W + x1) * 12;
    int off110 = corner_base + ((z1*H + y1)*W + x0) * 12;
    int off111 = corner_base + ((z1*H + y1)*W + x1) * 12;

    // interpolate and and affine in one loop
    #pragma unroll
    for (int ci = 0; ci < 12; ci++) {
        // fetch 8 corners
        auto v000 = bilagrid[off000 + ci];
        auto v001 = bilagrid[off001 + ci];
        auto v010 = bilagrid[off010 + ci];
        auto v011 = bilagrid[off011 + ci];
        auto v100 = bilagrid[off100 + ci];
        auto v101 = bilagrid[off101 + ci];
        auto v110 = bilagrid[off110 + ci];
        auto v111 = bilagrid[off111 + ci];

        // trilinear interp
        float c00 = v000*(1.0f-fx) + v001*fx;
        float c01 = v010*(1.0f-fx) + v011*fx;
        float c10 = v100*(1.0f-fx) + v101*fx;
        float c11 = v110*(1.0f-fx) + v111*fx;
        float c0 = c00*(1.0f-fy) + c01*fy;
        float c1 = c10*(1.0f-fy) + c11*fy;
        float val = c0*(1.0f-fz) + c1*fz;

        // affine transform
        int si = ci % 4;
        int di = ci / 4;
        (di == 0 ? dr : di == 1 ? dg : db) += val *
            (si==0 ? sr : si==1 ? sg : si==2 ? sb : 1.0f);
    }

    output[g_offset+0] = isfinite(dr) ? dr : 0.5f;
    output[g_offset+1] = isfinite(dg) ? dg : 0.5f;
    output[g_offset+2] = isfinite(db) ? db : 0.5f;
}
