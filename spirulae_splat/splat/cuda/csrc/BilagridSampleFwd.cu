#include "BilagridConfig.cuh"


__global__ void bilagrid_sample_forward_kernel(
    const float* __restrict__ bilagrid, // [N,12,L,H,W]
    const float* __restrict__ coords,  // [N,m,h,w,2]
    const float* __restrict__ rgb,  // [N,m,h,w,3]
    float* __restrict__ output,  // [N,m,h,w,3]
    int N, int L, int H, int W,
    int m, int h, int w
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int total = N * m * h * w;
    if (idx >= total) return;

    int tmp = idx;
    int wi = tmp % w; tmp /= w;
    int hi = tmp % h; tmp /= h;
    int mi = tmp % m; tmp /= m;
    int ni = tmp;

    // load colors
    int g_offset = (((ni * m + mi) * h + hi) * w + wi);
    float sr = rgb[3*g_offset+0];
    float sg = rgb[3*g_offset+1];
    float sb = rgb[3*g_offset+2];

    // read coords
    float gx = coords[2*g_offset+0];
    float gy = coords[2*g_offset+1];
    float gz = kC2G_r * sr + kC2G_g * sg + kC2G_b * sb;
    float x = gx * (W - 1);
    float y = gy * (H - 1);
    float z = gz * (L - 1);

    // find corner indices
    int x0 = (int)floorf(x);
    int y0 = (int)floorf(y);
    int z0 = (int)floorf(z);
    int x1 = x0 + 1;
    int y1 = y0 + 1;
    int z1 = z0 + 1;
    x0 = min(max(x0, 0), W-1);
    x1 = min(max(x1, 0), W-1);
    y0 = min(max(y0, 0), H-1);
    y1 = min(max(y1, 0), H-1);
    z0 = min(max(z0, 0), L-1);
    z1 = min(max(z1, 0), L-1);

    // interpolation parameters
    float fx = x - (float)x0;
    float fy = y - (float)y0;
    float fz = z - (float)z0;

    // output colors
    float dr = 0.0, dg = 0.0, db = 0.0;

    // interpolate and and affine in one loop
    #pragma unroll
    for (int ci = 0; ci < 12; ci++) {
        // base pointer for this volume
        const float* vol = &bilagrid[((ni*12 + ci)*L*H*W)];

        // fetch 8 corners
        auto v000 = vol[(z0*H+y0)*W+x0];
        auto v001 = vol[(z0*H+y0)*W+x1];
        auto v010 = vol[(z0*H+y1)*W+x0];
        auto v011 = vol[(z0*H+y1)*W+x1];
        auto v100 = vol[(z1*H+y0)*W+x0];
        auto v101 = vol[(z1*H+y0)*W+x1];
        auto v110 = vol[(z1*H+y1)*W+x0];
        auto v111 = vol[(z1*H+y1)*W+x1];

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

    output[3*g_offset+0] = dr;
    output[3*g_offset+1] = dg;
    output[3*g_offset+2] = db;
}


void bilagrid_sample_forward(
    const float* bilagrid,
    const float* coords,
    const float* rgb,
    float* output,
    int N, int L, int H, int W,
    int m, int h, int w,
    cudaStream_t stream
) {
    int total = N * m * h * w;
    int threads = 256;
    int blocks = (total + threads - 1) / threads;
    bilagrid_sample_forward_kernel<<<blocks, threads, 0, stream>>>(
        bilagrid, coords, rgb, output,
        N, L, H, W, m, h, w
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
    // cudaDeviceSynchronize();
}
