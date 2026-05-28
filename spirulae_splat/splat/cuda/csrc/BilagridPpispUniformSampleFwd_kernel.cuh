#include "BilagridConfig.cuh"

#include "BilagridPpispMath.cuh"

using namespace bilagrid_ppisp;


#ifdef PATCHED
__global__ void bilagrid_ppisp_patched_sample_forward_kernel(
#else
__global__ void bilagrid_ppisp_uniform_sample_forward_kernel(
#endif
    const float* __restrict__ bilagrid, // [N,9,L,H,W]
    const float* __restrict__ rgb_in,  // [N,m,h,w,3]
    float* __restrict__ rgb_out,  // [N,m,h,w,3]
    int N, int L, int H, int W,
    int m, int h, int w
#ifdef PATCHED
    , int h0, int w0,
    const int* __restrict__ offsets  // [N,m,2]
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

    // read coords
    int g_offset = (((ni * m + mi) * h + hi) * w + wi) * 3;

    // input and output colors
    float sr = rgb_in[g_offset+0];
    float sg = rgb_in[g_offset+1];
    float sb = rgb_in[g_offset+2];

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

    float exposure_param;
    ColorPPISPParams color_params;

    // interpolate and and affine in one loop
    #pragma unroll
    for (int ci = 0; ci < 9; ci++) {
        // base pointer for this volume
        int base = (ni*9 + ci)*L*H*W;

        // fetch 8 corners
        auto v000 = bilagrid[base+(z0*H+y0)*W+x0];
        auto v001 = bilagrid[base+(z0*H+y0)*W+x1];
        auto v010 = bilagrid[base+(z0*H+y1)*W+x0];
        auto v011 = bilagrid[base+(z0*H+y1)*W+x1];
        auto v100 = bilagrid[base+(z1*H+y0)*W+x0];
        auto v101 = bilagrid[base+(z1*H+y0)*W+x1];
        auto v110 = bilagrid[base+(z1*H+y1)*W+x0];
        auto v111 = bilagrid[base+(z1*H+y1)*W+x1];

        // trilinear interp
        float c00 = v000*(1.0f-fx) + v001*fx;
        float c01 = v010*(1.0f-fx) + v011*fx;
        float c10 = v100*(1.0f-fx) + v101*fx;
        float c11 = v110*(1.0f-fx) + v111*fx;
        float c0 = c00*(1.0f-fy) + c01*fy;
        float c1 = c10*(1.0f-fy) + c11*fy;
        float val = c0*(1.0f-fz) + c1*fz;

        (ci == 0 ? exposure_param :
            ci == 1 ? color_params.b.x :
            ci == 2 ? color_params.b.y :
            ci == 3 ? color_params.r.x :
            ci == 4 ? color_params.r.y :
            ci == 5 ? color_params.g.x :
            ci == 6 ? color_params.g.y :
            ci == 7 ? color_params.n.x :
            color_params.n.y) = val;
    }

    // apply PPISP
    float3 rgb = {sr, sg, sb};
    apply_exposure(rgb, exposure_param, rgb);
    apply_color_correction_ppisp(rgb, &color_params, rgb);

    // save results
    rgb_out[g_offset+0] = isfinite(rgb.x) ? rgb.x : 0.5f;
    rgb_out[g_offset+1] = isfinite(rgb.y) ? rgb.y : 0.5f;
    rgb_out[g_offset+2] = isfinite(rgb.z) ? rgb.z : 0.5f;
}
