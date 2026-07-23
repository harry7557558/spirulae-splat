#include "kernels/bilagrid/BilagridConfig.cuh"

#include "kernels/bilagrid/BilagridPpispMath.cuh"

using namespace bilagrid_ppisp;


#ifdef PACKED
__global__ void bilagrid_ppisp_packed_sample_forward_kernel(
    BilagridReader bilagrid, // [N,L,H,W,9]
    const int64_t* __restrict__ image_indices,  // [nnz]
    const float* __restrict__ coords,  // [nnz,2]
    const float* __restrict__ rgb_in,  // [nnz,3]
    float* __restrict__ rgb_out,  // [nnz,3]
    int N, int L, int H, int W,
    int nnz
#else
__global__ void bilagrid_ppisp_sample_forward_kernel(
    BilagridReader bilagrid, // [N,L,H,W,9]
    const float* __restrict__ coords,  // [N,m,h,w,2]
    const float* __restrict__ rgb_in,  // [N,m,h,w,3]
    float* __restrict__ rgb_out,  // [N,m,h,w,3]
    int N, int L, int H, int W,
    int m, int h, int w
#endif
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
#ifdef PACKED
    if (idx >= nnz) return;
    int ni = image_indices[idx];
    int g_offset = idx;
#else
    if (idx >= N * m * h * w) return;
    int tmp = idx;
    int wi = tmp % w; tmp /= w;
    int hi = tmp % h; tmp /= h;
    int mi = tmp % m; tmp /= m;
    int ni = tmp;
    int g_offset = (((ni * m + mi) * h + hi) * w + wi);
#endif

    // load colors
    float sr = rgb_in[3*g_offset+0];
    float sg = rgb_in[3*g_offset+1];
    float sb = rgb_in[3*g_offset+2];

    // read coords
    float gx = coords[2*g_offset+0];
    float gy = coords[2*g_offset+1];
    float gz = kC2G_r * sr + kC2G_g * sg + kC2G_b * sb;
    // Clamp gz to [0,1] -- equivalent to clipping z to [0, L-1]. Matches the
    // backward branch (which clamps), so fz here is in [0,1] same as bwd's.
    gz = fminf(fmaxf(gz, 0.0f), 1.0f);
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

    // interpolate

    float exposure_param;
    ColorPPISPParams color_params;

    // Channel-last: the 9 channels of one corner are contiguous in memory.
    int corner_base = ni * L * H * W * 9;
    int off000 = corner_base + ((z0*H + y0)*W + x0) * 9;
    int off001 = corner_base + ((z0*H + y0)*W + x1) * 9;
    int off010 = corner_base + ((z0*H + y1)*W + x0) * 9;
    int off011 = corner_base + ((z0*H + y1)*W + x1) * 9;
    int off100 = corner_base + ((z1*H + y0)*W + x0) * 9;
    int off101 = corner_base + ((z1*H + y0)*W + x1) * 9;
    int off110 = corner_base + ((z1*H + y1)*W + x0) * 9;
    int off111 = corner_base + ((z1*H + y1)*W + x1) * 9;

    #pragma unroll
    for (int ci = 0; ci < 9; ci++) {
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

    rgb_out[3*g_offset+0] = isfinite(rgb.x) ? rgb.x : 0.5f;
    rgb_out[3*g_offset+1] = isfinite(rgb.y) ? rgb.y : 0.5f;
    rgb_out[3*g_offset+2] = isfinite(rgb.z) ? rgb.z : 0.5f;
}

