#include "kernels/bilagrid/BilagridConfig.cuh"

#include "kernels/bilagrid/BilagridAxisAngle.cuh"


#ifdef PATCHED
__global__ void bilagrid_normal_patched_sample_forward_kernel(
#else
__global__ void bilagrid_normal_uniform_sample_forward_kernel(
#endif
    BilagridReader bilagrid, // [N,L,H,W,3]
#ifdef PATCHED
    const float* __restrict__ normal_in,  // [N,m,h,w,3]
    float* __restrict__ normal_out,  // [N,m,h,w,3]
#else
    const float* __restrict__ normal_in,  // [N,h,w,3]
    float* __restrict__ normal_out,  // [N,h,w,3]
#endif
    int N, int L, int H, int W,
#ifdef PATCHED
    int m,
#endif
    int h, int w
#ifdef PATCHED
    , int h0, int w0,
    const int* __restrict__ offsets  // [N,m,2]
#endif
#ifndef PATCHED
    , const int* __restrict__ grid_indices  // [N], or nullptr -> identity
#endif
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
#ifdef PATCHED
    int total = N * m * h * w;
#else
    int total = N * h * w;
#endif
    if (idx >= total) return;

    int tmp = idx;
    int wi = tmp % w; tmp /= w;
    int hi = tmp % h; tmp /= h;
#ifdef PATCHED
    int mi = tmp % m; tmp /= m;
#endif
    int ni = tmp;
#ifndef PATCHED
    int g_id = grid_indices ? grid_indices[ni] : ni;
#else
    int g_id = ni;
#endif

    // read coords
#ifdef PATCHED
    int g_offset = (((ni * m + mi) * h + hi) * w + wi) * 3;
#else
    int g_offset = ((ni * h + hi) * w + wi) * 3;
#endif

    // input and output colors
    float sr = normal_in[g_offset+0];
    float sg = normal_in[g_offset+1];
    float sb = normal_in[g_offset+2];

    // "No normal here" (core/Interpolation.cuh's gt_normal_valid) must survive
    // the grid: normalizing it would hand the loss a plausible direction to
    // train the sky and the unwarped corners against.
    if (!(sr + sg + sb > -2.366f && sr*sr + sg*sg + sb*sb > 0.25f)) {
        normal_out[g_offset+0] = sr;
        normal_out[g_offset+1] = sg;
        normal_out[g_offset+2] = sb;
        return;
    }

    float inv_norm = rsqrtf(sr*sr + sg*sg + sb*sb + 1e-20f);

    // grid coords
#ifdef PATCHED
    offsets += (ni * m + mi) * 2;
    float gx = (float)(wi + offsets[0]) / (float)(w0-1);
    float gy = (float)(hi + offsets[1]) / (float)(h0-1);
#else
    float gx = (float)wi / (float)(w-1);
    float gy = (float)hi / (float)(h-1);
#endif
    // float gz = acosf(fminf(fmaxf(sb * inv_norm, -1.0f), 1.0f)) * (1.0f / (float)M_PI);
    float gz = 0.5f + 0.5f * sb * inv_norm;
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
    int x1 = min(x0+1, W-1);
    int y1 = min(y0+1, H-1);
    int z1 = min(z0+1, L-1);

    // interpolation parameters
    float fx = x - (float)x0;
    float fy = y - (float)y0;
    float fz = z - (float)z0;

    float3 axis_angle;

    // Channel-last: the 3 channels of one corner are contiguous in memory.
    int corner_base = g_id * L * H * W * 3;
    int off000 = corner_base + ((z0*H + y0)*W + x0) * 3;
    int off001 = corner_base + ((z0*H + y0)*W + x1) * 3;
    int off010 = corner_base + ((z0*H + y1)*W + x0) * 3;
    int off011 = corner_base + ((z0*H + y1)*W + x1) * 3;
    int off100 = corner_base + ((z1*H + y0)*W + x0) * 3;
    int off101 = corner_base + ((z1*H + y0)*W + x1) * 3;
    int off110 = corner_base + ((z1*H + y1)*W + x0) * 3;
    int off111 = corner_base + ((z1*H + y1)*W + x1) * 3;

    // interpolate and and affine in one loop
    #pragma unroll
    for (int ci = 0; ci < 3; ci++) {
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

        (ci == 0 ? axis_angle.x
            : ci == 1 ? axis_angle.y
            : axis_angle.z) = val;
    }

    // apply normal
    float3 normal = {sr*inv_norm, sg*inv_norm, sb*inv_norm};
    normal = axis_angle_rotate_fwd(axis_angle, normal);

    // save results
    normal_out[g_offset+0] = isfinite(normal.x) ? normal.x : 0.5f;
    normal_out[g_offset+1] = isfinite(normal.y) ? normal.y : 0.5f;
    normal_out[g_offset+2] = isfinite(normal.z) ? normal.z : 0.5f;
}
