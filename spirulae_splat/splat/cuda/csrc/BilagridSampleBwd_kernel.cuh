#include "BilagridConfig.cuh"


// template<int compute_coords_grad>
#ifdef COMPUTE_COORDS_GRAD
__global__ void bilagrid_sample_backward_kernel_cg(
#else
__global__ void bilagrid_sample_backward_kernel(
#endif
    const float* __restrict__ bilagrid,  // [N,L,H,W,12]
    const float* __restrict__ coords,  // [N,m,h,w,2]
    const float* __restrict__ rgb,  // [N,m,h,w,3]
    const float* __restrict__ v_output,  // [N,m,h,w,3]
    float* __restrict__ v_bilagrid,  // [N,L,H,W,12]
    #ifdef COMPUTE_COORDS_GRAD
    float* __restrict__ v_coords,  // [N,m,h,w,2]
    #endif
    float* __restrict__ v_rgb,  // [N,m,h,w,3]
    int N, int L, int H, int W,
    int m, int h, int w
) {
    #if 0
    // faster when coords are random
    int wi = blockIdx.x * blockDim.x + threadIdx.x;
    int hi = blockIdx.y * blockDim.y + threadIdx.y;
    #else
    // faster when coords is a regular grid
    // reduces number of threads writing to the same address at a time in atomicAdd
    int wi = threadIdx.x * ((w+blockDim.x-1) / blockDim.x) + blockIdx.x;
    int hi = threadIdx.y * ((h+blockDim.y-1) / blockDim.y) + blockIdx.y;
    #endif

    int idx = blockIdx.z * blockDim.z + threadIdx.z;
    bool inside = (wi < w && hi < h && idx < (N*m));
    if (!inside) return;
    int mi = idx % m;
    int ni = idx / m;

    // grid coords
    int g_off = (((ni*m + mi)*h + hi)*w + wi);
    float sr = rgb[3*g_off+0], sg = rgb[3*g_off+1], sb = rgb[3*g_off+2];
    float gx = coords[2*g_off+0];
    float gy = coords[2*g_off+1];
    float gz = kC2G_r * sr + kC2G_g * sg + kC2G_b * sb;
    float x = gx * (W - 1);
    float y = gy * (H - 1);
    float z = gz * (L - 1);

    // floor + ceil, clamped
    int x0 = floorf(x), y0 = floorf(y), z0 = floorf(z);
    int x1 = x0 + 1, y1 = y0 + 1, z1 = z0 + 1;
    x0 = min(max(x0,0), W-1); x1 = min(max(x1,0), W-1);
    y0 = min(max(y0,0), H-1); y1 = min(max(y1,0), H-1);
    z0 = min(max(z0,0), L-1); z1 = min(max(z1,0), L-1);

    // fractional parts
    float fx = x - x0, fy = y - y0, fz = z - z0;
    float f000 = (1-fx)*(1-fy)*(1-fz);
    float f001 = fx*(1-fy)*(1-fz);
    float f010 = (1-fx)*fy*(1-fz);
    float f011 = fx*fy*(1-fz);
    float f100 = (1-fx)*(1-fy)*fz;
    float f101 = fx*(1-fy)*fz;
    float f110 = (1-fx)*fy*fz;
    float f111 = fx*fy*fz;

    // read rgb coeffs and upstream gradient
    float dr = v_output[3*g_off+0];
    float dg = v_output[3*g_off+1];
    float db = v_output[3*g_off+2];
    float vr = 0.0, vg = 0.0, vb = 0.0;

    // Channel-last corner offsets shared across the C-channel loop.
    int corner_base = ni * L * H * W * 12;
    int off000 = corner_base + ((z0*H + y0)*W + x0) * 12;
    int off001 = corner_base + ((z0*H + y0)*W + x1) * 12;
    int off010 = corner_base + ((z0*H + y1)*W + x0) * 12;
    int off011 = corner_base + ((z0*H + y1)*W + x1) * 12;
    int off100 = corner_base + ((z1*H + y0)*W + x0) * 12;
    int off101 = corner_base + ((z1*H + y0)*W + x1) * 12;
    int off110 = corner_base + ((z1*H + y1)*W + x0) * 12;
    int off111 = corner_base + ((z1*H + y1)*W + x1) * 12;

    // accumulate bilagrid gradient over 12 channels
    #pragma unroll
    for (int ci = 0; ci < 12; ++ci) {
        // weight from rgb channel
        int si = ci % 4, di = ci / 4;
        float r_coeff = (si==0 ? sr : si==1 ? sg : si==2 ? sb : 1.f);
        float gout = (di==0 ? dr : di==1 ? dg : db);
        float grad_weight = r_coeff * gout;

        // scatter back into the eight corners
        // accounts for >90% of run time, needs optimization
        atomicAdd(v_bilagrid + off000 + ci, f000 * grad_weight);
        atomicAdd(v_bilagrid + off001 + ci, f001 * grad_weight);
        atomicAdd(v_bilagrid + off010 + ci, f010 * grad_weight);
        atomicAdd(v_bilagrid + off011 + ci, f011 * grad_weight);
        atomicAdd(v_bilagrid + off100 + ci, f100 * grad_weight);
        atomicAdd(v_bilagrid + off101 + ci, f101 * grad_weight);
        atomicAdd(v_bilagrid + off110 + ci, f110 * grad_weight);
        atomicAdd(v_bilagrid + off111 + ci, f111 * grad_weight);

        // gradient w.r.t. rgb coefficients
        if (si < 3) {
            float val =
                ( ( (bilagrid[off000 + ci]*(1-fx) + bilagrid[off001 + ci]*fx)*(1-fy)
                    + (bilagrid[off010 + ci]*(1-fx) + bilagrid[off011 + ci]*fx)*fy )*(1-fz)
                + ( (bilagrid[off100 + ci]*(1-fx) + bilagrid[off101 + ci]*fx)*(1-fy)
                    + (bilagrid[off110 + ci]*(1-fx) + bilagrid[off111 + ci]*fx)*fy )*fz
                );
            (si == 0 ? vr : si == 1 ? vg : vb) += val * gout;
        }
    }

    // spatial derivatives for coords
    // dw000/dx = -(1-fy)*(1-fz), dw001/dx = +(1-fy)*(1-fz), etc...
#ifdef COMPUTE_COORDS_GRAD
    float dwdx[8] = {
        -(1-fy)*(1-fz),  (1-fy)*(1-fz),
        -fy*(1-fz),      fy*(1-fz),
        -(1-fy)*fz,      (1-fy)*fz,
        -fy*fz,          fy*fz
    };
    float dwdy[8] = {
        -(1-fx)*(1-fz), -fx*(1-fz),
         (1-fx)*(1-fz),  fx*(1-fz),
        -(1-fx)*fz,     -fx*fz,
         (1-fx)*fz,      fx*fz
    };
#endif
    float dwdz[8] = {
        -(1-fx)*(1-fy), -fx*(1-fy),
        -(1-fx)*fy,     -fx*fy,
         (1-fx)*(1-fy),  fx*(1-fy),
         (1-fx)*fy,      fx*fy
    };

    const int corner_offs[8] = {
        off000, off001, off010, off011,
        off100, off101, off110, off111
    };

    // accumulate gradient into coords (chain through bilagrid values and rgb)
    #ifdef COMPUTE_COORDS_GRAD
    float gx_grad = 0.f, gy_grad = 0.f;
    #endif
    float gz_grad = 0.f;
    #pragma unroll
    for (int corner = 0; corner < 8; ++corner) {
        int base = corner_offs[corner];
        float trilerp = 0.f;
        // gather the corresponding bilagrid value for each of the 12 channels
        #pragma unroll
        for (int ci = 0; ci < 12; ++ci) {
            float v = bilagrid[base + ci];
            int si = ci % 4, di = ci / 4;
            float r_coeff = (si==0 ? sr : si==1 ? sg : si==2 ? sb : 1.f);
            float gout = (di==0 ? dr : di==1 ? dg : db);
            trilerp += v * r_coeff * gout;
        }
        #ifdef COMPUTE_COORDS_GRAD
        gx_grad += dwdx[corner] * (W-1) * trilerp;
        gy_grad += dwdy[corner] * (H-1) * trilerp;
        #endif
        gz_grad += dwdz[corner] * (L-1) * trilerp;
    }
    // save gradient, with discontinuty masking
    #ifdef COMPUTE_COORDS_GRAD
    v_coords[2*g_off+0] = gx_grad * (float)(x0 != x && x1 != x);
    v_coords[2*g_off+1] = gy_grad * (float)(y0 != y && y1 != y);
    #endif
    gz_grad *= (float)(z0 != z && z1 != z);
    v_rgb[3*g_off+0] = vr + kC2G_r * gz_grad;
    v_rgb[3*g_off+1] = vg + kC2G_g * gz_grad;;
    v_rgb[3*g_off+2] = vb + kC2G_b * gz_grad;;
}
