#include "BilagridConfig.cuh"

#include "BilagridPpispMathBwd.cuh"

using namespace bilagrid_ppisp;


// template<int compute_coords_grad>
#ifdef PACKED
#ifdef COMPUTE_COORDS_GRAD
__global__ void bilagrid_ppisp_packed_sample_backward_kernel_cg(
#else
__global__ void bilagrid_ppisp_packed_sample_backward_kernel(
#endif  // COMPUTE_COORDS_GRAD
    BilagridReader bilagrid,  // [N,L,H,W,9]
    const int64_t* __restrict__ image_indices,  // [nnz]
    const float* __restrict__ coords,  // [nnz,2]
    const float* __restrict__ rgb_in,  // [nnz,3]
    const float* __restrict__ v_rgb_out,  // [nnz,3]
    float* __restrict__ v_bilagrid,  // [N,L,H,W,9]
    #ifdef COMPUTE_COORDS_GRAD
    float* __restrict__ v_coords,  // [nnz,2]
    #endif
    float* __restrict__ v_rgb_in,  // [nnz,3]
    int N, int L, int H, int W,
    int nnz
#else  // PACKED
#ifdef COMPUTE_COORDS_GRAD
__global__ void bilagrid_ppisp_sample_backward_kernel_cg(
#else
__global__ void bilagrid_ppisp_sample_backward_kernel(
#endif  // COMPUTE_COORDS_GRAD
    BilagridReader bilagrid,  // [N,L,H,W,9]
    const float* __restrict__ coords,  // [N,m,h,w,2]
    const float* __restrict__ rgb_in,  // [N,m,h,w,3]
    const float* __restrict__ v_rgb_out,  // [N,m,h,w,3]
    float* __restrict__ v_bilagrid,  // [N,L,H,W,9]
    #ifdef COMPUTE_COORDS_GRAD
    float* __restrict__ v_coords,  // [N,m,h,w,2]
    #endif
    float* __restrict__ v_rgb_in,  // [N,m,h,w,3]
    int N, int L, int H, int W,
    int m, int h, int w
#endif  // PACKED
) {
#ifdef PACKED
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    bool inside = (idx < nnz);
    if (!inside) return;
    int ni = image_indices[idx];

    // grid coords
    int g_off = idx;
#else  // PACKED
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
#endif  // PACKED
    float sr = rgb_in[3*g_off+0], sg = rgb_in[3*g_off+1], sb = rgb_in[3*g_off+2];
    float gx = coords[2*g_off+0];
    float gy = coords[2*g_off+1];
    // Clamp gz to [0,1] -- matches forward and V1 bilagrid-grad bwd branch.
    float gz_raw = kC2G_r * sr + kC2G_g * sg + kC2G_b * sb;
    float gz = fminf(fmaxf(gz_raw, 0.0f), 1.0f);
    bool  gz_in_range = (gz_raw >= 0.0f && gz_raw <= 1.0f);
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

    // Channel-last corner offsets shared across the C-channel loop.
    int corner_base = ni * L * H * W * 9;
    int off000 = corner_base + ((z0*H + y0)*W + x0) * 9;
    int off001 = corner_base + ((z0*H + y0)*W + x1) * 9;
    int off010 = corner_base + ((z0*H + y1)*W + x0) * 9;
    int off011 = corner_base + ((z0*H + y1)*W + x1) * 9;
    int off100 = corner_base + ((z1*H + y0)*W + x0) * 9;
    int off101 = corner_base + ((z1*H + y0)*W + x1) * 9;
    int off110 = corner_base + ((z1*H + y1)*W + x0) * 9;
    int off111 = corner_base + ((z1*H + y1)*W + x1) * 9;

    // load PPISP params
    float exposure_param;
    ColorPPISPParams color_params;
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

    // read rgb coeffs and upstream gradient
    float dr = v_rgb_out[3*g_off+0];
    float dg = v_rgb_out[3*g_off+1];
    float db = v_rgb_out[3*g_off+2];

    // apply PPISP
    float3 rgb = {sr, sg, sb};
    float3 grad_rgb = {dr, dg, db};
    float grad_exposure_param = 0.0f;
    ColorPPISPParams grad_color_params = {{0.0f, 0.0f}, {0.0f, 0.0f}, {0.0f, 0.0f}, {0.0f, 0.0f}};
    float3 rgb_after_exp;
    apply_exposure(rgb, exposure_param, rgb_after_exp);
    apply_color_correction_ppisp_bwd(rgb_after_exp, &color_params, grad_rgb, grad_rgb, &grad_color_params);
    apply_exposure_bwd(rgb, exposure_param, grad_rgb, grad_rgb, grad_exposure_param);

    float vr = grad_rgb.x, vg = grad_rgb.y, vb = grad_rgb.z;

    // accumulate bilagrid gradient over 9 parameters
    #pragma unroll
    for (int ci = 0; ci < 9; ++ci) {
        // weight from rgb channel
        float grad_weight = (ci == 0 ? grad_exposure_param :
            ci == 1 ? grad_color_params.b.x :
            ci == 2 ? grad_color_params.b.y :
            ci == 3 ? grad_color_params.r.x :
            ci == 4 ? grad_color_params.r.y :
            ci == 5 ? grad_color_params.g.x :
            ci == 6 ? grad_color_params.g.y :
            ci == 7 ? grad_color_params.n.x :
            grad_color_params.n.y);

        // scatter back into the eight corners
        atomicAdd(v_bilagrid + off000 + ci, f000 * grad_weight);
        atomicAdd(v_bilagrid + off001 + ci, f001 * grad_weight);
        atomicAdd(v_bilagrid + off010 + ci, f010 * grad_weight);
        atomicAdd(v_bilagrid + off011 + ci, f011 * grad_weight);
        atomicAdd(v_bilagrid + off100 + ci, f100 * grad_weight);
        atomicAdd(v_bilagrid + off101 + ci, f101 * grad_weight);
        atomicAdd(v_bilagrid + off110 + ci, f110 * grad_weight);
        atomicAdd(v_bilagrid + off111 + ci, f111 * grad_weight);

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
        // gather the corresponding bilagrid value for each of the 9 parameters
        #pragma unroll
        for (int ci = 0; ci < 9; ++ci) {
            float v = bilagrid[base + ci];
            float grad_weight = (ci == 0 ? grad_exposure_param :
                ci == 1 ? grad_color_params.b.x :
                ci == 2 ? grad_color_params.b.y :
                ci == 3 ? grad_color_params.r.x :
                ci == 4 ? grad_color_params.r.y :
                ci == 5 ? grad_color_params.g.x :
                ci == 6 ? grad_color_params.g.y :
                ci == 7 ? grad_color_params.n.x :
                grad_color_params.n.y);
            trilerp += v * grad_weight;
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
    // Zero gz_grad outside the [0,1] clamp range -- the clamp's vjp.
    if (!gz_in_range) gz_grad = 0.0f;
    v_rgb_in[3*g_off+0] = vr + kC2G_r * gz_grad;
    v_rgb_in[3*g_off+1] = vg + kC2G_g * gz_grad;;
    v_rgb_in[3*g_off+2] = vb + kC2G_b * gz_grad;;
}
