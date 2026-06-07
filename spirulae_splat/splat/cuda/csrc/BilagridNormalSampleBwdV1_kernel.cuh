#include "BilagridConfig.cuh"

#include "BilagridAxisAngle.cuh"

#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>

namespace cg = cooperative_groups;

#ifdef PATCHED
__global__ void bilagrid_normal_patched_sample_backward_v1_kernel_bilagrid(
#else
__global__ void bilagrid_normal_uniform_sample_backward_v1_kernel_bilagrid(
#endif
    BilagridReader bilagrid,  // [N,L,H,W,3]
#ifdef PATCHED
    const float* __restrict__ normal_in,  // [N,m,h,w,3]
    const float* __restrict__ v_normal_out,  // [N,m,h,w,3]
#else
    const float* __restrict__ normal_in,  // [N,h,w,3]
    const float* __restrict__ v_normal_out,  // [N,h,w,3]
#endif
    float* __restrict__ v_bilagrid,  // [N,L,H,W,3]
    int N, int L, int H, int W,
#ifdef PATCHED
    int m,
#endif
    int h, int w,
#ifdef PATCHED
    int h0, int w0,
    const int* __restrict__ offsets,  // [N,m,2]
#endif
    int mult_x, int mult_y
#ifdef PATCHED
    , int m_batch_stride
#endif
#ifndef PATCHED
    , const int* __restrict__ grid_indices  // [N], or nullptr -> identity
#endif
) {
#ifdef PATCHED
    int idx = blockIdx.z * blockDim.z + threadIdx.z;
    bool inside = (idx < (N*m*L));
    int zi = idx % L; idx /= L;
    int m_batch_i = idx % m_batch_stride; idx /= m_batch_stride;
    int ni = idx;

    // offsets += (ni * m + mi) * 2;
    // int2 offset = {offsets[0], offsets[1]};
    // int x_base = max((offset.x * W) / w0 - 1, 0);
    // int y_base = max((offset.y * H) / h0 - 1, 0);
    int x_base = 0, y_base = 0;

    int x_idx = blockIdx.x * blockDim.x + threadIdx.x;
    int y_idx = blockIdx.y * blockDim.y + threadIdx.y;
    int xi = x_idx / mult_x + x_base, xf = x_idx % mult_x;
    int yi = y_idx / mult_y + y_base, yf = y_idx % mult_y;
    // printf("x_idx=%d y_idx=%d  xi=%d xf=%d yi=%d yf=%d\n", x_idx, y_idx, xi, xf, yi, yf);

    inside &= (xi >= 0 && xi < W && yi >= 0 && yi < H);
#else
    int x_idx = blockIdx.x * blockDim.x + threadIdx.x;
    int y_idx = blockIdx.y * blockDim.y + threadIdx.y;
    int idx = blockIdx.z * blockDim.z + threadIdx.z;
    int xi = x_idx / mult_x, xf = x_idx % mult_x;
    int yi = y_idx / mult_y, yf = y_idx % mult_y;
    bool inside = (xi < W && yi < H && idx < (N*L));
    int zi = idx % L; idx /= L;
    int ni = idx;
#endif
    if (!inside && (
        mult_x*mult_y == 1 ||
        (mult_x % blockDim.x != 0 || mult_y % blockDim.y != 0)
    )) return;

#ifndef PATCHED
    bool use_indirect = (grid_indices != nullptr);
    int g_id = use_indirect ? grid_indices[ni] : ni;
#else
    int g_id = ni;
#endif

    // Loop bounds
#ifndef PATCHED
    float sw = float(w-1)/float(W-1);
    int block_wi0 = max((int)ceil((xi-1)*sw), 0);  // same for all threads in the block
    int block_wi1 = min((int)floor((xi+1)*sw), w-1) + 1;
    float sh = float(h-1)/float(H-1);
    int block_hi0 = max((int)ceil((yi-1)*sh), 0);
    int block_hi1 = min((int)floor((yi+1)*sh), h-1) + 1;
    int x_step = (block_wi1-block_wi0+mult_x-1)/mult_x;
    int y_step = (block_hi1-block_hi0+mult_y-1)/mult_y;

    int wi0 = block_wi0+xf*x_step;
    int hi0 = block_hi0+yf*y_step;
    // int wi1 = min(block_wi0+(xf+1)*x_step, w);
    // int hi1 = min(block_hi0+(yf+1)*y_step, h);
    int wi1 = min(block_wi0+(xf+1)*x_step, block_wi1);
    int hi1 = min(block_hi0+(yf+1)*y_step, block_hi1);
#endif

    // Result for each affine mat channel
    float accum[3] = {0.0f, 0.0f, 0.0f};

    // Loop over all samples for this batch
    if (inside)
#ifdef PATCHED
    for (int mi = m_batch_i; mi < m; mi += m_batch_stride)
#endif
    {
    #ifdef PATCHED
        int o_off = (ni * m + mi) * 2;
        int2 offset = {offsets[o_off+0], offsets[o_off+1]};
        float sw = float(w0-1)/float(W-1);
        int block_wi0 = max((int)ceil((xi-1)*sw), offset.x);  // same for all threads in the block
        int block_wi1 = min((int)floor((xi+1)*sw), min(offset.x+w,w0)-1) + 1;
        float sh = float(h0-1)/float(H-1);
        int block_hi0 = max((int)ceil((yi-1)*sh), offset.y);
        int block_hi1 = min((int)floor((yi+1)*sh), min(offset.y+h,h0)-1) + 1;
        int x_step = (block_wi1-block_wi0+mult_x-1)/mult_x;
        int y_step = (block_hi1-block_hi0+mult_y-1)/mult_y;
        // block_wi1 = max(block_wi1, block_wi0);
        // block_hi1 = max(block_hi1, block_hi0);
        if (!(block_wi1 > block_wi0 && block_hi1 > block_hi0) && (
            mult_x*mult_y == 1 ||
            (mult_x % blockDim.x != 0 || mult_y % blockDim.y != 0)
        )) continue;
        int wi0 = block_wi0+xf*x_step;
        int hi0 = block_hi0+yf*y_step;
        int wi1 = min(block_wi0+(xf+1)*x_step, block_wi1);
        int hi1 = min(block_hi0+(yf+1)*y_step, block_hi1);
    #endif

        for (int hi = hi0; hi < hi1; hi++)
        for (int wi = wi0; wi < wi1; wi++) {

        #ifdef PATCHED
            int g_off = (((ni*m + mi)*h + (hi-offset.y))*w + (wi-offset.x))*3;
        #else
            int g_off = ((ni*h + hi)*w + wi)*3;
        #endif
            float sr = normal_in[g_off+0];
            float sg = normal_in[g_off+1];
            float sb = normal_in[g_off+2];
            float inv_norm = rsqrtf(sr*sr + sg*sg + sb*sb + 1e-20f);

        #ifdef PATCHED
            float x = (float)wi / (float)(w0-1) * (float)(W-1);
            float y = (float)hi / (float)(h0-1) * (float)(H-1);
        #else
            float x = (float)wi / (float)(w-1) * (float)(W-1);
            float y = (float)hi / (float)(h-1) * (float)(H-1);
        #endif
            // float z = acosf(fminf(fmaxf(sb*inv_norm, -1.0f), 1.0f)) * (1.0f / (float)M_PI);
            float z = 0.5f + 0.5f * sb*inv_norm;
            z = min(max(z, 0.0f), 1.0f) * (float)(L-1);

            int x0 = floorf(x), y0 = floorf(y), z0 = floorf(z);
            int x1 = min(x0+1, W-1);
            int y1 = min(y0+1, H-1);
            int z1 = min(z0+1, L-1);
            if (zi != z0 && zi != z1) continue;

            float fx = x-x0, fy = y-y0, fz = z-z0;

            // Channel-last corner offsets shared across the C-channel loop.
            int corner_base = g_id * L * H * W * 3;
            int off000 = corner_base + ((z0*H + y0)*W + x0) * 3;
            int off001 = corner_base + ((z0*H + y0)*W + x1) * 3;
            int off010 = corner_base + ((z0*H + y1)*W + x0) * 3;
            int off011 = corner_base + ((z0*H + y1)*W + x1) * 3;
            int off100 = corner_base + ((z1*H + y0)*W + x0) * 3;
            int off101 = corner_base + ((z1*H + y0)*W + x1) * 3;
            int off110 = corner_base + ((z1*H + y1)*W + x0) * 3;
            int off111 = corner_base + ((z1*H + y1)*W + x1) * 3;

            // load normal params
            float3 axis_angle;
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

                (ci == 0 ? axis_angle.x : ci == 1 ? axis_angle.y : axis_angle.z) = val;
            }

            float dr = v_normal_out[g_off+0];
            float dg = v_normal_out[g_off+1];
            float db = v_normal_out[g_off+2];

            // apply normal
            float3 normal = {sr*inv_norm, sg*inv_norm, sb*inv_norm};
            float3 grad_axis_angle = {0.0f, 0.0f, 0.0f};
            float3 grad_normal = {0.0f, 0.0f, 0.0f};
            axis_angle_rotate_bwd(axis_angle, normal, {dr, dg, db}, grad_axis_angle, grad_normal);

            // accumulate gradient
            float accum_t = 0.0;
            if (zi == z0) {
                if (xi == x0 && yi == y0) accum_t += (1-fx)*(1-fy)*(1-fz);
                if (xi == x1 && yi == y0) accum_t += fx*(1-fy)*(1-fz);
                if (xi == x0 && yi == y1) accum_t += (1-fx)*fy*(1-fz);
                if (xi == x1 && yi == y1) accum_t += fx*fy*(1-fz);
            }
            if (zi == z1) {
                if (xi == x0 && yi == y0) accum_t += (1-fx)*(1-fy)*fz;
                if (xi == x1 && yi == y0) accum_t += fx*(1-fy)*fz;
                if (xi == x0 && yi == y1) accum_t += (1-fx)*fy*fz;
                if (xi == x1 && yi == y1) accum_t += fx*fy*fz;
            }

            #pragma unroll
            for (int ci = 0; ci < 3; ci++) {
                float grad_weight = (ci == 0 ? grad_axis_angle.x :
                    ci == 1 ? grad_axis_angle.y : grad_axis_angle.z);
                accum[ci] += accum_t * grad_weight;
            }

        }
    }

    // Output indexed by `ni` (sparse batch slot); reads use `g_id`.
    // Channel-last write: 3 channels for the cell are contiguous in memory.
    int out_idx_start = (((ni*L + zi)*H + yi)*W + xi) * 3;

    if (mult_x*mult_y == 1) {
        #pragma unroll
        for (int ci = 0; ci < 3; ci++) {
            int out_idx = out_idx_start + ci;
            if (isfinite(accum[ci]) && accum[ci] != 0.0f)
                atomicAdd(v_bilagrid + out_idx, accum[ci]);
        }
        return;
    }

    // out_idx can be different for each thread, fall back to global atomicAdd
    if (mult_x % blockDim.x != 0 || mult_y % blockDim.y != 0) {
        #pragma unroll
        for (int ci = 0; ci < 3; ci++) {
            int out_idx = out_idx_start + ci;
            if (isfinite(accum[ci]) && accum[ci] != 0.0f)
                atomicAdd(v_bilagrid + out_idx, accum[ci]);
        }
        return;
    }

    // fast atomicAdd

    __shared__ float sharedData[64];

    int blockSize = blockDim.x * blockDim.y * blockDim.z;
    int tid = (threadIdx.z * blockDim.y + threadIdx.y) * blockDim.x + threadIdx.x;

    #pragma unroll
    for (int ci = 0; ci < 3; ci++) {
        int out_idx = out_idx_start + ci;

        sharedData[tid] = isfinite(accum[ci]) ? accum[ci] : 0.0f;
        __syncthreads();

        for (int s = blockSize / 2; s > 0; s >>= 1) {
            if (tid < s)
                sharedData[tid] += sharedData[tid + s];
            __syncthreads();
        }

        if (tid == 0)
            atomicAdd(v_bilagrid + out_idx, sharedData[0]);
    }

}


#ifdef PATCHED
__global__ void bilagrid_normal_patched_sample_backward_v1_kernel_normal(
#else
__global__ void bilagrid_normal_uniform_sample_backward_v1_kernel_normal(
#endif
    BilagridReader bilagrid,  // [N,L,H,W,3]
#ifdef PATCHED
    const float* __restrict__ normal_in,  // [N,m,h,w,3]
    const float* __restrict__ v_normal_out,  // [N,m,h,w,3]
    float* __restrict__ v_normal_in,  // [N,m,h,w,3]
#else
    const float* __restrict__ normal_in,  // [N,h,w,3]
    const float* __restrict__ v_normal_out,  // [N,h,w,3]
    float* __restrict__ v_normal_in,  // [N,h,w,3]
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

    // input and output colors
#ifdef PATCHED
    int g_off = (((ni * m + mi) * h + hi) * w + wi) * 3;
#else
    int g_off = ((ni * h + hi) * w + wi) * 3;
#endif
    float sr = normal_in[g_off+0];
    float sg = normal_in[g_off+1];
    float sb = normal_in[g_off+2];
    float inv_norm = rsqrtf(sr*sr + sg*sg + sb*sb + 1e-20f);
    sr *= inv_norm, sg *= inv_norm, sb *= inv_norm;
    float dr = v_normal_out[g_off+0];
    float dg = v_normal_out[g_off+1];
    float db = v_normal_out[g_off+2];

    // grid coords
#ifdef PATCHED
    offsets += (ni * m + mi) * 2;
    float x = (float)(wi + offsets[0]) / (float)(w0-1) * (float)(W-1);
    float y = (float)(hi + offsets[1]) / (float)(h0-1) * (float)(H-1);
#else
    float x = (float)wi / (float)(w-1) * (float)(W-1);
    float y = (float)hi / (float)(h-1) * (float)(H-1);
#endif
    // float z = (acosf(fminf(fmaxf(sb, -1.0f), 1.0f)) * (1.0f / (float)M_PI)) * (L-1);
    float z = (0.5f + 0.5f * sb) * (L-1);
    int x0 = floorf(x), y0 = floorf(y), z0 = floorf(z);
    int x1 = min(x0+1, W-1);
    int y1 = min(y0+1, H-1);
    int z1 = z0 + 1;
    z0 = min(max(z0,0), L-1); z1 = min(max(z1,0), L-1);

    float fx = x-x0, fy = y-y0, fz = z-z0;

    // Channel-last corner offsets: 3 channels of each corner are contiguous.
    int grid_base = g_id * L * H * W * 3;
    int off000 = grid_base + ((z0*H + y0)*W + x0) * 3;
    int off001 = grid_base + ((z0*H + y0)*W + x1) * 3;
    int off010 = grid_base + ((z0*H + y1)*W + x0) * 3;
    int off011 = grid_base + ((z0*H + y1)*W + x1) * 3;
    int off100 = grid_base + ((z1*H + y0)*W + x0) * 3;
    int off101 = grid_base + ((z1*H + y0)*W + x1) * 3;
    int off110 = grid_base + ((z1*H + y1)*W + x0) * 3;
    int off111 = grid_base + ((z1*H + y1)*W + x1) * 3;

    // load normal params
    float3 axis_angle;
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

        (ci == 0 ? axis_angle.x : ci == 1 ? axis_angle.y : axis_angle.z) = val;
    }

    // apply normal
    float3 normal = {sr, sg, sb};
    float3 grad_axis_angle = {0.0f, 0.0f, 0.0f};
    float3 grad_normal = {0.0f, 0.0f, 0.0f};
    axis_angle_rotate_bwd(axis_angle, normal, {dr, dg, db}, grad_axis_angle, grad_normal);

    // spatial derivatives for coords
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

    // accumulate gradient into coords (chain through bilagrid values and normal)
    float gz_grad = 0.f;
    #pragma unroll
    for (int corner = 0; corner < 8; ++corner) {
        int base = corner_offs[corner];
        float trilerp = 0.f;
        // gather the corresponding bilagrid value for each of the 3 channels
        #pragma unroll
        for (int ci = 0; ci < 3; ++ci) {
            float v = bilagrid[base + ci];
            trilerp += v * (ci == 0 ? grad_axis_angle.x :
                ci == 1 ? grad_axis_angle.y : grad_axis_angle.z);
        }
        gz_grad += dwdz[corner] * (L-1) * trilerp;
    }
    // grad_normal.z += gz_grad * -rsqrtf(fmaxf(1.0f - sb*sb, 1e-20f)) * (1.0f / (float)M_PI);
    grad_normal.z += 0.5f * gz_grad;

    grad_normal = mul3(add3(grad_normal, mul3(normal, -dot3(normal, grad_normal))), inv_norm);
    v_normal_in[g_off+0] = isfinite(grad_normal.x) ? grad_normal.x : 0.0f;
    v_normal_in[g_off+1] = isfinite(grad_normal.y) ? grad_normal.y : 0.0f;
    v_normal_in[g_off+2] = isfinite(grad_normal.z) ? grad_normal.z : 0.0f;
}

