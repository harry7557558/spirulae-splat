#include "BilagridConfig.cuh"

#ifndef BILAGRID_DEPTH_BWD_V1_SHMEM_CONSTANTS_DEFINED
#define BILAGRID_DEPTH_BWD_V1_SHMEM_CONSTANTS_DEFINED
constexpr int kDepthShmemFootX  = (int)kBilagridBwdV1BlockX + 2;
constexpr int kDepthShmemFootY  = (int)kBilagridBwdV1BlockY + 2;
constexpr int kDepthShmemFootZ  = 3;
constexpr int kDepthShmemChans  = 2;
constexpr int kDepthShmemCells  = kDepthShmemFootX * kDepthShmemFootY * kDepthShmemFootZ;
constexpr int kDepthBlockSize   = (int)kBilagridBwdV1BlockX
                                * (int)kBilagridBwdV1BlockY
                                * (int)kBilagridBwdV1BlockZ;
#endif

#ifdef PATCHED
__global__ void bilagrid_depth_patched_sample_backward_v1_kernel_bilagrid(
#else
__global__ void bilagrid_depth_uniform_sample_backward_v1_kernel_bilagrid(
#endif
    BilagridReader bilagrid,  // [N,L,H,W,2]
#ifdef PATCHED
    const float* __restrict__ depth,     // [N,m,h,w,1]
#else
    const float* __restrict__ depth,     // [N,h,w,1]
#endif
    const float* __restrict__ scalars,   // [N]
#ifdef PATCHED
    const float* __restrict__ v_output,  // [N,m,h,w,1]
#else
    const float* __restrict__ v_output,  // [N,h,w,1]
#endif
    float* __restrict__ v_bilagrid,      // [N,L,H,W,2]
    int N, int L, int H, int W,
#ifdef PATCHED
    int m,
#endif
    int h, int w,
#ifdef PATCHED
    int h0, int w0,
    const int* __restrict__ offsets,     // [N,m,2]
#endif
    int mult_x, int mult_y
#ifdef PATCHED
    , int m_batch_stride
#endif
#ifndef PATCHED
    , const int* __restrict__ grid_indices  // [N], or nullptr -> identity
#endif
) {
    // 1. Block-uniform decode.
#ifdef PATCHED
    const int block_z_total = N * m_batch_stride * L;
#else
    const int block_z_total = N * L;
#endif
    if ((int)blockIdx.z >= block_z_total) return;

    int idx_block = (int)blockIdx.z;
    const int zi = idx_block % L; idx_block /= L;
#ifdef PATCHED
    const int m_batch_i = idx_block % m_batch_stride; idx_block /= m_batch_stride;
#endif
    const int ni = idx_block;

#ifndef PATCHED
    const int g_id = (grid_indices != nullptr) ? grid_indices[ni] : ni;
    const float scalar = scalars[g_id];
#else
    const int g_id = ni;
    const float scalar = scalars[0];
#endif
    const int corner_base = g_id * L * H * W * kDepthShmemChans;

    // 2. Cooperative preload.
    const int base_xi = (int)((blockIdx.x * kBilagridBwdV1BlockX)) / mult_x - 1;
    const int base_yi = (int)((blockIdx.y * kBilagridBwdV1BlockY)) / mult_y - 1;
    const int base_zi = zi - 1;

    __shared__ float shmem_corners
        [kDepthShmemFootZ][kDepthShmemFootY][kDepthShmemFootX][kDepthShmemChans];

    const int linear_tid =
        (threadIdx.z * (int)kBilagridBwdV1BlockY + threadIdx.y)
        * (int)kBilagridBwdV1BlockX + threadIdx.x;

    #pragma unroll
    for (int local = linear_tid; local < kDepthShmemCells; local += kDepthBlockSize) {
        const int lx = local % kDepthShmemFootX;
        const int ly = (local / kDepthShmemFootX) % kDepthShmemFootY;
        const int lz = local / (kDepthShmemFootX * kDepthShmemFootY);
        const int gx = min(max(base_xi + lx, 0), W - 1);
        const int gy = min(max(base_yi + ly, 0), H - 1);
        const int gz = min(max(base_zi + lz, 0), L - 1);
        const int g_off = corner_base + ((gz * H + gy) * W + gx) * kDepthShmemChans;
        #pragma unroll
        for (int ci = 0; ci < kDepthShmemChans; ci++) {
            shmem_corners[lz][ly][lx][ci] = bilagrid[g_off + ci];
        }
    }
    __syncthreads();

    // 3. Per-thread cell coords.
#ifdef PATCHED
    int x_base = 0, y_base = 0;
    int x_idx = blockIdx.x * kBilagridBwdV1BlockX + threadIdx.x;
    int y_idx = blockIdx.y * kBilagridBwdV1BlockY + threadIdx.y;
    int xi = x_idx / mult_x + x_base, xf = x_idx % mult_x;
    int yi = y_idx / mult_y + y_base, yf = y_idx % mult_y;
    bool inside = (xi >= 0 && xi < W && yi >= 0 && yi < H);
#else
    int x_idx = blockIdx.x * kBilagridBwdV1BlockX + threadIdx.x;
    int y_idx = blockIdx.y * kBilagridBwdV1BlockY + threadIdx.y;
    int xi = x_idx / mult_x, xf = x_idx % mult_x;
    int yi = y_idx / mult_y, yf = y_idx % mult_y;
    bool inside = (xi < W && yi < H);
#endif
    if (!inside && (
        mult_x * mult_y == 1 ||
        (mult_x % kBilagridBwdV1BlockX != 0 || mult_y % kBilagridBwdV1BlockY != 0)
    )) return;

#ifndef PATCHED
    float sw = float(w - 1) / float(W - 1);
    int block_wi0 = max((int)ceil((xi - 1) * sw), 0);
    int block_wi1 = min((int)floor((xi + 1) * sw), w - 1) + 1;
    float sh = float(h - 1) / float(H - 1);
    int block_hi0 = max((int)ceil((yi - 1) * sh), 0);
    int block_hi1 = min((int)floor((yi + 1) * sh), h - 1) + 1;
    int x_step = (block_wi1 - block_wi0 + mult_x - 1) / mult_x;
    int y_step = (block_hi1 - block_hi0 + mult_y - 1) / mult_y;
    int wi0 = block_wi0 + xf * x_step;
    int hi0 = block_hi0 + yf * y_step;
    int wi1 = min(block_wi0 + (xf + 1) * x_step, block_wi1);
    int hi1 = min(block_hi0 + (yf + 1) * y_step, block_hi1);
#endif

    float accum[kDepthShmemChans] = {0.0f, 0.0f};

    // 4. Pixel loop.
    if (inside)
#ifdef PATCHED
    for (int mi = m_batch_i; mi < m; mi += m_batch_stride)
#endif
    {
    #ifdef PATCHED
        int o_off = (ni * m + mi) * 2;
        int2 offset = {offsets[o_off + 0], offsets[o_off + 1]};
        float sw = float(w0 - 1) / float(W - 1);
        int block_wi0 = max((int)ceil((xi - 1) * sw), offset.x);
        int block_wi1 = min((int)floor((xi + 1) * sw), min(offset.x + w, w0) - 1) + 1;
        float sh = float(h0 - 1) / float(H - 1);
        int block_hi0 = max((int)ceil((yi - 1) * sh), offset.y);
        int block_hi1 = min((int)floor((yi + 1) * sh), min(offset.y + h, h0) - 1) + 1;
        int x_step = (block_wi1 - block_wi0 + mult_x - 1) / mult_x;
        int y_step = (block_hi1 - block_hi0 + mult_y - 1) / mult_y;
        if (!(block_wi1 > block_wi0 && block_hi1 > block_hi0) && (
            mult_x * mult_y == 1 ||
            (mult_x % kBilagridBwdV1BlockX != 0 || mult_y % kBilagridBwdV1BlockY != 0)
        )) continue;
        int wi0 = block_wi0 + xf * x_step;
        int hi0 = block_hi0 + yf * y_step;
        int wi1 = min(block_wi0 + (xf + 1) * x_step, block_wi1);
        int hi1 = min(block_hi0 + (yf + 1) * y_step, block_hi1);
    #endif

        for (int hi = hi0; hi < hi1; hi++)
        for (int wi = wi0; wi < wi1; wi++) {
        #ifdef PATCHED
            int g_off = (((ni * m + mi) * h + (hi - offset.y)) * w + (wi - offset.x));
        #else
            int g_off = ((ni * h + hi) * w + wi);
        #endif
            float sr = depth[g_off] * scalar;

        #ifdef PATCHED
            float x = (float)wi / (float)(w0 - 1) * (float)(W - 1);
            float y = (float)hi / (float)(h0 - 1) * (float)(H - 1);
        #else
            float x = (float)wi / (float)(w - 1) * (float)(W - 1);
            float y = (float)hi / (float)(h - 1) * (float)(H - 1);
        #endif
            float z = sr / (sr + 1.0f);
            z = min(max(z, 0.0f), 1.0f) * (float)(L - 1);

            int x0 = floorf(x), y0 = floorf(y), z0 = floorf(z);
            int x1 = min(x0 + 1, W - 1);
            int y1 = min(y0 + 1, H - 1);
            int z1 = min(z0 + 1, L - 1);
            if (zi != z0 && zi != z1) continue;

            float fx = x - x0, fy = y - y0, fz = z - z0;
            float accum_t = 0.0f;
            if (zi == z0) {
                if (xi == x0 && yi == y0) accum_t += (1 - fx) * (1 - fy) * (1 - fz);
                if (xi == x1 && yi == y0) accum_t += fx * (1 - fy) * (1 - fz);
                if (xi == x0 && yi == y1) accum_t += (1 - fx) * fy * (1 - fz);
                if (xi == x1 && yi == y1) accum_t += fx * fy * (1 - fz);
            }
            if (zi == z1) {
                if (xi == x0 && yi == y0) accum_t += (1 - fx) * (1 - fy) * fz;
                if (xi == x1 && yi == y0) accum_t += fx * (1 - fy) * fz;
                if (xi == x0 && yi == y1) accum_t += (1 - fx) * fy * fz;
                if (xi == x1 && yi == y1) accum_t += fx * fy * fz;
            }

            // Clamp defensively: see BilagridPpispUniformSampleBwdV1_kernel.cuh
            // for the fast-math integer-boundary rationale.
            const int sx0 = max(x0 - base_xi, 0);
            const int sx1 = min(x1 - base_xi, kDepthShmemFootX - 1);
            const int sy0 = max(y0 - base_yi, 0);
            const int sy1 = min(y1 - base_yi, kDepthShmemFootY - 1);
            const int sz0 = max(z0 - base_zi, 0);
            const int sz1 = min(z1 - base_zi, kDepthShmemFootZ - 1);

            float affine[kDepthShmemChans];
            #pragma unroll
            for (int ci = 0; ci < kDepthShmemChans; ci++) {
                float v000 = shmem_corners[sz0][sy0][sx0][ci];
                float v001 = shmem_corners[sz0][sy0][sx1][ci];
                float v010 = shmem_corners[sz0][sy1][sx0][ci];
                float v011 = shmem_corners[sz0][sy1][sx1][ci];
                float v100 = shmem_corners[sz1][sy0][sx0][ci];
                float v101 = shmem_corners[sz1][sy0][sx1][ci];
                float v110 = shmem_corners[sz1][sy1][sx0][ci];
                float v111 = shmem_corners[sz1][sy1][sx1][ci];

                float c00 = v000 * (1.0f - fx) + v001 * fx;
                float c01 = v010 * (1.0f - fx) + v011 * fx;
                float c10 = v100 * (1.0f - fx) + v101 * fx;
                float c11 = v110 * (1.0f - fx) + v111 * fx;
                float c0  = c00 * (1.0f - fy) + c01 * fy;
                float c1  = c10 * (1.0f - fy) + c11 * fy;
                float val = c0 * (1.0f - fz) + c1 * fz;

                affine[ci] = (ci == 0 ? __expf(val) : val);
            }
            sr = __logf(sr);
            float output = __expf(affine[0] * sr + affine[1]);

            float dr = v_output[g_off] * output;

            #pragma unroll
            for (int ci = 0; ci < kDepthShmemChans; ci++) {
                float r_coeff = (ci == 0 ? sr : 1.0f);
                float gout = dr;
                float grad_weight = r_coeff * gout;
                if (ci == 0) grad_weight *= affine[ci];
                accum[ci] += accum_t * grad_weight;
            }
        }
    }

    // 5. Writeback.
    int out_idx_start = (((ni * L + zi) * H + yi) * W + xi) * kDepthShmemChans;

    if (mult_x * mult_y == 1) {
        #pragma unroll
        for (int ci = 0; ci < kDepthShmemChans; ci++) {
            int out_idx = out_idx_start + ci;
            if (isfinite(accum[ci]) && accum[ci] != 0.0f)
                atomicAdd(v_bilagrid + out_idx, accum[ci]);
        }
        return;
    }

    if (mult_x % kBilagridBwdV1BlockX != 0 || mult_y % kBilagridBwdV1BlockY != 0) {
        #pragma unroll
        for (int ci = 0; ci < kDepthShmemChans; ci++) {
            int out_idx = out_idx_start + ci;
            if (isfinite(accum[ci]) && accum[ci] != 0.0f)
                atomicAdd(v_bilagrid + out_idx, accum[ci]);
        }
        return;
    }

    #pragma unroll
    for (int ci = 0; ci < kDepthShmemChans; ci++) {
        bilagrid_bwd_v1_block_atomic_add(v_bilagrid + out_idx_start + ci, accum[ci]);
        __syncthreads();
    }
}


#ifdef PATCHED
__global__ void bilagrid_depth_patched_sample_backward_v1_kernel_depth(
#else
__global__ void bilagrid_depth_uniform_sample_backward_v1_kernel_depth(
#endif
    BilagridReader bilagrid,  // [N,L,H,W,2]
#ifdef PATCHED
    const float* __restrict__ depth,  // [N,m,h,w,1]
#else
    const float* __restrict__ depth,  // [N,h,w,1]
#endif
    const float* __restrict__ scalars,  // [N]
#ifdef PATCHED
    const float* __restrict__ v_output,  // [N,m,h,w,1]
    float* __restrict__ v_depth,  // [N,m,h,w,1]
#else
    const float* __restrict__ v_output,  // [N,h,w,1]
    float* __restrict__ v_depth,  // [N,h,w,1]
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
    int idx = blockIdx.x * kBilagridBwdV1RgbThreads + threadIdx.x;
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

#ifndef PATCHED
    const float scalar = scalars[g_id];
#else
    const float scalar = scalars[0];
#endif

    // input and output colors
#ifdef PATCHED
    int g_off = (((ni * m + mi) * h + hi) * w + wi);
#else
    int g_off = ((ni * h + hi) * w + wi);
#endif
    float sr = depth[g_off] * scalar;
    float dr = v_output[g_off];
    float vr = 0.0;

    // grid coords
#ifdef PATCHED
    offsets += (ni * m + mi) * 2;
    float x = (float)(wi + offsets[0]) / (float)(w0-1) * (float)(W-1);
    float y = (float)(hi + offsets[1]) / (float)(h0-1) * (float)(H-1);
#else
    float x = (float)wi / (float)(w-1) * (float)(W-1);
    float y = (float)hi / (float)(h-1) * (float)(H-1);
#endif
    // Clamp gz to [0,1] -- matches forward and bilagrid-grad bwd branch.
    float gz_raw = sr / (sr + 1.0f);
    float gz = fminf(fmaxf(gz_raw, 0.0f), 1.0f);
    bool  gz_in_range = (gz_raw >= 0.0f && gz_raw <= 1.0f);
    float z = gz * (L-1);
    int x0 = floorf(x), y0 = floorf(y), z0 = floorf(z);
    int x1 = min(x0+1, W-1);
    int y1 = min(y0+1, H-1);
    int z1 = min(z0+1, L-1);

    // fractional parts
    float fx = x - x0, fy = y - y0, fz = z - z0;
    float w000 = (1-fx)*(1-fy)*(1-fz);
    float w001 = fx*(1-fy)*(1-fz);
    float w010 = (1-fx)*fy*(1-fz);
    float w011 = fx*fy*(1-fz);
    float w100 = (1-fx)*(1-fy)*fz;
    float w101 = fx*(1-fy)*fz;
    float w110 = (1-fx)*fy*fz;
    float w111 = fx*fy*fz;

    // Channel-last corner offsets: 2 channels of each corner are contiguous.
    int grid_base = g_id * L * H * W * 2;
    int off000 = grid_base + ((z0*H + y0)*W + x0) * 2;
    int off001 = grid_base + ((z0*H + y0)*W + x1) * 2;
    int off010 = grid_base + ((z0*H + y1)*W + x0) * 2;
    int off011 = grid_base + ((z0*H + y1)*W + x1) * 2;
    int off100 = grid_base + ((z1*H + y0)*W + x0) * 2;
    int off101 = grid_base + ((z1*H + y0)*W + x1) * 2;
    int off110 = grid_base + ((z1*H + y1)*W + x0) * 2;
    int off111 = grid_base + ((z1*H + y1)*W + x1) * 2;

    // accumulate bilagrid gradient over 2 channels
    float affine[2];
    #pragma unroll
    for (int ci = 0; ci < 2; ci++) {
        float val =
            bilagrid[off000 + ci] * w000 +
            bilagrid[off001 + ci] * w001 +
            bilagrid[off010 + ci] * w010 +
            bilagrid[off011 + ci] * w011 +
            bilagrid[off100 + ci] * w100 +
            bilagrid[off101 + ci] * w101 +
            bilagrid[off110 + ci] * w110 +
            bilagrid[off111 + ci] * w111;
        if (ci == 0)
            val = __expf(val);
        affine[ci] = val;
    }
    float log_sr = __logf(sr);
    dr *= __expf(affine[0] * log_sr + affine[1]);
    vr += affine[0] * dr / sr;

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

    // accumulate gradient into coords (chain through bilagrid values and depth)
    float gz_grad = 0.f;
    #pragma unroll
    for (int corner = 0; corner < 8; ++corner) {
        int base = corner_offs[corner];
        float trilerp = 0.f;
        // gather the corresponding bilagrid value for each of the 2 channels
        #pragma unroll
        for (int ci = 0; ci < 2; ++ci) {
            float v = bilagrid[base + ci];
            float r_coeff = (ci==0 ? log_sr * affine[ci] : 1.f);
            float gout = dr;
            trilerp += v * r_coeff * gout;
        }
        gz_grad += dwdz[corner] * (L-1) * trilerp;
    }
    // Zero gz_grad outside the [0,1] clamp range -- the clamp's vjp.
    if (!gz_in_range) gz_grad = 0.0f;
    vr += gz_grad / ((sr+1.0f) * (sr+1.0f));
    vr *= scalar;
    v_depth[g_off] = isfinite(vr) ? vr : 0.0f;
}

