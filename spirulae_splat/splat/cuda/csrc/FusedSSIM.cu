// Based on https://github.com/MrNeRF/optimized-fused-ssim/ by Rahul Goel and MrNeRF
// MIT License: https://github.com/MrNeRF/optimized-fused-ssim/blob/main/LICENSE

// Changes:
//  - Output scalar loss only, instead of full SSIM map and reduce on Python side
//  - Hard-core float3 with channel-last memory layout
//  - Hard-core padding=same

#include "FusedSSIM.cuh"

#include <ATen/ops/empty.h>
#include <ATen/ops/empty_like.h>
#include <ATen/ops/zeros.h>

#include <cooperative_groups.h>
#include <algorithm>
#include <iostream>
#include <c10/cuda/CUDAGuard.h>

namespace cg = cooperative_groups;

__forceinline__ __device__ float3 fmaxf(float3 v, float k) {
    return {
        fmaxf(v.x, k),
        fmaxf(v.y, k),
        fmaxf(v.z, k)
    };
}
__forceinline__ __device__ float3 sqrtf(float3 v) {
    return {
        sqrtf(fmaxf(v.x, 0.0f)),
        sqrtf(fmaxf(v.y, 0.0f)),
        sqrtf(fmaxf(v.z, 0.0f))
    };
}
__forceinline__ __device__ float3 fabsf(float3 v) {
    return {
        fabsf(v.x), fabsf(v.y), fabsf(v.z)
    };
}
__forceinline__ __device__ float3 fsignf(float3 v) {
    return {
        copysignf(1.0f, v.x),
        copysignf(1.0f, v.y),
        copysignf(1.0f, v.z)
    };
}

// ------------------------------------------
// Block and Shared Memory Dimensions
// ------------------------------------------
#define BLOCK_X 16
#define BLOCK_Y 16
#define HALO    5

#define SHARED_X (BLOCK_X + 2 * HALO)
#define SHARED_Y (BLOCK_Y + 2 * HALO)

// For partial results after horizontal pass
#define CONV_X BLOCK_X
#define CONV_Y SHARED_Y

// ------------------------------------------
// Constant Memory for Gaussian Coefficients
// ------------------------------------------
__constant__ float cGauss[HALO+1] = {
    0.001028380123898387f,
    0.0075987582094967365f,
    0.036000773310661316f,
    0.10936068743467331f,
    0.21300552785396576f,
    0.26601171493530273f,
    // 0.21300552785396576f,
    // 0.10936068743467331f,
    // 0.036000773310661316f,
    // 0.0075987582094967365f,
    // 0.001028380123898387f
};

// ------------------------------------------
// Utility: Safe pixel fetch w/ zero padding
// ------------------------------------------
__device__ __forceinline__ float3 get_pix_value(
    const float3* img, 
    int b, int y, int x,
    int H, int W
) {
    if (x < 0 || x >= W || y < 0 || y >= H) {
        return make_float3(0);
    }
    return img[b * H * W + y * W + x];
}

__device__ __forceinline__ float get_pix_value(
    const float3* img, 
    int b, int y, int x, int c,
    int H, int W
) {
    if (x < 0 || x >= W || y < 0 || y >= H) {
        return 0.0f;
    }
    return ((float*)img)[(b * H * W + y * W + x) * 3 + c];
}

// ------------------------------------------
// Forward Kernel: Fused SSIM
//  - Two-pass convolution to get mu1, mu2,
//    sigma1_sq, sigma2_sq, sigma12, etc.
//  - Writes final SSIM map to ssim_loss_map
//  - Optionally writes partial derivatives
//    to dm_dmu1, dm_dsigma1_sq, dm_dsigma12
// ------------------------------------------
template<bool is_l1, bool inplace>
__global__ void ssim_forward_kernel(
    int B, int H, int W,
    const float3* __restrict__ img1,
    const float3* __restrict__ img2,
    float* __restrict__ ssim,
    float ssim_loss_map_weight,
    float* __restrict__ ssim_loss_map,
    float3* __restrict__ dm_dmu1,
    float3* __restrict__ dm_dsigma1_sq,
    float3* __restrict__ dm_dsigma12
) {
    auto block = cg::this_thread_block();
    const int bIdx   = block.group_index().z;  // batch index
    const int pix_y  = block.group_index().y * BLOCK_Y + block.thread_index().y;
    const int pix_x  = block.group_index().x * BLOCK_X + block.thread_index().x;
    const int pix_id = pix_y * W + pix_x;
    const int num_pix = H * W;

    // Shared memory for the tile (img1, img2)
    __shared__ float3 sTile[SHARED_Y][SHARED_X][2];
    // After horizontal pass, store partial sums here
    // xconv[y][x] -> (sumX, sumX^2, sumY, sumY^2, sumXY)
    __shared__ float3 xconv[CONV_Y][CONV_X][5];

    // Each block processes B x C sub-batches. We loop over channels:

    // ------------------------------------------------------------
    // 1) Load (img1, img2) tile + halo into shared memory
    // ------------------------------------------------------------
    {
        const int tileSize = SHARED_Y * SHARED_X;
        const int threads = BLOCK_X * BLOCK_Y;
        const int steps = (tileSize + threads - 1) / threads;

        const int tileStartY = block.group_index().y * BLOCK_Y;
        const int tileStartX = block.group_index().x * BLOCK_X;

        for (int s = 0; s < steps; ++s) {
            int tid = s * threads + block.thread_rank();
            if (tid < tileSize) {
                int local_y = tid / SHARED_X;
                int local_x = tid % SHARED_X;
                int gy = tileStartY + local_y - HALO;
                int gx = tileStartX + local_x - HALO;

                float3 X = get_pix_value(img1, bIdx, gy, gx, H, W);
                float3 Y = get_pix_value(img2, bIdx, gy, gx, H, W);

                sTile[local_y][local_x][0] = X;
                sTile[local_y][local_x][1] = Y;
            }
        }
    }
    block.sync();

    // ------------------------------------------------------------
    // 2) Horizontal convolution (11x1) in shared memory
    //    We'll accumulate symmetrical pairs around center.
    // ------------------------------------------------------------
    {
        int ly = threadIdx.y;
        int lx = threadIdx.x + HALO;  // skip left halo

        float3 sumX   = make_float3(0);
        float3 sumX2  = make_float3(0);
        float3 sumY   = make_float3(0);
        float3 sumY2  = make_float3(0);
        float3 sumXY  = make_float3(0);

        // #pragma unroll for those 5 pairs
        #pragma unroll
        for (int d = 1; d <= HALO; ++d) {
            float w = cGauss[HALO - d];
            float3 Xleft  = sTile[ly][lx - d][0];
            float3 Yleft  = sTile[ly][lx - d][1];
            float3 Xright = sTile[ly][lx + d][0];
            float3 Yright = sTile[ly][lx + d][1];

            sumX  += (Xleft + Xright) * w;
            sumX2 += ((Xleft * Xleft) + (Xright * Xright)) * w;
            sumY  += (Yleft + Yright) * w;
            sumY2 += ((Yleft * Yleft) + (Yright * Yright)) * w;
            sumXY += ((Xleft * Yleft) + (Xright * Yright)) * w;
        }
        // center
        {
            float3 centerX = sTile[ly][lx][0];
            float3 centerY = sTile[ly][lx][1];
            float wc = cGauss[HALO];
            sumX  += centerX * wc;
            sumX2 += (centerX * centerX) * wc;
            sumY  += centerY * wc;
            sumY2 += (centerY * centerY) * wc;
            sumXY += (centerX * centerY) * wc;
        }

        // Write out partial sums
        xconv[ly][threadIdx.x][0] = sumX;
        xconv[ly][threadIdx.x][1] = sumX2;
        xconv[ly][threadIdx.x][2] = sumY;
        xconv[ly][threadIdx.x][3] = sumY2;
        xconv[ly][threadIdx.x][4] = sumXY;

        // Possibly handle second row in same warp
        int ly2 = ly + BLOCK_Y;
        if (ly2 < CONV_Y) {
            sumX   = make_float3(0); sumX2  = make_float3(0);
            sumY   = make_float3(0); sumY2  = make_float3(0);
            sumXY  = make_float3(0);

            #pragma unroll
            for (int d = 1; d <= HALO; ++d) {
                float w = cGauss[HALO - d];
                float3 Xleft  = sTile[ly2][lx - d][0];
                float3 Yleft  = sTile[ly2][lx - d][1];
                float3 Xright = sTile[ly2][lx + d][0];
                float3 Yright = sTile[ly2][lx + d][1];

                sumX  += (Xleft + Xright) * w;
                sumX2 += ((Xleft * Xleft) + (Xright * Xright)) * w;
                sumY  += (Yleft + Yright) * w;
                sumY2 += ((Yleft * Yleft) + (Yright * Yright)) * w;
                sumXY += ((Xleft * Yleft) + (Xright * Yright)) * w;
            }
            // center
            {
                float3 cx = sTile[ly2][lx][0];
                float3 cy = sTile[ly2][lx][1];
                float wc = cGauss[HALO];
                sumX  += cx * wc;
                sumX2 += (cx * cx) * wc;
                sumY  += cy * wc;
                sumY2 += (cy * cy) * wc;
                sumXY += (cx * cy) * wc;
            }
            xconv[ly2][threadIdx.x][0] = sumX;
            xconv[ly2][threadIdx.x][1] = sumX2;
            xconv[ly2][threadIdx.x][2] = sumY;
            xconv[ly2][threadIdx.x][3] = sumY2;
            xconv[ly2][threadIdx.x][4] = sumXY;
        }
    }
    block.sync();

    // ------------------------------------------------------------
    // 3) Vertical convolution (1x11) + final SSIM
    // ------------------------------------------------------------
    {
        int ly = threadIdx.y + HALO;
        int lx = threadIdx.x;

        float3 out0 = make_float3(0), out1 = make_float3(0), out2 = make_float3(0), out3 = make_float3(0), out4 = make_float3(0);

        #pragma unroll
        for (int d = 1; d <= HALO; ++d) {
            float w = cGauss[HALO - d];
            float3* top = xconv[ly - d][lx];
            float3* bot = xconv[ly + d][lx];

            out0 += (top[0] + bot[0]) * w;
            out1 += (top[1] + bot[1]) * w;
            out2 += (top[2] + bot[2]) * w;
            out3 += (top[3] + bot[3]) * w;
            out4 += (top[4] + bot[4]) * w;
        }
        // center
        {
            float wC = cGauss[HALO];
            float3* ctr = xconv[ly][lx];
            out0 += ctr[0] * wC;
            out1 += ctr[1] * wC;
            out2 += ctr[2] * wC;
            out3 += ctr[3] * wC;
            out4 += ctr[4] * wC;
        }

        float val = 0.0;
        int global_idx = bIdx * num_pix + pix_id;
        if (pix_x < W && pix_y < H) {
            static constexpr float kC1 = 0.01f * 0.01f;
            static constexpr float kC2 = 0.03f * 0.03f;

            float3 mu1 = out0;
            float3 mu2 = out2;
            float3 mu1_sq = mu1 * mu1;
            float3 mu2_sq = mu2 * mu2;

            float3 sigma1_sq = out1 - mu1_sq;
            float3 sigma2_sq = out3 - mu2_sq;
            float3 sigma12   = out4 - mu1 * mu2;

        if (!is_l1) {
            // standard SSIM

            float3 A = mu1_sq + mu2_sq + kC1;
            float3 B = sigma1_sq + sigma2_sq + kC2;
            float3 C_ = 2.f * mu1 * mu2 + kC1;
            float3 D_ = 2.f * sigma12 + kC2;

            float3 val3 = (C_ * D_) / (A * B);
            val = (val3.x + val3.y + val3.z) * (1.0f/3.0f);

            if (dm_dmu1) {
                // partial derivatives
                float3 d_m_dmu1 = (
                    (mu2 * 2.f * D_) / (A * B)
                    - (mu2 * 2.f * C_) / (A * B)
                    - (mu1 * 2.f * C_ * D_) / (A * A * B)
                    + (mu1 * 2.f * C_ * D_) / (A * B * B)
                );
                float3 d_m_dsigma1_sq = (-C_ * D_) / (A * B * B);
                float3 d_m_dsigma12   = (2.f * C_) / (A * B);

                dm_dmu1[global_idx]       = d_m_dmu1;
                dm_dsigma1_sq[global_idx] = d_m_dsigma1_sq;
                dm_dsigma12[global_idx]   = d_m_dsigma12;
                if (ssim_loss_map != nullptr) {
                    if (inplace)
                        ssim_loss_map[global_idx] += ssim_loss_map_weight * (1.0f - val);
                    else
                        ssim_loss_map[global_idx] = 1.0f - val;
                }
            }

        } else {
            // L1 version of SSIM
            float3 A = mu1_sq + mu2_sq + kC1;
            float3 B = sigma1_sq + sigma2_sq + kC2;
            float3 C_ = 2.f * mu1 * mu2 + kC1;
            float3 D_ = 2.f * sigma12 + kC2;

            float3 val3 = (C_ * D_) / (A * B);
            val3 = 1.0f - sqrtf(1.0f - val3);
            val = (val3.x + val3.y + val3.z) * (1.0f/3.0f);

            if (dm_dmu1) {
                // partial derivatives
                float3 d_m_dmu1 = (
                    (mu2 * 2.f * D_) / (A * B)
                    - (mu2 * 2.f * C_) / (A * B)
                    - (mu1 * 2.f * C_ * D_) / (A * A * B)
                    + (mu1 * 2.f * C_ * D_) / (A * B * B)
                );
                float3 d_m_dsigma1_sq = (-C_ * D_) / (A * B * B);
                float3 d_m_dsigma12   = (2.f * C_) / (A * B);

                float3 grad_val3 = 0.5f / fmaxf(1.0f - val3, 1e-5f);

                dm_dmu1[global_idx]       = d_m_dmu1 * grad_val3;
                dm_dsigma1_sq[global_idx] = d_m_dsigma1_sq * grad_val3;
                dm_dsigma12[global_idx]   = d_m_dsigma12 * grad_val3;
                if (ssim_loss_map != nullptr) {
                    if (inplace)
                        ssim_loss_map[global_idx] += ssim_loss_map_weight * (1.0f - val);
                    else
                        ssim_loss_map[global_idx] = 1.0f - val;
                }
            }
            }
        }

        // TODO: block reduce
        // Note that this value does not affect gradient; for display only
        cg::thread_block_tile<WARP_SIZE> warp = cg::tiled_partition<WARP_SIZE>(block);
        val = cg::reduce(warp, val, cg::plus<float>());
        if (warp.thread_rank() == 0)
            atomicAdd(ssim, val/(B*W*H));
    }

}

// ------------------------------------------
// Backward Kernel: Apply chain rule to get
//    dL/d(img1) from partial derivatives
//    (dm_dmu1, dm_dsigma1_sq, dm_dsigma12)
//    and dL/dmap (the gradient from above).
// ------------------------------------------
template<bool inplace>
__global__ void ssim_backward_kernel(
    int B, int H, int W,
    const float3* __restrict__ img1,
    const float3* __restrict__ img2,
    const float dL_dmap,
    float3* __restrict__ dL_dimg1,
    const float3* __restrict__ dm_dmu1,
    const float3* __restrict__ dm_dsigma1_sq,
    const float3* __restrict__ dm_dsigma12
) {
    auto block = cg::this_thread_block();

    float grad = dL_dmap / (3.0f*B*W*H);

    const int pix_y  = block.group_index().y * BLOCK_Y + block.thread_index().y;
    const int pix_x  = block.group_index().x * BLOCK_X + block.thread_index().x;
    const int pix_id = pix_y * W + pix_x;
    const int num_pix = H * W;
    const int bIdx   = block.group_index().z;

    // Shared memory for the fused data:
    // [0]: dm_dmu1*dL, [1]: dm_dsigma1_sq*dL, [2]: dm_dsigma12*dL
    __shared__ float3 sData[3][SHARED_Y][SHARED_X];
    __shared__ float3 sScratch[CONV_Y][CONV_X][3];

    float3 p1 = make_float3(0), p2 = make_float3(0);
    if (pix_x < W && pix_y < H) {
        p1 = get_pix_value(img1, bIdx, pix_y, pix_x, H, W);
        p2 = get_pix_value(img2, bIdx, pix_y, pix_x, H, W);
    }

    // (1) Load + fuse multiplication
    {
        const int start_y = block.group_index().y * BLOCK_Y;
        const int start_x = block.group_index().x * BLOCK_X;

        int tid = threadIdx.y * blockDim.x + threadIdx.x;
        int warp_id = tid / 32;
        int lane_id = tid % 32;
        int totalThreads = BLOCK_X * BLOCK_Y;
        int num_warps = (totalThreads + 31) / 32;

        for (int row = warp_id; row < SHARED_Y; row += num_warps) {
            int gy = start_y + row - HALO;
            for (int col = lane_id; col < SHARED_X; col += 32) {
                int gx = start_x + col - HALO;

                float3 vmu   = get_pix_value(dm_dmu1,      bIdx, gy, gx, H, W);
                float3 vs1   = get_pix_value(dm_dsigma1_sq,bIdx, gy, gx, H, W);
                float3 vs12  = get_pix_value(dm_dsigma12,  bIdx, gy, gx, H, W);

                sData[0][row][col] = vmu  * grad;
                sData[1][row][col] = vs1  * grad;
                sData[2][row][col] = vs12 * grad;
            }
        }
    }
    block.sync();

    // (2) Horizontal pass
    {
        int ly = threadIdx.y;
        int lx = threadIdx.x + HALO;

        for (int pass = 0; pass < 2; ++pass) {
            int yy = ly + pass * BLOCK_Y;
            if (yy < CONV_Y) {
                float3 accum0 = make_float3(0), accum1 = make_float3(0), accum2 = make_float3(0);

                #pragma unroll
                for (int d = 1; d <= HALO; ++d) {
                    float w = cGauss[HALO - d];
                    float3 left0  = sData[0][yy][lx - d];
                    float3 left1  = sData[1][yy][lx - d];
                    float3 left2  = sData[2][yy][lx - d];

                    float3 right0 = sData[0][yy][lx + d];
                    float3 right1 = sData[1][yy][lx + d];
                    float3 right2 = sData[2][yy][lx + d];

                    accum0 += (left0 + right0) * w;
                    accum1 += (left1 + right1) * w;
                    accum2 += (left2 + right2) * w;
                }
                // center
                {
                    float wc = cGauss[HALO];
                    float3 c0 = sData[0][yy][lx];
                    float3 c1 = sData[1][yy][lx];
                    float3 c2 = sData[2][yy][lx];
                    accum0 += c0 * wc;
                    accum1 += c1 * wc;
                    accum2 += c2 * wc;
                }

                sScratch[yy][threadIdx.x][0] = accum0;
                sScratch[yy][threadIdx.x][1] = accum1;
                sScratch[yy][threadIdx.x][2] = accum2;
            }
        }
    }
    block.sync();

    // (3) Vertical pass -> finalize dL/d(img1)
    if (pix_x < W && pix_y < H) {
        int ly = threadIdx.y + HALO;
        int lx = threadIdx.x;

        float3 sum0 = make_float3(0), sum1 = make_float3(0), sum2 = make_float3(0);

        #pragma unroll
        for (int d = 1; d <= HALO; ++d) {
            float w = cGauss[HALO - d];
            float3* top = sScratch[ly - d][lx];
            float3* bot = sScratch[ly + d][lx];

            sum0 += (top[0] + bot[0]) * w;
            sum1 += (top[1] + bot[1]) * w;
            sum2 += (top[2] + bot[2]) * w;
        }
        // center
        {
            float wc = cGauss[HALO];
            float3* ctr = sScratch[ly][lx];
            sum0 += ctr[0] * wc;
            sum1 += ctr[1] * wc;
            sum2 += ctr[2] * wc;
        }

        // final accumulation
        float3 dL_dpix = sum0 + (2.f * p1) * sum1 + (p2) * sum2;

        int out_idx = bIdx * num_pix + pix_id;
        if (inplace)
            dL_dimg1[out_idx] += dL_dpix;
        else
            dL_dimg1[out_idx] = dL_dpix;
    }

}

#define BLOCK_X_ME 24
#define BLOCK_Y_ME 24

#define SHARED_X_ME (BLOCK_X_ME + 2 * HALO)
#define SHARED_Y_ME (BLOCK_Y_ME + 2 * HALO)

#define CONV_X_ME BLOCK_X_ME
#define CONV_Y_ME SHARED_Y_ME

#define LOAD_X_ME (BLOCK_X_ME + 4 * HALO)
#define LOAD_Y_ME (BLOCK_Y_ME + 4 * HALO)

template<bool inplace>
__global__ void memory_efficient_ssim_backward_kernel(
    int B, int H, int W,
    const float3* __restrict__ img1,   // [B, H, W, 3]
    const float3* __restrict__ img2,   // [B, H, W, 3]
    const float dL_dmap, // [1]
    float3* __restrict__ dL_dimg1      // [B, H, W, 3]
) {
    auto block = cg::this_thread_block();

    const int bIdx   = block.group_index().z;
    const int pix_y  = block.group_index().y * BLOCK_Y_ME + block.thread_index().y;
    const int pix_x  = block.group_index().x * BLOCK_X_ME + block.thread_index().x;
    const int pix_id = pix_y * W + pix_x;
    const int num_pix = H * W;

    const float grad = dL_dmap / (3.0f * B * W * H);

    // __shared__ float sTile[LOAD_Y][LOAD_X][2];
    // __shared__ float xconv[LOAD_Y][SHARED_X][5];
    // __shared__ float sData[3][SHARED_Y][SHARED_X];
    // __shared__ float sScratch[CONV_Y][CONV_X][3];

    static constexpr int sTile_size = LOAD_Y_ME*LOAD_X_ME*2;
    static constexpr int xconv_size = LOAD_Y_ME*SHARED_X_ME*5;
    static constexpr int sData_size = 3*SHARED_Y_ME*SHARED_X_ME;
    static constexpr int sScratch_size = CONV_Y_ME*CONV_X_ME*3;
    __shared__ float _shared_buffer1[sTile_size>sData_size?sTile_size:sData_size];
    __shared__ float _shared_buffer2[xconv_size>sScratch_size?xconv_size:sScratch_size];
    #define sTile(i, j, k) _shared_buffer1[((i)*LOAD_X_ME+(j))*2+(k)]
    #define xconv(i, j, k) _shared_buffer2[((i)*SHARED_X_ME+(j))*5+(k)]
    #define sData(i, j, k) _shared_buffer1[((i)*SHARED_Y_ME+(j))*SHARED_X_ME+(k)]
    #define sScratch(i, j, k) _shared_buffer2[((i)*CONV_X_ME+(j))*3+(k)]

    // TODO: some shared memory can be reused to improve occupancy

    #pragma unroll
    for (int ci = 0; ci < 3; ci++) {
        if (ci != 0)
            block.sync();

    // ------------------------------------------------------------------
    // 1) Load (img1, img2) tile + halo into shared memory
    // ------------------------------------------------------------------
    {
        const int tileSize = LOAD_X_ME * LOAD_Y_ME;
        const int threads = BLOCK_X_ME * BLOCK_Y_ME;
        const int steps = (tileSize + threads - 1) / threads;

        const int tileStartY = block.group_index().y * BLOCK_Y_ME;
        const int tileStartX = block.group_index().x * BLOCK_X_ME;

        for (int s = 0; s < steps; ++s) {
            int tid = s * threads + block.thread_rank();
            if (tid < tileSize) {
                int local_y = tid / LOAD_X_ME;
                int local_x = tid % LOAD_X_ME;
                int gy = tileStartY + local_y - 2 * HALO;
                int gx = tileStartX + local_x - 2 * HALO;

                float X = get_pix_value(img1, bIdx, gy, gx, ci, H, W);
                float Y = get_pix_value(img2, bIdx, gy, gx, ci, H, W);

                sTile(local_y, local_x, 0) = X;
                sTile(local_y, local_x, 1) = Y;
            }
        }
    }
    block.sync();

    // ------------------------------------------------------------------
    // 2) Horizontal convolution (11x1) in shared memory
    //    We'll accumulate symmetrical pairs around center.
    // ------------------------------------------------------------------
    #if 0
    for (int ly = threadIdx.y; ly < LOAD_Y_ME; ly += BLOCK_Y_ME)
    for (int lx = threadIdx.x + HALO; lx + HALO < LOAD_X_ME; lx += BLOCK_X_ME) {
    #else
    for (int l = block.thread_rank(); l < SHARED_X_ME*LOAD_Y_ME; l += BLOCK_X_ME*BLOCK_Y_ME) {
        int lx = l % SHARED_X_ME + HALO;
        int ly = l / SHARED_Y_ME;
    #endif

        float sumX   = 0.0f;
        float sumX2  = 0.0f;
        float sumY   = 0.0f;
        float sumY2  = 0.0f;
        float sumXY  = 0.0f;

        // #pragma unroll for those 5 pairs
        #pragma unroll
        for (int d = 1; d <= HALO; ++d) {
            float w = cGauss[HALO - d];
            float Xleft  = sTile(ly, lx - d, 0);
            float Yleft  = sTile(ly, lx - d, 1);
            float Xright = sTile(ly, lx + d, 0);
            float Yright = sTile(ly, lx + d, 1);

            sumX  += (Xleft + Xright) * w;
            sumX2 += ((Xleft * Xleft) + (Xright * Xright)) * w;
            sumY  += (Yleft + Yright) * w;
            sumY2 += ((Yleft * Yleft) + (Yright * Yright)) * w;
            sumXY += ((Xleft * Yleft) + (Xright * Yright)) * w;
        }
        // center
        {
            float centerX = sTile(ly, lx, 0);
            float centerY = sTile(ly, lx, 1);
            float wc = cGauss[HALO];
            sumX  += centerX * wc;
            sumX2 += (centerX * centerX) * wc;
            sumY  += centerY * wc;
            sumY2 += (centerY * centerY) * wc;
            sumXY += (centerX * centerY) * wc;
        }

        // Write out partial sums
        xconv(ly, lx-HALO, 0) = sumX;
        xconv(ly, lx-HALO, 1) = sumX2;
        xconv(ly, lx-HALO, 2) = sumY;
        xconv(ly, lx-HALO, 3) = sumY2;
        xconv(ly, lx-HALO, 4) = sumXY;

    }
    block.sync();

    // ------------------------------------------------------------
    // 3) Vertical convolution (1x11) + final SSIM
    // ------------------------------------------------------------
    #if 0
    for (int ly = threadIdx.y + HALO; ly + HALO < LOAD_Y_ME; ly += BLOCK_Y_ME)
    for (int lx = threadIdx.x; lx < SHARED_X_ME; lx += BLOCK_X_ME) {
    #else
    for (int l = block.thread_rank(); l < SHARED_X_ME*SHARED_Y_ME; l += BLOCK_X_ME*BLOCK_Y_ME) {
        int lx = l % SHARED_X_ME;
        int ly = l / SHARED_Y_ME + HALO;
    #endif

        float out0 = 0.0f, out1 = 0.0f, out2 = 0.0f, out3 = 0.0f, out4 = 0.0f;

        #pragma unroll
        for (int d = 1; d <= HALO; ++d) {
            float w = cGauss[HALO - d];

            out0 += (xconv(ly-d, lx, 0) + xconv(ly+d, lx, 0)) * w;
            out1 += (xconv(ly-d, lx, 1) + xconv(ly+d, lx, 1)) * w;
            out2 += (xconv(ly-d, lx, 2) + xconv(ly+d, lx, 2)) * w;
            out3 += (xconv(ly-d, lx, 3) + xconv(ly+d, lx, 3)) * w;
            out4 += (xconv(ly-d, lx, 4) + xconv(ly+d, lx, 4)) * w;
        }
        // center
        {
            float wC = cGauss[HALO];
            out0 += xconv(ly, lx, 0) * wC;
            out1 += xconv(ly, lx, 1) * wC;
            out2 += xconv(ly, lx, 2) * wC;
            out3 += xconv(ly, lx, 3) * wC;
            out4 += xconv(ly, lx, 4) * wC;
        }

        static constexpr float kC1 = 0.01f * 0.01f;
        static constexpr float kC2 = 0.03f * 0.03f;

        float mu1 = out0;
        float mu2 = out2;
        float mu1_sq = mu1 * mu1;
        float mu2_sq = mu2 * mu2;

        float sigma1_sq = out1 - mu1_sq;
        float sigma2_sq = out3 - mu2_sq;
        float sigma12   = out4 - mu1 * mu2;

        float A = mu1_sq + mu2_sq + kC1;
        float B = sigma1_sq + sigma2_sq + kC2;
        float C_ = 2.f * mu1 * mu2 + kC1;
        float D_ = 2.f * sigma12 + kC2;

        // partial derivatives
        float d_m_dmu1 = (
            (mu2 * 2.f * D_) / (A * B)
            - (mu2 * 2.f * C_) / (A * B)
            - (mu1 * 2.f * C_ * D_) / (A * A * B)
            + (mu1 * 2.f * C_ * D_) / (A * B * B)
        );
        float d_m_dsigma1_sq = (-C_ * D_) / (A * B * B);
        float d_m_dsigma12   = (2.f * C_) / (A * B);

        const int pix_y  = block.group_index().y * BLOCK_Y_ME + ly-HALO;
        const int pix_x  = block.group_index().x * BLOCK_X_ME + lx;
        float masked_grad = grad * float(
            pix_x >= HALO && pix_y >= HALO && pix_x < W+HALO && pix_y < H+HALO
        );

        sData(0, ly-HALO, lx) = d_m_dmu1 * masked_grad;
        sData(1, ly-HALO, lx) = d_m_dsigma1_sq * masked_grad;
        sData(2, ly-HALO, lx) = d_m_dsigma12 * masked_grad;
    }
    block.sync();

    // (2) Horizontal pass
    {
        int ly = threadIdx.y;
        int lx = threadIdx.x + HALO;

        for (int pass = 0; pass < 2; ++pass) {
            int yy = ly + pass * BLOCK_Y_ME;
            if (yy < CONV_Y_ME) {
                float accum0 = 0.0f, accum1 = 0.0f, accum2 = 0.0f;

                #pragma unroll
                for (int d = 1; d <= HALO; ++d) {
                    float w = cGauss[HALO - d];
                    float left0  = sData(0, yy, lx - d);
                    float left1  = sData(1, yy, lx - d);
                    float left2  = sData(2, yy, lx - d);

                    float right0 = sData(0, yy, lx + d);
                    float right1 = sData(1, yy, lx + d);
                    float right2 = sData(2, yy, lx + d);

                    accum0 += (left0 + right0) * w;
                    accum1 += (left1 + right1) * w;
                    accum2 += (left2 + right2) * w;
                }
                // center
                {
                    float wc = cGauss[HALO];
                    float c0 = sData(0, yy, lx);
                    float c1 = sData(1, yy, lx);
                    float c2 = sData(2, yy, lx);
                    accum0 += c0 * wc;
                    accum1 += c1 * wc;
                    accum2 += c2 * wc;
                }

                sScratch(yy, threadIdx.x, 0) = accum0;
                sScratch(yy, threadIdx.x, 1) = accum1;
                sScratch(yy, threadIdx.x, 2) = accum2;
            }
        }
    }
    block.sync();

    // (3) Vertical pass -> finalize dL/d(img1)
    if (pix_x < W && pix_y < H) {
        int ly = threadIdx.y + HALO;
        int lx = threadIdx.x;

        float sum0 = 0.0f, sum1 = 0.0f, sum2 = 0.0f;

        #pragma unroll
        for (int d = 1; d <= HALO; ++d) {
            float w = cGauss[HALO - d];
            sum0 += (sScratch(ly - d, lx, 0) + sScratch(ly + d, lx, 0)) * w;
            sum1 += (sScratch(ly - d, lx, 1) + sScratch(ly + d, lx, 1)) * w;
            sum2 += (sScratch(ly - d, lx, 2) + sScratch(ly + d, lx, 2)) * w;
        }
        // center
        {
            float wc = cGauss[HALO];
            sum0 += sScratch(ly, lx, 0) * wc;
            sum1 += sScratch(ly, lx, 1) * wc;
            sum2 += sScratch(ly, lx, 2) * wc;
        }

        // final accumulation
        float p1 = get_pix_value(img1, bIdx, pix_y, pix_x, ci, H, W);
        float p2 = get_pix_value(img2, bIdx, pix_y, pix_x, ci, H, W);
        float dL_dpix = sum0 + (2.f * p1) * sum1 + (p2) * sum2;

        int out_idx = bIdx * num_pix + pix_id;
        if (inplace)
            ((float*)dL_dimg1)[out_idx*3+ci] += dL_dpix;
        else
            ((float*)dL_dimg1)[out_idx*3+ci] = dL_dpix;
    }

    }  // for (int ci = 0; ci < 3; ci++)

    #undef sTile
    #undef xconv
    #undef sData
    #undef sScratch
}

// ------------------------------------------
// PyTorch Interface (Forward)
//   Returns (ssim, ssim_loss_map (single channel non differentiable), dm_dmu1, dm_dsigma1_sq, dm_dsigma12).
//   If train=false, derivative Tensors are empty.
// ------------------------------------------
std::tuple<at::Tensor, std::optional<at::Tensor>, at::Tensor, at::Tensor, at::Tensor>
fused_ssim_forward(
    at::Tensor &img1,
    at::Tensor &img2,
    bool train,
    bool return_ssim_loss_map,
    bool is_l1
) {
    DEVICE_GUARD(img1);
    CHECK_INPUT(img1);
    CHECK_INPUT(img2);

    int B  = img1.size(0);
    int H  = img1.size(1);
    int W  = img1.size(2);
    int CH = img1.size(3);
    if (CH != 3)
        throw std::runtime_error("Image must be (B, H, W, 3)");

    // Output SSIM map
    at::Tensor ssim = at::zeros({}, img1.options());
    std::optional<at::Tensor> ssim_loss_map;
    if (return_ssim_loss_map) {
        ssim_loss_map = at::empty({B, H, W, 1}, img1.options());
        set_zero<float>(ssim_loss_map.value());
    }

    // Optionally allocate derivative Tensors
    at::Tensor dm_dmu1       = train ? at::empty_like(img1) : at::empty({0}, img1.options());
    at::Tensor dm_dsigma1_sq = train ? at::empty_like(img1) : at::empty({0}, img1.options());
    at::Tensor dm_dsigma12   = train ? at::empty_like(img1) : at::empty({0}, img1.options());

    (is_l1 ? ssim_forward_kernel<true, false> : ssim_forward_kernel<false, false>)
    <<<_LAUNCH_ARGS_3D(W, H, B, BLOCK_X, BLOCK_Y, 1)>>>(
        B, H, W,
        (float3*)img1.data_ptr<float>(),
        (float3*)img2.data_ptr<float>(),
        ssim.data_ptr<float>(),
        1.0f,
        return_ssim_loss_map ? ssim_loss_map.value().data_ptr<float>() : nullptr,
        train ? (float3*)dm_dmu1.data_ptr<float>()       : nullptr,
        train ? (float3*)dm_dsigma1_sq.data_ptr<float>() : nullptr,
        train ? (float3*)dm_dsigma12.data_ptr<float>()   : nullptr
    );

    return std::make_tuple(ssim, ssim_loss_map, dm_dmu1, dm_dsigma1_sq, dm_dsigma12);
}

std::tuple<at::Tensor, at::Tensor, at::Tensor, at::Tensor>
fused_ssim_forward_inplace(
    at::Tensor &img1,
    at::Tensor &img2,
    bool train,
    float ssim_loss_map_weight,
    at::Tensor &ssim_loss_map,
    bool is_l1
) {
    DEVICE_GUARD(img1);
    CHECK_INPUT(img1);
    CHECK_INPUT(img2);
    CHECK_INPUT(ssim_loss_map);

    int B  = img1.size(0);
    int H  = img1.size(1);
    int W  = img1.size(2);
    int CH = img1.size(3);
    if (CH != 3)
        throw std::runtime_error("Image must be (B, H, W, 3)");

    // Output SSIM map
    at::Tensor ssim = at::zeros({}, img1.options());

    // Optionally allocate derivative Tensors
    at::Tensor dm_dmu1       = train ? at::empty_like(img1) : at::empty({0}, img1.options());
    at::Tensor dm_dsigma1_sq = train ? at::empty_like(img1) : at::empty({0}, img1.options());
    at::Tensor dm_dsigma12   = train ? at::empty_like(img1) : at::empty({0}, img1.options());

    (is_l1 ? ssim_forward_kernel<true, true> : ssim_forward_kernel<false, true>)
    <<<_LAUNCH_ARGS_3D(W, H, B, BLOCK_X, BLOCK_Y, 1)>>>(
        B, H, W,
        (float3*)img1.data_ptr<float>(),
        (float3*)img2.data_ptr<float>(),
        ssim.data_ptr<float>(),
        ssim_loss_map_weight,
        ssim_loss_map.data_ptr<float>(),
        train ? (float3*)dm_dmu1.data_ptr<float>()       : nullptr,
        train ? (float3*)dm_dsigma1_sq.data_ptr<float>() : nullptr,
        train ? (float3*)dm_dsigma12.data_ptr<float>()   : nullptr
    );

    return std::make_tuple(ssim, dm_dmu1, dm_dsigma1_sq, dm_dsigma12);
}

// ------------------------------------------
// PyTorch Interface (Backward)
//   Takes the gradient wrt the SSIM map and
//   the partial derivatives from forward;
//   returns dL/d(img1).
// ------------------------------------------
at::Tensor
fused_ssim_backward(
    at::Tensor &img1,
    at::Tensor &img2,
    const float dL_dmap,
    std::optional<at::Tensor> &dm_dmu1,
    std::optional<at::Tensor> &dm_dsigma1_sq,
    std::optional<at::Tensor> &dm_dsigma12
) {
    DEVICE_GUARD(img1);
    CHECK_INPUT(img1);
    CHECK_INPUT(img2);

    int B  = img1.size(0);
    int H  = img1.size(1);
    int W  = img1.size(2);
    int CH = img1.size(3);

    auto dL_dimg1 = at::empty_like(img1);

    if (dm_dmu1.has_value() && dm_dsigma1_sq.has_value() && dm_dsigma12.has_value()) {
        CHECK_INPUT(dm_dmu1.value());
        CHECK_INPUT(dm_dsigma1_sq.value());
        CHECK_INPUT(dm_dsigma12.value());
        ssim_backward_kernel<false><<<_LAUNCH_ARGS_3D(W, H, B, BLOCK_X, BLOCK_Y, 1)>>>(
            B, H, W,
            (float3*)img1.data_ptr<float>(),
            (float3*)img2.data_ptr<float>(),
            dL_dmap,
            (float3*)dL_dimg1.data_ptr<float>(),
            (float3*)dm_dmu1.value().data_ptr<float>(),
            (float3*)dm_dsigma1_sq.value().data_ptr<float>(),
            (float3*)dm_dsigma12.value().data_ptr<float>()
        );
    }
    else {
        memory_efficient_ssim_backward_kernel<false><<<_LAUNCH_ARGS_3D(W, H, B, BLOCK_X_ME, BLOCK_Y_ME, 1)>>>(
            B, H, W,
            (float3*)img1.data_ptr<float>(),
            (float3*)img2.data_ptr<float>(),
            dL_dmap,
            (float3*)dL_dimg1.data_ptr<float>()
        );
    }

    return dL_dimg1;
}

void fused_ssim_backward_inplace(
    at::Tensor &img1,
    at::Tensor &img2,
    const float dL_dmap,
    std::optional<at::Tensor> &dm_dmu1,
    std::optional<at::Tensor> &dm_dsigma1_sq,
    std::optional<at::Tensor> &dm_dsigma12,
    at::Tensor &dL_dimg1
) {
    DEVICE_GUARD(img1);
    CHECK_INPUT(img1);
    CHECK_INPUT(img2);
    CHECK_INPUT(dL_dimg1);

    int B  = img1.size(0);
    int H  = img1.size(1);
    int W  = img1.size(2);
    int CH = img1.size(3);

    if (dm_dmu1.has_value() && dm_dsigma1_sq.has_value() && dm_dsigma12.has_value()) {
        CHECK_INPUT(dm_dmu1.value());
        CHECK_INPUT(dm_dsigma1_sq.value());
        CHECK_INPUT(dm_dsigma12.value());
        ssim_backward_kernel<true><<<_LAUNCH_ARGS_3D(W, H, B, BLOCK_X, BLOCK_Y, 1)>>>(
            B, H, W,
            (float3*)img1.data_ptr<float>(),
            (float3*)img2.data_ptr<float>(),
            dL_dmap,
            (float3*)dL_dimg1.data_ptr<float>(),
            (float3*)dm_dmu1.value().data_ptr<float>(),
            (float3*)dm_dsigma1_sq.value().data_ptr<float>(),
            (float3*)dm_dsigma12.value().data_ptr<float>()
        );
    }
    else {
        memory_efficient_ssim_backward_kernel<true><<<_LAUNCH_ARGS_3D(W, H, B, BLOCK_X_ME, BLOCK_Y_ME, 1)>>>(
            B, H, W,
            (float3*)img1.data_ptr<float>(),
            (float3*)img2.data_ptr<float>(),
            dL_dmap,
            (float3*)dL_dimg1.data_ptr<float>()
        );
    }

}
