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

// ------------------------------------------
// Constant Memory for Gaussian Coefficients
// ------------------------------------------
__constant__ float cGauss[11] = {
    0.001028380123898387f,
    0.0075987582094967365f,
    0.036000773310661316f,
    0.10936068743467331f,
    0.21300552785396576f,
    0.26601171493530273f,
    0.21300552785396576f,
    0.10936068743467331f,
    0.036000773310661316f,
    0.0075987582094967365f,
    0.001028380123898387f
};

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

// ------------------------------------------
// Forward Kernel: Fused SSIM
//  - Two-pass convolution to get mu1, mu2,
//    sigma1_sq, sigma2_sq, sigma12, etc.
//  - Writes final SSIM map to ssim_map
//  - Optionally writes partial derivatives
//    to dm_dmu1, dm_dsigma1_sq, dm_dsigma12
// ------------------------------------------
__global__ void fusedssimCUDA(
    int B, int H, int W,
    const float3* __restrict__ img1,
    const float3* __restrict__ img2,
    float* __restrict__ ssim_map,
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
            }
        }

        // TODO: block reduce
        // Note that this value does not affect gradient; for display only
        cg::thread_block_tile<WARP_SIZE> warp = cg::tiled_partition<WARP_SIZE>(block);
        val = cg::reduce(warp, val, cg::plus<float>());
        if (warp.thread_rank() == 0)
            atomicAdd(ssim_map, val/(B*W*H));
    }

}

// ------------------------------------------
// Backward Kernel: Apply chain rule to get
//    dL/d(img1) from partial derivatives
//    (dm_dmu1, dm_dsigma1_sq, dm_dsigma12)
//    and dL/dmap (the gradient from above).
// ------------------------------------------
__global__ void fusedssim_backwardCUDA(
    int B, int H, int W,
    const float3* __restrict__ img1,
    const float3* __restrict__ img2,
    const float* __restrict__ dL_dmap,
    float3* __restrict__ dL_dimg1,
    const float3* __restrict__ dm_dmu1,
    const float3* __restrict__ dm_dsigma1_sq,
    const float3* __restrict__ dm_dsigma12
) {
    auto block = cg::this_thread_block();

    float grad = dL_dmap[0] / (3.0f*B*W*H);

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
        dL_dimg1[out_idx] = dL_dpix;
    }

}

// ------------------------------------------
// PyTorch Interface (Forward)
//   Returns (ssim_map, dm_dmu1, dm_dsigma1_sq, dm_dsigma12).
//   If train=false, derivative Tensors are empty.
// ------------------------------------------
std::tuple<at::Tensor, at::Tensor, at::Tensor, at::Tensor>
fused_ssim_forward(
    at::Tensor &img1,
    at::Tensor &img2,
    bool train
) {
    const at::cuda::OptionalCUDAGuard device_guard(device_of(img1));
    int B  = img1.size(0);
    int H  = img1.size(1);
    int W  = img1.size(2);
    int CH = img1.size(3);
    if (CH != 3)
        throw std::runtime_error("Image must be (B, H, W, 3)");

    // Launch config
    dim3 grid((W + BLOCK_X - 1) / BLOCK_X,
              (H + BLOCK_Y - 1) / BLOCK_Y,
              B);
    dim3 block(BLOCK_X, BLOCK_Y);

    // Output SSIM map
    // auto ssim_map = at::empty_like(img1, img1.options());
    auto ssim_map = at::zeros({}, img1.options());

    // Optionally allocate derivative Tensors
    auto dm_dmu1       = train ? at::empty_like(img1) : at::empty({0}, img1.options());
    auto dm_dsigma1_sq = train ? at::empty_like(img1) : at::empty({0}, img1.options());
    auto dm_dsigma12   = train ? at::empty_like(img1) : at::empty({0}, img1.options());

    fusedssimCUDA<<<grid, block>>>(
        B, H, W,
        (float3*)img1.data_ptr<float>(),
        (float3*)img2.data_ptr<float>(),
        // (float3*)ssim_map.data_ptr<float>(),
        ssim_map.data_ptr<float>(),
        train ? (float3*)dm_dmu1.data_ptr<float>()       : nullptr,
        train ? (float3*)dm_dsigma1_sq.data_ptr<float>() : nullptr,
        train ? (float3*)dm_dsigma12.data_ptr<float>()   : nullptr
    );

    return std::make_tuple(ssim_map, dm_dmu1, dm_dsigma1_sq, dm_dsigma12);
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
    at::Tensor &dL_dmap,
    at::Tensor &dm_dmu1,
    at::Tensor &dm_dsigma1_sq,
    at::Tensor &dm_dsigma12
) {
    const at::cuda::OptionalCUDAGuard device_guard(device_of(img1));
    int B  = img1.size(0);
    int H  = img1.size(1);
    int W  = img1.size(2);
    int CH = img1.size(3);

    auto dL_dimg1 = at::empty_like(img1);

    dim3 grid((W + BLOCK_X - 1) / BLOCK_X,
              (H + BLOCK_Y - 1) / BLOCK_Y,
              B);
    dim3 block(BLOCK_X, BLOCK_Y);

    fusedssim_backwardCUDA<<<grid, block>>>(
        B, H, W,
        (float3*)img1.data_ptr<float>(),
        (float3*)img2.data_ptr<float>(),
        dL_dmap.data_ptr<float>(),
        (float3*)dL_dimg1.data_ptr<float>(),
        (float3*)dm_dmu1.data_ptr<float>(),
        (float3*)dm_dsigma1_sq.data_ptr<float>(),
        (float3*)dm_dsigma12.data_ptr<float>()
    );

    return dL_dimg1;
}
