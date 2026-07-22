#include <cuda_runtime.h>
#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>

#include "kernels/bilagrid/BilagridConfig.cuh"

namespace cg = cooperative_groups;


template<int C>
__global__ void tv_loss_forward_kernel(
    BilagridReader bilagrid,  // [N,L,H,W,C] (channel-last)
    float* __restrict__ tv_loss,
    int N, int L, int H, int W
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int wi = blockIdx.y * blockDim.y + threadIdx.y;
    int hi = blockIdx.z * blockDim.z + threadIdx.z;
    // bool inside = (wi < W && hi < H && idx < (L*C*N));
    bool inside = (wi < W && hi < H && idx < (L*N));
    int li = idx % L; idx /= L;
    // int ci = idx % C; idx /= C;
    int ni = idx;

    // Channel-last strides: ci=1, wi=C, hi=W*C, li=H*W*C, ni=L*H*W*C.
    const int sw = C;
    const int sh = W * C;
    const int sl = H * W * C;
    const int sn = L * H * W * C;

    float tv_sum = 0.0f;

    if (inside) {
        const int cell_base = ni * sn + li * sl + hi * sh + wi * sw;
        #pragma unroll
        for (int ci = 0; ci < C; ci++) {

        int cell_idx = cell_base + ci;

        float val = bilagrid[cell_idx];

        if (wi > 0) {
            float val0 = bilagrid[cell_idx - sw];
            float l2 = (val-val0) * (val-val0);
            tv_sum += l2 / (L*H*(W-1));
        }
        if (hi > 0) {
            float val0 = bilagrid[cell_idx - sh];
            float l2 = (val-val0) * (val-val0);
            tv_sum += l2 / (L*(H-1)*W);
        }
        if (li > 0) {
            float val0 = bilagrid[cell_idx - sl];
            float l2 = (val-val0) * (val-val0);
            tv_sum += l2 / ((L-1)*H*W);
        }

        }  // ci
        tv_sum /= (C*N);
    }

#if 0
    auto block = cg::this_thread_block();
    cg::thread_block_tile<32> warp = cg::tiled_partition<32>(block);

    tv_sum = cg::reduce(warp, tv_sum, cg::plus<float>());
    if (warp.thread_rank() == 0)
        atomicAdd(tv_loss, tv_sum);
#else
    static constexpr int blockSize = 64;
    __shared__ float sharedData[blockSize];

    int tid = (threadIdx.z * blockDim.y + threadIdx.y) * blockDim.x + threadIdx.x;

    sharedData[tid] = tv_sum;
    __syncthreads();

    #pragma unroll
    for (int s = blockSize / 2; s > 0; s >>= 1) {
        if (tid < s)
            sharedData[tid] += sharedData[tid + s];
        __syncthreads();
    }

    if (tid == 0)
        atomicAdd(tv_loss, sharedData[0]);
#endif
}


void tv_loss_forward(
    BilagridReader bilagrid,
    float* tv_loss,
    int N, int C, int L, int H, int W,
    cudaStream_t stream
) {
    // TODO: optimize memory access pattern
    dim3 block = { 4, 4, 4 };
    dim3 bounds = {
        // (N*C*L +block.z-1)/block.z,
        (N*L +block.z-1)/block.z,
        (W +block.x-1)/block.x,
        (H +block.y-1)/block.y,
    };
    if (C == 12)
        tv_loss_forward_kernel<12><<<bounds, block, 0, stream>>>(
            bilagrid, tv_loss,
            N, L, H, W
        );
    else if (C == 9)
        tv_loss_forward_kernel<9><<<bounds, block, 0, stream>>>(
            bilagrid, tv_loss,
            N, L, H, W
        );
    else if (C == 2)
        tv_loss_forward_kernel<2><<<bounds, block, 0, stream>>>(
            bilagrid, tv_loss,
            N, L, H, W
        );
    else if (C == 3)
        tv_loss_forward_kernel<3><<<bounds, block, 0, stream>>>(
            bilagrid, tv_loss,
            N, L, H, W
        );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}



template<int C>
__global__ void channel_mean_forward_kernel(
    BilagridReader bilagrid,  // [N,L,H,W,C] (channel-last)
    float* __restrict__ channel_mean,  // [C]
    int N, int L, int H, int W
) {
    int wi = blockIdx.x * blockDim.x + threadIdx.x;
    int hi = blockIdx.y * blockDim.y + threadIdx.y;
    int idx = blockIdx.z * blockDim.z + threadIdx.z;
    // bool inside = (wi < W && hi < H && idx < (L*C*N));
    bool inside = (wi < W && hi < H && idx < (L*N));
    int li = idx % L; idx /= L;
    // int ci = idx % C; idx /= C;
    int ni = idx;

    static constexpr int blockSize = 64;
    __shared__ float sharedData[blockSize];
    int tid = (threadIdx.z * blockDim.y + threadIdx.y) * blockDim.x + threadIdx.x;

    const int cell_base = ((((ni * L) + li) * H + hi) * W + wi) * C;

    #pragma unroll
    for (int ci = 0; ci < C; ci++) {

        float val = 0.0f;

        if (inside) {
            val = bilagrid[cell_base + ci];
        }

        sharedData[tid] = val;
        __syncthreads();

        #pragma unroll
        for (int s = blockSize / 2; s > 0; s >>= 1) {
            if (tid < s)
                sharedData[tid] += sharedData[tid + s];
            __syncthreads();
        }

        if (tid == 0)
            atomicAdd(&channel_mean[ci], sharedData[0] / (float)(N*L*H*W));
    }

}


void channel_mean_forward(
    BilagridReader bilagrid,
    float* channel_mean,
    int N, int C, int L, int H, int W,
    cudaStream_t stream
) {
    // TODO: optimize memory access pattern
    dim3 block = { 4, 4, 4 };
    dim3 bounds = {
        (W +block.x-1)/block.x,
        (H +block.y-1)/block.y,
        // (N*C*L +block.z-1)/block.z
        (N*L +block.z-1)/block.z
    };
    if (C == 12)
        channel_mean_forward_kernel<12><<<bounds, block, 0, stream>>>(
            bilagrid, channel_mean,
            N, L, H, W
        );
    else if (C == 9)
        channel_mean_forward_kernel<9><<<bounds, block, 0, stream>>>(
            bilagrid, channel_mean,
            N, L, H, W
        );
    else if (C == 2)
        channel_mean_forward_kernel<2><<<bounds, block, 0, stream>>>(
            bilagrid, channel_mean,
            N, L, H, W
        );
    else if (C == 3)
        channel_mean_forward_kernel<3><<<bounds, block, 0, stream>>>(
            bilagrid, channel_mean,
            N, L, H, W
        );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}
