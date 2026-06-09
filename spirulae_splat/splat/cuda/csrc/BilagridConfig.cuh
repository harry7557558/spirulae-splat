#pragma once

// Ported from fused-bilagrid/fused_bilagrid/config.h
// Self-contained to avoid pulling in common_utils.cuh, which would conflict
// with vector ops/dot/cross defined in BilagridPpispMath.cuh.

// RGB to gray constants
#define kC2G_r 0.299f
#define kC2G_g 0.587f
#define kC2G_b 0.114f

#ifdef __CUDACC__

#include <stdio.h>
#include <cuda_runtime.h>

// Function-call-style macro (different from the original statement-style),
// matching the form used by the rest of spirulae_splat/Common.cuh.
#ifndef CHECK_DEVICE_ERROR
#define CHECK_DEVICE_ERROR(call)                                    \
do {                                                                \
    cudaError_t err = call;                                         \
    if (err != cudaSuccess) {                                       \
        fprintf(stderr, "\033[41mCUDA Error at %s:%d: %s\033[m\n",  \
                __FILE__, __LINE__, cudaGetErrorString(err));       \
        exit(EXIT_FAILURE);                                         \
    }                                                               \
} while (0)
#endif

#include "BilagridReader.cuh"

// Block dimensions baked into every bilagrid *_sample_backward_v1 kernel
// (affine RGB, PPISP, log-linear, depth, normal). Threaded as compile-time
// constants so blockDim.x/y/z folds inside the kernels (saves register
// pressure on x/y/z index math and lets `mult_x % BX == 0` style checks
// disappear at compile time when call sites pick block-multiple tile sizes).
// Host launchers use these to set `dim3 block`; do NOT change one without
// the other. 8x8x1 = 64 threads per block matches `sharedData[64]` reduction
// buffer in the bilagrid-grad kernels.
constexpr unsigned kBilagridBwdV1BlockX = 8;
constexpr unsigned kBilagridBwdV1BlockY = 8;
constexpr unsigned kBilagridBwdV1BlockZ = 1;

// 1D block dim for the companion *_v1_kernel_rgb (image-side rgb-grad) pass
// in the same 5 v1 variants. Launcher and kernel must agree.
constexpr unsigned kBilagridBwdV1RgbThreads = 256;


// Block-wide reduce + single global atomicAdd per call, used by the 5
// bilagrid sample backward v1 kernels (affine RGB, PPISP, log-linear, depth,
// normal) for their channel-loop writeback. Functionally equivalent to
// `atomicAddFVec<64>` from Common.cuh but declared here so the v1 bwd kernel
// headers can avoid pulling in Common.cuh -- which provides global vector
// `dot`/`operator*` overloads that collide with bilagrid_ppisp::dot inside
// BilagridPpispMath{,Bwd}.cuh.
//
// Hardcoded for the 8x8x1 = 64-thread (2-warp) block geometry used by all 5
// v1 bwd kernels (see kBilagridBwdV1Block{X,Y,Z} in BilagridConfig.cuh).
//
// Caller requirements:
//   - threadblock is exactly 2 warps (8x8x1, 64 threads).
//   - All threads must call together.
//   - A __syncthreads() must separate successive calls within the same
//     kernel invocation: the static __shared__ scratch is reused, and
//     warp 0's reads of scratch[1] from the previous call can race with
//     warp 1's writes of scratch[1] for the next call otherwise.

#include <cuda_runtime.h>

inline __device__ void bilagrid_bwd_v1_block_atomic_add(
    float* addr, float val
) {
    constexpr unsigned kBlockSize = kBilagridBwdV1BlockX * kBilagridBwdV1BlockY * kBilagridBwdV1BlockZ;
    constexpr unsigned kWarpSize = 32;

    val = isfinite(val) ? val : 0.0f;

    // Intra-warp shuffle reduce. Each warp ends with its sum in lane 0.
    #pragma unroll
    for (int stride = kWarpSize >> 1; stride > 0; stride >>= 1) {
        val += __shfl_down_sync(0xFFFFFFFFu, val, stride);
    }

    // Lane 0 of each warp hands off its partial sum to a 2-slot scratch.
    static __shared__ float scratch[kBlockSize / kWarpSize];

    uint32_t lane_id;
    asm("mov.u32 %0, %%laneid;" : "=r"(lane_id));
    // 2D block, linear order: warpId = linear_tid / 32. blockDim is the
    // launch-time (constexpr) 8x8x1, so this collapses to a constant 0 or 1.
    const uint32_t linear_tid =
        (threadIdx.z * kBilagridBwdV1BlockY + threadIdx.y) * kBilagridBwdV1BlockX + threadIdx.x;
    const uint32_t warp_id = linear_tid / kWarpSize;

    if (lane_id == 0u) scratch[warp_id] = val;
    __syncthreads();

    // Warp 0 reduces the 2 partials and emits a single global atomicAdd.
    if (warp_id == 0u) {
        float sum = (lane_id < (kBlockSize / kWarpSize)) ? scratch[lane_id] : 0.0f;
        sum += __shfl_down_sync(0xFFFFFFFFu, sum, 1);
        if (lane_id == 0u && sum != 0.0f) atomicAdd(addr, sum);
    }
}

#endif  // __CUDACC__
