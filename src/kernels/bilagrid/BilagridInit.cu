// Identity initialization for bilagrid grids used by Engine.cpp.
//
// For the "affine" bilagrid (12 channels per cell), each cell holds a 3x4
// affine matrix; the identity is [[1,0,0,0],[0,1,0,0],[0,0,1,0]].
// "ppisp" / "loglinear" / "depth" / "normal" grids are initialized to zero
// (which is the identity transform for those parameterizations) and use a
// plain cudaMemset on the launcher side, not this kernel.

#include <cuda_runtime.h>
#include <cstdint>

#include "kernels/bilagrid/BilagridConfig.cuh"


// Flat over the N*L*H*W cells, one thread per cell. A 3D grid with N*L on
// gridDim.z would cap the camera count at 65535/L (~8k images for L=8), so
// the cell axis is 1D -- gridDim.x reaches 2^31-1. Matches the Vulkan
// dispatch_flat geometry in backend/vulkan/kernels/Bilagrid.cpp.
__global__ void bilagrid_affine_identity_init_kernel(
    float* __restrict__ grids,   // [N, L, H, W, 12]
    int64_t total_cells          // N*L*H*W
) {
    // Identity values for the 12-channel affine bilagrid (row-major 3x4).
    static constexpr float kIdent[12] = {
        1.0f, 0.0f, 0.0f, 0.0f,
        0.0f, 1.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 1.0f, 0.0f
    };

    const int64_t idx = (int64_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= total_cells) return;

    // Channel-last: cell base in floats is the cell index * 12.
    const int64_t cell = idx * 12;

    #pragma unroll
    for (int ci = 0; ci < 12; ci++) {
        grids[cell + ci] = kIdent[ci];
    }
}


void bilagrid_affine_identity_init(
    float* grids, int N, int L, int H, int W, cudaStream_t stream
) {
    if (N <= 0 || L <= 0 || H <= 0 || W <= 0) return;
    constexpr int kThreads = 256;
    const int64_t total_cells = (int64_t)N * L * H * W;
    const int64_t blocks = (total_cells + kThreads - 1) / kThreads;
    bilagrid_affine_identity_init_kernel<<<(unsigned)blocks, kThreads, 0, stream>>>(
        grids, total_cells);
    CHECK_DEVICE_ERROR(cudaGetLastError());
}
