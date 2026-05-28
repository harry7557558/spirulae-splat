// Identity initialization for bilagrid grids used by Engine.cpp.
//
// For the "affine" bilagrid (12 channels per cell), each cell holds a 3x4
// affine matrix; the identity is [[1,0,0,0],[0,1,0,0],[0,0,1,0]].
// "ppisp" / "loglinear" / "depth" / "normal" grids are initialized to zero
// (which is the identity transform for those parameterizations) and use a
// plain cudaMemset on the launcher side, not this kernel.

#include <cuda_runtime.h>
#include <cstdint>

#include "BilagridConfig.cuh"


__global__ void bilagrid_affine_identity_init_kernel(
    float* __restrict__ grids,   // [N, 12, L, H, W]
    int N, int L, int H, int W
) {
    // Identity values for the 12-channel affine bilagrid (row-major 3x4).
    static constexpr float kIdent[12] = {
        1.0f, 0.0f, 0.0f, 0.0f,
        0.0f, 1.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 1.0f, 0.0f
    };

    int wi = blockIdx.x * blockDim.x + threadIdx.x;
    int hi = blockIdx.y * blockDim.y + threadIdx.y;
    int li_n = blockIdx.z * blockDim.z + threadIdx.z;
    if (wi >= W || hi >= H || li_n >= L * N) return;
    int li = li_n % L;
    int ni = li_n / L;

    int64_t base_n = (int64_t)ni * 12 * L * H * W;
    int64_t plane = (int64_t)L * H * W;
    int64_t lhw = (int64_t)li * H * W + (int64_t)hi * W + wi;

    #pragma unroll
    for (int ci = 0; ci < 12; ci++) {
        grids[base_n + (int64_t)ci * plane + lhw] = kIdent[ci];
    }
}


void bilagrid_affine_identity_init(
    float* grids, int N, int L, int H, int W, cudaStream_t stream
) {
    if (N <= 0 || L <= 0 || H <= 0 || W <= 0) return;
    dim3 block(16, 16, 1);
    dim3 bounds((W + 15) / 16, (H + 15) / 16, (int)(((int64_t)L * N + 0) / 1));
    bilagrid_affine_identity_init_kernel<<<bounds, block, 0, stream>>>(
        grids, N, L, H, W);
    CHECK_DEVICE_ERROR(cudaGetLastError());
}
