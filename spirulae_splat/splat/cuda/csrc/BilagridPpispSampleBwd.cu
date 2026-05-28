#include "BilagridPpispSampleBwd_kernel.cuh"

#define COMPUTE_COORDS_GRAD
#include "BilagridPpispSampleBwd_kernel.cuh"

#undef COMPUTE_COORDS_GRAD

#define PACKED

#include "BilagridPpispSampleBwd_kernel.cuh"

#define COMPUTE_COORDS_GRAD
#include "BilagridPpispSampleBwd_kernel.cuh"


void bilagrid_ppisp_sample_backward(
    const float* bilagrid,
    const float* coords,
    const float* rgb,
    const float* v_output,
    float* v_bilagrid,
    float* v_coords,
    float* v_rgb,
    int N, int L, int H, int W,
    int m, int h, int w,
    cudaStream_t stream
) {
    dim3 block = { 16, 16, 1 };
    dim3 bounds = {
        (w +block.x-1)/block.x,
        (h +block.y-1)/block.y,
        (N*m +block.z-1)/block.z
    };
    if (v_coords == nullptr) {
        bilagrid_ppisp_sample_backward_kernel<<<bounds, block, 0, stream>>>(
            bilagrid, coords, rgb, v_output,
            v_bilagrid, v_rgb,
            N, L, H, W, m, h, w
        );
    }
    else {
        bilagrid_ppisp_sample_backward_kernel_cg<<<bounds, block, 0, stream>>>(
            bilagrid, coords, rgb, v_output,
            v_bilagrid, v_coords, v_rgb,
            N, L, H, W, m, h, w
        );
    }
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


void bilagrid_ppisp_packed_sample_backward(
    const float* bilagrid,
    const int64_t* image_indices,
    const float* coords,
    const float* rgb,
    const float* v_output,
    float* v_bilagrid,
    float* v_coords,
    float* v_rgb,
    int N, int L, int H, int W,
    int nnz,
    cudaStream_t stream
) {
    int block = 256;
    int bounds = (nnz +block-1)/block;
    if (v_coords == nullptr) {
        bilagrid_ppisp_packed_sample_backward_kernel<<<bounds, block, 0, stream>>>(
            bilagrid, image_indices, coords, rgb, v_output,
            v_bilagrid, v_rgb,
            N, L, H, W, nnz
        );
    }
    else {
        bilagrid_ppisp_packed_sample_backward_kernel_cg<<<bounds, block, 0, stream>>>(
            bilagrid, image_indices, coords, rgb, v_output,
            v_bilagrid, v_coords, v_rgb,
            N, L, H, W, nnz
        );
    }
    CHECK_DEVICE_ERROR(cudaGetLastError());
}
