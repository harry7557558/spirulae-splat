#include "kernels/bilagrid/BilagridPpispSampleFwd_kernel.cuh"

#define PACKED
#include "kernels/bilagrid/BilagridPpispSampleFwd_kernel.cuh"


void bilagrid_ppisp_sample_forward(
    BilagridReader bilagrid,
    const float* coords,
    const float* rgb,
    float* output,
    int N, int L, int H, int W,
    int m, int h, int w,
    cudaStream_t stream
) {
    int total = N * m * h * w;
    int threads = 256;
    int blocks = (total + threads - 1) / threads;
    bilagrid_ppisp_sample_forward_kernel<<<blocks, threads, 0, stream>>>(
        bilagrid, coords, rgb, output,
        N, L, H, W, m, h, w
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}

void bilagrid_ppisp_packed_sample_forward(
    BilagridReader bilagrid,
    const int64_t* image_indices,
    const float* coords,
    const float* rgb,
    float* output,
    int N, int L, int H, int W,
    int nnz,
    cudaStream_t stream
) {
    int threads = 256;
    int blocks = (nnz + threads - 1) / threads;
    bilagrid_ppisp_packed_sample_forward_kernel<<<blocks, threads, 0, stream>>>(
        bilagrid, image_indices, coords, rgb, output,
        N, L, H, W, nnz
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}
