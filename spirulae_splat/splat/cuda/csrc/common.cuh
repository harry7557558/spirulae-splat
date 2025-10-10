#pragma once

#define BLOCK_WIDTH 16
#define MAX_BLOCK_SIZE ( 16 * 16 )
#define N_THREADS 256


#include <c10/cuda/CUDAGuard.h>
#include <torch/types.h>

#define CHECK_CUDA(x) TORCH_CHECK(x.is_cuda(), #x " must be a CUDA tensor")
#define CHECK_CONTIGUOUS(x)                                                    \
    TORCH_CHECK(x.is_contiguous(), #x " must be contiguous")
#define CHECK_INPUT(x)                                                         \
    CHECK_CUDA(x);                                                             \
    CHECK_CONTIGUOUS(x)
#define DEVICE_GUARD(_ten) \
    const at::cuda::OptionalCUDAGuard device_guard(device_of(_ten));

#define CHECK_DEVICE_ERROR(call)                                      \
do {                                                                \
    cudaError_t err = call;                                         \
    if (err != cudaSuccess) {                                       \
        fprintf(stderr, "CUDA Error at %s:%d: %s\n",                \
                __FILE__, __LINE__, cudaGetErrorString(err));       \
        exit(EXIT_FAILURE);                                         \
    }                                                               \
} while (0)


#define _LAUNGH_ARGS_1D(n) (n+N_THREADS-1)/N_THREADS,N_THREADS

//--------------
#define CUDA_CALL(x)                                                           \
    do {                                                                       \
        if ((x) != cudaSuccess) {                                              \
            printf(                                                            \
                "Error at %s:%d - %s\n",                                       \
                __FILE__,                                                      \
                __LINE__,                                                      \
                cudaGetErrorString(cudaGetLastError())                         \
            );                                                                 \
            exit(EXIT_FAILURE);                                                \
        }                                                                      \
    } while (0)



#define BLOCK_DIM3 dim3(BLOCK_WIDTH, BLOCK_WIDTH, 1)


inline __host__ float4 tuple2float4(std::tuple<float, float, float, float> v) {
    return {std::get<0>(v), std::get<1>(v), std::get<2>(v), std::get<3>(v)};
}

inline __host__ dim3 tuple2dim3(std::tuple<unsigned, unsigned, unsigned> v) {
    return {std::get<0>(v), std::get<1>(v), std::get<2>(v)};
}

inline __host__ dim3 whb2tb(unsigned width, unsigned height, unsigned block_width=BLOCK_WIDTH) {
    return {
        (width + block_width - 1) / block_width,
        (height + block_width - 1) / block_width,
        1
    };
}

#include "common_utils.cuh"

template<typename T, int ndim>
TensorView<T, ndim> tensor2view(torch::Tensor& tensor) {
    TensorView<T, ndim> view;
    view.data = tensor.data_ptr<T>();
    for (int i = 0; i < ndim; i++) {
        view.shape[i] = tensor.size(i);
        view.strides[i] = *(tensor.strides().begin() + i);
    }
    return view;
}

