#pragma once

inline constexpr int WARP_SIZE = 32;

inline constexpr int TILE_SIZE = 16;

inline constexpr float ALPHA_THRESHOLD = (1.f/255.f);

#ifndef NO_TORCH
#endif

// #define CHECK_CUDA(x) TORCH_CHECK(x.is_cuda(), #x " must be a CUDA tensor")
// #define CHECK_HOST(x) TORCH_CHECK(!x.is_cuda(), #x " must be a CPU tensor")
// #define CHECK_CONTIGUOUS(x) \
//     TORCH_CHECK(x.is_contiguous(), #x " must be contiguous")
// #define CHECK_INPUT(x) \
//     do { CHECK_CUDA(x); CHECK_CONTIGUOUS(x); } while (0)
// #define DEVICE_GUARD(_ten) \
//     const at::cuda::OptionalCUDAGuard device_guard(device_of(_ten));

#define CHECK_CUDA(x)
#define CHECK_HOST(x)
#define CHECK_CONTIGUOUS(x)
#define CHECK_INPUT(x)
#define DEVICE_GUARD(_ten)

#define CHECK_DEVICE_ERROR(call)                                    \
do {                                                                \
    cudaError_t err = call;                                         \
    if (err != cudaSuccess) {                                       \
        fprintf(stderr, "\033[41mCUDA Error at %s:%d: %s\033[m\n",  \
                __FILE__, __LINE__, cudaGetErrorString(err));       \
        exit(EXIT_FAILURE);                                         \
    }                                                               \
} while (0)

#define _CEIL_DIV(n,m) (((n)+(m)-1)/(m))

#define _LAUNCH_ARGS_1D(n,b) ((n)==0?1:_CEIL_DIV(n,b)),b,0,(cudaStream_t)0
#define _LAUNCH_ARGS_2D(nx,ny,bx,by) dim3(_CEIL_DIV(nx,bx),_CEIL_DIV(ny,by),1),dim3(bx,by),0,(cudaStream_t)0
#define _LAUNCH_ARGS_3D(nx,ny,nz,bx,by,bz) dim3(_CEIL_DIV(nx,bx),_CEIL_DIV(ny,by),_CEIL_DIV(nz,bz)),dim3(bx,by,bz),0,(cudaStream_t)0

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


#include "common_utils.cuh"

// #ifndef NO_TORCH

// template<typename T, int ndim>
// TensorView<T, ndim> tensor2view(at::Tensor& tensor) {
//     TensorView<T, ndim> view;
//     view.data = tensor.data_ptr<T>();
//     for (int i = 0; i < ndim; i++) {
//         view.shape[i] = tensor.size(i);
//         view.strides[i] = *(tensor.strides().begin() + i);
//     }
//     return view;
// }

// #endif

// #ifndef NO_TORCH


// template<typename T>
// inline at::Tensor zeros_like(const at::Tensor& x) {
//     at::Tensor y = at::empty_like(x);
//     cudaMemset(y.data_ptr<T>(), 0, y.numel() * sizeof(T));
//     return y;
// }

// template<typename T>
// inline void set_zero(at::Tensor& x) {
//     cudaMemset(x.data_ptr<T>(), 0, x.numel() * sizeof(T));
// }

// at::Tensor zeros_like_tensor(const at::Tensor& x);

// void set_zero_tensor(at::Tensor& x);

// inline at::TensorOptions kTensorOptionF32()
//     { return at::TensorOptions(at::kCUDA).dtype(at::kFloat); }

// inline at::TensorOptions kTensorOptionI32()
//     { return at::TensorOptions(at::kCUDA).dtype(at::kInt); }

// inline at::TensorOptions kTensorOptionI64()
//     { return at::TensorOptions(at::kCUDA).dtype(at::kLong); }

// inline at::TensorOptions kTensorOptionByte()
//     { return at::TensorOptions(at::kCUDA).dtype(at::kByte); }

// inline at::TensorOptions kTensorOptionBool()
//     { return at::TensorOptions(at::kCUDA).dtype(at::kBool); }

// #endif
