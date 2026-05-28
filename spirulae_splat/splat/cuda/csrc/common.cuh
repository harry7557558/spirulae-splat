#pragma once

inline constexpr int WARP_SIZE = 32;

inline constexpr int TILE_SIZE = 16;

inline constexpr float ALPHA_THRESHOLD = (1.f/255.f);

#ifndef NO_TORCH
#endif

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
