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

#endif  // __CUDACC__
