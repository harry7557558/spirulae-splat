#pragma once

#define MAX_BLOCK_SIZE ( 16 * 16 )
#define N_THREADS 256
// for per-pixel sorting
#define N_THREADS_PPS 64

// kernel to use
// 0: Gaussian, 3 std
// 1: max(1-r^2, 0)
#define SPLAT_KERNEL 1

// sorted splats list size
#define MAX_SORTED_SPLATS 32
#define SORTED_INDEX_INF 0x7fffffff

// max number of cylindrical harmonics coefficients
#define MAX_CH_FLOAT3 (21)

// threshold for median depth
// a number less than 0.5 pushes splats away from camera
#define DEPTH_REG_MEDIAN_TH 0.5f

// ## DEPRECATED
// depth regularization
// 01: pairwise L1 with center depth
// 02: pairwise L2 with intersected depth
// 11: L1 intersected to reference
// 12: L2 intersected to reference
#define DEPTH_REG_L 0

// CH abs grad, 0 for mean, 1 for max
#define CH_ABSGRAD_REDUCE 0

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


// mean vs median depth
enum class DepthMode {
	Mean,
    Median,
};
