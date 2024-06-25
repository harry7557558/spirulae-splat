#define MAX_BLOCK_SIZE ( 16 * 16 )
#define N_THREADS 256

// kernel to use
// 0: Gaussian, 3 std
// 1: max(1-r^2, 0)
#define SPLAT_KERNEL 1

#define MAX_REGISTER_CHANNELS 3

// max number of cylindrical harmonics coefficients
#define MAX_CH_FLOAT3 (21)

// depth regularization
// 01: pairwise L1 with center depth
// 02: pairwise L2 with intersected depth
// 11: L1 intersected to reference
// 12: L2 intersected to reference
#define DEPTH_REG_L 11

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
