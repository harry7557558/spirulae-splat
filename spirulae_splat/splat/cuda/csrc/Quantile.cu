// from https://raw.githubusercontent.com/harry7557558/fused-bilagrid/refs/heads/dev/fused_bilagrid/quantile.cu

#include <cuda_runtime.h>
#include <stdint.h>
#include <cmath>
#include <algorithm>

namespace {

constexpr int kThreadsPerBlock = 256;
constexpr int kTileElems = 4096;   // Tune this if you want more/fewer blocks.
constexpr uint32_t kPosInfBits = 0x7f800000u;

__device__ __forceinline__ bool positive_finite_bits(uint32_t u) {
    // Strictly positive and finite:
    //   +0 excluded  -> u > 0
    //   +inf / NaN excluded -> u < 0x7f800000
    // This also excludes all negative values.
    return (u > 0u) && (u < kPosInfBits);
}

__device__ __forceinline__ uint32_t block_reduce_sum_u32(uint32_t v) {
    __shared__ uint32_t s[kThreadsPerBlock];
    s[threadIdx.x] = v;
    __syncthreads();

    for (int offset = kThreadsPerBlock >> 1; offset > 0; offset >>= 1) {
        if (threadIdx.x < offset) {
            s[threadIdx.x] += s[threadIdx.x + offset];
        }
        __syncthreads();
    }
    return s[0];
}

__device__ __forceinline__ uint32_t block_reduce_max_u32(uint32_t v) {
    __shared__ uint32_t s[kThreadsPerBlock];
    s[threadIdx.x] = v;
    __syncthreads();

    for (int offset = kThreadsPerBlock >> 1; offset > 0; offset >>= 1) {
        if (threadIdx.x < offset) {
            s[threadIdx.x] = (s[threadIdx.x] > s[threadIdx.x + offset])
                               ? s[threadIdx.x]
                               : s[threadIdx.x + offset];
        }
        __syncthreads();
    }
    return s[0];
}

__global__ void find_row_max_kernel(
    const float* __restrict__ x,   // [B*N], row-major
    int B,
    int N,
    int tiles_per_row,
    uint32_t* __restrict__ row_max_bits // [B], init to 0
) {
    const int tile_id = static_cast<int>(blockIdx.x);
    const int row = tile_id / tiles_per_row;
    if (row >= B) return;

    const int tile_in_row = tile_id - row * tiles_per_row;
    const int start = tile_in_row * kTileElems;
    const int end = (start + kTileElems < N) ? (start + kTileElems) : N;

    const float* row_ptr = x + static_cast<size_t>(row) * static_cast<size_t>(N);

    uint32_t local_max = 0u;
    for (int i = start + threadIdx.x; i < end; i += blockDim.x) {
        uint32_t u = __float_as_uint(fabsf(row_ptr[i]));
        if (positive_finite_bits(u) && u > local_max) {
            local_max = u;
        }
    }

    uint32_t block_max = block_reduce_max_u32(local_max);
    if (threadIdx.x == 0 && block_max != 0u) {
        atomicMax(&row_max_bits[row], block_max);
    }
}

__global__ void hist_pass_kernel(
    const float* __restrict__ x,         // [B*N]
    int B,
    int N,
    int tiles_per_row,
    int pass,                            // 0..3
    const uint32_t* __restrict__ row_max_bits,   // [B]
    const uint32_t* __restrict__ row_prefix,      // [B], prefix selected so far
    const uint32_t* __restrict__ row_done,        // [B], 1 => skip
    uint32_t* __restrict__ row_max_count,         // [B], only used/updated in pass 0
    uint32_t* __restrict__ row_hist                // [B*256], zeroed before each pass
) {
    __shared__ uint32_t shist[256];

    for (int b = threadIdx.x; b < 256; b += blockDim.x) shist[b] = 0u;
    __syncthreads();

    const int tile_id = static_cast<int>(blockIdx.x);
    const int row = tile_id / tiles_per_row;
    if (row >= B) return;
    if (row_done[row]) return;

    const int tile_in_row = tile_id - row * tiles_per_row;
    const int start = tile_in_row * kTileElems;
    const int end = (start + kTileElems < N) ? (start + kTileElems) : N;

    const float* row_ptr = x + static_cast<size_t>(row) * static_cast<size_t>(N);

    const uint32_t max_bits = row_max_bits[row];

    if (pass == 0) {
        // Count all strictly positive finite values by top byte.
        // Also count exact max occurrences so we can exclude them later.
        uint32_t local_max_count = 0u;

        for (int i = start + threadIdx.x; i < end; i += blockDim.x) {
            uint32_t u = __float_as_uint(fabsf(row_ptr[i]));
            if (!positive_finite_bits(u)) continue;

            atomicAdd(&shist[(u >> 24) & 255u], 1u);
            if (u == max_bits) ++local_max_count;
        }

        uint32_t block_max_count = block_reduce_sum_u32(local_max_count);
        if (threadIdx.x == 0 && block_max_count != 0u) {
            atomicAdd(&row_max_count[row], block_max_count);
        }
    } else {
        const uint32_t prefix = row_prefix[row];
        const uint32_t mask = 0xFFFFFFFFu << (32 - 8 * pass);
        const uint32_t shift = 24 - 8 * pass;

        for (int i = start + threadIdx.x; i < end; i += blockDim.x) {
            uint32_t u = __float_as_uint(fabsf(row_ptr[i]));
            if (!positive_finite_bits(u)) continue;

            // Only keep candidates that match the already chosen high bytes.
            if ((u & mask) == prefix) {
                atomicAdd(&shist[(u >> shift) & 255u], 1u);
            }
        }
    }

    __syncthreads();

    for (int b = threadIdx.x; b < 256; b += blockDim.x) {
        atomicAdd(&row_hist[static_cast<size_t>(row) * 256u + static_cast<size_t>(b)], shist[b]);
    }
}

template<bool invert_quantile>
__global__ void select_pass_kernel(
    int B,
    int pass,                            // 0..3
    float q,
    const uint32_t* __restrict__ row_max_bits,   // [B]
    const uint32_t* __restrict__ row_max_count,  // [B]
    uint32_t* __restrict__ row_prefix,            // [B]
    uint32_t* __restrict__ row_rank,              // [B]
    uint32_t* __restrict__ row_done,              // [B]
    const uint32_t* __restrict__ row_hist,        // [B*256]
    float* __restrict__ out                       // [B]
) {
    const int row = blockIdx.x * blockDim.x + threadIdx.x;
    if (row >= B) return;
    if (row_done[row]) return;

    const uint32_t max_bits = row_max_bits[row];
    const uint32_t max_count = row_max_count[row];

    if (max_bits == 0u) {
        // No strictly positive finite values exist.
        if (threadIdx.x == 0) {} // no-op to keep style consistent
        out[row] = 1.0f;
        row_done[row] = 1u;
        return;
    }

    const uint32_t* h = row_hist + static_cast<size_t>(row) * 256u;

    if (pass == 0) {
        uint64_t total = 0;
        for (int b = 0; b < 256; ++b) total += h[b];

        const uint32_t masked_count = (total > max_count) ? static_cast<uint32_t>(total - max_count) : 0u;
        if (masked_count == 0u) {
            out[row] = invert_quantile
                ? 2.0f / __uint_as_float(max_bits)
                : __uint_as_float(max_bits) * 0.5f;
            row_done[row] = 1u;
            return;
        }

        // Lower-rank quantile convention.
        const uint32_t rank = (uint32_t)(fminf(1.0f, fmaxf(0.0f, q)) * (float)(masked_count - 1u));

        row_rank[row] = rank;
        row_prefix[row] = 0u;
    }

    uint32_t prefix = row_prefix[row];
    uint32_t rank = row_rank[row];

    const uint32_t shift = 24 - 8 * pass;
    const uint32_t mask_prefix = (pass == 0) ? 0u : (0xFFFFFFFFu << (32 - 8 * pass));
    const bool exclude_max = (pass == 0) ? true : ((prefix & mask_prefix) == (max_bits & mask_prefix));
    const uint32_t max_bucket = (max_bits >> shift) & 255u;

    uint32_t cum = 0u;
    uint32_t chosen = 255u;

    for (uint32_t b = 0; b < 256u; ++b) {
        uint32_t c = h[b];
        if (exclude_max && b == max_bucket) {
            c -= max_count;
        }
        if (rank < cum + c) {
            chosen = b;
            break;
        }
        cum += c;
    }

    if (pass == 3) {
        prefix |= chosen;
        out[row] = invert_quantile
            ? 1.0f / __uint_as_float(prefix)
            : __uint_as_float(prefix);
        row_done[row] = 1u;
    } else {
        prefix |= (chosen << shift);
        row_prefix[row] = prefix;
        row_rank[row] = rank - cum;
    }
}

struct QuantileWorkspace {
    uint32_t* d_row_max_bits = nullptr;   // [B]
    uint32_t* d_row_max_count = nullptr;  // [B]
    uint32_t* d_row_prefix = nullptr;     // [B]
    uint32_t* d_row_rank = nullptr;       // [B]
    uint32_t* d_row_done = nullptr;       // [B]
    uint32_t* d_row_hist = nullptr;       // [B*256]
    int B = 0;

    void allocate(uint32_t B_, uint32_t* buffer) {
        B = B_;
        d_row_max_bits = &buffer[0*B];
        d_row_max_count = &buffer[1*B];
        d_row_prefix = &buffer[2*B];
        d_row_rank = &buffer[3*B];
        d_row_done = &buffer[4*B];
        d_row_hist = &buffer[5*B];
    }

    void release() {
        d_row_hist = nullptr;
        d_row_done = nullptr;
        d_row_rank = nullptr;
        d_row_prefix = nullptr;
        d_row_max_count = nullptr;
        d_row_max_bits = nullptr;
        B = 0;
    }
};

} // namespace

template<bool invert_quantile>
cudaError_t batch_quantile_masked_radix_select(
    const float* d_x,   // [B*N], row-major
    int B,
    int N,
    float q,
    float* d_out,       // [B]
    uint32_t* temp,  // [(256+5)*B]
    cudaStream_t stream
) {
    if (!d_x || !d_out || B <= 0 || N <= 0) return cudaErrorInvalidValue;

    const int tiles_per_row = (N + kTileElems - 1) / kTileElems;
    const size_t hist_bytes = static_cast<size_t>(B) * 256u * sizeof(uint32_t);
    const size_t row_bytes  = static_cast<size_t>(B) * sizeof(uint32_t);

    QuantileWorkspace ws;
    ws.allocate(B, temp);

    // Init state.
    cudaError_t err;
    err = cudaMemsetAsync(ws.d_row_max_bits,  0, row_bytes, stream); if (err != cudaSuccess) { ws.release(); return err; }
    err = cudaMemsetAsync(ws.d_row_max_count, 0, row_bytes, stream); if (err != cudaSuccess) { ws.release(); return err; }
    err = cudaMemsetAsync(ws.d_row_prefix,    0, row_bytes, stream); if (err != cudaSuccess) { ws.release(); return err; }
    err = cudaMemsetAsync(ws.d_row_rank,      0, row_bytes, stream); if (err != cudaSuccess) { ws.release(); return err; }
    err = cudaMemsetAsync(ws.d_row_done,      0, row_bytes, stream); if (err != cudaSuccess) { ws.release(); return err; }

    // 1) Find the per-row maximum positive finite value.
    {
        const int total_blocks = B * tiles_per_row;
        find_row_max_kernel<<<total_blocks, kThreadsPerBlock, 0, stream>>>(
            d_x, B, N, tiles_per_row, ws.d_row_max_bits
        );
        err = cudaPeekAtLastError();
        if (err != cudaSuccess) { ws.release(); return err; }
    }

    // 2) Radix passes.
    for (int pass = 0; pass < 4; ++pass) {
        err = cudaMemsetAsync(ws.d_row_hist, 0, hist_bytes, stream);
        if (err != cudaSuccess) { ws.release(); return err; }

        if (pass == 0) {
            err = cudaMemsetAsync(ws.d_row_max_count, 0, row_bytes, stream);
            if (err != cudaSuccess) { ws.release(); return err; }
        }

        {
            const int total_blocks = B * tiles_per_row;
            hist_pass_kernel<<<total_blocks, kThreadsPerBlock, 0, stream>>>(
                d_x, B, N, tiles_per_row, pass,
                ws.d_row_max_bits,
                ws.d_row_prefix,
                ws.d_row_done,
                ws.d_row_max_count,
                ws.d_row_hist
            );
            err = cudaPeekAtLastError();
            if (err != cudaSuccess) { ws.release(); return err; }
        }

        {
            const int threads = 256;
            const int blocks = (B + threads - 1) / threads;
            select_pass_kernel<invert_quantile><<<blocks, threads, 0, stream>>>(
                B, pass, q,
                ws.d_row_max_bits,
                ws.d_row_max_count,
                ws.d_row_prefix,
                ws.d_row_rank,
                ws.d_row_done,
                ws.d_row_hist,
                d_out
            );
            err = cudaPeekAtLastError();
            if (err != cudaSuccess) { ws.release(); return err; }
        }
    }

    ws.release();
    return cudaSuccess;
}

template cudaError_t batch_quantile_masked_radix_select<false>(
    const float* d_x,   // [B*N], row-major
    int B,
    int N,
    float q,
    float* d_out,       // [B]
    uint32_t* temp,  // [(256+5)*B]
    cudaStream_t stream
);

template cudaError_t batch_quantile_masked_radix_select<true>(
    const float* d_x,   // [B*N], row-major
    int B,
    int N,
    float q,
    float* d_out,       // [B]
    uint32_t* temp,  // [(256+5)*B]
    cudaStream_t stream
);
