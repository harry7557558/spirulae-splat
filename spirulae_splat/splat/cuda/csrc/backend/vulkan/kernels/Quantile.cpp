// Vulkan implementation of the batched masked quantile radix-select
// (csrc/Quantile.cu). Device work: slang/vulkan/quantile.slang. The selection
// is pure integer arithmetic (histograms + bit-image radix passes), so
// outputs are bit-exact against CUDA.
//
// Referenced by the portable engine layer (BilagridBindings.cpp, <true>) and
// by the multi-scale loss's RobustEdgeAware path (<false>).

#include <Tensor.h>

#include "KernelCommon.h"

namespace {

// Mirrors QuantileParams in slang/vulkan/quantile.slang.
struct QuantileParams {
    uint64_t x, row_max_bits, row_max_count, row_prefix, row_rank, row_done,
        row_hist, out;
    int32_t B, N, tiles_per_row, pass;
    float q;
    uint32_t wgs_per_row;
};
static_assert(sizeof(QuantileParams) == 8 * 8 + 6 * 4, "layout");

// Fold a 1D block range across (gx, gy) under the 65535 grid cap.
void dispatch_blocks(const char* entry, const backend::vk::SpecList& spec,
                     int64_t blocks, QuantileParams& p) {
    if (blocks <= 0) return;
    uint32_t per_row = (uint32_t)std::min<int64_t>(blocks, 65535);
    uint32_t rows = (uint32_t)((blocks + per_row - 1) / per_row);
    p.wgs_per_row = per_row;
    vkk::dispatch(entry, spec, per_row, rows, 1, &p, sizeof(p));
}

}  // namespace

template <bool invert_quantile>
int batch_quantile_masked_radix_select(
    const float* d_x,   // [B*N], row-major
    int B,
    int N,
    float q,
    float* d_out,       // [B]
    uint32_t* temp,     // [(256+5)*B]
    backend::Stream stream
) {
    (void)stream;
    if (!d_x || !d_out || B <= 0 || N <= 0) return 1;

    const int tiles_per_row = (N + 4096 - 1) / 4096;
    const size_t row_bytes = (size_t)B * sizeof(uint32_t);

    // Workspace layout mirrors Quantile.cu's QuantileWorkspace.
    uint32_t* d_row_max_bits = &temp[0 * B];
    uint32_t* d_row_max_count = &temp[1 * B];
    uint32_t* d_row_prefix = &temp[2 * B];
    uint32_t* d_row_rank = &temp[3 * B];
    uint32_t* d_row_done = &temp[4 * B];
    uint32_t* d_row_hist = &temp[5 * B];

    backend::memset_async(d_row_max_bits, 0, 5 * row_bytes,
                          backend::kDefaultStream);

    const backend::vk::SpecList spec{invert_quantile ? 1u : 0u};

    QuantileParams p{};
    p.x = (uint64_t)d_x;
    p.row_max_bits = (uint64_t)d_row_max_bits;
    p.row_max_count = (uint64_t)d_row_max_count;
    p.row_prefix = (uint64_t)d_row_prefix;
    p.row_rank = (uint64_t)d_row_rank;
    p.row_done = (uint64_t)d_row_done;
    p.row_hist = (uint64_t)d_row_hist;
    p.out = (uint64_t)d_out;
    p.B = B;
    p.N = N;
    p.tiles_per_row = tiles_per_row;
    p.q = q;

    const int64_t total_blocks = (int64_t)B * tiles_per_row;
    p.pass = 0;
    dispatch_blocks("quantile.qt_find_row_max", spec, total_blocks, p);

    for (int pass = 0; pass < 4; pass++) {
        backend::memset_async(d_row_hist, 0, (size_t)B * 256 * sizeof(uint32_t),
                              backend::kDefaultStream);
        if (pass == 0)
            backend::memset_async(d_row_max_count, 0, row_bytes,
                                  backend::kDefaultStream);
        p.pass = pass;
        dispatch_blocks("quantile.qt_hist_pass", spec, total_blocks, p);
        dispatch_blocks("quantile.qt_select_pass", spec, (B + 255) / 256, p);
    }
    return 0;
}

template int batch_quantile_masked_radix_select<false>(
    const float*, int, int, float, float*, uint32_t*, backend::Stream);
template int batch_quantile_masked_radix_select<true>(
    const float*, int, int, float, float*, uint32_t*, backend::Stream);
