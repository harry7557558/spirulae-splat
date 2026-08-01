#include "nn/core/Error.h"

#include <algorithm>
#include "nn/Ops.h"
#include "nn/vk/Stream.h"

namespace nn {

namespace {

// Mirrors GemmParams in shaders/gemm.slang. Pointers first (8-byte aligned at
// natural C offsets, which is what Slang emits for a uniform block of raw
// pointers), then scalars.
struct GemmParams {
    uint64_t out;
    uint64_t x;
    uint64_t w;
    uint64_t bias;
    uint64_t residual;
    uint32_t M, N, K;
    float    alpha;
    uint32_t x_row_stride;
    uint32_t groups_per_row;
};
static_assert(sizeof(GemmParams) == 64, "GemmParams layout");

// Below this many rows the 64x64 tiling has nothing to chew on and a
// thread-per-K reduction wins. SAM 3 hits this constantly: every MLP head runs
// a single token, the hypernetwork mask heads run one each.
constexpr int64_t kThinRowLimit = 4;

// At or above this on both axes the 128x128 tile wins: it does four times the
// work per shared-memory read, and the edge waste is under 2% at SAM 3's
// shapes. Below it the wasted quarter-tiles cost more than the ratio buys.
constexpr int64_t kBigTileLimit = 128;

void dispatch_gemm(const Tensor& out, const Tensor& x_in, const Tensor& w_in,
                   const Tensor& bias, const Tensor& residual, Act act, float alpha,
                   int64_t x_row_stride) {
    // Weights keep their PyTorch rank in the store; a [Cout, Cin, kh, kw] conv
    // kernel is the same bytes as the [Cout, Cin*kh*kw] matrix this wants.
    const Tensor w = w_in.asMatrix();
    const Tensor x = x_in;
    const int64_t M = out.rows();
    const int64_t N = out.cols();
    const int64_t K = w.cols();

    NN_CHECK(out.dtype == DType::F32, "gemm output must be f32");
    NN_CHECK(w.rows() == N, "gemm: weight is [%lld, %lld] but out has %lld columns",
               (long long)w.rows(), (long long)K, (long long)N);
    NN_CHECK(x.cols() == K || x_row_stride > 0,
               "gemm: x has %lld columns but the weight expects %lld",
               (long long)x.cols(), (long long)K);
    NN_CHECK(x.rows() == M || x_row_stride > 0,
               "gemm: x has %lld rows but out has %lld", (long long)x.rows(),
               (long long)M);
    if (residual.valid())
        NN_CHECK(residual.numel() >= M * N, "gemm: residual is too small");
    if (bias.valid())
        NN_CHECK(bias.numel() >= N, "gemm: bias has %lld entries, need %lld",
                   (long long)bias.numel(), (long long)N);

    GemmParams p{};
    p.out = out.ptr;
    p.x = x.ptr;
    p.w = w.ptr;
    // Never leave an optional pointer null: a CPU Vulkan implementation can
    // speculate the load past its guard and fault (AGENTS.md).
    p.bias = vk::or_fallback(bias.ptr);
    p.residual = vk::or_fallback(residual.ptr);
    p.M = (uint32_t)M;
    p.N = (uint32_t)N;
    p.K = (uint32_t)K;
    p.alpha = alpha;
    p.x_row_stride = (uint32_t)(x_row_stride > 0 ? x_row_stride : K);

    vk::SpecList spec{(uint32_t)(w.dtype == DType::F16),
                      (uint32_t)(x.dtype == DType::F16),
                      (uint32_t)act,
                      (uint32_t)(bias.valid() ? 1 : 0),
                      (uint32_t)(residual.valid() ? 1 : 0)};

    if (M <= kThinRowLimit) {
        // One workgroup per output element, so N alone can exceed the 65535
        // grid cap; fold it across x and y and give the row to z.
        const uint32_t per_row = (uint32_t)std::min<int64_t>(N, 65535);
        const uint32_t rows = (uint32_t)((N + per_row - 1) / per_row);
        p.groups_per_row = per_row;
        vk::Stream::get().dispatch("gemm.gemm_nt_thin", spec, per_row, rows, (uint32_t)M, &p,
                                   sizeof(p));
    } else {
        // gemm_nt_big keeps the weight tile packed as fp16 in shared, so it has
        // nothing to offer an fp32 second operand (a matmul against another
        // activation) and does not handle one.
        const bool big =
            M >= kBigTileLimit && N >= kBigTileLimit && w.dtype == DType::F16;
        const int64_t tile = big ? 128 : 64;
        const uint32_t tm = (uint32_t)((M + tile - 1) / tile);
        const uint32_t tn = (uint32_t)((N + tile - 1) / tile);
        NN_CHECK(tm <= 65535 && tn <= 65535,
                   "gemm: %lldx%lld output exceeds the dispatch grid cap", (long long)M,
                   (long long)N);
        // gemm_nt_big wants N on x (see the comment there); gemm_nt wants M.
        vk::Stream::get().dispatch(big ? "gemm.gemm_nt_big" : "gemm.gemm_nt", spec,
                                   big ? tn : tm, big ? tm : tn, 1, &p, sizeof(p));
    }
}

}  // namespace

void linear(const Tensor& out, const Tensor& x, const Tensor& w, const LinearOpts& o) {
    dispatch_gemm(out, x, w, o.bias, o.residual, o.act, o.alpha, o.x_row_stride);
}

void matmul_nt(const Tensor& out, const Tensor& a, const Tensor& b, float alpha,
               Act act) {
    dispatch_gemm(out, a, b, {}, {}, act, alpha, 0);
}

}  // namespace nn
