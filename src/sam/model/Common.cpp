// Attention and MLP helpers shared by every transformer in the model.

#include "sam/Common.h"
#include "nn/core/Error.h"
#include "sam/model/Modules.h"

#include <algorithm>
#include <cmath>

namespace sam {
namespace model {

using nn::Act;
using nn::AttnBias;
using nn::Tensor;

void mha_fused(vk::Arena& arena, const Tensor& out, const Tensor& q_src,
               const Tensor& k_src, const Tensor& v_src, int64_t nq, int64_t nkv,
               const Tensor& in_w, const Tensor& in_b, const Tensor& out_w,
               const Tensor& out_b, int n_heads, AttnBias bias_mode, const Tensor& bias,
               const Tensor& residual) {
    const int64_t D = out_w.cols();
    NN_CHECK(in_w.rows() == 3 * D,
               "mha_fused: in_proj is [%lld, %lld]; expected 3x%lld rows",
               (long long)in_w.rows(), (long long)in_w.cols(), (long long)D);
    NN_CHECK(D % n_heads == 0, "mha_fused: %lld dims do not split into %d heads",
               (long long)D, n_heads);

    vk::ArenaScope scope(arena);
    Tensor q = nn::arena_tensor(arena, nn::DType::F32, nq, D);
    Tensor k = nn::arena_tensor(arena, nn::DType::F32, nkv, D);
    Tensor v = nn::arena_tensor(arena, nn::DType::F32, nkv, D);

    // nn.MultiheadAttention stacks the three projections into one [3D, D]
    // matrix; slicing rows gives each one back with no copy.
    auto project = [&](const Tensor& dst, const Tensor& src, int64_t which, int64_t n) {
        nn::LinearOpts lo;
        lo.bias = in_b.slice0(which * D, D);
        lo.x_row_stride = src.cols();
        nn::linear(dst, src.view(n, src.cols()), in_w.slice0(which * D, D), lo);
    };
    project(q, q_src, 0, nq);
    project(k, k_src, 1, nkv);
    project(v, v_src, 2, nkv);

    Tensor ctx = nn::arena_tensor(arena, nn::DType::F32, nq, D);
    nn::AttnOpts ao;
    ao.arena = &arena;
    ao.n_heads = n_heads;
    ao.head_dim = (int)(D / n_heads);
    ao.bias_mode = bias_mode;
    ao.bias = bias;
    nn::attention(ctx, q, k, v, nq, nkv, ao);

    nn::LinearOpts lo;
    lo.bias = out_b;
    lo.residual = residual;
    nn::linear(out, ctx, out_w, lo);
}

void mha_split(vk::Arena& arena, const Tensor& out, const Tensor& q_src,
               const Tensor& k_src, const Tensor& v_src, int64_t nq, int64_t nkv,
               const Tensor& qw, const Tensor& qb, const Tensor& kw, const Tensor& kb,
               const Tensor& vw, const Tensor& vb, const Tensor& ow, const Tensor& ob,
               int n_heads, const Tensor& residual) {
    // SAM's Attention takes embedding_dim in and projects to internal_dim,
    // which is embedding_dim / downsample_rate for the cross-attentions.
    const int64_t ID = qw.rows();
    NN_CHECK(ID % n_heads == 0, "mha_split: internal dim %lld / %d heads",
               (long long)ID, n_heads);

    vk::ArenaScope scope(arena);
    Tensor q = nn::arena_tensor(arena, nn::DType::F32, nq, ID);
    Tensor k = nn::arena_tensor(arena, nn::DType::F32, nkv, ID);
    Tensor v = nn::arena_tensor(arena, nn::DType::F32, nkv, ID);

    auto project = [&](const Tensor& dst, const Tensor& src, const Tensor& w,
                       const Tensor& b, int64_t n) {
        nn::LinearOpts lo;
        lo.bias = b;
        lo.x_row_stride = src.cols();
        nn::linear(dst, src.view(n, src.cols()), w, lo);
    };
    project(q, q_src, qw, qb, nq);
    project(k, k_src, kw, kb, nkv);
    project(v, v_src, vw, vb, nkv);

    Tensor ctx = nn::arena_tensor(arena, nn::DType::F32, nq, ID);
    nn::AttnOpts ao;
    ao.arena = &arena;
    ao.n_heads = n_heads;
    ao.head_dim = (int)(ID / n_heads);
    nn::attention(ctx, q, k, v, nq, nkv, ao);

    nn::LinearOpts lo;
    lo.bias = ob;
    lo.residual = residual;
    nn::linear(out, ctx, ow, lo);
}

void mlp_relu(vk::Arena& arena, const Tensor& out, const Tensor& x, const SamModel& m,
              const char* prefix, int n_layers, Act final_act) {
    vk::ArenaScope scope(arena);
    Tensor cur = x;
    const int64_t rows = out.rows();
    for (int i = 0; i < n_layers; ++i) {
        Tensor w = m.w().getf((std::string(prefix) + ".%d.weight").c_str(), i);
        Tensor b = m.w().getf((std::string(prefix) + ".%d.bias").c_str(), i);
        const bool last = (i == n_layers - 1);
        Tensor dst = last ? out
                          : nn::arena_tensor(arena, nn::DType::F32, rows, w.rows());
        nn::LinearOpts lo;
        lo.bias = b;
        lo.act = last ? final_act : Act::Relu;
        nn::linear(dst, cur, w, lo);
        cur = dst;
    }
}

}  // namespace model
}  // namespace sam
