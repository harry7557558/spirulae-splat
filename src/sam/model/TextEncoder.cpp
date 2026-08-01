// CLIP-style text encoder: 24 causal pre-norm transformer blocks over a
// 32-token context, then a linear "resizer" down to the 256-dim model width.
//
// The reference keeps the full token sequence and throws the pooled output
// away, so `text_projection` is not in the checkpoint and is not used here --
// only the resizer output feeds the fusion encoder and the DETR decoder.

#include "sam/Common.h"
#include "nn/core/Error.h"
#include "nn/core/Log.h"
#include "sam/model/Modules.h"
#include "nn/vk/Stream.h"

namespace sam {
namespace model {

using nn::Act;
using nn::Tensor;

void encode_text(const SamModel& m, vk::Arena& arena,
                 const std::vector<int32_t>& token_ids, const Tensor& out) {
    NN_SCOPED_TIMER("encode_text");
    const Hparams& h = m.hp();
    const WeightStore& W = m.w();
    const int64_t L = h.text_ctx_len;
    const int64_t D = h.text_width;
    const int NH = h.text_heads;
    NN_CHECK((int64_t)token_ids.size() == L,
               "encode_text: %zu ids given, the context length is %lld",
               token_ids.size(), (long long)L);
    NN_CHECK(out.rows() == L && out.cols() == h.text_out_dim,
               "encode_text: output must be [%lld, %d]", (long long)L, h.text_out_dim);

    vk::ArenaScope scope(arena);

    Tensor ids = nn::arena_tensor(arena, nn::DType::I32, L);
    vk::Stream::get().upload(ids.ptr, token_ids.data(), token_ids.size() * 4);

    Tensor x = nn::arena_tensor(arena, nn::DType::F32, L, D);
    nn::gather_rows(x, W.get("text.token_embed.weight").asMatrix(), ids);
    nn::add(x, x, W.get("text.pos_embed"));

    for (int i = 0; i < h.text_layers; ++i) {
        vk::ArenaScope blk(arena);
        char p[48];
        std::snprintf(p, sizeof p, "text.blocks.%d", i);
        auto wt = [&](const char* s) { return W.getf("%s.%s", p, s); };

        Tensor xn = nn::arena_tensor(arena, nn::DType::F32, L, D);
        nn::layer_norm(xn, x, wt("ln_1.weight"), wt("ln_1.bias"));
        // Causal mask: the text encoder is autoregressive even though we only
        // ever use the whole sequence, and the trained weights depend on it.
        mha_fused(arena, x, xn, xn, xn, L, L, wt("attn.in_proj.weight"),
                  wt("attn.in_proj.bias"), wt("attn.out_proj.weight"),
                  wt("attn.out_proj.bias"), NH, nn::AttnBias::Causal, {}, x);

        nn::layer_norm(xn, x, wt("ln_2.weight"), wt("ln_2.bias"));
        Tensor hid = nn::arena_tensor(arena, nn::DType::F32, L, D * 4);
        nn::LinearOpts l1;
        l1.bias = wt("mlp.fc1.bias");
        l1.act = Act::GeluErf;
        nn::linear(hid, xn, wt("mlp.fc1.weight"), l1);
        nn::LinearOpts l2;
        l2.bias = wt("mlp.fc2.bias");
        l2.residual = x;
        nn::linear(x, hid, wt("mlp.fc2.weight"), l2);
    }

    nn::layer_norm(x, x, W.get("text.ln_final.weight"), W.get("text.ln_final.bias"));

    nn::LinearOpts lo;
    lo.bias = W.get("text.resizer.bias");
    nn::linear(out, x, W.get("text.resizer.weight"), lo);
}

}  // namespace model
}  // namespace sam
