// The SAM mask decoder: a two-way transformer between output tokens and image
// features, then a learned upscaling to quarter resolution and a hypernetwork
// that turns each mask token into a 32-channel filter.
//
// Shared by the PVS path (interactive points/boxes) and by video propagation --
// in the latter case `image_embed` is the memory-conditioned feature map rather
// than the raw neck output.

#include "sam/Common.h"
#include "nn/core/Error.h"
#include "nn/core/Log.h"
#include "sam/model/Modules.h"

#include <cmath>

namespace sam {
namespace model {

using nn::Act;
using nn::Tensor;

namespace {

// TwoWayAttentionBlock: tokens attend to themselves, then to the image, an MLP
// on the tokens, then the image attends back to the tokens.
void two_way_block(const SamModel& m, vk::Arena& arena, const Tensor& queries,
                   const Tensor& keys, const Tensor& query_pe, const Tensor& key_pe,
                   int64_t nq, int64_t nkv, int block_idx, bool skip_first_layer_pe) {
    const WeightStore& W = m.w();
    const int64_t D = m.hp().sam_embed_dim;
    const int n_heads = 8;

    vk::ArenaScope scope(arena);
    char p[48];
    std::snprintf(p, sizeof p, "sam_dec.twoway.%d", block_idx);
    auto wt = [&](const char* s) { return W.getf("%s.%s", p, s); };
    auto attn = [&](const char* pfx, const Tensor& out, const Tensor& q, const Tensor& k,
                    const Tensor& v, int64_t n_q, int64_t n_kv, const Tensor& residual) {
        mha_split(arena, out, q, k, v, n_q, n_kv, W.getf("%s.%s.q_proj.weight", p, pfx),
                  W.getf("%s.%s.q_proj.bias", p, pfx),
                  W.getf("%s.%s.k_proj.weight", p, pfx),
                  W.getf("%s.%s.k_proj.bias", p, pfx),
                  W.getf("%s.%s.v_proj.weight", p, pfx),
                  W.getf("%s.%s.v_proj.bias", p, pfx),
                  W.getf("%s.%s.out_proj.weight", p, pfx),
                  W.getf("%s.%s.out_proj.bias", p, pfx), n_heads, residual);
    };

    Tensor qp = nn::arena_tensor(arena, nn::DType::F32, nq, D);
    Tensor kp = nn::arena_tensor(arena, nn::DType::F32, nkv, D);

    // 1. token self-attention. On the first block the reference skips the
    // positional term AND the residual -- the queries at that point ARE the
    // positional encoding, so adding them twice would double-count.
    if (skip_first_layer_pe) {
        attn("sa", queries, queries, queries, queries, nq, nq, {});
    } else {
        nn::add(qp, queries, query_pe);
        attn("sa", queries, qp, qp, queries, nq, nq, queries);
    }
    nn::layer_norm(queries, queries, wt("norm1.weight"), wt("norm1.bias"));

    // 2. tokens -> image
    nn::add(qp, queries, query_pe);
    nn::add(kp, keys, key_pe);
    attn("cross_attn_token_to_image", queries, qp, kp, keys, nq, nkv, queries);
    nn::layer_norm(queries, queries, wt("norm2.weight"), wt("norm2.bias"));

    // 3. token MLP
    {
        Tensor hid = nn::arena_tensor(arena, nn::DType::F32, nq, wt("mlp.lin1.weight").rows());
        nn::LinearOpts l1;
        l1.bias = wt("mlp.lin1.bias");
        l1.act = Act::Relu;
        nn::linear(hid, queries, wt("mlp.lin1.weight"), l1);
        nn::LinearOpts l2;
        l2.bias = wt("mlp.lin2.bias");
        l2.residual = queries;
        nn::linear(queries, hid, wt("mlp.lin2.weight"), l2);
        nn::layer_norm(queries, queries, wt("norm3.weight"), wt("norm3.bias"));
    }

    // 4. image -> tokens (q and k swap roles)
    nn::add(qp, queries, query_pe);
    nn::add(kp, keys, key_pe);
    attn("cross_attn_image_to_token", keys, kp, qp, queries, nkv, nq, keys);
    nn::layer_norm(keys, keys, wt("norm4.weight"), wt("norm4.bias"));
}

}  // namespace

void sam_mask_decoder(const SamModel& m, vk::Arena& arena, const Tensor& image_embed,
                      const Tensor& sparse, int64_t n_sparse, const Tensor& feat_s0,
                      const Tensor& feat_s1, SamDecoderOutputs& out) {
    NN_SCOPED_TIMER("sam_mask_decoder");
    const Hparams& h = m.hp();
    const WeightStore& W = m.w();
    const int64_t D = h.sam_embed_dim;
    const int64_t G = m.grid();
    const int64_t NKV = G * G;
    const int64_t NM = h.num_mask_tokens();
    const int64_t n_special = 2 + NM;  // obj score, IoU, then the mask tokens
    const int64_t NQ = n_special + n_sparse;
    const int64_t S = G * 4;

    // Outputs are allocated BEFORE the scope: everything below is scratch that
    // the rewind reclaims, but these belong to the caller.
    out.masks = nn::arena_tensor(arena, nn::DType::F32, NM, S * S);
    out.iou = nn::arena_tensor(arena, nn::DType::F32, NM);
    out.obj_score = nn::arena_tensor(arena, nn::DType::F32, 1);
    out.mask_tokens = nn::arena_tensor(arena, nn::DType::F32, NM, D);

    vk::ArenaScope scope(arena);

    // ---- output tokens + the sparse prompt ----
    Tensor tokens = nn::arena_tensor(arena, nn::DType::F32, NQ, D);
    nn::strided_copy(tokens.slice0(0, 1), W.get("sam_dec.obj_score_token.weight"), 1, D, D, D);
    nn::strided_copy(tokens.slice0(1, 1), W.get("sam_dec.iou_token.weight"), 1, D, D, D);
    nn::strided_copy(tokens.slice0(2, NM), W.get("sam_dec.mask_tokens.weight"), NM, D, D, D);
    if (n_sparse > 0)
        nn::strided_copy(tokens.slice0(n_special, n_sparse), sparse, n_sparse, D, D, D);

    // The initial token values double as the query positional encoding
    // throughout the transformer, so they must be kept unmodified.
    Tensor query_pe = nn::arena_tensor(arena, nn::DType::F32, NQ, D);
    nn::strided_copy(query_pe, tokens, NQ, D, D, D);

    // ---- image side ----
    Tensor keys = nn::arena_tensor(arena, nn::DType::F32, NKV, D);
    nn::add(keys, image_embed.view(NKV, D), m.denseNoMask().view(NKV, D));
    Tensor key_pe = m.densePe().view(NKV, D);

    for (int i = 0; i < h.sam_dec_depth; ++i)
        two_way_block(m, arena, tokens, keys, query_pe, key_pe, NQ, NKV, i, i == 0);

    // ---- final tokens -> image attention ----
    {
        vk::ArenaScope inner(arena);
        Tensor qp = nn::arena_tensor(arena, nn::DType::F32, NQ, D);
        Tensor kp = nn::arena_tensor(arena, nn::DType::F32, NKV, D);
        nn::add(qp, tokens, query_pe);
        nn::add(kp, keys, key_pe);
        mha_split(arena, tokens, qp, kp, keys, NQ, NKV,
                  W.get("sam_dec.final_attn.q_proj.weight"),
                  W.get("sam_dec.final_attn.q_proj.bias"),
                  W.get("sam_dec.final_attn.k_proj.weight"),
                  W.get("sam_dec.final_attn.k_proj.bias"),
                  W.get("sam_dec.final_attn.v_proj.weight"),
                  W.get("sam_dec.final_attn.v_proj.bias"),
                  W.get("sam_dec.final_attn.out_proj.weight"),
                  W.get("sam_dec.final_attn.out_proj.bias"), 8, tokens);
        nn::layer_norm(tokens, tokens, W.get("sam_dec.final_norm.weight"),
                       W.get("sam_dec.final_norm.bias"));
    }

    Tensor obj_token = tokens.slice0(0, 1);
    Tensor iou_token = tokens.slice0(1, 1);
    nn::strided_copy(out.mask_tokens, tokens.slice0(2, NM), NM, D, D, D);

    // ---- upscale the attended image features to 4x, folding in the two
    // high-resolution FPN levels ----
    Tensor up2 = nn::arena_tensor(arena, nn::DType::F32, G * 4, G * 4, 32);
    {
        vk::ArenaScope inner(arena);
        Tensor up1 = nn::arena_tensor(arena, nn::DType::F32, G * 2, G * 2, 64);
        nn::conv_transpose2x2(arena, up1, keys.view(G, G, D), m.samUpscale(0),
                              W.get("sam_dec.upscale.0.bias"), Act::None);
        {
            // conv_s1 projects the 2x FPN level into the same 64 channels; the
            // reference adds it BEFORE the LayerNorm, not after.
            Tensor hs1 = nn::arena_tensor(arena, nn::DType::F32, G * 2 * G * 2, 64);
            nn::LinearOpts lo;
            lo.bias = W.get("sam_dec.conv_s1.bias");
            nn::linear(hs1, feat_s1.view(G * 2 * G * 2, D), W.get("sam_dec.conv_s1.weight"),
                       lo);
            nn::add(up1.view(G * 2 * G * 2, 64), up1.view(G * 2 * G * 2, 64), hs1);
        }
        nn::layer_norm(up1.view(G * 2 * G * 2, 64), up1.view(G * 2 * G * 2, 64),
                       W.get("sam_dec.upscale.1.weight"), W.get("sam_dec.upscale.1.bias"));
        nn::unary(up1.view(G * 2 * G * 2, 64), up1.view(G * 2 * G * 2, 64), Act::GeluErf);

        nn::conv_transpose2x2(arena, up2, up1, m.samUpscale(1),
                              W.get("sam_dec.upscale.3.bias"), Act::None);
        {
            Tensor hs0 = nn::arena_tensor(arena, nn::DType::F32, G * 4 * G * 4, 32);
            nn::LinearOpts lo;
            lo.bias = W.get("sam_dec.conv_s0.bias");
            nn::linear(hs0, feat_s0.view(G * 4 * G * 4, D), W.get("sam_dec.conv_s0.weight"),
                       lo);
            // No LayerNorm on this one -- act2(dc2(x) + feat_s0) in the reference.
            nn::add(up2.view(G * 4 * G * 4, 32), up2.view(G * 4 * G * 4, 32), hs0,
                    1.0f, 1.0f, Act::GeluErf);
        }
    }

    // ---- hypernetwork masks ----
    {
        vk::ArenaScope inner(arena);
        Tensor hyper = nn::arena_tensor(arena, nn::DType::F32, NM, 32);
        for (int64_t k = 0; k < NM; ++k) {
            char pfx[64];
            std::snprintf(pfx, sizeof pfx, "sam_dec.hyper.%d.layers", (int)k);
            mlp_relu(arena, hyper.slice0(k, 1), out.mask_tokens.slice0(k, 1), m, pfx, 3);
        }
        nn::matmul_nt(out.masks, hyper, up2.view(S * S, 32));
    }

    mlp_relu(arena, out.iou.view(1, NM), iou_token, m, "sam_dec.iou_prediction_head.layers",
             3, Act::Sigmoid);

    mlp_relu(arena, out.obj_score.view(1, 1), obj_token, m,
             "sam_dec.pred_obj_score_head.layers", 3);
}

}  // namespace model
}  // namespace sam
