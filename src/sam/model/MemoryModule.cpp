// Video memory: the encoder that turns a predicted mask into a compact spatial
// memory, and the attention that conditions the next frame on the bank.
//
// A note on GELU: the ggml reference uses the tanh approximation here (and the
// erf form in the ViT), while the PyTorch source is `nn.GELU()` -- erf -- in
// both. We use erf throughout, which is the reference the weights were trained
// against; the two differ by under 0.003 absolute, so this is a small move
// toward the original, not away from it.

#include "sam/Common.h"
#include "nn/core/Error.h"
#include "nn/core/Log.h"
#include "sam/model/Modules.h"

#include <cmath>

namespace sam {
namespace model {

using nn::Act;
using nn::Tensor;

// ---------------------------------------------------------------------------
// Memory encoder
// ---------------------------------------------------------------------------

void encode_memory(const SamModel& m, vk::Arena& arena, const Tensor& mask_logits,
                   int64_t mh, int64_t mw, const Tensor& pix_feats, const Tensor& out) {
    const Hparams& h = m.hp();
    const WeightStore& W = m.w();
    const int64_t D = h.neck_dim;
    const int64_t MD = h.mem_out_dim;
    const int64_t G = m.grid();
    const int64_t hi_res = m.imgSize();
    const int64_t interp = G * 16;  // the mask downsampler's four stride-2 stages

    vk::ArenaScope scope(arena);

    // Mask preprocessing: up to full input resolution, sigmoid, scale by the
    // checkpoint's sigmoid_scale/bias (20 and -10 in every shipped config, so
    // into [-10, 10]), then to the downsampler's input size. The two-step
    // resize is the reference's; going straight to `interp` would sample a
    // different grid. SAM 2 has interp == img_size, so the second one is a
    // no-op and is skipped.
    // ndim is spelled out: a trailing size-1 channel axis would otherwise be
    // collapsed by arena_tensor's rank inference, and the resamplers require a
    // genuine [H, W, C] tensor.
    Tensor hires = nn::arena_tensor(arena, nn::DType::F32, hi_res, hi_res, 1, 1, 3);
    nn::resize_bilinear(hires, mask_logits.view(mh, mw, 1));
    nn::unary(hires, hires, Act::Sigmoid, 1.0f, 0.0f, h.mem_sigmoid_scale,
              h.mem_sigmoid_bias);

    Tensor cur = hires;
    if (interp != hi_res) {
        cur = nn::arena_tensor(arena, nn::DType::F32, interp, interp, 1, 1, 3);
        nn::resize_bilinear(cur, hires);
    }

    // Learned mask downsampler: 4 x (conv3x3 stride 2, LayerNorm2d, GELU), then
    // a 1x1. Channels 1 -> 4 -> 16 -> 64 -> 256, resolution /16.
    {
        const int ds_conv_idx[4] = {0, 3, 6, 9};
        const int ds_norm_idx[4] = {1, 4, 7, 10};
        const int64_t ds_ch[5] = {1, 4, 16, 64, D};
        int64_t size = interp;
        for (int s = 0; s < 4; ++s) {
            const int64_t next = size / 2;
            const int64_t ch = ds_ch[s + 1];
            // The stage output is allocated in the ENCLOSING scope: it becomes
            // the next stage's input, so it has to outlive the scratch below.
            // The four of them together are ~21 MiB.
            Tensor normed = nn::arena_tensor(arena, nn::DType::F32, next, next, ch);
            {
                vk::ArenaScope stage(arena);
                Tensor conv = nn::arena_tensor(arena, nn::DType::F32, next, next, ch);
                nn::ConvOpts co;
                co.stride_y = co.stride_x = 2;
                co.pad_y = co.pad_x = 1;
                co.bias = W.getf("mem_enc.ds.%d.bias", ds_conv_idx[s]);
                nn::conv2d(arena, conv, cur,
                           W.getf("mem_enc.ds.%d.weight", ds_conv_idx[s]), 3, 3, co);
                nn::layer_norm(normed.view(next * next, ch), conv.view(next * next, ch),
                               W.getf("mem_enc.ds.%d.weight", ds_norm_idx[s]),
                               W.getf("mem_enc.ds.%d.bias", ds_norm_idx[s]));
            }
            nn::unary(normed.view(next * next, ch), normed.view(next * next, ch),
                      Act::GeluErf);
            cur = normed;
            size = next;
        }
    }
    Tensor mask_feat = nn::arena_tensor(arena, nn::DType::F32, G * G, D);
    {
        nn::LinearOpts lo;
        lo.bias = W.get("mem_enc.ds.12.bias");
        nn::linear(mask_feat, cur.view(G * G, D), W.get("mem_enc.ds.12.weight"), lo);
    }

    // Fuse with the projected pixel features (an ADD, not a multiply).
    Tensor fused = nn::arena_tensor(arena, nn::DType::F32, G * G, D);
    {
        nn::LinearOpts lo;
        lo.bias = W.get("mem_enc.pix_feat_proj.bias");
        lo.residual = mask_feat;
        nn::linear(fused, pix_feats.view(G * G, D), W.get("mem_enc.pix_feat_proj.weight"),
                   lo);
    }

    // Two CXBlocks: depthwise 7x7, LayerNorm2d, pointwise MLP, LayerScale
    // residual.
    for (int i = 0; i < 2; ++i) {
        vk::ArenaScope blk(arena);
        char p[48];
        std::snprintf(p, sizeof p, "mem_enc.fuser.%d", i);
        auto wt = [&](const char* s) { return W.getf("%s.%s", p, s); };

        Tensor dw = nn::arena_tensor(arena, nn::DType::F32, G, G, D);
        nn::ConvOpts co;
        co.pad_y = co.pad_x = 3;
        co.bias = wt("dwconv.bias");
        nn::conv2d_depthwise(dw, fused.view(G, G, D), wt("dwconv.weight"), 7, 7, co);
        nn::layer_norm(dw.view(G * G, D), dw.view(G * G, D), wt("norm.weight"),
                       wt("norm.bias"));

        Tensor hid = nn::arena_tensor(arena, nn::DType::F32, G * G, wt("pwconv1.weight").rows());
        nn::LinearOpts l1;
        l1.bias = wt("pwconv1.bias");
        l1.act = Act::GeluErf;
        nn::linear(hid, dw.view(G * G, D), wt("pwconv1.weight"), l1);
        Tensor proj = nn::arena_tensor(arena, nn::DType::F32, G * G, D);
        nn::LinearOpts l2;
        l2.bias = wt("pwconv2.bias");
        nn::linear(proj, hid, wt("pwconv2.weight"), l2);

        // x = x + gamma * f(x), gamma broadcast per channel.
        nn::mul(proj, proj, wt("gamma"));
        nn::add(fused, fused, proj);
    }

    nn::LinearOpts lo;
    lo.bias = W.get("mem_enc.out_proj.bias");
    nn::linear(out.view(G * G, MD), fused, W.get("mem_enc.out_proj.weight"), lo);
}

// ---------------------------------------------------------------------------
// Memory attention
// ---------------------------------------------------------------------------

void memory_attention(const SamModel& m, vk::Arena& arena, const Tensor& curr,
                      const Tensor& curr_pos, const Tensor& prompt,
                      const Tensor& prompt_pos, int64_t n_mem, int n_ptr_tokens,
                      const Tensor& out) {
    NN_SCOPED_TIMER("memory_attention");
    const Hparams& h = m.hp();
    const WeightStore& W = m.w();
    const int64_t D = h.neck_dim;
    const int64_t MD = h.mem_out_dim;
    const int64_t G = m.grid();
    const int64_t N = G * G;
    const int64_t n_spatial = n_mem - n_ptr_tokens;
    NN_CHECK(n_spatial >= 0 && n_spatial % N == 0,
               "memory_attention: %lld spatial memory tokens is not a multiple of the "
               "%lld-token grid", (long long)n_spatial, (long long)N);
    const int64_t n_frames = n_spatial / N;

    vk::ArenaScope scope(arena);

    // pos_enc_at_input: the current frame's own sinusoidal PE enters at 0.1
    // weight, not 1.0.
    nn::add(out, curr, curr_pos, 1.0f, 0.1f);

    // K over spatial memory is rotated with the SAME grid table for every
    // stored frame (rope_k_repeat in the reference), and the object-pointer
    // tokens at the tail are not rotated at all.
    Tensor kv_pos = nn::arena_tensor(arena, nn::DType::F32, n_mem, MD);
    nn::add(kv_pos, prompt, prompt_pos);

    for (int i = 0; i < h.mem_attn_layers; ++i) {
        vk::ArenaScope blk(arena);
        char p[48];
        std::snprintf(p, sizeof p, "mem_attn.layers.%d", i);
        auto wt = [&](const char* s) { return W.getf("%s.%s", p, s); };

        // ---- self-attention with rotary position, one 256-dim head ----
        {
            Tensor xn = nn::arena_tensor(arena, nn::DType::F32, N, D);
            nn::layer_norm(xn, out, wt("norm1.weight"), wt("norm1.bias"));

            Tensor q = nn::arena_tensor(arena, nn::DType::F32, N, D);
            Tensor k = nn::arena_tensor(arena, nn::DType::F32, N, D);
            Tensor v = nn::arena_tensor(arena, nn::DType::F32, N, D);
            auto proj = [&](const Tensor& dst, const char* wn, const char* bn) {
                nn::LinearOpts lo;
                lo.bias = wt(bn);
                nn::linear(dst, xn, wt(wn), lo);
            };
            proj(q, "sa.q_proj.weight", "sa.q_proj.bias");
            proj(k, "sa.k_proj.weight", "sa.k_proj.bias");
            proj(v, "sa.v_proj.weight", "sa.v_proj.bias");
            nn::rope(q, m.memRope(), 1, (int)D, N);
            nn::rope(k, m.memRope(), 1, (int)D, N);

            Tensor ctx = nn::arena_tensor(arena, nn::DType::F32, N, D);
            nn::AttnOpts ao;
            ao.arena = &arena;
            ao.n_heads = 1;
            ao.head_dim = (int)D;
            nn::attention(ctx, q, k, v, N, N, ao);

            nn::LinearOpts lo;
            lo.bias = wt("sa.out_proj.bias");
            lo.residual = out;
            nn::linear(out, ctx, wt("sa.out_proj.weight"), lo);
        }

        // ---- cross-attention to the memory bank (kv_in_dim = 64) ----
        {
            Tensor xn = nn::arena_tensor(arena, nn::DType::F32, N, D);
            nn::layer_norm(xn, out, wt("norm2.weight"), wt("norm2.bias"));

            Tensor q = nn::arena_tensor(arena, nn::DType::F32, N, D);
            Tensor k = nn::arena_tensor(arena, nn::DType::F32, n_mem, D);
            Tensor v = nn::arena_tensor(arena, nn::DType::F32, n_mem, D);
            {
                nn::LinearOpts lo;
                lo.bias = wt("ca.q_proj.bias");
                nn::linear(q, xn, wt("ca.q_proj.weight"), lo);
                // Keys see the memory's positional encoding, values do not.
                lo.bias = wt("ca.k_proj.bias");
                nn::linear(k, kv_pos, wt("ca.k_proj.weight"), lo);
                lo.bias = wt("ca.v_proj.bias");
                nn::linear(v, prompt, wt("ca.v_proj.weight"), lo);
            }
            nn::rope(q, m.memRope(), 1, (int)D, N);
            // One rope call per stored frame: the table is the grid's, repeated,
            // so this avoids materializing n_frames copies of it.
            for (int64_t f = 0; f < n_frames; ++f)
                nn::rope(k.offsetElems(f * N * D), m.memRope(), 1, (int)D, N);

            Tensor ctx = nn::arena_tensor(arena, nn::DType::F32, N, D);
            nn::AttnOpts ao;
            ao.arena = &arena;
            ao.n_heads = 1;
            ao.head_dim = (int)D;
            nn::attention(ctx, q, k, v, N, n_mem, ao);

            nn::LinearOpts lo;
            lo.bias = wt("ca.out_proj.bias");
            lo.residual = out;
            nn::linear(out, ctx, wt("ca.out_proj.weight"), lo);
        }

        // ---- FFN ----
        {
            Tensor xn = nn::arena_tensor(arena, nn::DType::F32, N, D);
            nn::layer_norm(xn, out, wt("norm3.weight"), wt("norm3.bias"));
            Tensor hid = nn::arena_tensor(arena, nn::DType::F32, N,
                                          wt("linear1.weight").rows());
            nn::LinearOpts l1;
            l1.bias = wt("linear1.bias");
            l1.act = Act::Relu;
            nn::linear(hid, xn, wt("linear1.weight"), l1);
            nn::LinearOpts l2;
            l2.bias = wt("linear2.bias");
            l2.residual = out;
            nn::linear(out, hid, wt("linear2.weight"), l2);
        }
    }

    nn::layer_norm(out, out, W.get("mem_attn.norm.weight"), W.get("mem_attn.norm.bias"));
}

}  // namespace model
}  // namespace sam
