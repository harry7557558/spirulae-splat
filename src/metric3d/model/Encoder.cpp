// The DINOv2-with-registers encoder, read against ViT_DINO_reg.py.
//
// Standard pre-norm ViT with two details that are not:
//
//   * The positional embedding is stored for a 37x37 patch grid and bicubically
//     resampled to whatever grid the input produces. DINOv2 asks torch for a
//     scale of (grid + 0.1) / 37 -- its own guard against a float comparison
//     landing on the wrong side of an integer -- but the export emits a
//     SIZES-based Resize, and an ONNX Resize given sizes maps coordinates by
//     out/in. So the 0.1 decides the output size, which is `grid` either way,
//     and nothing else. Carrying it into the mapping shifts every sample by a
//     fraction of a source patch and costs ~4% on the embedding.
//   * The register tokens carry NO positional embedding. The order is
//     [cls, registers, patches], but the embedding is added to [cls, patches]
//     before the registers are spliced in.
//
// LayerScale is applied as its own multiply rather than folded into the
// projection weight it follows. Folding is exact in fp32 and would save a pass,
// but the weights go to the device as f16 and a trained gamma reaches 1e-4:
// the product of a small weight and a small gamma is subnormal there, so the
// fold would quietly zero layers.

#include "metric3d/model/Model.h"

#include "metric3d/Common.h"
#include "nn/Ops.h"
#include "nn/vk/Stream.h"

#include <cmath>
#include <vector>

namespace metric3d {
namespace {

using nn::Act;
using nn::DType;
using nn::LinearOpts;
using nn::Tensor;

// torch's bicubic with align_corners=False and torch's own -0.75 kernel, which
// is what F.interpolate(mode='bicubic') and the exported Resize both use.
float cubic_weight(float t, float a) {
    t = std::fabs(t);
    if (t <= 1.0f) return ((a + 2.0f) * t - (a + 3.0f)) * t * t + 1.0f;
    if (t < 2.0f) return ((a * t - 5.0f * a) * t + 8.0f * a) * t - 4.0f * a;
    return 0.0f;
}

// [gi, gi, C] -> [ho, wo, C], torch's bicubic with align_corners=False:
// src = (dst + 0.5) * in/out - 0.5, indices clamped to the edge and their
// weights kept (ONNX's exclude_outside=0).
std::vector<float> bicubic_resize(const std::vector<float>& src, int gi, int64_t C,
                                  int64_t ho, int64_t wo) {
    std::vector<float> dst((size_t)ho * wo * C);
    const float a = -0.75f;
    for (int64_t oy = 0; oy < ho; ++oy) {
        const double sy = ((double)oy + 0.5) * (double)gi / (double)ho - 0.5;
        const int64_t y0 = (int64_t)std::floor(sy);
        float wy[4];
        for (int k = 0; k < 4; ++k) wy[k] = cubic_weight((float)(sy - (y0 - 1 + k)), a);
        for (int64_t ox = 0; ox < wo; ++ox) {
            const double sx = ((double)ox + 0.5) * (double)gi / (double)wo - 0.5;
            const int64_t x0 = (int64_t)std::floor(sx);
            float wx[4];
            for (int k = 0; k < 4; ++k) wx[k] = cubic_weight((float)(sx - (x0 - 1 + k)), a);
            float* out = &dst[(size_t)(oy * wo + ox) * C];
            for (int ky = 0; ky < 4; ++ky) {
                const int64_t yy = std::min<int64_t>(std::max<int64_t>(y0 - 1 + ky, 0), gi - 1);
                for (int kx = 0; kx < 4; ++kx) {
                    const int64_t xx =
                        std::min<int64_t>(std::max<int64_t>(x0 - 1 + kx, 0), gi - 1);
                    const float w = wy[ky] * wx[kx];
                    if (w == 0.0f) continue;
                    const float* s = &src[(size_t)(yy * gi + xx) * C];
                    for (int64_t c = 0; c < C; ++c) out[c] += w * s[c];
                }
            }
        }
    }
    return dst;
}

}  // namespace

Model::~Model() {
    if (pos_blob) vk::device_free(pos_blob);
}

void Model::ensurePosEmbed(int64_t gh, int64_t gw) {
    if (pos_gh == gh && pos_gw == gw && pos_blob) return;
    const Hparams& h = hp();
    const int64_t D = h.embed_dim;
    const int gi = weights.posGrid();

    std::vector<float> patch = bicubic_resize(weights.patchPos(), gi, D, gh, gw);

    if (pos_blob) vk::device_free(pos_blob);
    const uint64_t bytes = (uint64_t)(patch.size() + (size_t)D) * 4;
    pos_blob = vk::device_alloc(bytes, "metric3d-pos-embed");
    pos_patch = Tensor(pos_blob, DType::F32, gh * gw, D);
    pos_cls = Tensor(pos_blob + (uint64_t)patch.size() * 4, DType::F32, D);
    nn::tensor_from_host(pos_patch, patch.data(), (int64_t)patch.size());
    nn::tensor_from_host(pos_cls, weights.clsPos().data(), D);
    vk::Stream::get().sync();
    pos_gh = gh;
    pos_gw = gw;
}

Features Model::encode(const Tensor& image) {
    const Hparams& h = hp();
    const int64_t D = h.embed_dim;
    const int64_t P = h.patch;
    const int64_t gh = image.shape[0] / P, gw = image.shape[1] / P;
    const int64_t np = gh * gw;
    const int64_t pre = h.prefix_tokens();
    const int64_t N = pre + np;
    ensurePosEmbed(gh, gw);

    Features out;
    out.gh = gh;
    out.gw = gw;

    // The token sequence is allocated OUTSIDE the block scopes: it is the one
    // thing that stays live from here to the decoder.
    Tensor x = nn::arena_tensor(arena, DType::F32, N, D);
    {
        vk::ArenaScope scope(arena);
        // patch_embed is a stride-14 convolution, so its [gh, gw, D] output is
        // already the patch token sequence in memory.
        Tensor patches = nn::arena_tensor(arena, DType::F32, gh, gw, D);
        nn::ConvOpts co;
        co.stride_y = co.stride_x = (int)P;
        co.bias = weights.get("encoder.patch_embed.proj.bias");
        nn::conv2d(arena, patches, image, weights.get("encoder.patch_embed.proj.weight"),
                   (int)P, (int)P, co);

        Tensor patch_rows(patches.ptr, DType::F32, np, D);
        nn::add(patch_rows, patch_rows, pos_patch);

        // [cls + pos] [registers] [patches + pos]
        Tensor cls_row(x.ptr, DType::F32, 1, D);
        nn::add(cls_row, weights.get("encoder.cls_token"), pos_cls);
        nn::strided_copy(Tensor(x.ptr + (uint64_t)D * 4, DType::F32, pre - 1, D),
                         weights.get("encoder.register_tokens"), pre - 1, D, D, D);
        nn::strided_copy(Tensor(x.ptr + (uint64_t)pre * D * 4, DType::F32, np, D),
                         patch_rows, np, D, D, D);
    }
    if (dump_enabled()) dump_tensor("tokens_in", x, {N, D});

    // A tapped block's output is read again after later blocks have overwritten
    // x, so the four of them are allocated here rather than inside a block's
    // scope -- which would rewind them at the end of that iteration.
    if (!h.final_norm)
        for (int i = 0; i < 4; ++i) out.tokens[i] = nn::arena_tensor(arena, DType::F32, N, D);

    int tap = 0;
    for (int b = 0; b < h.depth; ++b) {
        vk::ArenaScope scope(arena);
        Tensor t = nn::arena_tensor(arena, DType::F32, N, D);

        // x = x + ls1 * proj(attn(norm1(x)))
        nn::layer_norm(t, x, weights.getf("encoder.blocks.%d.norm1.weight", b),
                       weights.getf("encoder.blocks.%d.norm1.bias", b), 1e-6f);
        {
            vk::ArenaScope inner(arena);
            Tensor qkv = nn::arena_tensor(arena, DType::F32, N, 3 * D);
            LinearOpts lo;
            lo.bias = weights.getf("encoder.blocks.%d.attn.qkv.bias", b);
            nn::linear(qkv, t, weights.getf("encoder.blocks.%d.attn.qkv.weight", b), lo);

            Tensor attn = nn::arena_tensor(arena, DType::F32, N, D);
            nn::AttnOpts ao;
            ao.n_heads = h.num_heads;
            ao.head_dim = h.head_dim();
            ao.arena = &arena;
            ao.q_stride = ao.k_stride = ao.v_stride = 3 * D;
            nn::attention(attn, Tensor(qkv.ptr, DType::F32, N, D),
                          Tensor(qkv.ptr + (uint64_t)D * 4, DType::F32, N, D),
                          Tensor(qkv.ptr + (uint64_t)2 * D * 4, DType::F32, N, D), N, N, ao);

            LinearOpts po;
            po.bias = weights.getf("encoder.blocks.%d.attn.proj.bias", b);
            nn::linear(t, attn, weights.getf("encoder.blocks.%d.attn.proj.weight", b), po);
        }
        nn::mul(t, t, weights.getf("encoder.blocks.%d.ls1.gamma", b));
        nn::add(x, x, t);

        // x = x + ls2 * mlp(norm2(x))
        nn::layer_norm(t, x, weights.getf("encoder.blocks.%d.norm2.weight", b),
                       weights.getf("encoder.blocks.%d.norm2.bias", b), 1e-6f);
        {
            vk::ArenaScope inner(arena);
            if (h.swiglu) {
                // SwiGLUFFN: one projection to 2*hidden, split in halves,
                // silu(first) * second. The halves are contiguous column ranges
                // of the same row, so the gate is a strided elementwise pass and
                // nothing is copied.
                const int64_t Hh = h.mlp_hidden;
                Tensor w12 = nn::arena_tensor(arena, DType::F32, N, 2 * Hh);
                LinearOpts lo;
                lo.bias = weights.getf("encoder.blocks.%d.mlp.w12.bias", b);
                nn::linear(w12, t, weights.getf("encoder.blocks.%d.mlp.w12.weight", b), lo);

                Tensor gate = nn::arena_tensor(arena, DType::F32, N, Hh);
                nn::strided_copy(gate, Tensor(w12.ptr, DType::F32, N, Hh), N, Hh, 2 * Hh, Hh);
                nn::unary(gate, gate, Act::Silu);
                Tensor up = nn::arena_tensor(arena, DType::F32, N, Hh);
                nn::strided_copy(up, Tensor(w12.ptr + (uint64_t)Hh * 4, DType::F32, N, Hh),
                                 N, Hh, 2 * Hh, Hh);
                nn::mul(gate, gate, up);

                LinearOpts o3;
                o3.bias = weights.getf("encoder.blocks.%d.mlp.w3.bias", b);
                nn::linear(t, gate, weights.getf("encoder.blocks.%d.mlp.w3.weight", b), o3);
            } else {
                Tensor hid = nn::arena_tensor(arena, DType::F32, N, h.mlp_hidden);
                LinearOpts lo;
                lo.bias = weights.getf("encoder.blocks.%d.mlp.fc1.bias", b);
                lo.act = Act::GeluErf;
                nn::linear(hid, t, weights.getf("encoder.blocks.%d.mlp.fc1.weight", b), lo);

                LinearOpts o2;
                o2.bias = weights.getf("encoder.blocks.%d.mlp.fc2.bias", b);
                nn::linear(t, hid, weights.getf("encoder.blocks.%d.mlp.fc2.weight", b), o2);
            }
        }
        nn::mul(t, t, weights.getf("encoder.blocks.%d.ls2.gamma", b));
        nn::add(x, x, t);

        if (b == 0 && dump_enabled()) dump_tensor("tokens_block0", x, {N, D});
        if (!h.final_norm && tap < 4 && b == h.taps[tap]) {
            nn::copy(out.tokens[tap], x);
            if (dump_enabled())
                dump_tensor(("tokens_tap" + std::to_string(tap)).c_str(), out.tokens[tap],
                            {N, D});
            ++tap;
        }
    }

    if (h.final_norm) {
        Tensor normed = nn::arena_tensor(arena, DType::F32, N, D);
        nn::layer_norm(normed, x, weights.get("encoder.norm.weight"),
                       weights.get("encoder.norm.bias"), 1e-6f);
        for (int i = 0; i < 4; ++i) out.tokens[i] = normed;
        if (dump_enabled()) dump_tensor("tokens", normed, {N, D});
    } else {
        NN_CHECK(tap == 4, "encoder tapped %d of 4 intermediate features", tap);
    }
    return out;
}

}  // namespace metric3d
