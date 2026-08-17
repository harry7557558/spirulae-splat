// The DINOv2 encoder, read against DINOv2's own vision_transformer.py: a
// pre-norm ViT with LayerScale and no register tokens.
//
// LayerScale is applied as its own multiply rather than folded into the
// projection weight it follows. Folding is exact in fp32, but the weights go to
// the device as f16 and a trained gamma reaches 1e-4: the product of a small
// weight and a small gamma is subnormal there, so the fold would quietly zero
// layers.
//
// src/moge/README.md, "Four conventions", is the rest of it.

#include "moge/model/Model.h"

#include "moge/Common.h"
#include "nn/Ops.h"
#include "nn/vk/Stream.h"

#include <cmath>
#include <vector>

namespace moge {
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
    if (uv_blob) vk::device_free(uv_blob);
}

void Model::ensurePosEmbed(int64_t gh, int64_t gw) {
    if (pos_gh == gh && pos_gw == gw && pos_blob) return;
    const int64_t D = hp().embed_dim;

    std::vector<float> patch =
        bicubic_resize(weights.patchPos(), hp().pos_grid, D, gh, gw);

    if (pos_blob) vk::device_free(pos_blob);
    const uint64_t bytes = (uint64_t)(patch.size() + (size_t)D) * 4;
    pos_blob = vk::device_alloc(bytes, "moge-pos-embed");
    pos_patch = Tensor(pos_blob, DType::F32, gh * gw, D);
    pos_cls = Tensor(pos_blob + (uint64_t)patch.size() * 4, DType::F32, D);
    nn::tensor_from_host(pos_patch, patch.data(), (int64_t)patch.size());
    nn::tensor_from_host(pos_cls, weights.clsPos().data(), D);
    vk::Stream::get().sync();
    pos_gh = gh;
    pos_gw = gw;
}

// normalized_view_plane_uv: the image plane scaled so its DIAGONAL is 2, which
// is the frame MoGe's point map is affine in. `aspect` is the source image's,
// not the level's -- every level covers the same field of view.
void Model::ensureUv(int64_t gh, int64_t gw, double aspect) {
    if (uv_gh == gh && uv_gw == gw && uv_aspect == aspect && uv_blob) return;

    const double span_x = aspect / std::sqrt(1.0 + aspect * aspect);
    const double span_y = 1.0 / std::sqrt(1.0 + aspect * aspect);

    uint64_t total = 0;
    for (int l = 0; l < Hparams::kLevels; ++l) total += (uint64_t)(gh << l) * (gw << l) * 2;

    std::vector<float> host((size_t)total);
    uint64_t at = 0;
    for (int l = 0; l < Hparams::kLevels; ++l) {
        const int64_t h = gh << l, w = gw << l;
        float* dst = host.data() + at;
        for (int64_t y = 0; y < h; ++y) {
            const float v = (float)(span_y * (2.0 * ((double)y + 0.5) / (double)h - 1.0));
            for (int64_t x = 0; x < w; ++x) {
                dst[(size_t)(y * w + x) * 2 + 0] =
                    (float)(span_x * (2.0 * ((double)x + 0.5) / (double)w - 1.0));
                dst[(size_t)(y * w + x) * 2 + 1] = v;
            }
        }
        at += (uint64_t)h * w * 2;
    }

    if (uv_blob) vk::device_free(uv_blob);
    uv_blob = vk::device_alloc(total * 4, "moge-uv");
    at = 0;
    for (int l = 0; l < Hparams::kLevels; ++l) {
        const int64_t h = gh << l, w = gw << l;
        uv[l] = Tensor(uv_blob + at * 4, DType::F32, h, w, 2);
        at += (uint64_t)h * w * 2;
    }
    nn::tensor_from_host(Tensor(uv_blob, DType::F32, (int64_t)total), host.data(),
                         (int64_t)total);
    vk::Stream::get().sync();
    uv_gh = gh;
    uv_gw = gw;
    uv_aspect = aspect;
}

Features Model::encode(const Tensor& image) {
    const Hparams& h = hp();
    const int64_t D = h.embed_dim;
    const int64_t P = h.patch;
    const int64_t gh = image.shape[0] / P, gw = image.shape[1] / P;
    const int64_t np = gh * gw;
    const int64_t N = 1 + np;
    ensurePosEmbed(gh, gw);

    Features out;
    out.gh = gh;
    out.gw = gw;
    // Both outlive the block loop: the summed feature map is what the neck
    // reads and the class token is what the scale head reads.
    out.map = nn::arena_tensor(arena, DType::F32, gh, gw, D);
    out.cls = nn::arena_tensor(arena, DType::F32, 1, D);
    nn::fill(out.map, 0.0f);

    // Nothing else the encoder allocates outlives it, and the token sequence is
    // the largest thing in the pass after the head's feature maps.
    vk::ArenaScope tokens(arena);
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

        Tensor cls_row(x.ptr, DType::F32, 1, D);
        nn::add(cls_row, weights.get("encoder.cls_token"), pos_cls);
        nn::strided_copy(Tensor(x.ptr + (uint64_t)D * 4, DType::F32, np, D), patch_rows, np,
                         D, D, D);
    }
    if (dump_enabled()) dump_tensor("tokens_in", x, {N, D});

    Tensor normed = nn::arena_tensor(arena, DType::F32, N, D);
    int tap = 0;
    for (int b = 0; b < h.depth; ++b) {
        vk::ArenaScope scope(arena);
        Tensor t = nn::arena_tensor(arena, DType::F32, N, D);

        // x = x + ls1 * proj(attn(norm1(x)))
        nn::layer_norm(t, x, weights.getf("encoder.blocks.%d.norm1.weight", b),
                       weights.getf("encoder.blocks.%d.norm1.bias", b), Hparams::kNormEps);
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
                       weights.getf("encoder.blocks.%d.norm2.bias", b), Hparams::kNormEps);
        {
            vk::ArenaScope inner(arena);
            Tensor hid = nn::arena_tensor(arena, DType::F32, N, h.mlp_hidden);
            LinearOpts lo;
            lo.bias = weights.getf("encoder.blocks.%d.mlp.fc1.bias", b);
            lo.act = Act::GeluErf;
            nn::linear(hid, t, weights.getf("encoder.blocks.%d.mlp.fc1.weight", b), lo);

            LinearOpts o2;
            o2.bias = weights.getf("encoder.blocks.%d.mlp.fc2.bias", b);
            nn::linear(t, hid, weights.getf("encoder.blocks.%d.mlp.fc2.weight", b), o2);
        }
        nn::mul(t, t, weights.getf("encoder.blocks.%d.ls2.gamma", b));
        nn::add(x, x, t);

        if (tap >= h.n_taps || b != h.taps[tap]) continue;

        // Each tap goes through the SAME final LayerNorm, is projected to the
        // shared width and summed into the feature map. Only the last tap's
        // class token is kept.
        nn::layer_norm(normed, x, weights.get("encoder.norm.weight"),
                       weights.get("encoder.norm.bias"), Hparams::kNormEps);
        {
            LinearOpts lo;
            lo.bias = weights.getf("encoder.output_projections.%d.bias", tap);
            // The residual is the output, which nn::linear accumulates in place.
            lo.residual = out.map.view(np, D);
            nn::linear(out.map.view(np, D),
                       Tensor(normed.ptr + (uint64_t)D * 4, DType::F32, np, D),
                       weights.getf("encoder.output_projections.%d.weight", tap), lo);
        }
        if (tap == h.n_taps - 1) nn::copy(out.cls, Tensor(normed.ptr, DType::F32, 1, D));
        if (dump_enabled())
            dump_tensor(("tap" + std::to_string(tap)).c_str(), normed, {N, D});
        ++tap;
    }
    NN_CHECK(tap == h.n_taps, "the encoder tapped %d of %d blocks", tap, h.n_taps);
    if (dump_enabled()) dump_tensor("encoder_feat", out.map, {gh, gw, D});
    return out;
}

}  // namespace moge
