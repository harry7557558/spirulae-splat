// ViT backbone (Meta's Perception Encoder trunk) + SimpleFPN neck.
//
// This is ~90% of the model's FLOPs: 32 blocks over 5184 tokens at width 1024,
// with 4 global-attention blocks and 28 windowed ones, then two necks (one for
// the detector, one for the tracker) that upsample to 288x288.
//
// Differences from the ggml original, all of them consequences of the
// channel-last layout and the fused ops:
//   * no permutes anywhere -- the fused qkv projection lands in exactly the
//     layout the attention kernel indexes, and LayerNorm2d is LayerNorm;
//   * residual adds ride along on the projection that produces them;
//   * the neck's transposed convolutions run as GEMM + scatter.

#include "sam/Common.h"
#include "nn/core/Error.h"
#include "nn/core/Log.h"
#include "nn/core/Parallel.h"
#include "sam/model/Modules.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace sam {
namespace model {

using nn::Act;
using nn::Tensor;

// ---------------------------------------------------------------------------
// Preprocessing
// ---------------------------------------------------------------------------

std::vector<float> preprocess_image(const SamModel& m, const uint8_t* rgb, int width,
                                    int height) {
    const int img_size = m.imgSize();
    // SAM 3 normalizes to [-1, 1]; SAM 2 uses the ImageNet statistics its Hiera
    // trunk was pretrained with.
    const bool sam2 = m.family() == Family::Sam2;
    const float mean[3] = {sam2 ? 0.485f : 0.5f, sam2 ? 0.456f : 0.5f,
                           sam2 ? 0.406f : 0.5f};
    const float inv_std[3] = {1.0f / (sam2 ? 0.229f : 0.5f),
                              1.0f / (sam2 ? 0.224f : 0.5f),
                              1.0f / (sam2 ? 0.225f : 0.5f)};
    // Resize in uint8, in double, exactly as torchvision does:
    // F.interpolate(bilinear, align_corners=False) on the uint8 tensor followed
    // by a round to uint8. Doing the arithmetic in double makes the result
    // independent of FMA/SIMD/compiler choices, so a mask computed here matches
    // one computed by the Python reference pixel for pixel.
    //
    // This is the only step of the pipeline that stays on the host, and at
    // 1008x1008x3 it is ~3 M output pixels of double-precision bilinear -- tens
    // of milliseconds on one core, which is why it runs over rows in parallel.
    // Each row is independent and reads only source rows, so splitting changes
    // no arithmetic: the result is still bit-identical to the reference.
    NN_SCOPED_TIMER("preprocess_image");
    std::vector<uint8_t> resized;
    const uint8_t* pixels = rgb;
    if (width != img_size || height != img_size) {
        resized.resize((size_t)img_size * img_size * 3);
        const double sx = (double)width / img_size;
        const double sy = (double)height / img_size;
        parallel_for(img_size, [&](int64_t y_lo, int64_t y_hi) {
            for (int64_t y = y_lo; y < y_hi; ++y) {
                double fy = (y + 0.5) * sy - 0.5;
                if (fy < 0.0) fy = 0.0;
                const int y0 = (int)fy;
                const int y1 = (y0 < height - 1) ? y0 + 1 : y0;
                const double wy = fy - y0, wy0 = 1.0 - wy;
                for (int x = 0; x < img_size; ++x) {
                    double fx = (x + 0.5) * sx - 0.5;
                    if (fx < 0.0) fx = 0.0;
                    const int x0 = (int)fx;
                    const int x1 = (x0 < width - 1) ? x0 + 1 : x0;
                    const double wx = fx - x0, wx0 = 1.0 - wx;
                    for (int c = 0; c < 3; ++c) {
                        const double p00 = rgb[((size_t)y0 * width + x0) * 3 + c];
                        const double p01 = rgb[((size_t)y0 * width + x1) * 3 + c];
                        const double p10 = rgb[((size_t)y1 * width + x0) * 3 + c];
                        const double p11 = rgb[((size_t)y1 * width + x1) * 3 + c];
                        double v =
                            wy0 * (wx0 * p00 + wx * p01) + wy * (wx0 * p10 + wx * p11);
                        int iv = (int)(v + 0.5);
                        resized[((size_t)y * img_size + x) * 3 + c] =
                            (uint8_t)std::min(255, std::max(0, iv));
                    }
                }
            }
        }, 16);
        pixels = resized.data();
    }

    std::vector<float> out((size_t)img_size * img_size * 3);
    parallel_for((int64_t)out.size(), [&](int64_t lo, int64_t hi) {
        for (int64_t i = lo; i < hi; ++i) {
            const int c = (int)(i % 3);
            out[(size_t)i] = ((float)pixels[i] / 255.0f - mean[c]) * inv_std[c];
        }
    }, 1 << 16);
    return out;
}

// ---------------------------------------------------------------------------
// ViT
// ---------------------------------------------------------------------------

namespace {

// One transformer block, in place on `x` ([G*G, E]).
void vit_block(const SamModel& m, vk::Arena& arena, const Tensor& x, int block_idx) {
    const Hparams& h = m.hp();
    const WeightStore& W = m.w();
    const int G = m.grid();
    const int E = h.vit_embed_dim;
    const int NH = h.vit_num_heads;
    const int HD = h.vit_head_dim();
    const int WS = h.vit_window_size;
    const bool global = h.is_global_attn(block_idx);
    const int64_t N_all = (int64_t)G * G;

    vk::ArenaScope scope(arena);
    char p[64];
    std::snprintf(p, sizeof p, "vit.blocks.%d", block_idx);
    auto wt = [&](const char* suffix) {
        return W.getf("%s.%s", p, suffix);
    };

    // ---- attention ----
    {
        Tensor xn = nn::arena_tensor(arena, nn::DType::F32, N_all, E);
        nn::layer_norm(xn, x, wt("norm1.weight"), wt("norm1.bias"));

        // Windowed blocks attend inside ws x ws tiles; the partition makes each
        // tile an independent batch entry for the attention kernel.
        const int nwh = (G + WS - 1) / WS, nww = (G + WS - 1) / WS;
        const int64_t batch = global ? 1 : (int64_t)nwh * nww;
        const int64_t Ntok = global ? N_all : (int64_t)WS * WS;

        Tensor src = xn;
        if (!global) {
            Tensor xw = nn::arena_tensor(arena, nn::DType::F32, batch * Ntok, E);
            nn::window_partition(xw, xn, G, G, E, WS);
            src = xw;
        }

        // One GEMM for all three projections: q/k/v then live at strides of 3E
        // inside a single buffer, which the attention kernel reads directly.
        Tensor qkv = nn::arena_tensor(arena, nn::DType::F32, batch * Ntok, 3 * E);
        {
            nn::LinearOpts lo;
            lo.bias = wt("attn.qkv.bias");
            nn::linear(qkv, src, wt("attn.qkv.weight"), lo);
        }

        Tensor freqs = wt("attn.freqs_cis");
        Tensor q = qkv;
        Tensor k = qkv.offsetElems(E);
        Tensor v = qkv.offsetElems(2 * E);
        nn::rope(q, freqs, NH, HD, Ntok, (int)batch, 3 * E);
        nn::rope(k, freqs, NH, HD, Ntok, (int)batch, 3 * E);

        Tensor ctx = nn::arena_tensor(arena, nn::DType::F32, batch * Ntok, E);
        nn::AttnOpts ao;
        ao.arena = &arena;
        ao.n_heads = NH;
        ao.head_dim = HD;
        ao.batch = (int)batch;
        ao.q_stride = ao.k_stride = ao.v_stride = 3 * E;
        nn::attention(ctx, q, k, v, Ntok, Ntok, ao);

        nn::LinearOpts lo;
        lo.bias = wt("attn.proj.bias");
        if (global) {
            // Residual folds into the projection's epilogue; `x` aliases the
            // residual, which is safe (one thread per element).
            lo.residual = x;
            nn::linear(x, ctx, wt("attn.proj.weight"), lo);
        } else {
            Tensor proj = nn::arena_tensor(arena, nn::DType::F32, batch * Ntok, E);
            nn::linear(proj, ctx, wt("attn.proj.weight"), lo);
            Tensor merged = nn::arena_tensor(arena, nn::DType::F32, N_all, E);
            nn::window_unpartition(merged, proj, G, G, E, WS);
            nn::add(x, x, merged);
        }
    }

    // ---- MLP ----
    {
        Tensor xn = nn::arena_tensor(arena, nn::DType::F32, N_all, E);
        nn::layer_norm(xn, x, wt("norm2.weight"), wt("norm2.bias"));

        Tensor hid = nn::arena_tensor(arena, nn::DType::F32, N_all, h.vit_mlp_dim);
        nn::LinearOpts l1;
        l1.bias = wt("mlp.lin1.bias");
        l1.act = Act::GeluErf;
        nn::linear(hid, xn, wt("mlp.lin1.weight"), l1);

        nn::LinearOpts l2;
        l2.bias = wt("mlp.lin2.bias");
        l2.residual = x;
        nn::linear(x, hid, wt("mlp.lin2.weight"), l2);
    }
}

// ---------------------------------------------------------------------------
// SimpleFPN neck
// ---------------------------------------------------------------------------

void neck_forward(const SamModel& m, vk::Arena& arena, const Tensor& vit_out, bool det,
                  Tensor out[4]) {
    const Hparams& h = m.hp();
    const WeightStore& W = m.w();
    const int G = m.grid();
    const int E = h.vit_embed_dim;
    const int D = h.neck_dim;
    const std::string pre = det ? "neck.det." : "neck.trk.";
    const vk::PoolSlot slot = det ? vk::PoolSlot::NeckDet : vk::PoolSlot::NeckTrk;

    // conv1x1 (as a Linear over positions) then conv3x3, the tail every scale
    // shares. `src` is [S, S, Cin], `dst` is [S, S, D].
    auto conv_tail = [&](const Tensor& dst, const Tensor& src, int scale, int64_t S,
                         int64_t Cin) {
        vk::ArenaScope tail(arena);
        Tensor mid = nn::arena_tensor(arena, nn::DType::F32, S, S, D);
        nn::LinearOpts lo;
        lo.bias = W.getf("%s%d.conv_1x1.bias", pre.c_str(), scale);
        nn::linear(mid.view(S * S, D), src.view(S * S, Cin),
                   W.getf("%s%d.conv_1x1.weight", pre.c_str(), scale), lo);

        nn::ConvOpts co;
        co.pad_y = co.pad_x = 1;
        co.bias = W.getf("%s%d.conv_3x3.bias", pre.c_str(), scale);
        nn::conv2d(arena, dst, mid, W.getf("%s%d.conv_3x3.weight", pre.c_str(), scale), 3,
                   3, co);
    };

    // scale 0: two 2x deconvs (GELU between them) to 4x, then the tail.
    {
        vk::ArenaScope s0(arena);
        Tensor t1 = nn::arena_tensor(arena, nn::DType::F32, G * 2, G * 2, 512);
        nn::conv_transpose2x2(arena, t1, vit_out, m.neckDeconv(det, 0, 0),
                              W.getf("%s0.dconv_2x2_0.bias", pre.c_str()), Act::GeluErf);
        Tensor t2 = nn::arena_tensor(arena, nn::DType::F32, G * 4, G * 4, D);
        nn::conv_transpose2x2(arena, t2, t1, m.neckDeconv(det, 0, 1),
                              W.getf("%s0.dconv_2x2_1.bias", pre.c_str()), Act::None);
        out[0] = nn::pool_tensor(slot, 0, nn::DType::F32, G * 4, G * 4, D);
        conv_tail(out[0], t2, 0, G * 4, D);
    }
    // scale 1: one 2x deconv (no activation -- the reference has none here).
    {
        vk::ArenaScope s1(arena);
        Tensor t1 = nn::arena_tensor(arena, nn::DType::F32, G * 2, G * 2, 512);
        nn::conv_transpose2x2(arena, t1, vit_out, m.neckDeconv(det, 1, 0),
                              W.getf("%s1.dconv_2x2.bias", pre.c_str()), Act::None);
        out[1] = nn::pool_tensor(slot, 1, nn::DType::F32, G * 2, G * 2, D);
        conv_tail(out[1], t1, 1, G * 2, 512);
    }
    // scale 2: identity resolution.
    {
        out[2] = nn::pool_tensor(slot, 2, nn::DType::F32, G, G, D);
        conv_tail(out[2], vit_out, 2, G, E);
    }
    // scale 3: 2x2 max pool then the tail. Dropped by `scalp=1` in the
    // reference's detector, but the tracker's neck still produces it and the
    // cost is negligible, so keep it available.
    {
        vk::ArenaScope s3(arena);
        Tensor pooled = nn::arena_tensor(arena, nn::DType::F32, G / 2, G / 2, E);
        nn::maxpool2x2(pooled, vit_out);
        out[3] = nn::pool_tensor(slot, 3, nn::DType::F32, G / 2, G / 2, D);
        conv_tail(out[3], pooled, 3, G / 2, E);
    }
}

}  // namespace

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

void encode_image(const SamModel& m, vk::Arena& arena, const Tensor& image,
                  ImageFeatures& out) {
    if (m.family() == Family::Sam2) {
        encode_image_hiera(m, arena, image, out);
        return;
    }
    encode_image_vit(m, arena, image, out);
}

void encode_image_vit(const SamModel& m, vk::Arena& arena, const Tensor& image,
                      ImageFeatures& out) {
    NN_SCOPED_TIMER("encode_image");
    const Hparams& h = m.hp();
    const WeightStore& W = m.w();
    const int G = m.grid();
    const int E = h.vit_embed_dim;
    const int P = h.patch_size;
    const int64_t N = (int64_t)G * G;

    vk::ArenaScope scope(arena);

    // Patch embedding: a stride-14 kernel-14 convolution, i.e. a patch
    // extraction followed by one GEMM against the [1024, 3*14*14] weight (which
    // is the checkpoint's memory layout unchanged).
    Tensor x = nn::arena_tensor(arena, nn::DType::F32, N, E);
    {
        vk::ArenaScope inner(arena);
        Tensor patches = nn::arena_tensor(arena, nn::DType::F32, N, (int64_t)P * P * 3);
        nn::patchify(patches, image, P);
        nn::linear(x, patches, W.get("vit.patch_embed.proj.weight"));
    }

    // Absolute positional embedding: the checkpoint ships the 24x24 grid from
    // the 336-pixel pretraining and it is TILED (repeated), not interpolated.
    {
        Tensor pos = W.get("vit.pos_embed");
        const int th = (int)pos.dim(pos.ndim - 3);
        const int tw = (int)pos.dim(pos.ndim - 2);
        NN_CHECK(pos.dim(pos.ndim - 1) == E,
               "vit.pos_embed has %lld channels, expected %d",
               (long long)pos.dim(pos.ndim - 1), E);
        nn::add_tiled(x, x, pos, G, G, E, th, tw);
    }
    nn::layer_norm(x, x, W.get("vit.ln_pre.weight"), W.get("vit.ln_pre.bias"));

    for (int i = 0; i < h.vit_depth; ++i) {
        vit_block(m, arena, x, i);
        if (log_level() >= 3) {
            char label[48];
            std::snprintf(label, sizeof label, "vit.block[%d]", i);
            nn::tensor_debug_dump(label, x);
        }
    }

    Tensor vit_out = x.view(G, G, E);
    if (!m.visualOnly()) neck_forward(m, arena, vit_out, /*det=*/true, out.det);
    neck_forward(m, arena, vit_out, /*det=*/false, out.trk);
    out.grid = G;
}

}  // namespace model
}  // namespace sam
