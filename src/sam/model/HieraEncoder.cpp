// SAM 2's image backbone: the Hiera trunk plus an FPN neck.
//
// Hiera is a hierarchical ViT. Four stages run at strides 4, 8, 16 and 32, and
// each stage transition doubles the channel count while a 2x2 max pool on the
// *queries* halves the resolution -- so the residual has to be pooled too, and
// the block's window is the one from the stage it is leaving, not the one it is
// entering. Most blocks attend inside a window; three per model attend
// globally. The FPN then maps every stage to 256 channels, with a single
// top-down add into the stride-16 level.
//
// This is ~95% of SAM 2's per-frame cost, and it is roughly a thirtieth of
// SAM 3's ViT: 1024x1024 into 12 blocks at width 96 (tiny) or 48 at width 144
// (large), against 32 blocks at width 1024 over a 72x72 grid.
//
// Two details are worth stating because getting either wrong produces a model
// that runs and returns plausible garbage:
//
//   * window_partition pads with ZEROS and the padded tokens are real keys --
//     the reference does not mask them. So the partition has to happen before
//     the qkv projection, where a padded row becomes the projection's bias
//     rather than a zero.
//   * LayerNorm here uses eps = 1e-6, not the 1e-5 everything else in this
//     port defaults to.

#include "sam/Common.h"
#include "nn/core/Error.h"
#include "nn/core/Log.h"
#include "sam/model/Modules.h"

#include <algorithm>
#include <cstdio>
#include <vector>

namespace sam {
namespace model {

using nn::Act;
using nn::Tensor;

namespace {

constexpr float kHieraEps = 1e-6f;
// Row-chunk budget for the block MLP's hidden activation; see the note there.
constexpr int64_t kMlpChunkBytes = 32ll << 20;

// One MultiScaleBlock. `x` is [H, W, b.dim_in]; `out` is [Ho, Wo, b.dim_out],
// with the resolution halved exactly when the block pools.
void hiera_block(const SamModel& m, vk::Arena& arena, const Tensor& out, const Tensor& x,
                 int H, int W, const HieraBlock& b, int idx) {
    const WeightStore& WS = m.w();
    char p[48];
    std::snprintf(p, sizeof p, "hiera.blocks.%d", idx);
    auto wt = [&](const char* s) { return WS.getf("%s.%s", p, s); };

    const int Ci = b.dim_in, Co = b.dim_out;
    const int pool = b.q_pool ? 2 : 1;
    const int Ho = H / pool, Wo = W / pool;
    const int HD = Co / b.n_heads;
    const int64_t Nin = (int64_t)H * W, Nout = (int64_t)Ho * Wo;

    // Window geometry. ws = 0 is global attention, which is the same code path
    // with one window covering everything and no padding.
    const int ws = b.window;
    const int nwh = ws ? (H + ws - 1) / ws : 1;
    const int nww = ws ? (W + ws - 1) / ws : 1;
    const int64_t nwin = (int64_t)nwh * nww;
    const int64_t tk = ws ? (int64_t)ws * ws : Nin;             // key tokens per window
    const int64_t tq = b.q_pool ? tk / 4 : tk;                  // query tokens per window
    NN_CHECK(!b.q_pool || ws == 0 || (ws % 2) == 0,
               "hiera block %d pools a %d-wide window; the pool would straddle two "
               "windows", idx, ws);
    NN_CHECK(!b.q_pool || Ci != Co,
               "hiera block %d pools without changing dimension, so its residual would "
               "not be pooled either", idx);
    NN_CHECK(Co % b.n_heads == 0, "hiera block %d: %d channels over %d heads", idx, Co,
               b.n_heads);

    vk::ArenaScope scope(arena);

    // ---- residual, and the partitioned normalized input -------------------
    // These two are allocated before the LayerNorm scratch so that the rewind
    // below drops the full-resolution normalized map -- the block's
    // second-largest tensor -- before the projections allocate.
    Tensor shortcut = x.view(Nin, Ci);
    Tensor sc_buf;
    if (Ci != Co) sc_buf = nn::arena_tensor(arena, nn::DType::F32, Nout, Co);
    Tensor src = nn::arena_tensor(arena, nn::DType::F32, nwin * tk, Ci);

    // `proj` sees the NORMALIZED input: the reference rebinds x to norm1(x)
    // before taking the shortcut, which is easy to miss reading it.
    auto take_shortcut = [&](const Tensor& xn) {
        if (Ci == Co) return;
        vk::ArenaScope proj_scope(arena);
        Tensor projected =
            b.q_pool ? nn::arena_tensor(arena, nn::DType::F32, H, W, Co) : sc_buf;
        nn::LinearOpts lo;
        lo.bias = wt("proj.bias");
        nn::linear(projected.view(Nin, Co), xn, wt("proj.weight"), lo);
        if (b.q_pool) nn::maxpool2x2(sc_buf.view(Ho, Wo, Co), projected);
        shortcut = sc_buf;
    };

    if (ws) {
        vk::ArenaScope inner(arena);
        Tensor xn = nn::arena_tensor(arena, nn::DType::F32, Nin, Ci);
        nn::layer_norm(xn, x.view(Nin, Ci), wt("norm1.weight"), wt("norm1.bias"),
                       kHieraEps);
        take_shortcut(xn);
        nn::window_partition(src, xn.view(H, W, Ci), H, W, Ci, ws);
    } else {
        // Global attention: one "window" covering the map, so the normalized
        // input already is the partitioned input.
        nn::layer_norm(src, x.view(Nin, Ci), wt("norm1.weight"), wt("norm1.bias"),
                       kHieraEps);
        take_shortcut(src);
    }

    // ---- attention --------------------------------------------------------
    Tensor ctx = nn::arena_tensor(arena, nn::DType::F32, nwin * tq, Co);
    {
        const Tensor qkv_w = wt("attn.qkv.weight");   // [3*Co, Ci]
        const Tensor qkv_b = wt("attn.qkv.bias");
        nn::AttnOpts ao;
        ao.arena = &arena;
        ao.n_heads = b.n_heads;
        ao.head_dim = HD;
        ao.batch = (int)nwin;

        if (b.q_pool) {
            // Q is pooled inside its window, so it needs its own compact
            // buffer; K and V stay fused. Rows of the partitioned tensor run
            // (window, wy, wx) and ws is even, so viewing them as an
            // [nwin*ws, ws] map lets the ordinary 2x2 pool do the job without
            // ever crossing a window boundary.
            Tensor q = nn::arena_tensor(arena, nn::DType::F32, nwin * tq, Co);
            {
                vk::ArenaScope inner(arena);
                Tensor q_full = nn::arena_tensor(arena, nn::DType::F32, nwin * tk, Co);
                nn::LinearOpts lo;
                lo.bias = qkv_b.slice0(0, Co);
                nn::linear(q_full, src, qkv_w.slice0(0, Co), lo);
                const int64_t rows = ws ? nwin * ws : H;
                const int64_t cols = ws ? ws : W;
                nn::maxpool2x2(q.view(rows / 2, cols / 2, Co),
                               q_full.view(rows, cols, Co));
            }
            Tensor kv = nn::arena_tensor(arena, nn::DType::F32, nwin * tk, 2 * Co);
            {
                nn::LinearOpts lo;
                lo.bias = qkv_b.slice0(Co, 2 * Co);
                nn::linear(kv, src, qkv_w.slice0(Co, 2 * Co), lo);
            }
            ao.k_stride = ao.v_stride = 2 * Co;
            nn::attention(ctx, q, kv, kv.offsetElems(Co), tq, tk, ao);
        } else {
            Tensor qkv = nn::arena_tensor(arena, nn::DType::F32, nwin * tk, 3 * Co);
            nn::LinearOpts lo;
            lo.bias = qkv_b;
            nn::linear(qkv, src, qkv_w, lo);
            ao.q_stride = ao.k_stride = ao.v_stride = 3 * Co;
            nn::attention(ctx, qkv, qkv.offsetElems(Co), qkv.offsetElems(2 * Co), tk, tk,
                          ao);
        }
    }

    {
        vk::ArenaScope inner(arena);
        Tensor proj = nn::arena_tensor(arena, nn::DType::F32, nwin * tq, Co);
        nn::LinearOpts lo;
        lo.bias = wt("attn.proj.bias");
        nn::linear(proj, ctx, wt("attn.proj.weight"), lo);
        // After pooling the windows are half as wide, and the padded extent is
        // recomputed from the pooled map -- which is why the reference derives
        // pad_hw from the shortcut rather than reusing the pre-pool one.
        if (ws) nn::window_unpartition(out.view(Ho, Wo, Co), proj, Ho, Wo, Co, ws / pool);
        else    nn::copy(out.view(Nout, Co), proj);
        nn::add(out.view(Nout, Co), out.view(Nout, Co), shortcut.view(Nout, Co));
    }

    // ---- MLP --------------------------------------------------------------
    // Chunked over rows. The hidden layer is 4x the width over the whole
    // stride-4 map, which for Hiera-L is 151 MiB and was the encoder's arena
    // peak all by itself; the rows are independent, so bounding the chunk
    // costs a handful of extra dispatches and nothing else. Same trick, and
    // the same budget, as nn::conv2d's column buffer.
    {
        vk::ArenaScope inner(arena);
        const int64_t hidden = wt("mlp.fc1.weight").rows();
        const int64_t chunk =
            std::min(Nout, std::max<int64_t>(256, kMlpChunkBytes / (hidden * 4)));
        Tensor xn = nn::arena_tensor(arena, nn::DType::F32, chunk, Co);
        Tensor hid = nn::arena_tensor(arena, nn::DType::F32, chunk, hidden);
        for (int64_t r0 = 0; r0 < Nout; r0 += chunk) {
            const int64_t n = std::min(chunk, Nout - r0);
            Tensor rows = out.view(Nout, Co).slice0(r0, n);
            nn::layer_norm(xn.slice0(0, n), rows, wt("norm2.weight"), wt("norm2.bias"),
                           kHieraEps);
            nn::LinearOpts l1;
            l1.bias = wt("mlp.fc1.bias");
            l1.act = Act::GeluErf;
            nn::linear(hid.slice0(0, n), xn.slice0(0, n), wt("mlp.fc1.weight"), l1);
            nn::LinearOpts l2;
            l2.bias = wt("mlp.fc2.bias");
            l2.residual = rows;
            nn::linear(rows, hid.slice0(0, n), wt("mlp.fc2.weight"), l2);
        }
    }
}

// FpnNeck: a 1x1 conv per stage down to neck_dim, and a nearest-neighbour
// top-down add into the levels listed in fpn_top_down_levels (just the
// stride-16 one in every shipped config). `stage` is finest-first.
void fpn_neck(const SamModel& m, vk::Arena& arena, const Tensor stage[4],
              const int size[4], Tensor out[4]) {
    const Hparams& h = m.hp();
    const WeightStore& WS = m.w();
    const int D = h.neck_dim;

    auto is_top_down = [&](int level) {
        for (int i = 0; i < h.fpn_top_down_n && i < 4; ++i)
            if (h.fpn_top_down[i] == level) return true;
        return false;
    };

    // Levels the scalp discards still have to be computed: the top-down path
    // reads them. They live in the arena, the surviving ones in the pool.
    vk::ArenaScope scope(arena);
    const int keep = h.hiera_n_stages - h.scalp;
    Tensor prev;
    int prev_size = 0;
    for (int level = h.hiera_n_stages - 1; level >= 0; --level) {
        const int64_t S = size[level];
        const int64_t Cin = stage[level].cols();
        // convs[] is indexed coarsest-first, the reverse of the stage order.
        const int conv = h.hiera_n_stages - 1 - level;

        Tensor lateral = (level < keep)
                             ? nn::pool_tensor(vk::PoolSlot::NeckTrk, (uint32_t)level,
                                               nn::DType::F32, S, S, D)
                             : nn::arena_tensor(arena, nn::DType::F32, S, S, D);
        nn::LinearOpts lo;
        lo.bias = WS.getf("fpn.convs.%d.bias", conv);
        nn::linear(lateral.view(S * S, D), stage[level].view(S * S, Cin),
                   WS.getf("fpn.convs.%d.weight", conv), lo);

        if (prev.valid() && is_top_down(level)) {
            vk::ArenaScope inner(arena);
            NN_CHECK(S == 2 * prev_size,
                       "fpn: level %d is %lldx%lld, not twice level %d's %dx%d", level,
                       (long long)S, (long long)S, level + 1, prev_size, prev_size);
            Tensor up = nn::arena_tensor(arena, nn::DType::F32, S, S, D);
            nn::upsample_nearest2x(up, prev.view(prev_size, prev_size, D));
            nn::add(lateral.view(S * S, D), lateral.view(S * S, D), up.view(S * S, D));
        }
        prev = lateral;
        prev_size = (int)S;
        if (level < keep) out[level] = lateral;
    }
}

}  // namespace

void encode_image_hiera(const SamModel& m, vk::Arena& arena, const Tensor& image,
                        ImageFeatures& out) {
    NN_SCOPED_TIMER("encode_image_hiera");
    const Hparams& h = m.hp();
    const WeightStore& WS = m.w();
    const std::vector<HieraBlock>& blocks = m.hieraBlocks();
    const int S0 = m.imgSize() / 4;   // patch_embed: kernel 7, stride 4, padding 3

    vk::ArenaScope scope(arena);

    // Stage outputs live across the whole trunk, so they are allocated first
    // and never rewound until this scope exits.
    Tensor stage[4];
    int stage_size[4] = {0, 0, 0, 0};
    {
        int size = S0, dim = h.hiera_embed_dim, s = 0;
        for (size_t i = 0; i < blocks.size(); ++i) {
            if (blocks[i].q_pool) size /= 2;
            dim = blocks[i].dim_out;
            if (!blocks[i].stage_end) continue;
            NN_CHECK(s < 4, "hiera declares more than four stages");
            stage[s] = nn::arena_tensor(arena, nn::DType::F32, size, size, dim);
            stage_size[s] = size;
            ++s;
        }
        NN_CHECK(s == h.hiera_n_stages, "hiera block table produced %d stage ends, "
                   "expected %d", s, h.hiera_n_stages);
    }

    // Two slabs sized for the largest block tensor, alternating. A bump
    // allocator cannot free the older of two live buffers, and copying the
    // block output down each time would move ~1.8 GiB over a 48-block trunk.
    int64_t slab_elems = (int64_t)S0 * S0 * h.hiera_embed_dim;
    {
        int size = S0;
        for (const HieraBlock& b : blocks) {
            if (b.q_pool) size /= 2;
            slab_elems = std::max(slab_elems, (int64_t)size * size * b.dim_out);
        }
    }
    Tensor slab[2] = {nn::arena_tensor(arena, nn::DType::F32, slab_elems),
                      nn::arena_tensor(arena, nn::DType::F32, slab_elems)};

    // ---- patch embedding + absolute position ------------------------------
    int H = S0, W = S0, C = h.hiera_embed_dim;
    Tensor x = slab[0].view(H, W, C);
    {
        vk::ArenaScope inner(arena);
        nn::ConvOpts co;
        co.stride_y = co.stride_x = 4;
        co.pad_y = co.pad_x = 3;
        co.bias = WS.get("hiera.patch_embed.bias");
        nn::conv2d(arena, x, image, WS.get("hiera.patch_embed.weight"), 7, 7, co);
    }
    nn::add(x.view((int64_t)H * W, C), x.view((int64_t)H * W, C),
            m.hieraPosEmbed().view((int64_t)H * W, C));

    // ---- trunk ------------------------------------------------------------
    int cur = 0, s = 0;
    for (size_t i = 0; i < blocks.size(); ++i) {
        const HieraBlock& b = blocks[i];
        const int Ho = b.q_pool ? H / 2 : H, Wo = b.q_pool ? W / 2 : W;
        Tensor y = slab[cur ^ 1].view(Ho, Wo, b.dim_out);
        hiera_block(m, arena, y, x, H, W, b, (int)i);
        x = y;
        H = Ho;
        W = Wo;
        C = b.dim_out;
        cur ^= 1;
        if (b.stage_end) {
            nn::strided_copy(stage[s].view((int64_t)H * W, C), x.view((int64_t)H * W, C),
                             (int64_t)H * W, C, C, C);
            ++s;
        }
        if (log_level() >= 3) {
            char label[48];
            std::snprintf(label, sizeof label, "hiera.block[%zu]", i);
            nn::tensor_debug_dump(label, x.view((int64_t)H * W, C));
        }
    }

    fpn_neck(m, arena, stage, stage_size, out.trk);
    out.grid = m.grid();
    NN_CHECK(out.trk[2].valid() && stage_size[2] == m.grid(),
               "hiera produced a %dx%d stride-16 level, expected %dx%d", stage_size[2],
               stage_size[2], m.grid(), m.grid());
}

}  // namespace model
}  // namespace sam
