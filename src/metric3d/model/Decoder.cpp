// RAFTDepthNormalDPT5, read against RAFTDepthNormalDPTDecoder5.py.
//
// Three resolutions run through the whole file and the module names do not
// match any of them: what the source calls gru08 / gru16 / gru32 lives at 1/4,
// 1/7 and 1/14 of the input, because the ViT's stride is 14 and the DPT
// pyramid is built from it by factors of 7/2 and 2. The names are RAFT's,
// inherited from a stereo network with a stride-8 backbone.
//
// The refinement state is six channels -- depth, depth confidence, normal xyz,
// normal concentration -- carried at 1/4 and convex-upsampled to full
// resolution at the end. The source frames it as a RAFT flow field, with
// `coords1 - coords0` and an `initialize_flow`; coords_grid() returns zeros
// (read it), so coords0 is nothing, `flow` IS the state, and the loop is just
// `state += delta`. Nothing here reconstructs the grid.

#include "metric3d/model/Model.h"

#include "metric3d/Common.h"
#include "nn/Ops.h"
#include "nn/vk/Stream.h"

#include <cmath>
#include <string>
#include <vector>

namespace metric3d {
namespace {

using nn::Act;
using nn::DType;
using nn::LinearOpts;
using nn::Tensor;

// LayerNorm2d over a channel-last map is nn::layer_norm with no permute; the
// decoder's are torch's default epsilon, unlike the encoder's 1e-6.
constexpr float kNormEps = 1e-5f;

int64_t scaled(int64_t n, double s) { return (int64_t)std::floor((double)n * s); }

}  // namespace

// The forward pass carries a lot of state that every helper needs -- the
// weights, the arena, the working resolutions -- so the helpers are members of
// a local struct rather than free functions with eight arguments each.
namespace {

struct Decoder {
    Model&    m;
    vk::Arena& arena;
    const Hparams& h;

    Decoder(Model& model)
        : m(model), arena(model.arena), h(model.weights.hparams()) {}

    Tensor w(const std::string& n) const { return m.weights.get(n); }

    Tensor alloc(int64_t H, int64_t W, int64_t C) {
        return nn::arena_tensor(arena, DType::F32, H, W, C);
    }

    // ---- layer helpers ---------------------------------------------------

    void conv3x3(const Tensor& out, const Tensor& in, const std::string& p, Act act) {
        nn::ConvOpts co;
        co.pad_y = co.pad_x = 1;
        co.act = act;
        co.bias = w(p + ".bias");
        nn::conv2d(arena, out, in, w(p + ".weight"), 3, 3, co);
    }

    // A 1x1 convolution on a channel-last map is a Linear over its rows, and
    // the [Cout, Cin, 1, 1] weight is the same memory as [Cout, Cin] -- so this
    // skips the im2col a general conv2d would do to discover the same thing.
    void conv1x1(const Tensor& out, const Tensor& in, const std::string& p, Act act) {
        const Tensor wt = w(p + ".weight");
        LinearOpts lo;
        lo.bias = w(p + ".bias");
        lo.act = act;
        nn::linear(Tensor(out.ptr, DType::F32, out.rows(), out.shape[2]),
                   Tensor(in.ptr, DType::F32, in.rows(), in.shape[2]),
                   Tensor(wt.ptr, wt.dtype, wt.shape[0], wt.shape[1]), lo);
    }

    void layerNorm(const Tensor& out, const Tensor& in, const std::string& p) {
        nn::layer_norm(out, in, w(p + ".weight"), w(p + ".bias"), kNormEps);
    }

    // DPT's ConvBlock -- and it RECTIFIES `x` IN PLACE. That is not a liberty
    // taken here; it is what the reference does and what the rest of the
    // decoder depends on.
    //
    // `self.act` is nn.ReLU(inplace=True), so `out = self.act(x)` overwrites
    // the caller's tensor. Two things follow that the module source suggests
    // otherwise, and the exported graph settles both:
    //
    //   * the residual adds relu(x), not x -- way_trunk/Add takes
    //     act/Relu_output_0;
    //   * the rectification ESCAPES the block. DecoderFeature runs before
    //     ContextFeatureEncoder over the same feature list, so by the time the
    //     context branches read x1 and x2 the FuseBlocks have rectified them:
    //     outputs08's first convolution takes
    //     upconv_1/way_branch/act/Relu_output_0, not read_1's output. read_0
    //     never passes through a ConvBlock (upconv_0 is commented out
    //     upstream) and reaches outputs04 untouched, which is why missing this
    //     breaks two of the three context branches and not the third.
    void convBlock(const Tensor& out, const Tensor& x, const std::string& p) {
        const int64_t H = x.shape[0], W = x.shape[1], C = x.shape[2];
        vk::ArenaScope scope(arena);
        nn::unary(x, x, Act::Relu);
        Tensor b = alloc(H, W, C);
        conv3x3(b, x, p + ".conv1", Act::Relu);
        conv3x3(out, b, p + ".conv2", Act::None);
        nn::add(out, out, x);
    }

    // DPT's FuseBlock. `branch` is empty for the trunk-only first stage.
    // `scale` 0 means no upsample.
    Tensor fuseBlock(const Tensor& trunk, const Tensor& branch, const std::string& p,
                     double scale, int64_t out_ch) {
        const int64_t H = trunk.shape[0], W = trunk.shape[1], C = trunk.shape[2];
        // Without a branch the trunk goes straight into the ConvBlock, which
        // rectifies it in place exactly as the reference does. Nothing reads
        // read_3 again, so that is unobservable -- but copying to avoid it
        // would be a copy taken to hide a faithful mutation.
        Tensor fused = trunk;
        if (branch.valid()) {
            fused = alloc(H, W, C);
            vk::ArenaScope scope(arena);
            Tensor b = alloc(H, W, C);
            convBlock(b, branch, p + ".way_branch");
            nn::add(fused, trunk, b);
        }

        Tensor body = alloc(H, W, C);
        convBlock(body, fused, p + ".way_trunk");

        if (scale > 0.0) {
            const int64_t Ho = scaled(H, scale), Wo = scaled(W, scale);
            Tensor up = alloc(Ho, Wo, C);
            nn::resize_bilinear(up, body, /*align_corners=*/true);
            Tensor out = alloc(Ho, Wo, out_ch);
            conv1x1(out, up, p + ".out_conv", Act::None);
            return out;
        }
        Tensor out = alloc(H, W, out_ch);
        conv1x1(out, body, p + ".out_conv", Act::None);
        return out;
    }

    // RAFT's ResidualBlock with norm_fn='layer'. Always strided 1 here, and
    // always widening, so the downsample branch is always present.
    Tensor residualBlock(const Tensor& in, const std::string& p, int64_t out_ch) {
        const int64_t H = in.shape[0], W = in.shape[1];
        Tensor out = alloc(H, W, out_ch);
        vk::ArenaScope scope(arena);
        Tensor y = alloc(H, W, out_ch);
        Tensor t = alloc(H, W, out_ch);
        conv3x3(y, in, p + ".conv1", Act::None);
        layerNorm(t, y, p + ".norm1");
        nn::unary(t, t, Act::Relu);
        conv3x3(y, t, p + ".conv2", Act::None);
        layerNorm(t, y, p + ".norm2");
        nn::unary(t, t, Act::Relu);

        conv1x1(y, in, p + ".downsample.0", Act::None);
        layerNorm(out, y, p + ".norm3");
        nn::add(out, out, t, 1.0f, 1.0f, Act::Relu);
        return out;
    }

    // Readout + resample: one ViT token tensor to one feature map.
    Tensor token2feature(const Tensor& tokens, int idx, int64_t gh, int64_t gw,
                         double scale, int64_t out_ch) {
        const int64_t D = h.embed_dim, pre = h.prefix_tokens(), np = gh * gw;
        const std::string p = "decoder.token2feature.read_" + std::to_string(idx);

        Tensor feat = alloc(gh, gw, D);
        {
            vk::ArenaScope scope(arena);
            // The cls and register tokens are one flat row: project_learn is a
            // Linear over (1 + registers) * D, which is exactly their memory.
            Tensor learn = nn::arena_tensor(arena, DType::F32, 1, D);
            const Tensor wl = w(p + ".readoper.project_learn.weight");
            nn::linear(learn, Tensor(tokens.ptr, DType::F32, 1, pre * D), wl);

            Tensor patch = Tensor(tokens.ptr + (uint64_t)pre * D * 4, DType::F32, np, D);
            LinearOpts lo;
            lo.bias = w(p + ".readoper.project_patch.bias");
            Tensor flat(feat.ptr, DType::F32, np, D);
            nn::linear(flat, patch, w(p + ".readoper.project_patch.weight"), lo);
            // x_learn broadcasts over every patch, and the GELU is the Readout's.
            nn::add(flat, flat, learn, 1.0f, 1.0f, Act::GeluErf);
        }

        if (scale == 1.0) return feat;
        if (scale == 2.0) {
            Tensor out = alloc(gh * 2, gw * 2, out_ch);
            nn::conv_transpose2x2(arena, out, feat, w(p + ".sample.weight"),
                                  w(p + ".sample.bias"));
            return out;
        }
        // 7/2: nearest upsample then a 1x1 projection. torch was asked for the
        // scale, so the source pixel is floor(dst / 3.5) even though the output
        // is floor(in * 3.5) -- nn::resize_nearest takes the scale for that
        // reason.
        const int64_t Ho = scaled(gh, scale), Wo = scaled(gw, scale);
        Tensor out = alloc(Ho, Wo, out_ch);
        vk::ArenaScope scope(arena);
        Tensor up = alloc(Ho, Wo, D);
        nn::resize_nearest(up, feat, (float)scale, (float)scale);
        conv1x1(out, up, p + ".sample.0", Act::None);
        return out;
    }

    // ---- the ConvGRU -----------------------------------------------------

    // h <- (1-z)*h + z*tanh(convq([r*h, x]) + cq), in place.
    //
    // `hx` is built once with h in its leading columns and the concatenated
    // inputs after; the q gate reuses the same buffer with r*h written over
    // those leading columns, so the inputs are copied once rather than twice.
    void convGru(const Tensor& state, const std::string& p, const Tensor& cz,
                 const Tensor& cr, const Tensor& cq, const std::vector<Tensor>& xs) {
        const int64_t H = state.shape[0], W = state.shape[1], C = state.shape[2];
        const int64_t rows = H * W;
        int64_t Cx = 0;
        for (const Tensor& x : xs) Cx += x.shape[2];

        vk::ArenaScope scope(arena);
        Tensor hx = alloc(H, W, C + Cx);
        const int64_t stride = C + Cx;
        nn::strided_copy(hx, state, rows, C, C, stride);
        int64_t off = C;
        for (const Tensor& x : xs) {
            const int64_t xc = x.shape[2];
            nn::strided_copy(Tensor(hx.ptr + (uint64_t)off * 4, DType::F32, rows, xc), x,
                             rows, xc, xc, stride);
            off += xc;
        }

        Tensor z = alloc(H, W, C);
        conv3x3(z, hx, p + ".convz", Act::None);
        nn::add(z, z, cz, 1.0f, 1.0f, Act::Sigmoid);

        Tensor r = alloc(H, W, C);
        conv3x3(r, hx, p + ".convr", Act::None);
        nn::add(r, r, cr, 1.0f, 1.0f, Act::Sigmoid);

        nn::mul(r, r, state);
        nn::strided_copy(hx, r, rows, C, C, stride);

        Tensor q = alloc(H, W, C);
        conv3x3(q, hx, p + ".convq", Act::None);
        nn::add(q, q, cq, 1.0f, 1.0f, Act::Tanh);

        blend(state, z, q);
    }

    void blend(const Tensor& state, const Tensor& z, const Tensor& q) {
        struct { uint64_t h, z, q; uint32_t n, groups_per_row; } bp{};
        bp.h = state.ptr;
        bp.z = z.ptr;
        bp.q = q.ptr;
        bp.n = (uint32_t)state.numel();
        vk::SpecList spec{0u, 0u};
        vk::Stream::get().dispatchFlat("metric3d.gru_blend", spec, state.numel(), 256, &bp,
                                       sizeof(bp), &bp.groups_per_row);
    }

    Tensor pool2x(const Tensor& x) {
        const int64_t H = x.shape[0], W = x.shape[1], C = x.shape[2];
        Tensor out = alloc((H + 2 - 3) / 2 + 1, (W + 2 - 3) / 2 + 1, C);
        nn::avgpool(out, x, 3, 2, 1);
        return out;
    }

    Tensor interpTo(const Tensor& x, const Tensor& like) {
        Tensor out = alloc(like.shape[0], like.shape[1], x.shape[2]);
        nn::resize_bilinear(out, x, /*align_corners=*/true);
        return out;
    }
};

}  // namespace

void Model::decode(const Features& f, Tensor* out_depth, Tensor* out_normal,
                   Tensor* out_kappa) {
    Decoder d(*this);
    const Hparams& h = weights.hparams();
    const int64_t gh = f.gh, gw = f.gw;
    const int64_t g7h = gh * 2, g7w = gw * 2;
    const int64_t g4h = scaled(gh, 3.5), g4w = scaled(gw, 3.5);

    // ---- token2feature: 1/14, 1/14, 1/7, 1/4 -----------------------------
    Tensor x3 = d.token2feature(f.tokens[3], 3, gh, gw, 1.0, h.feat_ch[3]);
    Tensor x2 = d.token2feature(f.tokens[2], 2, gh, gw, 1.0, h.feat_ch[2]);
    Tensor x1 = d.token2feature(f.tokens[1], 1, gh, gw, 2.0, h.feat_ch[1]);
    Tensor x0 = d.token2feature(f.tokens[0], 0, gh, gw, 3.5, h.feat_ch[0]);
    if (dump_enabled()) {
        dump_tensor("read3", x3, {gh, gw, h.feat_ch[3]});
        dump_tensor("read2", x2, {gh, gw, h.feat_ch[2]});
        dump_tensor("read1", x1, {g7h, g7w, h.feat_ch[1]});
        dump_tensor("read0", x0, {g4h, g4w, h.feat_ch[0]});
    }

    // ---- decoder_mono: coarse to 1/4 -------------------------------------
    // Allocated before the scope so the three FuseBlocks' working maps -- which
    // are the peak of the whole decoder -- are rewound while this survives.
    Tensor ref = nn::arena_tensor(arena, DType::F32, g4h, g4w, h.dec_ch[1]);
    {
        vk::ArenaScope scope(arena);
        Tensor t = d.fuseBlock(x3, Tensor(), "decoder.decoder_mono.upconv_3", 0.0,
                               h.dec_ch[3]);
        if (dump_enabled()) dump_tensor("upconv3", t, {gh, gw, h.dec_ch[3]});
        t = d.fuseBlock(t, x2, "decoder.decoder_mono.upconv_2", 2.0, h.dec_ch[2]);
        if (dump_enabled()) dump_tensor("upconv2", t, {g7h, g7w, h.dec_ch[2]});
        t = d.fuseBlock(t, x1, "decoder.decoder_mono.upconv_1", 1.75, h.dec_ch[1]);
        nn::copy(ref, t);
    }
    if (dump_enabled()) dump_tensor("ref_feat", ref, {g4h, g4w, h.dec_ch[1]});

    // ---- the two heads over the shared feature map ------------------------
    const int64_t used = h.used_res_channel;
    const int64_t px4 = g4h * g4w;
    Tensor state = nn::arena_tensor(arena, DType::F32, g4h, g4w, 6);
    {
        vk::ArenaScope scope(arena);
        // ref's leading `used` channels are the shared map; the last two are
        // the depth and normal confidences, which pack_state reads in place.
        Tensor feat = d.alloc(g4h, g4w, used);
        nn::strided_copy(feat, ref, px4, used, h.dec_ch[1], used);

        Tensor prob = d.alloc(g4h, g4w, h.n_bins);
        {
            vk::ArenaScope inner(arena);
            Tensor t = d.alloc(g4h, g4w, h.n_bins);
            d.conv3x3(t, feat, "decoder.depth_regressor.0", Act::Relu);
            d.conv1x1(prob, t, "decoder.depth_regressor.2", Act::None);
        }
        nn::softmax_rows(Tensor(prob.ptr, DType::F32, px4, h.n_bins),
                         Tensor(prob.ptr, DType::F32, px4, h.n_bins));

        // The expectation over the anchors is a [px, bins] x [1, bins]^T
        // matmul, so the tuned GEMM does it and no reduction kernel is needed.
        Tensor bins = nn::arena_tensor(arena, DType::F32, 1, h.n_bins);
        nn::tensor_from_host(bins, weights.depthBins().data(), h.n_bins);
        Tensor expect = d.alloc(g4h, g4w, 1);
        nn::matmul_nt(Tensor(expect.ptr, DType::F32, px4, 1),
                      Tensor(prob.ptr, DType::F32, px4, h.n_bins), bins);

        Tensor normal = d.alloc(g4h, g4w, 3);
        {
            vk::ArenaScope inner(arena);
            Tensor a = d.alloc(g4h, g4w, 128);
            Tensor b = d.alloc(g4h, g4w, 128);
            d.conv3x3(a, feat, "decoder.normal_predictor.0", Act::Relu);
            d.conv1x1(b, a, "decoder.normal_predictor.2", Act::Relu);
            d.conv1x1(a, b, "decoder.normal_predictor.4", Act::Relu);
            d.conv1x1(normal, a, "decoder.normal_predictor.6", Act::None);
        }

        struct {
            uint64_t out, expect, conf, normal;
            uint32_t n, conf_ch;
            float min_depth, max_depth, regress_scale;
            uint32_t groups_per_row;
        } pp{};
        pp.out = state.ptr;
        pp.expect = expect.ptr;
        pp.conf = ref.ptr;
        pp.normal = normal.ptr;
        pp.n = (uint32_t)px4;
        pp.conf_ch = (uint32_t)h.dec_ch[1];
        pp.min_depth = h.min_depth;
        pp.max_depth = h.max_depth;
        pp.regress_scale = Hparams::kRegressScale;
        vk::SpecList spec{0u, 0u};
        vk::Stream::get().dispatchFlat("metric3d.pack_state", spec, px4, 256, &pp,
                                       sizeof(pp), &pp.groups_per_row);
    }

    if (dump_enabled()) dump_tensor("state0", state, {g4h, g4w, 6});

    // ---- context: hidden state and the pre-convolved GRU inputs -----------
    // The three scales are fed the feature maps in the reverse of the DPT's
    // order, so `net[0]` is the 1/4 map and `net[2]` the 1/14 one.
    const Tensor scales[3] = {x0, x1, x2};
    const char* names[3] = {"outputs04", "outputs08", "outputs16"};
    Tensor net[3], ctx[3][3];
    for (int i = 0; i < 3; ++i) {
        const int64_t H = scales[i].shape[0], W = scales[i].shape[1];
        const std::string base = "decoder.context_feature_encoder." + std::string(names[i]);
        net[i] = nn::arena_tensor(arena, DType::F32, H, W, h.hidden);
        for (int g = 0; g < 3; ++g)
            ctx[i][g] = nn::arena_tensor(arena, DType::F32, H, W, h.hidden);

        {
            vk::ArenaScope scope(arena);
            Tensor t = d.residualBlock(scales[i], base + ".0.0", h.hidden);
            d.conv3x3(net[i], t, base + ".0.1", Act::None);
        }
        nn::unary(net[i], net[i], Act::Tanh);
        if (dump_enabled())
            dump_tensor(("net" + std::to_string(i)).c_str(), net[i], {H, W, h.hidden});

        // One 3*hidden convolution per scale, split once into the three gate
        // biases -- the reference hoists it out of the loop for the same
        // reason, it is `iters` repetitions of work that does not change.
        vk::ArenaScope scope(arena);
        Tensor inp = d.alloc(H, W, h.hidden);
        {
            vk::ArenaScope inner(arena);
            Tensor t = d.residualBlock(scales[i], base + ".1.0", h.hidden);
            d.conv3x3(inp, t, base + ".1.1", Act::None);
        }
        nn::unary(inp, inp, Act::Relu);

        Tensor zqr = d.alloc(H, W, 3 * h.hidden);
        d.conv3x3(zqr, inp, "decoder.context_zqr_convs." + std::to_string(i), Act::None);
        for (int g = 0; g < 3; ++g)
            nn::strided_copy(ctx[i][g],
                             Tensor(zqr.ptr + (uint64_t)g * h.hidden * 4, DType::F32,
                                    H * W, h.hidden),
                             H * W, h.hidden, 3 * h.hidden, h.hidden);
    }

    // ---- refinement -------------------------------------------------------
    Tensor mask = nn::arena_tensor(arena, DType::F32, g4h, g4w, 9 * h.up_factor * h.up_factor);
    for (int it = 0; it < h.iters; ++it) {
        // The slow-fast schedule: the coarse GRUs run ahead of the fine one, so
        // one iteration is three calls and only the last produces an update.
        auto gru32 = [&] {
            vk::ArenaScope scope(arena);
            d.convGru(net[2], "decoder.update_block.gru32", ctx[2][0], ctx[2][1],
                      ctx[2][2], {d.pool2x(net[1])});
        };
        auto gru16 = [&] {
            vk::ArenaScope scope(arena);
            Tensor a = d.interpTo(d.pool2x(net[0]), net[1]);
            Tensor b = d.interpTo(net[2], net[1]);
            d.convGru(net[1], "decoder.update_block.gru16", ctx[1][0], ctx[1][1],
                      ctx[1][2], {a, b});
        };
        gru32();
        gru32();
        gru16();
        gru32();
        gru16();
        {
            vk::ArenaScope scope(arena);
            Tensor a = d.interpTo(net[1], net[0]);
            d.convGru(net[0], "decoder.update_block.gru08", ctx[0][0], ctx[0][1],
                      ctx[0][2], {state, a});
        }

        vk::ArenaScope scope(arena);
        {
            // FlowHead is two independent branches over the same hidden state,
            // concatenated as (depth, depth confidence | normal xyz, kappa) --
            // which is the channel order the six-channel state already has, so
            // the concat is two strided copies into it.
            vk::ArenaScope inner(arena);
            Tensor t = d.alloc(g4h, g4w, h.hidden);
            Tensor dd = d.alloc(g4h, g4w, 2);
            Tensor dn = d.alloc(g4h, g4w, 4);
            d.conv3x3(t, net[0], "decoder.update_block.flow_head.conv1d", Act::Relu);
            d.conv3x3(dd, t, "decoder.update_block.flow_head.conv2d", Act::None);
            d.conv3x3(t, net[0], "decoder.update_block.flow_head.conv1n", Act::Relu);
            d.conv3x3(dn, t, "decoder.update_block.flow_head.conv2n", Act::None);

            Tensor packed = d.alloc(g4h, g4w, 6);
            nn::strided_copy(packed, dd, px4, 2, 2, 6);
            nn::strided_copy(Tensor(packed.ptr + 2 * 4, DType::F32, px4, 4), dn, px4, 4,
                             4, 6);
            nn::add(state, state, packed);
        }

        {
            vk::ArenaScope inner(arena);
            Tensor t = d.alloc(g4h, g4w, h.hidden);
            d.conv3x3(t, net[0], "decoder.update_block.mask.0", Act::Relu);
            d.conv1x1(mask, t, "decoder.update_block.mask.2", Act::None);
            // The reference scales the mask before its softmax, and softmax is
            // not scale-invariant, so this is part of the arithmetic.
            nn::unary(mask, mask, Act::None, 0.25f);
        }
    }

    // ---- convex upsample and unpack --------------------------------------
    const int64_t f4 = h.up_factor;
    const int64_t Ho = g4h * f4, Wo = g4w * f4;
    Tensor up = nn::arena_tensor(arena, DType::F32, Ho, Wo, 6);
    {
        struct {
            uint64_t out, flow, mask;
            uint32_t H, W, D, groups_per_row;
        } cp{};
        cp.out = up.ptr;
        cp.flow = state.ptr;
        cp.mask = mask.ptr;
        cp.H = (uint32_t)g4h;
        cp.W = (uint32_t)g4w;
        cp.D = 6;
        vk::SpecList spec{0u, (uint32_t)f4};
        vk::Stream::get().dispatchFlat("metric3d.convex_upsample", spec, Ho * Wo, 256, &cp,
                                       sizeof(cp), &cp.groups_per_row);
    }

    if (out_depth) *out_depth = nn::arena_tensor(arena, DType::F32, Ho, Wo, 1);
    if (out_normal) *out_normal = nn::arena_tensor(arena, DType::F32, Ho, Wo, 3);
    if (out_kappa) *out_kappa = nn::arena_tensor(arena, DType::F32, Ho, Wo, 1);
    {
        // The kernel writes all three; a caller that wants only some still gets
        // one pass, and the scratch for the rest is the arena's cheapest kind.
        Tensor scratch_d, scratch_n, scratch_k;
        if (!out_depth) scratch_d = nn::arena_tensor(arena, DType::F32, Ho, Wo, 1);
        if (!out_normal) scratch_n = nn::arena_tensor(arena, DType::F32, Ho, Wo, 3);
        if (!out_kappa) scratch_k = nn::arena_tensor(arena, DType::F32, Ho, Wo, 1);
        struct {
            uint64_t depth, normal, kappa, state;
            uint32_t n;
            float min_depth, max_depth, regress_scale;
            uint32_t groups_per_row;
        } up_p{};
        up_p.depth = out_depth ? out_depth->ptr : scratch_d.ptr;
        up_p.normal = out_normal ? out_normal->ptr : scratch_n.ptr;
        up_p.kappa = out_kappa ? out_kappa->ptr : scratch_k.ptr;
        up_p.state = up.ptr;
        up_p.n = (uint32_t)(Ho * Wo);
        up_p.min_depth = h.min_depth;
        up_p.max_depth = h.max_depth;
        up_p.regress_scale = Hparams::kRegressScale;
        vk::SpecList spec{0u, 0u};
        vk::Stream::get().dispatchFlat("metric3d.unpack_state", spec, Ho * Wo, 256, &up_p,
                                       sizeof(up_p), &up_p.groups_per_row);
    }

    if (dump_enabled()) {
        if (out_depth) dump_tensor("depth", *out_depth, {Ho, Wo});
        if (out_normal) dump_tensor("normal", *out_normal, {Ho, Wo, 3});
        if (out_kappa) dump_tensor("normal_conf", *out_kappa, {Ho, Wo});
    }
}

}  // namespace metric3d
