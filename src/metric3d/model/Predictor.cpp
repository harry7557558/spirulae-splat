// The public surface: a checkpoint, an arena, and one image at a time.

#include "metric3d/Metric3D.h"

#include "metric3d/Common.h"
#include "metric3d/model/Fetch.h"
#include "metric3d/model/Model.h"
#include "nn/Ops.h"
#include "nn/vk/EmbeddedSpirv.h"
#include "nn/vk/Stream.h"

#include "core/Env.h"
#include "external/npy.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <vector>

// This library's SPIR-V blobs. Registered by an explicit call from load(), not
// by a static initializer: ss_metric3d is a static archive, and an object
// nothing references is not linked -- the kernels would come back "no shader
// module". See nn/vk/EmbeddedSpirv.h.
NN_DECLARE_EMBEDDED_MODULES(metric3d)

namespace metric3d {
namespace {

using nn::DType;
using nn::Tensor;

// Resolved and created once: every dump after the first would otherwise race
// to make the same directory, and npy::write_npy throws rather than creating
// one for itself.
const char* dump_dir() {
    static const char* dir = [] {
        const char* d = spirula::env("METRIC3D_DUMP");
        if (d) {
            std::error_code ec;
            std::filesystem::create_directories(d, ec);
        }
        return d;
    }();
    return dir;
}

void write_npy(const char* name, const float* data, const std::vector<int64_t>& shape) {
    npy::shape_t s;
    for (int64_t d : shape) s.push_back((unsigned long)d);
    npy::npy_data_ptr<float> d{data, s, false};
    npy::write_npy(std::string(dump_dir()) + "/" + name + ".npy", d);
}

}  // namespace

bool dump_enabled() { return dump_dir() != nullptr; }

void dump_tensor(const char* name, const Tensor& t, const std::vector<int64_t>& shape) {
    if (!dump_dir()) return;
    std::vector<float> host((size_t)t.numel());
    nn::tensor_to_host(t, host.data(), t.numel());
    write_npy(name, host.data(), shape);
}

struct Predictor::Impl {
    Model model;
};

Predictor::Predictor() : impl_(new Impl) {}
Predictor::~Predictor() { delete impl_; }

int Predictor::patchSize() { return 14; }
int Predictor::sizeGranularity() { return 28; }

void Predictor::load(const std::string& model) {
    NN_ENSURE_EMBEDDED_MODULES(metric3d);
    impl_->model.weights.load(resolve_model(model));
}

bool Predictor::loaded() const { return impl_->model.weights.loaded(); }

uint64_t Predictor::deviceBytes() const { return impl_->model.weights.deviceBytes(); }

uint64_t Model::planArenaBytes(int64_t H, int64_t W) const {
    const Hparams& h = hp();
    const int64_t gh = H / h.patch, gw = W / h.patch;
    const int64_t p14 = gh * gw;
    const int64_t p7 = 4 * p14;
    const int64_t p4 = (int64_t)std::floor(gh * 3.5) * (int64_t)std::floor(gw * 3.5);
    const int64_t pfull = p4 * h.up_factor * h.up_factor;
    const int64_t ntok = p14 + h.prefix_tokens();
    const int64_t D = h.embed_dim;

    // The token sequence outlives the encoder: it is what the decoder reads.
    // One tensor when the export ends in a LayerNorm (all four reads share it),
    // otherwise the running state plus four taps.
    const int64_t tokens = ntok * D * (h.final_norm ? 2 : 5);
    // The widest thing live inside a block: the fused qkv, the attention output
    // and the MLP's hidden layer.
    const int64_t enc_tmp =
        ntok * (3 * D + D + (h.swiglu ? 4 * (int64_t)h.mlp_hidden : (int64_t)h.mlp_hidden));

    const int64_t dec_live =
        2 * p14 * D + p7 * h.feat_ch[1] + p4 * h.feat_ch[0]   // the four reads
        + p4 * h.dec_ch[1]                                     // ref_feat
        + p4 * (6 + 9 * h.up_factor * h.up_factor)             // state + mask
        + 4 * (p4 + p7 + p14) * h.hidden                       // net + the three gates
        + pfull * (6 + 5);                                     // upsampled + outputs
    // The deepest transient is the FuseBlock chain's upsample stage; the depth
    // regressor's 256-bin map is the runner-up and they are not live together.
    const int64_t dec_tmp =
        std::max<int64_t>(4 * p7 * h.dec_ch[2] + p4 * (h.dec_ch[2] + h.dec_ch[1]),
                          2 * p4 * h.n_bins + p4 * h.used_res_channel);

    const int64_t peak = tokens + std::max(enc_tmp, dec_live + dec_tmp);
    // nn::conv2d chunks its column workspace to ~32 MiB, the GRU concatenations
    // hold one more copy of their inputs, and flash-decoding takes its partials
    // from the same arena.
    return (uint64_t)peak * 4 + (192ull << 20);
}

Prediction Predictor::predict(const float* rgb, int width, int height,
                              const PredictOptions& opts) {
    Model& m = impl_->model;
    NN_CHECK(m.weights.loaded(), "metric3d::Predictor::predict before load()");
    NN_CHECK(rgb != nullptr && width > 0 && height > 0,
             "metric3d::Predictor::predict: empty image");
    const Hparams& h = m.weights.hparams();

    const int64_t P = h.patch;
    const int64_t H = (height / P) * P, W = (width / P) * P;
    NN_CHECK(H >= 2 * P && W >= 2 * P,
             "%dx%d is too small for a stride-%lld patch embedding", width, height,
             (long long)P);

    m.arena.reserve(m.planArenaBytes(H, W));
    vk::ArenaScope root(m.arena);

    // The graph starts with (pixel_values - mean) / std over 0..255 values; the
    // caller works in 0..1, so the 255 rides along here. One pass over the
    // image next to a backbone that is hundreds of them.
    Tensor image = nn::arena_tensor(m.arena, DType::F32, H, W, 3);
    {
        std::vector<float> host((size_t)H * W * 3);
        const float sc[3] = {255.0f / h.rgb_std[0], 255.0f / h.rgb_std[1],
                             255.0f / h.rgb_std[2]};
        const float off[3] = {-h.rgb_mean[0] / h.rgb_std[0], -h.rgb_mean[1] / h.rgb_std[1],
                              -h.rgb_mean[2] / h.rgb_std[2]};
        nn::parallel_for(H, [&](int64_t y0, int64_t y1) {
            for (int64_t y = y0; y < y1; ++y) {
                const float* src = rgb + (size_t)y * width * 3;
                float* dst = host.data() + (size_t)y * W * 3;
                for (int64_t x = 0; x < W; ++x)
                    for (int c = 0; c < 3; ++c)
                        dst[x * 3 + c] = src[x * 3 + c] * sc[c] + off[c];
            }
        });
        nn::tensor_from_host(image, host.data(), (int64_t)host.size());
        if (dump_enabled()) {
            // `pixel_values` is BEFORE the normalization: the graph carries the
            // mean and standard deviation as its first two nodes, so handing
            // the reference what we hand the patch embedding would normalize
            // it twice.
            std::vector<float> chw((size_t)H * W * 3);
            for (int64_t y = 0; y < H; ++y)
                for (int64_t x = 0; x < W; ++x)
                    for (int c = 0; c < 3; ++c)
                        chw[((size_t)c * H + y) * W + x] =
                            rgb[((size_t)y * width + x) * 3 + c] * 255.0f;
            write_npy("input", chw.data(), {1, 3, H, W});
        }
    }

    const Features f = m.encode(image);
    Tensor depth, normal, kappa;
    m.decode(f, opts.want_depth ? &depth : nullptr, opts.want_normal ? &normal : nullptr,
             opts.want_normal ? &kappa : nullptr);

    Prediction out;
    const Tensor any = depth.valid() ? depth : normal;
    out.height = (int)any.shape[0];
    out.width = (int)any.shape[1];
    if (depth.valid()) {
        out.depth.resize((size_t)depth.numel());
        nn::tensor_to_host(depth, out.depth.data(), depth.numel());
    }
    if (normal.valid()) {
        out.normal.resize((size_t)normal.numel());
        nn::tensor_to_host(normal, out.normal.data(), normal.numel());
        out.confidence.resize((size_t)kappa.numel());
        nn::tensor_to_host(kappa, out.confidence.data(), kappa.numel());
    }
    return out;
}

}  // namespace metric3d
