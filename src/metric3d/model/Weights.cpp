#include "metric3d/model/Weights.h"

#include "metric3d/Common.h"

#include "core/Env.h"
#include "nn/vk/Memory.h"
#include "nn/vk/Stream.h"

#include <algorithm>
#include <cstdarg>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>

namespace metric3d {
namespace {

using nn::DType;
using nn::Tensor;

constexpr uint64_t kAlign = 256;
// NVIDIA's Windows driver refuses a single 1 GiB device-local allocation on an
// 8 GiB laptop part with 7 GiB free; giant2 is 2.75 GB of weights.
constexpr uint64_t kChunkBytes = 256ull << 20;
constexpr const char* kPrefix = "meta_arch.depth_model.";

uint64_t align_up(uint64_t v, uint64_t a) { return (v + a - 1) / a * a; }

bool starts_with(const std::string& s, const std::string& p) {
    return s.size() >= p.size() && s.compare(0, p.size(), p) == 0;
}
bool ends_with(const std::string& s, const std::string& p) {
    return s.size() >= p.size() && s.compare(s.size() - p.size(), p.size(), p) == 0;
}

bool replace_first(std::string& s, const char* from, const char* to) {
    const size_t at = s.find(from);
    if (at == std::string::npos) return false;
    s.replace(at, std::strlen(from), to);
    return true;
}

// `meta_arch.depth_model.encoder.blocks.0.23.norm1.weight`
//   -> `encoder.blocks.23.norm1.weight`
// `/depth_model/encoder/blocks.0/blocks.0.23/attn/qkv`
//   -> `encoder.blocks.23.attn.qkv`
//
// Flattening the BlockChunk is only sound because there is exactly one of them;
// load() checks that before calling this on anything.
std::string normalize(std::string s) {
    if (starts_with(s, kPrefix)) s.erase(0, std::strlen(kPrefix));
    if (starts_with(s, "/depth_model/")) s.erase(0, std::strlen("/depth_model/"));
    for (char& c : s)
        if (c == '/') c = '.';
    // A node path names the chunk twice ("blocks.0/blocks.0.23"), an
    // initializer once ("blocks.0.23"). Trying both in sequence would strip a
    // node path's real block index along with its chunk.
    if (!replace_first(s, "blocks.0.blocks.0.", "blocks."))
        replace_first(s, "blocks.0.", "blocks.");
    return s;
}

// A matrix headed for nn::linear or nn::conv2d goes to the device as f16 --
// that is the operand the tuned GEMM wants, and the published exports are
// already f16 so nothing is lost.  Everything else stays f32.
//
// SS_METRIC3D_F32_WEIGHTS=1 keeps them all f32. That is not a quality knob:
// it costs twice the VRAM and gives up the tensor-core GEMM for a difference
// no output shows. It exists because it is what makes
// tools/metric3d/compare_ort.py tight -- against an fp32 export the rounding
// here is the entire disagreement, and a real bug of the same size would hide
// under it.
bool wants_f16(const std::string& name, size_t rank) {
    static const bool force_f32 = [] {
        const char* v = spirula::env("METRIC3D_F32_WEIGHTS");
        return v && v[0] && v[0] != '0';
    }();
    return !force_f32 && rank >= 2 && ends_with(name, ".weight");
}

struct Staged {
    std::string          name;
    std::vector<int64_t> shape;
    std::vector<float>   data;
    bool                 f16 = false;
};

// [in, out] -> [out, in]. Every Linear in the graph arrives as the `W` of
// `x @ W`, which is the transpose of what nn::linear indexes.
void transpose2d(Staged& s) {
    const int64_t in = s.shape[0], out = s.shape[1];
    std::vector<float> t((size_t)in * out);
    for (int64_t i = 0; i < in; ++i)
        for (int64_t o = 0; o < out; ++o) t[(size_t)(o * in + i)] = s.data[(size_t)(i * out + o)];
    s.data.swap(t);
    s.shape = {out, in};
}

// ConvTranspose2d's [Cin, Cout, 2, 2] -> the [Cout*4, Cin] matrix
// nn::conv_transpose2x2 multiplies by, with the four kernel taps as four
// output-channel groups.
void repack_conv_transpose(Staged& s) {
    const int64_t cin = s.shape[0], cout = s.shape[1];
    NN_CHECK(s.shape[2] == 2 && s.shape[3] == 2,
             "'%s' is a %lldx%lld transposed convolution; only 2x2 stride 2 is "
             "implemented",
             s.name.c_str(), (long long)s.shape[2], (long long)s.shape[3]);
    std::vector<float> dst((size_t)cin * cout * 4);
    for (int64_t ci = 0; ci < cin; ++ci)
        for (int64_t co = 0; co < cout; ++co)
            for (int64_t i = 0; i < 2; ++i)
                for (int64_t j = 0; j < 2; ++j)
                    dst[(size_t)((co * 4 + i * 2 + j) * cin + ci)] =
                        s.data[(size_t)(((ci * cout + co) * 2 + i) * 2 + j)];
    s.data.swap(dst);
    s.shape = {cout * 4, cin};
}

}  // namespace

Weights::~Weights() {
    for (nn::DevicePtr p : chunks_) vk::device_free(p);
}

void Weights::load(const std::string& onnx_path) {
    NN_CHECK(!loaded_, "metric3d::Weights::load called twice");
    const char* path = onnx_path.c_str();

    // ---- pass 1: structure ------------------------------------------------
    // Two things live only in the graph and not in any tensor: which module
    // each anonymous MatMul weight belongs to, and how many refinement
    // iterations the export unrolled (the GRU weights are shared across them,
    // so the count is invisible in the weight list -- small runs 4 where large
    // and giant2 run 8).
    const OnnxFile graph = nn::read_onnx_structure(onnx_path);

    std::unordered_map<std::string, std::string> matmul_module;  // init name -> module
    int n_update_blocks = 0, max_chunk = 0;
    bool has_encoder_norm = false;
    {
        std::vector<std::string> update_groups;
        for (const OnnxNode& n : graph.nodes) {
            if (n.op_type == "MatMul" && n.inputs.size() == 2 &&
                ends_with(n.name, "/MatMul")) {
                std::string mod = n.name.substr(0, n.name.size() - 7);
                matmul_module[n.inputs[1]] = normalize(mod);
            }
            if (n.name.find("/encoder/norm/") != std::string::npos) has_encoder_norm = true;
            // "/depth_model/decoder/update_block_7/..." -- the trailing index
            // is the unrolled call, three per iteration (the slow-fast schedule
            // runs the 1/16 and 1/8 GRUs on their own before the full step).
            const size_t at = n.name.find("/decoder/update_block");
            if (at != std::string::npos) {
                const size_t beg = at + std::strlen("/decoder/");
                const size_t end = n.name.find('/', beg);
                if (end != std::string::npos) {
                    std::string g = n.name.substr(beg, end - beg);
                    if (std::find(update_groups.begin(), update_groups.end(), g) ==
                        update_groups.end())
                        update_groups.push_back(g);
                }
            }
            const size_t bat = n.name.find("/blocks.");
            if (bat != std::string::npos) {
                const int chunk = std::atoi(n.name.c_str() + bat + 8);
                max_chunk = std::max(max_chunk, chunk);
            }
        }
        n_update_blocks = (int)update_groups.size();
    }
    NN_CHECK(max_chunk == 0,
             "'%s' has %d encoder block chunks; this reader flattens a single "
             "chunk and every published Metric3D export has one",
             path, max_chunk + 1);
    NN_CHECK(n_update_blocks > 0 && n_update_blocks % 3 == 0,
             "'%s' unrolls %d update blocks; expected a multiple of 3 (one "
             "refinement iteration is three GRU calls)",
             path, n_update_blocks);
    hp_.iters = n_update_blocks / 3;
    hp_.final_norm = has_encoder_norm;

    // ---- pass 2: initializers, one at a time ------------------------------
    // Streamed rather than collected: giant2's f16 weights are 2.75 GB on disk
    // and would be 5.5 GB held as host floats all at once.
    std::vector<Staged> staged;
    std::vector<float> rgb_mean, rgb_std;
    int n_skipped = 0;
    nn::read_onnx_initializers(onnx_path, [&](OnnxTensor&& t) {
        const size_t rank = t.shape.size();
        if (t.name == "rgb_mean") { rgb_mean = std::move(t.data); return; }
        if (t.name == "rgb_std")  { rgb_std = std::move(t.data); return; }

        const bool named = starts_with(t.name, kPrefix);
        auto mm = matmul_module.find(t.name);

        if (!named && mm == matmul_module.end()) {
            // The graph's own constants. Three of them are weights in
            // everything but name; the rest are shapes, axes and scalars the
            // forward pass spells out itself.
            if (rank == 4 && t.shape[0] == 1 && t.shape[2] == t.shape[3] &&
                t.shape[2] > 1) {
                // pos_embed, already transposed to [1, C, g, g] for the graph's
                // own bicubic resize. Held as [g*g, C] to match token order.
                pos_grid_ = (int)t.shape[2];
                const int64_t C = t.shape[1], G = t.shape[2] * t.shape[3];
                patch_pos_.resize((size_t)G * C);
                for (int64_t c = 0; c < C; ++c)
                    for (int64_t g = 0; g < G; ++g)
                        patch_pos_[(size_t)(g * C + c)] = t.data[(size_t)(c * G + g)];
            } else if (rank == 3 && t.shape[0] == 1 && t.shape[1] == 1) {
                cls_pos_ = std::move(t.data);
            } else if (rank == 2 && t.shape[0] == 1 && t.shape[1] >= 64) {
                bins_ = std::move(t.data);
            } else {
                ++n_skipped;
            }
            return;
        }

        Staged s;
        s.name = (mm != matmul_module.end()) ? mm->second + ".weight" : normalize(t.name);
        s.shape = std::move(t.shape);
        s.data = std::move(t.data);
        if (mm != matmul_module.end()) {
            NN_CHECK(s.shape.size() == 2, "'%s': linear weight '%s' has rank %zu", path,
                     s.name.c_str(), s.shape.size());
            transpose2d(s);
        } else if (s.shape.size() == 4 && s.name.find(".sample.weight") != std::string::npos) {
            repack_conv_transpose(s);
        }
        s.f16 = wants_f16(s.name, s.shape.size());
        staged.push_back(std::move(s));
    });
    NN_CHECK(!staged.empty(), "'%s' holds no Metric3D weights", path);
    NN_CHECK(!patch_pos_.empty() && !cls_pos_.empty(),
             "'%s' has no positional embedding constant", path);
    NN_CHECK(!bins_.empty(), "'%s' has no depth anchor constant", path);
    for (size_t i = 1; i < bins_.size(); ++i)
        NN_CHECK(bins_[i] > bins_[i - 1], "'%s': depth anchors are not increasing", path);
    hp_.min_depth = bins_.front();
    hp_.max_depth = bins_.back();
    hp_.n_bins = (int)bins_.size();
    if (rgb_mean.size() == 3 && rgb_std.size() == 3)
        for (int i = 0; i < 3; ++i) {
            hp_.rgb_mean[i] = rgb_mean[(size_t)i];
            hp_.rgb_std[i] = rgb_std[(size_t)i];
        }

    // ---- upload, chunked --------------------------------------------------
    uint64_t cursor = 0, chunk_used = 0;
    std::vector<uint64_t> chunk_sizes;
    std::vector<std::pair<uint32_t, uint64_t>> placement(staged.size());
    for (size_t i = 0; i < staged.size(); ++i) {
        const uint64_t bytes =
            (uint64_t)staged[i].data.size() * (staged[i].f16 ? 2u : 4u);
        chunk_used = align_up(chunk_used, kAlign);
        if (chunk_used != 0 && chunk_used + bytes > kChunkBytes) {
            chunk_sizes.push_back(chunk_used);
            chunk_used = 0;
        }
        placement[i] = {(uint32_t)chunk_sizes.size(), chunk_used};
        chunk_used += bytes;
        cursor += bytes;
    }
    if (chunk_used) chunk_sizes.push_back(chunk_used);
    for (size_t i = 0; i < chunk_sizes.size(); ++i)
        chunks_.push_back(vk::device_alloc(chunk_sizes[i], "metric3d-weights"));
    device_bytes_ = cursor;

    for (size_t i = 0; i < staged.size(); ++i) {
        Staged& s = staged[i];
        NN_CHECK(s.shape.size() <= 4, "'%s': '%s' has rank %zu; nn::Tensor holds 4", path,
                 s.name.c_str(), s.shape.size());
        Tensor t;
        t.ptr = chunks_[placement[i].first] + placement[i].second;
        t.dtype = s.f16 ? DType::F16 : DType::F32;
        t.ndim = (int32_t)std::max<size_t>(s.shape.size(), 1);
        for (size_t d = 0; d < s.shape.size(); ++d) t.shape[d] = s.shape[d];
        if (s.shape.empty()) t.shape[0] = 1;
        nn::tensor_from_host(t, s.data.data(), (int64_t)s.data.size());
        tensors_[s.name] = t;
        std::vector<float>().swap(s.data);
    }
    vk::Stream::get().sync();

    // ---- hyperparameters, from shapes -------------------------------------
    {
        const Tensor pe = get("encoder.patch_embed.proj.weight");
        hp_.embed_dim = (int)pe.shape[0];
        hp_.patch = (int)pe.shape[2];
        NN_CHECK(pe.shape[2] == pe.shape[3], "'%s': non-square patch embedding", path);
        hp_.num_register_tokens = (int)get("encoder.register_tokens").shape[1];
        // DINOv2 keeps head_dim at 64 for every width (384/6, 1024/16,
        // 1536/24) and the export records the head count nowhere.
        NN_CHECK(hp_.embed_dim % 64 == 0, "'%s': embed_dim %d is not a multiple of 64",
                 path, hp_.embed_dim);
        hp_.num_heads = hp_.embed_dim / 64;
        hp_.swiglu = has("encoder.blocks.0.mlp.w12.weight");
        hp_.mlp_hidden = hp_.swiglu ? (int)get("encoder.blocks.0.mlp.w12.weight").shape[0] / 2
                                    : (int)get("encoder.blocks.0.mlp.fc1.weight").shape[0];

        hp_.depth = 0;
        while (has("encoder.blocks." + std::to_string(hp_.depth) + ".norm1.weight"))
            ++hp_.depth;
        NN_CHECK(hp_.depth > 0, "'%s' has no encoder blocks", path);
        if (!hp_.final_norm) {
            NN_CHECK(hp_.depth % 4 == 0,
                     "'%s' taps four intermediate blocks but has %d of them", path,
                     hp_.depth);
            for (int i = 0; i < 4; ++i) hp_.taps[i] = hp_.depth / 4 * (i + 1) - 1;
        }

        hp_.feat_ch[0] = (int)get("decoder.token2feature.read_0.sample.0.weight").shape[0];
        // read_1's ConvTranspose was repacked to [Cout*4, Cin] above.
        hp_.feat_ch[1] = (int)get("decoder.token2feature.read_1.sample.weight").shape[0] / 4;
        hp_.feat_ch[2] = hp_.feat_ch[3] = hp_.embed_dim;

        hp_.dec_ch[3] = (int)get("decoder.decoder_mono.upconv_3.out_conv.weight").shape[0];
        hp_.dec_ch[2] = (int)get("decoder.decoder_mono.upconv_2.out_conv.weight").shape[0];
        // upconv_1 emits the shared feature map plus one confidence channel for
        // each of depth and normal.
        hp_.dec_ch[1] = (int)get("decoder.decoder_mono.upconv_1.out_conv.weight").shape[0];
        hp_.used_res_channel = hp_.dec_ch[1] - 2;
        hp_.hidden = (int)get("decoder.context_zqr_convs.0.weight").shape[0] / 3;

        const int64_t mask_ch = get("decoder.update_block.mask.2.weight").shape[0];
        NN_CHECK(mask_ch % 9 == 0, "'%s': convex upsample mask has %lld channels", path,
                 (long long)mask_ch);
        const int64_t f2 = mask_ch / 9;
        hp_.up_factor = (int)std::lround(std::sqrt((double)f2));
        NN_CHECK((int64_t)hp_.up_factor * hp_.up_factor == f2,
                 "'%s': convex upsample factor is not square", path);

        NN_CHECK((int)get("decoder.depth_regressor.2.weight").shape[0] == hp_.n_bins,
                 "'%s': %d depth anchors but the regressor emits %lld", path, hp_.n_bins,
                 (long long)get("decoder.depth_regressor.2.weight").shape[0]);
        NN_CHECK((int)patch_pos_.size() == pos_grid_ * pos_grid_ * hp_.embed_dim,
                 "'%s': positional embedding is %zu for a %dx%d grid at width %d", path,
                 patch_pos_.size(), pos_grid_, pos_grid_, hp_.embed_dim);
    }

    path_ = onnx_path;
    loaded_ = true;
    NN_LOG_INFO("[metric3d] %s: dim=%d depth=%d heads=%d%s, hidden=%d iters=%d up=%d, "
                "%zu tensors, %.2f MB on device in %zu chunks (%d graph constants "
                "skipped)\n",
                onnx_path.c_str(), hp_.embed_dim, hp_.depth, hp_.num_heads,
                hp_.swiglu ? " swiglu" : "", hp_.hidden, hp_.iters, hp_.up_factor,
                tensors_.size(), (double)device_bytes_ / 1e6, chunks_.size(), n_skipped);
}

nn::Tensor Weights::get(const std::string& name) const {
    auto it = tensors_.find(name);
    if (it != tensors_.end()) return it->second;

    // Near misses: the failure is almost always a name that drifted, and
    // printing the neighbourhood turns a five-minute hunt into a glance.
    std::vector<std::string> near;
    const size_t dot = name.rfind('.');
    for (const auto& kv : tensors_)
        if (dot != std::string::npos && kv.first.compare(0, dot, name, 0, dot) == 0)
            near.push_back(kv.first);
    std::sort(near.begin(), near.end());
    std::string hint;
    for (size_t i = 0; i < near.size() && i < 8; ++i) hint += "\n    " + near[i];
    nn::fail("no weight named '%s'%s%s", name.c_str(),
             hint.empty() ? "" : " (did you mean one of these?)", hint.c_str());
}

nn::Tensor Weights::getf(const char* fmt, ...) const {
    char buf[256];
    va_list ap;
    va_start(ap, fmt);
    std::vsnprintf(buf, sizeof buf, fmt, ap);
    va_end(ap);
    return get(buf);
}

}  // namespace metric3d
