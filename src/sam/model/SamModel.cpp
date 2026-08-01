#include "sam/model/SamModel.h"

#include "sam/Common.h"
#include "nn/core/Error.h"
#include "nn/core/Log.h"
#include "nn/vk/Stream.h"

#define _USE_MATH_DEFINES
#include <cmath>
#include <cstring>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace sam {

// ================
// Host-side encodings
// ================

std::vector<float> sinusoidal_pe_2d(int H, int W, int d_model) {
    const int half = d_model / 2;
    const float scale = 2.0f * (float)M_PI;
    const float temperature = 10000.0f;
    std::vector<float> pe((size_t)d_model * H * W);
    for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x) {
            // The reference normalizes by (pos + 1) / max_pos, NOT (pos + 0.5).
            const float pos_y = ((float)(y + 1) / (float)H) * scale;
            const float pos_x = ((float)(x + 1) / (float)W) * scale;
            for (int i = 0; i < half; ++i) {
                // dim_t = temperature^(2*(i//2) / half): sin/cos alternate over
                // i while the frequency only changes every other channel.
                const int paired = i & ~1;
                const float dim_t = std::pow(temperature, (float)paired / (float)half);
                float vy, vx;
                if ((i % 2) == 0) {
                    vy = std::sin(pos_y / dim_t);
                    vx = std::sin(pos_x / dim_t);
                } else {
                    vy = std::cos(pos_y / dim_t);
                    vx = std::cos(pos_x / dim_t);
                }
                const size_t base = ((size_t)y * W + x) * d_model;
                pe[base + i] = vy;
                pe[base + half + i] = vx;
            }
        }
    return pe;
}

void sine_pe_1d(float* out, float pos, int dim, float temperature) {
    const int half = dim / 2;
    for (int i = 0; i < half; ++i) {
        const int paired = i & ~1;
        const float dim_t = std::pow(temperature, (float)paired / (float)half);
        const float v = pos / dim_t;
        out[i] = std::sin(v);
        out[half + i] = std::cos(v);
    }
}

std::vector<float> axial_rope_table(int dim, int end_x, int end_y, float theta) {
    // dim/4 complex frequencies per axis: the first half of the head's complex
    // slots rotate by x, the second half by y.
    const int quarter = dim / 4;
    std::vector<float> freqs((size_t)quarter);
    for (int i = 0; i < quarter; ++i)
        freqs[(size_t)i] = 1.0f / std::pow(theta, (float)(i * 4) / (float)dim);

    const int N = end_x * end_y;
    std::vector<float> out((size_t)N * dim);
    for (int idx = 0; idx < N; ++idx) {
        const float tx = (float)(idx % end_x);
        const float ty = (float)(idx / end_x);
        float* row = out.data() + (size_t)idx * dim;
        for (int i = 0; i < quarter; ++i) {
            const float ax = tx * freqs[(size_t)i];
            row[i * 2 + 0] = std::cos(ax);
            row[i * 2 + 1] = std::sin(ax);
        }
        for (int i = 0; i < quarter; ++i) {
            const float ay = ty * freqs[(size_t)i];
            row[quarter * 2 + i * 2 + 0] = std::cos(ay);
            row[quarter * 2 + i * 2 + 1] = std::sin(ay);
        }
    }
    return out;
}

void sine_encode_box(float* out, float cx, float cy, float w, float h, int num_pos_feats,
                     float temperature) {
    const float scale = 2.0f * (float)M_PI;
    const float x_embed = cx * scale;
    const float y_embed = cy * scale;
    for (int i = 0; i < num_pos_feats; ++i) {
        const int div_idx = 2 * (i / 2);
        const float dim_t = std::pow(temperature, (float)div_idx / (float)num_pos_feats);
        const float px = x_embed / dim_t;
        const float py = y_embed / dim_t;
        if ((i % 2) == 0) {
            out[i] = std::sin(py);
            out[num_pos_feats + i] = std::sin(px);
        } else {
            out[i] = std::cos(py);
            out[num_pos_feats + i] = std::cos(px);
        }
    }
    out[2 * num_pos_feats] = h;
    out[2 * num_pos_feats + 1] = w;
}

std::vector<float> bicubic_resize(const float* src, int Hi, int Wi, int Ho, int Wo) {
    // PyTorch's upsample_bicubic2d: Keys cubic convolution with A = -0.75,
    // align_corners=False source mapping, and clamped edge access. No
    // antialiasing -- F.interpolate does not apply any unless asked.
    auto weights = [](float t, float* w) {
        // Keys' A, declared inside: MSVC refuses to implicitly capture a
        // constexpr from the enclosing scope in a captureless lambda.
        constexpr float A = -0.75f;
        // ((A+2)x - (A+3))x^2 + 1 on [0,1]; ((Ax - 5A)x + 8A)x - 4A on [1,2].
        const float t1 = t, t2 = 1.0f - t;
        w[1] = ((A + 2.0f) * t1 - (A + 3.0f)) * t1 * t1 + 1.0f;
        w[2] = ((A + 2.0f) * t2 - (A + 3.0f)) * t2 * t2 + 1.0f;
        const float s1 = t + 1.0f, s2 = 2.0f - t;
        w[0] = ((A * s1 - 5.0f * A) * s1 + 8.0f * A) * s1 - 4.0f * A;
        w[3] = ((A * s2 - 5.0f * A) * s2 + 8.0f * A) * s2 - 4.0f * A;
    };
    auto clamp = [](int v, int lo, int hi) { return v < lo ? lo : (v > hi ? hi : v); };

    const float sy = (float)Hi / (float)Ho, sx = (float)Wi / (float)Wo;
    std::vector<float> out((size_t)Ho * Wo);
    for (int oy = 0; oy < Ho; ++oy) {
        const float fy = ((float)oy + 0.5f) * sy - 0.5f;
        const int y0 = (int)std::floor(fy);
        float wy[4];
        weights(fy - (float)y0, wy);
        for (int ox = 0; ox < Wo; ++ox) {
            const float fx = ((float)ox + 0.5f) * sx - 0.5f;
            const int x0 = (int)std::floor(fx);
            float wx[4];
            weights(fx - (float)x0, wx);
            float acc = 0.0f;
            for (int j = 0; j < 4; ++j) {
                const int sy_i = clamp(y0 - 1 + j, 0, Hi - 1);
                float row = 0.0f;
                for (int i = 0; i < 4; ++i)
                    row += wx[i] * src[(size_t)sy_i * Wi + clamp(x0 - 1 + i, 0, Wi - 1)];
                acc += wy[j] * row;
            }
            out[(size_t)oy * Wo + ox] = acc;
        }
    }
    return out;
}

void pe_encode_coord(float* out, float x_norm, float y_norm, const float* gaussian,
                     int num_pos_feats) {
    const float c[2] = {2.0f * x_norm - 1.0f, 2.0f * y_norm - 1.0f};
    for (int i = 0; i < num_pos_feats; ++i) {
        // `gaussian` is [2, num_pos_feats] row-major.
        float dot = c[0] * gaussian[i] + c[1] * gaussian[num_pos_feats + i];
        dot *= 2.0f * (float)M_PI;
        out[i] = std::sin(dot);
        out[i + num_pos_feats] = std::cos(dot);
    }
}

// ================
// Model build
// ================

namespace {

// Small pool-backed f32 tensor uploaded from host data.
nn::Tensor upload_table(vk::PoolSlot slot, uint32_t sub, const std::vector<float>& v,
                        int64_t d0, int64_t d1 = 1, int64_t d2 = 1) {
    nn::Tensor t = nn::pool_tensor(slot, sub, nn::DType::F32, d0, d1, d2);
    NN_CHECK(t.numel() == (int64_t)v.size(),
               "upload_table: shape holds %lld elements, data has %zu",
               (long long)t.numel(), v.size());
    nn::tensor_from_host(t, v.data(), (int64_t)v.size());
    return t;
}

// Absent tensors come back empty, not zero-filled: callers test emptiness to
// decide whether the feature exists at all (SAM 2.0 ships neither
// no_obj_embed_spatial nor obj_ptr_tpos_proj).
std::vector<float> download(const nn::Tensor& t) {
    if (!t.valid()) return {};
    std::vector<float> v((size_t)t.numel());
    nn::tensor_to_host(t, v.data(), t.numel());
    return v;
}

}  // namespace

void SamModel::load(const std::string& path, int img_size_override) {
    store_.load(path);
    const Hparams& h = hp();

    img_size_ = img_size_override > 0 ? img_size_override : h.img_size;
    NN_CHECK(img_size_ % h.patch_size == 0,
               "--img-size %d is not a multiple of the model's feature stride %d",
               img_size_, h.patch_size);
    grid_ = img_size_ / h.patch_size;
    if (img_size_ != h.img_size) {
        if (h.family == Family::Sam3) {
            // The ViT's RoPE tables ship in the checkpoint, sized for the
            // native token grid. Running another resolution needs them
            // regenerated, which the exporter does not do -- fail loudly rather
            // than silently attending with the wrong frequencies.
            fail("this checkpoint's rotary tables are exported for a %dx%d input "
                 "(%dx%d tokens); --img-size %d would need them regenerated",
                 h.img_size, h.img_size, h.grid(), h.grid(), img_size_);
        }
        // Hiera has no baked tables: the position embedding is interpolated
        // below at whatever grid we ask for. It does need all four stage
        // resolutions to stay even.
        NN_CHECK(img_size_ % 32 == 0,
                   "--img-size %d is not a multiple of 32; Hiera halves the resolution "
                   "at each of its three stage transitions", img_size_);
    }

    buildTables();
    buildPromptEncoderCache();
}

void SamModel::buildTables() {
    const Hparams& h = hp();
    const int D = h.neck_dim;
    const int G = grid_;

    {
        auto pe = sinusoidal_pe_2d(G, G, D);
        neck_pe_ = upload_table(vk::PoolSlot::NeckPe, 0, pe, G, G, D);
    }
    {
        auto pe = sinusoidal_pe_2d(G, G, h.mem_out_dim);
        mem_pe_ = upload_table(vk::PoolSlot::PosEncTable, 0, pe, G, G, h.mem_out_dim);
    }
    {
        // Memory attention runs one 256-dim head, so the axial table covers the
        // full model dim rather than a per-head slice.
        auto rope = axial_rope_table(D, G, G);
        mem_rope_ = upload_table(vk::PoolSlot::RopeTable, 0, rope, (int64_t)G * G, D / 2, 2);
    }
    {
        std::vector<float> dim_t(64);
        for (int i = 0; i < 64; ++i)
            dim_t[(size_t)i] =
                2.0f * (float)M_PI / std::pow(10000.0f, 2.0f * (float)i / 128.0f);
        sine_dim_t_ = upload_table(vk::PoolSlot::PosEncTable, 1, dim_t, 64);
    }
    rpb_coords_.resize((size_t)G);
    for (int i = 0; i < G; ++i) rpb_coords_[(size_t)i] = (float)i / (float)G;

    if (h.family == Family::Sam2) buildHieraTables();
}

// Hiera._get_pos_embed, done once instead of per frame: stretch the 7x7
// background embedding to the patch-embed grid with bicubic interpolation, then
// add the window embedding tiled over it. The result is what every frame adds
// to the patch-embed output, so it is stored channel-last, ready for nn::add.
void SamModel::buildHieraTables() {
    const Hparams& h = hp();
    hiera_blocks_ = hiera_block_table(h);

    const int C = h.hiera_embed_dim;
    const int S = img_size_ / 4;   // patch_embed is kernel 7, stride 4, padding 3

    const nn::Tensor pe = store_.get("hiera.pos_embed");            // [1, C, ph, pw]
    const nn::Tensor win = store_.get("hiera.pos_embed_window");    // [1, C, ws, ws]
    const int ph = (int)pe.dim(pe.ndim - 2), pw = (int)pe.dim(pe.ndim - 1);
    const int wh = (int)win.dim(win.ndim - 2), ww = (int)win.dim(win.ndim - 1);
    NN_CHECK(pe.numel() == (int64_t)C * ph * pw && win.numel() == (int64_t)C * wh * ww,
               "hiera.pos_embed%s has an unexpected shape for embed dim %d",
               pe.numel() == (int64_t)C * ph * pw ? "_window" : "", C);
    NN_CHECK(S % wh == 0 && S % ww == 0,
               "the %dx%d Hiera window embedding does not tile a %dx%d grid", wh, ww, S, S);

    std::vector<float> pe_host((size_t)pe.numel()), win_host((size_t)win.numel());
    nn::tensor_to_host(pe, pe_host.data(), pe.numel());
    nn::tensor_to_host(win, win_host.data(), win.numel());

    std::vector<float> out((size_t)S * S * C);
    for (int c = 0; c < C; ++c) {
        const std::vector<float> plane =
            bicubic_resize(pe_host.data() + (size_t)c * ph * pw, ph, pw, S, S);
        const float* wplane = win_host.data() + (size_t)c * wh * ww;
        for (int y = 0; y < S; ++y)
            for (int x = 0; x < S; ++x)
                out[((size_t)y * S + x) * C + c] =
                    plane[(size_t)y * S + x] + wplane[(size_t)(y % wh) * ww + (x % ww)];
    }
    hiera_pos_embed_ = upload_table(vk::PoolSlot::PosEncTable, 2, out, S, S, C);
}

void SamModel::buildPromptEncoderCache() {
    const Hparams& h = hp();
    const int D = h.sam_embed_dim;
    const int G = grid_;
    const int num_pos_feats = D / 2;

    // ---- transposed-conv weight repacking ----
    // Pool subs 1.. are reserved for these (sub 0 is the weight blob itself).
    uint32_t sub = 1;
    auto repack = [&](const char* name, int64_t cin, int64_t cout) {
        return repack_conv_transpose2x2(store_.get(name), cin, cout, sub++);
    };
    // The SimpleFPN necks belong to SAM 3's ViT; Hiera's FpnNeck is 1x1 convs
    // and a nearest-neighbour top-down path, with nothing to repack.
    if (h.family == Family::Sam3) {
        const int E = h.vit_embed_dim;
        for (int path = 0; path < 2; ++path) {
            const char* prefix = (path == 0) ? "neck.det." : "neck.trk.";
            if (path == 0 && store_.visualOnly()) continue;
            deconv_[path][0][0] =
                repack((std::string(prefix) + "0.dconv_2x2_0.weight").c_str(), E, 512);
            deconv_[path][0][1] =
                repack((std::string(prefix) + "0.dconv_2x2_1.weight").c_str(), 512,
                       h.neck_dim);
            deconv_[path][1][0] =
                repack((std::string(prefix) + "1.dconv_2x2.weight").c_str(), E, 512);
        }
    }
    sam_upscale_[0] = repack("sam_dec.upscale.0.weight", h.neck_dim, 64);
    sam_upscale_[1] = repack("sam_dec.upscale.3.weight", 64, 32);
    if (!store_.visualOnly())
        geom_box_pool_w_ = repack_conv_to_hwc(store_.get("geom.boxes_pool_project.weight"),
                                              h.neck_dim, h.neck_dim, 7, sub++);

    // ---- SAM prompt encoder: small learned vectors, kept on the host ----
    pe_gaussian_ = download(store_.get("sam_pe.pe_gaussian"));
    for (int i = 0; i < 4; ++i)
        point_embed_[i] =
            download(store_.getf("sam_pe.point_embeddings.%d.weight", i));
    not_a_point_ = download(store_.get("sam_pe.not_a_point_embed.weight"));
    const std::vector<float> no_mask = download(store_.get("sam_pe.no_mask_embed.weight"));

    // ---- dense grids: the positional encoding the mask decoder attends with,
    // and the "no mask given" dense embedding, both broadcast over the grid ----
    {
        std::vector<float> pe((size_t)G * G * D);
        std::vector<float> nomask((size_t)G * G * D);
        for (int row = 0; row < G; ++row)
            for (int col = 0; col < G; ++col) {
                // Pixel centers, matching PositionEmbeddingRandom.forward.
                const float xn = ((float)col + 0.5f) / (float)G;
                const float yn = ((float)row + 0.5f) / (float)G;
                pe_encode_coord(pe.data() + ((size_t)row * G + col) * D, xn, yn,
                                pe_gaussian_.data(), num_pos_feats);
                std::memcpy(nomask.data() + ((size_t)row * G + col) * D, no_mask.data(),
                            (size_t)D * sizeof(float));
            }
        dense_pe_ = upload_table(vk::PoolSlot::PromptEncCache, 0, pe, G, G, D);
        dense_nomask_ = upload_table(vk::PoolSlot::PromptEncCache, 1, nomask, G, G, D);
    }

    // ---- tracker odds and ends kept on the host ----
    for (int i = 0; i < 3; ++i) {
        obj_ptr_w_[i] = download(store_.getf("obj_ptr_proj.layers.%d.weight", i));
        obj_ptr_b_[i] = download(store_.getf("obj_ptr_proj.layers.%d.bias", i));
    }
    no_obj_ptr_ = download(store_.get("no_obj_ptr"));
    // SAM 2.0 sets add_tpos_enc_to_obj_ptrs = false and therefore has no
    // projection; its object pointers carry no temporal encoding.
    obj_ptr_tpos_w_ = download(store_.tryGet("obj_ptr_tpos_proj.weight"));
    obj_ptr_tpos_b_ = download(store_.tryGet("obj_ptr_tpos_proj.bias"));
    maskmem_tpos_ = download(store_.get("mem_enc.tpos_enc"));
    no_mem_embed_ = download(store_.get("no_mem_embed"));
    no_obj_embed_spatial_ = download(store_.tryGet("no_obj_embed_spatial"));

    NN_LOG_INFO("[model] derived tables built (grid %dx%d, mask %dx%d)\n", G, G, G * 4,
                  G * 4);
}

}  // namespace sam
