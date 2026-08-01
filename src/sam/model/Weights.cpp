#include "sam/model/Weights.h"

#include "sam/Common.h"
#include "nn/core/Error.h"
#include "nn/core/Half.h"
#include "nn/core/Log.h"
#include "nn/vk/Stream.h"

#include <algorithm>
#include <cstdarg>
#include <cstring>
#include <fstream>

namespace sam {

namespace {

// Tensors smaller than this stay fp32 on device: learned tokens, reference
// points, label/type embeddings, temporal PE tables. Together they are a few
// megabytes, and several of them feed coordinate arithmetic where fp16's
// 11-bit mantissa is visible in the output boxes.
constexpr int64_t kSmallTensorElems = 65536;
// Chunk size for the weight blob. Kept well under the ~1 GiB single-allocation
// ceiling NVIDIA's Windows driver enforces regardless of free VRAM, and large
// enough that a 1.6 GiB model is seven allocations rather than 1465.
constexpr uint64_t kWeightChunkBytes = 256ull << 20;

bool force_f32(const std::string& name) {
    // The RoPE tables are read as raw f32 by misc.slang::rope_apply.
    return name.find("freqs_cis") != std::string::npos;
}

template <typename T>
bool read_pod(std::istream& in, T& v) {
    in.read(reinterpret_cast<char*>(&v), sizeof(T));
    return (bool)in;
}

const char* file_dtype_name(FileDType t) {
    switch (t) {
        case FileDType::F32:  return "f32";
        case FileDType::F16:  return "f16";
        case FileDType::Q4_0: return "q4_0";
        case FileDType::Q4_1: return "q4_1";
        case FileDType::Q8_0: return "q8_0";
    }
    return "?";
}

// Bytes a tensor of `n` elements occupies in the file.
uint64_t file_bytes_for(FileDType t, int64_t n, int64_t row_len) {
    switch (t) {
        case FileDType::F32: return (uint64_t)n * 4;
        case FileDType::F16: return (uint64_t)n * 2;
        // ggml block quantizations: 32 elements per block, along the row.
        case FileDType::Q4_0: return (uint64_t)(n / 32) * 18;
        case FileDType::Q4_1: return (uint64_t)(n / 32) * 20;
        case FileDType::Q8_0: return (uint64_t)(n / 32) * 34;
    }
    (void)row_len;
    return 0;
}

// ---------------------------------------------------------------------------
// ggml block dequantization
//
// These are the exact layouts ggml writes; the nibble ordering of q4 in
// particular is not the obvious one (low nibbles hold elements 0..15 and high
// nibbles elements 16..31 of the same block), and getting it wrong produces a
// model that loads, runs, and returns noise.
// ---------------------------------------------------------------------------

void dequant_q4_0(const uint8_t* src, float* dst, int64_t n) {
    const int64_t nb = n / 32;
    for (int64_t b = 0; b < nb; ++b) {
        const uint8_t* p = src + b * 18;
        uint16_t dh;
        std::memcpy(&dh, p, 2);
        const float d = half_to_float(dh);
        const uint8_t* qs = p + 2;
        float* out = dst + b * 32;
        for (int j = 0; j < 16; ++j) {
            out[j] = ((int)(qs[j] & 0x0F) - 8) * d;
            out[j + 16] = ((int)(qs[j] >> 4) - 8) * d;
        }
    }
}

void dequant_q4_1(const uint8_t* src, float* dst, int64_t n) {
    const int64_t nb = n / 32;
    for (int64_t b = 0; b < nb; ++b) {
        const uint8_t* p = src + b * 20;
        uint16_t dh, mh;
        std::memcpy(&dh, p, 2);
        std::memcpy(&mh, p + 2, 2);
        const float d = half_to_float(dh), m = half_to_float(mh);
        const uint8_t* qs = p + 4;
        float* out = dst + b * 32;
        for (int j = 0; j < 16; ++j) {
            out[j] = (float)(qs[j] & 0x0F) * d + m;
            out[j + 16] = (float)(qs[j] >> 4) * d + m;
        }
    }
}

void dequant_q8_0(const uint8_t* src, float* dst, int64_t n) {
    const int64_t nb = n / 32;
    for (int64_t b = 0; b < nb; ++b) {
        const uint8_t* p = src + b * 34;
        uint16_t dh;
        std::memcpy(&dh, p, 2);
        const float d = half_to_float(dh);
        const int8_t* qs = reinterpret_cast<const int8_t*>(p + 2);
        float* out = dst + b * 32;
        for (int j = 0; j < 32; ++j) out[j] = (float)qs[j] * d;
    }
}

}  // namespace

// ================
// Loading
// ================

void WeightStore::load(const std::string& path) {
    std::ifstream fin(path, std::ios::binary);
    NN_CHECK((bool)fin, "cannot open model file '%s'", path.c_str());

    uint32_t magic = 0;
    int32_t version = 0, ftype = 0, n_tensors = 0;
    NN_CHECK(read_pod(fin, magic) && read_pod(fin, version) && read_pod(fin, ftype) &&
                   read_pod(fin, n_tensors),
               "'%s' is too short to be a model file", path.c_str());
    NN_CHECK(magic == kSam3Magic || magic == kSam2Magic,
               "'%s' has magic 0x%08x; expected 0x%08x (\"sam3\") or 0x%08x (\"sam2\")",
               path.c_str(), magic, kSam3Magic, kSam2Magic);
    const bool sam2 = magic == kSam2Magic;
    const int32_t want_version = sam2 ? kSam2FileVersion : kSam3FileVersion;
    NN_CHECK(version == want_version,
               "'%s' is format version %d, this build reads version %d for %s files",
               path.c_str(), version, want_version, sam2 ? "sam2" : "sam3");
    NN_CHECK(n_tensors > 0 && n_tensors < 100000, "'%s' declares %d tensors",
               path.c_str(), n_tensors);

    // ---- hyperparameters (fixed field order, see Hparams.h) ----
    {
        Hparams& h = hparams_;
        auto rd = [&](int32_t& v) { read_pod(fin, v); };
        auto skip = [&](int n) {
            int32_t junk;
            for (int i = 0; i < n; ++i) rd(junk);
        };
        if (sam2) {
            h.family = Family::Sam2;
            int32_t backbone_type = 0;
            rd(h.img_size);
            rd(backbone_type);
            NN_CHECK(backbone_type == 1,
                       "'%s' declares backbone type %d; only 1 (Hiera) is implemented",
                       path.c_str(), backbone_type);
            rd(h.hiera_embed_dim); rd(h.hiera_num_heads); rd(h.hiera_n_stages);
            for (int i = 0; i < 4; ++i) rd(h.hiera_stages[i]);
            rd(h.n_global_attn);
            for (int i = 0; i < 8; ++i) rd(h.global_attn_idx[i]);
            rd(h.hiera_q_pool);
            for (int i = 0; i < 4; ++i) rd(h.hiera_window[i]);
            rd(h.hiera_pe_bkg_h); rd(h.hiera_pe_bkg_w);
            rd(h.scalp);
            rd(h.neck_dim);
            rd(h.fpn_top_down_n);
            for (int i = 0; i < 4; ++i) rd(h.fpn_top_down[i]);
            rd(h.sam_embed_dim); rd(h.sam_dec_depth); rd(h.sam_n_multimask);
            rd(h.sam_iou_head_depth);
            rd(h.mem_out_dim); rd(h.mem_attn_layers); rd(h.num_maskmem);
            rd(h.max_obj_ptrs);
            int32_t sig_scale = 0, sig_bias = 0;
            rd(sig_scale); rd(sig_bias);
            h.mem_sigmoid_scale = (float)sig_scale / 100.0f;
            h.mem_sigmoid_bias = (float)sig_bias / 100.0f;
            // Flags, in the converter's order. Only the ones this port acts on
            // are kept; the rest describe defaults it already follows.
            skip(1);                             // use_high_res_features
            skip(1);                             // use_obj_ptrs_in_encoder
            skip(1);                             // pred_obj_scores
            skip(1);                             // use_multimask_token_for_obj_ptr
            skip(1);                             // directly_add_no_mem_embed
            skip(1);                             // non_overlap_masks_for_mem_enc
            skip(1);                             // binarize_mask_from_pts
            rd(h.multimask_for_tracking);
            rd(h.multimask_min_pts); rd(h.multimask_max_pts);
            skip(1);                             // fixed_no_obj_ptr
            skip(1);                             // iou_prediction_use_sigmoid
            skip(1);                             // use_mask_input_as_output
            rd(h.multimask_in_sam);
            rd(h.is_sam2_1);

            // Hiera's stride-16 stage is the grid every head below the neck
            // works on, so it plays the ViT's patch size.
            h.patch_size = 16;
            h.visual_only = 1;                   // no text tower, no detector
            NN_CHECK(h.hiera_n_stages == 4,
                       "'%s' declares %d Hiera stages; the FPN neck assumes 4",
                       path.c_str(), h.hiera_n_stages);
        } else {
            rd(h.img_size); rd(h.patch_size); rd(h.vit_embed_dim); rd(h.vit_depth);
            rd(h.vit_num_heads);
            int32_t mlp_ratio_x1000 = 0;
            rd(mlp_ratio_x1000);
            h.vit_mlp_dim = (int32_t)((float)h.vit_embed_dim * (mlp_ratio_x1000 / 1000.0f));
            rd(h.vit_window_size);
            rd(h.n_global_attn);
            NN_CHECK(h.n_global_attn >= 0 && h.n_global_attn <= 8,
                       "model declares %d global-attention blocks", h.n_global_attn);
            for (int i = 0; i < h.n_global_attn; ++i) rd(h.global_attn_idx[i]);
            rd(h.text_width); rd(h.text_heads); rd(h.text_layers); rd(h.text_ctx_len);
            rd(h.text_vocab_size); rd(h.text_out_dim);
            rd(h.neck_dim);
            rd(h.fenc_layers); rd(h.fenc_heads); rd(h.fenc_ffn_dim);
            rd(h.ddec_layers); rd(h.ddec_heads); rd(h.ddec_ffn_dim);
            rd(h.ddec_num_queries);
            rd(h.geom_layers); rd(h.n_presence_tokens); rd(h.n_geom_queries);
            rd(h.sam_embed_dim); rd(h.sam_dec_depth); rd(h.sam_n_multimask);
            rd(h.sam_iou_head_depth);
            rd(h.mem_out_dim); rd(h.mem_attn_layers); rd(h.num_maskmem);
            rd(h.max_obj_ptrs);
            rd(h.n_amb_experts); rd(h.visual_only);
        }
        NN_CHECK(h.n_global_attn >= 0 && h.n_global_attn <= 8,
                   "model declares %d global-attention blocks", h.n_global_attn);
        NN_CHECK((bool)fin, "'%s' ended inside the hyperparameter block", path.c_str());
    }

    if (sam2)
        NN_LOG_INFO("[model] %s: SAM %s, %d tensors, ftype %s, %dx%d input, "
                      "hiera dim %d depth %d\n",
                      path.c_str(), hparams_.is_sam2_1 ? "2.1" : "2.0", n_tensors,
                      file_dtype_name((FileDType)ftype), hparams_.img_size,
                      hparams_.img_size, hparams_.hiera_embed_dim, hparams_.hiera_depth());
    else
        NN_LOG_INFO("[model] %s: SAM 3, %d tensors, ftype %s, %dx%d input, %s\n",
                      path.c_str(), n_tensors, file_dtype_name((FileDType)ftype),
                      hparams_.img_size, hparams_.img_size,
                      hparams_.visual_only ? "visual-only" : "full (text + detector)");

    // ---- pass 1: read the records, plan the device layout ----
    std::vector<Record> records;
    records.reserve((size_t)n_tensors);
    uint64_t device_cursor = 0;
    uint64_t chunk_cursor = 0;
    std::vector<uint64_t> chunk_sizes;

    for (int32_t t = 0; t < n_tensors; ++t) {
        int32_t n_dims = 0, name_len = 0, dtype = 0;
        NN_CHECK(read_pod(fin, n_dims) && read_pod(fin, name_len) && read_pod(fin, dtype),
                   "truncated tensor record %d in '%s'", t, path.c_str());
        NN_CHECK(n_dims >= 1 && n_dims <= 4, "tensor %d has %d dimensions", t, n_dims);
        NN_CHECK(name_len > 0 && name_len < 1024, "tensor %d has a %d-byte name", t,
                   name_len);

        Record r;
        r.dtype = (FileDType)dtype;
        // The file stores ggml's `ne[]`: fastest-varying dimension first, which
        // is the reverse of PyTorch's shape (convert_sam3_to_ggml.py writes
        // `for dim in reversed(data.shape)`). Undo it here, once, so every
        // shape above this line is PyTorch order and also memory order.
        r.shape.resize((size_t)n_dims);
        for (int i = 0; i < n_dims; ++i) {
            int32_t d = 0;
            read_pod(fin, d);
            r.shape[(size_t)(n_dims - 1 - i)] = d;
        }
        r.name.resize((size_t)name_len);
        fin.read(&r.name[0], name_len);

        // Tensor payloads are 32-byte aligned in the file.
        uint64_t pos = (uint64_t)fin.tellg();
        uint64_t pad = (32 - pos % 32) % 32;
        if (pad) fin.seekg((std::streamoff)pad, std::ios::cur);

        int64_t n = 1;
        for (int64_t d : r.shape) n *= d;
        r.file_offset = (uint64_t)fin.tellg();
        r.file_bytes = file_bytes_for(r.dtype, n, r.shape.back());
        NN_CHECK(r.file_bytes > 0, "tensor '%s' has unsupported dtype %d",
                   r.name.c_str(), dtype);

        const bool small = n < kSmallTensorElems;
        r.device_dtype = (n_dims == 1 || small || force_f32(r.name)) ? nn::DType::F32
                                                                    : nn::DType::F16;
        const uint64_t bytes = ((uint64_t)n * nn::dtype_size(r.device_dtype) + 255) & ~255ull;
        if (chunk_cursor != 0 && chunk_cursor + bytes > kWeightChunkBytes) {
            chunk_sizes.push_back(chunk_cursor);
            chunk_cursor = 0;
        }
        r.chunk = (uint32_t)chunk_sizes.size();
        r.device_offset = chunk_cursor;
        chunk_cursor += bytes;
        device_cursor += bytes;

        fin.seekg((std::streamoff)r.file_bytes, std::ios::cur);
        NN_CHECK((bool)fin, "'%s' ended inside tensor '%s'", path.c_str(),
                   r.name.c_str());
        records.push_back(std::move(r));
    }

    // ---- pass 2: one allocation per chunk, then upload ----
    if (chunk_cursor) chunk_sizes.push_back(chunk_cursor);
    device_bytes_ = device_cursor;
    chunks_.clear();
    for (size_t i = 0; i < chunk_sizes.size(); ++i)
        chunks_.push_back(
            vk::VramPool::get().acquire(vk::PoolSlot::Weights, (uint32_t)i, chunk_sizes[i]));
    NN_LOG_INFO("[model] weights: %.1f MiB on device in %zu chunks (%d tensors)\n",
                  device_cursor / 1048576.0, chunk_sizes.size(), n_tensors);

    std::vector<uint8_t> raw;
    std::vector<float> f32;
    std::vector<uint16_t> f16;
    const double t0 = now_ms();

    for (const Record& r : records) {
        int64_t n = 1;
        for (int64_t d : r.shape) n *= d;

        fin.clear();
        fin.seekg((std::streamoff)r.file_offset, std::ios::beg);
        raw.resize((size_t)r.file_bytes);
        fin.read((char*)raw.data(), (std::streamsize)r.file_bytes);
        NN_CHECK((bool)fin, "failed to read tensor '%s'", r.name.c_str());

        const nn::DevicePtr dst = chunks_[r.chunk] + r.device_offset;

        // Fast path: the file already holds what the device wants.
        if ((r.dtype == FileDType::F16 && r.device_dtype == nn::DType::F16) ||
            (r.dtype == FileDType::F32 && r.device_dtype == nn::DType::F32)) {
            vk::Stream::get().upload(dst, raw.data(), r.file_bytes);
        } else {
            // Everything else goes through f32 as the common form.
            f32.resize((size_t)n);
            switch (r.dtype) {
                case FileDType::F32:
                    std::memcpy(f32.data(), raw.data(), (size_t)n * 4);
                    break;
                case FileDType::F16:
                    for (int64_t i = 0; i < n; ++i) {
                        uint16_t h;
                        std::memcpy(&h, raw.data() + i * 2, 2);
                        f32[(size_t)i] = half_to_float(h);
                    }
                    break;
                case FileDType::Q4_0: dequant_q4_0(raw.data(), f32.data(), n); break;
                case FileDType::Q4_1: dequant_q4_1(raw.data(), f32.data(), n); break;
                case FileDType::Q8_0: dequant_q8_0(raw.data(), f32.data(), n); break;
            }
            if (r.device_dtype == nn::DType::F32) {
                vk::Stream::get().upload(dst, f32.data(), (uint64_t)n * 4);
            } else {
                f16.resize((size_t)n);
                for (int64_t i = 0; i < n; ++i) f16[(size_t)i] = float_to_half(f32[(size_t)i]);
                vk::Stream::get().upload(dst, f16.data(), (uint64_t)n * 2);
            }
        }

        nn::Tensor tensor;
        tensor.ptr = dst;
        tensor.dtype = r.device_dtype;
        tensor.ndim = (int32_t)r.shape.size();
        for (size_t i = 0; i < r.shape.size(); ++i) tensor.shape[i] = r.shape[i];
        tensors_.emplace(r.name, tensor);
    }
    vk::Stream::get().sync();
    NN_LOG_INFO("[model] upload: %.0f ms\n", now_ms() - t0);

    // ---- embedded tokenizer (full models only) ----
    if (!hparams_.visual_only) {
        if (!tokenizer_.loadFromStream(fin))
            NN_LOG_WARN("[model] no usable tokenizer in '%s'; text prompts will be "
                          "rejected\n", path.c_str());
    } else {
        NN_LOG_INFO("[model] visual-only checkpoint: no tokenizer, no detector path\n");
    }
}

// ================
// Lookup
// ================

nn::Tensor WeightStore::tryGet(const std::string& name) const {
    auto it = tensors_.find(name);
    return it == tensors_.end() ? nn::Tensor{} : it->second;
}

nn::Tensor WeightStore::get(const std::string& name) const {
    auto it = tensors_.find(name);
    if (it != tensors_.end()) return it->second;

    // Name typos are otherwise a null tensor that faults somewhere unrelated.
    // Offer the closest-looking keys by shared prefix.
    std::string best;
    size_t best_len = 0;
    for (const auto& kv : tensors_) {
        size_t i = 0;
        while (i < name.size() && i < kv.first.size() && name[i] == kv.first[i]) ++i;
        if (i > best_len) { best_len = i; best = kv.first; }
    }
    fail("model has no tensor '%s'%s%s", name.c_str(),
         best_len > 4 ? "; closest name is " : "", best_len > 4 ? best.c_str() : "");
}

nn::Tensor WeightStore::getf(const char* fmt, ...) const {
    char buf[512];
    va_list ap;
    va_start(ap, fmt);
    std::vsnprintf(buf, sizeof buf, fmt, ap);
    va_end(ap);
    return get(buf);
}

nn::Tensor WeightStore::tryGetf(const char* fmt, ...) const {
    char buf[512];
    va_list ap;
    va_start(ap, fmt);
    std::vsnprintf(buf, sizeof buf, fmt, ap);
    va_end(ap);
    return tryGet(buf);
}

// ================
// Post-processing
// ================

nn::Tensor repack_conv_transpose2x2(const nn::Tensor& w, int64_t cin, int64_t cout,
                                    uint32_t pool_sub) {
    NN_CHECK(w.numel() == cin * cout * 4,
               "repack_conv_transpose2x2: weight has %lld elements, expected %lld",
               (long long)w.numel(), (long long)(cin * cout * 4));
    // Small enough to shuffle on the host (the largest is 1024x512x2x2 = 4 MiB)
    // and it happens once at model build, not per frame.
    std::vector<float> src((size_t)w.numel());
    nn::tensor_to_host(w, src.data(), w.numel());

    std::vector<float> dst((size_t)w.numel());
    for (int64_t ci = 0; ci < cin; ++ci)
        for (int64_t co = 0; co < cout; ++co)
            for (int64_t i = 0; i < 2; ++i)
                for (int64_t j = 0; j < 2; ++j)
                    dst[(size_t)((co * 4 + i * 2 + j) * cin + ci)] =
                        src[(size_t)(((ci * cout + co) * 2 + i) * 2 + j)];

    nn::Tensor out = nn::pool_tensor(vk::PoolSlot::Weights, pool_sub, w.dtype, cout * 4, cin);
    nn::tensor_from_host(out, dst.data(), (int64_t)dst.size());
    return out;
}

nn::Tensor repack_conv_to_hwc(const nn::Tensor& w, int64_t cout, int64_t cin, int64_t k,
                              uint32_t pool_sub) {
    NN_CHECK(w.numel() == cout * cin * k * k,
               "repack_conv_to_hwc: weight has %lld elements, expected %lld",
               (long long)w.numel(), (long long)(cout * cin * k * k));
    std::vector<float> src((size_t)w.numel());
    nn::tensor_to_host(w, src.data(), w.numel());
    std::vector<float> dst((size_t)w.numel());
    for (int64_t co = 0; co < cout; ++co)
        for (int64_t ci = 0; ci < cin; ++ci)
            for (int64_t ky = 0; ky < k; ++ky)
                for (int64_t kx = 0; kx < k; ++kx)
                    dst[(size_t)(co * k * k * cin + (ky * k + kx) * cin + ci)] =
                        src[(size_t)(((co * cin + ci) * k + ky) * k + kx)];
    nn::Tensor out =
        nn::pool_tensor(vk::PoolSlot::Weights, pool_sub, w.dtype, cout, k * k * cin);
    nn::tensor_from_host(out, dst.data(), (int64_t)dst.size());
    return out;
}

}  // namespace sam
