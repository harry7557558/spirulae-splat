#pragma once
// Hyperparameters, read from the checkpoint header.
//
// Two file formats are read, both sam3.cpp's: a magic + version + ftype +
// tensor count, then a fixed sequence of int32 fields, then the tensor records,
// then (SAM 3 only) an embedded BPE vocabulary. Keeping the formats
// byte-compatible means `convert_sam3_to_ggml.py` and `convert_sam2_to_ggml.py`
// load here unchanged -- there is no second conversion step to keep in sync,
// and no reason to invent one.
//
// The two families share everything from the neck down: SAM 3's tracker head is
// SAM 2's, and the converters give both the same tensor names (`sam_dec.*`,
// `sam_pe.*`, `mem_attn.*`, `mem_enc.*`). They differ in the image backbone
// (a 32-block ViT at stride 14 vs. a 4-stage Hiera at stride 16), in whether
// there is a text encoder and detector at all, and in a handful of tracking
// behaviours flagged below.

#include <cstdint>
#include <vector>

namespace sam {

// "sam3" / "sam2" / "tok\0"
constexpr uint32_t kSam3Magic = 0x73616D33u;
constexpr uint32_t kSam2Magic = 0x73616D32u;
constexpr uint32_t kTokMagic = 0x746F6B00u;
constexpr int32_t  kSam3FileVersion = 3;
constexpr int32_t  kSam2FileVersion = 1;

enum class Family : int32_t {
    Sam3 = 0,   // Perception Encoder ViT + text tower + DETR detector
    Sam2 = 1,   // Hiera + FPN, visual prompts only
};

// ggml type ids, as written into each tensor record.
enum class FileDType : int32_t {
    F32 = 0,
    F16 = 1,
    Q4_0 = 2,
    Q4_1 = 3,
    Q8_0 = 8,
};

struct Hparams {
    Family  family = Family::Sam3;

    // ---- ViT backbone (Perception Encoder; SAM 3 only) ----
    int32_t img_size = 1008;
    // Stride of the backbone's primary feature level, which for the ViT is
    // literally the patch size and for Hiera is the stride-16 stage. `grid()`
    // is img_size / this either way.
    int32_t patch_size = 14;
    int32_t vit_embed_dim = 1024;
    int32_t vit_depth = 32;
    int32_t vit_num_heads = 16;
    int32_t vit_mlp_dim = 4736;      // embed_dim * 4.625
    int32_t vit_window_size = 24;
    int32_t n_global_attn = 4;
    int32_t global_attn_idx[8] = {7, 15, 23, 31, 0, 0, 0, 0};

    // ---- text encoder ----
    int32_t text_width = 1024;
    int32_t text_heads = 16;
    int32_t text_layers = 24;
    int32_t text_ctx_len = 32;
    int32_t text_vocab_size = 49408;
    int32_t text_out_dim = 256;

    // ---- shared model dim ----
    int32_t neck_dim = 256;

    // ---- detector ----
    int32_t fenc_layers = 6;
    int32_t fenc_heads = 8;
    int32_t fenc_ffn_dim = 2048;
    int32_t ddec_layers = 6;
    int32_t ddec_heads = 8;
    int32_t ddec_ffn_dim = 2048;
    int32_t ddec_num_queries = 200;
    int32_t geom_layers = 3;
    int32_t n_presence_tokens = 1;
    int32_t n_geom_queries = 4;

    // ---- SAM (tracker) head ----
    int32_t sam_embed_dim = 256;
    int32_t sam_dec_depth = 2;
    int32_t sam_n_multimask = 3;
    int32_t sam_iou_head_depth = 3;

    // ---- memory / video ----
    int32_t mem_out_dim = 64;
    int32_t mem_attn_layers = 4;
    int32_t num_maskmem = 7;
    int32_t max_obj_ptrs = 16;

    int32_t n_amb_experts = 2;
    int32_t visual_only = 0;   // 1 = no text encoder / detector path

    // ---- Hiera backbone + FPN neck (SAM 2 only) ----
    int32_t hiera_embed_dim = 0;
    int32_t hiera_num_heads = 0;
    int32_t hiera_n_stages = 0;
    int32_t hiera_stages[4] = {0, 0, 0, 0};   // blocks per stage
    int32_t hiera_q_pool = 0;                 // stage transitions that pool Q
    int32_t hiera_window[4] = {0, 0, 0, 0};   // window_spec, per stage
    // window_pos_embed_bkg_spatial_size: the shape of `hiera.pos_embed`, which
    // is bicubically stretched to the token grid at load.
    int32_t hiera_pe_bkg_h = 0, hiera_pe_bkg_w = 0;
    int32_t scalp = 0;                        // FPN levels dropped from the top
    int32_t fpn_top_down_n = 0;
    int32_t fpn_top_down[4] = {0, 0, 0, 0};   // levels that take a top-down add

    // ---- tracking behaviour that differs between the families ----
    // The memory encoder maps the predicted mask through
    // sigmoid(x) * scale + bias before downsampling.
    float   mem_sigmoid_scale = 20.0f;
    float   mem_sigmoid_bias = -10.0f;
    // SAM 2 decodes three ambiguity masks and keeps the one with the best
    // predicted IoU whenever the prompt is ambiguous enough (SAM2Base._use_multimask:
    // an unprompted propagation step counts as zero points, so tracking uses
    // it too). SAM 3 always takes the single-mask slot.
    int32_t multimask_in_sam = 0;
    int32_t multimask_for_tracking = 0;
    int32_t multimask_min_pts = 0;
    int32_t multimask_max_pts = 0;
    // SAM 2.0 has no obj_ptr_tpos_proj: object pointers carry no temporal
    // encoding at all. SAM 2.1 and SAM 3 do.
    int32_t is_sam2_1 = 0;

    // ---- derived ----
    int32_t grid() const { return img_size / patch_size; }             // 72
    int32_t n_img_tokens() const { return grid() * grid(); }           // 5184
    int32_t vit_head_dim() const { return vit_embed_dim / vit_num_heads; }  // 64
    int32_t mask_size() const { return grid() * 4; }                   // 288
    int32_t num_mask_tokens() const { return sam_n_multimask + 1; }    // 4

    bool is_global_attn(int layer) const {
        for (int i = 0; i < n_global_attn && i < 8; ++i)
            if (global_attn_idx[i] == layer) return true;
        return false;
    }

    // Number of tokens the ViT's attention sees in a given block: the whole
    // grid for a global block, one window otherwise.
    int32_t block_tokens(int layer) const {
        return is_global_attn(layer) ? n_img_tokens()
                                     : vit_window_size * vit_window_size;
    }

    // SAM2Base._use_multimask. `n_points` counts every prompt token the SAM
    // decoder sees (a box contributes its two corners); propagation passes 0.
    bool use_multimask(int n_points, bool init_frame) const {
        if (!multimask_in_sam) return false;
        if (!init_frame && !multimask_for_tracking) return false;
        return n_points >= multimask_min_pts && n_points <= multimask_max_pts;
    }

    int32_t hiera_depth() const {
        int32_t n = 0;
        for (int i = 0; i < hiera_n_stages && i < 4; ++i) n += hiera_stages[i];
        return n;
    }
};

// ---------------------------------------------------------------------------
// Hiera block table
// ---------------------------------------------------------------------------

// What one MultiScaleBlock does. Derived entirely from the hyperparameters;
// see Hiera.__init__ in hieradet.py, whose ordering has two traps worth
// spelling out:
//
//   * `window` is read from the PREVIOUS stage's window_spec on the first
//     block of each new stage, because the spec lookup happens before
//     cur_stage is bumped. That is not a typo in the reference -- the block
//     partitions at the pre-pool resolution, so it needs the pre-pool window.
//   * dim and head count double on the first block of a new stage, and that is
//     also exactly where Q pooling happens, so `proj` and `q_pool` always come
//     as a pair.
struct HieraBlock {
    int dim_in = 0;
    int dim_out = 0;      // != dim_in on the first block of a stage (has `proj`)
    int n_heads = 0;
    int window = 0;       // 0 = global attention over the whole map
    bool q_pool = false;  // 2x2 max pool on Q (and on the residual)
    bool stage_end = false;
};

inline std::vector<HieraBlock> hiera_block_table(const Hparams& h) {
    std::vector<HieraBlock> out;
    const int depth = h.hiera_depth();
    out.reserve((size_t)depth);

    int stage_end[4] = {0, 0, 0, 0};
    int acc = 0;
    for (int s = 0; s < h.hiera_n_stages && s < 4; ++s) {
        acc += h.hiera_stages[s];
        stage_end[s] = acc - 1;
    }
    auto is_stage_end = [&](int i) {
        for (int s = 0; s < h.hiera_n_stages && s < 4; ++s)
            if (stage_end[s] == i) return true;
        return false;
    };

    int dim = h.hiera_embed_dim;
    int heads = h.hiera_num_heads;
    int cur_stage = 1;
    int pooled = 0;
    for (int i = 0; i < depth; ++i) {
        HieraBlock b;
        b.dim_in = dim;
        b.dim_out = dim;
        b.window = h.hiera_window[cur_stage - 1];
        if (h.is_global_attn(i)) b.window = 0;
        if (i > 0 && is_stage_end(i - 1)) {
            b.dim_out = dim * 2;
            heads *= 2;
            ++cur_stage;
            if (pooled < h.hiera_q_pool) {
                b.q_pool = true;
                ++pooled;
            }
        }
        b.n_heads = heads;
        b.stage_end = is_stage_end(i);
        out.push_back(b);
        dim = b.dim_out;
    }
    return out;
}

}  // namespace sam
