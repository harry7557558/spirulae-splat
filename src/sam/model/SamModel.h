#pragma once
// The loaded model: weights plus everything derivable from them that is worth
// computing once.
//
// The derived tables (sinusoidal position encodings, the memory-attention RoPE
// grid, the prompt encoder's dense grids, the repacked transposed-conv weights)
// are pure functions of the hyperparameters and the checkpoint. sam3.cpp
// recomputes several of them per frame and per prompt; here they are built at
// load and live in the VRAM pool.

#include "sam/model/Weights.h"
#include "nn/Ops.h"
#include "nn/vk/Memory.h"

#include <memory>
#include <string>
#include <vector>

namespace sam {

class SamModel {
public:
    // Loads the checkpoint and builds the derived tables. `img_size` overrides
    // the checkpoint's native input resolution (0 = keep it). SAM 3's ViT is
    // resolution-agnostic apart from its RoPE tables, which ship in the
    // checkpoint sized for one grid, so an override is rejected there; Hiera
    // interpolates its position embedding at load and accepts any multiple of
    // 32.
    void load(const std::string& path, int img_size_override = 0);

    const Hparams&      hp() const { return store_.hparams(); }
    const WeightStore&  w() const { return store_; }
    const BpeTokenizer& tokenizer() const { return store_.tokenizer(); }

    Family family() const { return store_.hparams().family; }
    bool visualOnly() const { return store_.visualOnly(); }
    bool canTextPrompt() const { return !visualOnly() && store_.hasTokenizer(); }

    int imgSize() const { return img_size_; }
    int grid() const { return grid_; }
    int maskSize() const { return grid_ * 4; }

    // ---- derived tables (device resident) ----

    // Sinusoidal PositionEmbeddingSine on the token grid, [G, G, neck_dim]:
    // memory attention's `curr_pos` and the fusion encoder's image PE. The
    // other FPN levels have one too in the reference, but nothing downstream
    // reads them, and at 288x288x256 each they are not worth 100 MiB of VRAM.
    const nn::Tensor& neckPe() const { return neck_pe_; }
    // Hiera's absolute position embedding, already bicubically stretched to the
    // patch-embed grid and summed with the tiled window embedding:
    // [img_size/4, img_size/4, hiera_embed_dim]. Empty for SAM 3.
    const nn::Tensor& hieraPosEmbed() const { return hiera_pos_embed_; }
    // Per-block dims / heads / window size, derived from the hyperparameters.
    const std::vector<HieraBlock>& hieraBlocks() const { return hiera_blocks_; }
    // The same encoding at mem_out_dim (64) on the token grid -- the spatial PE
    // stored with every memory-bank slot.
    const nn::Tensor& memPe() const { return mem_pe_; }
    // 2-D axial RoPE for memory attention: [grid*grid, neck_dim/2, 2].
    const nn::Tensor& memRope() const { return mem_rope_; }

    // SAM prompt encoder, precomputed once: the dense positional grid and the
    // tiled no-mask embedding, both [grid, grid, 256].
    const nn::Tensor& densePe() const { return dense_pe_; }
    const nn::Tensor& denseNoMask() const { return dense_nomask_; }
    // Host copies of the small prompt embeddings; sparse prompt tokens are
    // assembled on the host (a handful of 256-float rows) and uploaded whole.
    const float* pointEmbed(int label) const { return point_embed_[label].data(); }
    const float* notAPointEmbed() const { return not_a_point_.data(); }
    const float* peGaussian() const { return pe_gaussian_.data(); }

    // DETR: the [64] angle multipliers for gen_sineembed_for_position, and the
    // [grid] normalized coordinates the box-relative bias is built from.
    const nn::Tensor& sineDimT() const { return sine_dim_t_; }
    const std::vector<float>& rpbCoords() const { return rpb_coords_; }

    // Transposed-conv weights repacked to [Cout*4, Cin] (see nn::conv_transpose2x2).
    const nn::Tensor& neckDeconv(bool det, int scale, int which) const {
        return deconv_[(det ? 0 : 1)][scale][which];
    }
    const nn::Tensor& samUpscale(int which) const { return sam_upscale_[which]; }

    // Geometry encoder's Conv2d(D, D, 7) box pooling projection, reordered to
    // [D, 7*7*D] so it multiplies the channel-last ROI-Align patch directly.
    const nn::Tensor& geomBoxPoolWeight() const { return geom_box_pool_w_; }

    // Host copies of the object-pointer projection, used by the tracker's
    // per-instance pointer extraction (one 256-vector per instance per frame --
    // three tiny matmuls that are not worth a dispatch each).
    const std::vector<float>& objPtrProjW(int layer) const { return obj_ptr_w_[layer]; }
    const std::vector<float>& objPtrProjB(int layer) const { return obj_ptr_b_[layer]; }
    const std::vector<float>& noObjPtr() const { return no_obj_ptr_; }
    const std::vector<float>& objPtrTposW() const { return obj_ptr_tpos_w_; }
    const std::vector<float>& objPtrTposB() const { return obj_ptr_tpos_b_; }
    const std::vector<float>& maskMemTpos() const { return maskmem_tpos_; }
    const std::vector<float>& noMemEmbed() const { return no_mem_embed_; }
    const std::vector<float>& noObjEmbedSpatial() const { return no_obj_embed_spatial_; }

private:
    void buildTables();
    void buildHieraTables();
    void buildPromptEncoderCache();

    WeightStore store_;
    int img_size_ = 0;
    int grid_ = 0;

    nn::Tensor neck_pe_;
    nn::Tensor hiera_pos_embed_;
    std::vector<HieraBlock> hiera_blocks_;
    nn::Tensor mem_pe_;
    nn::Tensor mem_rope_;
    nn::Tensor dense_pe_;
    nn::Tensor dense_nomask_;
    nn::Tensor sine_dim_t_;
    std::vector<float> rpb_coords_;

    nn::Tensor deconv_[2][2][2];   // [det/trk][scale 0..1][deconv 0..1]
    nn::Tensor sam_upscale_[2];
    nn::Tensor geom_box_pool_w_;

    std::vector<float> point_embed_[4];
    std::vector<float> not_a_point_;
    std::vector<float> pe_gaussian_;
    std::vector<float> obj_ptr_w_[3], obj_ptr_b_[3];
    std::vector<float> no_obj_ptr_;
    std::vector<float> obj_ptr_tpos_w_, obj_ptr_tpos_b_;
    std::vector<float> maskmem_tpos_;
    std::vector<float> no_mem_embed_;
    std::vector<float> no_obj_embed_spatial_;
};

// ---------------------------------------------------------------------------
// Host-side encodings shared by several modules
// ---------------------------------------------------------------------------

// torch.nn.functional.interpolate(mode="bicubic", align_corners=False) on one
// [Hi, Wi] plane, with PyTorch's A = -0.75 Keys kernel and its edge clamping.
// Hiera's position embedding ships at 7x7 and is stretched to the token grid.
std::vector<float> bicubic_resize(const float* src, int Hi, int Wi, int Ho, int Wo);

// PositionEmbeddingSine(num_pos_feats = d_model/2, temperature 10000,
// normalize=True, scale=2*pi) laid out channel-last as [H, W, d_model]. The
// first half of the channels encode y, the second half x -- the order the
// reference concatenates them in.
std::vector<float> sinusoidal_pe_2d(int H, int W, int d_model);

// get_1d_sine_pe(pos, dim): first dim/2 sin, then dim/2 cos.
void sine_pe_1d(float* out, float pos, int dim, float temperature = 10000.0f);

// compute_axial_cis(dim, end_x, end_y): [N, dim/2, 2] (cos, sin) pairs, x
// frequencies in the first quarter-dim complex slots and y in the second.
std::vector<float> axial_rope_table(int dim, int end_x, int end_y,
                                    float theta = 10000.0f);

// PositionEmbeddingSine.encode_boxes -> [2*num_pos_feats + 2] = pos_y, pos_x, h, w
void sine_encode_box(float* out, float cx, float cy, float w, float h, int num_pos_feats,
                     float temperature = 10000.0f);

// The SAM prompt encoder's random-Fourier coordinate encoding:
// out[256] = [sin(2*pi * c @ G), cos(2*pi * c @ G)] with c mapped from [0,1] to
// [-1,1] and G the learned [2, 128] Gaussian matrix.
void pe_encode_coord(float* out, float x_norm, float y_norm, const float* gaussian,
                     int num_pos_feats);

}  // namespace sam
