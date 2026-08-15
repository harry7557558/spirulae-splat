#pragma once
// The shape of a Metric3D v2 checkpoint, read off the file rather than keyed
// to the model id.
//
// The three published variants differ in more than width -- giant2 swaps the
// MLP for a SwiGLU, drops the encoder's final LayerNorm and taps four
// intermediate blocks instead of reusing one output, and small runs half as
// many refinement iterations -- and every one of those differences is visible
// in the file. Deriving them means one code path, and a checkpoint that is not
// Metric3D fails while being read instead of somewhere in the forward pass.

#include <cstdint>

namespace metric3d {

struct Hparams {
    // ---- encoder (DINOv2 with registers) --------------------------------
    int embed_dim = 0;
    int depth = 0;             // transformer blocks
    int num_heads = 0;         // embed_dim / 64; DINOv2 keeps head_dim 64 at
                               // every width, and nothing in the file says so
    int patch = 14;
    int num_register_tokens = 0;
    int mlp_hidden = 0;
    bool swiglu = false;       // giant2's SwiGLUFFN in place of Mlp
    // With a final LayerNorm the export ends in `encoder.norm` and all four
    // decoder reads share that one tensor; without it they tap `taps[]`.
    bool final_norm = false;
    int taps[4] = {-1, -1, -1, -1};

    // ---- decoder (RAFTDepthNormalDPT5) ----------------------------------
    // Token2Feature output widths at 1/4, 1/7, 1/14, 1/14 of the input.
    int feat_ch[4] = {0, 0, 0, 0};
    int dec_ch[4] = {0, 0, 0, 0};   // FuseBlock trunk widths, coarse to fine
    int used_res_channel = 0;       // the shared depth/normal feature map
    int hidden = 0;                 // ConvGRU hidden and context width
    int n_bins = 0;                 // depth regressor anchors
    int iters = 0;                  // refinement iterations
    int up_factor = 0;              // convex upsample, sqrt(mask_ch / 9)

    // The depth anchors are a graph constant, so min/max come from the file;
    // `regress_scale` does not appear as a named tensor and is Metric3D's own
    // hard-coded 100.0 (RAFTDepthNormalDPTDecoder5.py, self.regress_scale).
    float min_depth = 0.0f;
    float max_depth = 0.0f;
    static constexpr float kRegressScale = 100.0f;
    // norm_normalize's floor on the concentration channel, and the epsilon it
    // ADDS to the norm rather than clamping with.
    static constexpr float kMinKappa = 0.01f;
    static constexpr float kNormEps = 1e-10f;
    // ImageNet statistics, in 0..255 units. Baked into the graph as `rgb_mean`
    // and `rgb_std`; read from the file and these are the fallback.
    float rgb_mean[3] = {123.675f, 116.28f, 103.53f};
    float rgb_std[3] = {58.395f, 57.12f, 57.375f};

    int head_dim() const { return embed_dim / num_heads; }
    // Tokens the transformer sees: cls + registers + patches.
    int64_t prefix_tokens() const { return 1 + num_register_tokens; }
};

}  // namespace metric3d
