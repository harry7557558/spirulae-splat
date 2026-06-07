#pragma once

// EngineConfig -- per-step config structs accepted by the engine_*_step
// entrypoints. These collapse what used to be ~40-arg signatures into a
// handful of typed bundles, while keeping all scheduling on the Python side
// (every numeric field below is the already-resolved value for this step).

#include "PerPixelLoss.cuh"   // LossWeightIndex
#include "PixelWise.cuh"      // PPISPRegLossIndex

#include <array>


// engine_compute_loss_backward keeps its (step, weights, w_ssim, scales, map)
// signature; LossConfig exists to bundle those four scalars when called
// transitively from engine_train_step.
struct LossConfig {
    std::array<float, (int)LossWeightIndex::length> weights{};
    float w_ssim          = 0.0f;
    int   num_loss_scales = 1;
    bool  compute_loss_map = false;
    // When true, the densification loss map is filled with the SSIM structure
    // term only (s(x,y) = (2*sigma12 + C2) / (2*sigma1*sigma2 + C2)) and the
    // per-pixel L1/L2 + auxiliary supervisory terms are NOT accumulated into
    // the map. Affects loss_map output only; gradients, scalar loss values
    // and the SSIM display scalar are unchanged.
    bool  structure_only_loss_map = false;
    // Image-space overexposure regularization weight. When non-zero, a
    // dedicated kernel adds dL/dx of L = w * mean(max(-x, x-1, 0)^2) directly
    // into v_render_rgb (in the pre-bilagrid / pre-PPISP / pre-color-space
    // working space, which is what raster bwd consumes). The scalar loss is
    // never computed.
    float overexposure_reg_weight = 0.0f;
};


// Per-param learning rates + per-splat regularization weights + flags.
// All LR values are *already* multiplied by train_frame_scale / alpha on the
// Python side where needed; the engine consumes them verbatim.
struct OptimConfig {
    float lr_means        = 0.0f;
    float lr_quats        = 0.0f;
    float lr_scales       = 0.0f;
    float lr_opacities    = 0.0f;
    float lr_features_dc  = 0.0f;
    float lr_features_sh  = 0.0f;

    float max_gauss_ratio              = 0.0f;
    float scale_regularization_weight  = 0.0f;
    float mcmc_opacity_reg_weight      = 0.0f;
    float mcmc_scale_reg_weight        = 0.0f;
    float erank_reg_weight             = 0.0f;
    float erank_reg_weight_s3          = 0.0f;
    float quat_norm_reg_weight         = 0.0f;
    float sh_reg_weight                = 0.0f;

    bool  use_scale_agnostic_mean         = false;
    // SH-Adam optimizer-state quantization bit depth. 32 = no quantization (full
    // fp32 g1/g2). 4 or 8 = QuantizedAdamState<BITS, 256> with one float4 per
    // 256-cell block of (u, sqrt_g2) bounds. Other values are rejected.
    int   sh_optim_bits                   = 32;
    bool  use_per_splat_bias_correction   = false;

    // When true, fold projection-backward and Adam-based per-splat optim into
    // a single fused kernel (FusedProjectionBwdOptim). The engine then skips
    // the standalone projection_*_backward + engine_optim_step calls. SH
    // quantization is not yet supported by the fused path; the MVP allocates
    // full fp32 g1/g2 momentum (no `sh_quant_state`) in this mode.
    bool  use_fused_proj_bwd_optim        = false;

    // Trust-region color-space Adam.  Mirrors the
    // fused_adamtr_(linear_)rgb_(sh_)optim variants used in Python when
    // splat_color_is_linear or splat_color_gamut is set: the DC and SH color
    // updates are clipped to +/-kSh0*sqrt(4*eps_tr*c/opac) per step so the
    // working-color-space update stays inside the model's confidence radius.
    bool  use_color_trust_region          = false;
    bool  color_is_linear                 = false;   // gradient gets divided by linear->sRGB Jacobian
    float eps_tr                          = 1e-6f;
};


// MCMC / revised-relocate densification controls. max_world_size / noise_lr*
// are pre-scaled by alpha on the Python side.
struct DensifyConfig {
    int   refine_start_iter             = 0;
    int   refine_stop_num_iter          = 0;
    int   refine_every                  = 0;
    float growth_factor                 = 1.0f;
    float min_opacity                   = 0.0f;
    float max_screen_size               = 0.0f;
    float max_screen_size_clip_hardness = 0.0f;
    float max_world_size                = 0.0f;
    float noise_lr                      = 0.0f;
    float noise_lr_final                = 0.0f;
    bool use_revised_densification      = true;
    int  score_mode                     = 0;    // 0=mean, 1=max, 2=median, 3=geom
};


// Per-type Adam LR + TV regularization weight. lr <= 0 disables the channel
// for the current step (so a single config covers "enabled but skipped" too).
//
// When the engine was init'd with use_adagrad=true for a given type, those
// (lr_*) values are the AdaGrad LRs (typically larger than Adam's, e.g. 1e-1
// for RGB, since AdaGrad's effective per-parameter LR shrinks as the
// accumulator grows).
struct BilagridStepConfig {
    float lr_rgb         = 0.0f;
    float lr_depth       = 0.0f;
    float lr_normal      = 0.0f;
    float tv_weight_rgb  = 0.0f;
    float tv_weight_depth = 0.0f;
    float tv_weight_normal = 0.0f;
    // RGB color-shift regularizer (design 1, BilagridShiftReg.cu). Penalizes
    // the dataset-wide mean sign(post - pre); zero disables the kernel launch
    // entirely. beta is the EMA decay (e.g. 1 - 1/T where T = #batches per
    // epoch); ignored when weight is 0. Applied during loss backward, BEFORE
    // bilagrid bwd consumes v_render_rgb.
    float shift_reg_weight_rgb = 0.0f;
    float shift_reg_beta_rgb   = 0.0f;
};


// Background blending step config. When mode=Noise, only `randomize_weight`
// Adam over the per-band SH coefficients (index 0 = DC color, 1+ = bands).
struct BackgroundStepConfig {
    float    lr_dc           = 0.0f;
    float    lr_sh           = 0.0f;
    float    randomize_weight = 0.0f;
    uint32_t seed             = 0;
};


// PPISP Adam LR + 6-component regularization weights (indexed by PPISPRegLossIndex).
struct PpispStepConfig {
    float lr = 0.0f;
    std::array<float, (int)PPISPRegLossIndex::length> reg_weights{};
};


// Bundle passed to engine_train_step covering all per-step config groups.
struct EngineStepConfig {
    LossConfig           loss;
    OptimConfig          optim;
    DensifyConfig        densify;
    BilagridStepConfig   bilagrid;
    PpispStepConfig      ppisp;
    BackgroundStepConfig background;
};
