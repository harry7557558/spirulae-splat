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
    // Selects what the densification loss map is filled with. See
    // DensifyLossMapMode in PerPixelLoss.cuh:
    //   0 = None (no map populated; densification falls back to uniform alpha*T)
    //   1 = LossFull (per-pixel L1/L2/aux + full SSIM LCS)
    //   2 = SsimFull (full SSIM LCS only)
    //   3 = SsimContrastStruct (D_/B, no luminance)
    //   4 = SsimStructure (structure-only, D_/(2*sqrt(sig1*sig2)+C2))
    //   5 = EdgeAware (canny edge magnitude of GT rgb)
    //   6 = RobustEdgeAware (Tukey biweight on |render - GT| + canny;
    //       cutoff at the per-image robust_edge_aware_quantile of |r|)
    // Affects loss_map output only; gradients, scalar loss values and the
    // SSIM display scalar are unchanged.
    int   loss_map_mode = 0;
    // Per-image quantile of the BT.601-luma residual used as the Tukey
    // cutoff in RobustEdgeAware mode. Lower values are more aggressive
    // outlier rejection (more pixels treated as distractors). Ignored
    // unless loss_map_mode == 6. Typical values: 0.8-0.95; default 0.9.
    float robust_edge_aware_quantile = 0.9f;
    // Image-space overexposure regularization weight. When non-zero, a
    // dedicated kernel adds dL/dx of L = w * mean(max(-x, x-1, 0)^2) directly
    // into v_render_rgb (in the pre-bilagrid / pre-PPISP / pre-color-space
    // working space, which is what raster bwd consumes). The scalar loss is
    // never computed.
    float overexposure_reg_weight = 0.0f;
    // Color-shift regularizer (design 1, ColorShiftReg.cu). Penalizes the
    // dataset-wide mean sign(post - pre) across the COMBINED bilagrid + PPISP
    // color transform: `pre` is the input to whichever of the two runs first
    // and `post` is the final post-both image fed to the photometric loss.
    // Zero disables the kernel launch entirely. beta is the EMA decay (e.g.
    // 1 - 1/T where T = #batches per epoch); ignored when weight is 0.
    // Applied during loss backward, BEFORE either bilagrid or PPISP bwd hooks
    // consume v_render_rgb. Active when at least one of bilagrid_rgb / PPISP
    // is enabled.
    float color_shift_reg_weight = 0.0f;
    float color_shift_reg_beta   = 0.0f;
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
    // SH PARAMETER (value) quantization bit depth. 32 = no quantization (full
    // fp32 features_sh). 8 or 16 = QuantizedTensor<BITS, 256> with one float2
    // per 256-cell block of (min, max) bounds; canonical storage is the packed
    // buffer, and every consumer kernel dequantizes on load via the codec.
    // Other values are rejected. NOTE: end-to-end kernel threading is
    // multi-turn work; the engine currently allocates storage and exposes the
    // flag but the read/write paths through Slang harmonics + FPBO + densify
    // are not yet plumbed. Setting this != 32 at training time will throw at
    // optimizer state init until that work lands.
    int   sh_value_bits                   = 8;
    // Non-SH Adam-state quantization bit depth (means, quats, scales, opacities,
    // features_dc). 32 = no quantization (full fp32 g1/g2). 16 =
    // QuantizedAdamState<16, 256> -- joint (u, log_s) at 16-bit per primitive
    // (4 B / cell). FPBO-only; the non-FPBO Adam kernel throws when this is
    // non-32. Coupled to quantization_level on the Python side: level 1
    // sets this to 16 automatically.
    int   non_sh_optim_bits               = 32;
    // Single SH quantization level. Collapses the FPBO dispatch's
    // (sh_optim_bits, sh_value_bits) axes down to 2 instantiations:
    //   0 = off    : fp32 everywhere
    //   1 = light  : 16-bit SH value, 8x2-bit SH optim, 16x2-bit non-SH optim
    // Python is expected to ALSO set sh_optim_bits / sh_value_bits /
    // non_sh_optim_bits to the values implied by this level (the FPBO
    // dispatcher only uses the level; non-FPBO paths and EngineState still
    // read the individual bits).
    int   quantization_level           = 0;
    bool  use_per_splat_bias_correction   = false;

    // When true, fold projection-backward and Adam-based per-splat optim into
    // a single fused kernel (FusedProjectionBwdOptim). The engine then skips
    // the standalone projection_*_backward + engine_optim_step calls. SH
    // quantization is not yet supported by the fused path; the MVP allocates
    // full fp32 g1/g2 momentum (no `sh_quant_state`) in this mode.
    bool  use_fused_proj_bwd_optim        = false;

    // When true, split the camera batch into one-camera sub-batches inside
    // engine_train_step. Forward + bilagrid/PPISP fwd + loss + raster/proj
    // bwd run once per sub-batch and atomicAdd into the per-splat grad
    // accumulators; a single optimizer + bilagrid/PPISP optim + densify pass
    // runs at the end. Inside the splat Adam kernels, the accumulated data
    // gradient is scaled by `1/B` (B = full batch size) before adding the
    // per-splat regularization terms, so per-image grad magnitude vs reg
    // weight is batch-size invariant. Frees the per-sub-batch screen-space
    // buffers (splats_s, aabb, depths, isect/flatten ids, render_Ts, raster
    // bwd v_splats_s) between sub-batches -- peak VRAM scales by ~1/B.
    //
    // Not compatible with use_fused_proj_bwd_optim or use_color_trust_region
    // in this turn; the engine throws when either is also set. Warped
    // training_step path also throws (would need per-input-image splitting).
    bool  split_batch      = false;

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
//
// run_before_bilagrid controls the forward order when both PPISP and bilagrid
// RGB are enabled:
//   false (default): render -> bilagrid -> PPISP -> loss.
//   true           : render -> PPISP    -> bilagrid -> loss.
// Backward hooks invert this automatically (the inner hook runs first).
// Ignored when only one of the two is enabled.
struct PpispStepConfig {
    float lr = 0.0f;
    std::array<float, (int)PPISPRegLossIndex::length> reg_weights{};
    bool  run_before_bilagrid = false;
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
