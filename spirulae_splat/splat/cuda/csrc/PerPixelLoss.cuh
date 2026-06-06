#pragma once

#include <Tensor.h>

enum class RawLossIndex {
    RgbLoss,
    RgbL2,
    DepthSupX,
    DepthSupY,
    DepthSupXX,
    DepthSupYY,
    DepthSupXY,
    RenderNormalSup,
    DepthNormalSup,
    AlphaSup,
    AlphaSupUnder,
    NormalReg,
    AlphaReg,
    RgbDistReg,
    DepthDistReg,
    NormalDistReg,
    PixelsTotal,
    MaskTotal,
    DepthMaskTotal,
    RenderNormalMaskTotal,
    DepthNormalMaskTotal,
    NormalRegMaskTotal,
    AlphaMaskTotal,
    length
};

enum class LossWeightIndex {
    RgbSupL1,
    RgbSupL2,
    // YUV (BT.601) per-pixel supervision; folded into RgbLoss in the slang
    // kernel so it shares the MaskTotal divisor. Keep the order in sync with
    // per_pixel_losses.slang::LossWeightIndex.
    YSupL1,
    YSupL2,
    USupL2,
    VSupL2,
    DepthSup,
    NormalSup,
    AlphaSup,
    AlphaSupUnder,
    NormalReg,
    AlphaReg,
    RgbDistReg,
    DepthDistReg,
    NormalDistReg,
    length
};

enum class LossIndex {
    RgbLoss,
    RgbPSNR,
    DepthSup,
    NormalSup,
    AlphaSup,
    NormalReg,
    AlphaReg,
    RgbDistReg,
    DepthDistReg,
    NormalDistReg,
    length
};

struct PerPixelGrads {
    TorchTensorView v_render_rgb, v_ref_rgb, v_render_depth, v_ref_depth;
    TorchTensorView v_render_normal, v_depth_normal, v_ref_normal;
    TorchTensorView v_render_Ts, v_rgb_dist, v_depth_dist, v_normal_dist;
};

// Per-step loss values returned by compute_multi_scale_per_pixel_losses.
// Fields mirror LossIndex (raw weighted sums) plus ssim.
struct LossValues {
    float rgb_loss;
    float rgb_psnr;
    float depth_sup;
    float normal_sup;
    float alpha_sup;
    float normal_reg;
    float alpha_reg;
    float rgb_dist_reg;
    float depth_dist_reg;
    float normal_dist_reg;
    float ssim;
};


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



LossValues compute_multi_scale_per_pixel_losses(
    int num_loss_scales,
    TorchTensorView render_rgb,
    TorchTensorView ref_rgb,
    TorchTensorView render_depth,
    TorchTensorView ref_depth,
    TorchTensorView render_normal,
    TorchTensorView depth_normal,
    TorchTensorView ref_normal,
    TorchTensorView render_Ts,
    TorchTensorView rgb_dist,
    TorchTensorView depth_dist,
    TorchTensorView normal_dist,
    TorchTensorView ref_alpha,
    bool has_mask,
    const std::array<float, (int)LossWeightIndex::length> loss_weights_0,
    const float w_ssim,
    TorchTensorView v_losses,
    std::vector<bool> needs_input_grad,
    long num_train_images,
    TorchTensorView camera_indices,
    TorchTensorView loss_map_out,
    bool structure_only_loss_map,
    PerPixelGrads& grads_out
);
