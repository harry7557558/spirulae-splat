#pragma once

#include <torch/types.h>


enum class RawLossIndex {
    RgbL1,
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
    length
};

enum class LossWeightIndex {
    RgbSup,
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
    RgbL1,
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


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



std::tuple<torch::Tensor, torch::Tensor>
compute_per_pixel_losses_forward_tensor(
    std::optional<at::Tensor> render_rgb,
    std::optional<at::Tensor> ref_rgb,
    std::optional<at::Tensor> render_depth,
    std::optional<at::Tensor> ref_depth,
    std::optional<at::Tensor> render_normal,
    std::optional<at::Tensor> depth_normal,
    std::optional<at::Tensor> ref_normal,
    std::optional<at::Tensor> render_alpha,
    std::optional<at::Tensor> rgb_dist,
    std::optional<at::Tensor> depth_dist,
    std::optional<at::Tensor> normal_dist,
    std::optional<at::Tensor> ref_alpha,
    std::optional<at::Tensor> mask,
    std::optional<at::Tensor> depth_mask,
    std::optional<at::Tensor> normal_mask,
    const std::array<float, (uint)LossIndex::length> loss_weights_0,
    long num_train_images,
    std::optional<at::Tensor> camera_indices
);


std::tuple<
    std::optional<at::Tensor>, // render_rgb
    std::optional<at::Tensor>, // ref_rgb
    std::optional<at::Tensor>, // render_depth
    std::optional<at::Tensor>, // ref_depth
    std::optional<at::Tensor>, // render_normal
    std::optional<at::Tensor>, // depth_normal
    std::optional<at::Tensor>, // ref_normal
    std::optional<at::Tensor>, // render_alpha
    std::optional<at::Tensor>, // rgb_dist
    std::optional<at::Tensor>, // depth_dist
    std::optional<at::Tensor> // normal_dist
> compute_per_pixel_losses_backward_tensor(
    std::optional<at::Tensor> render_rgb,
    std::optional<at::Tensor> ref_rgb,
    std::optional<at::Tensor> render_depth,
    std::optional<at::Tensor> ref_depth,
    std::optional<at::Tensor> render_normal,
    std::optional<at::Tensor> depth_normal,
    std::optional<at::Tensor> ref_normal,
    std::optional<at::Tensor> render_alpha,
    std::optional<at::Tensor> rgb_dist,
    std::optional<at::Tensor> depth_dist,
    std::optional<at::Tensor> normal_dist,
    std::optional<at::Tensor> ref_alpha,
    std::optional<at::Tensor> mask,
    std::optional<at::Tensor> depth_mask,
    std::optional<at::Tensor> normal_mask,
    at::Tensor raw_losses,
    const std::array<float, (uint)LossIndex::length> loss_weights_0,
    at::Tensor v_losses,
    long num_train_images,
    std::optional<at::Tensor> camera_indices
);
