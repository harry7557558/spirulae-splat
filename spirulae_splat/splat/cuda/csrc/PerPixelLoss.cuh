#pragma once

#include <ATen/Tensor.h>
#include <c10/core/TensorOptions.h>
#include <ATen/Device.h>


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


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



std::tuple<at::Tensor, at::Tensor, std::optional<at::Tensor>>
compute_per_pixel_losses_forward_tensor(
    std::optional<at::Tensor> render_rgb,
    std::optional<at::Tensor> ref_rgb,
    std::optional<at::Tensor> render_depth,
    std::optional<at::Tensor> ref_depth,
    std::optional<at::Tensor> render_normal,
    std::optional<at::Tensor> depth_normal,
    std::optional<at::Tensor> ref_normal,
    std::optional<at::Tensor> render_Ts,
    std::optional<at::Tensor> rgb_dist,
    std::optional<at::Tensor> depth_dist,
    std::optional<at::Tensor> normal_dist,
    std::optional<at::Tensor> ref_alpha,
    std::optional<at::Tensor> mask,
    std::optional<at::Tensor> depth_mask,
    std::optional<at::Tensor> normal_mask,
    std::optional<at::Tensor> alpha_mask,
    const std::array<float, (int)LossWeightIndex::length> loss_weights_0,
    long num_train_images,
    std::optional<at::Tensor> camera_indices,
    bool return_loss_map
);


std::tuple<  // returns gradients
    std::optional<at::Tensor>, // render_rgb
    std::optional<at::Tensor>, // ref_rgb
    std::optional<at::Tensor>, // render_depth
    std::optional<at::Tensor>, // ref_depth
    std::optional<at::Tensor>, // render_normal
    std::optional<at::Tensor>, // depth_normal
    std::optional<at::Tensor>, // ref_normal
    std::optional<at::Tensor>, // render_Ts
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
    std::optional<at::Tensor> render_Ts,
    std::optional<at::Tensor> rgb_dist,
    std::optional<at::Tensor> depth_dist,
    std::optional<at::Tensor> normal_dist,
    std::optional<at::Tensor> ref_alpha,
    std::optional<at::Tensor> mask,
    std::optional<at::Tensor> depth_mask,
    std::optional<at::Tensor> normal_mask,
    std::optional<at::Tensor> alpha_mask,
    at::Tensor raw_losses,
    const std::array<float, (int)LossWeightIndex::length> loss_weights_0,
    at::Tensor v_losses,
    std::vector<bool> needs_input_grad,
    long num_train_images,
    std::optional<at::Tensor> camera_indices
);


std::tuple<
    at::Tensor,  // mean losses
    std::optional<at::Tensor>,  // loss map
    std::optional<std::tuple<
        std::optional<at::Tensor>, // render_rgb
        std::optional<at::Tensor>, // ref_rgb
        std::optional<at::Tensor>, // render_depth
        std::optional<at::Tensor>, // ref_depth
        std::optional<at::Tensor>, // render_normal
        std::optional<at::Tensor>, // depth_normal
        std::optional<at::Tensor>, // ref_normal
        std::optional<at::Tensor>, // render_Ts
        std::optional<at::Tensor>, // rgb_dist
        std::optional<at::Tensor>, // depth_dist
        std::optional<at::Tensor> // normal_dist
    >>,
    std::tuple<
        float,  // psnr value
        float  // ssim value
    >
> compute_multi_scale_per_pixel_losses_tensor(
    int num_loss_scales,
    std::optional<at::Tensor> render_rgb,
    std::optional<at::Tensor> ref_rgb,
    std::optional<at::Tensor> render_depth,
    std::optional<at::Tensor> ref_depth,
    std::optional<at::Tensor> render_normal,
    std::optional<at::Tensor> depth_normal,
    std::optional<at::Tensor> ref_normal,
    std::optional<at::Tensor> render_Ts,
    std::optional<at::Tensor> rgb_dist,
    std::optional<at::Tensor> depth_dist,
    std::optional<at::Tensor> normal_dist,
    std::optional<at::Tensor> ref_alpha,
    std::optional<at::Tensor> mask,
    std::optional<at::Tensor> depth_mask,
    std::optional<at::Tensor> normal_mask,
    std::optional<at::Tensor> alpha_mask,
    const std::array<float, (int)LossWeightIndex::length> loss_weights_0,
    const float w_ssim,
    std::optional<at::Tensor> v_losses,
    std::vector<bool> needs_input_grad,
    long num_train_images,
    std::optional<at::Tensor> camera_indices,
    bool return_loss_map
);
