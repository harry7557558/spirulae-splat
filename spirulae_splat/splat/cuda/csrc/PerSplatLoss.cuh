#pragma once

#include <torch/types.h>


static constexpr uint kNumPerSplatLosses = 5;


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



torch::Tensor compute_per_splat_losses_forward_tensor(
    torch::Tensor &scales,  // [N, 3] or [N, 2]
    torch::Tensor &opacities,  // [N, 1]
    torch::Tensor &quats,  // [N, 4]
    float mcmc_opacity_reg_weight,
    float mcmc_scale_reg_weight,
    float max_gauss_ratio,
    float scale_regularization_weight,
    float erank_reg_weight,
    float erank_reg_weight_s3,
    float quat_norm_reg_weight
);


std::tuple<torch::Tensor, torch::Tensor, torch::Tensor>
compute_per_splat_losses_backward_tensor(
    torch::Tensor &scales,  // [N, 3] or [N, 2]
    torch::Tensor &opacities,  // [N, 1]
    torch::Tensor &quats,  // [N, 4]
    torch::Tensor &v_losses,  // [kNumPerSplatLosses]
    float mcmc_opacity_reg_weight,
    float mcmc_scale_reg_weight,
    float max_gauss_ratio,
    float scale_regularization_weight,
    float erank_reg_weight,
    float erank_reg_weight_s3,
    float quat_norm_reg_weight
);
