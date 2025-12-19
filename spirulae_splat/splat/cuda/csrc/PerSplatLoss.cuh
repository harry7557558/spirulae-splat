#pragma once

#include <ATen/Tensor.h>


static constexpr uint kNumPerSplatLosses = 5;


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



at::Tensor compute_per_splat_losses_forward_tensor(
    at::Tensor &scales,  // [N, 3] or [N, 2]
    at::Tensor &opacities,  // [N, 1]
    at::Tensor &quats,  // [N, 4]
    float mcmc_opacity_reg_weight,
    float mcmc_scale_reg_weight,
    float max_gauss_ratio,
    float scale_regularization_weight,
    float erank_reg_weight,
    float erank_reg_weight_s3,
    float quat_norm_reg_weight
);


std::tuple<at::Tensor, at::Tensor, at::Tensor>
compute_per_splat_losses_backward_tensor(
    at::Tensor &scales,  // [N, 3] or [N, 2]
    at::Tensor &opacities,  // [N, 1]
    at::Tensor &quats,  // [N, 4]
    at::Tensor &v_losses,  // [kNumPerSplatLosses]
    float mcmc_opacity_reg_weight,
    float mcmc_scale_reg_weight,
    float max_gauss_ratio,
    float scale_regularization_weight,
    float erank_reg_weight,
    float erank_reg_weight_s3,
    float quat_norm_reg_weight
);


void mcmc_add_noise_3dgs_tensor(
    std::string primitive,
    float scaler, float min_opacity,
    at::Tensor &means,
    at::Tensor &scales,
    at::Tensor &quats,
    at::Tensor &opacs
);
