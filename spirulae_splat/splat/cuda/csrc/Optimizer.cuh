#pragma once

#include <ATen/Tensor.h>
#include <c10/core/TensorOptions.h>
#include <ATen/Device.h>


#include "types.cuh"


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



void fused_adam(
    at::Tensor param,
    at::Tensor grad,
    at::Tensor exp_avg,
    at::Tensor exp_avg_sq,
    float lr,
    float beta1,
    float beta2,
    float eps,
    int step
);


void offloaded_adam(
    at::Tensor param,      // Device
    at::Tensor grad,       // Device
    at::Tensor exp_avg,    // Host
    at::Tensor exp_avg_sq, // Host
    float lr, float beta1, float beta2, float eps, int step
);


void semi_offloaded_adam(
    at::Tensor param,      // Device
    at::Tensor grad,       // Device
    at::Tensor exp_avg,    // Host
    at::Tensor exp_avg_sq, // Host
    float lr, float beta1, float beta2, float eps, int step
);


void fused_adam_multi(
    std::vector<at::Tensor> params,
    std::vector<at::Tensor> grads,
    std::vector<at::Tensor> exp_avgs,
    std::vector<at::Tensor> exp_avg_sqs,
    float lr,
    float beta1,
    float beta2,
    float eps,
    int step
);


void fused_adam_riemannian_quat(
    at::Tensor param,
    at::Tensor grad,
    at::Tensor exp_avg,
    at::Tensor exp_avg_sq,
    float lr,
    float beta1,
    float beta2,
    float eps,
    int step
);


void fused_newton(
    at::Tensor param,
    at::Tensor grad,
    at::Tensor hess_diag,
    at::Tensor exp_avg,
    at::Tensor exp_avg_sq,
    float lr,
    float beta1,
    float beta2,
    float eps,
    int step1,
    int step2
);


void fused_newton_multi(
    std::vector<at::Tensor> params,
    std::vector<at::Tensor> grads,
    std::vector<at::Tensor> hess_diags,
    std::vector<at::Tensor> exp_avgs,
    std::vector<at::Tensor> exp_avg_sqs,
    float lr,
    float beta1,
    float beta2,
    float eps,
    int step1,
    int step2
);


void fused_adam_scale_agnostic_mean(
    at::Tensor param,
    at::Tensor grad,
    at::Tensor exp_avg,
    at::Tensor exp_avg_sq,
    at::Tensor scales,
    at::Tensor quats,
    at::Tensor opacities,
    at::Tensor radii,
    float lr,
    float beta1,
    float beta2,
    float eps,
    float eps_tr,
    int step
);


void fused_optim_3dgs_geometry(
    at::Tensor means,
    at::Tensor v_means,
    at::Tensor g1_means,
    at::Tensor g2_means,
    at::Tensor quats,
    at::Tensor v_quats,
    at::Tensor g1_quats,
    at::Tensor g2_quats,
    at::Tensor scales,
    at::Tensor v_scales,
    at::Tensor g1_scales,
    at::Tensor g2_scales,
    at::Tensor opacities,
    at::Tensor v_opacities,
    at::Tensor g1_opacities,
    at::Tensor g2_opacities,
    at::Tensor radii,
    const float lr_means,
    const float lr_quats,
    const float lr_scales,
    const float lr_opacs,
    const float mcmc_noise_lr,
    const float min_opacity,
    const float max_gauss_ratio,
    const float scale_regularization_weight,
    const float mcmc_opacity_reg_weight,
    const float mcmc_scale_reg_weight,
    const float erank_reg_weight,
    const float erank_reg_weight_s3,
    const float quat_norm_reg_weight,
    bool use_scale_agnostic_mean,
    std::variant<int32_t, at::Tensor> step
);


void fused_adam_with_steps_tensor(
    at::Tensor param,
    at::Tensor grad,
    at::Tensor exp_avg,
    at::Tensor exp_avg_sq,
    float lr,
    std::variant<int32_t, at::Tensor> step
);


void fused_3dgs2tr_mean_optim(
    at::Tensor means,
    at::Tensor vr_means,
    at::Tensor h_means,
    at::Tensor scales,
    at::Tensor quats,
    at::Tensor opacities,
    at::Tensor exp_avg_means,
    at::Tensor exp_avg_sq_means,
    float lr,
    float beta1,
    float beta2,
    float eps,
    float eps_tr,
    int step1,
    int step2
);


void fused_3dgs2tr_scale_optim(
    at::Tensor scales,
    at::Tensor vr_scales,
    at::Tensor h_scales,
    at::Tensor opacities,
    at::Tensor exp_avg_scales,
    at::Tensor exp_avg_sq_scales,
    float lr,
    float beta1,
    float beta2,
    float eps,
    float eps_tr,
    int step1,
    int step2
);


void fused_3dgs2tr_color_optim(
    at::Tensor colors,
    at::Tensor vr_colors,
    at::Tensor h_colors,
    at::Tensor opacities,
    at::Tensor exp_avg_colors,
    at::Tensor exp_avg_sq_colors,
    float lr,
    float beta1,
    float beta2,
    float eps,
    float eps_tr,
    int step1,
    int step2
);


void fused_3dgs2tr_opacity_optim(
    at::Tensor opacities,
    at::Tensor vr_opacities,
    at::Tensor h_opacities,
    at::Tensor exp_avg_opacities,
    at::Tensor exp_avg_sq_opacities,
    float lr,
    float beta1,
    float beta2,
    float eps,
    float eps_tr,
    int step1,
    int step2
);


void fused_3dgs2tr_quat_optim(
    at::Tensor quats,
    at::Tensor vr_quats,
    at::Tensor h_quats,
    at::Tensor scales,
    at::Tensor opacities,
    at::Tensor exp_avg_quats,
    at::Tensor exp_avg_sq_quats,
    float lr,
    float beta1,
    float beta2,
    float eps,
    float eps_tr,
    int step1,
    int step2
);


void fused_adam_linear_rgb_optim(
    at::Tensor param,
    at::Tensor grad,
    at::Tensor exp_avg,
    at::Tensor exp_avg_sq,
    float lr,
    float beta1,
    float beta2,
    float eps,
    int step
);


void fused_adamtr_linear_rgb_optim(
    at::Tensor param,
    at::Tensor grad,
    at::Tensor exp_avg,
    at::Tensor exp_avg_sq,
    at::Tensor opacities,
    float lr,
    float beta1,
    float beta2,
    float eps,
    float eps_tr,
    int step
);


void fused_adamtr_rgb_optim(
    at::Tensor param,
    at::Tensor grad,
    at::Tensor exp_avg,
    at::Tensor exp_avg_sq,
    at::Tensor opacities,
    float lr,
    float beta1,
    float beta2,
    float eps,
    float eps_tr,
    int step
);


void fused_adamtr_linear_rgb_sh_optim(
    at::Tensor param,
    at::Tensor grad,
    at::Tensor exp_avg,
    at::Tensor exp_avg_sq,
    at::Tensor colors,
    at::Tensor opacities,
    float lr,
    float beta1,
    float beta2,
    float eps,
    float eps_tr,
    int step
);


void fused_adamtr_rgb_sh_optim(
    at::Tensor param,
    at::Tensor grad,
    at::Tensor exp_avg,
    at::Tensor exp_avg_sq,
    at::Tensor colors,
    at::Tensor opacities,
    float lr,
    float beta1,
    float beta2,
    float eps,
    float eps_tr,
    int step
);
