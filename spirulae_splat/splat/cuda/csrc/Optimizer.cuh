#pragma once

#include <Tensor.h>
#include "types.cuh"


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



void set_zero_tensor(TorchTensorView x);


void fused_adam(
    TorchTensorView param,
    TorchTensorView grad,
    TorchTensorView exp_avg,
    TorchTensorView exp_avg_sq,
    float lr,
    float beta1,
    float beta2,
    float eps,
    int step
);


void offloaded_adam(
    TorchTensorView param,      // Device
    TorchTensorView grad,       // Device
    TorchTensorView exp_avg,    // Host
    TorchTensorView exp_avg_sq, // Host
    float lr, float beta1, float beta2, float eps, int step
);


void semi_offloaded_adam(
    TorchTensorView param,      // Device
    TorchTensorView grad,       // Device
    TorchTensorView exp_avg,    // Host
    TorchTensorView exp_avg_sq, // Host
    float lr, float beta1, float beta2, float eps, int step
);


void fused_adam_multi(
    std::vector<TorchTensorView> params,
    std::vector<TorchTensorView> grads,
    std::vector<TorchTensorView> exp_avgs,
    std::vector<TorchTensorView> exp_avg_sqs,
    float lr,
    float beta1,
    float beta2,
    float eps,
    int step
);


void fused_adam_riemannian_quat(
    TorchTensorView param,
    TorchTensorView grad,
    TorchTensorView exp_avg,
    TorchTensorView exp_avg_sq,
    float lr,
    float beta1,
    float beta2,
    float eps,
    int step
);


void fused_newton(
    TorchTensorView param,
    TorchTensorView grad,
    TorchTensorView hess_diag,
    TorchTensorView exp_avg,
    TorchTensorView exp_avg_sq,
    float lr,
    float beta1,
    float beta2,
    float eps,
    int step1,
    int step2
);


void fused_newton_multi(
    std::vector<TorchTensorView> params,
    std::vector<TorchTensorView> grads,
    std::vector<TorchTensorView> hess_diags,
    std::vector<TorchTensorView> exp_avgs,
    std::vector<TorchTensorView> exp_avg_sqs,
    float lr,
    float beta1,
    float beta2,
    float eps,
    int step1,
    int step2
);


void fused_adam_scale_agnostic_mean(
    TorchTensorView param,
    TorchTensorView grad,
    TorchTensorView exp_avg,
    TorchTensorView exp_avg_sq,
    TorchTensorView scales,
    TorchTensorView quats,
    TorchTensorView opacities,
    TorchTensorView radii,
    float lr,
    float beta1,
    float beta2,
    float eps,
    float eps_tr,
    int step
);


void fused_optim_3dgs_geometry(
    int64_t num_splats,
    TorchTensorView means,
    TorchTensorView v_means,
    TorchTensorView g1_means,
    TorchTensorView g2_means,
    TorchTensorView quats,
    TorchTensorView v_quats,
    TorchTensorView g1_quats,
    TorchTensorView g2_quats,
    TorchTensorView scales,
    TorchTensorView v_scales,
    TorchTensorView g1_scales,
    TorchTensorView g2_scales,
    TorchTensorView opacities,
    TorchTensorView v_opacities,
    TorchTensorView g1_opacities,
    TorchTensorView g2_opacities,
    TorchTensorView radii,
    const float lr_means,
    const float lr_quats,
    const float lr_scales,
    const float lr_opacs,
    const float max_gauss_ratio,
    const float scale_regularization_weight,
    const float mcmc_opacity_reg_weight,
    const float mcmc_scale_reg_weight,
    const float erank_reg_weight,
    const float erank_reg_weight_s3,
    const float quat_norm_reg_weight,
    bool use_scale_agnostic_mean,
    std::variant<int32_t, TorchTensorView> step
);


void fused_adam_with_steps_tensor(
    uint64_t num_splats,
    TorchTensorView param,
    TorchTensorView grad,
    TorchTensorView exp_avg,
    TorchTensorView exp_avg_sq,
    float lr,
    std::variant<int32_t, TorchTensorView> step,
    float l2_reg,
    float l2_reg_offset
);


void fused_adam_with_steps_8bit_tensor(
    uint64_t num_splats,
    TorchTensorView param,
    TorchTensorView grad,
    TorchTensorView exp_avg,
    TorchTensorView exp_avg_sq,
    TorchTensorView quant_bounds,
    float lr,
    std::variant<int32_t, TorchTensorView> step,
    float l2_reg,
    float l2_reg_offset
);


void fused_3dgs2tr_mean_optim(
    TorchTensorView means,
    TorchTensorView vr_means,
    TorchTensorView h_means,
    TorchTensorView scales,
    TorchTensorView quats,
    TorchTensorView opacities,
    TorchTensorView exp_avg_means,
    TorchTensorView exp_avg_sq_means,
    float lr,
    float beta1,
    float beta2,
    float eps,
    float eps_tr,
    int step1,
    int step2
);


void fused_3dgs2tr_scale_optim(
    TorchTensorView scales,
    TorchTensorView vr_scales,
    TorchTensorView h_scales,
    TorchTensorView opacities,
    TorchTensorView exp_avg_scales,
    TorchTensorView exp_avg_sq_scales,
    float lr,
    float beta1,
    float beta2,
    float eps,
    float eps_tr,
    int step1,
    int step2
);


void fused_3dgs2tr_color_optim(
    TorchTensorView colors,
    TorchTensorView vr_colors,
    TorchTensorView h_colors,
    TorchTensorView opacities,
    TorchTensorView exp_avg_colors,
    TorchTensorView exp_avg_sq_colors,
    float lr,
    float beta1,
    float beta2,
    float eps,
    float eps_tr,
    int step1,
    int step2
);


void fused_3dgs2tr_opacity_optim(
    TorchTensorView opacities,
    TorchTensorView vr_opacities,
    TorchTensorView h_opacities,
    TorchTensorView exp_avg_opacities,
    TorchTensorView exp_avg_sq_opacities,
    float lr,
    float beta1,
    float beta2,
    float eps,
    float eps_tr,
    int step1,
    int step2
);


void fused_3dgs2tr_quat_optim(
    TorchTensorView quats,
    TorchTensorView vr_quats,
    TorchTensorView h_quats,
    TorchTensorView scales,
    TorchTensorView opacities,
    TorchTensorView exp_avg_quats,
    TorchTensorView exp_avg_sq_quats,
    float lr,
    float beta1,
    float beta2,
    float eps,
    float eps_tr,
    int step1,
    int step2
);


void fused_adam_linear_rgb_optim(
    TorchTensorView param,
    TorchTensorView grad,
    TorchTensorView exp_avg,
    TorchTensorView exp_avg_sq,
    float lr,
    float beta1,
    float beta2,
    float eps,
    int step
);


void fused_adamtr_linear_rgb_optim(
    TorchTensorView param,
    TorchTensorView grad,
    TorchTensorView exp_avg,
    TorchTensorView exp_avg_sq,
    TorchTensorView opacities,
    float lr,
    float beta1,
    float beta2,
    float eps,
    float eps_tr,
    int step
);


void fused_adamtr_rgb_optim(
    TorchTensorView param,
    TorchTensorView grad,
    TorchTensorView exp_avg,
    TorchTensorView exp_avg_sq,
    TorchTensorView opacities,
    float lr,
    float beta1,
    float beta2,
    float eps,
    float eps_tr,
    int step
);


void fused_adamtr_linear_rgb_sh_optim(
    TorchTensorView param,
    TorchTensorView grad,
    TorchTensorView exp_avg,
    TorchTensorView exp_avg_sq,
    TorchTensorView colors,
    TorchTensorView opacities,
    float lr,
    float beta1,
    float beta2,
    float eps,
    float eps_tr,
    int step
);


void fused_adamtr_rgb_sh_optim(
    TorchTensorView param,
    TorchTensorView grad,
    TorchTensorView exp_avg,
    TorchTensorView exp_avg_sq,
    TorchTensorView colors,
    TorchTensorView opacities,
    float lr,
    float beta1,
    float beta2,
    float eps,
    float eps_tr,
    int step
);


void engine_optim_3dgs_geometry(
    int64_t num_splats,
    TorchTensorView means, TorchTensorView v_means, TorchTensorView g1_means, TorchTensorView g2_means,
    TorchTensorView quats, TorchTensorView v_quats, TorchTensorView g1_quats, TorchTensorView g2_quats,
    TorchTensorView scales, TorchTensorView v_scales, TorchTensorView g1_scales, TorchTensorView g2_scales,
    TorchTensorView opacities, TorchTensorView v_opacities, TorchTensorView g1_opacities, TorchTensorView g2_opacities,
    TorchTensorView radii,
    float lr_means, float lr_quats, float lr_scales, float lr_opacs,
    float max_gauss_ratio, float scale_regularization_weight,
    float mcmc_opacity_reg_weight, float mcmc_scale_reg_weight,
    float erank_reg_weight, float erank_reg_weight_s3, float quat_norm_reg_weight,
    bool use_scale_agnostic_mean,
    int32_t step
);


void engine_adam_step(
    int64_t num_splats,
    TorchTensorView param,
    TorchTensorView grad,
    TorchTensorView exp_avg,
    TorchTensorView exp_avg_sq,
    float lr,
    int32_t step,
    float l2_reg,
    float l2_reg_offset
);


void engine_adam_step_8bit(
    int64_t num_splats,
    TorchTensorView param,
    TorchTensorView grad,
    TorchTensorView exp_avg,
    TorchTensorView exp_avg_sq,
    TorchTensorView quant_bounds,
    float lr,
    int32_t step,
    float l2_reg,
    float l2_reg_offset
);
