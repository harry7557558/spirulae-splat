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
