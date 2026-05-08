#pragma once

#include <ATen/Tensor.h>
#include <ATen/Device.h>

#include "common.cuh"


std::tuple<at::Tensor, std::optional<at::Tensor>, at::Tensor, at::Tensor, at::Tensor>
fused_ssim_forward(
    at::Tensor &img1,
    at::Tensor &img2,
    bool train,
    bool return_ssim_map,
    bool is_l1
);

std::tuple<at::Tensor, at::Tensor, at::Tensor, at::Tensor>
fused_ssim_forward_inplace(
    at::Tensor &img1,
    at::Tensor &img2,
    bool train,
    float ssim_loss_map_weight,
    at::Tensor &ssim_loss_map,
    bool is_l1
);

at::Tensor
fused_ssim_backward(
    at::Tensor &img1,
    at::Tensor &img2,
    const float dL_dmap,
    std::optional<at::Tensor> &dm_dmu1,
    std::optional<at::Tensor> &dm_dsigma1_sq,
    std::optional<at::Tensor> &dm_dsigma12
);

void fused_ssim_backward_inplace(
    at::Tensor &img1,
    at::Tensor &img2,
    const float dL_dmap,
    std::optional<at::Tensor> &dm_dmu1,
    std::optional<at::Tensor> &dm_dsigma1_sq,
    std::optional<at::Tensor> &dm_dsigma12,
    at::Tensor &dL_dimg1
);

float fused_ssim_inplace(
    at::Tensor &img1,
    at::Tensor &img2,
    const float dL_dmap,
    at::Tensor &dL_dimg1,
    bool return_ssim_val,
    std::optional<at::Tensor> ssim_loss_map,
    float ssim_loss_map_weight
);
