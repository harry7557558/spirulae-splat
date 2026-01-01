#pragma once

#include <ATen/Tensor.h>
#include <ATen/Device.h>

#include "common.cuh"


std::tuple<at::Tensor,at::Tensor,at::Tensor,at::Tensor>
fused_ssim_forward(
    at::Tensor &img1,
    at::Tensor &img2,
    bool train
);

at::Tensor
fused_ssim_backward(
    at::Tensor &img1,
    at::Tensor &img2,
    at::Tensor &dL_dmap,
    at::Tensor &dm_dmu1,
    at::Tensor &dm_dsigma1_sq,
    at::Tensor &dm_dsigma12
);
