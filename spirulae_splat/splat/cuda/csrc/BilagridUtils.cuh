#pragma once

#include <ATen/Tensor.h>
#include <ATen/Device.h>

#include "common.cuh"


std::tuple<at::Tensor, std::optional<at::Tensor>, at::Tensor, at::Tensor, at::Tensor>
fused_ssim_forward(
    at::Tensor &img1,
    at::Tensor &img2,
    bool train,
    bool return_ssim_map
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


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



at::Tensor dct3d_type1_ortho_tensor(at::Tensor input);
