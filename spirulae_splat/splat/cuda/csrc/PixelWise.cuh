#pragma once

#include <torch/types.h>


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



torch::Tensor blend_background_forward_tensor(
    torch::Tensor &rgb,  // [H, W, 3]
    torch::Tensor &alpha,  // [H, W, 1]
    torch::Tensor &background  // [H, W, 3]
);


std::tuple<torch::Tensor, torch::Tensor, torch::Tensor>
blend_background_backward_tensor(
    torch::Tensor &rgb,  // [H, W, 3]
    torch::Tensor &alpha,  // [H, W, 1]
    torch::Tensor &background,  // [H, W, 3]
    torch::Tensor &v_out_rgb  // [H, W, 3]
);
