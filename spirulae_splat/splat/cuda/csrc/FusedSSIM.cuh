#pragma once

#include "common.cuh"
#include <Tensor.h>


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



float fused_ssim_forward(
    TorchTensorView img1,           // [B, H, W, 3]
    TorchTensorView img2,           // [B, H, W, 3]
    TorchTensorView dm_dmu1,        // [B, H, W, 3] output, or null
    TorchTensorView dm_dsigma1_sq,  // [B, H, W, 3] output, or null
    TorchTensorView dm_dsigma12,    // [B, H, W, 3] output, or null
    TorchTensorView ssim_loss_map,  // [B, H, W, 1] output, or null
    float ssim_loss_map_weight,
    bool inplace,
    bool is_l1
);


void fused_ssim_backward(
    TorchTensorView img1,           // [B, H, W, 3]
    TorchTensorView img2,           // [B, H, W, 3]
    float dL_dmap,
    TorchTensorView dm_dmu1,        // [B, H, W, 3] input, or null
    TorchTensorView dm_dsigma1_sq,  // [B, H, W, 3] input, or null
    TorchTensorView dm_dsigma12,    // [B, H, W, 3] input, or null
    TorchTensorView dL_dimg1,       // [B, H, W, 3] output
    bool inplace
);


float fused_ssim_inplace(
    TorchTensorView img1,           // [B, H, W, 3]
    TorchTensorView img2,           // [B, H, W, 3]
    TorchTensorView mask,           // [B, H, W] bool, or null
    float dL_dmap,
    TorchTensorView dL_dimg1,       // [B, H, W, 3] output (accumulated)
    bool return_ssim_val,
    TorchTensorView ssim_loss_map,  // [B, H, W, 1] output, or null
    float ssim_loss_map_weight
);


float fused_ssim_inplace_async(
    TorchTensorView img1,
    TorchTensorView img2,
    TorchTensorView mask,
    float dL_dmap,
    TorchTensorView dL_dimg1,
    TorchTensorView ssim_loss_map,
    float ssim_loss_map_weight,
    AsyncReadout<float>& readout
);
