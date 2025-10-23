#pragma once

#include <torch/types.h>

#include <gsplat/Common.h>


typedef std::tuple<
    std::optional<at::Tensor>,
    std::optional<at::Tensor>,
    std::optional<at::Tensor>
> CameraDistortionCoeffsTensor;

struct CameraDistortionCoeffsBuffer {
    float4* __restrict__ radial_coeffs;
    float2* __restrict__ tangential_coeffs;
    float2* __restrict__ thin_prism_coeffs;

    CameraDistortionCoeffsBuffer(const CameraDistortionCoeffsTensor &tensors);
};


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


torch::Tensor depth_to_normal_forward_tensor(
    gsplat::CameraModelType camera_model,
    torch::Tensor Ks,  // [B, 3, 3]
    CameraDistortionCoeffsTensor dist_coeffs,
    bool is_ray_depth,
    torch::Tensor depths  // [B, H, W, 1]
);


torch::Tensor depth_to_normal_backward_tensor(
    gsplat::CameraModelType camera_model,
    torch::Tensor Ks,  // [B, 3, 3]
    CameraDistortionCoeffsTensor dist_coeffs,
    bool is_ray_depth,
    torch::Tensor depths,  // [B, H, W, 1]
    torch::Tensor v_normals  // [B, H, W, 3]
);
