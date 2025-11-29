#pragma once

#include <torch/types.h>

#ifdef __CUDACC__
#define TensorView _Slang_TensorView
#include "generated/slang_all.cu"
#undef TensorView
#endif

#include "types.cuh"


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



torch::Tensor blend_background_forward_tensor(
    torch::Tensor &rgb,  // [B, H, W, 3]
    torch::Tensor &alpha,  // [B, H, W, 1]
    torch::Tensor &background  // [B, H, W, 3]
);


std::tuple<torch::Tensor, torch::Tensor, torch::Tensor>
blend_background_backward_tensor(
    torch::Tensor &rgb,  // [B, H, W, 3]
    torch::Tensor &alpha,  // [B, H, W, 1]
    torch::Tensor &background,  // [B, H, W, 3]
    torch::Tensor &v_out_rgb  // [B, H, W, 3]
);


torch::Tensor log_map_image_forward_tensor(
    torch::Tensor &rgb,  // [B, H, W, 3]
    float t
);


torch::Tensor log_map_image_backward_tensor(
    torch::Tensor &rgb,  // [B, H, W, 3]
    float t,
    torch::Tensor &v_out_rgb  // [B, H, W, 3]
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


torch::Tensor ray_depth_to_linear_depth_forward_tensor(
    gsplat::CameraModelType camera_model,
    torch::Tensor Ks,  // [B, 3, 3]
    CameraDistortionCoeffsTensor dist_coeffs,
    torch::Tensor depths  // [B, H, W, 1]
);


torch::Tensor ray_depth_to_linear_depth_backward_tensor(
    gsplat::CameraModelType camera_model,
    torch::Tensor Ks,  // [B, 3, 3]
    CameraDistortionCoeffsTensor dist_coeffs,
    torch::Tensor v_out_depths  // [B, H, W, 1]
);


torch::Tensor distort_image_tensor(
    gsplat::CameraModelType camera_model,
    torch::Tensor Ks,  // [B, 3, 3]
    CameraDistortionCoeffsTensor dist_coeffs,
    torch::Tensor in_image  // [B, H, W, C]
);


torch::Tensor undistort_image_tensor(
    gsplat::CameraModelType camera_model,
    torch::Tensor Ks,  // [B, 3, 3]
    CameraDistortionCoeffsTensor dist_coeffs,
    torch::Tensor in_image  // [B, H, W, C]
);
