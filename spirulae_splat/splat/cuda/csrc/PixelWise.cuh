#pragma once

#include <ATen/Tensor.h>

#ifdef __CUDACC__
#define TensorView _Slang_TensorView
#include "generated/slang_all.cu"
#undef TensorView
#endif

#include "types.cuh"


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



at::Tensor blend_background_forward_tensor(
    at::Tensor &rgb,  // [B, H, W, 3]
    at::Tensor &alpha,  // [B, H, W, 1]
    at::Tensor &background  // [B, H, W, 3]
);


std::tuple<at::Tensor, at::Tensor, at::Tensor>
blend_background_backward_tensor(
    at::Tensor &rgb,  // [B, H, W, 3]
    at::Tensor &alpha,  // [B, H, W, 1]
    at::Tensor &background,  // [B, H, W, 3]
    at::Tensor &v_out_rgb  // [B, H, W, 3]
);


at::Tensor log_map_image_forward_tensor(
    at::Tensor &rgb,  // [B, H, W, 3]
    float t
);


at::Tensor log_map_image_backward_tensor(
    at::Tensor &rgb,  // [B, H, W, 3]
    float t,
    at::Tensor &v_out_rgb  // [B, H, W, 3]
);


at::Tensor depth_to_normal_forward_tensor(
    gsplat::CameraModelType camera_model,
    at::Tensor Ks,  // [B, 3, 3]
    CameraDistortionCoeffsTensor dist_coeffs,
    bool is_ray_depth,
    at::Tensor depths  // [B, H, W, 1]
);


at::Tensor depth_to_normal_backward_tensor(
    gsplat::CameraModelType camera_model,
    at::Tensor Ks,  // [B, 3, 3]
    CameraDistortionCoeffsTensor dist_coeffs,
    bool is_ray_depth,
    at::Tensor depths,  // [B, H, W, 1]
    at::Tensor v_normals  // [B, H, W, 3]
);


at::Tensor ray_depth_to_linear_depth_forward_tensor(
    gsplat::CameraModelType camera_model,
    at::Tensor Ks,  // [B, 3, 3]
    CameraDistortionCoeffsTensor dist_coeffs,
    at::Tensor depths  // [B, H, W, 1]
);


at::Tensor ray_depth_to_linear_depth_backward_tensor(
    gsplat::CameraModelType camera_model,
    at::Tensor Ks,  // [B, 3, 3]
    CameraDistortionCoeffsTensor dist_coeffs,
    at::Tensor v_out_depths  // [B, H, W, 1]
);


at::Tensor distort_image_tensor(
    gsplat::CameraModelType camera_model,
    at::Tensor Ks,  // [B, 3, 3]
    CameraDistortionCoeffsTensor dist_coeffs,
    at::Tensor in_image  // [B, H, W, C]
);


at::Tensor undistort_image_tensor(
    gsplat::CameraModelType camera_model,
    at::Tensor Ks,  // [B, 3, 3]
    CameraDistortionCoeffsTensor dist_coeffs,
    at::Tensor in_image  // [B, H, W, C]
);
