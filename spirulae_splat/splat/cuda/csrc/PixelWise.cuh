#pragma once

#include <ATen/Tensor.h>

#include "types.cuh"

enum class RawPPISPRegLossIndex {
    SumExposure,
    SumVignettingCrSquared,
    SumVignettingAlpha0Relu,
    SumVignettingAlpha1Relu,
    SumVignettingAlpha2Relu,
    SumVignettingCxChannelVariance,
    SumVignettingCyChannelVariance,
    SumVignettingAlpha0ChannelVariance,
    SumVignettingAlpha1ChannelVariance,
    SumVignettingAlpha2ChannelVariance,
    SumColorBx,
    SumColorBy,
    SumColorRx,
    SumColorRy,
    SumColorGx,
    SumColorGy,
    SumColorNx,
    SumColorNy,
    SumCRFToeChannelVariance,
    SumCRFShoulderChannelVariance,
    SumCRFGammaChannelVariance,
    SumCRFCenterChannelVariance,
    length
};

enum class RawPPISPRegLossIndexRQS {
    SumExposure,
    SumVignettingCrSquared,
    SumVignettingAlpha0Relu,
    SumVignettingAlpha1Relu,
    SumVignettingAlpha2Relu,
    SumVignettingCxChannelVariance,
    SumVignettingCyChannelVariance,
    SumVignettingAlpha0ChannelVariance,
    SumVignettingAlpha1ChannelVariance,
    SumVignettingAlpha2ChannelVariance,
    SumColorBx,
    SumColorBy,
    SumColorRx,
    SumColorRy,
    SumColorGx,
    SumColorGy,
    SumColorNx,
    SumColorNy,
    SumCRFG0ChannelVariance,
    SumCRFG1ChannelVariance,
    SumCRFX0ChannelVariance,
    SumCRFY0ChannelVariance,
    SumCRFGcChannelVariance,
    length
};

enum class PPISPRegLossIndex {
    ExposureMean,
    VignettingCenter,
    VignettingNonPositivity,
    VignettingChannelVariance,
    ColorMean,
    CRFChannelVariance,
    length
};


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



at::Tensor uint8_image_to_float_tensor(
    at::Tensor &img_in  // [B, H, W, C]
);


at::Tensor uint16_image_to_float_tensor(
    at::Tensor &img_in  // [B, H, W, C]
);


at::Tensor rendered_depth_to_expected_depth_forward_tensor(
    at::Tensor &depth,  // [B, H, W, 1]
    at::Tensor &transmittance  // [B, H, W, 1]
);


std::tuple<at::Tensor, at::Tensor>
rendered_depth_to_expected_depth_backward_tensor(
    at::Tensor &depth,  // [B, H, W, 1]
    at::Tensor &transmittance,  // [B, H, W, 1]
    at::Tensor &v_out_depth  // [B, H, W, 1]
);


at::Tensor blend_background_forward_tensor(
    at::Tensor &rgb,  // [B, H, W, 3]
    at::Tensor &transmittance,  // [B, H, W, 1]
    at::Tensor &background  // [B, H, W, 3]
);


std::tuple<at::Tensor, at::Tensor, at::Tensor>
blend_background_backward_tensor(
    at::Tensor &rgb,  // [B, H, W, 3]
    at::Tensor &transmittance,  // [B, H, W, 1]
    at::Tensor &background,  // [B, H, W, 3]
    at::Tensor &v_out_rgb  // [B, H, W, 3]
);


at::Tensor blend_background_noise_forward_tensor(
    bool is_linear,
    at::Tensor &rgb,  // [B, H, W, 3]
    at::Tensor &transmittance,  // [B, H, W, 1]
    float randomize_weight,
    uint32_t seed
);


std::tuple<at::Tensor, at::Tensor>
blend_background_noise_backward_tensor(
    bool is_linear,
    at::Tensor &rgb,  // [B, H, W, 3]
    at::Tensor &transmittance,  // [B, H, W, 1]
    float randomize_weight,
    uint32_t seed,
    at::Tensor &v_out_rgb  // [B, H, W, 3]
);


at::Tensor rgb_to_srgb_forward_tensor(
    bool is_input_linear,
    at::Tensor &rgb,  // [B, H, W, 3]
    at::Tensor &color_matrix   // [3, 3]
);


at::Tensor rgb_to_srgb_backward_tensor(
    bool is_input_linear,
    at::Tensor &rgb,  // [B, H, W, 3]
    at::Tensor &color_matrix,   // [3, 3]
    at::Tensor &v_out_rgb  // [B, H, W, 3]
);


at::Tensor depth_to_points_forward_tensor(
    std::string camera_model,
    at::Tensor intrins,  // fx, fy, cx, cy
    CameraDistortionCoeffsTensor dist_coeffs,
    bool is_ray_depth,
    at::Tensor depths  // [B, H, W, 1]
);


at::Tensor depth_to_points_backward_tensor(
    std::string camera_model,
    at::Tensor intrins,  // fx, fy, cx, cy
    CameraDistortionCoeffsTensor dist_coeffs,
    bool is_ray_depth,
    at::Tensor in_depths,  // [B, H, W, 1]
    at::Tensor v_out_points  // [B, H, W, 3]
);


at::Tensor depth_to_normal_forward_tensor(
    std::string camera_model,
    at::Tensor intrins,  // fx, fy, cx, cy
    CameraDistortionCoeffsTensor dist_coeffs,
    bool is_ray_depth,
    at::Tensor depths  // [B, H, W, 1]
);


at::Tensor depth_to_normal_backward_tensor(
    std::string camera_model,
    at::Tensor intrins,  // fx, fy, cx, cy
    CameraDistortionCoeffsTensor dist_coeffs,
    bool is_ray_depth,
    at::Tensor depths,  // [B, H, W, 1]
    at::Tensor v_normals  // [B, H, W, 3]
);


at::Tensor depth_normal_loss_forward_tensor(
    std::string camera_model,
    at::Tensor intrins,  // fx, fy, cx, cy
    CameraDistortionCoeffsTensor dist_coeffs,
    bool is_ray_depth,
    at::Tensor depths,  // [B, H, W, 1]
    at::Tensor gt_normals  // [B, H, W, 3]
);


std::tuple<at::Tensor, at::Tensor> depth_normal_loss_backward_tensor(
    std::string camera_model,
    at::Tensor intrins,  // fx, fy, cx, cy
    CameraDistortionCoeffsTensor dist_coeffs,
    bool is_ray_depth,
    at::Tensor depths,  // [B, H, W, 1]
    at::Tensor gt_normals,  // [B, H, W, 3]
    at::Tensor v_losses  // [B, H, W, 1]
);


at::Tensor ray_depth_to_linear_depth_forward_tensor(
    std::string camera_model,
    at::Tensor intrins,  // fx, fy, cx, cy
    CameraDistortionCoeffsTensor dist_coeffs,
    at::Tensor depths  // [B, H, W, 1]
);


at::Tensor ray_depth_to_linear_depth_backward_tensor(
    std::string camera_model,
    at::Tensor intrins,  // fx, fy, cx, cy
    CameraDistortionCoeffsTensor dist_coeffs,
    at::Tensor v_out_depths  // [B, H, W, 1]
);


at::Tensor distort_image_tensor(
    std::string camera_model,
    at::Tensor intrins,  // fx, fy, cx, cy
    CameraDistortionCoeffsTensor dist_coeffs,
    at::Tensor in_image  // [B, H, W, C]
);


at::Tensor undistort_image_tensor(
    std::string camera_model,
    at::Tensor intrins,  // fx, fy, cx, cy
    CameraDistortionCoeffsTensor dist_coeffs,
    at::Tensor in_image  // [B, H, W, C]
);


at::Tensor warp_image_wide_to_pinhole_tensor(
    std::string camera_model,
    at::Tensor intrins,  // fx, fy, cx, cy
    CameraDistortionCoeffsTensor dist_coeffs,
    at::Tensor wide_image,  // [B, H, W, C]
    at::Tensor axes,  // [K, 3, 3]
    int out_w, int out_h
);


at::Tensor warp_image_equirectangular_to_pinhole_tensor(
    at::Tensor equirectangular_image,  // [B, H, W, C]
    at::Tensor axes,  // [K, 3, 3]
    int out_w, int out_h
);


at::Tensor warp_image_pinhole_to_wide_tensor(
    std::string camera_model,
    at::Tensor intrins,  // fx, fy, cx, cy
    CameraDistortionCoeffsTensor dist_coeffs,
    at::Tensor pinhole_images,  // [B, K, H, W, C]
    at::Tensor axes,  // [K, 3, 3]
    int out_w, int out_h
);


at::Tensor warp_linear_depth_pinhole_to_wide_tensor(
    std::string camera_model,
    at::Tensor intrins,  // fx, fy, cx, cy
    CameraDistortionCoeffsTensor dist_coeffs,
    at::Tensor pinhole_images,  // [B, K, H, W, C]
    at::Tensor axes,  // [K, 3, 3]
    int out_w, int out_h
);


at::Tensor warp_ray_depth_pinhole_to_wide_tensor(
    std::string camera_model,
    at::Tensor intrins,  // fx, fy, cx, cy
    CameraDistortionCoeffsTensor dist_coeffs,
    at::Tensor pinhole_images,  // [B, K, H, W, C]
    at::Tensor axes,  // [K, 3, 3]
    int out_w, int out_h
);


at::Tensor warp_points_pinhole_to_wide_tensor(
    std::string camera_model,
    at::Tensor intrins,  // fx, fy, cx, cy
    CameraDistortionCoeffsTensor dist_coeffs,
    at::Tensor pinhole_images,  // [B, K, H, W, C]
    at::Tensor axes,  // [K, 3, 3]
    int out_w, int out_h
);


at::Tensor warp_depth_pinhole_to_wide_scale_matrix_tensor(
    std::string camera_model,
    at::Tensor intrins,  // fx, fy, cx, cy
    CameraDistortionCoeffsTensor dist_coeffs,
    at::Tensor pinhole_images,  // [B, K, H, W, C]
    at::Tensor axes,  // [K, 3, 3]
    int out_w, int out_h
);


at::Tensor ppisp_forward_tensor(
    at::Tensor &in_image,  // [B, H, W, C]
    at::Tensor &ppisp_params,  // [B, PPISP_NUM_PARAMS]
    at::Tensor &intrins,  // [B, 4]
    const float actual_image_width,
    const float actual_image_height,
    std::string param_type
);


std::tuple<at::Tensor, at::Tensor> ppisp_backward_tensor(
    at::Tensor &in_image,  // [B, H, W, C]
    at::Tensor &ppisp_params,  // [B, PPISP_NUM_PARAMS]
    at::Tensor &intrins,  // [B, 4]
    const float actual_image_width,
    const float actual_image_height,
    at::Tensor &v_out_image,  // [B, H, W, C]
    std::string param_type
);


std::tuple<at::Tensor, at::Tensor>
compute_ppsip_regularization_forward_tensor(
    at::Tensor &ppisp_params,  // [B, PPISP_NUM_PARAMS]
    const std::array<float, (int)PPISPRegLossIndex::length> loss_weights_0,
    std::string param_type
);


at::Tensor compute_ppsip_regularization_backward_tensor(
    at::Tensor &ppisp_params,  // [B, PPISP_NUM_PARAMS]
    const std::array<float, (int)PPISPRegLossIndex::length> loss_weights_0,
    at::Tensor &raw_losses,  // [B+1, RawPPISPRegLossIndex::length]
    at::Tensor &v_losses,  // [PPISPRegLossIndex::length]
    std::string param_type
);
