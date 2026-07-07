#pragma once

#include <cuda_runtime.h>
#include <cstdint>

#include <Tensor.h>

#include <Common.cuh>

#include "Primitive3DGS.cuh"
#include "Primitive3DGUT.cuh"

#include "NonShQuantState.h"


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



void fused_projection_bwd_optimizer_3dgs(
    const int64_t num_splats,
    const int max_sh_degree,
    std::vector<DeviceTensorFloatND> splats_world,
    TorchTensorView viewmats,
    TorchTensorView intrins,
    const uint32_t image_width,
    const uint32_t image_height,
    const std::string camera_model,
    const TorchTensorView dist_coeffs,
    const DeviceVector<int32_t> camera_ids,
    const DeviceVector<int32_t> gaussian_ids,
    DeviceTensorFloatND aabb,
    const std::vector<DeviceTensorFloatND> v_splats_world,
    const std::vector<DeviceTensorFloatND> v_splats_screen,
    const std::vector<DeviceTensorFloatND> g1_splats_world,
    const std::vector<DeviceTensorFloatND> g2_splats_world,
    const std::optional<TorchTensorView> sh_packed,
    const std::optional<TorchTensorView> sh_quant_bounds,
    const std::optional<TorchTensorView> sh_value_packed,
    const std::optional<TorchTensorView> sh_value_bounds,
    NonShQuantState non_sh,
    DeviceVector<float> radii,
    DeviceVector<float> densify_score,
    const float lr_means,
    const float lr_quats,
    const float lr_scales,
    const float lr_opacs,
    const float lr_features_dc,
    const float lr_features_sh,
    const float max_gauss_ratio,
    const float scale_regularization_weight,
    const float mcmc_opacity_reg_weight,
    const float mcmc_scale_reg_weight,
    const float erank_reg_weight,
    const float erank_reg_weight_s3,
    const float quat_norm_reg_weight,
    const float sh_reg_weight,
    bool use_scale_agnostic_mean,
    bool color_trust_linear,
    float eps_tr,
    std::variant<int32_t, TorchTensorView> step,
    int quantization_level
);


void fused_projection_bwd_optimizer_mip(
    const int64_t num_splats,
    const int max_sh_degree,
    std::vector<DeviceTensorFloatND> splats_world,
    TorchTensorView viewmats,
    TorchTensorView intrins,
    const uint32_t image_width,
    const uint32_t image_height,
    const std::string camera_model,
    const TorchTensorView dist_coeffs,
    const DeviceVector<int32_t> camera_ids,
    const DeviceVector<int32_t> gaussian_ids,
    DeviceTensorFloatND aabb,
    const std::vector<DeviceTensorFloatND> v_splats_world,
    const std::vector<DeviceTensorFloatND> v_splats_screen,
    const std::vector<DeviceTensorFloatND> g1_splats_world,
    const std::vector<DeviceTensorFloatND> g2_splats_world,
    const std::optional<TorchTensorView> sh_packed,
    const std::optional<TorchTensorView> sh_quant_bounds,
    const std::optional<TorchTensorView> sh_value_packed,
    const std::optional<TorchTensorView> sh_value_bounds,
    NonShQuantState non_sh,
    DeviceVector<float> radii,
    DeviceVector<float> densify_score,
    const float lr_means,
    const float lr_quats,
    const float lr_scales,
    const float lr_opacs,
    const float lr_features_dc,
    const float lr_features_sh,
    const float max_gauss_ratio,
    const float scale_regularization_weight,
    const float mcmc_opacity_reg_weight,
    const float mcmc_scale_reg_weight,
    const float erank_reg_weight,
    const float erank_reg_weight_s3,
    const float quat_norm_reg_weight,
    const float sh_reg_weight,
    bool use_scale_agnostic_mean,
    bool color_trust_linear,
    float eps_tr,
    std::variant<int32_t, TorchTensorView> step,
    int quantization_level
);


void fused_projection_bwd_optimizer_3dgut(
    const int64_t num_splats,
    const int max_sh_degree,
    std::vector<DeviceTensorFloatND> splats_world,
    TorchTensorView viewmats,
    TorchTensorView intrins,
    const uint32_t image_width,
    const uint32_t image_height,
    const std::string camera_model,
    const TorchTensorView dist_coeffs,
    const DeviceVector<int32_t> camera_ids,
    const DeviceVector<int32_t> gaussian_ids,
    DeviceTensorFloatND aabb,
    const std::vector<DeviceTensorFloatND> v_splats_world,
    const std::vector<DeviceTensorFloatND> v_splats_screen,
    const std::vector<DeviceTensorFloatND> g1_splats_world,
    const std::vector<DeviceTensorFloatND> g2_splats_world,
    const std::optional<TorchTensorView> sh_packed,
    const std::optional<TorchTensorView> sh_quant_bounds,
    const std::optional<TorchTensorView> sh_value_packed,
    const std::optional<TorchTensorView> sh_value_bounds,
    NonShQuantState non_sh,
    DeviceVector<float> radii,
    DeviceVector<float> densify_score,
    const float lr_means,
    const float lr_quats,
    const float lr_scales,
    const float lr_opacs,
    const float lr_features_dc,
    const float lr_features_sh,
    const float max_gauss_ratio,
    const float scale_regularization_weight,
    const float mcmc_opacity_reg_weight,
    const float mcmc_scale_reg_weight,
    const float erank_reg_weight,
    const float erank_reg_weight_s3,
    const float quat_norm_reg_weight,
    const float sh_reg_weight,
    bool use_scale_agnostic_mean,
    bool color_trust_linear,
    float eps_tr,
    std::variant<int32_t, TorchTensorView> step,
    int quantization_level
);
