#pragma once

#include <cuda_runtime.h>
#include <cstdint>

#include <Tensor.h>

#include <gsplat/Common.h>

#include "Primitive3DGS.cuh"
#include "Primitive3DGUT.cuh"
#include "Primitive3DGUT_SV.cuh"
#include "PrimitiveOpaqueTriangle.cuh"
#include "PrimitiveVoxel.cuh"

#include "types.cuh"
#include "common.cuh"


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



void fused_projection_bwd_optimizer_3dgut(
    // fwd inputs
    const int64_t num_splats,
    std::vector<DeviceTensorFloatND> splats_world,
    TorchTensorView viewmats,
    TorchTensorView intrins,
    const uint32_t image_width,
    const uint32_t image_height,
    const std::string camera_model,
    const TorchTensorView dist_coeffs,
    // fwd outputs
    const DeviceVector<int32_t> camera_ids,
    const DeviceVector<int32_t> gaussian_ids,
    DeviceTensorFloatND aabb,
    // grad outputs
    const std::vector<DeviceTensorFloatND> v_splats_world,
    const std::optional<std::vector<DeviceTensorFloatND>> vr_splats_world,
    const std::optional<std::vector<DeviceTensorFloatND>> h_splats_world,
    const std::vector<DeviceTensorFloatND> v_splats_screen,
    const std::optional<std::vector<DeviceTensorFloatND>> vr_splats_screen,
    const std::optional<std::vector<DeviceTensorFloatND>> h_splats_screen,
    // optimizer states
    const std::vector<DeviceTensorFloatND> g1_splats_world,
    const std::vector<DeviceTensorFloatND> g2_splats_world,
    const std::optional<TorchTensorView> sh_packed,         // AoS packed SH state
    const std::optional<TorchTensorView> sh_quant_bounds,
    // optimizer params
    DeviceVector<float> radii,
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
    std::variant<int32_t, TorchTensorView> step
);
