#pragma once

#include <cuda_runtime.h>
#include <cstdint>

#include <Tensor.h>

#include <gsplat/Common.h>

#include "Primitive3DGS.cuh"
#include "Primitive3DGUT.cuh"

#include "common.cuh"
#include "types.cuh"


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



std::tuple<
    DeviceVector<int32_t>, DeviceVector<int32_t>, DeviceVector<float4>,
    DeviceVector<float>, std::vector<DeviceTensorFloatND>
> projection_3dgs_packed_forward(
    const int64_t num_splats,
    const int max_sh_degree,
    const std::vector<DeviceTensorFloatND> &in_splats,
    TorchTensorView viewmats,
    TorchTensorView intrins,
    const uint32_t image_width,
    const uint32_t image_height,
    const std::string camera_model,
    const TorchTensorView dist_coeffs,
    DeviceVector<float> radii
);


std::tuple<
    DeviceVector<int32_t>, DeviceVector<int32_t>, DeviceVector<float4>,
    DeviceVector<float>, std::vector<DeviceTensorFloatND>
> projection_mip_packed_forward(
    const int64_t num_splats,
    const int max_sh_degree,
    const std::vector<DeviceTensorFloatND> &in_splats,
    TorchTensorView viewmats,
    TorchTensorView intrins,
    const uint32_t image_width,
    const uint32_t image_height,
    const std::string camera_model,
    const TorchTensorView dist_coeffs,
    DeviceVector<float> radii
);


std::tuple<
    DeviceVector<int32_t>, DeviceVector<int32_t>, DeviceVector<float4>,
    DeviceVector<float>, std::vector<DeviceTensorFloatND>
> projection_3dgut_packed_forward(
    const int64_t num_splats,
    const int max_sh_degree,
    const std::vector<DeviceTensorFloatND> &in_splats,
    TorchTensorView viewmats,
    TorchTensorView intrins,
    const uint32_t image_width,
    const uint32_t image_height,
    const std::string camera_model,
    const TorchTensorView dist_coeffs,
    DeviceVector<float> radii
);
