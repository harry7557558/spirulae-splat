#pragma once

#include <cuda_runtime.h>
#include <cstdint>

#include <Tensor.h>

#include <Common.cuh>

#include "Primitive3DGS.cuh"
#include "Primitive3DGUT.cuh"


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



std::tuple<
    DeviceTensor2D<float4>, DeviceTensor2D<float>, std::vector<DeviceTensorFloatND>
> projection_3dgs_forward(
    const int64_t num_splats, const int max_sh_degree,
    const std::vector<DeviceTensorFloatND> &in_splats,
    TorchTensorView viewmats, TorchTensorView intrins,
    const uint32_t image_width, const uint32_t image_height,
    const std::string camera_model, const TorchTensorView dist_coeffs,
    DeviceVector<float> radii,
    const std::optional<TorchTensorView> sh_value_packed,
    const std::optional<TorchTensorView> sh_value_bounds,
    const uint32_t num_sh_buffer,
    const int sh_value_bits,
    // SH value-bounds cell stride. 0 = FPBO per-splat-block layout
    // (256 * 3 * num_sh_buffer cells/bound). 256 = non-FPBO per-cell-block
    // layout. Plumbed unchanged from EngineForward.cpp.
    const int64_t sh_bounds_stride
);


std::tuple<
    DeviceTensor2D<float4>, DeviceTensor2D<float>, std::vector<DeviceTensorFloatND>
> projection_mip_forward(
    const int64_t num_splats, const int max_sh_degree,
    const std::vector<DeviceTensorFloatND> &in_splats,
    TorchTensorView viewmats, TorchTensorView intrins,
    const uint32_t image_width, const uint32_t image_height,
    const std::string camera_model, const TorchTensorView dist_coeffs,
    DeviceVector<float> radii,
    const std::optional<TorchTensorView> sh_value_packed,
    const std::optional<TorchTensorView> sh_value_bounds,
    const uint32_t num_sh_buffer,
    const int sh_value_bits,
    // SH value-bounds cell stride. 0 = FPBO per-splat-block layout
    // (256 * 3 * num_sh_buffer cells/bound). 256 = non-FPBO per-cell-block
    // layout. Plumbed unchanged from EngineForward.cpp.
    const int64_t sh_bounds_stride
);


std::tuple<
    DeviceTensor2D<float4>, DeviceTensor2D<float>, std::vector<DeviceTensorFloatND>
> projection_3dgut_forward(
    const int64_t num_splats, const int max_sh_degree,
    const std::vector<DeviceTensorFloatND> &in_splats,
    TorchTensorView viewmats, TorchTensorView intrins,
    const uint32_t image_width, const uint32_t image_height,
    const std::string camera_model, const TorchTensorView dist_coeffs,
    DeviceVector<float> radii,
    const std::optional<TorchTensorView> sh_value_packed,
    const std::optional<TorchTensorView> sh_value_bounds,
    const uint32_t num_sh_buffer,
    const int sh_value_bits,
    // SH value-bounds cell stride. 0 = FPBO per-splat-block layout
    // (256 * 3 * num_sh_buffer cells/bound). 256 = non-FPBO per-cell-block
    // layout. Plumbed unchanged from EngineForward.cpp.
    const int64_t sh_bounds_stride
);
