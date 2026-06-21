#pragma once

#include <cuda.h>
#include <cuda_runtime.h>
#include <cstdint>

#include <Tensor.h>

#include "Primitive3DGS.cuh"
// #include "Primitive3DGUT.cuh"

#include "Common.cuh"


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



std::tuple<
    RenderOutput::TensorTuple,  // renders
    DeviceTensor3D<float>,  // transmittances
    DeviceTensor3D<int32_t>,  // last_ids
    RenderOutput::TensorTuple,  // renders2, optional
    RenderOutput::TensorTuple,  // distortions, optional
    DeviceTensor3D<float>,  // median depth, optional
    DeviceTensor3D<float2>  // median near-anchor (z1, v1), optional
> rasterize_to_pixels_3dgs_fwd(
    // Gaussian parameters
    int64_t num_splats,
    std::vector<DeviceTensorFloatND> splats_w,
    std::vector<DeviceTensorFloatND> splats_s,
    DeviceVector<int32_t> gaussian_ids,
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    // intersections
    const DeviceTensor3D<int32_t> tile_offsets,
    const DeviceVector<int32_t> flatten_ids,
    bool output_distortion,
    bool output_median
);


std::tuple<
    RenderOutput::TensorTuple,  // renders
    DeviceTensor3D<float>,  // transmittances
    DeviceTensor3D<int32_t>,  // last_ids
    RenderOutput::TensorTuple,  // renders2, optional
    RenderOutput::TensorTuple,  // distortions, optional
    DeviceTensor3D<float>,  // median depth, optional
    DeviceTensor3D<float2>  // median near-anchor (z1, v1), optional
> rasterize_to_pixels_mip_fwd(
    // Gaussian parameters
    int64_t num_splats,
    std::vector<DeviceTensorFloatND> splats_w,
    std::vector<DeviceTensorFloatND> splats_s,
    DeviceVector<int32_t> gaussian_ids,
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    // intersections
    const DeviceTensor3D<int32_t> tile_offsets,
    const DeviceVector<int32_t> flatten_ids,
    bool output_distortion,
    bool output_median
);
