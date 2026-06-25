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
    std::vector<DeviceTensorFloatND>, std::vector<DeviceTensorFloatND>,  // gradient
    DeviceVector<float>  // accum_weight
> rasterize_to_pixels_3dgs_bwd(
    // Gaussian parameters
    int64_t num_splats,
    std::vector<DeviceTensorFloatND> splats_w,
    std::vector<DeviceTensorFloatND> splats_s,
    DeviceVector<int32_t> gaussian_ids,
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    // intersections
    const DeviceTensor3D<int32_t> tile_offsets, // [I, tile_height, tile_width]
    const DeviceVector<int32_t> flatten_ids,    // [n_isects]
    // forward outputs
    const DeviceTensor3D<float> render_Ts,  // [I, image_height, image_width]
    const DeviceTensor3D<int32_t> last_ids, // [I, image_height, image_width]
    RenderOutput::TensorTuple render_outputs_tuple,
    std::optional<RenderOutput::TensorTuple> distortion_fwd_outputs,  // forward D
    DistortionType dist_type,  // distortion channel set (None/D/RGB_D)
    DeviceTensor3D<float> accum_weight_map,  // [I, H, W]
    // gradients of outputs
    RenderOutput::TensorTuple v_render_outputs,
    const DeviceTensor3D<float> v_render_Ts,
    const DeviceTensor3D<float> v_median,  // [I, H, W], optional
    std::optional<RenderOutput::TensorTuple> v_distortion_outputs,
    std::optional<std::vector<DeviceTensorFloatND>> v_splats_w,
    std::optional<std::vector<DeviceTensorFloatND>> v_splats_s
);
