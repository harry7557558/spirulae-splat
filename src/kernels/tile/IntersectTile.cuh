#pragma once

#include <cstdint>

#include <core/Tensor.h>

#include "primitives/Primitive3DGS.cuh"
#include "primitives/Primitive3DGUT.cuh"


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



int64_t intersect_tile_count(int width, int height);


void compute_tile_active(
    TorchTensorView mask,   // [I, H_mask, W_mask] bool
    int I, int width, int height,
    int32_t* tile_active    // [I, tile_h, tile_w]
);


std::tuple<
    DeviceVector<int64_t>,    // isect_ids [n_isects]
    DeviceVector<int32_t>,    // flatten_ids [n_isects]
    DeviceTensor3D<int32_t>   // offsets [I, tile_h, tile_w]
> do_intersect_tile_generic(
    DeviceTensorFloatND aabb,     // [*N, 4] float32
    DeviceTensorFloatND depths,   // [*N] float32
    DeviceTensorFloatND* proj_xy,    // null for AABB mode
    DeviceTensorFloatND* proj_conic, // non-null enables ellipse mode
    DeviceTensorFloatND* proj_opac,  // null for AABB mode
    const uint32_t I,
    TorchTensorView intrins,      // [I, 4]
    const uint32_t image_width,
    const uint32_t image_height,
    DeviceVector<int32_t>* image_ids, // null for non-packed
    const int32_t* tile_active        // [I, tile_h, tile_w]; null = all live
);
