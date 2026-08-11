#pragma once

#include "backend/api/BackendTypes.h"

#include <core/Tensor.h>

#include <core/Common.cuh>

#include "primitives/Primitive3DGS.cuh"


template<CameraModelType camera_model, CameraDistortionType distortion>
struct TileBuffers {
    long size;
    float width, height;
    // Row-major 4x4 view matrices, stride 16 floats per tile.
    const float* __restrict__ viewmats;
    const float4* __restrict__ intrins;  // fx, fy, cx, cy
    CameraDistortionCoeffsBuffer dist_coeffs;

    TileBuffers(
        unsigned width,
        unsigned height,
        const float* viewmats_ptr,   // [B*16], row-major 4x4
        const float4* intrins_ptr,   // [B], fx, fy, cx, cy
        float* dist_coeffs_ptr,      // [B*8] or nullptr
        long size
    )
        : size(size),
          width((float)width), height((float)height),
          viewmats(viewmats_ptr), intrins(intrins_ptr),
          dist_coeffs(dist_coeffs_ptr) {}
};


template<typename Primitive, CameraModelType camera_model, CameraDistortionType distortion>
struct SplatTileIntersector {

    typename Primitive::WorldBuffer splats;
    long numSplats;
    float rel_scale;
    TileBuffers<camera_model, distortion> tiles;

    SplatTileIntersector(
        const typename Primitive::WorldBuffer &splats,
        const TileBuffers<camera_model, distortion> &tiles,
        float rel_scale
    );

    std::tuple<DeviceVector<int32_t>, DeviceVector<int32_t>> getIntersections_brute();

    std::tuple<DeviceVector<int32_t>, DeviceVector<int32_t>> getIntersections_lbvh();

};


std::tuple<DeviceVector<int32_t>, DeviceVector<int32_t>>
intersect_splat_tile_3dgs(
    std::vector<DeviceTensorFloatND> splats_tuple,
    unsigned width,
    unsigned height,
    DeviceVector<float> viewmats,       // [C*16], row-major 4x4 per camera
    DeviceVector<float4> intrins,       // [C]
    const std::string& camera_model,
    const std::string& distortion,
    const DeviceTensor2D<float>& dist_coeffs,  // [C, 8], may be null
    float rel_scale
);
