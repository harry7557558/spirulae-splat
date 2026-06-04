#pragma once

#include <cuda_runtime.h>

#include <Tensor.h>

#include <gsplat/Common.h>

#include "Primitive3DGS.cuh"

#include "common.cuh"
#include "types.cuh"



#if 0

template<ssplat::CameraModelType camera_model>
struct TileBuffers {
    long size;
    float width, height;
    const glm::mat4* __restrict__ viewmats;
    const float4* __restrict__ intrins;  // fx, fy, cx, cy
    CameraDistortionCoeffsBuffer dist_coeffs;

    TileBuffers(
        unsigned width,
        unsigned height,
        TorchTensorView viewmats,  // [B, 4, 4]
        TorchTensorView intrins,  // [B, 4], fx, fy, cx, cy
        const TorchTensorView& dist_coeffs
    ) : width((float)width), height((float)height), dist_coeffs(dist_coeffs) {
        // TODO: check dimension and shape

        static_assert(sizeof(glm::mat4) == 16*sizeof(float));
        static_assert(sizeof(glm::mat3) == 9*sizeof(float));

        this->size = std::get<2>(viewmats)[0];
        this->viewmats = (const glm::mat4*)std::get<0>(viewmats);
        this->intrins = (const float4*)std::get<0>(intrins);
    }
};


template<typename Primitive, ssplat::CameraModelType camera_model>
struct SplatTileIntersector {

    typename Primitive::WorldBuffer splats;
    long numSplats;
    float rel_scale;
    TileBuffers<camera_model> tiles;

    SplatTileIntersector(
        const typename Primitive::WorldBuffer &splats,
        const TileBuffers<camera_model> &tiles,
        float rel_scale
    );

    std::tuple<at::Tensor, at::Tensor> getIntersections_brute();

    std::tuple<at::Tensor, at::Tensor> getIntersections_lbvh();

};


std::tuple<at::Tensor, at::Tensor>
intersect_splat_tile_3dgs(
    std::vector<DeviceTensorFloatND> splats_tuple,
    unsigned width,
    unsigned height,
    TorchTensorView viewmats,
    TorchTensorView intrins,
    const std::string& camera_model,
    const TorchTensorView& dist_coeffs,
    float rel_scale
);

#endif
