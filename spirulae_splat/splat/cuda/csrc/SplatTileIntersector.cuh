#pragma once

#include <cuda_runtime.h>

#include <ATen/Tensor.h>

#include <gsplat/Common.h>

#include "Primitive3DGS.cuh"
#include "PrimitiveOpaqueTriangle.cuh"
#include "PrimitiveVoxel.cuh"

#include "common.cuh"
#include "types.cuh"



template<gsplat::CameraModelType camera_model>
struct TileBuffers {
    long size;
    float width, height;
    const glm::mat4* __restrict__ viewmats;
    const float4* __restrict__ intrins;  // fx, fy, cx, cy
    CameraDistortionCoeffsBuffer dist_coeffs;

    TileBuffers(
        unsigned width,
        unsigned height,
        const at::Tensor& viewmats,  // [B, 4, 4]
        const at::Tensor& intrins,  // [B, 4], fx, fy, cx, cy
        const CameraDistortionCoeffsTensor& dist_coeffs
    ) : width((float)width), height((float)height), dist_coeffs(dist_coeffs) {
        DEVICE_GUARD(viewmats);
        CHECK_INPUT(viewmats);
        CHECK_INPUT(intrins);
        // TODO: check dimension and shape

        static_assert(sizeof(glm::mat4) == 16*sizeof(float));
        static_assert(sizeof(glm::mat3) == 9*sizeof(float));

        this->size = viewmats.size(0);
        this->viewmats = (glm::mat4*)viewmats.data_ptr<float>();
        this->intrins = (float4*)intrins.data_ptr<float>();
    }
};


template<typename Primitive, gsplat::CameraModelType camera_model>
struct SplatTileIntersector {

    c10::TensorOptions tensorF32;
    c10::TensorOptions tensorI32, tensorI16, tensorI64, tensorU8;
    typename Primitive::World::Buffer splats;
    long numSplats;
    float rel_scale;
    TileBuffers<camera_model> tiles;

    SplatTileIntersector(
        const typename Primitive::World::Tensor &splats,
        const TileBuffers<camera_model> &tiles,
        float rel_scale
    );

    std::tuple<at::Tensor, at::Tensor> getIntersections_brute();

    std::tuple<at::Tensor, at::Tensor> getIntersections_lbvh();

};


std::tuple<at::Tensor, at::Tensor>
intersect_splat_tile_3dgs(
    Vanilla3DGS::World::TensorTuple splats_tuple,
    unsigned width,
    unsigned height,
    const at::Tensor& viewmats,
    const at::Tensor& intrins,
    const gsplat::CameraModelType& camera_model,
    const CameraDistortionCoeffsTensor& dist_coeffs,
    float rel_scale
);

std::tuple<at::Tensor, at::Tensor>
intersect_splat_tile_opaque_triangle(
    OpaqueTriangle::World::TensorTuple splats_tuple,
    unsigned width,
    unsigned height,
    const at::Tensor& viewmats,
    const at::Tensor& intrins,
    const gsplat::CameraModelType& camera_model,
    const CameraDistortionCoeffsTensor& dist_coeffs,
    float rel_scale
);

std::tuple<at::Tensor, at::Tensor>
intersect_splat_tile_voxel(
    VoxelPrimitive::World::TensorTuple splats_tuple,
    unsigned width,
    unsigned height,
    const at::Tensor& viewmats,
    const at::Tensor& intrins,
    const gsplat::CameraModelType& camera_model,
    const CameraDistortionCoeffsTensor& dist_coeffs,
    float rel_scale
);
