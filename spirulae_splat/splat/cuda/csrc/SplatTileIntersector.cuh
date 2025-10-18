#pragma once

#include "common.cuh"

struct SplatBuffers {
    long size;
    const glm::vec3* __restrict__ means;
    const glm::vec3* __restrict__ scales;
    const float* __restrict__ opacs;
    const float4* __restrict__ quats;

    SplatBuffers(
        torch::Tensor& means,
        torch::Tensor& scales,
        torch::Tensor& opacs,
        torch::Tensor& quats
    ) {
        DEVICE_GUARD(means);
        CHECK_INPUT(means);
        CHECK_INPUT(scales);
        CHECK_INPUT(opacs);
        CHECK_INPUT(quats);
        // TODO: check dimension and shape

        this->size = means.size(0);
        this->means = (glm::vec3*)means.data_ptr<float>();
        this->scales = (glm::vec3*)scales.data_ptr<float>();
        this->opacs = opacs.data_ptr<float>();
        this->quats = (float4*)quats.data_ptr<float>();
    }
};



struct TileBuffers {
    long size;
    float width, height;
    const glm::mat4* __restrict__ viewmats;
    const glm::mat3* __restrict__ Ks;  // TODO: make aligned

    TileBuffers(
        unsigned width,
        unsigned height,
        torch::Tensor& viewmats,  // [B, 4, 4]
        torch::Tensor& Ks  // [B, 3, 3]
    ) : width((float)width), height((float)height) {
        DEVICE_GUARD(viewmats);
        CHECK_INPUT(viewmats);
        CHECK_INPUT(Ks);
        // TODO: check dimension and shape

        static_assert(sizeof(glm::mat4) == 16*sizeof(float));
        static_assert(sizeof(glm::mat3) == 9*sizeof(float));

        this->size = viewmats.size(0);
        this->viewmats = (glm::mat4*)viewmats.data_ptr<float>();
        this->Ks = (glm::mat3*)Ks.data_ptr<float>();
    }
};


struct SplatTileIntersector {

    float3 rootAABBMin, rootAABBMax;

    c10::TensorOptions tensorF32;
    c10::TensorOptions tensorI32, tensorI16, tensorI64, tensorU8;
    SplatBuffers splats;
    TileBuffers tiles;

    SplatTileIntersector(
        c10::TensorOptions tensorOptions,
        const SplatBuffers &splats,
        const TileBuffers &tiles
    );

    std::tuple<torch::Tensor, torch::Tensor> getIntersections_brute();

    std::tuple<torch::Tensor, torch::Tensor> getIntersections_octree();

    std::tuple<torch::Tensor, torch::Tensor> getIntersections_lbvh();

    static std::tuple<torch::Tensor, torch::Tensor>
    intersect_splat_tile(
        torch::Tensor& means,
        torch::Tensor& scales,
        torch::Tensor& opacs,
        torch::Tensor& quats,
        unsigned width,
        unsigned height,
        torch::Tensor& viewmats,
        torch::Tensor& Ks
    ) {
        SplatBuffers splat_buffers = {means, scales, opacs, quats};
        TileBuffers tile_buffers = {width, height, viewmats, Ks};

        return SplatTileIntersector(
            means.options(),
            splat_buffers,
            tile_buffers
        ).getIntersections_lbvh();
    }

};

