#pragma once

#include <torch/types.h>
#include "glm/glm/glm.hpp"

#define CHECK_CUDA(x) TORCH_CHECK(x.is_cuda(), #x " must be a CUDA tensor")
#define CHECK_CONTIGUOUS(x)                                                    \
    TORCH_CHECK(x.is_contiguous(), #x " must be contiguous")
#define CHECK_INPUT(x)                                                         \
    CHECK_CUDA(x);                                                             \
    CHECK_CONTIGUOUS(x)
#define DEVICE_GUARD(_ten) \
    const at::cuda::OptionalCUDAGuard device_guard(device_of(_ten));

#define CHECK_DEVICE_ERROR(call)                                      \
do {                                                                \
    cudaError_t err = call;                                         \
    if (err != cudaSuccess) {                                       \
        fprintf(stderr, "CUDA Error at %s:%d: %s\n",                \
                __FILE__, __LINE__, cudaGetErrorString(err));       \
        exit(EXIT_FAILURE);                                         \
    }                                                               \
} while (0)


struct SplatBuffers {
    long size;
    glm::vec3* __restrict__ means;
    glm::vec3* __restrict__ scales;
    float* __restrict__ opacs;
    float4* __restrict__ quats;

    SplatBuffers(
        torch::Tensor& means,
        torch::Tensor& scales,
        torch::Tensor& opacs,
        torch::Tensor& quats
    ) {
        CHECK_CUDA(means);
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
    glm::mat4* __restrict__ viewmats;
    glm::mat3* __restrict__ Ks;  // TODO: make aligned

    TileBuffers(
        torch::Tensor& viewmats,  // [B, 4, 4]
        torch::Tensor& Ks  // [B, 3, 3]
    ) {
        CHECK_CUDA(viewmats);
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
        TileBuffers tiles
    );

    std::tuple<torch::Tensor, torch::Tensor> getIntersections_brute();

    std::tuple<torch::Tensor, torch::Tensor> getIntersections_octree();

    std::tuple<torch::Tensor, torch::Tensor> getIntersections_lbvh();

};
