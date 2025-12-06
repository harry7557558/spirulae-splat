#pragma once

#include <cuda_runtime.h>
#include <torch/types.h>

#include <gsplat/Common.h>

#include "PrimitiveVoxel.cuh"

#include "common.cuh"
#include "types.cuh"


struct SVHashTensor {

    at::Tensor vertList;  // [numVerts, 3], int32
    at::Tensor cellList;  // [numCells, 4], int32

    at::Tensor vertIndexMap;  // [#, 2], uint32, hash table
    at::Tensor cellIndexMap;  // [#, 2], uint32, hash table

    float3 origin;
    float scale;

    inline long numVerts() const { return vertList.size(-2); }
    inline long numCells() const { return cellList.size(-2); }

    typedef std::tuple<
        at::Tensor, at::Tensor,
        at::Tensor, at::Tensor,
        std::tuple<float, float, float, float>
    > TensorTuple;

    TensorTuple tuple() {
        return std::make_tuple(
            vertList, cellList,
            vertIndexMap, cellIndexMap,
            std::make_tuple(origin.x, origin.y, origin.z, scale)
        );
    }

    SVHashTensor() {};

    SVHashTensor(TensorTuple& tp) {
        vertList = std::get<0>(tp);
        cellList = std::get<1>(tp);
        vertIndexMap = std::get<2>(tp);
        cellIndexMap = std::get<3>(tp);
        DEVICE_GUARD(vertList);
        CHECK_INPUT(vertList);
        CHECK_INPUT(cellList);
        CHECK_INPUT(vertIndexMap);
        CHECK_INPUT(cellIndexMap);
        auto cs = std::get<4>(tp);
        origin = { std::get<0>(cs), std::get<1>(cs), std::get<2>(cs) };
        scale = std::get<3>(cs);
    };

    SVHashTensor(long num_verts, long num_cells);

};


SVHashTensor::TensorTuple svhashCreateInitialVolume(
    std::tuple<float, float, float> p0_,
    std::tuple<float, float, float> p1_,
    std::tuple<unsigned, unsigned, unsigned> res_
);

std::tuple<at::Tensor, at::Tensor>
svhashGetVoxels(SVHashTensor::TensorTuple tensorTuples);

std::tuple<SVHashTensor::TensorTuple, at::Tensor, at::Tensor, at::Tensor>
svhashSplitVoxels(
    SVHashTensor::TensorTuple tensorTuples,
    at::Tensor split_mask
);
