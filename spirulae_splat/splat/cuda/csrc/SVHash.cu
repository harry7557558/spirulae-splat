#include "SVHash.cuh"


// hash functions - https://github.com/aappleby/smhasher/wiki/MurmurHash3

__forceinline__ __device__ uint32_t _hash(uint32_t k) {
    k ^= k >> 16;
    k *= 0x85ebca6b;
    k ^= k >> 13;
    k *= 0xc2b2ae35;
    k ^= k >> 16;
    return k;
}

__forceinline__ __device__ uint64_t _hash(uint64_t k) {
    k ^= k >> 33;
    k *= 0xff51afd7ed558ccd;
    k ^= k >> 33;
    k *= 0xc4ceb9fe1a85ec53;
    k ^= k >> 33;
    return k;
}

__forceinline__ __device__ uint32_t _probe(uint32_t k) {
    return k + 1;
}

__forceinline__ __device__ uint64_t _probe(uint64_t k) {
    return k + 1;
}


// device hash table

template<typename Value>
struct _SVHashTable {

    static constexpr uint64_t kEmpty = ~(uint64_t)0;

    struct KeyValue {
        uint64_t key;
        Value value;
    };
    #ifdef __CUDACC__
    static_assert(sizeof(KeyValue) == sizeof(uint64_t) + sizeof(Value));  // packed
    static_assert(sizeof(KeyValue) % sizeof(uint64_t) == 0);  // required for memory aligned atomic
    #endif

    int64_t capacity;
    KeyValue* __restrict__ data;

    _SVHashTable() : capacity(0), data(nullptr) {}

    _SVHashTable(at::Tensor &tensor) {
        DEVICE_GUARD(tensor);
        CHECK_INPUT(tensor);
        capacity = tensor.numel() / (sizeof(KeyValue) / sizeof(int64_t));
        data = (KeyValue*)tensor.data_ptr<int64_t>();
    }

    __device__ bool lookup(uint64_t key, Value &out_value) const {
        uint64_t slot = _hash(key) % capacity;
        while (true) {
            if (data[slot].key == key) {
                out_value = data[slot].value;
                return true;
            }
            if (data[slot].key == kEmpty)
                return false;
            slot = _probe(slot) % capacity;
        }
    }

    template<bool overwrite>
    __device__ bool insert(uint64_t key, Value value) {
        uint64_t slot = _hash(key) % capacity;
        while (true) {
            uint64_t prev = atomicCAS((uint64_t*)&data[slot].key, (uint64_t)kEmpty, (uint64_t)key);
            if (prev == kEmpty || (overwrite && prev == key)) {
                data[slot].value = value;
                return true;
            } else if (prev == key) {
                return false;
            }
            slot = _probe(slot) % capacity;
        }
    }

};


// device hash grid

struct SVHash {

    static constexpr int kNumResolutions = 18;  // inclusive
    static constexpr int kUnitResolution = 9;
    struct Cell {
        int3 cell;  // xyz
        int resolution;  // 0 is unit cube, negative is larger, positive is smaller
    };

    int64_t numVerts;
    int64_t numCells;

    int3* __restrict__ vertList;
    int4* __restrict__ cellList;

    _SVHashTable<int64_t> vertIndexMap;  // map xyz to index
    _SVHashTable<int64_t> cellIndexMap;  // map cell to index

    float3 origin;
    float scale;

    SVHash() {}

    SVHash(SVHashTensor tensor) :
        vertIndexMap(tensor.vertIndexMap),
        cellIndexMap(tensor.cellIndexMap)
    {
        DEVICE_GUARD(tensor.vertList);
        CHECK_INPUT(tensor.vertList);
        CHECK_INPUT(tensor.cellList);
        this->numVerts = tensor.vertList.size(0);
        this->numCells = tensor.cellList.size(0);
        vertList = (int3*)tensor.vertList.data_ptr<int>();
        cellList = (int4*)tensor.cellList.data_ptr<int>();
        this->origin = tensor.origin;
        this->scale = tensor.scale;
    }

    static __forceinline__ __device__ uint64_t getKey(int3 vert) {
        static constexpr uint64_t mult = (uint64_t)1 << (uint64_t)(kNumResolutions + 1);
        uint64_t z = (uint64_t)vert.z & (mult-1);
        uint64_t y = (uint64_t)vert.y & (mult-1);
        uint64_t x = (uint64_t)vert.x & (mult-1);
        return (z * mult + y) * mult + x;
    }
    static __forceinline__ __device__ uint64_t getKey(int4 cell) {
        uint64_t key = getKey(make_int3(cell.x, cell.y, cell.z));
        return key | ((uint64_t)cell.w << (uint64_t)(3*(kNumResolutions+1)));
    }
};


SVHashTensor::SVHashTensor(int64_t num_verts, int64_t num_cells) {
    auto kOptI32 = at::TensorOptions().dtype(torch::kInt32).device(torch::kCUDA);
    vertList = at::empty({num_verts, 3}, kOptI32);
    cellList = at::empty({num_cells, 4}, kOptI32);
    auto kOptI64 = at::TensorOptions().dtype(torch::kInt64).device(torch::kCUDA);
    vertIndexMap = at::empty({2*num_verts, 2}, kOptI64);
    cellIndexMap = at::empty({2*num_cells, 2}, kOptI64);
    cudaMemset(vertIndexMap.data_ptr<int64_t>(), 0xff, vertIndexMap.numel()*sizeof(int64_t));
    cudaMemset(cellIndexMap.data_ptr<int64_t>(), 0xff, cellIndexMap.numel()*sizeof(int64_t));
}



template<typename dtype>
__global__ void computeIndexMap_kernel(
    int64_t num_elements,
    const dtype* __restrict__ list,
    _SVHashTable<int64_t> indexMap
) {
    uint tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_elements)
        return;

    dtype key = list[tid];
    indexMap.insert<true>(SVHash::getKey(key), tid);
}

template<typename dtype>
at::Tensor computeIndexMap(const at::Tensor& list) {
    int64_t num_elements = list.numel() / (sizeof(dtype) / sizeof(int32_t));
    int64_t hash_table_capacity = (int)ceil(num_elements * 1.5);
    at::Tensor indexMap = at::empty({hash_table_capacity, 2}, list.options().dtype(torch::kInt64));
    cudaMemset(indexMap.data_ptr<int64_t>(), 0xff, indexMap.numel()*sizeof(int64_t));

    computeIndexMap_kernel<<<_LAUNCH_ARGS_1D(num_elements, 256)>>>(
        num_elements,
        (dtype*)list.data_ptr<int32_t>(),
        _SVHashTable<int64_t>(indexMap)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return indexMap;
}




/* Initialization */

__global__ void svhashCreateInitialVolume_kernel(
    SVHash svhash,
    float3 p0, float3 p1, uint3 res,
    int level0, int level_sc
) {
    uint tid = blockIdx.x * blockDim.x + threadIdx.x;
    uint xi = tid % (res.x + 1); tid /= res.x + 1;
    uint yi = tid % (res.y + 1); tid /= res.y + 1;
    uint zi = tid;
    if (zi >= res.z + 1)
        return;

    // insert vert
    int4 key = {
        ((int)xi - (int)(res.x + 1) / 2) * level_sc,
        ((int)yi - (int)(res.y + 1) / 2) * level_sc,
        ((int)zi - (int)(res.z + 1) / 2) * level_sc,
        level0
    };
    uint vert_idx = (zi * (res.y+1) + yi) * (res.x+1) + xi;
    svhash.vertList[vert_idx] = {key.x, key.y, key.z};
    svhash.vertIndexMap.insert<true>(svhash.getKey(make_int3(key.x, key.y, key.z)), vert_idx);

    if (xi == res.x || yi == res.y || zi == res.z)
        return;

    // insert cell
    uint cell_idx = (zi * res.y + yi) * res.x + xi;
    svhash.cellList[cell_idx] = key;
    svhash.cellIndexMap.insert<true>(svhash.getKey(key), cell_idx);

}


SVHashTensor::TensorTuple svhashCreateInitialVolume(
    std::tuple<float, float, float> p0_,
    std::tuple<float, float, float> p1_,
    std::tuple<unsigned, unsigned, unsigned> res_
) {
    float3 p0 = { std::get<0>(p0_), std::get<1>(p0_), std::get<2>(p0_) };
    float3 p1 = { std::get<0>(p1_), std::get<1>(p1_), std::get<2>(p1_) };
    uint3 res = { std::get<0>(res_), std::get<1>(res_), std::get<2>(res_) };

    int64_t num_verts = (res.x+1)*(res.y+1)*(res.z+1);
    int64_t num_cells = res.x*res.y*res.z;

    SVHashTensor result(num_verts, num_cells);

    int level0 = max(SVHash::kUnitResolution + (int)round(-1.0f/6.0f * log2(res.x*res.y*res.z)), 0);
    int level_sc = 1 << level0;

    svhashCreateInitialVolume_kernel<<<_LAUNCH_ARGS_1D(num_verts, 256)>>>(
        SVHash(result), p0, p1, res, level0, level_sc
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    // TODO: exact alignment
    result.origin = { 0.5f*(p0.x+p1.x), 0.5f*(p0.y+p1.y), 0.5f*(p0.z+p1.z) };
    result.scale = fmax(fmax((p1.x-p0.x)/res.x, (p1.y-p0.y)/res.y), (p1.z-p0.z)/res.z) / (float)level_sc;

    return result.tuple();
}




/* Access elements */

__global__ void svhashGetVoxels_kernel(
    const SVHash svhash,
    float4* __restrict__ out_voxels,
    FixedArray<int, 8>* __restrict__ out_vert_indices
) {
    uint tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= svhash.numCells)
        return;

    int4 voxel = svhash.cellList[tid];
    int stride = 1 << voxel.w;
    out_voxels[tid] = {
        (float)voxel.x * svhash.scale + svhash.origin.x,
        (float)voxel.y * svhash.scale + svhash.origin.x,
        (float)voxel.z * svhash.scale + svhash.origin.x,
        (float)stride * svhash.scale
    };

    FixedArray<int, 8> vert_indices = {-1, -1, -1, -1, -1, -1, -1, -1};
    #pragma unroll
    for (int i = 0; i < 8; i++) {
        int3 vi = {
            voxel.x + stride * (i & 1),
            voxel.y + stride * ((i >> 1) & 1),
            voxel.z + stride * (i >> 2)
        };
        int64_t vid;
        if (svhash.vertIndexMap.lookup(svhash.getKey(vi), vid))
            vert_indices[i] = (int)vid;
    }
    out_vert_indices[tid] = vert_indices;

}

std::tuple<at::Tensor, at::Tensor>
svhashGetVoxels(SVHashTensor::TensorTuple tensorTuples) {
    SVHashTensor tensors(tensorTuples);

    auto kOptI32 = at::TensorOptions().dtype(torch::kInt32).device(torch::kCUDA);
    auto kOptF32 = at::TensorOptions().dtype(torch::kFloat32).device(torch::kCUDA);

    at::Tensor pos_sizes = at::empty({tensors.numCells(), 4}, kOptF32);
    at::Tensor vert_indices = at::empty({tensors.numCells(), 8}, kOptI32);

    svhashGetVoxels_kernel<<<_LAUNCH_ARGS_1D(tensors.numCells(), 256)>>>(
        SVHash(tensors),
        (float4*)pos_sizes.data_ptr<float>(),
        (FixedArray<int, 8>*)vert_indices.data_ptr<int>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return std::make_tuple(pos_sizes, vert_indices);
}




/* Densification - Split */

__global__ void svhashSplitVoxels_expandMask(
    const SVHash svhash,
    const bool* __restrict__ in_split_mask,
    bool* __restrict__ out_split_mask
) {
    uint tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= svhash.numCells)
        return;

    bool mask = in_split_mask[tid];
    int4 cell = svhash.cellList[tid];
    mask &= (cell.w > 0);

    // TODO: expand so split factor with neighbor cell is small

    out_split_mask[tid] = mask;
}

__global__ void svhashSplitVoxels_fillNewCells(
    const int64_t num_cells,
    const int64_t num_cells_split,
    const int64_t* __restrict__ split_indices,
    int4* __restrict__ new_cell_list
) {
    uint tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_cells_split)
        return;

    int64_t src_idx = split_indices[tid];
    int4 cell = new_cell_list[src_idx];
    int stride = 1 << (cell.w - 1);

    new_cell_list[src_idx] = {cell.x, cell.y, cell.z, cell.w-1};
    #pragma unroll
    for (int i = 1; i < 8; i++) {
        int4 new_cell = {
            cell.x + stride * (i & 1),
            cell.y + stride * ((i >> 1) & 1),
            cell.z + stride * (i >> 2),
            cell.w-1
        };
        int64_t dst_idx = num_cells + 7 * tid + (i-1);
        new_cell_list[dst_idx] = new_cell;
    }
}

template<bool isCountingPass>
__global__ void svhashSplitVoxels_fillNewVerts(
    const int64_t num_cells_split,
    const int64_t* __restrict__ split_indices,
    const int4* __restrict__ new_cell_list,
    _SVHashTable<int64_t> verts_set,
    int* __restrict__ counts,  // value or PSA
    _SVHashTable<int64_t> old_verts_set = _SVHashTable<int64_t>(),
    int3* __restrict__ vert_list = nullptr,
    FixedArray<uint32_t, 8>* __restrict__ interpolateIndices = nullptr,
    FixedArray<float, 8>* __restrict__ interpolateWeights = nullptr
) {
    uint tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_cells_split)
        return;

    int64_t src_idx = split_indices[tid];
    int4 cell = new_cell_list[src_idx];
    int stride = 1 << cell.w;

    int count = 0;
    if (!isCountingPass)
        count = (tid == 0 ? 0 : counts[tid-1]);

    FixedArray<uint32_t, 8> interpIndices;
    if (!isCountingPass) {
        #pragma unroll
        for (int i = 0; i < 8; i++) {
            int3 offset = { i & 1, (i >> 1) & 1, i >> 2 };
            int3 vert = {
                cell.x + 2 * stride * offset.x,
                cell.y + 2 * stride * offset.y,
                cell.z + 2 * stride * offset.z
            };
            int64_t idx = -1;
            old_verts_set.lookup(SVHash::getKey(vert), idx);
            interpIndices[i] = (uint32_t)idx;
        }
    }

    #pragma unroll
    for (int i = 0; i < 27; i++) {
        int3 offset = { i % 3, (i / 3) % 3, i / (3 * 3) };
        if ((offset.x & 1) + (offset.y & 1) + (offset.z & 1) == 0)
            continue;  // skip corners

        int3 new_vert = {
            cell.x + stride * offset.x,
            cell.y + stride * offset.y,
            cell.z + stride * offset.z
        };
        
        if (isCountingPass) {
            bool written = verts_set.insert<false>(SVHash::getKey(new_vert), (int64_t)tid);
            count += (int)written;
        }
        else {
            int64_t written_idx;
            if (verts_set.lookup(SVHash::getKey(new_vert), written_idx) && written_idx == (int64_t)tid) {
                vert_list[count] = new_vert;
                interpolateIndices[count] = interpIndices;

                // TODO: fused kernel to save VRAM
                float3 w = make_float3(offset.x, offset.y, offset.z) * 0.5f;
                FixedArray<float, 8> weights;
                weights[0] = (1.0f-w.z)*(1.0f-w.y)*(1.0f-w.x);
                weights[1] = (1.0f-w.z)*(1.0f-w.y)*w.x;
                weights[2] = (1.0f-w.z)*w.y*(1.0f-w.x);
                weights[3] = (1.0f-w.z)*w.y*w.x;
                weights[4] = w.z*(1.0f-w.y)*(1.0f-w.x);
                weights[5] = w.z*(1.0f-w.y)*w.x;
                weights[6] = w.z*w.y*(1.0f-w.x);
                weights[7] = w.z*w.y*w.x;
                interpolateWeights[count] = weights;

                ++count;
            }
        }
    }

    if (isCountingPass) {
        counts[tid] = count;
    }
}

std::tuple<SVHashTensor::TensorTuple, at::Tensor, at::Tensor, at::Tensor>
svhashSplitVoxels(
    SVHashTensor::TensorTuple tensorTuples,
    at::Tensor split_mask
) {
    SVHashTensor svhash(tensorTuples);
    DEVICE_GUARD(split_mask);
    CHECK_INPUT(split_mask);
    if (split_mask.numel() != svhash.numCells())
        AT_ERROR("Mismatch in number of cells between voxel grid and split mask");

    int64_t num_cells = svhash.numCells();
    int64_t num_verts = svhash.numVerts();

    auto optI64 = svhash.vertIndexMap.options();
    auto optI32 = optI64.dtype(torch::kInt32);
    auto optF32 = optI32.dtype(torch::kFloat32);

    // expand split mask
    at::Tensor split_mask_expanded = at::zeros_like(split_mask);
    svhashSplitVoxels_expandMask<<<_LAUNCH_ARGS_1D(num_verts, 256)>>>(
        SVHash(svhash),
        split_mask.data_ptr<bool>(),
        split_mask_expanded.data_ptr<bool>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    // add new cells
    at::Tensor split_indices = at::where(split_mask_expanded)[0];
    int64_t num_cells_split = split_indices.numel();
    if (num_cells_split == 0)
        return std::make_tuple(
            svhash.tuple(),
            split_indices.repeat_interleave(7),
            at::empty({0, 8}, optI32),
            at::empty({0, 8}, optF32)
        );
    int64_t new_num_cells = num_cells + 7 * num_cells_split;
    at::Tensor newCellList = at::empty({new_num_cells, 4}, optI32);
    cudaMemcpy(
        newCellList.data_ptr<int32_t>(), svhash.cellList.data_ptr<int32_t>(),
        4*num_cells*sizeof(int32_t), cudaMemcpyDeviceToDevice
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
    svhashSplitVoxels_fillNewCells<<<_LAUNCH_ARGS_1D(num_cells_split, 256)>>>(
        num_cells,
        num_cells_split,
        split_indices.data_ptr<int64_t>(),
        (int4*)newCellList.data_ptr<int32_t>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    // add new vertices - counting pass
    int64_t vert_hash_table_size = 24 * num_cells_split;  // at least 19x for 19 non-corner vertices
    at::Tensor newVertSet = at::empty({vert_hash_table_size, 2}, optI64);
    cudaMemset(newVertSet.data_ptr<int64_t>(), 0xff, newVertSet.numel()*sizeof(int64_t));
    at::Tensor newVertCounts = at::empty({num_cells_split}, optI32);
    svhashSplitVoxels_fillNewVerts<true><<<_LAUNCH_ARGS_1D(num_cells_split, 256)>>>(
        num_cells_split,
        split_indices.data_ptr<int64_t>(),
        (int4*)newCellList.data_ptr<int32_t>(),
        _SVHashTable<int64_t>(newVertSet),
        newVertCounts.data_ptr<int32_t>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    // std::cout << newVertCounts << std::endl;
    // std::cout << newVertSet << std::endl;

    // add new vertices - filling pass
    newVertCounts = at::cumsum(newVertCounts, 0).to(torch::kInt32);
    int64_t num_verts_added = (int64_t)newVertCounts[-1].item<int>();
    int64_t new_num_verts = num_verts + num_verts_added;
    at::Tensor newVertList = at::empty({new_num_verts, 3}, optI32);
    cudaMemcpy(
        newVertList.data_ptr<int32_t>(), svhash.vertList.data_ptr<int32_t>(),
        3*num_verts*sizeof(int32_t), cudaMemcpyDeviceToDevice
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
    at::Tensor newVertInterpIndices = at::empty({num_verts_added, 8}, optI32);
    at::Tensor newVertInterpWeights = at::empty({num_verts_added, 8}, optF32);
    svhashSplitVoxels_fillNewVerts<false><<<_LAUNCH_ARGS_1D(num_cells_split, 256)>>>(
        num_cells_split,
        split_indices.data_ptr<int64_t>(),
        (int4*)newCellList.data_ptr<int32_t>(),
        _SVHashTable<int64_t>(newVertSet),
        newVertCounts.data_ptr<int32_t>(),
        _SVHashTable<int64_t>(svhash.vertIndexMap),
        (int3*)newVertList.data_ptr<int32_t>() + num_verts,
        (FixedArray<uint32_t, 8>*)newVertInterpIndices.data_ptr<int32_t>(),
        (FixedArray<float, 8>*)newVertInterpWeights.data_ptr<float>()
    );

    // recompute hash tables
    svhash.vertList = newVertList;
    svhash.cellList = newCellList;
    svhash.vertIndexMap = computeIndexMap<int3>(newVertList);
    svhash.cellIndexMap= computeIndexMap<int4>(newCellList);
    return std::make_tuple(
        svhash.tuple(),
        split_indices.repeat_interleave(7),
        newVertInterpIndices,
        newVertInterpWeights
    );
}
