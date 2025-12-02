#include "SVHash.cuh"


__forceinline__ __device__ uint32_t _hash(uint32_t k) {
    // 32 bit Murmur3 hash
    k ^= k >> 16;
    k *= 0x85ebca6b;
    k ^= k >> 13;
    k *= 0xc2b2ae35;
    k ^= k >> 16;
    return k;
}

__forceinline__ __device__ uint32_t _probe(uint32_t k) {
    return k + 1;
}


// convert value to hash key (a hash will be applied afterwards anyway)

inline __device__ uint32_t mixi2(int2 k) {
    return _hash((uint)k.x) ^ (uint)k.y;
}

inline __device__ uint32_t mixi3(int3 k) {
    return _hash(_hash((uint)k.x) ^ (uint)k.y) ^ (uint)k.z;
}

inline __device__ uint32_t mixi4(int4 k) {
    return _hash(_hash(_hash((uint)k.x) ^ (uint)k.y) ^ (uint)k.z) ^ (uint)k.w;
}

inline __device__ uint32_t mixu2(uint2 k) {
    return _hash(k.x) ^ k.y;
}

inline __device__ uint32_t mixu3(uint3 k) {
    return _hash(_hash(k.x) ^ k.y) ^ k.z;
}

inline __device__ uint32_t mixu4(uint4 k) {
    return _hash(_hash(_hash(k.x) ^ k.y) ^ k.z) ^ k.w;
}



template<typename Value>
struct _SVHashTable {

    static constexpr uint32_t kEmpty = 0xffffffff;

    struct KeyValue {
        uint32_t key;
        Value value;
    };

    uint32_t capacity;
    KeyValue* __restrict__ data;

    _SVHashTable() : capacity(0) {}

    _SVHashTable(at::Tensor &tensor) {
        DEVICE_GUARD(tensor);
        CHECK_INPUT(tensor);
        static_assert(sizeof(Value) == sizeof(int32_t));
        capacity = tensor.numel() / 2;
        data = (KeyValue*)tensor.data_ptr<int32_t>();
    }

    __device__ bool lookup(uint32_t key, Value &out_value) const {
        uint32_t slot = _hash(key) % capacity;
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

    __device__ void insert(uint32_t key, Value value) {
        uint32_t slot = _hash(key) % capacity;
        while (true) {
            uint32_t prev = atomicCAS(&data[slot].key, kEmpty, key);
            if (prev == kEmpty || prev == key) {
                data[slot].value = value;
                return;
            }
            slot = _probe(slot) % capacity;
        }
    }

};

struct SVHash {

    static constexpr int kNumResolutions = 18;  // inclusive
    static constexpr int kUnitResolution = 9;
    struct Cell {
        int3 cell;  // xyz
        int resolution;  // 0 is unit cube, negative is larger, positive is smaller
    };

    uint32_t numVerts;
    uint32_t numCells;

    int3* __restrict__ vertList;
    int4* __restrict__ cellList;

    _SVHashTable<uint32_t> vertIndexMap;  // map xyz to index
    _SVHashTable<uint32_t> cellIndexMap;  // map cell to index

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
};


SVHashTensor::SVHashTensor(long num_verts, long num_cells) {
    auto kOptI32 = at::TensorOptions().dtype(torch::kInt32).device(torch::kCUDA);
    vertList = at::empty({num_verts, 3}, kOptI32);
    cellList = at::empty({num_cells, 4}, kOptI32);
    vertIndexMap = at::empty({2*num_verts, 2}, kOptI32);
    cellIndexMap = at::empty({2*num_cells, 2}, kOptI32);
    cudaMemset(vertIndexMap.data_ptr<int32_t>(), 0xff, vertIndexMap.numel()*sizeof(int32_t));
    cudaMemset(cellIndexMap.data_ptr<int32_t>(), 0xff, cellIndexMap.numel()*sizeof(int32_t));
}


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
    svhash.vertIndexMap.insert(mixi3({key.x, key.y, key.z}), vert_idx);

    if (xi == res.x || yi == res.y || zi == res.z)
        return;

    // insert cell
    uint cell_idx = (zi * res.y + yi) * res.x + xi;
    svhash.cellList[cell_idx] = key;
    svhash.cellIndexMap.insert(mixi4(key), cell_idx);

}


SVHashTensor::TensorTuple svhashCreateInitialVolume(
    std::tuple<float, float, float> p0_,
    std::tuple<float, float, float> p1_,
    std::tuple<unsigned, unsigned, unsigned> res_
) {
    float3 p0 = { std::get<0>(p0_), std::get<1>(p0_), std::get<2>(p0_) };
    float3 p1 = { std::get<0>(p1_), std::get<1>(p1_), std::get<2>(p1_) };
    uint3 res = { std::get<0>(res_), std::get<1>(res_), std::get<2>(res_) };

    long num_verts = (res.x+1)*(res.y+1)*(res.z+1);
    long num_cells = res.x*res.y*res.z;

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
        uint32_t vid;
        if (svhash.vertIndexMap.lookup(mixi3(vi), vid))
            vert_indices[i] = vid;
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
