#include "SplatTileIntersector.cuh"

#include <cuda_runtime.h>
#include <cub/cub.cuh>
#include <thrust/sort.h>
#include <thrust/device_ptr.h>
#include <vector>

#include <c10/cuda/CUDAGuard.h>

#include "helpers.cuh"

// #define DEBUG


struct Splat {
    float3 mean;
    glm::mat3 cov;

    __device__ __forceinline__ void getAABB(float3 &aabb_min, float3 &aabb_max) const {
        float3 extend = make_float3(sqrtf(cov[0][0]), sqrtf(cov[1][1]), sqrtf(cov[2][2]));
        if (isfinite(dot(mean, extend))) {
            aabb_min = mean - extend;
            aabb_max = mean + extend;
        } else {
            aabb_min = aabb_max = make_float3(0.0f);
        }
    }
};

__device__ __forceinline__ Splat loadSplat(unsigned splatIdx, const SplatBuffers& buffers) {

    glm::vec3 mean = buffers.means[splatIdx];
    glm::vec3 scales = buffers.scales[splatIdx];
    float opac = buffers.opacs[splatIdx];
    float4 quat = buffers.quats[splatIdx];

    // opac = 1.0f / (1.0f+__expf(-opac));
    // scales = { __expf(scales.x), __expf(scales.y), __expf(scales.z) };
    quat = normalize(quat);

    float extend = fmin(3.33f, sqrt(2.0f * __logf(opac / ALPHA_THRESHOLD)));

    glm::mat3 S = {
        scales.x, 0.0f, 0.0f,
        0.0f, scales.y, 0.0f,
        0.0f, 0.0f, scales.z
    };
    glm::mat3 R = quat_to_rotmat(quat);
    glm::mat3 M = extend * R * S;
    return {
        {mean.x, mean.y, mean.z},
        M * glm::transpose(M)
    };
}


struct Ray {
    float3 ro;
    float3 rd;

    __device__ __forceinline__ bool isOverlap(float3 aabb_min, float3 aabb_max) const {
        float3 aabb_center = 0.5f*(aabb_min+aabb_max);
        float3 aabb_size = 0.5f*(aabb_max-aabb_min);
        float3 m = 1.0f / rd;
        float3 n = m * (ro - aabb_center);
        float3 k = fabs(m) * aabb_size;
        float3 t1 = -n - k;
        float3 t2 = -n + k;
        float tn = fmax(fmax(t1.x, t1.y), t1.z);
        float tf = fmin(fmin(t2.x, t2.y), t2.z);
        return (tn < tf && tf > 0.0f);
    }

    // return negative if no overlap, strictly positive for sorting ID
    __device__ __forceinline__ float isOverlap(const Splat &splat) const {
        float3 aabb_min, aabb_max;
        splat.getAABB(aabb_min, aabb_max);
        if (!isOverlap(aabb_min, aabb_max))
            return -1.0f;
        return dot(splat.mean-ro, rd);  // negative if center is behind
    }

};

struct Tile {
    float x0, x1, y0, y1;
    glm::mat4x3 view;
    glm::vec3 ro, rd;
    glm::vec3 n0, n1, n2, n3;

    __device__ __forceinline__ void precompute() {
        glm::mat3 R = glm::transpose(glm::mat3(view));
        ro = -R * glm::vec3(view[3]);
        rd = R[2];
        glm::vec3 e0 = R * glm::vec3(x0, y0, 1.0f);
        glm::vec3 e1 = R * glm::vec3(x0, y1, 1.0f);
        glm::vec3 e2 = R * glm::vec3(x1, y1, 1.0f);
        glm::vec3 e3 = R * glm::vec3(x1, y0, 1.0f);
        n0 = glm::cross(e0, e1);
        n1 = glm::cross(e1, e2);
        n2 = glm::cross(e2, e3);
        n3 = glm::cross(e3, e0);
    }

    __device__ __forceinline__ bool isOverlap(float3 aabb_min, float3 aabb_max) const {

        float3 c_ = 0.5f*(aabb_min+aabb_max);
        float3 r_ = 0.5f*(aabb_max-aabb_min);
        glm::vec3 c = {c_.x, c_.y, c_.z};
        glm::vec3 r = {r_.x, r_.y, r_.z};

        // intersection test using separating axis theorem
        // has false positive, may be good enough in practice
        glm::vec3 roc = c - ro;
        float s0 = glm::dot(n0, roc) - glm::dot(r, glm::abs(n0));
        float s1 = glm::dot(n1, roc) - glm::dot(r, glm::abs(n1));
        float s2 = glm::dot(n2, roc) - glm::dot(r, glm::abs(n2));
        float s3 = glm::dot(n3, roc) - glm::dot(r, glm::abs(n3));
        float s = fmax(fmax(s0, s1), fmax(s2, s3));
        float sz = -glm::dot(rd, roc) - glm::dot(r, glm::abs(rd));
        return fmax(s, sz) < 0.0f;
    }

    // return negative if no overlap, strictly positive for sorting ID
    __device__ __forceinline__ float isOverlap(const Splat &splat) const {
        // TODO
        float3 aabb_min, aabb_max;
        splat.getAABB(aabb_min, aabb_max);
        if (!isOverlap(aabb_min, aabb_max))
            return -1.0f;
        glm::vec3 mean = {splat.mean.x, splat.mean.y, splat.mean.z};
        return glm::dot(mean-ro, rd);  // negative if center is behind
    }

};

__device__ __forceinline__ Tile loadTile(unsigned tileIdx, const TileBuffers buffers) {
    static_assert(sizeof(glm::mat4) == 16*sizeof(float));
    static_assert(sizeof(glm::mat3) == 9*sizeof(float));

    glm::mat4 view = glm::transpose(buffers.viewmats[tileIdx]);
    glm::mat3 intrins = buffers.Ks[tileIdx];  // take it as row major

    float fx = intrins[0][0];
    float fy = intrins[1][1];
    float cx = intrins[0][2];
    float cy = intrins[1][2];
    // printf("%f %f %f %f\n", fx, fy, cx, cy);

    Tile res = {
        -cx / fx, (buffers.width - cx) / fx,
        -cy / fy, (buffers.height - cy) / fy,
        glm::mat4x3(glm::vec3(view[0]), glm::vec3(view[1]), glm::vec3(view[2]), glm::vec3(view[3]))
    };
    res.precompute();
    return res;
}


__global__ void computeSplatAABB(
    const SplatBuffers& splats,
    float3* __restrict__ aabb,
    float3* __restrict__ aabb_reduced
) {
    unsigned gid = blockIdx.x * blockDim.x + threadIdx.x;
    if (gid >= splats.size)
        return;

    Splat g = loadSplat(gid, splats);
    float3 aabb_min, aabb_max;
    g.getAABB(aabb_min, aabb_max);
    if (aabb != nullptr) {
        aabb[2*gid+0] = aabb_min;
        aabb[2*gid+1] = aabb_max;
    }

    auto block = cg::this_thread_block();
    cg::thread_block_tile<WARP_SIZE> warp = cg::tiled_partition<WARP_SIZE>(block);

    aabb_min.x = cg::reduce(warp, aabb_min.x, cg::less<float>());
    if (warp.thread_rank() == 0) atomicMin((float*)aabb_reduced+0, aabb_min.x);
    aabb_min.y = cg::reduce(warp, aabb_min.y, cg::less<float>());
    if (warp.thread_rank() == 0) atomicMin((float*)aabb_reduced+1, aabb_min.y);
    aabb_min.z = cg::reduce(warp, aabb_min.z, cg::less<float>());
    if (warp.thread_rank() == 0) atomicMin((float*)aabb_reduced+2, aabb_min.z);

    aabb_max.x = cg::reduce(warp, aabb_max.x, cg::greater<float>());
    if (warp.thread_rank() == 0) atomicMax((float*)aabb_reduced+3, aabb_max.x);
    aabb_max.y = cg::reduce(warp, aabb_max.y, cg::greater<float>());
    if (warp.thread_rank() == 0) atomicMax((float*)aabb_reduced+4, aabb_max.y);
    aabb_max.z = cg::reduce(warp, aabb_max.z, cg::greater<float>());
    if (warp.thread_rank() == 0) atomicMax((float*)aabb_reduced+5, aabb_max.z);

}


template<uint BRANCH_FACTOR>
__device__ unsigned getLevel(
    float3 aabb_min, float3 aabb_max,
    float3 root_min, float3 root_max,
    unsigned num_levels
) {
    float3 size = aabb_max - aabb_min;
    float max_size = fmaxf(size.x, fmaxf(size.y, size.z));
    float3 root_size = root_max - root_min;
    float root_max_size = fmaxf(root_size.x, fmaxf(root_size.y, root_size.z));

    // will overlap with max 8 cells if root is cube
    float ratio = fmaxf(root_max_size / max_size, 1.0f);
    float level = __logf(ratio) / __logf(BRANCH_FACTOR);
    return min(max((unsigned)level, (unsigned)0), num_levels-1);
}

template<uint BRANCH_FACTOR>
__device__ __forceinline__ uint getSubcellOffset(uint3 subcell) {
    uint3 i = subcell % BRANCH_FACTOR;
    return (i.z * BRANCH_FACTOR + i.y) * BRANCH_FACTOR + i.x;
}

template<uint BRANCH_FACTOR>
__global__ void countCellOverlaps(
    const SplatBuffers& splats,
    float3 root_min, float3 root_max,
    unsigned num_levels,
    unsigned* __restrict__ overlap_counts
) {
    unsigned tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= splats.size) return;
    
    Splat g = loadSplat(tid, splats);
    float3 aabb_min, aabb_max;
    g.getAABB(aabb_min, aabb_max);
    
    unsigned level = getLevel<BRANCH_FACTOR>(aabb_min, aabb_max, root_min, root_max, num_levels);
    
    unsigned count_splat = 1;
    unsigned count_subcell = level;
    overlap_counts[tid] = count_splat + count_subcell;
}

__device__ __forceinline__ uint64_t insert_2_zeros_between_bits(uint64_t x) {
    x = (x | (x << 32)) & (uint64_t)0xFFFF00000000FFFFULL;
    x = (x | (x << 16)) & (uint64_t)0x00FF0000FF0000FFULL;
    x = (x | (x << 8))  & (uint64_t)0x100F00F00F00F00FULL;
    x = (x | (x << 4))  & (uint64_t)0x10C30C30C30C30C3ULL;
    x = (x | (x << 2))  & (uint64_t)0x1249249249249249ULL;
    return x;
}

template<uint BRANCH_FACTOR>
__device__ __forceinline__ uint64_t getCellKey(
    unsigned level,
    uint3 pos, bool isSplat
) {
    static_assert(BRANCH_FACTOR == 2);
    // 6 bit level, 1 bit splat vs cell, (64-6-1)/3=19 Morton bits in each dimension
    constexpr unsigned kMortonBitsPerDim = 19;
    uint64_t x = (uint64_t)(pos.x & ((1<<level)-1)) << (kMortonBitsPerDim - level);
    uint64_t y = (uint64_t)(pos.y & ((1<<level)-1)) << (kMortonBitsPerDim - level);
    uint64_t z = (uint64_t)(pos.z & ((1<<level)-1)) << (kMortonBitsPerDim - level);
  #if 0
    x = insert_2_zeros_between_bits(x) & (((uint64_t)1<<(3*kMortonBitsPerDim))-1);
    y = insert_2_zeros_between_bits(y) & (((uint64_t)1<<(3*kMortonBitsPerDim))-1);
    z = insert_2_zeros_between_bits(z) & (((uint64_t)1<<(3*kMortonBitsPerDim))-1);
    uint64_t morton = (x * 2 + y) * 2 + z;
  #else
    uint64_t morton = (((x << kMortonBitsPerDim) | y) << kMortonBitsPerDim) | z;
  #endif
    // printf("%u  %d %d %d  %llx\n", level, pos.x, pos.y, pos.z, morton);
    return ((
        (((uint64_t)level << (3*kMortonBitsPerDim)) | morton) << 1
    ) | (uint64_t)isSplat);
}

template<uint BRANCH_FACTOR>
__global__ void fillCellOverlaps(
    const SplatBuffers& splats,
    float3 root_min, float3 root_max,
    unsigned num_levels,
    unsigned* __restrict__ overlap_offsets,
    uint64_t* __restrict__ cell_keys,
    int32_t* __restrict__ splat_ids,
    uint8_t* __restrict__ subcell_masks,
    int32_t* __restrict__ subcell_ids,
    float3* __restrict__ subcell_aabb
) {
    unsigned tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= splats.size) return;
    
    Splat g = loadSplat(tid, splats);
    float3 aabb_min, aabb_max;
    g.getAABB(aabb_min, aabb_max);
    float3 aabb_center = (aabb_min + aabb_max) * 0.5f;

    unsigned level = getLevel<BRANCH_FACTOR>(aabb_min, aabb_max, root_min, root_max, num_levels);
    unsigned offset = overlap_offsets[tid];
    unsigned num_cells = overlap_offsets[tid+1] - offset;
    if (num_cells == 0)
        return;
    unsigned idx = 0;

    // Fill splat cells
    float scale = powf(BRANCH_FACTOR, level);
    float3 cell_size = (root_max - root_min) / scale;
    
  #if 0
    uint3 min_cell = make_uint3((aabb_min - root_min) / cell_size);
    uint3 max_cell = make_uint3((aabb_max - root_min) / cell_size);
    min_cell = clamp(min_cell, 0, (int)(scale-0.99f));
    max_cell = clamp(max_cell, 0, (int)(scale-0.99f));
    max_cell.x = min(max_cell.x, min_cell.x + 1);
    max_cell.y = min(max_cell.y, min_cell.y + 1);
    max_cell.z = min(max_cell.z, min_cell.z + 1);
    uint3 cells[8];  // max this number of cells guaranteed by getLevel
    int cellsRef[8];
  #else
    uint3 min_cell = make_uint3((aabb_center - root_min) / cell_size + 0.5f);
    uint3 max_cell = min_cell;
    uint3 cells[1];
    int cellsRef[1];
  #endif
    
    uint cellCount = 0;
    for (uint z = min_cell.z; z <= max_cell.z; z++) {
        for (uint y = min_cell.y; y <= max_cell.y; y++) {
            for (uint x = min_cell.x; x <= max_cell.x; x++) {
                cell_keys[offset+idx] = getCellKey<BRANCH_FACTOR>(level, {x, y, z}, true);
                splat_ids[offset+idx] = tid;
                subcell_masks[offset+idx] = (uint8_t)0;
                subcell_ids[offset+idx] = -1;
                subcell_aabb[2*(offset+idx)+0] = aabb_min;
                subcell_aabb[2*(offset+idx)+1] = aabb_max;
                cells[cellCount] = {x, y, z};
                cellsRef[cellCount] = offset+idx;
                if (++idx >= num_cells)
                    return;
                ++cellCount;
            }
        }
    }
    
    // Fill parent cells
    while (level--) {
        // reduce cell list while writing grid
        // notice cells are sorted by z, then y, then x
        uint i0 = 0;
        for (uint i = 0; i < cellCount; i++) {
            // uint8_t mask = (uint8_t)1 << getSubcellOffset<BRANCH_FACTOR>(cells[i]);
            uint8_t mask = (uint8_t)(1 + getSubcellOffset<BRANCH_FACTOR>(cells[i]));
            cells[i0] = (cells[i] >> 1);
            cell_keys[offset+idx] = getCellKey<BRANCH_FACTOR>(level, cells[i0], false);;
            splat_ids[offset+idx] = -1;
            subcell_masks[offset+idx] = mask;
            subcell_ids[offset+idx] = cellsRef[i];
            subcell_aabb[2*(offset+idx)+0] = aabb_min;
            subcell_aabb[2*(offset+idx)+1] = aabb_max;
            if (i0 == 0 || cells[i0-1] != cells[i0]) {
                cellsRef[i0] = offset+idx;
                ++i0;
            }
            if (++idx >= num_cells)
                return;
        }
        cellCount = i0;
        if (level == 0) break;
    }

    // set the rest empty
    while (idx < num_cells) {
        cell_keys[offset+idx] = (~((uint64_t)0)) >> 1;
        splat_ids[offset+idx] = -1;
        subcell_masks[offset+idx] = (uint8_t)0;
        subcell_ids[offset+idx] = -1;
        idx++;
    }
}



template<typename T>
__global__ void invertPermutation(
    size_t size,
    const T* __restrict__ perm,
    T* __restrict__ inverse
) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= size) return;
    T p = perm[idx];
    inverse[p] = (T)idx;
}

__global__ void gatherAndRemap(
    size_t size,
    const int32_t* __restrict__ subcell,
    const int32_t* __restrict__ perm,
    const int32_t* __restrict__ perm_inverse,
    int32_t* __restrict__ subcell_out
) {
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= size) return;
    int32_t temp = subcell[perm[i]];
    subcell_out[i] = (temp == -1) ? -1 : perm_inverse[temp];
}

__global__ void getCellDifferential(
    unsigned num_elem,
    const uint64_t* __restrict__ keys,
    int32_t* __restrict__ cell_id_differentials,
    int32_t* __restrict__ splat_differentials
) {
    unsigned tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid > num_elem) return;

    uint64_t key = keys[tid];

    cell_id_differentials[tid] = (tid == num_elem) ? 0 :
        (int)((key >> 1) != (keys[tid+1] >> 1));

    splat_differentials[tid] = (tid == num_elem) ? 0 :
        (int)(key & 1);
}

__global__ void fillTreeSplats(
    unsigned num_overlaps,
    const uint64_t* __restrict__ cell_keys,
    const int32_t* __restrict__ splat_ids,
    const unsigned* __restrict__ splat_idx_map,
    const unsigned* __restrict__ cell_id_map,
    unsigned* __restrict__ splatRanges,
    unsigned* __restrict__ splatIndices
) {
    unsigned tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_overlaps) return;

    // filter only splats (not cells)
    if ((int)cell_keys[tid] & 1 == 1) {

        // update splat idx
        unsigned splat_idx = splat_idx_map[tid];
        splatIndices[splat_idx] = splat_ids[tid];

        // update range
        if (tid == 0 || cell_keys[tid] != cell_keys[tid-1])
            splatRanges[cell_id_map[tid]*2] = splat_idx;
        if (tid == num_overlaps-1 || cell_keys[tid] != cell_keys[tid+1])
            splatRanges[cell_id_map[tid]*2+1] = splat_idx + 1;
    }
}


template<typename T>
__device__ __forceinline__ void lower_upper_bounds(
    const T *arr, unsigned n, T value, unsigned &lo, unsigned &hi
) {
    unsigned left = 0, right = n;

    // Find lower bound (first index >= value)
    while (left < right) {
        int mid = left + (right - left) / 2;
        if (arr[mid] < value)
            left = mid + 1;
        else
            right = mid;
    }
    lo = left;

    // Find upper bound (first index > value)
    right = n;
    while (left < right) {
        int mid = left + (right - left) / 2;
        if (arr[mid] <= value)
            left = mid + 1;
        else
            right = mid;
    }
    hi = left;
}



template<uint BRANCH_FACTOR>
__global__ void fillTreeSubcells_perCell(
    unsigned num_overlaps,
    unsigned num_cells,
    const unsigned* __restrict__ cell_id_map,
    const int32_t* __restrict__ subcell_ids,
    const uint8_t* __restrict__ subcell_masks,
    const float3* __restrict__ subcell_aabb,
    int32_t* __restrict__ children,
    float3* __restrict__ treeAABB
) {
    unsigned tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_cells) return;

    // binary search bounds [i0, i1)
    unsigned i0, i1;
    lower_upper_bounds(cell_id_map, num_overlaps, tid, i0, i1);
    // printf("%u [%u, %u)\n", tid, i0, i1);

    float3 aabbMin, aabbMax;

    // fill
    constexpr unsigned B3 = BRANCH_FACTOR * BRANCH_FACTOR * BRANCH_FACTOR;
    for (unsigned i = i0; i < i1; i++) {

        // update children
        int mask = (int)subcell_masks[i];
        if (mask != 0) {
            // unsigned cid = tid * B3 + (31 - __clz(mask));
            unsigned cid = tid * B3 + (mask - 1);
            unsigned sid = subcell_ids[i];
            children[cid] = cell_id_map[sid];
        }

        // update AABB
        if (i == i0) {
            aabbMin = subcell_aabb[2*i+0];
            aabbMax = subcell_aabb[2*i+1];
        } else {
            aabbMin = fmin(aabbMin, subcell_aabb[2*i+0]);
            aabbMax = fmax(aabbMax, subcell_aabb[2*i+1]);
        }
    }

    treeAABB[2*tid+0] = aabbMin;
    treeAABB[2*tid+1] = aabbMax;
}


__global__ void fillTreeSubcells_initAABB(
    unsigned num_cells,
    float3* __restrict__ treeAABB
) {
    unsigned tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_cells) return;
    treeAABB[2*tid+0] = make_float3(1e10);
    treeAABB[2*tid+1] = -make_float3(1e10);
}

template<uint BRANCH_FACTOR>
__global__ void fillTreeSubcells_perOverlap(
    unsigned num_overlaps,
    unsigned num_cells,
    const unsigned* __restrict__ cell_id_map,
    const int32_t* __restrict__ subcell_ids,
    const uint8_t* __restrict__ subcell_masks,
    const float3* __restrict__ subcell_aabb,
    int32_t* __restrict__ children,
    float3* __restrict__ treeAABB
) {
    unsigned oid = blockIdx.x * blockDim.x + threadIdx.x;
    if (oid >= num_overlaps) return;
    unsigned tid = cell_id_map[oid];

    constexpr unsigned B3 = BRANCH_FACTOR * BRANCH_FACTOR * BRANCH_FACTOR;

    // update children
    int mask = (int)subcell_masks[oid];
    if (mask != 0) {
        unsigned cid = tid * B3 + (mask - 1);
        unsigned sid = subcell_ids[oid];
        atomicMax(&children[cid], cell_id_map[sid]);
    }

    // update AABB
    float3 aabbMin = subcell_aabb[2*oid+0];
    atomicMin(&treeAABB[2*tid+0].x, aabbMin.x);
    atomicMin(&treeAABB[2*tid+0].y, aabbMin.y);
    atomicMin(&treeAABB[2*tid+0].z, aabbMin.z);
    float3 aabbMax = subcell_aabb[2*oid+1];
    atomicMax(&treeAABB[2*tid+1].x, aabbMax.x);
    atomicMax(&treeAABB[2*tid+1].y, aabbMax.y);
    atomicMax(&treeAABB[2*tid+1].z, aabbMax.z);
}


template<uint MAX_NUM_LEVELS, uint BRANCH_FACTOR>
__global__ void getTileSplatIntersections_octree(
    const TileBuffers tiles, const SplatBuffers& splats,
    const float3 rootAABBMin, const float3 rootAABBMax,
    const int32_t* __restrict__ children,
    const float3* __restrict__ treeAABB,
    const uint32_t* __restrict__ splatRanges,
    const uint32_t* __restrict__ splatIndices,
    uint32_t* __restrict__ intersect_counts,  // to be filled or exclusive scan
    uint32_t* __restrict__ intersectionSplatID  // nullptr or to be filled
) {
    static_assert(BRANCH_FACTOR == 2);
    constexpr uint kNumSubtree = BRANCH_FACTOR*BRANCH_FACTOR*BRANCH_FACTOR;  // 8
    constexpr uint kThreadsPerTile = BRANCH_FACTOR*BRANCH_FACTOR*BRANCH_FACTOR;  // 8
    constexpr uint kTilesPerWarp = 32 / kThreadsPerTile;  // 4

    unsigned tid = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned globalTileIdx = tid / kThreadsPerTile;

    unsigned warpTileIdx = threadIdx.x / kThreadsPerTile;
    unsigned tileThreadIdx = threadIdx.x % kThreadsPerTile;
    bool isActive = (globalTileIdx < tiles.size);

    bool isCountingPass = (intersectionSplatID == nullptr);
    uint32_t intersectGlobalOffset = 0, intersectGlobalOffsetMax = 0;
    if (!isCountingPass) {
        intersectGlobalOffset = intersect_counts[globalTileIdx];
        intersectGlobalOffsetMax = intersect_counts[globalTileIdx+1];
    }

    Tile tile;
    if (isActive)
        tile = loadTile(globalTileIdx, tiles);

    struct StackElem {
        uint32_t cellIdx;
        uint32_t level;
        uint3 offset;
    };
    __shared__ StackElem stack[kTilesPerWarp][(MAX_NUM_LEVELS+1)*kThreadsPerTile];
    uint stackSize = 0;
    if (isActive && tile.isOverlap(treeAABB[0], treeAABB[1])) {
        if (tileThreadIdx == 0)
            stack[warpTileIdx][stackSize] = {0, 0, make_uint3(0)};
        stackSize++;
    }
    __shared__ uint numSplatIntersects[kTilesPerWarp];
    if (tileThreadIdx == 0)
        numSplatIntersects[warpTileIdx] = 0;
    __syncwarp();

    for (uint _num_steps = 0; _num_steps < 65536; _num_steps++) {

        if (stackSize == 0)
            isActive = false;
        if (__ballot_sync(~0u, isActive) == 0)
            break;

        --stackSize;
        StackElem elem;
        if (isActive)
            elem = stack[warpTileIdx][stackSize];

        // Process splats
        if (isActive) {
            uint gi0 = splatRanges[2*elem.cellIdx+0], gi1 = splatRanges[2*elem.cellIdx+1];

            for (uint gi_ = gi0; gi_ < gi1; gi_ += kThreadsPerTile) {
                uint gi = gi_ + tileThreadIdx;
                if (gi < gi1) {
                    uint splatIdx = splatIndices[gi];
                    Splat splat = loadSplat(splatIdx, splats);
                    float overlap = tile.isOverlap(splat);
                    if (overlap > 0.0) {
                        uint idx = atomicAdd(&numSplatIntersects[warpTileIdx], 1) + intersectGlobalOffset;
                        if (idx < intersectGlobalOffsetMax) {
                            intersectionSplatID[idx] = splatIdx;
                        }
                    }
                }
            }
        }
        __syncwarp();
        if (!isCountingPass && numSplatIntersects[warpTileIdx] >= intersectGlobalOffsetMax-intersectGlobalOffset)
            isActive = false;
        if (__ballot_sync(~0u, isActive) == 0)
            break;

        // Process subcells
        #pragma unroll
        for (uint si_ = 0; si_ < kNumSubtree; si_ += kThreadsPerTile) {
            uint si = si_ + tileThreadIdx;
            int childIdx = isActive && si < kNumSubtree ?
                children[kNumSubtree*elem.cellIdx + si] : -1;
            int isActiveChild = int(childIdx >= 0);
            if (isActiveChild)
                isActiveChild &= tile.isOverlap(treeAABB[2*childIdx+0], treeAABB[2*childIdx+1]);

            int inclusiveActiveSum = isActiveChild;
            #pragma unroll
            for (unsigned offset = 1; offset < kThreadsPerTile; offset <<= 1) {
                int temp = __shfl_up_sync(~0u, inclusiveActiveSum, offset, kThreadsPerTile);
                if (tileThreadIdx >= offset)
                    inclusiveActiveSum += temp;
            }
            int exclusiveActiveSum = inclusiveActiveSum - isActiveChild;

            int last_lane = warpTileIdx * kThreadsPerTile + (kThreadsPerTile - 1);
            int activeSum = __shfl_sync(~0u, inclusiveActiveSum, last_lane, kThreadsPerTile);

            if (isActiveChild != 0) {
                uint3 delta = {
                    si / (MAX_NUM_LEVELS * MAX_NUM_LEVELS),
                    (si / MAX_NUM_LEVELS) % MAX_NUM_LEVELS,
                    si % MAX_NUM_LEVELS
                };
                stack[warpTileIdx][stackSize+exclusiveActiveSum] = {
                    (uint)childIdx,
                    elem.level + 1,
                    elem.offset * make_uint3(MAX_NUM_LEVELS) + delta
                };
            }
            stackSize += activeSum;
        }
        __syncwarp();
    }

    if (tileThreadIdx == 0) {
        if (isCountingPass)
            intersect_counts[globalTileIdx] = numSplatIntersects[warpTileIdx];
        else {
            uint32_t idx = numSplatIntersects[warpTileIdx] + intersectGlobalOffset;
            while (idx < intersectGlobalOffsetMax) {
                intersectionSplatID[idx] = 0;
                ++idx;
            }
        }
    }
}


__global__ void getTileSplatIntersections_brute(
    const TileBuffers tiles, const SplatBuffers& splats,
    uint32_t* __restrict__ intersect_counts,  // to be filled or exclusive scan
    uint32_t* __restrict__ intersectionSplatID  // nullptr or to be filled
) {
    unsigned tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= tiles.size)
        return;

    bool isCountingPass = (intersectionSplatID == nullptr);
    uint32_t intersectGlobalOffset = 0, intersectGlobalOffsetMax = 0;
    if (!isCountingPass) {
        intersectGlobalOffset = intersect_counts[tid];
        intersectGlobalOffsetMax = intersect_counts[tid+1];
    }
    uint32_t intersectCount = 0;

    Tile tile = loadTile(tid, tiles);

    for (uint32_t sid = 0; sid < splats.size; sid++) {
        Splat splat = loadSplat(sid, splats);
        float overlap = tile.isOverlap(splat);
        if (overlap > 0.0) {
            uint32_t idx = intersectGlobalOffset + intersectCount;
            intersectCount += 1;
            if (idx < intersectGlobalOffsetMax) {
                intersectionSplatID[idx] = sid;
            }
        }
    }

    if (isCountingPass)
        intersect_counts[tid] = intersectCount;
    else {
        uint32_t idx = intersectCount + intersectGlobalOffset;
        while (idx < intersectGlobalOffsetMax) {
            intersectionSplatID[idx] = 0;
            ++idx;
        }
    }
}


__device__ __forceinline__ uint64_t getSplatSortingKey(
    uint level, uint3 pos
) {
    // 7 bit level, (64-7-1)/3=19 Morton bits in each dimension
    constexpr uint kMortonBitsPerDim = 19;
    uint64_t x = (uint64_t)(pos.x & ((1<<level)-1)) << (kMortonBitsPerDim - level);
    uint64_t y = (uint64_t)(pos.y & ((1<<level)-1)) << (kMortonBitsPerDim - level);
    uint64_t z = (uint64_t)(pos.z & ((1<<level)-1)) << (kMortonBitsPerDim - level);
    x = insert_2_zeros_between_bits(x) & (((uint64_t)1<<(3*kMortonBitsPerDim))-1);
    y = insert_2_zeros_between_bits(y) & (((uint64_t)1<<(3*kMortonBitsPerDim))-1);
    z = insert_2_zeros_between_bits(z) & (((uint64_t)1<<(3*kMortonBitsPerDim))-1);
    uint64_t morton = (x * 2 + y) * 2 + z;
    return ((uint64_t)level << (3*kMortonBitsPerDim)) | morton;
}

__global__ void fillSplatSortingKeys(
    const SplatBuffers& splats,
    float3 root_min, float3 root_max,
    unsigned num_levels,
    uint64_t* __restrict__ splat_keys
) {
    unsigned tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= splats.size) return;
    
    Splat g = loadSplat(tid, splats);
    float3 aabb_min, aabb_max;
    g.getAABB(aabb_min, aabb_max);
    float3 aabb_center = (aabb_min + aabb_max) * 0.5f;

    // unsigned level = getLevel<2>(aabb_min, aabb_max, root_min, root_max, num_levels);
    unsigned level = num_levels-1;

    float scale = exp2f(level);
    float3 cell_size = (root_max - root_min) / scale;
    uint3 cell = make_uint3((aabb_center - root_min) / cell_size + 0.5f);

    uint64_t key = getSplatSortingKey(level, cell);
    splat_keys[tid] = key;
}

__global__ void fillLbvhInternalNodes(
    unsigned num_splats,
    const uint64_t* __restrict__ morton,
    const int32_t* __restrict__ splat_idx,
    int2* __restrict__ internal_nodes,
    int32_t* __restrict__ parent_nodes
) {
    // https://developer.nvidia.com/blog/parallelforall/wp-content/uploads/2012/11/karras2012hpg_paper.pdf
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= num_splats-1) return;

    #define delta(i, j) \
        (((j)<0 || (j)>=num_splats) ? -1 : morton[i] == morton[j] ? 64 + __clz(i ^ j) : __clzll(morton[i] ^ morton[j]))

    // Determine direction of the range (+1 or -1)
    int d = delta(i,i+1) - delta(i,i-1);
    d = d > 0 ? 1 : d < 0 ? -1 : 0;

    // Compute upper bound for the length of the range
    int delta_min = delta(i, i-d);
    int lmax = 2;
    while (delta(i, i+lmax*d) > delta_min)
        lmax <<= 1;

    // Find the other end using binary search
    int l = 0;
    for (int t = lmax>>1; t >= 1; t >>= 1)
        if (delta(i, i+(l+t)*d) > delta_min)
            l += t;
    int j = i + l * d;

    // Find the split position using binary search
    int delta_node = delta(i, j);
    int s = 0;
    for (int tf = 2, t; (t = (l+tf-1)/tf) >= 1; tf <<= 1)
        if (delta(i, i+(s+t)*d) > delta_node)
            s += t;
    int gamma = i + s*d + min(d, 0);

    // Output child pointers
    // regular for internal node, bit flip for leaf
    int left = min(i,j) == gamma ? ~splat_idx[gamma] : gamma;
    int right = max(i,j) == gamma+1 ? ~splat_idx[gamma+1] : gamma+1;
    internal_nodes[i] = make_int2(left, right);

    // Output parent pointers
    if (left >= 0)
        atomicMax(&parent_nodes[left], i);
    if (right >= 0)
        atomicMax(&parent_nodes[right], i);

    #undef delta
}

__global__ void computeLbvhAABB(
    const SplatBuffers& splats,
    const int2* __restrict__ internal_nodes,
    const int32_t* __restrict__ parent_nodes,
    float3* __restrict__ treeAABB
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= splats.size-1) return;

    int2 children = internal_nodes[i];
    if (children.x >= 0 && children.y >= 0)
        return;

    // find splat AABB
    float3 aabb_min, aabb_max;
    if (children.x < 0)
        loadSplat(~children.x, splats).getAABB(aabb_min, aabb_max);
    if (children.y < 0) {
        float3 aabb_min1, aabb_max1;
        loadSplat(~children.y, splats).getAABB(aabb_min1, aabb_max1);
        if (children.x < 0)
            aabb_min = fmin(aabb_min, aabb_min1),
            aabb_max = fmax(aabb_max, aabb_max1);
        else
            aabb_min = aabb_min1, aabb_max = aabb_max1;
    }

    // fill parent AABB
    do {
        #if 0
        atomicMin(&treeAABB[2*i].x, aabb_min.x);
        atomicMin(&treeAABB[2*i].y, aabb_min.y);
        atomicMin(&treeAABB[2*i].z, aabb_min.z);
        atomicMax(&treeAABB[2*i+1].x, aabb_max.x);
        atomicMax(&treeAABB[2*i+1].y, aabb_max.y);
        atomicMax(&treeAABB[2*i+1].z, aabb_max.z);
        #else
        if (atomicMin(&treeAABB[2*i].x, aabb_min.x) < aabb_min.x &
            atomicMin(&treeAABB[2*i].y, aabb_min.y) < aabb_min.y &
            atomicMin(&treeAABB[2*i].z, aabb_min.z) < aabb_min.z &
            atomicMax(&treeAABB[2*i+1].x, aabb_max.x) > aabb_max.x &
            atomicMax(&treeAABB[2*i+1].y, aabb_max.y) > aabb_max.y &
            atomicMax(&treeAABB[2*i+1].z, aabb_max.z) > aabb_max.z
        ) break;
        #endif
    } while ((i = parent_nodes[i]) >= 0);

}


__global__ void getTileSplatIntersections_lbvh(
    const TileBuffers tiles, const SplatBuffers& splats,
    const int2* __restrict__ internal_nodes,
    float3* __restrict__ treeAABB,
    uint32_t* __restrict__ intersect_counts,  // to be filled or exclusive scan
    uint32_t* __restrict__ intersectionSplatID  // nullptr or to be filled
) {
    // one thread per tile
    unsigned tileIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if (tileIdx >= tiles.size)
        return;
    uint laneIdx = tileIdx % WARP_SIZE;

    bool isCountingPass = (intersectionSplatID == nullptr);
    uint32_t intersectGlobalOffset = 0, intersectGlobalOffsetMax = 0;
    if (!isCountingPass) {
        intersectGlobalOffset = intersect_counts[tileIdx];
        intersectGlobalOffsetMax = intersect_counts[tileIdx+1];
    }

    Tile tile = loadTile(tileIdx, tiles);

    struct StackElem {
        uint32_t nodeIdx;
    };
    constexpr uint MAX_STACK_SIZE = 8*sizeof(int32_t)+1;
    __shared__ StackElem stack[WARP_SIZE][MAX_STACK_SIZE];
    uint stackSize = 0;
    if (tile.isOverlap(treeAABB[0], treeAABB[1])) {
        stack[laneIdx][stackSize] = { 0 };
        stackSize++;
    }
    uint numSplatIntersects = 0;

    for (uint _num_steps = 0; _num_steps < 65536; _num_steps++) {
        if (stackSize == 0)
            break;

        --stackSize;
        StackElem elem = stack[laneIdx][stackSize];
        int2 node = internal_nodes[elem.nodeIdx];
        // printf("[%u] stack %u - node %d %d\n", _num_steps, stackSize, node.x, node.y);

        for (uint ci = 0; ci < 2; ci++) {
            int childIdx = ci == 0 ? node.x : node.y;
            // splat
            if (childIdx < 0) {
                int splatIdx = ~childIdx;
                Splat splat = loadSplat(splatIdx, splats);
                float overlap = tile.isOverlap(splat);
                if (overlap > 0.0) {
                    uint idx = numSplatIntersects + intersectGlobalOffset;
                    if (idx < intersectGlobalOffsetMax) {
                        intersectionSplatID[idx] = splatIdx;
                    }
                    numSplatIntersects += 1;
                }
            }
            // node
            else if (tile.isOverlap(treeAABB[2*childIdx+0], treeAABB[2*childIdx+1])
                    && stackSize < MAX_STACK_SIZE) {
                stack[laneIdx][stackSize] = { (uint)childIdx };
                stackSize += 1;
            }
        }

    }

    if (isCountingPass)
        intersect_counts[tileIdx] = numSplatIntersects;
    else {
        uint32_t idx = numSplatIntersects + intersectGlobalOffset;
        while (idx < intersectGlobalOffsetMax) {
            intersectionSplatID[idx] = 0;
            ++idx;
        }
    }
}


__global__ void getTileSplatIntersections_lbvh_warp(
    const TileBuffers tiles, const SplatBuffers& splats,
    const int2* __restrict__ internal_nodes,
    float3* __restrict__ treeAABB,
    uint32_t* __restrict__ intersect_counts,  // to be filled or exclusive scan
    uint32_t* __restrict__ intersectionSplatID  // nullptr or to be filled
) {
    // one warp per tile, blockDim.x must be warp size
    unsigned tileIdx = blockIdx.x;
    if (tileIdx >= tiles.size)
        return;
    unsigned laneIdx = threadIdx.x % WARP_SIZE;

    bool isCountingPass = (intersectionSplatID == nullptr);
    uint32_t intersectGlobalOffset = 0, intersectGlobalOffsetMax = 0;
    if (!isCountingPass) {
        intersectGlobalOffset = intersect_counts[tileIdx];
        intersectGlobalOffsetMax = intersect_counts[tileIdx+1];
    }

    Tile tile = loadTile(tileIdx, tiles);

    struct StackElem {
        uint32_t nodeIdx;
    };
    constexpr uint MAX_STACK_SIZE = (8*sizeof(int32_t)+1)*WARP_SIZE;
    __shared__ StackElem stack[MAX_STACK_SIZE];
    uint stackSize = 0;
    if (tile.isOverlap(treeAABB[0], treeAABB[1])) {
        if (laneIdx == 0)
            stack[stackSize] = { 0 };
        stackSize++;
    }
    __syncwarp();
    uint numSplatIntersects = 0;

    auto block = cg::this_thread_block();
    cg::thread_block_tile<WARP_SIZE> warp = cg::tiled_partition<WARP_SIZE>(block);

    for (uint _num_steps = 0; _num_steps < 65536; _num_steps++) {
        if (stackSize == 0)
            break;

        int elemIdx = (int)stackSize - (int)laneIdx - 1;
        bool active = elemIdx >= 0;

        StackElem elem;
        if (active) elem = stack[elemIdx];

        stackSize -= __syncthreads_count(active);

        #pragma unroll
        for (uint ci = 0; ci < 2; ci++) {

            bool hasSplat = false;
            bool hasNode = false;
            int splatIdx = 0;
            int nodeIdx = 0;

            if (active) {
                int2 node = internal_nodes[elem.nodeIdx];

                int childIdx = ci == 0 ? node.x : node.y;
                // splat
                if (childIdx < 0) {
                    splatIdx = ~childIdx;
                    Splat splat = loadSplat(splatIdx, splats);
                    float overlap = tile.isOverlap(splat);
                    if (overlap > 0.0)
                        hasSplat = true;
                }
                // node
                else if (tile.isOverlap(treeAABB[2*childIdx+0], treeAABB[2*childIdx+1])) {
                    nodeIdx = childIdx;
                    hasNode = true;
                }
            }

            // fill splats
            if (!isCountingPass) {
                int inclusiveScan = (int)hasSplat;
                #pragma unroll
                for (unsigned offset = 1; offset < WARP_SIZE; offset <<= 1) {
                    int temp = __shfl_up_sync(~0u, inclusiveScan, offset, WARP_SIZE);
                    if (laneIdx >= offset)
                        inclusiveScan += temp;
                }
                int exclusiveScan = inclusiveScan - (int)hasSplat;

                uint idx = (numSplatIntersects + exclusiveScan) + intersectGlobalOffset;
                if (hasSplat && idx < intersectGlobalOffsetMax)
                    intersectionSplatID[idx] = splatIdx;
            }
            numSplatIntersects += __syncthreads_count(hasSplat);

            // fill nodes
            {
                int inclusiveScan = (int)hasNode;
                #pragma unroll
                for (unsigned offset = 1; offset < WARP_SIZE; offset <<= 1) {
                    int temp = __shfl_up_sync(~0u, inclusiveScan, offset, WARP_SIZE);
                    if (laneIdx >= offset)
                        inclusiveScan += temp;
                }
                int exclusiveScan = inclusiveScan - (int)hasNode;

                uint idx = stackSize + exclusiveScan;
                if (hasNode && idx < MAX_STACK_SIZE)
                    stack[idx] = { (uint)nodeIdx };
            }
            stackSize += __syncthreads_count(hasNode);
            stackSize = min(stackSize, MAX_STACK_SIZE);
        }

    }

    if (isCountingPass)
        intersect_counts[tileIdx] = numSplatIntersects;
    else {
        uint32_t idx = numSplatIntersects + intersectGlobalOffset;
        while (idx < intersectGlobalOffsetMax) {
            if (idx + laneIdx < intersectGlobalOffsetMax)
                intersectionSplatID[idx + laneIdx] = 0;
            idx += WARP_SIZE;
        }
    }
}


__forceinline__ torch::Tensor exclusiveScan(torch::Tensor &tensor) {
    torch::Tensor result = torch::empty_like(tensor);
    size_t temp_storage_bytes = 0;
    cub::DeviceScan::ExclusiveSum(nullptr, temp_storage_bytes,
        (unsigned*)tensor.data_ptr<int32_t>(),
        (unsigned*)result.data_ptr<int32_t>(),
        tensor.size(0));
    torch::Tensor temp_storage = torch::empty({(long)temp_storage_bytes}, tensor.options().dtype(torch::kUInt8));
    cub::DeviceScan::ExclusiveSum(temp_storage.data_ptr<uint8_t>(), temp_storage_bytes,
        (unsigned*)tensor.data_ptr<int32_t>(),
        (unsigned*)result.data_ptr<int32_t>(),
        tensor.size(0));
    return result;
}

__forceinline__ torch::Tensor invertPermutation(torch::Tensor &tensor) {
    constexpr uint block = 256;
    torch::Tensor result = torch::empty_like(tensor);
    invertPermutation<int32_t><<<(tensor.size(0)+block-1)/block, block>>>(
        tensor.size(0),
        tensor.data_ptr<int32_t>(),
        result.data_ptr<int32_t>()
    );
    return result;
}

template<typename T>
void _print_tensor(std::string name, torch::Tensor tensor) {
    cudaDeviceSynchronize();
    printf("%s ", name.c_str());
    // printf("\n"); return;
    tensor = tensor.cpu();
    if (tensor.ndimension() == 1) {
        for (unsigned i = 0; i < tensor.size(0); i++)
            printf("%lg ", (double)tensor.data_ptr<T>()[i]);
    }
    else if (tensor.ndimension() == 2) {
        for (unsigned i = 0; i < tensor.size(0); i++) {
            for (unsigned j = 0; j < tensor.size(1); j++)
                printf("%lg ", (double)tensor.data_ptr<T>()[i*tensor.size(1)+j]);
            printf(" ");
        }
    }
    printf("\n");
}

// #define print_tensor(dtype, tensor) _print_tensor<dtype>(#tensor, tensor)
#define print_tensor(dtype, tensor) ;


#ifdef DEBUG
void printAABB(float3 p0, float3 p1) {
    float3 c = -0.5f*(p0+p1), r = -0.5f*(p1-p0);
    printf("\\max\\left(\\left|x%+f\\right|%+f,\\left|y%+f\\right|%+f,\\left|z%+f\\right|%+f\\right)=0\n", c.x, r.x, c.y, r.y, c.z, r.z);
}

void printAABB_wireframe(float3 p0, float3 p1) {
    // B\left(x_{0},y_{0},z_{0},x_{1},y_{1},z_{1},t\right)=\left[\left(x_{0},y_{0},z_{0}\right),\left(x_{1},y_{0},z_{0}\right),\left(x_{1},y_{1},z_{0}\right),\left(x_{0},y_{1},z_{0}\right),\left(x_{0},y_{0},z_{0}\right),\left(x_{0},y_{1},z_{0}\right),\left(x_{0},y_{1},z_{1}\right),\left(x_{0},y_{0},z_{1}\right),\left(x_{0},y_{0},z_{0}\right),\left(x_{1},y_{0},z_{0}\right),\left(x_{1},y_{0},z_{1}\right),\left(x_{0},y_{0},z_{1}\right)\right]\left(1-t\right)+\left[\left(x_{0},y_{0},z_{1}\right),\left(x_{1},y_{0},z_{1}\right),\left(x_{1},y_{1},z_{1}\right),\left(x_{0},y_{1},z_{1}\right),\left(x_{1},y_{0},z_{0}\right),\left(x_{1},y_{1},z_{0}\right),\left(x_{1},y_{1},z_{1}\right),\left(x_{1},y_{0},z_{1}\right),\left(x_{0},y_{1},z_{0}\right),\left(x_{1},y_{1},z_{0}\right),\left(x_{1},y_{1},z_{1}\right),\left(x_{0},y_{1},z_{1}\right)\right]t
    printf("B\\left(%f,%f,%f,%f,%f,%f,t\\right)\n", p0.x, p0.y, p0.z, p1.x, p1.y, p1.z);
}

void printTile(glm::vec3 ro, glm::vec3 rd) {
    rd = glm::normalize(rd);
    printf("\\left(%f,%f,%f\\right)+\\left(%f,%f,%f\\right)100t\n",
        ro.x, ro.y, ro.z, rd.x, rd.y, rd.z);
}

void printCells(torch::Tensor cell_keys, torch::Tensor splat_ids, torch::Tensor subcell_masks, torch::Tensor subcell_ids) {
    int num_cells = cell_keys.size(0);
    printf("%d overlaps\n", num_cells);
    for (int i = 0; i < num_cells; i++) {
        uint64_t cell_key = cell_keys[i].item().toUInt64();
        int32_t splat_id = splat_ids[i].item().toInt();
        uint8_t subcell_mask = subcell_masks[i].item().toByte();
        int32_t subcell_id = subcell_ids[i].item().toInt();
        printf("[%d] layer %d, offset %llx, alias %b, splat %d, subcell %03b, %d\n", i, (int)(cell_key>>58), cell_key>>1, (unsigned)(cell_key&1), splat_id, (unsigned)(subcell_mask==0?0:31-__builtin_clz((unsigned)subcell_mask)), subcell_id);
    }
}

void checkMatch(std::string name, const torch::Tensor &a, const torch::Tensor &b) {
    if (a.size(0) != b.size(0)) {
        printf("%s: shape mismatch (%d != %d)\n", name.c_str(), (int)a.size(0), (int)b.size(0));
        return;
    }
    int numdiff = torch::abs(a - b).clip(0, 1).sum().item<int>();
    printf("%s: %d / %d mismatch\n", name.c_str(), numdiff, (int)a.size(0));
}
#endif


SplatTileIntersector::SplatTileIntersector(
    c10::TensorOptions tensorOptions,
    const SplatBuffers &splats,
    TileBuffers tiles
) : splats(splats), tiles(tiles)
{
    tensorF32 = tensorOptions.dtype(torch::kFloat32);
    tensorI32 = tensorOptions.dtype(torch::kInt32);
    tensorI16 = tensorOptions.dtype(torch::kInt16);
    tensorI64 = tensorOptions.dtype(torch::kInt64);
    tensorU8 = tensorOptions.dtype(torch::kUInt8);

    #ifdef DEBUG
    std::chrono::system_clock::time_point t0, t1;

    for (int i = 0; i < 1000; i++) {
        t0 = std::chrono::high_resolution_clock::now();
        // auto [icm1, sid1] = getIntersections_octree<12, 2>();
        auto [icm1, sid1] = getIntersections_lbvh();
        t1 = std::chrono::high_resolution_clock::now();
        printf("tree: %.2f ms\n", std::chrono::duration<float>(t1-t0).count()*1e3f);

        // continue;

        t0 = std::chrono::high_resolution_clock::now();
        auto [icm0, sid0] = getIntersections_brute();
        t1 = std::chrono::high_resolution_clock::now();
        printf("brute: %.2f ms\n", std::chrono::duration<float>(t1-t0).count()*1e3f);

        icm0 = icm0.cpu();
        icm1 = icm1.cpu();
        sid0 = sid0.cpu();
        sid1 = sid1.cpu();
        for (int k = 0; k+1 < icm0.size(0); k++) {
            int i0 = icm0[k].item<int>();
            int i1 = icm0[k+1].item<int>();
            std::sort(sid0.data_ptr<int>()+i0, sid0.data_ptr<int>()+i1);
            i0 = icm1[k].item<int>();
            i1 = icm1[k+1].item<int>();
            std::sort(sid1.data_ptr<int>()+i0, sid1.data_ptr<int>()+i1);
        }
        checkMatch("icm", icm1, icm0);
        checkMatch("sid", sid1, sid0);
    }
    exit(0);
    #endif
}

std::tuple<torch::Tensor, torch::Tensor> SplatTileIntersector::getIntersections_brute() {
    constexpr unsigned warp = 32;

    torch::Tensor intersection_count = torch::zeros({tiles.size+1}, tensorI32);
    getTileSplatIntersections_brute<<<(tiles.size+warp-1)/warp, warp>>>(
        tiles, splats,
        (uint32_t*)intersection_count.data_ptr<int32_t>(),
        nullptr
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
    cudaDeviceSynchronize();
    print_tensor(int, intersection_count);

    torch::Tensor intersection_count_map = exclusiveScan(intersection_count);
    print_tensor(int, intersection_count_map);
    unsigned total_intersections = (unsigned)intersection_count_map[tiles.size].item<int32_t>();

    torch::Tensor intersectionSplatID = torch::empty({total_intersections}, tensorI32);
    getTileSplatIntersections_brute<<<(tiles.size+warp-1)/warp, warp>>>(
        tiles, splats,
        (uint32_t*)intersection_count_map.data_ptr<int32_t>(),
        (uint32_t*)intersectionSplatID.data_ptr<int32_t>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
    print_tensor(int32_t, intersectionSplatID);

    return std::make_tuple(intersection_count_map, intersectionSplatID);
}


std::tuple<torch::Tensor, torch::Tensor> SplatTileIntersector::getIntersections_octree() {
    constexpr uint MAX_NUM_LEVELS = 12;
    constexpr uint BRANCH_FACTOR = 2;

    static_assert(MAX_NUM_LEVELS < 16);
    static_assert(BRANCH_FACTOR == 2 || BRANCH_FACTOR == 3 || BRANCH_FACTOR == 4);

    constexpr unsigned block = 256;
    constexpr unsigned warp = 32;
    constexpr int kFloatPInfByte = 0x7f;  // 0x7f7f7f7f -> 3.39615e+38
    constexpr int kFloatNInfByte = 0xfe;  // 0xfefefefe -> -1.69474e+38

    // find splat AABB
    // torch::Tensor splat_aabb = torch::empty({splats.size, 2, 3}, tensorF32);
    torch::Tensor root_aabb_tensor = torch::empty({2, 3}, tensorF32);
    cudaMemset(root_aabb_tensor.data_ptr<float>()+0, kFloatPInfByte, 3*sizeof(float));
    cudaMemset(root_aabb_tensor.data_ptr<float>()+3, kFloatNInfByte, 3*sizeof(float));
    computeSplatAABB<<<(splats.size+block-1)/block, block>>>(
        splats,
        // (float3*)splat_aabb.data_ptr<float>(),
        nullptr,
        (float3*)root_aabb_tensor.data_ptr<float>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
    #if 0
    {
        torch::Tensor splat_aabb_cpu = splat_aabb.cpu();
        for (int i = 0; i < splats.size; i++) {
            float3* ps = (float3*)splat_aabb_cpu.data_ptr<float>() + 2*i;
            printAABB(ps[0], ps[1]);
        }
    }
    #endif

    // find root AABB, pad them to cubes
    {
        root_aabb_tensor = root_aabb_tensor.cpu();
        float3* root_aabb = (float3*)root_aabb_tensor.data_ptr<float>();
        rootAABBMin = root_aabb[0];
        rootAABBMax = root_aabb[1];
        float3 center = 0.5f * (rootAABBMax + rootAABBMin);
        float3 extend = 0.5f * (rootAABBMax - rootAABBMin);
        float max_size = 1.01f * fmax(extend.x, fmax(extend.y, extend.z));
        rootAABBMin = center - make_float3(max_size);
        rootAABBMax = center + make_float3(max_size);        
    }
    // printAABB(rootAABBMin, rootAABBMax);
    // printf("%f %f %f  %f %f %f\n", rootAABBMin.x, rootAABBMin.y, rootAABBMin.z, rootAABBMax.x, rootAABBMax.y, rootAABBMax.z);

    // determine number of levels
    unsigned numLevels = MAX_NUM_LEVELS;
    
    // count number of cell overlaps for each splat
    torch::Tensor splat_cell_overlap_counts = torch::zeros({splats.size+1}, tensorI32);
    countCellOverlaps<BRANCH_FACTOR><<<(splats.size+block-1)/block, block>>>(
        splats, rootAABBMin, rootAABBMax, numLevels,
        (unsigned*)splat_cell_overlap_counts.data_ptr<int32_t>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
    print_tensor(int32_t, splat_cell_overlap_counts);

    // prefix sum to get offsets
    torch::Tensor splat_cell_overlap_offsets = exclusiveScan(splat_cell_overlap_counts);
    print_tensor(int32_t, splat_cell_overlap_offsets);
    unsigned total_overlaps = (unsigned)splat_cell_overlap_offsets[splats.size].item<int32_t>();
    // printf("%d\n", (int)total_overlaps);

    // fill overlap data
    torch::Tensor cell_keys = torch::zeros({total_overlaps}, tensorI64);
    torch::Tensor splat_ids = torch::zeros({total_overlaps}, tensorI32);
    torch::Tensor subcell_masks = torch::zeros({total_overlaps}, tensorU8);
    torch::Tensor subcell_ids = torch::zeros({total_overlaps}, tensorI32);
    torch::Tensor subcell_aabb = torch::zeros({total_overlaps, 2, 3}, tensorF32);
    fillCellOverlaps<BRANCH_FACTOR><<<(splats.size+block-1)/block, block>>>(
        splats, rootAABBMin, rootAABBMax, numLevels,
        (unsigned*)splat_cell_overlap_offsets.data_ptr<int32_t>(),
        (uint64_t*)cell_keys.data_ptr<int64_t>(),
        splat_ids.data_ptr<int32_t>(),
        subcell_masks.data_ptr<uint8_t>(),
        subcell_ids.data_ptr<int32_t>(),
        (float3*)subcell_aabb.data_ptr<float>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
    // printCells(cell_keys, splat_ids, subcell_masks, subcell_ids);

    // sort cells by keys
    // thrust::device_ptr<uint64_t> keys_ptr((uint64_t*)cell_keys.data_ptr<int64_t>());
    // thrust::device_ptr<int32_t> vals1_ptr(splat_ids.data_ptr<int32_t>());
    // thrust::device_ptr<uint64_t> vals2_ptr((uint64_t*)subcell_masks.data_ptr<int64_t>());
    // thrust::device_ptr<int32_t> vals3_ptr(subcell_ids.data_ptr<int32_t>());

    // auto first = thrust::make_zip_iterator(thrust::make_tuple(vals1_ptr, vals2_ptr, vals3_ptr));
    // thrust::sort_by_key(keys_ptr, keys_ptr + total_overlaps, first);

    auto [cell_keys_sorted, cell_keys_argsort] = torch::sort(cell_keys);
    cell_keys_argsort = cell_keys_argsort.to(torch::kInt32);
    cell_keys = cell_keys_sorted;
    splat_ids = torch::index(splat_ids, {cell_keys_argsort});
    subcell_masks = torch::index(subcell_masks, {cell_keys_argsort});
    subcell_aabb = torch::index(subcell_aabb, {cell_keys_argsort});
    #if 0
    subcell_ids = torch::index(subcell_ids, {cell_keys_argsort});
    subcell_ids = torch::where(subcell_ids == -1, subcell_ids, torch::index(cell_keys_argsort.argsort(), {subcell_ids})).contiguous().to(torch::kInt32);
    #else
    torch::Tensor cell_keys_argsort_argsort = invertPermutation(cell_keys_argsort);
    CHECK_DEVICE_ERROR(cudaGetLastError());
    torch::Tensor new_subcell_ids = torch::empty_like(subcell_ids);
    gatherAndRemap<<<(total_overlaps+block-1)/block, block>>>(
        total_overlaps,
        subcell_ids.data_ptr<int32_t>(),
        cell_keys_argsort.data_ptr<int32_t>(),
        cell_keys_argsort_argsort.data_ptr<int32_t>(),
        new_subcell_ids.data_ptr<int32_t>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
    subcell_ids = new_subcell_ids;
    #endif
    
    print_tensor(long, cell_keys_argsort);
    print_tensor(long, cell_keys_argsort.argsort());
    // printCells(cell_keys, splat_ids, subcell_masks, subcell_ids);

    // Get number of cells and splats as well as index map
    torch::Tensor cell_id_differential_map = torch::empty({total_overlaps+1}, tensorI32);
    torch::Tensor is_splat_map = torch::empty({total_overlaps+1}, tensorI32);
    getCellDifferential<<<(total_overlaps+block-1)/block, block>>>(
        total_overlaps,
        (uint64_t*)cell_keys.data_ptr<int64_t>(),
        cell_id_differential_map.data_ptr<int32_t>(),
        is_splat_map.data_ptr<int32_t>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    torch::Tensor cell_id_map = exclusiveScan(cell_id_differential_map);
    torch::Tensor splat_idx_map = exclusiveScan(is_splat_map);
    unsigned total_cells = (unsigned)cell_id_map[total_overlaps].item<int32_t>();
    unsigned total_splat_overlaps = (unsigned)splat_idx_map[total_overlaps].item<int32_t>();
    // print_tensor(int, cell_id_differential_map);
    print_tensor(int, cell_id_map);
    // print_tensor(int, is_splat_map);
    print_tensor(int, splat_idx_map);
    // printf("%d cells, %d splat overlaps\n", (int)total_cells, (int)total_splat_overlaps);

    // Fill splats
    torch::Tensor splatRanges = torch::zeros({total_cells, 2}, tensorI32);
    torch::Tensor splatIndices = torch::empty({total_splat_overlaps}, tensorI32);

    torch::Tensor treeAABB = torch::empty({total_cells, 2, 3}, tensorF32);
    fillTreeSplats<<<(total_overlaps+block-1)/block, block>>>(
        total_overlaps,
        (uint64_t*)cell_keys.data_ptr<int64_t>(),
        splat_ids.data_ptr<int32_t>(),
        (unsigned*)splat_idx_map.data_ptr<int32_t>(),
        (unsigned*)cell_id_map.data_ptr<int32_t>(),
        (unsigned*)splatRanges.data_ptr<int32_t>(),
        (unsigned*)splatIndices.data_ptr<int32_t>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
    print_tensor(int, splatRanges);
    print_tensor(int, splatIndices);
    print_tensor(uint8_t, subcell_masks);

    // Fill subcell map
    constexpr unsigned B3 = BRANCH_FACTOR * BRANCH_FACTOR * BRANCH_FACTOR;
    torch::Tensor children = torch::empty({total_cells, B3}, tensorI32);
    cudaMemset(children.data_ptr<int32_t>(), 0xff, total_cells*B3*sizeof(int32_t));
    
    print_tensor(int, cell_id_map);
    // print_tensor(int, torch::arange(cell_id_map.size(0)).cuda().to(torch::kInt32) - cell_id_map);
    print_tensor(int, subcell_ids);
    #if 0
    fillTreeSubcells_perCell<BRANCH_FACTOR><<<(total_cells+block-1)/block, block>>>(
        total_overlaps, total_cells,
        (unsigned*)cell_id_map.data_ptr<int32_t>(),
        subcell_ids.data_ptr<int32_t>(),
        subcell_masks.data_ptr<uint8_t>(),
        (float3*)subcell_aabb.data_ptr<float>(),
        children.data_ptr<int32_t>(),
        (float3*)treeAABB.data_ptr<float>()
    );
    #else
    fillTreeSubcells_initAABB<<<(total_cells+block-1)/block, block>>>(
        total_cells,
        (float3*)treeAABB.data_ptr<float>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
    fillTreeSubcells_perOverlap<BRANCH_FACTOR><<<(total_overlaps+block-1)/block, block>>>(
        total_overlaps, total_cells,
        (unsigned*)cell_id_map.data_ptr<int32_t>(),
        subcell_ids.data_ptr<int32_t>(),
        subcell_masks.data_ptr<uint8_t>(),
        (float3*)subcell_aabb.data_ptr<float>(),
        children.data_ptr<int32_t>(),
        (float3*)treeAABB.data_ptr<float>()
    );
    #endif
    CHECK_DEVICE_ERROR(cudaGetLastError());
    print_tensor(int, children);
    cudaDeviceSynchronize();
    // std::cout << treeAABB << std::endl;
    #ifdef DEBUG
    if (0) {
        torch::Tensor treeAABB_cpu = treeAABB.cpu();
        for (int i = 0; i < total_cells; i++) {
            float3* p = (float3*)treeAABB_cpu.data_ptr<float>() + 2*i;
            printAABB_wireframe(p[0], p[1]);
        }
    }
    #endif
    // return;

    // Traverse tree - get counts
    torch::Tensor intersection_count = torch::zeros({tiles.size+1}, tensorI32);
    print_tensor(int, intersection_count);
    getTileSplatIntersections_octree<MAX_NUM_LEVELS, BRANCH_FACTOR>
    <<<(tiles.size*B3+warp-1)/warp, warp>>>(
        tiles, splats, rootAABBMin, rootAABBMax,
        children.data_ptr<int32_t>(),
        (float3*)treeAABB.data_ptr<float>(),
        (unsigned*)splatRanges.data_ptr<int32_t>(),
        (unsigned*)splatIndices.data_ptr<int32_t>(),
        (uint32_t*)intersection_count.data_ptr<int32_t>(),
        nullptr
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
    cudaDeviceSynchronize();
    print_tensor(int, intersection_count);

    // Traverse tree - get offsets
    torch::Tensor intersection_count_map = exclusiveScan(intersection_count);
    print_tensor(int, intersection_count_map);
    unsigned total_intersections = (unsigned)intersection_count_map[tiles.size].item<int32_t>();

    // Traverse tree - write data
    torch::Tensor intersectionSplatID = torch::empty({total_intersections}, tensorI32);
    getTileSplatIntersections_octree<MAX_NUM_LEVELS, BRANCH_FACTOR>
    <<<(tiles.size*B3+warp-1)/warp, warp>>>(
        tiles, splats, rootAABBMin, rootAABBMax,
        children.data_ptr<int32_t>(),
        (float3*)treeAABB.data_ptr<float>(),
        (unsigned*)splatRanges.data_ptr<int32_t>(),
        (unsigned*)splatIndices.data_ptr<int32_t>(),
        (uint32_t*)intersection_count_map.data_ptr<int32_t>(),
        (uint32_t*)intersectionSplatID.data_ptr<int32_t>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
    print_tensor(int32_t, intersectionSplatID);

    return std::make_tuple(intersection_count_map, intersectionSplatID);
}


std::tuple<torch::Tensor, torch::Tensor> SplatTileIntersector::getIntersections_lbvh() {
    constexpr uint MAX_NUM_LEVELS = 12;

    static_assert(MAX_NUM_LEVELS < 32);

    constexpr unsigned block = 256;
    constexpr unsigned warp = 32;
    constexpr int kFloatPInfByte = 0x7f;  // 0x7f7f7f7f -> 3.39615e+38
    constexpr int kFloatNInfByte = 0xfe;  // 0xfefefefe -> -1.69474e+38

    // find splat AABB
    torch::Tensor splat_aabb = torch::empty({splats.size, 2, 3}, tensorF32);
    torch::Tensor root_aabb_tensor = torch::empty({2, 3}, tensorF32);
    cudaMemset(root_aabb_tensor.data_ptr<float>()+0, kFloatPInfByte, 3*sizeof(float));
    cudaMemset(root_aabb_tensor.data_ptr<float>()+3, kFloatNInfByte, 3*sizeof(float));
    computeSplatAABB<<<(splats.size+block-1)/block, block>>>(
        splats,
        (float3*)splat_aabb.data_ptr<float>(),
        (float3*)root_aabb_tensor.data_ptr<float>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
    #ifdef DEBUG
    if (0) {
        torch::Tensor splat_aabb_cpu = splat_aabb.cpu();
        for (int i = 0; i < splats.size; i++) {
            float3* ps = (float3*)splat_aabb_cpu.data_ptr<float>() + 2*i;
            printAABB(ps[0], ps[1]);
        }
    }
    #endif

    // find root AABB, pad them to cubes
    {
        root_aabb_tensor = root_aabb_tensor.cpu();
        float3* root_aabb = (float3*)root_aabb_tensor.data_ptr<float>();
        rootAABBMin = root_aabb[0];
        rootAABBMax = root_aabb[1];
        float3 center = 0.5f * (rootAABBMax + rootAABBMin);
        float3 extend = 0.5f * (rootAABBMax - rootAABBMin);
        float max_size = 1.01f * fmax(extend.x, fmax(extend.y, extend.z));
        rootAABBMin = center - make_float3(max_size);
        rootAABBMax = center + make_float3(max_size);
    }
    // printAABB_wireframe(rootAABBMin, rootAABBMax);

    // compute sorting keys (level and Morton code)
    torch::Tensor morton = torch::empty({splats.size}, tensorI64);
    fillSplatSortingKeys<<<(splats.size+block-1)/block, block>>>(
        splats, rootAABBMin, rootAABBMax, MAX_NUM_LEVELS,
        (uint64_t*)morton.data_ptr<int64_t>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
    cudaDeviceSynchronize();

    auto [sorted_morton, splat_argsort] = torch::sort(morton);
    splat_argsort = splat_argsort.to(torch::kInt32);

    #ifdef DEBUG
    if (0) {
        splat_aabb = torch::index(splat_aabb, {splat_argsort});
        // torch::Tensor splat_idx = invertPermutation(splat_argsort);
        torch::Tensor splat_aabb_cpu = splat_aabb.cpu();
        float3* aabb = (float3*)splat_aabb_cpu.data_ptr<float>();
        printf("\\left[");
        for (int i = 0; i < splats.size; i++) {
            float3 p = 0.5f*(aabb[2*i]+aabb[2*i+1]);
            printf("\\left(%f,%f,%f\\right),", p.x, p.y, p.z);
        }
        printf("\b\\right]\n");
    }
    #endif

    // Build tree
    torch::Tensor internal_nodes = torch::empty({splats.size-1, 2}, tensorI32);
    torch::Tensor parent_nodes = torch::empty({splats.size-1}, tensorI32);
    cudaMemset(parent_nodes.data_ptr<int32_t>(), 0xff, (splats.size-1)*sizeof(int32_t));
    fillLbvhInternalNodes<<<((splats.size-1)+block-1)/block, block>>>(
        splats.size,
        (uint64_t*)sorted_morton.data_ptr<int64_t>(),
        splat_argsort.data_ptr<int32_t>(),
        (int2*)internal_nodes.data_ptr<int32_t>(),
        parent_nodes.data_ptr<int32_t>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
    cudaDeviceSynchronize();
    print_tensor(int, internal_nodes);
    print_tensor(int, parent_nodes);

    // Compute AABB
    torch::Tensor treeAABB = torch::empty({splats.size, 2, 3}, tensorF32);
    fillTreeSubcells_initAABB<<<((splats.size-1)+block-1)/block, block>>>(
        splats.size-1,
        (float3*)treeAABB.data_ptr<float>()
    );
    computeLbvhAABB<<<((splats.size-1)+block-1)/block, block>>>(
        splats,
        (int2*)internal_nodes.data_ptr<int32_t>(),
        parent_nodes.data_ptr<int32_t>(),
        (float3*)treeAABB.data_ptr<float>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
    cudaDeviceSynchronize();
    #ifdef DEBUG
    if (0) {
        torch::Tensor treeAABB_cpu = treeAABB.cpu();
        for (int i = 0; i < splats.size-1; i++) {
            float3* p = (float3*)treeAABB_cpu.data_ptr<float>() + 2*i;
            printAABB_wireframe(p[0], p[1]);
        }
    }
    #endif

    // Traverse to find intersections
    torch::Tensor intersection_count = torch::zeros({tiles.size+1}, tensorI32);
    // getTileSplatIntersections_lbvh<<<(tiles.size+warp-1)/warp, warp>>>(
    getTileSplatIntersections_lbvh_warp<<<tiles.size, warp>>>(
        tiles, splats,
        (int2*)internal_nodes.data_ptr<int32_t>(),
        (float3*)treeAABB.data_ptr<float>(),
        (uint32_t*)intersection_count.data_ptr<int32_t>(),
        nullptr
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
    cudaDeviceSynchronize();
    print_tensor(int, intersection_count);

    torch::Tensor intersection_count_map = exclusiveScan(intersection_count);
    cudaDeviceSynchronize();
    print_tensor(int, intersection_count_map);
    unsigned total_intersections = (unsigned)intersection_count_map[tiles.size].item<int32_t>();

    torch::Tensor intersectionSplatID = torch::empty({total_intersections}, tensorI32);
    // getTileSplatIntersections_lbvh<<<(tiles.size+warp-1)/warp, warp>>>(
    getTileSplatIntersections_lbvh_warp<<<tiles.size, warp>>>(
        tiles, splats,
        (int2*)internal_nodes.data_ptr<int32_t>(),
        (float3*)treeAABB.data_ptr<float>(),
        (uint32_t*)intersection_count_map.data_ptr<int32_t>(),
        (uint32_t*)intersectionSplatID.data_ptr<int32_t>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
    cudaDeviceSynchronize();
    print_tensor(int32_t, intersectionSplatID);

    return std::make_tuple(intersection_count_map, intersectionSplatID);
}


#ifdef DEBUG

#include <stdio.h>

#include <regex>
#include <filesystem>
#include <fstream>

torch::Tensor loadBinaryToTensor(const std::string& filepath) {
    namespace fs = std::filesystem;
    
    // Extract filename from path
    fs::path p(filepath);
    std::string filename = p.filename().string();
    
    // Remove extension if present
    size_t dot_pos = filename.rfind('.');
    if (dot_pos != std::string::npos) {
        filename = filename.substr(0, dot_pos);
    }
    
    // Parse dimensions from filename using regex
    // Pattern: alphanumeric_name followed by numbers separated by underscores
    std::regex pattern(R"(^[a-zA-Z0-9]+_(.+)$)");
    std::smatch match;
    
    if (!std::regex_match(filename, match, pattern)) {
        throw std::runtime_error("Filename does not match expected pattern");
    }
    
    std::string dims_str = match[1].str();
    
    // Parse dimension values
    std::vector<int64_t> shape;
    std::stringstream ss(dims_str);
    std::string token;
    
    while (std::getline(ss, token, '_')) {
        if (!token.empty()) {
            try {
                shape.push_back(std::stoll(token));
            } catch (const std::exception& e) {
                throw std::runtime_error("Failed to parse dimension: " + token);
            }
        }
    }
    
    if (shape.empty()) {
        throw std::runtime_error("No dimensions found in filename");
    }
    
    // Calculate total number of elements
    int64_t total_elements = 1;
    for (int64_t dim : shape) {
        total_elements *= dim;
    }
    
    // Read binary file
    std::ifstream file(filepath, std::ios::binary | std::ios::ate);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file: " + filepath);
    }
    
    std::streamsize file_size = file.tellg();
    file.seekg(0, std::ios::beg);
    
    int64_t expected_bytes = total_elements * sizeof(float);
    if (file_size != expected_bytes) {
        throw std::runtime_error(
            "File size mismatch. Expected " + std::to_string(expected_bytes) +
            " bytes, but file contains " + std::to_string(file_size) + " bytes"
        );
    }
    
    // Allocate buffer and read data
    auto buffer = std::make_unique<float[]>(total_elements);
    file.read(reinterpret_cast<char*>(buffer.get()), file_size);
    file.close();
    
    if (!file) {
        throw std::runtime_error("Failed to read file completely");
    }
    
    // Create torch tensor on CUDA device
    torch::Tensor tensor = torch::from_blob(
        buffer.get(),
        shape,
        torch::kFloat32
    ).clone();
    
    return tensor.to(torch::kCUDA);
}

int main(int argc, char** argv) {
    auto seed = std::stoll(argv[3]);
    torch::manual_seed(seed);

    // unsigned num_splat = std::stod(argv[1]);
    // torch::Tensor means = torch::randn({num_splat, 3}).cuda();
    // torch::Tensor scales = 0.2*torch::randn({num_splat, 3}).cuda() - 0.4f*logf(num_splat) - 0.0f;
    // torch::Tensor opacs = torch::randn({num_splat, 1}).cuda();
    // torch::Tensor quats = torch::randn({num_splat, 4}).cuda();
    // unsigned num_tiles = std::stod(argv[2]);
    // auto [viewmats, Ks] = generate_random_camera_poses(num_tiles, seed);

    // if (0) {
    //     torch::Tensor tile_ro_cpu = tile_apex.cpu();
    //     torch::Tensor tile_rd_cpu = tile_dirs.cpu();
    //     for (unsigned i = 0; i < num_tiles; i++)
    //         printTile(((glm::vec3*)tile_ro_cpu.data_ptr<float>())[i],
    //         ((glm::vec3*)tile_rd_cpu.data_ptr<float>())[i]);
    // }

    // scales = torch::exp(scales);
    // opacs = torch::sigmoid(opacs);

    torch::Tensor means = loadBinaryToTensor("means_765390_3.bin");
    torch::Tensor scales = loadBinaryToTensor("scales_765390_3.bin");
    torch::Tensor opacs = loadBinaryToTensor("opacities_765390.bin");
    torch::Tensor quats = loadBinaryToTensor("quats_765390_4.bin");
    torch::Tensor viewmats = loadBinaryToTensor("viewmats_672_4_4.bin");
    torch::Tensor Ks = loadBinaryToTensor("Ks_672_3_3.bin");

    SplatTileIntersector::intersect_splat_tile(
        means, scales, opacs, quats,
        TILE_SIZE, TILE_SIZE,
        viewmats, Ks
    );
    return 0;
}

#endif


// /usr/local/cuda-12.8/bin/nvcc -I/media/harry/d/gs/spirulae-splat/spirulae_splat/splat/cuda/csrc/glm -I/home/harry/.venv/base/lib/python3.12/site-packages/torch/include -I/home/harry/.venv/base/lib/python3.12/site-packages/torch/include/torch/csrc/api/include -I/usr/local/cuda-12.8/include /media/harry/d/gs/spirulae-splat/spirulae_splat/splat/cuda/csrc/SplatTileIntersector.cu -o ./temp -D__CUDA_NO_HALF_OPERATORS__ -D__CUDA_NO_HALF_CONVERSIONS__ -D__CUDA_NO_BFLOAT16_CONVERSIONS__ -D__CUDA_NO_HALF2_OPERATORS__ -DDEBUG --expt-relaxed-constexpr --compiler-options ''"'"'-fPIC'"'"'' -O3 --use_fast_math --expt-relaxed-constexpr -Xcudafe=--diag_suppress=20012 -Xcudafe=--diag_suppress=550 -DTORCH_API_INCLUDE_EXTENSION_H -gencode=arch=compute_120,code=compute_120 -gencode=arch=compute_120,code=sm_120 -std=c++17 -L/home/harry/.venv/base/lib/python3.12/site-packages/torch/lib -L/usr/local/cuda-12.8/lib64 -L/usr/lib/x86_64-linux-gnu -lc10 -ltorch -ltorch_cpu -lcudart -lc10_cuda -ltorch_cuda

