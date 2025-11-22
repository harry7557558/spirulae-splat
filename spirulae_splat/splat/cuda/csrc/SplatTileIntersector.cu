#include "SplatTileIntersector.cuh"

#include <cuda_runtime.h>
#include <cub/cub.cuh>
#include <thrust/sort.h>
#include <thrust/device_ptr.h>
#include <vector>

#include <c10/cuda/CUDAGuard.h>

#include <gsplat/Utils.cuh>

#include "helpers.cuh"

// #define DEBUG


__device__ __forceinline__ float remapFunction(float x) {
    // identity near origin, proportional to x^(1/k) for large x
    constexpr float k = 2.5f;
    return k * sinhf( (1.0f/k) * asinhf( x ) );
}

__device__ __forceinline__ float3 remapAABB(float3 b) {
    return {
        remapFunction(b.x),
        remapFunction(b.y),
        remapFunction(b.z),
    };
}

template<bool remap>
__device__ bool getAABB(
    const Vanilla3DGS::World::Buffer& splatBuffer, long idx,
    float3 &aabb_min, float3 &aabb_max
) {
    float3 mean = splatBuffer.means[idx];
    float3 scales = splatBuffer.scales[idx];
    float opac = splatBuffer.opacities[idx];
    float4 quat = splatBuffer.quats[idx];

    opac = 1.0f / (1.0f+__expf(-opac));
    scales = { __expf(scales.x), __expf(scales.y), __expf(scales.z) };
    quat = normalize(quat);

    float3x3 covar;
    quat_scale_to_covar(quat, scales, &covar);

    float extend = fmin(3.33f, sqrt(2.0f * __logf(fmaxf(opac / ALPHA_THRESHOLD, 1.0f))));

    float3 bound = extend * make_float3(
        sqrtf(covar[0].x), sqrtf(covar[1].y), sqrtf(covar[2].z)
    );
    aabb_min = mean - bound;
    aabb_max = mean + bound;

    if (remap) {
        aabb_min = remapAABB(aabb_min);
        aabb_max = remapAABB(aabb_max);
    }
    return opac > ALPHA_THRESHOLD && isfinite(dot(aabb_min, aabb_max));
}

template<bool remap>
__device__ bool getAABB(
    const OpaqueTriangle::World::Buffer& splatBuffer, long idx,
    float3 &aabb_min, float3 &aabb_max
) {
    float3 mean = splatBuffer.means[idx];
    float3 scales = splatBuffer.scales[idx];
    float4 quat = splatBuffer.quats[idx];

    float3 vert0, vert1, vert2;
    map_opaque_triangle(mean, quat, scales, &vert0, &vert1, &vert2);

    aabb_min = {
        fmin(fmin(vert0.x, vert1.x), vert2.x),
        fmin(fmin(vert0.y, vert1.y), vert2.y),
        fmin(fmin(vert0.z, vert1.z), vert2.z),
    };
    aabb_max = {
        fmax(fmax(vert0.x, vert1.x), vert2.x),
        fmax(fmax(vert0.y, vert1.y), vert2.y),
        fmax(fmax(vert0.z, vert1.z), vert2.z),
    };

    if (remap) {
        aabb_min = remapAABB(aabb_min);
        aabb_max = remapAABB(aabb_max);
    }
    return isfinite(dot(aabb_min, aabb_max));
}


template<gsplat::CameraModelType camera_model>
struct Tile {
    glm::vec3 ro, rd;
    glm::vec3 n0, n1, n2, n3;

    __device__ bool init(
        CameraDistortionCoeffs dist_coeffs,
        float x0, float x1, float y0, float y1,
        glm::mat4x3 view
    ) {
        bool is_fisheye = (camera_model == gsplat::CameraModelType::FISHEYE);

        glm::mat3 R = glm::transpose(glm::mat3(view));
        ro = -R * glm::vec3(view[3]);
        rd = R[2];

        // TODO: better way to handle this in nonlinear and partially invalid case
        // May not matter in training with small tiles, but obvious artifact when rendering >180deg fisheye
        float3 e0_, e1_, e2_, e3_;
        bool valid0 = unproject_point({x0, y0}, is_fisheye, &dist_coeffs, &e0_);
        bool valid1 = unproject_point({x0, y1}, is_fisheye, &dist_coeffs, &e1_);
        bool valid2 = unproject_point({x1, y1}, is_fisheye, &dist_coeffs, &e2_);
        bool valid3 = unproject_point({x1, y0}, is_fisheye, &dist_coeffs, &e3_);
        if (!valid0 && valid3 && valid1 && valid2)
            e0_ = e3_ + e1_ - e2_, valid0 = true;
        if (!valid1 && valid0 && valid2 && valid3)
            e1_ = e0_ + e2_ - e3_, valid1 = true;
        if (!valid2 && valid1 && valid3 && valid0)
            e2_ = e1_ + e3_ - e0_, valid2 = true;
        if (!valid3 && valid2 && valid0 && valid1)
            e3_ = e2_ + e0_ - e1_, valid3 = true;
        glm::vec3 e0 = R * glm::vec3(e0_.x, e0_.y, e0_.z);
        glm::vec3 e1 = R * glm::vec3(e1_.x, e1_.y, e1_.z);
        glm::vec3 e2 = R * glm::vec3(e2_.x, e2_.y, e2_.z);
        glm::vec3 e3 = R * glm::vec3(e3_.x, e3_.y, e3_.z);

        n0 = glm::normalize(glm::cross(e0, e1));
        n1 = glm::normalize(glm::cross(e1, e2));
        n2 = glm::normalize(glm::cross(e2, e3));
        n3 = glm::normalize(glm::cross(e3, e0));
        return (valid0 && valid1 && valid2 && valid3 &&
            isfinite(glm::length(n0+n1+n2+n3)));
    }

    __device__ __forceinline__ bool isOverlap(float3 aabb_min, float3 aabb_max) const {

        float3 c_ = 0.5f*(aabb_min+aabb_max);
        float3 r_ = 0.5f*(aabb_max-aabb_min);
        glm::vec3 c = {c_.x, c_.y, c_.z};
        glm::vec3 r = {r_.x, r_.y, r_.z};

        // intersection test using separating axis theorem
        // has false positive; TODO: tighter one may help performance, most latency is global memory load
        glm::vec3 roc = c - ro;
        float s0 = glm::dot(n0, roc) - glm::dot(r, glm::abs(n0));
        float s1 = glm::dot(n1, roc) - glm::dot(r, glm::abs(n1));
        float s2 = glm::dot(n2, roc) - glm::dot(r, glm::abs(n2));
        float s3 = glm::dot(n3, roc) - glm::dot(r, glm::abs(n3));
        float s = fmax(fmax(s0, s1), fmax(s2, s3));
        if (camera_model != gsplat::CameraModelType::PINHOLE)
            return s < 0.0f;
        float sz = -glm::dot(rd, roc) - glm::dot(r, glm::abs(rd));
        return fmax(s, sz) < 0.0f;
    }

    // return negative if no overlap, strictly positive for sorting ID
    __device__ __forceinline__ float isOverlap(const Vanilla3DGS::World::Buffer& splatBuffer, long idx) const {
        // TODO: primitive aware version with less false positives
        float3 aabb_min, aabb_max;
        bool valid_aabb = getAABB<false>(splatBuffer, idx, aabb_min, aabb_max);
        if (!valid_aabb || !isOverlap(aabb_min, aabb_max))
            return -1.0f;
        float3 mean_ = 0.5f * (aabb_min + aabb_max);
        glm::vec3 mean(mean_.x, mean_.y, mean_.z);
        return camera_model != gsplat::CameraModelType::PINHOLE ?
            glm::length(mean - ro) :
            glm::dot(mean - ro, rd);  // negative if center is behind
    }

    __device__ __forceinline__ float isOverlap(const OpaqueTriangle::World::Buffer& splatBuffer, long idx) const {
        // TODO: primitive aware version with less false positives
        float3 aabb_min, aabb_max;
        bool valid_aabb = getAABB<false>(splatBuffer, idx, aabb_min, aabb_max);
        if (!valid_aabb || !isOverlap(aabb_min, aabb_max))
            return -1.0f;
        float3 mean_ = 0.5f * (aabb_min + aabb_max);
        glm::vec3 mean(mean_.x, mean_.y, mean_.z);
        return camera_model != gsplat::CameraModelType::PINHOLE ?
            glm::length(mean - ro) :
            glm::dot(mean - ro, rd);  // negative if center is behind
    }

};


template<gsplat::CameraModelType camera_model>
__device__ __forceinline__ Tile<camera_model>
loadTile(unsigned tileIdx, const TileBuffers<camera_model> buffers, bool& isActive) {
    static_assert(sizeof(glm::mat4) == 16*sizeof(float));
    static_assert(sizeof(glm::mat3) == 9*sizeof(float));

    glm::mat4 view = glm::transpose(buffers.viewmats[tileIdx]);
    glm::mat3 intrins = buffers.Ks[tileIdx];  // take it as row major

    float fx = intrins[0][0];
    float fy = intrins[1][1];
    float cx = intrins[0][2];
    float cy = intrins[1][2];
    // printf("%f %f %f %f\n", fx, fy, cx, cy);

    Tile<camera_model> res;
    isActive &= res.init(
        buffers.dist_coeffs.load(tileIdx),
        -cx / fx, (buffers.width - cx) / fx,
        -cy / fy, (buffers.height - cy) / fy,
        glm::mat4x3(glm::vec3(view[0]), glm::vec3(view[1]), glm::vec3(view[2]), glm::vec3(view[3]))
    );
    return res;
}


template<typename Primitive>
__global__ void computeSplatAABB(
    long numSplats,
    const typename Primitive::World::Buffer splatBuffer,
    float3* __restrict__ aabb,
    float3* __restrict__ aabb_reduced
) {
    unsigned gid = blockIdx.x * blockDim.x + threadIdx.x;
    if (gid >= numSplats)
        return;

    float3 aabb_min, aabb_max;
    bool valid_aabb = getAABB<true>(splatBuffer, gid, aabb_min, aabb_max);
    if (!valid_aabb)
        aabb_min = aabb_max = make_float3(0.0f);  // TODO: better way
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


__device__ unsigned getLevel(
    float3 aabb_min, float3 aabb_max,
    float3 root_min, float3 root_max,
    unsigned num_levels, float branch_factor
) {
    float3 size = aabb_max - aabb_min;
    float max_size = fmaxf(size.x, fmaxf(size.y, size.z));
    float3 root_size = root_max - root_min;
    float root_max_size = fmaxf(root_size.x, fmaxf(root_size.y, root_size.z));

    // will overlap with max 8 cells if root is cube
    float ratio = fmaxf(root_max_size / max_size, 1.0f);
    float level = __logf(ratio) / __logf(branch_factor);
    return min(max((unsigned)level, (unsigned)0), num_levels-1);
}

__device__ __forceinline__ uint64_t insert_2_zeros_between_bits(uint64_t x) {
    x = (x | (x << 32)) & (uint64_t)0xFFFF00000000FFFFULL;
    x = (x | (x << 16)) & (uint64_t)0x00FF0000FF0000FFULL;
    x = (x | (x << 8))  & (uint64_t)0x100F00F00F00F00FULL;
    x = (x | (x << 4))  & (uint64_t)0x10C30C30C30C30C3ULL;
    x = (x | (x << 2))  & (uint64_t)0x1249249249249249ULL;
    return x;
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



template<typename Primitive, gsplat::CameraModelType camera_model>
__global__ void getTileSplatIntersections_brute(
    const long numSplats,
    const TileBuffers<camera_model> tiles,
    const typename Primitive::World::Buffer splatBuffer,
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

    bool isActive = true;
    Tile<camera_model> tile = loadTile(tid, tiles, isActive);
    if (!isActive) {
        if (isCountingPass)
            intersect_counts[tid] = 0;
        return;
    }

    for (uint32_t sid = 0; sid < numSplats; sid++) {
        float overlap = tile.isOverlap(splatBuffer, sid);
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


// 7 bit level, (64-7-1)/3=19 Morton bits in each dimension
inline constexpr uint kMortonBitsPerDim = 19;

__device__ __forceinline__ uint64_t getSplatSortingKey(
    uint level, float3 pos
) {
    constexpr uint64_t mask_comp = ((uint64_t)1 << (uint64_t)kMortonBitsPerDim) - 1;
    uint64_t x = uint64_t(pos.x * exp2f(kMortonBitsPerDim - level) + 0.5f) & mask_comp;
    uint64_t y = uint64_t(pos.y * exp2f(kMortonBitsPerDim - level) + 0.5f) & mask_comp;
    uint64_t z = uint64_t(pos.z * exp2f(kMortonBitsPerDim - level) + 0.5f) & mask_comp;
    constexpr uint64_t mask_full = (((uint64_t)1 << (uint64_t)(3*kMortonBitsPerDim)) - 1);
    x = insert_2_zeros_between_bits(x) & mask_full;
    y = insert_2_zeros_between_bits(y) & mask_full;
    z = insert_2_zeros_between_bits(z) & mask_full;
    uint64_t morton = (x * 2 + y) * 2 + z;
    return ((uint64_t)level << (3*kMortonBitsPerDim)) | morton;
}

template<typename Primitive>
__global__ void fillSplatSortingKeys(
    const long numSplats,
    const typename Primitive::World::Buffer splatBuffer,
    float3 root_min, float3 root_max,
    unsigned num_levels, float branch_factor,
    uint64_t* __restrict__ splat_keys
) {
    unsigned tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= numSplats) return;
    
    float3 aabb_min, aabb_max;
    bool valid_aabb = getAABB<true>(splatBuffer, tid, aabb_min, aabb_max);
    if (!valid_aabb) {
        splat_keys[tid] = 0;
        return;
    }
    float3 aabb_center = (aabb_min + aabb_max) * 0.5f;

    unsigned level = getLevel(aabb_min, aabb_max, root_min, root_max, num_levels, branch_factor);
    // unsigned level = num_levels-1;

    float scale = exp2f(level);
    float3 cell_size = (root_max - root_min) / scale;
    float3 cell = (aabb_center - root_min) / cell_size;

    uint64_t key = getSplatSortingKey(level, cell);
    splat_keys[tid] = key;
}

__global__ void fillLbvhTreeRanges(
    unsigned num_levels,
    unsigned num_nodes,
    const uint64_t* __restrict__ keys,
    uint2* __restrict__ trees_ranges
) {
    unsigned tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_nodes) return;

    uint level0 = (tid > 0 ? keys[tid-1] : ~(uint64_t)0) >> (3*kMortonBitsPerDim);
    uint level1 = keys[tid] >> (3*kMortonBitsPerDim);
    if (level0 != level1)
        trees_ranges[level1].x = tid;

    level0 = level1;
    level1 = (tid+1 < num_nodes ? keys[tid+1] : ~(uint64_t)0) >> (3*kMortonBitsPerDim);
    if (level0 != level1)
        trees_ranges[level0].y = tid+1;
}

__global__ void sanitizeLbvhTreeRanges(
    unsigned num_levels,
    uint2* __restrict__ trees_ranges
) {
    unsigned level = threadIdx.x;  // <= warp size
    bool active = (level < num_levels);

    uint2 range = active ? trees_ranges[level] : make_uint2(~0u);

    // move nonempty ranges to beginning

    bool isNonEmpty = (range.y != range.x);
    int inclusiveScan = (int)isNonEmpty;
    #pragma unroll
    for (unsigned offset = 1; offset < WARP_SIZE; offset <<= 1) {
        int temp = __shfl_up_sync(~0u, inclusiveScan, offset, WARP_SIZE);
        if (threadIdx.x >= offset)
            inclusiveScan += temp;
    }
    int exclusiveScan = inclusiveScan - (int)isNonEmpty;

    __shared__ uint2 ranges[32];
    if (isNonEmpty)
        ranges[exclusiveScan] = range;
    __syncwarp();
    int numNonEmpty = __syncthreads_count(isNonEmpty);
    range = threadIdx.x < numNonEmpty ? ranges[threadIdx.x] : make_uint2(~0u);

    // merge ranges with size 1 to left
    // ranges with size 1 can be a problem in tree traversal

    for (uint iter = 0; iter < num_levels+1; iter++) {
        bool isInvalid = (range.y - range.x == 1);
        // printf("%u %u %d\n", iter, threadIdx.x, (int)isInvalid);
        unsigned invalidMask = __ballot_sync(~0u, isInvalid);
        if (invalidMask == 0)
            break;
        int elim = 31 - __clz(invalidMask);
        if (elim != 0)
            elim--;
        uint2 nextRange = { __shfl_down_sync(~0u, range.x, 1), __shfl_down_sync(~0u, range.y, 1) };
        range = (threadIdx.x < elim) ? range :
            (threadIdx.x == elim) ? make_uint2(range.x, nextRange.y) :
            (threadIdx.x < numNonEmpty-1) ? nextRange :
            make_uint2(~0u);
        numNonEmpty--;
    }

    if (threadIdx.x < num_levels)
        trees_ranges[threadIdx.x] = range;

}

__global__ void fillLbvhInternalNodes(
    unsigned num_levels,
    const uint2* __restrict__ trees_ranges,
    const uint64_t* __restrict__ morton,
    const int32_t* __restrict__ splat_idx,
    int2* __restrict__ internal_nodes,
    int32_t* __restrict__ parent_nodes
) {
    // https://developer.nvidia.com/blog/parallelforall/wp-content/uploads/2012/11/karras2012hpg_paper.pdf
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    uint2 range = {0, 0};
    for (uint l = 0; l < num_levels; l++) {
        range = trees_ranges[l];
        if (i >= range.x && i < range.y) {
            i -= range.x;
            morton += range.x;
            splat_idx += range.x;
            internal_nodes += range.x;
            parent_nodes += range.x;
            break;
        }
    }
    int num_splats = (int)(range.y - range.x);
    if (i >= num_splats-1)
        return;

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

template<typename Primitive>
__global__ void computeLbvhAABB(
    const long numSplats,
    const typename Primitive::World::Buffer splatBuffer,
    unsigned num_levels,
    const uint2* __restrict__ trees_ranges,
    const int2* __restrict__ internal_nodes,
    const int32_t* __restrict__ parent_nodes,
    float3* __restrict__ treeAABB
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numSplats-1) return;

    uint2 range = {0, 0};
    for (uint l = 0; l < num_levels; l++) {
        range = trees_ranges[l];
        if (range.x == range.y)
            continue;
        if (i >= range.x && i < range.y) {
            i -= range.x;
            internal_nodes += range.x;
            parent_nodes += range.x;
            treeAABB += 2 * range.x;
            break;
        }
    }
    if (i >= (int)(range.y-range.x)-1)
        return;

    int2 children = internal_nodes[i];
    if (children.x >= 0 && children.y >= 0)
        return;

    // find splat AABB
    float3 aabb_min, aabb_max;
    bool valid_aabb = true;
    if (children.x < 0)
        valid_aabb = getAABB<false>(splatBuffer, ~children.x, aabb_min, aabb_max);
    if (children.y < 0) {
        float3 aabb_min1, aabb_max1;
        valid_aabb = getAABB<false>(splatBuffer, ~children.y, aabb_min1, aabb_max1);
        if (children.x < 0)
            aabb_min = fmin(aabb_min, aabb_min1),
            aabb_max = fmax(aabb_max, aabb_max1);
        else
            aabb_min = aabb_min1, aabb_max = aabb_max1;
    }
    if (!valid_aabb)
        return;

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


template<typename Primitive, gsplat::CameraModelType camera_model>
__global__ void getTileSplatIntersections_lbvh(
    const TileBuffers<camera_model> tiles,
    const typename Primitive::World::Buffer splatBuffer,
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

    bool isActive = true;
    Tile<camera_model> tile = loadTile(tileIdx, tiles, isActive);
    if (!isActive) {
        if (isCountingPass)
            intersect_counts[tileIdx] = 0;
        return;
    }

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
                float overlap = tile.isOverlap(splatBuffer, splatIdx);
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


template<typename Primitive, gsplat::CameraModelType camera_model>
__global__ void getTileSplatIntersections_lbvh_warp(
    const TileBuffers<camera_model> tiles,
    const typename Primitive::World::Buffer splatBuffer,
    unsigned num_levels,
    const uint2* __restrict__ trees_ranges,
    const int2* __restrict__ internal_nodes_0,
    float3* __restrict__ treeAABB_0,
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

    bool isActive = true;
    Tile<camera_model> tile = loadTile(tileIdx, tiles, isActive);
    if (__ballot_sync(~0u, isActive) != ~0u) {
        if (isCountingPass) {
            if (laneIdx == 0)
                intersect_counts[tileIdx] = 0;
        } else {
            uint32_t idx = intersectGlobalOffset;
            while (idx < intersectGlobalOffsetMax) {
                if (idx + laneIdx < intersectGlobalOffsetMax)
                    intersectionSplatID[idx + laneIdx] = 0;
                idx += WARP_SIZE;
            }
        }
        return;
    }

    struct StackElem {
        uint32_t nodeIdx;
    };
    constexpr uint MAX_STACK_SIZE = (8*sizeof(int32_t)+1)*WARP_SIZE;
    __shared__ StackElem stack[MAX_STACK_SIZE];

    uint numSplatIntersects = 0;

  for (uint level = 0; level < num_levels; level++) {
    uint2 range = trees_ranges[level];
    if (range.y-range.x == 0)
        continue;
    float3* treeAABB = treeAABB_0 + 2*range.x;
    const int2* internal_nodes = internal_nodes_0 + range.x;

    // handle this case where treeAABB may be uninitialized
    if (range.y-range.x == 1) {
        // if (laneIdx == 0)
        // printf("range.y-range.x == 1: %u %u %u\n", level, range.x, range.y);
        continue;
    }

    uint stackSize = 0;
    if (tile.isOverlap(treeAABB[0], treeAABB[1])) {
        if (laneIdx == 0)
            stack[stackSize] = { 0 };
        stackSize++;
    }
    __syncwarp();

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
                    float overlap = tile.isOverlap(splatBuffer, splatIdx);
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

  }  // for (uint level = 0; level < num_levels; level++)

    if (isCountingPass) {
        if (laneIdx == 0)
            intersect_counts[tileIdx] = numSplatIntersects;
    }
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

void clearL2Cache() {
    size_t l2_cache_size = 128 << 20;
    void* temp_buffer;
    cudaMalloc(&temp_buffer, l2_cache_size);
    cudaMemset(temp_buffer, 0, l2_cache_size);
    cudaFree(temp_buffer);
}
#endif


template<typename Primitive, gsplat::CameraModelType camera_model>
SplatTileIntersector<Primitive, camera_model>::SplatTileIntersector(
    const typename Primitive::World::Tensor &splats,
    const TileBuffers<camera_model> &tiles
) : tiles(tiles), splats(splats)
{
    if (splats.batchSize() != 1)
        AT_ERROR("Patched mode only supports splat batch size 1");
    this->numSplats = splats.size();
    
    this->tensorF32 = splats.options().dtype(torch::kFloat32);
    this->tensorI32 = splats.options().dtype(torch::kInt32);
    this->tensorI16 = splats.options().dtype(torch::kInt16);
    this->tensorI64 = splats.options().dtype(torch::kInt64);
    this->tensorU8 = splats.options().dtype(torch::kUInt8);

    #ifdef DEBUG
    std::chrono::system_clock::time_point t0, t1;

    for (int i = 0; i < 10; i++) {
        clearL2Cache();
        t0 = std::chrono::high_resolution_clock::now();
        // auto [icm1, sid1] = getIntersections_octree<12, 2>();
        auto [icm1, sid1] = getIntersections_lbvh();
        t1 = std::chrono::high_resolution_clock::now();
        printf("tree: %.2f ms\n", std::chrono::duration<float>(t1-t0).count()*1e3f);

        // continue;

        clearL2Cache();
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

template<typename Primitive, gsplat::CameraModelType camera_model>
std::tuple<torch::Tensor, torch::Tensor> SplatTileIntersector<Primitive, camera_model>::getIntersections_brute() {
    constexpr unsigned warp = 32;

    torch::Tensor intersection_count = torch::zeros({tiles.size+1}, tensorI32);
    getTileSplatIntersections_brute<<<(tiles.size+warp-1)/warp, warp>>>(
        tiles, splats,
        (uint32_t*)intersection_count.data_ptr<int32_t>(),
        nullptr
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
    print_tensor(int, intersection_count);

    torch::Tensor intersection_count_map = exclusiveScan(intersection_count);
    print_tensor(int, intersection_count_map);
    unsigned total_intersections = (unsigned)intersection_count_map[(long)tiles.size].item<int32_t>();

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


template<typename Primitive, gsplat::CameraModelType camera_model>
std::tuple<torch::Tensor, torch::Tensor> SplatTileIntersector<Primitive, camera_model>::getIntersections_lbvh() {
    // TODO: use a separate rotated AABB aligned with (1,1,1) for thin off-diagnoal Gaussians?
    constexpr uint MAX_NUM_LEVELS = 28;
    constexpr float BRANCH_FACTOR = 2.0f;
    static_assert(MAX_NUM_LEVELS < 32);

    constexpr unsigned block = 256;
    constexpr unsigned warp = 32;
    constexpr int kFloatPInfByte = 0x7f;  // 0x7f7f7f7f -> 3.39615e+38
    constexpr int kFloatNInfByte = 0xfe;  // 0xfefefefe -> -1.69474e+38

    cudaDeviceSynchronize();

    // find splat AABB
    torch::Tensor splat_aabb = torch::empty({numSplats, 2, 3}, tensorF32);
    torch::Tensor root_aabb_tensor = torch::empty({2, 3}, tensorF32);
    cudaMemset(root_aabb_tensor.data_ptr<float>()+0, kFloatPInfByte, 3*sizeof(float));
    cudaMemset(root_aabb_tensor.data_ptr<float>()+3, kFloatNInfByte, 3*sizeof(float));
    computeSplatAABB<Primitive><<<(numSplats+block-1)/block, block>>>(
        numSplats, splats,
        (float3*)splat_aabb.data_ptr<float>(),
        (float3*)root_aabb_tensor.data_ptr<float>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
    #ifdef DEBUG
    if (0) {
        torch::Tensor splat_aabb_cpu = splat_aabb.cpu();
        for (int i = 0; i < numSplats; i++) {
            float3* ps = (float3*)splat_aabb_cpu.data_ptr<float>() + 2*i;
            printAABB(ps[0], ps[1]);
        }
    }
    #endif

    // find root AABB, pad them to cubes
    glm::vec3 rootAABBMin, rootAABBMax;
    {
        root_aabb_tensor = root_aabb_tensor.cpu();
        glm::vec3* root_aabb = (glm::vec3*)root_aabb_tensor.data_ptr<float>();
        rootAABBMin = root_aabb[0];
        rootAABBMax = root_aabb[1];
        glm::vec3 center = 0.5f * (rootAABBMax + rootAABBMin);
        glm::vec3 extend = 0.5f * (rootAABBMax - rootAABBMin);
        float max_size = 1.01f * fmax(extend.x, fmax(extend.y, extend.z));
        rootAABBMin = center - glm::vec3(max_size);
        rootAABBMax = center + glm::vec3(max_size);
    }
    // printf("AABB: %f %f %f  %f %f %f\n", rootAABBMin.x, rootAABBMin.y, rootAABBMin.z, rootAABBMax.x, rootAABBMax.y, rootAABBMax.z);
    // printAABB_wireframe(rootAABBMin, rootAABBMax);

    // compute sorting keys (level and Morton code)
    torch::Tensor morton = torch::empty({numSplats}, tensorI64);
    fillSplatSortingKeys<Primitive><<<(numSplats+block-1)/block, block>>>(
        numSplats, splats,
        *(float3*)&rootAABBMin, *(float3*)&rootAABBMax, MAX_NUM_LEVELS, BRANCH_FACTOR,
        (uint64_t*)morton.data_ptr<int64_t>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    auto [sorted_morton, splat_argsort] = torch::sort(morton);
    splat_argsort = splat_argsort.to(torch::kInt32);

    torch::Tensor tree_ranges = torch::empty({MAX_NUM_LEVELS, 2}, tensorI32);
    cudaMemset(tree_ranges.data_ptr<int32_t>(), 0xff, (2*MAX_NUM_LEVELS)*sizeof(int32_t));
    fillLbvhTreeRanges<<<(numSplats+block-1)/block, block>>>(
        MAX_NUM_LEVELS, numSplats,
        (uint64_t*)sorted_morton.data_ptr<int64_t>(),
        (uint2*)tree_ranges.data_ptr<int32_t>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
    // std::cout << tree_ranges << std::endl;
    sanitizeLbvhTreeRanges<<<1, WARP_SIZE>>>(
        MAX_NUM_LEVELS,
        (uint2*)tree_ranges.data_ptr<int32_t>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
    // std::cout << tree_ranges << std::endl;
    // cudaDeviceSynchronize();
    // exit(0);

    #ifdef DEBUG
    if (0) {
        splat_aabb = torch::index(splat_aabb, {splat_argsort});
        // torch::Tensor splat_idx = invertPermutation(splat_argsort);
        torch::Tensor splat_aabb_cpu = splat_aabb.cpu();
        float3* aabb = (float3*)splat_aabb_cpu.data_ptr<float>();
        printf("\\left[");
        for (int i = 0; i < numSplats; i++) {
            float3 p = 0.5f*(aabb[2*i]+aabb[2*i+1]);
            printf("\\left(%f,%f,%f\\right),", p.x, p.y, p.z);
        }
        printf("\b\\right]\n");
    }
    #endif

    // auto uniques = torch::unique_consecutive(sorted_morton, true, true);
    // int64_t num_unique = std::get<0>(uniques).numel();
    // int64_t max_collision_count = std::get<2>(uniques).amax().item<int64_t>();
    // printf("%d/%d collisions, %d max\n", (int)(numSplats-num_unique), (int)numSplats, (int)max_collision_count);

    // Build tree
    torch::Tensor internal_nodes = torch::empty({numSplats-1, 2}, tensorI32);
    torch::Tensor parent_nodes = torch::empty({numSplats-1}, tensorI32);
    cudaMemset(parent_nodes.data_ptr<int32_t>(), 0xff, (numSplats-1)*sizeof(int32_t));
    CHECK_DEVICE_ERROR(cudaGetLastError());
    fillLbvhInternalNodes<<<((numSplats-1)+block-1)/block, block>>>(
        MAX_NUM_LEVELS,
        (uint2*)tree_ranges.data_ptr<int32_t>(),
        (uint64_t*)sorted_morton.data_ptr<int64_t>(),
        splat_argsort.data_ptr<int32_t>(),
        (int2*)internal_nodes.data_ptr<int32_t>(),
        parent_nodes.data_ptr<int32_t>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
    print_tensor(int, internal_nodes);
    print_tensor(int, parent_nodes);

    // Compute AABB
    torch::Tensor treeAABB = torch::empty({numSplats, 2, 3}, tensorF32);
    fillTreeSubcells_initAABB<<<((numSplats-1)+block-1)/block, block>>>(
        numSplats-1,
        (float3*)treeAABB.data_ptr<float>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
    computeLbvhAABB<Primitive><<<((numSplats-1)+block-1)/block, block>>>(
        numSplats, splats, MAX_NUM_LEVELS,
        (uint2*)tree_ranges.data_ptr<int32_t>(),
        (int2*)internal_nodes.data_ptr<int32_t>(),
        parent_nodes.data_ptr<int32_t>(),
        (float3*)treeAABB.data_ptr<float>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
    #ifdef DEBUG
    if (0) {
        torch::Tensor treeAABB_cpu = treeAABB.cpu();
        for (int i = 0; i < numSplats-1; i++) {
            float3* p = (float3*)treeAABB_cpu.data_ptr<float>() + 2*i;
            printAABB_wireframe(p[0], p[1]);
        }
    }
    #endif

    // Traverse to find intersections
    torch::Tensor intersection_count = torch::zeros({tiles.size+1}, tensorI32);
    // cudaDeviceSynchronize();
    // getTileSplatIntersections_lbvh<<<(tiles.size+warp-1)/warp, warp>>>(
    getTileSplatIntersections_lbvh_warp<Primitive, camera_model><<<tiles.size, warp>>>(
        tiles, splats, MAX_NUM_LEVELS,
        (uint2*)tree_ranges.data_ptr<int32_t>(),
        (int2*)internal_nodes.data_ptr<int32_t>(),
        (float3*)treeAABB.data_ptr<float>(),
        (uint32_t*)intersection_count.data_ptr<int32_t>(),
        nullptr
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
    print_tensor(int, intersection_count);

    torch::Tensor intersection_count_map = exclusiveScan(intersection_count);
    print_tensor(int, intersection_count_map);
    unsigned total_intersections = (unsigned)intersection_count_map[(long)tiles.size].item<int32_t>();

    torch::Tensor intersectionSplatID = torch::empty({total_intersections}, tensorI32);
    // getTileSplatIntersections_lbvh<<<(tiles.size+warp-1)/warp, warp>>>(
    getTileSplatIntersections_lbvh_warp<Primitive, camera_model><<<tiles.size, warp>>>(
        tiles, splats, MAX_NUM_LEVELS,
        (uint2*)tree_ranges.data_ptr<int32_t>(),
        (int2*)internal_nodes.data_ptr<int32_t>(),
        (float3*)treeAABB.data_ptr<float>(),
        (uint32_t*)intersection_count_map.data_ptr<int32_t>(),
        (uint32_t*)intersectionSplatID.data_ptr<int32_t>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
    print_tensor(int32_t, intersectionSplatID);

    cudaDeviceSynchronize();
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
    // auto seed = std::stoll(argv[3]);
    // torch::manual_seed(seed);

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

    torch::Tensor means = loadBinaryToTensor("means_1000000_3.bin");
    torch::Tensor scales = loadBinaryToTensor("scales_1000000_3.bin");
    torch::Tensor opacs = loadBinaryToTensor("opacities_1000000.bin");
    torch::Tensor quats = loadBinaryToTensor("quats_1000000_4.bin");
    torch::Tensor viewmats = loadBinaryToTensor("viewmats_256_4_4.bin");
    torch::Tensor Ks = loadBinaryToTensor("Ks_256_3_3.bin");

    SplatTileIntersector::intersect_splat_tile(
        means, scales, opacs, quats,
        64, 64,
        viewmats, Ks
    );
    return 0;
}

// /usr/local/cuda-12.8/bin/nvcc -I/media/harry/d/gs/spirulae-splat/spirulae_splat/splat/cuda/csrc/glm -I/home/harry/.venv/base/lib/python3.12/site-packages/torch/include -I/home/harry/.venv/base/lib/python3.12/site-packages/torch/include/torch/csrc/api/include -I/usr/local/cuda-12.8/include /media/harry/d/gs/spirulae-splat/spirulae_splat/splat/cuda/csrc/SplatTileIntersector.cu -o ./temp -D__CUDA_NO_HALF_OPERATORS__ -D__CUDA_NO_HALF_CONVERSIONS__ -D__CUDA_NO_BFLOAT16_CONVERSIONS__ -D__CUDA_NO_HALF2_OPERATORS__ -DDEBUG --expt-relaxed-constexpr --compiler-options ''"'"'-fPIC'"'"'' -O3 --use_fast_math -lineinfo --generate-line-info --source-in-ptx --expt-relaxed-constexpr -Xcudafe=--diag_suppress=20012 -Xcudafe=--diag_suppress=550 -DTORCH_API_INCLUDE_EXTENSION_H -gencode=arch=compute_120,code=compute_120 -gencode=arch=compute_120,code=sm_120 -std=c++17 -L/home/harry/.venv/base/lib/python3.12/site-packages/torch/lib -L/usr/local/cuda-12.8/lib64 -L/usr/lib/x86_64-linux-gnu -lc10 -ltorch -ltorch_cpu -lcudart -lc10_cuda -ltorch_cuda


#endif



std::tuple<torch::Tensor, torch::Tensor>
intersect_splat_tile_3dgs(
    Vanilla3DGS::World::TensorTuple splats_tuple,
    unsigned width,
    unsigned height,
    const torch::Tensor& viewmats,
    const torch::Tensor& Ks,
    const gsplat::CameraModelType& camera_model,
    const CameraDistortionCoeffsTensor& dist_coeffs
) {
    Vanilla3DGS::World::Tensor splats_tensor(splats_tuple);

    if (camera_model == gsplat::CameraModelType::PINHOLE) {
        TileBuffers<gsplat::CameraModelType::PINHOLE> tile_buffers =
            {width, height, viewmats, Ks, dist_coeffs};
        return SplatTileIntersector<Vanilla3DGS, gsplat::CameraModelType::PINHOLE>
            (splats_tensor, tile_buffers).getIntersections_lbvh();
    }
    else if (camera_model == gsplat::CameraModelType::FISHEYE) {
        TileBuffers<gsplat::CameraModelType::FISHEYE> tile_buffers =
            {width, height, viewmats, Ks, dist_coeffs};
        return SplatTileIntersector<Vanilla3DGS, gsplat::CameraModelType::FISHEYE>
            (splats_tensor, tile_buffers).getIntersections_lbvh();
    }
    else
        throw std::runtime_error("Unsupported camera model");

}

std::tuple<torch::Tensor, torch::Tensor>
intersect_splat_tile_opaque_triangle(
    OpaqueTriangle::World::TensorTuple splats_tuple,
    unsigned width,
    unsigned height,
    const torch::Tensor& viewmats,
    const torch::Tensor& Ks,
    const gsplat::CameraModelType& camera_model,
    const CameraDistortionCoeffsTensor& dist_coeffs
) {
    OpaqueTriangle::World::Tensor splats_tensor(splats_tuple);

    if (camera_model == gsplat::CameraModelType::PINHOLE) {
        TileBuffers<gsplat::CameraModelType::PINHOLE> tile_buffers =
            {width, height, viewmats, Ks, dist_coeffs};
        return SplatTileIntersector<OpaqueTriangle, gsplat::CameraModelType::PINHOLE>
            (splats_tensor, tile_buffers).getIntersections_lbvh();
    }
    else if (camera_model == gsplat::CameraModelType::FISHEYE) {
        TileBuffers<gsplat::CameraModelType::FISHEYE> tile_buffers =
            {width, height, viewmats, Ks, dist_coeffs};
        return SplatTileIntersector<OpaqueTriangle, gsplat::CameraModelType::FISHEYE>
            (splats_tensor, tile_buffers).getIntersections_lbvh();
    }
    else
        throw std::runtime_error("Unsupported camera model");

}
