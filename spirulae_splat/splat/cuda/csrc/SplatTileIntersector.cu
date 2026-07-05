#include "SplatTileIntersector.cuh"

#include <cuda_runtime.h>
#include <cub/cub.cuh>


#include "generated/slang.cuh"
namespace SlangProjectionUtils {
#include "generated/set_namespace.cuh"
#include "generated/projection_utils.cuh"
}



__device__ __forceinline__ float remapFunction(float x, float rel_scale) {
    // identity near origin, proportional to x^(1/k) for large x
    constexpr float k = 2.5f;
    return k * sinhf( (1.0f/k) * asinhf( x / rel_scale ) ) * rel_scale;
}

__device__ __forceinline__ float3 remapAABB(float3 b, float rel_scale) {
    return {
        remapFunction(b.x, rel_scale),
        remapFunction(b.y, rel_scale),
        remapFunction(b.z, rel_scale),
    };
}

template<bool remap>
__device__ bool getAABB(
    const Vanilla3DGS<0>::WorldBuffer& splatBuffer, long idx,
    float3 &aabb_min, float3 &aabb_max, float rel_scale = 1.0f
) {
    float3 mean = splatBuffer.means(idx);
    float3 scales = splatBuffer.scales(idx);
    float opac = splatBuffer.opacities(idx);
    float4 quat = splatBuffer.quats(idx);

    opac = 1.0f / (1.0f+__expf(-opac));
    scales = { __expf(scales.x), __expf(scales.y), __expf(scales.z) };
    quat = normalize(quat);

    float3x3 covar;
    SlangProjectionUtils::quat_scale_to_covar(quat, scales, &covar);

    float extend = fmin(3.33f, sqrt(2.0f * __logf(fmaxf(opac / ALPHA_THRESHOLD, 1.0f))));

    float3 bound = extend * make_float3(
        sqrtf(covar[0].x), sqrtf(covar[1].y), sqrtf(covar[2].z)
    );
    aabb_min = mean - bound;
    aabb_max = mean + bound;

    if (remap) {
        aabb_min = remapAABB(aabb_min, rel_scale);
        aabb_max = remapAABB(aabb_max, rel_scale);
    }
    return opac > ALPHA_THRESHOLD && isfinite(dot(aabb_min, aabb_max));
}


template<CameraModelType camera_model>
struct Tile {
    float3 ro, rd;
    float3 n0, n1, n2, n3;

    __device__ bool init(
        CameraDistortionCoeffs dist_coeffs,
        float x0, float x1, float y0, float y1,
        const float3 R0, const float3 R1, const float3 R2,  // columns of R (world rotation)
        const float3 t                                       // translation
    ) {
        // Camera origin in world: ro = -R * t.  R has columns (R0, R1, R2).
        ro = -(R0 * t.x + R1 * t.y + R2 * t.z);
        // Camera forward (z-axis) in world.
        rd = R2;

        // TODO: better way to handle this in nonlinear and partially invalid case
        // May not matter in training with small tiles, but obvious artifact when rendering >180deg fisheye
        float3 e0_, e1_, e2_, e3_;
        bool valid0 = SlangProjectionUtils::unproject_point({x0, y0}, (int)camera_model, dist_coeffs, &e0_);
        bool valid1 = SlangProjectionUtils::unproject_point({x0, y1}, (int)camera_model, dist_coeffs, &e1_);
        bool valid2 = SlangProjectionUtils::unproject_point({x1, y1}, (int)camera_model, dist_coeffs, &e2_);
        bool valid3 = SlangProjectionUtils::unproject_point({x1, y0}, (int)camera_model, dist_coeffs, &e3_);
        if (!valid0 && valid3 && valid1 && valid2)
            e0_ = e3_ + e1_ - e2_, valid0 = true;
        if (!valid1 && valid0 && valid2 && valid3)
            e1_ = e0_ + e2_ - e3_, valid1 = true;
        if (!valid2 && valid1 && valid3 && valid0)
            e2_ = e1_ + e3_ - e0_, valid2 = true;
        if (!valid3 && valid2 && valid0 && valid1)
            e3_ = e2_ + e0_ - e1_, valid3 = true;
        // Camera-space rays -> world-space rays: e = R * e_cam.
        float3 e0 = R0 * e0_.x + R1 * e0_.y + R2 * e0_.z;
        float3 e1 = R0 * e1_.x + R1 * e1_.y + R2 * e1_.z;
        float3 e2 = R0 * e2_.x + R1 * e2_.y + R2 * e2_.z;
        float3 e3 = R0 * e3_.x + R1 * e3_.y + R2 * e3_.z;

        n0 = normalize(cross(e0, e1));
        n1 = normalize(cross(e1, e2));
        n2 = normalize(cross(e2, e3));
        n3 = normalize(cross(e3, e0));
        return (valid0 && valid1 && valid2 && valid3 &&
            isfinite(length(n0+n1+n2+n3)));
    }

    __device__ __forceinline__ bool isOverlap(float3 aabb_min, float3 aabb_max) const {

        float3 c = 0.5f*(aabb_min+aabb_max);
        float3 r = 0.5f*(aabb_max-aabb_min);

        // intersection test using separating axis theorem
        // has false positive; TODO: tighter one may help performance, most latency is global memory load
        float3 roc = c - ro;
        float s0 = dot(n0, roc) - dot(r, fabs(n0));
        float s1 = dot(n1, roc) - dot(r, fabs(n1));
        float s2 = dot(n2, roc) - dot(r, fabs(n2));
        float s3 = dot(n3, roc) - dot(r, fabs(n3));
        float s = fmaxf(fmaxf(s0, s1), fmaxf(s2, s3));
        if (camera_model != CameraModelType::PINHOLE)
            return s < 0.0f;
        float sz = -dot(rd, roc) - dot(r, fabs(rd));
        return fmaxf(s, sz) < 0.0f;
    }

    // return negative if no overlap, strictly positive for sorting ID
    __device__ __forceinline__ float isOverlap(const Vanilla3DGS<0>::WorldBuffer& splatBuffer, long idx) const {
        // TODO: primitive aware version with less false positives
        float3 aabb_min, aabb_max;
        bool valid_aabb = getAABB<false>(splatBuffer, idx, aabb_min, aabb_max);
        if (!valid_aabb || !isOverlap(aabb_min, aabb_max))
            return -1.0f;
        float3 mean = 0.5f * (aabb_min + aabb_max);
        return camera_model != CameraModelType::PINHOLE ?
            length(mean - ro) :
            dot(mean - ro, rd);  // negative if center is behind
    }

};


template<CameraModelType camera_model>
__device__ __forceinline__ Tile<camera_model>
loadTile(unsigned tileIdx, const TileBuffers<camera_model> buffers, bool& isActive) {
    // Row-major 4x4 view matrix: m[4*row + col].  Top-left 3x3 is the world rotation R,
    // right column is translation t.  Columns of R are read across the row stride.
    const float* m = buffers.viewmats + 16 * tileIdx;
    const float3 R0 = { m[ 0], m[ 4], m[ 8] };
    const float3 R1 = { m[ 1], m[ 5], m[ 9] };
    const float3 R2 = { m[ 2], m[ 6], m[10] };
    const float3 t  = { m[ 3], m[ 7], m[11] };

    float4 intrin = buffers.intrins[tileIdx];  // fx, fy, cx, cy
    float fx = intrin.x;
    float fy = intrin.y;
    float cx = intrin.z;
    float cy = intrin.w;

    Tile<camera_model> res;
    isActive &= res.init(
        buffers.dist_coeffs.load(tileIdx),
        -cx / fx, (buffers.width - cx) / fx,
        -cy / fy, (buffers.height - cy) / fy,
        R0, R1, R2, t
    );
    return res;
}


template<typename Primitive>
__global__ void computeSplatAABB(
    long numSplats,
    const typename Primitive::WorldBuffer splatBuffer,
    float rel_scale,
    float3* __restrict__ aabb,
    float3* __restrict__ aabb_reduced
) {
    unsigned gid = blockIdx.x * blockDim.x + threadIdx.x;
    if (gid >= numSplats)
        return;

    float3 aabb_min, aabb_max;
    bool valid_aabb = getAABB<true>(splatBuffer, gid, aabb_min, aabb_max, rel_scale);
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


__global__ void fillIota(int32_t* __restrict__ out, int64_t n) {
    int64_t idx = (int64_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    out[idx] = (int32_t)idx;
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


template<typename Primitive, CameraModelType camera_model>
__global__ void getTileSplatIntersections_brute(
    const long numSplats,
    const TileBuffers<camera_model> tiles,
    const typename Primitive::WorldBuffer splatBuffer,
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
    const typename Primitive::WorldBuffer splatBuffer,
    float3 root_min, float3 root_max,
    unsigned num_levels, float branch_factor,
    float rel_scale,
    uint64_t* __restrict__ splat_keys
) {
    unsigned tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= numSplats) return;

    float3 aabb_min, aabb_max;
    bool valid_aabb = getAABB<true>(splatBuffer, tid, aabb_min, aabb_max, rel_scale);
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
    const typename Primitive::WorldBuffer splatBuffer,
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


template<typename Primitive, CameraModelType camera_model>
__global__ void getTileSplatIntersections_lbvh_warp(
    const TileBuffers<camera_model> tiles,
    const typename Primitive::WorldBuffer splatBuffer,
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


__forceinline__ DeviceVector<int32_t> exclusiveScan(
    PoolSlot key, const DeviceVector<int32_t>& input
) {
    DeviceVector<int32_t> result;
    result.resize(key, input.size());
    CUB_WRAPPER(cub::DeviceScan::ExclusiveSum,
        (unsigned*)input.data_ptr(),
        (unsigned*)result.data_ptr(),
        (int)input.size());
    return result;
}



template<typename Primitive, CameraModelType camera_model>
SplatTileIntersector<Primitive, camera_model>::SplatTileIntersector(
    const typename Primitive::WorldBuffer &splats,
    const TileBuffers<camera_model> &tiles,
    float rel_scale
) : tiles(tiles), splats(splats), rel_scale(rel_scale)
{
    this->numSplats = splats.size();
}

template<typename Primitive, CameraModelType camera_model>
std::tuple<DeviceVector<int32_t>, DeviceVector<int32_t>>
SplatTileIntersector<Primitive, camera_model>::getIntersections_brute() {
    constexpr unsigned warp = 32;

    DeviceVector<int32_t> intersection_count;
    intersection_count.resize(PoolSlot::StiBruteCount, tiles.size+1);
    intersection_count.zero();
    getTileSplatIntersections_brute<Primitive, camera_model><<<_LAUNCH_ARGS_1D(tiles.size, warp)>>>(
        numSplats,
        tiles, splats,
        (uint32_t*)intersection_count.data_ptr(),
        nullptr
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    DeviceVector<int32_t> intersection_count_map = exclusiveScan(PoolSlot::StiBruteCountMap, intersection_count);
    int32_t total_intersections = 0;
    cudaMemcpy(&total_intersections,
        intersection_count_map.data_ptr() + tiles.size,
        sizeof(int32_t), cudaMemcpyDeviceToHost);

    DeviceVector<int32_t> intersectionSplatID;
    intersectionSplatID.resize(PoolSlot::StiBruteIds, (int64_t)total_intersections);
    getTileSplatIntersections_brute<Primitive, camera_model><<<_LAUNCH_ARGS_1D(tiles.size, warp)>>>(
        numSplats,
        tiles, splats,
        (uint32_t*)intersection_count_map.data_ptr(),
        (uint32_t*)intersectionSplatID.data_ptr()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return std::make_tuple(intersection_count_map, intersectionSplatID);
}


template<typename Primitive, CameraModelType camera_model>
std::tuple<DeviceVector<int32_t>, DeviceVector<int32_t>>
SplatTileIntersector<Primitive, camera_model>::getIntersections_lbvh() {
    // TODO: use a separate rotated AABB aligned with (1,1,1) for thin off-diagnoal Gaussians?
    constexpr uint MAX_NUM_LEVELS = 28;
    constexpr float BRANCH_FACTOR = 2.0f;
    static_assert(MAX_NUM_LEVELS < 32);

    constexpr unsigned block = 256;
    constexpr unsigned warp = 32;
    constexpr int kFloatPInfByte = 0x7f;  // 0x7f7f7f7f -> 3.39615e+38
    constexpr int kFloatNInfByte = 0xfe;  // 0xfefefefe -> -1.69474e+38

    // find splat AABB
    DeviceVector<float3> splat_aabb;
    splat_aabb.resize(PoolSlot::StiLbvhSplatAabb, numSplats * 2);
    DeviceVector<float3> root_aabb;
    root_aabb.resize(PoolSlot::StiLbvhRootAabb, 2);
    cudaMemset((float*)root_aabb.data_ptr() + 0, kFloatPInfByte, 3*sizeof(float));
    cudaMemset((float*)root_aabb.data_ptr() + 3, kFloatNInfByte, 3*sizeof(float));
    computeSplatAABB<Primitive><<<_LAUNCH_ARGS_1D(numSplats, block)>>>(
        numSplats, splats, rel_scale,
        splat_aabb.data_ptr(),
        root_aabb.data_ptr()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    // find root AABB on host, pad to a cube
    float3 rootAABBMin, rootAABBMax;
    {
        float root_aabb_host[6];
        cudaMemcpy(root_aabb_host, root_aabb.data_ptr(), 6*sizeof(float), cudaMemcpyDeviceToHost);
        rootAABBMin = make_float3(root_aabb_host[0], root_aabb_host[1], root_aabb_host[2]);
        rootAABBMax = make_float3(root_aabb_host[3], root_aabb_host[4], root_aabb_host[5]);
        float3 center = 0.5f * (rootAABBMax + rootAABBMin);
        float3 extend = 0.5f * (rootAABBMax - rootAABBMin);
        float max_size = 1.01f * fmaxf(extend.x, fmaxf(extend.y, extend.z));
        float3 max_size3 = {max_size, max_size, max_size};
        rootAABBMin = center - max_size3;
        rootAABBMax = center + max_size3;
    }

    // compute sorting keys (level and Morton code)
    DeviceVector<int64_t> morton;
    morton.resize(PoolSlot::StiLbvhMorton, numSplats);
    fillSplatSortingKeys<Primitive><<<_LAUNCH_ARGS_1D(numSplats, block)>>>(
        numSplats, splats,
        rootAABBMin, rootAABBMax, MAX_NUM_LEVELS, BRANCH_FACTOR, rel_scale,
        (uint64_t*)morton.data_ptr()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    // Sort morton keys with splat indices (iota) as values via CUB radix sort.
    DeviceVector<int32_t> splat_argsort_in;
    splat_argsort_in.resize(PoolSlot::StiLbvhArgsortIn, numSplats);
    fillIota<<<_LAUNCH_ARGS_1D(numSplats, block)>>>(splat_argsort_in.data_ptr(), numSplats);
    CHECK_DEVICE_ERROR(cudaGetLastError());

    DeviceVector<int64_t> sorted_morton;
    sorted_morton.resize(PoolSlot::StiLbvhSortedMorton, numSplats);
    DeviceVector<int32_t> splat_argsort;
    splat_argsort.resize(PoolSlot::StiLbvhArgsort, numSplats);
    CUB_WRAPPER(cub::DeviceRadixSort::SortPairs,
        (const uint64_t*)morton.data_ptr(),
        (uint64_t*)sorted_morton.data_ptr(),
        splat_argsort_in.data_ptr(),
        splat_argsort.data_ptr(),
        (int)numSplats);
    CHECK_DEVICE_ERROR(cudaGetLastError());

    DeviceVector<int32_t> tree_ranges;
    tree_ranges.resize(PoolSlot::StiLbvhTreeRanges, MAX_NUM_LEVELS * 2);
    cudaMemset(tree_ranges.data_ptr(), 0xff, (2*MAX_NUM_LEVELS)*sizeof(int32_t));
    fillLbvhTreeRanges<<<_LAUNCH_ARGS_1D(numSplats, block)>>>(
        MAX_NUM_LEVELS, numSplats,
        (uint64_t*)sorted_morton.data_ptr(),
        (uint2*)tree_ranges.data_ptr()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
    sanitizeLbvhTreeRanges<<<_LAUNCH_ARGS_1D(MAX_NUM_LEVELS, WARP_SIZE)>>>(
        MAX_NUM_LEVELS,
        (uint2*)tree_ranges.data_ptr()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    // Build tree
    DeviceVector<int32_t> internal_nodes;
    internal_nodes.resize(PoolSlot::StiLbvhInternalNodes, (numSplats-1) * 2);
    DeviceVector<int32_t> parent_nodes;
    parent_nodes.resize(PoolSlot::StiLbvhParentNodes, numSplats-1);
    cudaMemset(parent_nodes.data_ptr(), 0xff, (numSplats-1)*sizeof(int32_t));
    CHECK_DEVICE_ERROR(cudaGetLastError());
    fillLbvhInternalNodes<<<_LAUNCH_ARGS_1D(numSplats-1, block)>>>(
        MAX_NUM_LEVELS,
        (uint2*)tree_ranges.data_ptr(),
        (uint64_t*)sorted_morton.data_ptr(),
        splat_argsort.data_ptr(),
        (int2*)internal_nodes.data_ptr(),
        parent_nodes.data_ptr()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    // Compute AABB
    DeviceVector<float3> treeAABB;
    treeAABB.resize(PoolSlot::StiLbvhTreeAabb, numSplats * 2);
    fillTreeSubcells_initAABB<<<_LAUNCH_ARGS_1D(numSplats-1, block)>>>(
        numSplats-1,
        treeAABB.data_ptr()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
    computeLbvhAABB<Primitive><<<_LAUNCH_ARGS_1D(numSplats-1, block)>>>(
        numSplats, splats, MAX_NUM_LEVELS,
        (uint2*)tree_ranges.data_ptr(),
        (int2*)internal_nodes.data_ptr(),
        parent_nodes.data_ptr(),
        treeAABB.data_ptr()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    // Traverse to find intersections
    DeviceVector<int32_t> intersection_count;
    intersection_count.resize(PoolSlot::StiLbvhCount, tiles.size+1);
    intersection_count.zero();
    getTileSplatIntersections_lbvh_warp<Primitive, camera_model><<<tiles.size, warp>>>(
        tiles, splats, MAX_NUM_LEVELS,
        (uint2*)tree_ranges.data_ptr(),
        (int2*)internal_nodes.data_ptr(),
        treeAABB.data_ptr(),
        (uint32_t*)intersection_count.data_ptr(),
        nullptr
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    DeviceVector<int32_t> intersection_count_map = exclusiveScan(PoolSlot::StiLbvhCountMap, intersection_count);
    int32_t total_intersections = 0;
    cudaMemcpy(&total_intersections,
        intersection_count_map.data_ptr() + tiles.size,
        sizeof(int32_t), cudaMemcpyDeviceToHost);

    DeviceVector<int32_t> intersectionSplatID;
    intersectionSplatID.resize(PoolSlot::StiLbvhIds, (int64_t)total_intersections);
    getTileSplatIntersections_lbvh_warp<Primitive, camera_model><<<tiles.size, warp>>>(
        tiles, splats, MAX_NUM_LEVELS,
        (uint2*)tree_ranges.data_ptr(),
        (int2*)internal_nodes.data_ptr(),
        treeAABB.data_ptr(),
        (uint32_t*)intersection_count_map.data_ptr(),
        (uint32_t*)intersectionSplatID.data_ptr()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return std::make_tuple(intersection_count_map, intersectionSplatID);
}



std::tuple<DeviceVector<int32_t>, DeviceVector<int32_t>>
intersect_splat_tile_3dgs(
    std::vector<DeviceTensorFloatND> splats_tuple,
    unsigned width,
    unsigned height,
    DeviceVector<float> viewmats,        // [C*16] row-major 4x4
    DeviceVector<float4> intrins,        // [C]
    const std::string& camera_model,
    const DeviceTensor2D<float>& dist_coeffs,  // [C, 10] or null
    float rel_scale
) {
    Vanilla3DGS<0>::WorldBuffer splats_tensor(splats_tuple);

    const float*  viewmats_ptr    = viewmats.data_ptr();
    const float4* intrins_ptr     = intrins.data_ptr();
    float*        dist_coeffs_ptr = dist_coeffs.data_ptr();
    const long    num_cams        = intrins.size();

    if (cmt(camera_model) == CameraModelType::PINHOLE) {
        TileBuffers<CameraModelType::PINHOLE> tile_buffers(
            width, height, viewmats_ptr, intrins_ptr, dist_coeffs_ptr, num_cams);
        return SplatTileIntersector<Vanilla3DGS<0>, CameraModelType::PINHOLE>
            (splats_tensor, tile_buffers, rel_scale).getIntersections_lbvh();
    }
    else if (cmt(camera_model) == CameraModelType::FISHEYE) {
        TileBuffers<CameraModelType::FISHEYE> tile_buffers(
            width, height, viewmats_ptr, intrins_ptr, dist_coeffs_ptr, num_cams);
        return SplatTileIntersector<Vanilla3DGS<0>, CameraModelType::FISHEYE>
            (splats_tensor, tile_buffers, rel_scale).getIntersections_lbvh();
    }
    else if (cmt(camera_model) == CameraModelType::EQUISOLID) {
        TileBuffers<CameraModelType::EQUISOLID> tile_buffers(
            width, height, viewmats_ptr, intrins_ptr, dist_coeffs_ptr, num_cams);
        return SplatTileIntersector<Vanilla3DGS<0>, CameraModelType::EQUISOLID>
            (splats_tensor, tile_buffers, rel_scale).getIntersections_lbvh();
    }
    else if (cmt(camera_model) == CameraModelType::EQUIRECTANGULAR) {
        TileBuffers<CameraModelType::EQUIRECTANGULAR> tile_buffers(
            width, height, viewmats_ptr, intrins_ptr, dist_coeffs_ptr, num_cams);
        return SplatTileIntersector<Vanilla3DGS<0>, CameraModelType::EQUIRECTANGULAR>
            (splats_tensor, tile_buffers, rel_scale).getIntersections_lbvh();
    }
    else
        throw std::runtime_error("Unsupported camera model");
}
