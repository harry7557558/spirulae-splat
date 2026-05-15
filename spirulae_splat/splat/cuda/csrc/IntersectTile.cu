#include "IntersectTile.cuh"

#include "generated/slang.cuh"
namespace SlangProjectionUtils {
#include "generated/set_namespace.cuh"
#include "generated/projection_utils.cuh"
}

#include "common.cuh"

#include <gsplat/Common.h>
#include <gsplat/Utils.cuh>

#include <cooperative_groups.h>

#include <ATen/ops/empty_like.h>
#include <ATen/ops/empty.h>
#include <ATen/ops/zeros.h>
#include <ATen/ops/cumsum.h>

#include <c10/cuda/CUDACachingAllocator.h>
#include <cub/cub.cuh>


namespace cg = cooperative_groups;


// inline constexpr int SAMPLE_TILE_WIDTH = 8;
// inline constexpr int SAMPLE_TILE_HEIGHT = 4;
inline constexpr int SAMPLE_TILE_WIDTH = TILE_SIZE;
inline constexpr int SAMPLE_TILE_HEIGHT = TILE_SIZE;

inline constexpr float kTransmitThreshold = 1e-4f;


// https://github.com/harry7557558/vksplat/blob/main/vksplat/slang/utils.slang

template<typename T>
inline __device__ T clamp(T x, T a, T b)
    { return x < a ? a : x > b ? b : x; }

inline __device__ float2 ellipse_range_bound(
    float3 inv_cov, float y0, float y1
) {
    // find x_min and x_max of an ellipse clipped to y0 < y <= y1 
    // ellipse: centered at origin, defined by inv_cov and r=1

    float a = inv_cov.x, b = inv_cov.y, c = inv_cov.z;
    float ym = -b/c * sqrtf(c/(a*c-b*b));  // can be pre-computed

    float v0 = clamp(-ym, y0, y1);
    float v1 = clamp(ym, y0, y1);

    float bv0 = -b*v0, bv1 = -b*v1;

    float inv_a = 1.0f / a;
    float x0 = inv_a * (bv0 - sqrtf(bv0*bv0 - a*(c*v0*v0-1.f)));
    float x1 = inv_a * (bv1 + sqrtf(bv1*bv1 - a*(c*v1*v1-1.f)));

    return {x0, x1};
}

inline __device__ int count_ellipse_grid_overlaps(
    float2 xy,
    float3 inv_cov,
    int grid_xmin, int grid_xmax,
    int grid_ymin, int grid_ymax
) {
    // count the number grid cells that overlap with an ellipse

    int n_tiles = 0;

    if (grid_ymax-grid_ymin <= grid_xmax-grid_xmin) {
        for (int y = grid_ymin; y < grid_ymax; y++) {
            float y0 = y * TILE_SIZE - xy.y;
            float y1 = y0 + TILE_SIZE;
            float2 bound = ellipse_range_bound(inv_cov, y0, y1);
            int x0 = int(floor((bound.x + xy.x) / TILE_SIZE));
            int x1 = int(ceil((bound.y + xy.x) / TILE_SIZE));
            x0 = clamp(x0, grid_xmin, grid_xmax);
            x1 = clamp(x1, grid_xmin, grid_xmax);
            n_tiles += x1-x0;
        }
    } else {
        inv_cov = {inv_cov.z, inv_cov.y, inv_cov.x};
        for (int x = grid_xmin; x < grid_xmax; x++) {
            float x0 = x * TILE_SIZE - xy.x;
            float x1 = x0 + TILE_SIZE;
            float2 bound = ellipse_range_bound(inv_cov, x0, x1);
            int y0 = int(floor((bound.x + xy.y) / TILE_SIZE));
            int y1 = int(ceil((bound.y + xy.y) / TILE_SIZE));
            y0 = clamp(y0, grid_ymin, grid_ymax);
            y1 = clamp(y1, grid_ymin, grid_ymax);
            n_tiles += y1-y0;
        }
    }

    return n_tiles;
}


template<bool is_counting_pass, bool is_ellipse>
__global__ void intersect_tile_kernel(
    const uint32_t I,  // or 1 in packed mode
    const uint32_t N,  // or nnz in packed mode
    const int32_t *__restrict__ image_ids,  // [nnz], packed mode only
    const float4* __restrict__ intrins,
    const float4 *__restrict__ aabb_buffer,  // [..., N, 4], int32, xyxy in pixels
    const float *__restrict__ depths_buffer,  // [..., N]
    const float2 *__restrict__ proj_xy,  // [..., N, 2], float, optional
    const float3 *__restrict__ proj_conic,  // [..., N, 3], float, optional
    const float *__restrict__ proj_opac,  // [..., N, 1], float, optional
    const int64_t *__restrict__ cum_tiles_per_splat, // [..., N], optional for counting pass
    const uint32_t tile_width,
    const uint32_t tile_height,
    int64_t *__restrict__ tiles_per_splat, // [..., N]
    int64_t *__restrict__ isect_ids,  // [n_isects]
    int32_t *__restrict__ flatten_ids,  // [n_isects]
    float *__restrict__ radii  // [N]
) {
    uint32_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= I * N) {
        return;
    }

    float4 aabb = aabb_buffer[idx];
    float xmin = aabb.x, ymin = aabb.y;
    float xmax = aabb.z, ymax = aabb.w;

    if (xmax <= xmin || ymax <= ymin) {
        if (is_counting_pass) {
            tiles_per_splat[idx] = 0;
        }
        return;
    }

    float2 xy = is_ellipse ? (proj_xy ? proj_xy[idx] : float2{0.5f*(xmin+xmax), 0.5f*(ymin+ymax)}) : float2{};
    float3 conic = is_ellipse ? proj_conic[idx] : float3{};
    if constexpr (is_ellipse)
        conic = conic * (0.5f / __logf(proj_opac[idx] / ALPHA_THRESHOLD));

    uint2 tile_min, tile_max;
    tile_min.x = (uint32_t)min(max(0, (int)floorf((xmin + 0.5f) / TILE_SIZE)), (int)tile_width);
    tile_min.y = (uint32_t)min(max(0, (int)floorf((ymin + 0.5f) / TILE_SIZE)), (int)tile_height);
    tile_max.x = (uint32_t)min(max(0, (int)ceilf((xmax + 0.5f) / TILE_SIZE)), (int)tile_width);
    tile_max.y = (uint32_t)min(max(0, (int)ceilf((ymax + 0.5f) / TILE_SIZE)), (int)tile_height);
    if constexpr (is_counting_pass) {
        // counting pass only writes out tiles_per_splat
        if constexpr (is_ellipse)
            tiles_per_splat[idx] = count_ellipse_grid_overlaps(
                xy, conic,
                tile_min.x, tile_max.x, tile_min.y, tile_max.y
            );
        else
            tiles_per_splat[idx] = static_cast<int32_t>(
                (tile_max.y - tile_min.y) * (tile_max.x - tile_min.x)
            );
        return;
    }

    int64_t iid = (image_ids ? image_ids[idx] : idx / N);

    float depth_f32 = depths_buffer[idx];
    depth_f32 = fabsf(depth_f32);
    uint32_t depth_u32 = __float_as_uint(depth_f32);
    if (depth_u32 >> 31)  // negative
        depth_u32 = ~depth_u32;
    else  // positive
        depth_u32 ^= (1u << 31u);
    
    int64_t cur_idx = (idx == 0) ? 0 : cum_tiles_per_splat[idx - 1];
    int64_t max_idx = cum_tiles_per_splat[idx];
    if constexpr (is_ellipse) {
        if (tile_max.y-tile_min.y <= tile_max.x-tile_min.x) {
            for (int y = tile_min.y; y < tile_max.y; y++) {
                float y0 = y * TILE_SIZE - xy.y;
                float y1 = y0 + TILE_SIZE;

                float2 bound = ellipse_range_bound(conic, y0, y1);
                int min_x = clamp(
                    (unsigned)floor((bound.x + xy.x) / TILE_SIZE),
                    tile_min.x, tile_max.x
                );
                int max_x = clamp(
                    (unsigned)ceil((bound.y + xy.x) / TILE_SIZE),
                    tile_min.x, tile_max.x
                );
                for (int x = min_x; x < max_x && cur_idx < max_idx; x++) {
                    int64_t tile_id = iid * tile_width * tile_height + int64_t(y) * tile_width + int64_t(x);
                    isect_ids[cur_idx] = (tile_id << 32) | (int64_t)depth_u32;
                    flatten_ids[cur_idx] = static_cast<int32_t>(idx);
                    ++cur_idx;
                }
            }
        } else {
            conic = {conic.z, conic.y, conic.x};
            for (int x = tile_min.x; x < tile_max.x; x++) {
                float x0 = x * TILE_SIZE - xy.x;
                float x1 = x0 + TILE_SIZE;

                float2 bound = ellipse_range_bound(conic, x0, x1);
                int min_y = clamp(
                    (unsigned)floor((bound.x + xy.y) / TILE_SIZE),
                    tile_min.y, tile_max.y
                );
                int max_y = clamp(
                    (unsigned)ceil((bound.y + xy.y) / TILE_SIZE),
                    tile_min.y, tile_max.y
                );
                for (int y = min_y; y < max_y && cur_idx < max_idx; y++) {
                    int64_t tile_id = iid * tile_width * tile_height + int64_t(y) * tile_width + int64_t(x);
                    isect_ids[cur_idx] = (tile_id << 32) | (int64_t)depth_u32;
                    flatten_ids[cur_idx] = static_cast<int32_t>(idx);
                    ++cur_idx;
                }
            }
        }
    }
    else {
        for (int32_t i = tile_min.y; i < tile_max.y; ++i) {
            for (int32_t j = tile_min.x; j < tile_max.x; ++j) {
                if (cur_idx >= max_idx) break;
                int64_t tile_id = iid * tile_width * tile_height + i * tile_width + j;
                isect_ids[cur_idx] = (tile_id << 32) | (int64_t)depth_u32;
                flatten_ids[cur_idx] = static_cast<int32_t>(idx);
                ++cur_idx;
            }
        }
    }
    // this can happen with floating point optimization, make sure it doesn't introduce invalid ID
    while (cur_idx < max_idx) {
        int64_t tile_id = iid * tile_width * tile_height;
        isect_ids[cur_idx] = (tile_id << 32) | (int64_t)0xffffffff;
        flatten_ids[cur_idx] = static_cast<int32_t>(idx);
        ++cur_idx;
    }

    // save radii
    float4 intrin = intrins[iid];
    // float radius = 0.5f * fmaxf((float)(xmax - xmin) / intrin.x, (float)(ymax - ymin) / intrin.y);
    const float image_width = (float)(tile_width*TILE_SIZE);
    const float image_height = (float)(tile_height*TILE_SIZE);
    float radius = 0.5f * fmaxf((float)(xmax - xmin) / fminf(intrin.x, image_width), (float)(ymax - ymin) / fminf(intrin.y, image_height));
    atomicMax(&radii[idx % N], radius);
}



__global__ void intersect_offset_kernel(
    const uint32_t n_isects,
    const int64_t *__restrict__ isect_ids,
    const uint32_t I,
    const uint32_t tile_width,
    const uint32_t tile_height,
    int32_t *__restrict__ offsets // [I, tile_height, tile_width]
) {
    uint32_t idx = cg::this_grid().thread_rank();
    if (idx >= n_isects)
        return;

    int64_t tile_id = isect_ids[idx] >> 32;

    if (idx == 0) {
        // write out the offsets until the first valid tile (inclusive)
        for (uint32_t i = 0; i < tile_id + 1; ++i)
            offsets[i] = static_cast<int32_t>(idx);
    }
    if (idx == n_isects - 1) {
        // write out the rest of the offsets
        for (uint32_t i = tile_id + 1; i < I * tile_width * tile_height; ++i)
            offsets[i] = static_cast<int32_t>(n_isects);
    }

    if (idx > 0) {
        // visit the current and previous isect_id and check if the (bid, cid,
        // tile_id) tuple changes.
        int64_t tile_id_prev = isect_ids[idx - 1] >> 32;
        if (tile_id_prev == tile_id)
            return;

        // write out the offsets between the previous and current tiles
        for (uint32_t i = tile_id_prev + 1; i < tile_id + 1; ++i)
            offsets[i] = static_cast<int32_t>(idx);
    }
}

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    at::Tensor,  // isect_ids, [n_isects], int64
    at::Tensor,  // flatten_ids, [n_isects], int32
    at::Tensor,  // offsets, [I * n_tiles], int32
    at::Tensor  // radii, [N], float32
> do_intersect_tile_generic(
    at::Tensor aabb,  // [..., N, 4], float32, xyxy in pixels
    at::Tensor depths,  // [..., N], float32
    std::optional<std::tuple<at::optional<at::Tensor>, at::Tensor, at::Tensor>> splats_proj,
    const uint32_t I,
    at::Tensor intrins,
    const uint32_t image_width,
    const uint32_t image_height,
    std::optional<at::Tensor> image_ids
) {
    DEVICE_GUARD(aabb);
    CHECK_INPUT(aabb);
    CHECK_INPUT(depths);
    if (image_ids.has_value())
        CHECK_INPUT(image_ids.value());
    if (splats_proj.has_value()) {
        if (std::get<0>(splats_proj.value()).has_value())
            CHECK_INPUT(std::get<0>(splats_proj.value()).value());
        CHECK_INPUT(std::get<1>(splats_proj.value()));
        CHECK_INPUT(std::get<2>(splats_proj.value()));
    }

    if (aabb.size(-1) != 4)
        AT_ERROR("aabb must be of shape [..., N, 4]");
    bool packed = image_ids.has_value();
    const uint32_t N = aabb.size(-2);
    if (depths.numel() != (packed ? N : I * N))
        AT_ERROR("depths must be of shape [..., N]");

    uint32_t tile_width = _CEIL_DIV(image_width, TILE_SIZE);
    uint32_t tile_height = _CEIL_DIV(image_height, TILE_SIZE);
    uint32_t n_tiles = tile_width * tile_height * I;

    /* For each splat, count tiles intersected by its AABB */
    at::Tensor tiles_per_splat = at::empty(
        {packed ? N : I*N},
        at::TensorOptions().device(aabb.device()).dtype(at::kLong)
    );
    (splats_proj.has_value() ?
        intersect_tile_kernel<true, true> :
        intersect_tile_kernel<true, false>
    )<<<_LAUNCH_ARGS_1D(I*N, 256)>>>(
        packed ? 1 : I,
        N,
        nullptr,  // image_ids
        nullptr,  // intrins
        reinterpret_cast<const float4 *>(aabb.data_ptr<float>()),
        depths.data_ptr<float>(),
        splats_proj.has_value() &&  std::get<0>(splats_proj.value()).has_value() ?
            (float2*)std::get<0>(splats_proj.value()).value().data_ptr<float>() : (float2*)nullptr,
        splats_proj.has_value() ? (float3*)std::get<1>(splats_proj.value()).data_ptr<float>() : (float3*)nullptr,
        splats_proj.has_value() ? std::get<2>(splats_proj.value()).data_ptr<float>() : (float*)nullptr,
        nullptr,  // cum_tiles_per_splat
        tile_width,
        tile_height,
        tiles_per_splat.data_ptr<int64_t>(),
        nullptr,  // isect_ids
        nullptr,  // flatten_ids
        nullptr  // radii
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    /* Prefix sum the count to get offsets */
    at::Tensor cum_tiles_per_splat = at::cumsum(tiles_per_splat, /*dim=*/-1);

    /* For each splat, write out keys (intersected tile ID in higher bits, and depth in lower bits) at the correct offsets */
    int64_t n_isects = cum_tiles_per_splat.numel() > 0 ?
        cum_tiles_per_splat[-1].item<int64_t>() : 0;
    at::Tensor isect_ids = at::empty(
        {n_isects},
        at::TensorOptions().device(aabb.device()).dtype(at::kLong)
    );
    at::Tensor flatten_ids = at::empty(
        {n_isects},
        at::TensorOptions().device(aabb.device()).dtype(at::kInt)
    );
    at::Tensor offsets = at::empty(
        {I, tile_height, tile_width},
        at::TensorOptions().device(aabb.device()).dtype(at::kInt)
    );
    at::Tensor radii = at::empty(
        {N,},
        at::TensorOptions().device(aabb.device()).dtype(at::kFloat)
    );
    set_zero_tensor(radii);
    if (n_isects == 0)
        return std::make_tuple(
            isect_ids,
            flatten_ids,
            offsets.zero_(),
            radii
        );
    (splats_proj.has_value() ?
        intersect_tile_kernel<false, true> :
        intersect_tile_kernel<false, false>
    )<<<_LAUNCH_ARGS_1D(I*N, 256)>>>(
        packed ? 1 : I,
        N,
        image_ids.has_value() ? image_ids.value().data_ptr<int32_t>() : nullptr,
        (float4*)intrins.data_ptr<float>(),
        reinterpret_cast<const float4 *>(aabb.data_ptr<float>()),
        depths.data_ptr<float>(),
        splats_proj.has_value() &&  std::get<0>(splats_proj.value()).has_value() ?
            (float2*)std::get<0>(splats_proj.value()).value().data_ptr<float>() : (float2*)nullptr,
        splats_proj.has_value() ? (float3*)std::get<1>(splats_proj.value()).data_ptr<float>() : (float3*)nullptr,
        splats_proj.has_value() ? std::get<2>(splats_proj.value()).data_ptr<float>() : (float*)nullptr,
        reinterpret_cast<const int64_t *>(cum_tiles_per_splat.data_ptr<int64_t>()),
        tile_width,
        tile_height,
        nullptr,  // tiles_per_splat
        reinterpret_cast<int64_t *>(isect_ids.data_ptr<int64_t>()),
        reinterpret_cast<int32_t *>(flatten_ids.data_ptr<int32_t>()),
        radii.data_ptr<float>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    /* Sort tiles by keys */
    at::Tensor isect_ids_sorted = at::empty_like(isect_ids);
    at::Tensor flatten_ids_sorted = at::empty_like(flatten_ids);
    cub::DoubleBuffer<int64_t> d_keys(
        isect_ids.data_ptr<int64_t>(), isect_ids_sorted.data_ptr<int64_t>()
    );
    cub::DoubleBuffer<int32_t> d_values(
        flatten_ids.data_ptr<int32_t>(), flatten_ids_sorted.data_ptr<int32_t>()
    );
    int tile_n_bits = 0;
    while ((1U << tile_n_bits) <= n_tiles)
        ++tile_n_bits;
    CUB_WRAPPER(
        cub::DeviceRadixSort::SortPairs,
        d_keys,
        d_values,
        n_isects,
        0,
        32 + tile_n_bits,
        at::cuda::getCurrentCUDAStream()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
    switch (d_keys.selector) {
    case 0: // sorted items are stored in isect_ids
        break;
    case 1: // sorted items are stored in isect_ids_sorted
        isect_ids = isect_ids_sorted;
        break;
    }
    switch (d_values.selector) {
    case 0: // sorted items are stored in flatten_ids
        break;
    case 1: // sorted items are stored in flatten_ids_sorted
        flatten_ids = flatten_ids_sorted;
        break;
    }

    /* Compute offsets for each tile */
    intersect_offset_kernel<<<_LAUNCH_ARGS_1D(n_isects, 256)>>>(
        n_isects,
        isect_ids.data_ptr<int64_t>(),
        I, tile_width, tile_height,
        offsets.data_ptr<int32_t>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return std::make_tuple(
        isect_ids,
        flatten_ids,
        offsets,
        radii
    );
}


std::tuple<
    at::Tensor,  // isect_ids, [n_isects], int64
    at::Tensor,  // flatten_ids, [n_isects], int32
    at::Tensor,  // offsets, [I * n_tiles], int32
    at::Tensor  // radii, [N], float32
> do_intersect_tile_post(
    at::Tensor isect_ids,
    at::Tensor flatten_ids,
    at::Tensor offsets,
    at::Tensor mask,
    at::Tensor radii,
    const uint32_t I,
    const uint32_t image_width,
    const uint32_t image_height
) {
    uint32_t tile_width = _CEIL_DIV(image_width, TILE_SIZE);
    uint32_t tile_height = _CEIL_DIV(image_height, TILE_SIZE);

    /* Remove not masked tiles */
    isect_ids = isect_ids.masked_select(mask);
    flatten_ids = flatten_ids.masked_select(mask);
    uint32_t n_isects = isect_ids.numel();

    /* Update offsets for each tile */
    offsets = at::empty(
        {I, tile_height, tile_width},
        at::TensorOptions().device(mask.device()).dtype(at::kInt)
    );
    intersect_offset_kernel<<<_LAUNCH_ARGS_1D(n_isects, 256)>>>(
        n_isects,
        isect_ids.data_ptr<int64_t>(),
        I, tile_width, tile_height,
        offsets.data_ptr<int32_t>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return std::make_tuple(
        isect_ids,
        flatten_ids,
        offsets,
        radii
    );
}
