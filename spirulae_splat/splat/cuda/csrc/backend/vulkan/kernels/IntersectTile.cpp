// Vulkan implementation of the tile-intersection launch API
// (csrc/IntersectTile.cuh). Mirrors do_intersect_tile_generic in
// IntersectTile.cu with the CUB spine swapped for backend:: SortScan:
// count -> inclusive_sum<int64> -> n_isects readback -> key write ->
// sort_pairs<int64,int32> (begin_bit 0, end_bit 32 + tile bits) -> offsets.
// Device work runs slang/vulkan/intersect_tile.slang.

#include <IntersectTile.cuh>

#include "../../common/SortScan.h"
#include "../VulkanInternal.h"
#include "../VulkanPipelines.h"

#include <stdexcept>

namespace {

inline constexpr int TILE_SIZE_IX = TILE_SIZE_X * MACRO_TILE_SIZE_X;
inline constexpr int TILE_SIZE_IY = TILE_SIZE_Y * MACRO_TILE_SIZE_Y;

// Mirror the params structs in slang/vulkan/intersect_tile.slang.
struct IsectCountParams {
    uint64_t aabb, proj_xy, proj_conic, proj_opac, tiles_per_splat;
    uint32_t tile_width, tile_height, total, wgs_per_row;
};
static_assert(sizeof(IsectCountParams) == 5 * 8 + 4 * 4,
              "params layout must match the slang struct");

struct IsectWriteParams {
    uint64_t image_ids, aabb, depths, proj_xy, proj_conic, proj_opac;
    uint64_t cum_tiles, isect_ids, flatten_ids;
    uint32_t N, tile_width, tile_height, total, wgs_per_row;
};
static_assert(sizeof(IsectWriteParams) == 9 * 8 + 5 * 4 + 4 /*pad*/,
              "params layout must match the slang struct");

struct IsectOffsetParams {
    uint64_t isect_ids, offsets;
    uint32_t n_isects, n_offsets, wgs_per_row;
};
static_assert(sizeof(IsectOffsetParams) == 2 * 8 + 3 * 4 + 4 /*pad*/,
              "params layout must match the slang struct");

struct Grid {
    uint32_t per_row, rows;
};
Grid fold_1d(uint64_t threads, uint32_t block) {
    uint32_t wgs = (uint32_t)((threads + block - 1) / block);
    uint32_t per_row = std::min(std::max(wgs, 1u), 65535u);
    return {per_row, (wgs + per_row - 1) / per_row};
}

}  // namespace

/* API definition matching csrc/IntersectTile.cuh */

std::tuple<
    DeviceVector<int64_t>,    // isect_ids [n_isects]
    DeviceVector<int32_t>,    // flatten_ids [n_isects]
    DeviceTensor3D<int32_t>   // offsets [I, tile_h, tile_w]
> do_intersect_tile_generic(
    DeviceTensorFloatND aabb,     // [*N, 4] float32
    DeviceTensorFloatND depths,   // [*N] float32
    DeviceTensorFloatND* proj_xy,    // null for AABB mode
    DeviceTensorFloatND* proj_conic, // non-null enables ellipse mode
    DeviceTensorFloatND* proj_opac,  // null for AABB mode
    const uint32_t I,
    TorchTensorView intrins,      // [I, 4] (unused by the kernels, as in CUDA)
    const uint32_t image_width,
    const uint32_t image_height,
    DeviceVector<int32_t>* image_ids  // null for non-packed
) {
    (void)intrins;
    bool packed = image_ids != nullptr;
    // depths is always [*N] float32 (numel = N or nnz), while aabb is [*N, 4].
    const uint32_t total_count = (uint32_t)depths.numel();
    const uint32_t N = packed ? total_count : total_count / I;

    uint32_t tile_width = _CEIL_DIV(image_width, TILE_SIZE_IX);
    uint32_t tile_height = _CEIL_DIV(image_height, TILE_SIZE_IY);
    uint32_t n_tiles = tile_width * tile_height * I;

    const uint64_t p_xy = proj_xy ? (uint64_t)proj_xy->data_ptr() : 0;
    const uint64_t p_conic = proj_conic ? (uint64_t)proj_conic->data_ptr() : 0;
    const uint64_t p_opac = proj_opac ? (uint64_t)proj_opac->data_ptr() : 0;

    /* Count tiles intersected per splat */
    DeviceVector<int64_t> tiles_per_splat;
    tiles_per_splat.resize(PoolSlot::IsectTilesPerSplat, total_count);
    {
        Grid g = fold_1d(total_count, 256);
        IsectCountParams cp{};
        cp.aabb = (uint64_t)aabb.data_ptr();
        cp.proj_xy = p_xy;
        cp.proj_conic = p_conic;
        cp.proj_opac = p_opac;
        cp.tiles_per_splat = (uint64_t)tiles_per_splat.data_ptr();
        cp.tile_width = tile_width;
        cp.tile_height = tile_height;
        cp.total = total_count;
        cp.wgs_per_row = g.per_row;
        if (!backend::vk::dispatch(backend::kDefaultStream,
                                   "intersect_tile.intersect_tile_count", {},
                                   g.per_row, g.rows, 1, &cp, sizeof(cp)))
            throw std::runtime_error(
                "Vulkan backend: intersect count dispatch failed");
    }

    /* Inclusive prefix sum -> cumulative tile counts */
    DeviceVector<int64_t> cum_tiles_per_splat;
    cum_tiles_per_splat.resize(PoolSlot::IsectCumTiles, total_count);
    backend::inclusive_sum<int64_t>(tiles_per_splat.data_ptr(),
                                    cum_tiles_per_splat.data_ptr(),
                                    total_count);

    /* Read total intersection count from the last element */
    int64_t n_isects = 0;
    if (total_count > 0)
        backend::memcpy_sync(&n_isects,
                             cum_tiles_per_splat.data_ptr() + (total_count - 1),
                             sizeof(int64_t),
                             backend::MemcpyKind::DeviceToHost);

    DeviceTensor3D<int32_t> offsets_out;
    offsets_out.resize(PoolSlot::IsectOffsets, I, tile_height, tile_width);

    if (n_isects == 0) {
        offsets_out.zero();
        return std::make_tuple(DeviceVector<int64_t>{},
                               DeviceVector<int32_t>{}, offsets_out);
    }

    /* Write isect keys into double-buffer pair */
    DeviceVector<int64_t> isect_ids_a, isect_ids_b;
    isect_ids_a.resize(PoolSlot::IsectIdsA, n_isects);
    isect_ids_b.resize(PoolSlot::IsectIdsB, n_isects);
    DeviceVector<int32_t> flatten_ids_a, flatten_ids_b;
    flatten_ids_a.resize(PoolSlot::IsectFlatA, n_isects);
    flatten_ids_b.resize(PoolSlot::IsectFlatB, n_isects);

    {
        Grid g = fold_1d(total_count, 256);
        IsectWriteParams wp{};
        wp.image_ids = packed ? (uint64_t)image_ids->data_ptr() : 0;
        wp.aabb = (uint64_t)aabb.data_ptr();
        wp.depths = (uint64_t)depths.data_ptr();
        wp.proj_xy = p_xy;
        wp.proj_conic = p_conic;
        wp.proj_opac = p_opac;
        wp.cum_tiles = (uint64_t)cum_tiles_per_splat.data_ptr();
        wp.isect_ids = (uint64_t)isect_ids_a.data_ptr();
        wp.flatten_ids = (uint64_t)flatten_ids_a.data_ptr();
        wp.N = N;
        wp.tile_width = tile_width;
        wp.tile_height = tile_height;
        wp.total = total_count;
        wp.wgs_per_row = g.per_row;
        if (!backend::vk::dispatch(backend::kDefaultStream,
                                   "intersect_tile.intersect_tile_write", {},
                                   g.per_row, g.rows, 1, &wp, sizeof(wp)))
            throw std::runtime_error(
                "Vulkan backend: intersect write dispatch failed");
    }

    /* Sort by (tile_id << 32 | depth) key */
    backend::DoubleBuffer<int64_t> d_keys(isect_ids_a.data_ptr(),
                                          isect_ids_b.data_ptr());
    backend::DoubleBuffer<int32_t> d_values(flatten_ids_a.data_ptr(),
                                            flatten_ids_b.data_ptr());
    int tile_n_bits = 0;
    while ((1U << tile_n_bits) <= n_tiles) ++tile_n_bits;
    backend::sort_pairs(d_keys, d_values, n_isects, 0, 32 + tile_n_bits);

    /* Pick whichever buffer the sort left the result in */
    DeviceVector<int64_t> isect_ids_out =
        d_keys.selector ? isect_ids_b : isect_ids_a;
    DeviceVector<int32_t> flatten_ids_out =
        d_values.selector ? flatten_ids_b : flatten_ids_a;

    /* Compute per-tile start offsets */
    {
        Grid g = fold_1d(n_isects, 256);
        IsectOffsetParams op{};
        op.isect_ids = (uint64_t)isect_ids_out.data_ptr();
        op.offsets = (uint64_t)offsets_out.data_ptr();
        op.n_isects = (uint32_t)n_isects;
        op.n_offsets = n_tiles;
        op.wgs_per_row = g.per_row;
        if (!backend::vk::dispatch(backend::kDefaultStream,
                                   "intersect_tile.intersect_offset", {},
                                   g.per_row, g.rows, 1, &op, sizeof(op)))
            throw std::runtime_error(
                "Vulkan backend: intersect offset dispatch failed");
    }

    return std::make_tuple(isect_ids_out, flatten_ids_out, offsets_out);
}
