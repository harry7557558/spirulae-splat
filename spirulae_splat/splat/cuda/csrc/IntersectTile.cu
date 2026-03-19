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


template<bool is_counting_pass>
__global__ void intersect_tile_kernel(
    const uint32_t I,  // or 1 in packed mode
    const uint32_t N,  // or nnz in packed mode
    const int32_t *__restrict__ image_ids,  // [nnz], packed mode only
    const float4 *__restrict__ aabb_buffer,  // [..., N, 4], int32, xyxy in pixels
    const float *__restrict__ depths_buffer,  // [..., N]
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

    uint2 tile_min, tile_max;
    tile_min.x = (uint32_t)min(max(0, (int)floorf((xmin + 0.5f) / TILE_SIZE)), (int)tile_width);
    tile_min.y = (uint32_t)min(max(0, (int)floorf((ymin + 0.5f) / TILE_SIZE)), (int)tile_height);
    tile_max.x = (uint32_t)min(max(0, (int)ceilf((xmax + 0.5f) / TILE_SIZE)), (int)tile_width);
    tile_max.y = (uint32_t)min(max(0, (int)ceilf((ymax + 0.5f) / TILE_SIZE)), (int)tile_height);
    if (is_counting_pass) {
        // counting pass only writes out tiles_per_splat
        tiles_per_splat[idx] = static_cast<int32_t>(
            (tile_max.y - tile_min.y) * (tile_max.x - tile_min.x)
        );
        return;
    }

    int64_t iid = (image_ids ? image_ids[idx] : idx / N);

    int32_t depth_i32 = __float_as_int(depths_buffer[idx]);
    
    int64_t cur_idx = (idx == 0) ? 0 : cum_tiles_per_splat[idx - 1];
    int64_t max_idx = cum_tiles_per_splat[idx];
    for (int32_t i = tile_min.y; i < tile_max.y; ++i) {
        for (int32_t j = tile_min.x; j < tile_max.x; ++j) {
            if (cur_idx >= max_idx) break;
            int64_t tile_id = iid * tile_width * tile_height + i * tile_width + j;
            isect_ids[cur_idx] = (tile_id << 32) | (int64_t)depth_i32;
            flatten_ids[cur_idx] = static_cast<int32_t>(idx);
            ++cur_idx;
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
    float radius = 0.5f * fmaxf((float)(xmax - xmin), (float)(ymax - ymin));
    atomicMax(&radii[idx % N], radius);
}


template <typename SplatPrimitive>
__global__ void intersect_mask_eval3d_kernel(
    const uint32_t n_isects,
    const uint32_t I,
    const uint32_t tile_width,
    const uint32_t tile_height,
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t N,
    const int32_t *__restrict__ tile_offsets,  // [I, tile_height, tile_width]
    const int32_t *__restrict__ flatten_ids,  // [n_isects]
    const typename SplatPrimitive::Screen::Buffer splat_buffer,
    const float *__restrict__ viewmats, // [B, I, 4, 4]
    const float4 *__restrict__ intrins,  // [B, I, 4], fx, fy, cx, cy
    ssplat::CameraModelType camera_model,
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    bool *__restrict__ mask  // [n_isects]
) {
    auto block = cg::this_thread_block();
    int32_t image_id = block.group_index().x;
    int32_t tile_id =
        block.group_index().y * tile_width + block.group_index().z;
    float py = (float)(block.group_index().y * TILE_SIZE) +
        ((float)block.thread_index().y + 0.5f) * (float)TILE_SIZE / (float)SAMPLE_TILE_HEIGHT;
    float px = (float)(block.group_index().z * TILE_SIZE) +
        ((float)block.thread_index().x + 0.5f) * (float)TILE_SIZE / (float)SAMPLE_TILE_WIDTH;

    // Load camera
    viewmats += image_id * 16;
    float4 intrin = intrins[image_id];
    float3x3 R = {
        viewmats[0], viewmats[1], viewmats[2],  // 1st row
        viewmats[4], viewmats[5], viewmats[6],  // 2nd row
        viewmats[8], viewmats[9], viewmats[10],  // 3rd row
    };
    float3 t = { viewmats[3], viewmats[7], viewmats[11] };
    float fx = intrin.x, fy = intrin.y, cx = intrin.z, cy = intrin.w;
    CameraDistortionCoeffs dist_coeffs = dist_coeffs_buffer.load(image_id);

    bool inside = (py < image_height && px < image_width);

    float3 raydir;
    inside &= SlangProjectionUtils::generate_ray(
        {(px-cx)/fx, (py-cy)/fy},
        camera_model == ssplat::CameraModelType::FISHEYE, dist_coeffs,
        &raydir
    );
    float3 ray_o = SlangProjectionUtils::transform_ray_o(R, t);
    float3 ray_d = SlangProjectionUtils::transform_ray_d(R, raydir);
    inside &= (dot(ray_d, ray_d) > 0.0f);

    // have all threads in tile process the same gaussians in batches
    bool done = !inside;
    float T = 1.0f;
    tile_offsets += image_id * tile_height * tile_width;
    int32_t range_start = tile_offsets[tile_id];
    int32_t range_end =
        (image_id == I - 1) && (tile_id == tile_width * tile_height - 1)
            ? n_isects
            : tile_offsets[tile_id + 1];
    static constexpr uint32_t BLOCK_SIZE = SAMPLE_TILE_WIDTH * SAMPLE_TILE_HEIGHT;
    uint32_t num_batches =
        (range_end - range_start + BLOCK_SIZE - 1) / BLOCK_SIZE;

    __shared__ typename SplatPrimitive::Screen splat_batch[BLOCK_SIZE];
    __shared__ bool mask_batch[BLOCK_SIZE];

    uint32_t tr = block.thread_rank();
    for (uint32_t b = 0; b < num_batches; ++b) {
        // resync all threads before beginning next batch
        // end early if entire tile is done
        done |= (T <= kTransmitThreshold);
        if (__syncthreads_count(done) >= BLOCK_SIZE) {
            break;
        }

        // each thread fetch 1 gaussian from front to back
        // index of gaussian to load
        uint32_t batch_start = range_start + BLOCK_SIZE * b;
        uint32_t idx = batch_start + tr;
        if (idx < range_end) {
            int32_t g = flatten_ids[idx]; // flatten index in [I * N] or [nnz]
            splat_batch[tr] = SplatPrimitive::Screen::loadWithPrecompute(splat_buffer, g, nullptr);
            mask_batch[tr] = false;
        }
        block.sync();

        // process gaussians in the current batch for this pixel
        uint32_t batch_size = min(BLOCK_SIZE, range_end - batch_start);
        for (uint32_t t = 0; t < batch_size; ++t) {
            typename SplatPrimitive::Screen splat = splat_batch[t];
            float alpha = splat.evaluate_alpha(ray_o, ray_d);
            bool is_visible = inside && (alpha > ALPHA_THRESHOLD) &&
                ((T *= 1.0f - alpha) > kTransmitThreshold);

            // atomic OR is_visible to mask[idx]
            bool is_any_visible = (__ballot_sync(~0u, is_visible) != 0);
            if (is_any_visible && tr % WARP_SIZE == 0)
                mask_batch[t] = true;  // no atomic needed since write same value
        }
        block.sync();

        if (idx < range_end) {
            mask[idx] = mask_batch[tr];
        }
    }
}


template <typename SplatPrimitive>
__global__ void intersect_mask_kernel(
    const uint32_t n_isects,
    const uint32_t I,
    const uint32_t tile_width,
    const uint32_t tile_height,
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t N,
    const int32_t *__restrict__ tile_offsets,  // [I, tile_height, tile_width]
    const int32_t *__restrict__ flatten_ids,  // [n_isects]
    const typename SplatPrimitive::Screen::Buffer splat_buffer,
    bool *__restrict__ mask  // [n_isects]
) {
    auto block = cg::this_thread_block();
    int32_t image_id = block.group_index().x;
    int32_t tile_id =
        block.group_index().y * tile_width + block.group_index().z;
    float py = (float)(block.group_index().y * TILE_SIZE) +
        ((float)block.thread_index().y + 0.5f) * (float)TILE_SIZE / (float)SAMPLE_TILE_HEIGHT;
    float px = (float)(block.group_index().z * TILE_SIZE) +
        ((float)block.thread_index().x + 0.5f) * (float)TILE_SIZE / (float)SAMPLE_TILE_WIDTH;

    bool inside = (py < image_height && px < image_width);

    // have all threads in tile process the same gaussians in batches
    bool done = !inside;
    float T = 1.0f;
    tile_offsets += image_id * tile_height * tile_width;
    int32_t range_start = tile_offsets[tile_id];
    int32_t range_end =
        (image_id == I - 1) && (tile_id == tile_width * tile_height - 1)
            ? n_isects
            : tile_offsets[tile_id + 1];
    static constexpr uint32_t BLOCK_SIZE = SAMPLE_TILE_WIDTH * SAMPLE_TILE_HEIGHT;
    uint32_t num_batches =
        (range_end - range_start + BLOCK_SIZE - 1) / BLOCK_SIZE;

    __shared__ typename SplatPrimitive::Screen splat_batch[BLOCK_SIZE];
    __shared__ bool mask_batch[BLOCK_SIZE];

    uint32_t tr = block.thread_rank();
    for (uint32_t b = 0; b < num_batches; ++b) {
        // resync all threads before beginning next batch
        // end early if entire tile is done
        done |= (T <= kTransmitThreshold);
        if (__syncthreads_count(done) >= BLOCK_SIZE) {
            break;
        }

        // each thread fetch 1 gaussian from front to back
        // index of gaussian to load
        uint32_t batch_start = range_start + BLOCK_SIZE * b;
        uint32_t idx = batch_start + tr;
        if (idx < range_end) {
            int32_t g = flatten_ids[idx]; // flatten index in [I * N] or [nnz]
            splat_batch[tr] = SplatPrimitive::Screen::loadWithPrecompute(splat_buffer, g, nullptr);
            mask_batch[tr] = false;
        }
        block.sync();

        // process gaussians in the current batch for this pixel
        uint32_t batch_size = min(BLOCK_SIZE, range_end - batch_start);
        for (uint32_t t = 0; t < batch_size; ++t) {
            typename SplatPrimitive::Screen splat = splat_batch[t];
            float alpha = splat.evaluate_alpha(px, py);
            bool is_visible = inside && (alpha > ALPHA_THRESHOLD) &&
                ((T *= 1.0f - alpha) > kTransmitThreshold);

            // atomic OR is_visible to mask[idx]
            bool is_any_visible = (__ballot_sync(~0u, is_visible) != 0);
            if (is_any_visible && tr % WARP_SIZE == 0)
                mask_batch[t] = true;  // no atomic needed since write same value
        }
        block.sync();

        if (idx < range_end) {
            mask[idx] = mask_batch[tr];
        }
    }
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
    const uint32_t I,
    const uint32_t image_width,
    const uint32_t image_height,
    std::optional<at::Tensor> image_ids
) {
    DEVICE_GUARD(aabb);
    CHECK_INPUT(aabb);
    CHECK_INPUT(depths);
    if (image_ids.has_value())
        CHECK_INPUT(image_ids.value());

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
    intersect_tile_kernel<true><<<_LAUNCH_ARGS_1D(I*N, 256)>>>(
        packed ? 1 : I,
        N,
        nullptr,  // image_ids
        reinterpret_cast<const float4 *>(aabb.data_ptr<float>()),
        depths.data_ptr<float>(),
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
    int64_t n_isects = cum_tiles_per_splat[-1].item<int64_t>();
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
    set_zero<float>(radii);
    if (n_isects == 0)
        return std::make_tuple(
            isect_ids,
            flatten_ids,
            offsets.zero_(),
            radii
        );
    intersect_tile_kernel<false><<<_LAUNCH_ARGS_1D(I*N, 256)>>>(
        packed ? 1 : I,
        N,
        image_ids.has_value() ? image_ids.value().data_ptr<int32_t>() : nullptr,
        reinterpret_cast<const float4 *>(aabb.data_ptr<float>()),
        depths.data_ptr<float>(),
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

template <typename SplatPrimitive>
std::tuple<
    at::Tensor,  // isect_ids, [n_isects], int64
    at::Tensor,  // flatten_ids, [n_isects], int32
    at::Tensor,  // offsets, [I * n_tiles], int32
    at::Tensor  // radii, [N], float32
> do_intersect_tile_eval3d(
    at::Tensor aabb,  // [..., N, 4], int32, xyxy in pixels
    at::Tensor depths,  // [..., N], float32
    const uint32_t I,
    const uint32_t image_width,
    const uint32_t image_height,
    typename SplatPrimitive::Screen::TensorTuple splats_tuple,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs
) {
    auto [isect_ids, flatten_ids, offsets, radii] = do_intersect_tile_generic(
        aabb,
        depths,
        I,
        image_width,
        image_height,
        std::nullopt
    );
    uint32_t n_isects = isect_ids.numel();

    uint32_t tile_width = _CEIL_DIV(image_width, TILE_SIZE);
    uint32_t tile_height = _CEIL_DIV(image_height, TILE_SIZE);
    uint32_t N = aabb.size(-2);

    /* Compute mask of tiles intersected by AABB but not by the splat itself */
    typename SplatPrimitive::Screen::Tensor splats(splats_tuple);
    at::Tensor mask = at::zeros(
        {n_isects},
        at::TensorOptions().device(aabb.device()).dtype(at::kBool)
    );
    dim3 threads = {SAMPLE_TILE_WIDTH, SAMPLE_TILE_HEIGHT, 1};
    dim3 grid = {I, tile_height, tile_width};
    intersect_mask_eval3d_kernel<SplatPrimitive>
    <<<grid, threads, 0, at::cuda::getCurrentCUDAStream()>>>(
        n_isects,
        I,
        tile_width,
        tile_height,
        image_width,
        image_height,
        N,
        reinterpret_cast<const int32_t *>(offsets.data_ptr<int32_t>()),
        reinterpret_cast<const int32_t *>(flatten_ids.data_ptr<int32_t>()),
        splats.buffer(),
        viewmats.data_ptr<float>(),
        reinterpret_cast<const float4 *>(intrins.data_ptr<float>()),
        cmt(camera_model),
        dist_coeffs,
        mask.data_ptr<bool>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return do_intersect_tile_post(
        isect_ids,
        flatten_ids,
        offsets,
        mask,
        radii,
        I,
        image_width,
        image_height
    );
}

template <typename SplatPrimitive>
std::tuple<
    at::Tensor,  // isect_ids, [n_isects], int64
    at::Tensor,  // flatten_ids, [n_isects], int32
    at::Tensor,  // offsets, [I * n_tiles], int32
    at::Tensor  // radii, [N], float32
> do_intersect_tile(
    at::Tensor aabb,  // [..., N, 4], int32, xyxy in pixels
    at::Tensor depths,  // [..., N], float32
    const uint32_t I,
    const uint32_t image_width,
    const uint32_t image_height,
    typename SplatPrimitive::Screen::TensorTuple splats_tuple
) {
    auto [isect_ids, flatten_ids, offsets, radii] = do_intersect_tile_generic(
        aabb,
        depths,
        I,
        image_width,
        image_height,
        std::nullopt
    );
    uint32_t n_isects = isect_ids.numel();

    uint32_t tile_width = _CEIL_DIV(image_width, TILE_SIZE);
    uint32_t tile_height = _CEIL_DIV(image_height, TILE_SIZE);
    uint32_t N = aabb.size(-2);

    /* Compute mask of tiles intersected by AABB but not by the splat itself */
    typename SplatPrimitive::Screen::Tensor splats(splats_tuple);
    at::Tensor mask = at::zeros(
        {n_isects},
        at::TensorOptions().device(aabb.device()).dtype(at::kBool)
    );
    dim3 threads = {SAMPLE_TILE_WIDTH, SAMPLE_TILE_HEIGHT, 1};
    dim3 grid = {I, tile_height, tile_width};
    intersect_mask_kernel<SplatPrimitive>
    <<<grid, threads, 0, at::cuda::getCurrentCUDAStream()>>>(
        n_isects,
        I,
        tile_width,
        tile_height,
        image_width,
        image_height,
        N,
        reinterpret_cast<const int32_t *>(offsets.data_ptr<int32_t>()),
        reinterpret_cast<const int32_t *>(flatten_ids.data_ptr<int32_t>()),
        splats.buffer(),
        mask.data_ptr<bool>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return do_intersect_tile_post(
        isect_ids,
        flatten_ids,
        offsets,
        mask,
        radii,
        I,
        image_width,
        image_height
    );
}

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    at::Tensor,  // isect_ids, [n_isects], int64
    at::Tensor,  // flatten_ids, [n_isects], int32
    at::Tensor,  // offsets, [I * n_tiles], int32
    at::Tensor  // radii, [N], float32
> intersect_tile_3dgs_tensor(
    at::Tensor aabb,  // [..., N, 4], float, xyxy in pixels
    at::Tensor depths,  // [..., N], float32
    const uint32_t I,
    const uint32_t image_width,
    const uint32_t image_height,
    typename Vanilla3DGS::Screen::TensorTuple splats,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs
) {
    return do_intersect_tile<Vanilla3DGS>(
        aabb,
        depths,
        I,
        image_width,
        image_height,
        splats
    );
}

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    at::Tensor,  // isect_ids, [n_isects], int64
    at::Tensor,  // flatten_ids, [n_isects], int32
    at::Tensor,  // offsets, [I * n_tiles], int32
    at::Tensor  // radii, [N], float32
> intersect_tile_mip_tensor(
    at::Tensor aabb,  // [..., N, 4], float, xyxy in pixels
    at::Tensor depths,  // [..., N], float32
    const uint32_t I,
    const uint32_t image_width,
    const uint32_t image_height,
    typename MipSplatting::Screen::TensorTuple splats,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs
) {
    return do_intersect_tile<MipSplatting>(
        aabb,
        depths,
        I,
        image_width,
        image_height,
        splats
    );
}

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    at::Tensor,  // isect_ids, [n_isects], int64
    at::Tensor,  // flatten_ids, [n_isects], int32
    at::Tensor,  // offsets, [I * n_tiles], int32
    at::Tensor  // radii, [N], float32
> intersect_tile_3dgut_tensor(
    at::Tensor aabb,  // [..., N, 4], float, xyxy in pixels
    at::Tensor depths,  // [..., N], float32
    const uint32_t I,
    const uint32_t image_width,
    const uint32_t image_height,
    typename Vanilla3DGUT::Screen::TensorTuple splats,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs
) {
    return do_intersect_tile_eval3d<Vanilla3DGUT>(
        aabb,
        depths,
        I,
        image_width,
        image_height,
        splats,
        viewmats,
        intrins,
        camera_model,
        dist_coeffs
    );
}

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    at::Tensor,  // isect_ids, [n_isects], int64
    at::Tensor,  // flatten_ids, [n_isects], int32
    at::Tensor,  // offsets, [I * n_tiles], int32
    at::Tensor  // radii, [N], float32
> intersect_tile_3dgut_sv_tensor(
    at::Tensor aabb,  // [..., N, 4], float, xyxy in pixels
    at::Tensor depths,  // [..., N], float32
    const uint32_t I,
    const uint32_t image_width,
    const uint32_t image_height,
    typename SphericalVoronoi3DGUT_Default::Screen::TensorTuple splats,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs
) {
    return do_intersect_tile_eval3d<SphericalVoronoi3DGUT_Default>(
        aabb,
        depths,
        I,
        image_width,
        image_height,
        splats,
        viewmats,
        intrins,
        camera_model,
        dist_coeffs
    );
}

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    at::Tensor,  // isect_ids, [n_isects], int64
    at::Tensor,  // flatten_ids, [n_isects], int32
    at::Tensor,  // offsets, [I * n_tiles], int32
    at::Tensor  // radii, [N], float32
> intersect_tile_opaque_triangle_tensor(
    at::Tensor aabb,  // [..., N, 4], float, xyxy in pixels
    at::Tensor depths,  // [..., N], float32
    const uint32_t I,
    const uint32_t image_width,
    const uint32_t image_height,
    typename OpaqueTriangle::Screen::TensorTuple splats,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs
) {
    return do_intersect_tile_eval3d<OpaqueTriangle>(
        aabb,
        depths,
        I,
        image_width,
        image_height,
        splats,
        viewmats,
        intrins,
        camera_model,
        dist_coeffs
    );
}

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    at::Tensor,  // isect_ids, [n_isects], int64
    at::Tensor,  // flatten_ids, [n_isects], int32
    at::Tensor,  // offsets, [I * n_tiles], int32
    at::Tensor  // radii, [N], float32
> intersect_tile_voxel_tensor(
    at::Tensor aabb,  // [..., N, 4], float, xyxy in pixels
    at::Tensor depths,  // [..., N], float32
    const uint32_t I,
    const uint32_t image_width,
    const uint32_t image_height,
    typename VoxelPrimitive::Screen::TensorTuple splats,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs
) {
    return do_intersect_tile_eval3d<VoxelPrimitive>(
        aabb,
        depths,
        I,
        image_width,
        image_height,
        splats,
        viewmats,
        intrins,
        camera_model,
        dist_coeffs
    );
}
