#include "PerPixelLoss.cuh"
#include "FusedSSIM.cuh"

#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>
namespace cg = cooperative_groups;

#include "generated/slang.cuh"
namespace SlangPerPixelLosses {
#include "generated/set_namespace.cuh"
#include "generated/per_pixel_losses.cuh"
}

#include "common.cuh"
#include <Tensor.h>

template<typename T>
static inline TensorView<T, 4> _make_view4(T* data, long B, long H, long W, long C) {
    TensorView<T, 4> v;
    v.data = data;
    v.shape[0] = B; v.shape[1] = H; v.shape[2] = W; v.shape[3] = C;
    v.strides[0] = H*W*C; v.strides[1] = W*C; v.strides[2] = C; v.strides[3] = 1;
    return v;
}

static inline float* _fptr(const TorchTensorView& tv) { return (float*)std::get<0>(tv); }
static inline float3* _f3ptr(const TorchTensorView& tv) { return (float3*)std::get<0>(tv); }
static inline bool* _bptr(const TorchTensorView& tv) { return (bool*)std::get<0>(tv); }
static inline int64_t* _i64ptr(const TorchTensorView& tv) { return (int64_t*)std::get<0>(tv); }
static inline bool _has(const TorchTensorView& tv) { return std::get<0>(tv) != 0; }

static inline TorchTensorView _pool_alloc_f(const std::string& key, long B, long H, long W, long C) {
    float* p = DevicePool::global().acquire<float>(key, (size_t)(B * H * W * C));
    return TorchTensorView((uint64_t)p, 4, {B, H, W, C});
}
static inline TorchTensorView _pool_alloc_f_zero(const std::string& key, long B, long H, long W, long C) {
    float* p = DevicePool::global().acquire<float>(key, (size_t)(B * H * W * C));
    cudaMemset(p, 0, B * H * W * C * sizeof(float));
    return TorchTensorView((uint64_t)p, 4, {B, H, W, C});
}
static inline TorchTensorView _pool_alloc_b(const std::string& key, long B, long H, long W) {
    bool* p = DevicePool::global().acquire<bool>(key, (size_t)(B * H * W));
    return TorchTensorView((uint64_t)p, 1, {B, H, W, 1});
}

static inline TensorView<float, 4> _tv_view4(const TorchTensorView& tv) {
    const auto& s = std::get<2>(tv);
    return _make_view4(_fptr(tv), s[0], s[1], s[2], s[3]);
}
static inline TensorView<uint8_t, 4> _tv_view4_u8(const TorchTensorView& tv) {
    const auto& s = std::get<2>(tv);
    return _make_view4((uint8_t*)std::get<0>(tv), s[0], s[1], s[2], s[3]);
}

template<typename T, uint size>
inline __device__ FixedArray<T, size> loadFixedArray(const T* p, long idx) {
    FixedArray<T, size> arr;
    #pragma unroll
    for (int i = 0; i < size; i++)
        arr[i] = p[idx * size + i];
    return arr;
}

template<typename T, uint size>
inline __device__ void saveFixedArray(T* __restrict__ p, long idx, const FixedArray<T, size> &arr) {
    #pragma unroll
    for (int i = 0; i < size; i++)
        p[idx * size + i] = arr[i];
}


inline __device__ int atomicAddWhich(RawLossIndex idx) {
    if (
        idx == RawLossIndex::DepthSupX ||
        idx == RawLossIndex::DepthSupY ||
        idx == RawLossIndex::DepthSupXX ||
        idx == RawLossIndex::DepthSupYY ||
        idx == RawLossIndex::DepthSupXY
    ) return 1;
    if (
        idx == RawLossIndex::DepthMaskTotal
    ) return 2;
    return 0;
}

inline __device__ bool isTotalReduce(LossIndex idx) {
    return idx != LossIndex::DepthSup;
}


__global__ void per_pixel_losses_forward_kernel(
    const size_t batch_size,
    const size_t pixels_per_image,
    const int64_t* __restrict__ camera_indices,
    const float3* __restrict__ render_rgb,
    const float3* __restrict__ ref_rgb,
    const float* __restrict__ render_depth,
    const float* __restrict__ ref_depth,
    const float3* __restrict__ render_normal,
    const float3* __restrict__ depth_normal,
    const float3* __restrict__ ref_normal,
    const float* __restrict__ render_Ts,
    const float3* __restrict__ rgb_dist,
    const float* __restrict__ depth_dist,
    const float3* __restrict__ normal_dist,
    const bool* __restrict__ ref_alpha,
    const bool* __restrict__ mask,
    const bool* __restrict__ depth_mask,
    const bool* __restrict__ normal_mask,
    const bool* __restrict__ alpha_mask,
    FixedArray<float, (uint)LossWeightIndex::length> loss_weights,
    float* __restrict__ out_loss_map,  // non differentiable
    float* __restrict__ out_losses
) {
    size_t pixel_idx = blockIdx.x * blockDim.x + threadIdx.x;
    size_t batch_idx = blockIdx.y * blockDim.y + threadIdx.y;
    if (batch_idx >= batch_size)
        return;
    size_t idx = batch_idx * pixels_per_image + pixel_idx;

    FixedArray<float, (uint)RawLossIndex::length> losses;

    bool inside = pixel_idx < pixels_per_image;
    if (inside) {
        SlangPerPixelLosses::per_pixel_losses(
            render_rgb ? render_rgb[idx] : make_float3(0),
            ref_rgb ? ref_rgb[idx] : make_float3(0),
            render_depth ? render_depth[idx] : 1.f,
            ref_depth ? ref_depth[idx] : 1.f,
            render_normal ? render_normal[idx] : make_float3(0),
            depth_normal ? depth_normal[idx] : make_float3(0),
            ref_normal ? ref_normal[idx] : make_float3(0),
            render_Ts ? render_Ts[idx] : 1.f,
            rgb_dist ? rgb_dist[idx] : make_float3(0),
            depth_dist ? depth_dist[idx] : 0.f,
            normal_dist ? normal_dist[idx] : make_float3(0),
            ref_alpha ? ref_alpha[idx] : true,
            mask ? mask[idx] : true,
            depth_mask ? depth_mask[idx] : true,
            normal_mask ? normal_mask[idx] : true,
            alpha_mask ? alpha_mask[idx] : true,
            loss_weights,
            &losses
        );

        if (out_loss_map != nullptr) {
            out_loss_map[idx] = fabsf(
                losses[(int)RawLossIndex::RgbLoss] +
                losses[(int)RawLossIndex::RenderNormalSup] +
                losses[(int)RawLossIndex::DepthNormalSup] +
                losses[(int)RawLossIndex::AlphaSup] +
                losses[(int)RawLossIndex::AlphaSupUnder] +
                losses[(int)RawLossIndex::NormalReg] +
                losses[(int)RawLossIndex::AlphaReg] +
                losses[(int)RawLossIndex::RgbDistReg] +
                losses[(int)RawLossIndex::DepthDistReg] +
                losses[(int)RawLossIndex::NormalDistReg]
            );
            // TODO: more accurate version
        }
    }

    auto block = cg::this_thread_block();
    cg::thread_block_tile<WARP_SIZE> warp = cg::tiled_partition<WARP_SIZE>(block);
    uint warp_idx = block.thread_rank() / WARP_SIZE;

    __shared__ float atomic_reduce[WARP_SIZE];

    if (camera_indices != nullptr)
        batch_idx = camera_indices[batch_idx];
    for (uint i = 0; i < (uint)RawLossIndex::length; i++) {
        float loss = inside ? losses[i] : 0.0f;
        loss = isfinite(loss) ? loss : 0.0f;
        float loss_reduced = cg::reduce(warp, loss, cg::plus<float>());
        if (warp.thread_rank() == 0)
            atomic_reduce[warp_idx] = loss_reduced;
        __syncthreads();
        loss = (warp_idx == 0) ? atomic_reduce[warp.thread_rank()] : 0.0f;
        if (__ballot_sync(~0u, loss != 0.0f) == 0)
            continue;
        loss_reduced = cg::reduce(warp, loss, cg::plus<float>());
        if (block.thread_rank() == 0 && loss_reduced != 0.0f && isfinite(loss_reduced)) {
            int which = atomicAddWhich((RawLossIndex)i);
            if (which == 0 || which == 2)
                atomicAdd(out_losses + i, loss_reduced);
            if (which == 1 || which == 2)
                atomicAdd(out_losses + (batch_idx+1)*(size_t)RawLossIndex::length + i, loss_reduced);
        }
    }
}

__global__ void per_pixel_losses_backward_kernel(
    const size_t batch_size,
    const size_t pixels_per_image,
    const int64_t* __restrict__ camera_indices,
    const float3* __restrict__ render_rgb,
    const float3* __restrict__ ref_rgb,
    const float* __restrict__ render_depth,
    const float* __restrict__ ref_depth,
    const float3* __restrict__ render_normal,
    const float3* __restrict__ depth_normal,
    const float3* __restrict__ ref_normal,
    const float* __restrict__ render_Ts,
    const float3* __restrict__ rgb_dist,
    const float* __restrict__ depth_dist,
    const float3* __restrict__ normal_dist,
    const bool* __restrict__ ref_alpha,
    const bool* __restrict__ mask,
    const bool* __restrict__ depth_mask,
    const bool* __restrict__ normal_mask,
    const bool* __restrict__ alpha_mask,
    FixedArray<float, (uint)LossWeightIndex::length> loss_weights,
    const float* __restrict__ v_out_losses,
    float3* __restrict__ v_render_rgb,
    float3* __restrict__ v_ref_rgb,
    float* __restrict__ v_render_depth,
    float* __restrict__ v_ref_depth,
    float3* __restrict__ v_render_normal,
    float3* __restrict__ v_depth_normal,
    float3* __restrict__ v_ref_normal,
    float* __restrict__ v_render_Ts,
    float3* __restrict__ v_rgb_dist,
    float* __restrict__ v_depth_dist,
    float3* __restrict__ v_normal_dist
) {
    size_t pixel_idx = blockIdx.x * blockDim.x + threadIdx.x;
    size_t batch_idx = blockIdx.y * blockDim.y + threadIdx.y;
    if (batch_idx >= batch_size)
        return;
    size_t idx = batch_idx * pixels_per_image + pixel_idx;

    bool inside = pixel_idx < pixels_per_image;
    if (!inside) return;

    if (camera_indices != nullptr)
        batch_idx = camera_indices[batch_idx];
    FixedArray<float, (uint)RawLossIndex::length> v_losses;
    for (uint i = 0; i < (uint)RawLossIndex::length; i++) {
        float temp_loss = 0.0f;
        int which = atomicAddWhich((RawLossIndex)i);
        if (which == 0 || which == 2)
            temp_loss += v_out_losses[i];
        if (which == 1 || which == 2)
            temp_loss += v_out_losses[(batch_idx+1)*(size_t)RawLossIndex::length + i];
        v_losses[i] = temp_loss;
    }

    float3 temp_v_render_rgb;
    float3 temp_v_ref_rgb;
    float temp_v_render_depth;
    float temp_v_ref_depth;
    float3 temp_v_render_normal;
    float3 temp_v_depth_normal;
    float3 temp_v_ref_normal;
    float temp_v_render_Ts;
    float3 temp_v_rgb_dist;
    float temp_v_depth_dist;
    float3 temp_v_normal_dist;

    SlangPerPixelLosses::per_pixel_losses_bwd(
        render_rgb ? render_rgb[idx] : make_float3(0),
        ref_rgb ? ref_rgb[idx] : make_float3(0),
        render_depth ? render_depth[idx] : 1.f,
        ref_depth ? ref_depth[idx] : 1.f,
        render_normal ? render_normal[idx] : make_float3(0),
        depth_normal ? depth_normal[idx] : make_float3(0),
        ref_normal ? ref_normal[idx] : make_float3(0),
        render_Ts ? render_Ts[idx] : 1.f,
        rgb_dist ? rgb_dist[idx] : make_float3(0),
        depth_dist ? depth_dist[idx] : 0.f,
        normal_dist ? normal_dist[idx] : make_float3(0),
        ref_alpha ? ref_alpha[idx] : true,
        mask ? mask[idx] : true,
        depth_mask ? depth_mask[idx] : true,
        normal_mask ? normal_mask[idx] : true,
        alpha_mask ? alpha_mask[idx] : true,
        loss_weights,
        v_losses,
        &temp_v_render_rgb,
        &temp_v_ref_rgb,
        &temp_v_render_depth,
        &temp_v_ref_depth,
        &temp_v_render_normal,
        &temp_v_depth_normal,
        &temp_v_ref_normal,
        &temp_v_render_Ts,
        &temp_v_rgb_dist,
        &temp_v_depth_dist,
        &temp_v_normal_dist
    );

    if (v_render_rgb) v_render_rgb[idx] = temp_v_render_rgb;
    if (v_ref_rgb) v_ref_rgb[idx] = temp_v_ref_rgb;
    if (v_render_depth) v_render_depth[idx] = temp_v_render_depth;
    if (v_ref_depth) v_ref_depth[idx] = temp_v_ref_depth;
    if (v_render_normal) v_render_normal[idx] = temp_v_render_normal;
    if (v_depth_normal) v_depth_normal[idx] = temp_v_depth_normal;
    if (v_ref_normal) v_ref_normal[idx] = temp_v_ref_normal;
    if (v_render_Ts) v_render_Ts[idx] = temp_v_render_Ts;
    if (v_rgb_dist) v_rgb_dist[idx] = temp_v_rgb_dist;
    if (v_depth_dist) v_depth_dist[idx] = temp_v_depth_dist;
    if (v_normal_dist) v_normal_dist[idx] = temp_v_normal_dist;
}

__global__ void per_pixel_losses_reduce_forward_kernel(
    const size_t batch_size,
    const float* __restrict__ raw_losses,  // [B, RawLossIndex::length]
    FixedArray<float, (uint)LossWeightIndex::length> loss_weights,
    float* __restrict__ losses  // [LossIndex::length]
) {
    size_t batch_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (batch_idx > batch_size)
        return;

    FixedArray<float, (uint)RawLossIndex::length> local_raw_losses =
        loadFixedArray<float, (uint)RawLossIndex::length>(raw_losses, batch_idx);

    FixedArray<float, (uint)LossIndex::length> local_losses;
    SlangPerPixelLosses::per_pixel_losses_reduce(
        local_raw_losses, loss_weights,
        &local_losses
    );
    #pragma unroll
    for (int i = 0; i < (uint)LossIndex::length; i++)
        if (local_losses[i] != 0.0f && isfinite(local_losses[i])) {
            if (isTotalReduce((LossIndex)i)) {
                if (batch_idx == 0)
                    losses[i] = local_losses[i];
            }
            else if (batch_idx != 0)
                atomicAdd(&losses[i], local_losses[i] / (float)batch_size);
        }
}

__global__ void per_pixel_losses_reduce_backward_kernel(
    const size_t batch_size,
    const float* __restrict__ raw_losses,  // [B, RawLossIndex::length]
    FixedArray<float, (uint)LossWeightIndex::length> loss_weights,
    const float* __restrict__ v_losses,  // [LossIndex::length]
    float* __restrict__ v_raw_losses  // [B, RawLossIndex::length]
) {
    size_t batch_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (batch_idx > batch_size)
        return;

    FixedArray<float, (uint)LossIndex::length> local_v_losses;
    #pragma unroll
    for (int i = 0; i < (uint)LossIndex::length; i++) {
        float temp_v_loss = 0.0f;
        if (isTotalReduce((LossIndex)i)) {
            if (batch_idx == 0)
                temp_v_loss = v_losses[i];
        }
        else if (batch_idx != 0)
            temp_v_loss = v_losses[i] / (float)batch_size;
        local_v_losses[i] = temp_v_loss;
    }

    FixedArray<float, (uint)RawLossIndex::length> local_raw_losses =
        loadFixedArray<float, (uint)RawLossIndex::length>(raw_losses, batch_idx);
    FixedArray<float, (uint)RawLossIndex::length> local_v_raw_losses;

    SlangPerPixelLosses::per_pixel_losses_reduce_bwd(
        local_raw_losses, loss_weights,
        local_v_losses, &local_v_raw_losses
    );

    saveFixedArray<float, (uint)RawLossIndex::length>(v_raw_losses, batch_idx, local_v_raw_losses);
}


static void _compute_per_pixel_losses_forward(
    long B, long pixels_per_image,
    TorchTensorView render_rgb,
    TorchTensorView ref_rgb,
    TorchTensorView render_depth,
    TorchTensorView ref_depth,
    TorchTensorView render_normal,
    TorchTensorView depth_normal,
    TorchTensorView ref_normal,
    TorchTensorView render_Ts,
    TorchTensorView rgb_dist,
    TorchTensorView depth_dist,
    TorchTensorView normal_dist,
    TorchTensorView ref_alpha,
    TorchTensorView mask,
    TorchTensorView depth_mask,
    TorchTensorView normal_mask,
    TorchTensorView alpha_mask,
    FixedArray<float, (uint)LossWeightIndex::length> loss_weights,
    long num_train_images,
    TorchTensorView camera_indices,
    float* loss_map_ptr,
    float* raw_losses_ptr,
    float* losses_ptr
) {
    per_pixel_losses_forward_kernel<<<_LAUNCH_ARGS_2D(pixels_per_image, B, WARP_SIZE*WARP_SIZE, 1)>>>(
        B, pixels_per_image,
        _i64ptr(camera_indices),
        _f3ptr(render_rgb),
        _f3ptr(ref_rgb),
        _fptr(render_depth),
        _fptr(ref_depth),
        _f3ptr(render_normal),
        _f3ptr(depth_normal),
        _f3ptr(ref_normal),
        _fptr(render_Ts),
        _f3ptr(rgb_dist),
        _fptr(depth_dist),
        _f3ptr(normal_dist),
        _bptr(ref_alpha),
        _bptr(mask),
        _bptr(depth_mask),
        _bptr(normal_mask),
        _bptr(alpha_mask),
        loss_weights,
        loss_map_ptr,
        raw_losses_ptr
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    per_pixel_losses_reduce_forward_kernel
    <<<_LAUNCH_ARGS_1D(num_train_images+1, WARP_SIZE)>>>(
        num_train_images,
        raw_losses_ptr,
        loss_weights,
        losses_ptr
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


static void _compute_per_pixel_losses_backward(
    long B, long pixels_per_image,
    TorchTensorView render_rgb,
    TorchTensorView ref_rgb,
    TorchTensorView render_depth,
    TorchTensorView ref_depth,
    TorchTensorView render_normal,
    TorchTensorView depth_normal,
    TorchTensorView ref_normal,
    TorchTensorView render_Ts,
    TorchTensorView rgb_dist,
    TorchTensorView depth_dist,
    TorchTensorView normal_dist,
    TorchTensorView ref_alpha,
    TorchTensorView mask,
    TorchTensorView depth_mask,
    TorchTensorView normal_mask,
    TorchTensorView alpha_mask,
    float* raw_losses_ptr,
    FixedArray<float, (uint)LossWeightIndex::length> loss_weights,
    float* v_losses_ptr,
    long num_train_images,
    TorchTensorView camera_indices,
    PerPixelGrads& grads
) {
    float* v_raw_losses = DevicePool::global().acquire<float>(
        "ppl.v_raw_losses", (size_t)(num_train_images+1) * (uint)RawLossIndex::length);

    per_pixel_losses_reduce_backward_kernel<<<_LAUNCH_ARGS_1D(num_train_images+1, WARP_SIZE)>>>(
        num_train_images,
        raw_losses_ptr,
        loss_weights,
        v_losses_ptr,
        v_raw_losses
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    per_pixel_losses_backward_kernel<<<_LAUNCH_ARGS_2D(pixels_per_image, B, WARP_SIZE*WARP_SIZE, 1)>>>(
        B, pixels_per_image,
        _i64ptr(camera_indices),
        _f3ptr(render_rgb),
        _f3ptr(ref_rgb),
        _fptr(render_depth),
        _fptr(ref_depth),
        _f3ptr(render_normal),
        _f3ptr(depth_normal),
        _f3ptr(ref_normal),
        _fptr(render_Ts),
        _f3ptr(rgb_dist),
        _fptr(depth_dist),
        _f3ptr(normal_dist),
        _bptr(ref_alpha),
        _bptr(mask),
        _bptr(depth_mask),
        _bptr(normal_mask),
        _bptr(alpha_mask),
        loss_weights,
        v_raw_losses,
        _f3ptr(grads.v_render_rgb),
        _f3ptr(grads.v_ref_rgb),
        _fptr(grads.v_render_depth),
        _fptr(grads.v_ref_depth),
        _f3ptr(grads.v_render_normal),
        _f3ptr(grads.v_depth_normal),
        _f3ptr(grads.v_ref_normal),
        _fptr(grads.v_render_Ts),
        _f3ptr(grads.v_rgb_dist),
        _fptr(grads.v_depth_dist),
        _f3ptr(grads.v_normal_dist)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


__global__ void avg_pool_downsample_float_kernel(
    const TensorView<float, 4> image_hs,
    TensorView<float, 4> image_ls
) {
    uint32_t xid = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t yid = blockIdx.y * blockDim.y + threadIdx.y;
    uint32_t bid = blockIdx.z * blockDim.z + threadIdx.z;
    if (yid >= image_ls.shape[1] || xid >= image_ls.shape[2])
        return;
    for (int c = 0; c < image_ls.shape[3]; ++c) {
        float v =
            image_hs.at(bid, 2*yid+0, 2*xid+0, c) +
            image_hs.at(bid, 2*yid+0, 2*xid+1, c) +
            image_hs.at(bid, 2*yid+1, 2*xid+0, c) +
            image_hs.at(bid, 2*yid+1, 2*xid+1, c);
        image_ls.at(bid, yid, xid, c) = 0.25f*v;
    }
}

template<typename uintx_t>
__global__ void avg_pool_downsample_integral_kernel(
    const TensorView<uintx_t, 4> image_hs,
    TensorView<uintx_t, 4> image_ls
) {
    uint32_t xid = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t yid = blockIdx.y * blockDim.y + threadIdx.y;
    uint32_t bid = blockIdx.z * blockDim.z + threadIdx.z;
    if (yid >= image_ls.shape[1] || xid >= image_ls.shape[2])
        return;
    for (int c = 0; c < image_ls.shape[3]; ++c) {
        float v =
            (float)image_hs.at(bid, 2*yid+0, 2*xid+0, c) +
            (float)image_hs.at(bid, 2*yid+0, 2*xid+1, c) +
            (float)image_hs.at(bid, 2*yid+1, 2*xid+0, c) +
            (float)image_hs.at(bid, 2*yid+1, 2*xid+1, c);
        image_ls.at(bid, yid, xid, c) = (uintx_t)(0.25f*v + 0.5f);
    }
}

__global__ void avg_pool_downsample_bool_kernel(
    const TensorView<uint8_t, 4> image_hs,
    TensorView<uint8_t, 4> image_ls
) {
    uint32_t xid = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t yid = blockIdx.y * blockDim.y + threadIdx.y;
    uint32_t bid = blockIdx.z * blockDim.z + threadIdx.z;
    if (yid >= image_ls.shape[1] || xid >= image_ls.shape[2])
        return;
    for (int c = 0; c < image_ls.shape[3]; ++c) {
        uint8_t v =
            image_hs.at(bid, 2*yid+0, 2*xid+0, c) +
            image_hs.at(bid, 2*yid+0, 2*xid+1, c) +
            image_hs.at(bid, 2*yid+1, 2*xid+0, c) +
            image_hs.at(bid, 2*yid+1, 2*xid+1, c);
        image_ls.at(bid, yid, xid, c) = (uint8_t)(v >= 2);
    }
}


static void _avg_pool_downsample_float(const TorchTensorView& src, const TorchTensorView& dst) {
    const auto& s = std::get<2>(dst);
    long b = s[0], h = s[1], w = s[2], c = s[3];
    avg_pool_downsample_float_kernel<<<_LAUNCH_ARGS_3D(w, h, b, 16, 16, 1)>>>(
        _tv_view4(src), _tv_view4(dst));
    CHECK_DEVICE_ERROR(cudaGetLastError());
}

static void _avg_pool_downsample_bool(const TorchTensorView& src, const TorchTensorView& dst) {
    const auto& s = std::get<2>(dst);
    long b = s[0], h = s[1], w = s[2];
    avg_pool_downsample_bool_kernel<<<_LAUNCH_ARGS_3D(w, h, b, 16, 16, 1)>>>(
        _tv_view4_u8(src), _tv_view4_u8(dst));
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


__global__ void vector_add_scaled_kernel(
    float* __restrict__ dst,
    const float* __restrict__ src,
    float scale,
    int n
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) dst[i] += src[i] * scale;
}

__global__ void avg_pool_upsample_float_kernel(
    TensorView<float, 4> image_hs,
    const TensorView<float, 4> image_ls,
    int scale,
    float a, float b
) {
    uint32_t xid = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t yid = blockIdx.y * blockDim.y + threadIdx.y;
    uint32_t bid = blockIdx.z * blockDim.z + threadIdx.z;
    if (yid >= image_hs.shape[1] || xid >= image_hs.shape[2])
        return;
    for (int c = 0; c < image_hs.shape[3]; ++c) {
        float v = (a == 0.0f) ? 0.0f :
            a * image_hs.at(bid, yid, xid, c);
        if (yid/scale < image_ls.shape[1] && xid/scale < image_ls.shape[2])
            v += b * image_ls.at(bid, yid/scale, xid/scale, c);
        image_hs.at(bid, yid, xid, c) = v;
    }
}


/*[AutoHeaderGeneratorExport]*/
std::tuple<float, float> compute_multi_scale_per_pixel_losses(
    int num_loss_scales,
    TorchTensorView render_rgb,
    TorchTensorView ref_rgb,
    TorchTensorView render_depth,
    TorchTensorView ref_depth,
    TorchTensorView render_normal,
    TorchTensorView depth_normal,
    TorchTensorView ref_normal,
    TorchTensorView render_Ts,
    TorchTensorView rgb_dist,
    TorchTensorView depth_dist,
    TorchTensorView normal_dist,
    TorchTensorView ref_alpha,
    TorchTensorView mask,
    TorchTensorView depth_mask,
    TorchTensorView normal_mask,
    TorchTensorView alpha_mask,
    const std::array<float, (int)LossWeightIndex::length> loss_weights_0,
    const float w_ssim,
    TorchTensorView v_losses,
    std::vector<bool> needs_input_grad,
    long num_train_images,
    TorchTensorView camera_indices,
    TorchTensorView loss_map_out,
    TorchTensorView total_losses_out,
    PerPixelGrads& grads_out
) {
    const auto& s = std::get<2>(render_rgb);
    long B = s[0], H = s[1], W = s[2];

    FixedArray<float, (uint)LossWeightIndex::length> loss_weights =
        *reinterpret_cast<const FixedArray<float, (uint)LossWeightIndex::length>*>(loss_weights_0.data());

    if (!_has(camera_indices))
        num_train_images = B;

    // Scale arrays: [scale] -> TorchTensorView (null if not used)
    constexpr int MAX_SCALES = 4;
    if (num_loss_scales > MAX_SCALES)
        throw std::runtime_error("num_loss_scales > MAX_SCALES");

    TorchTensorView render_rgb_s[MAX_SCALES], ref_rgb_s[MAX_SCALES];
    TorchTensorView render_depth_s[MAX_SCALES], ref_depth_s[MAX_SCALES];
    TorchTensorView render_normal_s[MAX_SCALES], depth_normal_s[MAX_SCALES], ref_normal_s[MAX_SCALES];
    TorchTensorView render_Ts_s[MAX_SCALES];
    TorchTensorView rgb_dist_s[MAX_SCALES], depth_dist_s[MAX_SCALES], normal_dist_s[MAX_SCALES];
    TorchTensorView ref_alpha_s[MAX_SCALES], mask_s[MAX_SCALES];
    TorchTensorView depth_mask_s[MAX_SCALES], normal_mask_s[MAX_SCALES], alpha_mask_s[MAX_SCALES];

    render_rgb_s[0] = render_rgb; ref_rgb_s[0] = ref_rgb;
    render_depth_s[0] = render_depth; ref_depth_s[0] = ref_depth;
    render_normal_s[0] = render_normal; depth_normal_s[0] = depth_normal; ref_normal_s[0] = ref_normal;
    render_Ts_s[0] = render_Ts;
    rgb_dist_s[0] = rgb_dist; depth_dist_s[0] = depth_dist; normal_dist_s[0] = normal_dist;
    ref_alpha_s[0] = ref_alpha; mask_s[0] = mask;
    depth_mask_s[0] = depth_mask; normal_mask_s[0] = normal_mask; alpha_mask_s[0] = alpha_mask;

    // Downsample to create scales
    for (int sc = 1; sc < num_loss_scales; ++sc) {
        const auto& ps = std::get<2>(render_rgb_s[sc-1]);
        long Hp = ps[1], Wp = ps[2];
        long Hn = Hp / 2, Wn = Wp / 2;
        std::string pfx = "ppl.s" + std::to_string(sc) + ".";

        auto ds_f = [&](TorchTensorView& prev, TorchTensorView& curr, const std::string& name, int C) {
            if (_has(prev)) {
                curr = _pool_alloc_f(pfx + name, B, Hn, Wn, C);
                _avg_pool_downsample_float(prev, curr);
            }
        };
        auto ds_b = [&](TorchTensorView& prev, TorchTensorView& curr, const std::string& name) {
            if (_has(prev)) {
                curr = _pool_alloc_b(pfx + name, B, Hn, Wn);
                _avg_pool_downsample_bool(prev, curr);
            }
        };

        ds_f(render_rgb_s[sc-1], render_rgb_s[sc], "rrgb", 3);
        ds_f(ref_rgb_s[sc-1], ref_rgb_s[sc], "frgb", 3);
        ds_f(render_depth_s[sc-1], render_depth_s[sc], "rd", 1);
        ds_f(ref_depth_s[sc-1], ref_depth_s[sc], "fd", 1);
        ds_f(render_normal_s[sc-1], render_normal_s[sc], "rn", 3);
        ds_f(depth_normal_s[sc-1], depth_normal_s[sc], "dn", 3);
        ds_f(ref_normal_s[sc-1], ref_normal_s[sc], "fn", 3);
        ds_f(render_Ts_s[sc-1], render_Ts_s[sc], "rT", 1);
        ds_f(rgb_dist_s[sc-1], rgb_dist_s[sc], "rgbd", 3);
        ds_f(depth_dist_s[sc-1], depth_dist_s[sc], "dd", 1);
        ds_f(normal_dist_s[sc-1], normal_dist_s[sc], "nd", 3);
        ds_b(ref_alpha_s[sc-1], ref_alpha_s[sc], "ra");
        ds_b(mask_s[sc-1], mask_s[sc], "m");
        ds_b(depth_mask_s[sc-1], depth_mask_s[sc], "dm");
        ds_b(normal_mask_s[sc-1], normal_mask_s[sc], "nm");
        ds_b(alpha_mask_s[sc-1], alpha_mask_s[sc], "am");
    }

    // Total losses accumulator
    float* total_losses_ptr = _fptr(total_losses_out);
    cudaMemset(total_losses_ptr, 0, (uint)LossIndex::length * sizeof(float));

    float psnr_val = 0.0f, ssim_val = 0.0f;

    for (int scale = 0; scale < num_loss_scales; ++scale) {
        const auto& ss = std::get<2>(render_rgb_s[scale]);
        long Hs = ss[1], Ws = ss[2];
        size_t ppi = (size_t)(Hs * Ws);

        // Forward losses
        float* raw_losses_ptr = DevicePool::global().acquire<float>(
            "ppl.raw_losses", (size_t)(num_train_images+1) * (uint)RawLossIndex::length);
        cudaMemset(raw_losses_ptr, 0, (size_t)(num_train_images+1) * (uint)RawLossIndex::length * sizeof(float));

        float* losses_ptr = DevicePool::global().acquire<float>(
            "ppl.losses", (uint)LossIndex::length);
        cudaMemset(losses_ptr, 0, (uint)LossIndex::length * sizeof(float));

        float* loss_map_ptr = nullptr;
        TorchTensorView loss_map_scale = {};
        if (_has(loss_map_out)) {
            loss_map_scale = _pool_alloc_f_zero("ppl.loss_map_scale", B, Hs, Ws, 1);
            loss_map_ptr = _fptr(loss_map_scale);
        }

        _compute_per_pixel_losses_forward(
            B, ppi,
            render_rgb_s[scale], ref_rgb_s[scale], render_depth_s[scale], ref_depth_s[scale],
            render_normal_s[scale], depth_normal_s[scale], ref_normal_s[scale], render_Ts_s[scale],
            rgb_dist_s[scale], depth_dist_s[scale], normal_dist_s[scale],
            ref_alpha_s[scale], mask_s[scale], depth_mask_s[scale], normal_mask_s[scale], alpha_mask_s[scale],
            loss_weights, num_train_images, camera_indices,
            loss_map_ptr, raw_losses_ptr, losses_ptr
        );

        // Backward losses
        PerPixelGrads scale_grads = {};
        auto alloc_grad_f = [&](TorchTensorView& out, bool need, TorchTensorView& input, const std::string& name, int C) {
            if (need && _has(input))
                out = _pool_alloc_f("ppl.g." + name, B, Hs, Ws, C);
        };
        alloc_grad_f(scale_grads.v_render_rgb, needs_input_grad[0], render_rgb_s[scale], "vrgb", 3);
        alloc_grad_f(scale_grads.v_ref_rgb, needs_input_grad[1], ref_rgb_s[scale], "vfrgb", 3);
        alloc_grad_f(scale_grads.v_render_depth, needs_input_grad[2], render_depth_s[scale], "vrd", 1);
        alloc_grad_f(scale_grads.v_ref_depth, needs_input_grad[3], ref_depth_s[scale], "vfd", 1);
        alloc_grad_f(scale_grads.v_render_normal, needs_input_grad[4], render_normal_s[scale], "vrn", 3);
        alloc_grad_f(scale_grads.v_depth_normal, needs_input_grad[5], depth_normal_s[scale], "vdn", 3);
        alloc_grad_f(scale_grads.v_ref_normal, needs_input_grad[6], ref_normal_s[scale], "vfn", 3);
        alloc_grad_f(scale_grads.v_render_Ts, needs_input_grad[7], render_Ts_s[scale], "vrT", 1);
        alloc_grad_f(scale_grads.v_rgb_dist, needs_input_grad[8], rgb_dist_s[scale], "vrgbd", 3);
        alloc_grad_f(scale_grads.v_depth_dist, needs_input_grad[9], depth_dist_s[scale], "vdd", 1);
        alloc_grad_f(scale_grads.v_normal_dist, needs_input_grad[10], normal_dist_s[scale], "vnd", 3);

        _compute_per_pixel_losses_backward(
            B, ppi,
            render_rgb_s[scale], ref_rgb_s[scale], render_depth_s[scale], ref_depth_s[scale],
            render_normal_s[scale], depth_normal_s[scale], ref_normal_s[scale], render_Ts_s[scale],
            rgb_dist_s[scale], depth_dist_s[scale], normal_dist_s[scale],
            ref_alpha_s[scale], mask_s[scale], depth_mask_s[scale], normal_mask_s[scale], alpha_mask_s[scale],
            raw_losses_ptr, loss_weights, _fptr(v_losses),
            num_train_images, camera_indices, scale_grads
        );

        // SSIM backward (fused forward+backward)
        float ssim = fused_ssim_inplace(
            render_rgb_s[scale], ref_rgb_s[scale], mask_s[scale],
            -w_ssim,
            scale_grads.v_render_rgb,
            scale == 0,
            loss_map_scale,
            w_ssim
        );

        if (scale == 0)
            ssim_val = ssim;

        // Accumulate losses on device: total_losses += losses
        constexpr int N_LOSSES = (int)LossIndex::length;
        vector_add_scaled_kernel<<<1, N_LOSSES>>>(
            total_losses_ptr, losses_ptr, 1.0f, N_LOSSES);
        CHECK_DEVICE_ERROR(cudaGetLastError());

        // Upsample loss map
        if (_has(loss_map_out) && loss_map_ptr) {
            if (scale == 0) {
                cudaMemcpy(_fptr(loss_map_out), loss_map_ptr, B * H * W * sizeof(float), cudaMemcpyDeviceToDevice);
            } else {
                avg_pool_upsample_float_kernel<<<_LAUNCH_ARGS_3D(W, H, B, 16, 16, 1)>>>(
                    _make_view4(_fptr(loss_map_out), B, H, W, 1L),
                    _make_view4(loss_map_ptr, B, Hs, Ws, 1L),
                    1 << scale,
                    scale == 1 ? 1.0f / num_loss_scales : 1.0f,
                    1.0f / num_loss_scales
                );
                CHECK_DEVICE_ERROR(cudaGetLastError());
            }
        }

        // Upsample gradients and accumulate into grads_out
        auto upsample_grad = [&](TorchTensorView& grad_scale, TorchTensorView& grad_acc, int C) {
            if (_has(grad_acc) && scale == 0) {
                cudaMemcpy(_fptr(grad_acc), _fptr(grad_scale), B * H * W * C * sizeof(float), cudaMemcpyDeviceToDevice);
                return;
            }
            if (_has(grad_scale) && _has(grad_acc)) {
                float a = (scale == 1 ? 1.0f / num_loss_scales : 1.0f);
                float b = powf(0.25f, (float)scale) / num_loss_scales;
                avg_pool_upsample_float_kernel<<<_LAUNCH_ARGS_3D(W, H, B, 16, 16, 1)>>>(
                    _make_view4(_fptr(grad_acc), B, H, W, (long)C),
                    _make_view4(_fptr(grad_scale), B, Hs, Ws, (long)C),
                    1 << scale, a, b
                );
                CHECK_DEVICE_ERROR(cudaGetLastError());
            }
        };

        upsample_grad(scale_grads.v_render_rgb, grads_out.v_render_rgb, 3);
        upsample_grad(scale_grads.v_ref_rgb, grads_out.v_ref_rgb, 3);
        upsample_grad(scale_grads.v_render_depth, grads_out.v_render_depth, 1);
        upsample_grad(scale_grads.v_ref_depth, grads_out.v_ref_depth, 1);
        upsample_grad(scale_grads.v_render_normal, grads_out.v_render_normal, 3);
        upsample_grad(scale_grads.v_depth_normal, grads_out.v_depth_normal, 3);
        upsample_grad(scale_grads.v_ref_normal, grads_out.v_ref_normal, 3);
        upsample_grad(scale_grads.v_render_Ts, grads_out.v_render_Ts, 1);
        upsample_grad(scale_grads.v_rgb_dist, grads_out.v_rgb_dist, 3);
        upsample_grad(scale_grads.v_depth_dist, grads_out.v_depth_dist, 1);
        upsample_grad(scale_grads.v_normal_dist, grads_out.v_normal_dist, 3);
    }

    // Scale total losses by 1/num_scales (device-side)
    if (num_loss_scales > 1) {
        constexpr int N_LOSSES = (int)LossIndex::length;
        vector_add_scaled_kernel<<<1, N_LOSSES>>>(
            total_losses_ptr, total_losses_ptr,
            1.0f / (float)num_loss_scales - 1.0f, N_LOSSES);
        CHECK_DEVICE_ERROR(cudaGetLastError());
    }

    // psnr_val is not read here; caller reads it from total_losses_out after D2H
    return std::make_tuple(psnr_val, ssim_val);
}
