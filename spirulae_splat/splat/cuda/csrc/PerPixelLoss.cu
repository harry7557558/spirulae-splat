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

#include <ATen/ops/empty.h>
#include <ATen/ops/empty_like.h>
#include <ATen/ops/zeros.h>


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


/*[AutoHeaderGeneratorExport]*/
std::tuple<at::Tensor, at::Tensor, std::optional<at::Tensor>>
compute_per_pixel_losses_forward_tensor(
    std::optional<at::Tensor> render_rgb,
    std::optional<at::Tensor> ref_rgb,
    std::optional<at::Tensor> render_depth,
    std::optional<at::Tensor> ref_depth,
    std::optional<at::Tensor> render_normal,
    std::optional<at::Tensor> depth_normal,
    std::optional<at::Tensor> ref_normal,
    std::optional<at::Tensor> render_Ts,
    std::optional<at::Tensor> rgb_dist,
    std::optional<at::Tensor> depth_dist,
    std::optional<at::Tensor> normal_dist,
    std::optional<at::Tensor> ref_alpha,
    std::optional<at::Tensor> mask,
    std::optional<at::Tensor> depth_mask,
    std::optional<at::Tensor> normal_mask,
    std::optional<at::Tensor> alpha_mask,
    const std::array<float, (int)LossWeightIndex::length> loss_weights_0,
    long num_train_images,
    std::optional<at::Tensor> camera_indices,
    bool return_loss_map
) {
    DEVICE_GUARD(render_rgb.value());

    long B = -1, H = -1, W = -1;
    auto check_generic = [&](std::string name, const at::Tensor& tensor) {
        CHECK_CUDA(tensor);
        if (tensor.ndimension() != 4)
            AT_ERROR(name + " must be (B, H, W, C)");
        if (B == -1)
            B = tensor.size(0), H = tensor.size(1), W = tensor.size(2);
        else if (B != tensor.size(0) || H != tensor.size(1) || W != tensor.size(2))
            AT_ERROR("Tensor shape mismatch with render_rgb (" + name + ")");
    };
    auto check_float = [&](std::string name, const std::optional<at::Tensor>& tensor, int ncomp) {
        if (tensor.has_value()) {
            check_generic(name, tensor.value());
            if (tensor.value().size(-1) != ncomp)
                AT_ERROR("Last dimension of " + name + " must be " + std::to_string(ncomp));
        }
    };

    if (!render_rgb.has_value())
        AT_ERROR("render_rgb must be provided");
    check_float("render_rgb", render_rgb, 3);
    check_float("ref_rgb", ref_rgb, 3);
    check_float("render_depth", render_depth, 1);
    check_float("ref_depth", ref_depth, 1);
    check_float("render_normal", render_normal, 3);
    check_float("depth_normal", depth_normal, 3);
    check_float("ref_normal", ref_normal, 3);
    check_float("render_Ts", render_Ts, 1);
    check_float("rgb_dist", rgb_dist, 3);
    check_float("depth_dist", depth_dist, 1);
    check_float("normal_dist", normal_dist, 3);
    check_float("ref_alpha", ref_alpha, 1);
    check_float("mask", mask, 1);
    check_float("depth_mask", depth_mask, 1);
    check_float("normal_mask", normal_mask, 1);
    check_float("alpha_mask", alpha_mask, 1);

    size_t pixels_per_image = render_rgb.value().numel() / (3 * B);

    FixedArray<float, (uint)LossWeightIndex::length> loss_weights =
        *reinterpret_cast<const FixedArray<float, (uint)LossWeightIndex::length>*>(loss_weights_0.data());

    if (!camera_indices.has_value())
        num_train_images = B;
    at::Tensor raw_losses = at::zeros({num_train_images+1, (uint)RawLossIndex::length}, render_rgb.value().options());
    at::Tensor losses = at::zeros({(uint)LossIndex::length}, render_rgb.value().options());
    std::optional<at::Tensor> loss_map;
    if (return_loss_map) {
        loss_map = at::empty({B, H, W, 1}, render_rgb.value().options());
        set_zero_tensor(loss_map.value());
    }

    per_pixel_losses_forward_kernel<<<_LAUNCH_ARGS_2D(pixels_per_image, B, WARP_SIZE*WARP_SIZE, 1)>>>(
        B, pixels_per_image,
        camera_indices.has_value() ? camera_indices.value().contiguous().data_ptr<int64_t>() : nullptr,
        render_rgb.has_value() ? (float3*)render_rgb.value().contiguous().data_ptr<float>() : nullptr,
        ref_rgb.has_value() ? (float3*)ref_rgb.value().contiguous().data_ptr<float>() : nullptr,
        render_depth.has_value() ? (float*)render_depth.value().contiguous().data_ptr<float>() : nullptr,
        ref_depth.has_value() ? (float*)ref_depth.value().contiguous().data_ptr<float>() : nullptr,
        render_normal.has_value() ? (float3*)render_normal.value().contiguous().data_ptr<float>() : nullptr,
        depth_normal.has_value() ? (float3*)depth_normal.value().contiguous().data_ptr<float>() : nullptr,
        ref_normal.has_value() ? (float3*)ref_normal.value().contiguous().data_ptr<float>() : nullptr,
        render_Ts.has_value() ? (float*)render_Ts.value().contiguous().data_ptr<float>() : nullptr,
        rgb_dist.has_value() ? (float3*)rgb_dist.value().contiguous().data_ptr<float>() : nullptr,
        depth_dist.has_value() ? (float*)depth_dist.value().contiguous().data_ptr<float>() : nullptr,
        normal_dist.has_value() ? (float3*)normal_dist.value().contiguous().data_ptr<float>() : nullptr,
        ref_alpha.has_value() ? (bool*)ref_alpha.value().contiguous().data_ptr<bool>() : nullptr,
        mask.has_value() ? (bool*)mask.value().contiguous().data_ptr<bool>() : nullptr,
        depth_mask.has_value() ? (bool*)depth_mask.value().contiguous().data_ptr<bool>() : nullptr,
        normal_mask.has_value() ? (bool*)normal_mask.value().contiguous().data_ptr<bool>() : nullptr,
        alpha_mask.has_value() ? (bool*)alpha_mask.value().contiguous().data_ptr<bool>() : nullptr,
        loss_weights,
        return_loss_map ? loss_map.value().data_ptr<float>() : nullptr,
        raw_losses.data_ptr<float>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    per_pixel_losses_reduce_forward_kernel
    <<<_LAUNCH_ARGS_1D(num_train_images+1, WARP_SIZE)>>>(
        num_train_images,
        raw_losses.data_ptr<float>(),
        loss_weights,
        losses.data_ptr<float>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return std::make_tuple(losses, raw_losses, loss_map);
}


/*[AutoHeaderGeneratorExport]*/
std::tuple<  // returns gradients
    std::optional<at::Tensor>, // render_rgb
    std::optional<at::Tensor>, // ref_rgb
    std::optional<at::Tensor>, // render_depth
    std::optional<at::Tensor>, // ref_depth
    std::optional<at::Tensor>, // render_normal
    std::optional<at::Tensor>, // depth_normal
    std::optional<at::Tensor>, // ref_normal
    std::optional<at::Tensor>, // render_Ts
    std::optional<at::Tensor>, // rgb_dist
    std::optional<at::Tensor>, // depth_dist
    std::optional<at::Tensor> // normal_dist
> compute_per_pixel_losses_backward_tensor(
    std::optional<at::Tensor> render_rgb,
    std::optional<at::Tensor> ref_rgb,
    std::optional<at::Tensor> render_depth,
    std::optional<at::Tensor> ref_depth,
    std::optional<at::Tensor> render_normal,
    std::optional<at::Tensor> depth_normal,
    std::optional<at::Tensor> ref_normal,
    std::optional<at::Tensor> render_Ts,
    std::optional<at::Tensor> rgb_dist,
    std::optional<at::Tensor> depth_dist,
    std::optional<at::Tensor> normal_dist,
    std::optional<at::Tensor> ref_alpha,
    std::optional<at::Tensor> mask,
    std::optional<at::Tensor> depth_mask,
    std::optional<at::Tensor> normal_mask,
    std::optional<at::Tensor> alpha_mask,
    at::Tensor raw_losses,
    const std::array<float, (int)LossWeightIndex::length> loss_weights_0,
    at::Tensor v_losses,
    std::vector<bool> needs_input_grad,
    long num_train_images,
    std::optional<at::Tensor> camera_indices
) {
    DEVICE_GUARD(render_rgb.value());

    long B = render_rgb.value().size(0);
    size_t pixels_per_image = render_rgb.value().numel() / (3 * B);

    std::optional<at::Tensor> v_render_rgb = needs_input_grad[0] && render_rgb.has_value() ? (std::optional<at::Tensor>)at::empty_like(render_rgb.value()) : std::nullopt;
    std::optional<at::Tensor> v_ref_rgb = needs_input_grad[1] && ref_rgb.has_value() ? (std::optional<at::Tensor>)at::empty_like(ref_rgb.value()) : std::nullopt;
    std::optional<at::Tensor> v_render_depth = needs_input_grad[2] && render_depth.has_value() ? (std::optional<at::Tensor>)at::empty_like(render_depth.value()) : std::nullopt;
    std::optional<at::Tensor> v_ref_depth = needs_input_grad[3] && ref_depth.has_value() ? (std::optional<at::Tensor>)at::empty_like(ref_depth.value()) : std::nullopt;
    std::optional<at::Tensor> v_render_normal = needs_input_grad[4] && render_normal.has_value() ? (std::optional<at::Tensor>)at::empty_like(render_normal.value()) : std::nullopt;
    std::optional<at::Tensor> v_depth_normal = needs_input_grad[5] && depth_normal.has_value() ? (std::optional<at::Tensor>)at::empty_like(depth_normal.value()) : std::nullopt;
    std::optional<at::Tensor> v_ref_normal = needs_input_grad[6] && ref_normal.has_value() ? (std::optional<at::Tensor>)at::empty_like(ref_normal.value()) : std::nullopt;
    std::optional<at::Tensor> v_render_Ts = needs_input_grad[7] && render_Ts.has_value() ? (std::optional<at::Tensor>)at::empty_like(render_Ts.value()) : std::nullopt;
    std::optional<at::Tensor> v_rgb_dist = needs_input_grad[8] && rgb_dist.has_value() ? (std::optional<at::Tensor>)at::empty_like(rgb_dist.value()) : std::nullopt;
    std::optional<at::Tensor> v_depth_dist = needs_input_grad[9] && depth_dist.has_value() ? (std::optional<at::Tensor>)at::empty_like(depth_dist.value()) : std::nullopt;
    std::optional<at::Tensor> v_normal_dist = needs_input_grad[10] && normal_dist.has_value() ? (std::optional<at::Tensor>)at::empty_like(normal_dist.value()) : std::nullopt;

    FixedArray<float, (uint)LossWeightIndex::length> loss_weights =
        *reinterpret_cast<const FixedArray<float, (uint)LossWeightIndex::length>*>(loss_weights_0.data());

    if (!camera_indices.has_value())
        num_train_images = B;
    at::Tensor v_raw_losses = at::empty({num_train_images+1, (uint)RawLossIndex::length}, render_rgb.value().options());

    per_pixel_losses_reduce_backward_kernel<<<_LAUNCH_ARGS_1D(num_train_images+1, WARP_SIZE)>>>(
        num_train_images,
        raw_losses.data_ptr<float>(),
        loss_weights,
        v_losses.data_ptr<float>(),
        v_raw_losses.data_ptr<float>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    per_pixel_losses_backward_kernel<<<_LAUNCH_ARGS_2D(pixels_per_image, B, WARP_SIZE*WARP_SIZE, 1)>>>(
        B, pixels_per_image,
        camera_indices.has_value() ? camera_indices.value().contiguous().data_ptr<int64_t>() : nullptr,
        render_rgb.has_value() ? (float3*)render_rgb.value().contiguous().data_ptr<float>() : nullptr,
        ref_rgb.has_value() ? (float3*)ref_rgb.value().contiguous().data_ptr<float>() : nullptr,
        render_depth.has_value() ? (float*)render_depth.value().contiguous().data_ptr<float>() : nullptr,
        ref_depth.has_value() ? (float*)ref_depth.value().contiguous().data_ptr<float>() : nullptr,
        render_normal.has_value() ? (float3*)render_normal.value().contiguous().data_ptr<float>() : nullptr,
        depth_normal.has_value() ? (float3*)depth_normal.value().contiguous().data_ptr<float>() : nullptr,
        ref_normal.has_value() ? (float3*)ref_normal.value().contiguous().data_ptr<float>() : nullptr,
        render_Ts.has_value() ? (float*)render_Ts.value().contiguous().data_ptr<float>() : nullptr,
        rgb_dist.has_value() ? (float3*)rgb_dist.value().contiguous().data_ptr<float>() : nullptr,
        depth_dist.has_value() ? (float*)depth_dist.value().contiguous().data_ptr<float>() : nullptr,
        normal_dist.has_value() ? (float3*)normal_dist.value().contiguous().data_ptr<float>() : nullptr,
        ref_alpha.has_value() ? (bool*)ref_alpha.value().contiguous().data_ptr<bool>() : nullptr,
        mask.has_value() ? (bool*)mask.value().contiguous().data_ptr<bool>() : nullptr,
        depth_mask.has_value() ? (bool*)depth_mask.value().contiguous().data_ptr<bool>() : nullptr,
        normal_mask.has_value() ? (bool*)normal_mask.value().contiguous().data_ptr<bool>() : nullptr,
        alpha_mask.has_value() ? (bool*)alpha_mask.value().contiguous().data_ptr<bool>() : nullptr,
        loss_weights,
        v_raw_losses.data_ptr<float>(),
        v_render_rgb.has_value() ? (float3*)v_render_rgb.value().contiguous().data_ptr<float>() : nullptr,
        v_ref_rgb.has_value() ? (float3*)v_ref_rgb.value().contiguous().data_ptr<float>() : nullptr,
        v_render_depth.has_value() ? (float*)v_render_depth.value().contiguous().data_ptr<float>() : nullptr,
        v_ref_depth.has_value() ? (float*)v_ref_depth.value().contiguous().data_ptr<float>() : nullptr,
        v_render_normal.has_value() ? (float3*)v_render_normal.value().contiguous().data_ptr<float>() : nullptr,
        v_depth_normal.has_value() ? (float3*)v_depth_normal.value().contiguous().data_ptr<float>() : nullptr,
        v_ref_normal.has_value() ? (float3*)v_ref_normal.value().contiguous().data_ptr<float>() : nullptr,
        v_render_Ts.has_value() ? (float*)v_render_Ts.value().contiguous().data_ptr<float>() : nullptr,
        v_rgb_dist.has_value() ? (float3*)v_rgb_dist.value().contiguous().data_ptr<float>() : nullptr,
        v_depth_dist.has_value() ? (float*)v_depth_dist.value().contiguous().data_ptr<float>() : nullptr,
        v_normal_dist.has_value() ? (float3*)v_normal_dist.value().contiguous().data_ptr<float>() : nullptr
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    // TODO: investigate why program crashes without this
    // if (camera_indices.has_value())
    //     cudaDeviceSynchronize();

    return std::make_tuple(
        v_render_rgb,
        v_ref_rgb,
        v_render_depth,
        v_ref_depth,
        v_render_normal,
        v_depth_normal,
        v_ref_normal,
        v_render_Ts,
        v_rgb_dist,
        v_depth_dist,
        v_normal_dist
    );
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


/*[AutoHeaderGeneratorExport]*/
at::Tensor avg_pool_downsample_tensor(
    at::Tensor tensor
) {
    DEVICE_GUARD(tensor);
    CHECK_INPUT(tensor);

    auto h = tensor.size(1) / 2, w = tensor.size(2) / 2, b = tensor.size(0);
    at::Tensor result = at::empty({b, h, w, tensor.size(3)}, tensor.options());

    if (tensor.dtype() == at::kFloat)
        avg_pool_downsample_float_kernel<<<_LAUNCH_ARGS_3D(w, h, b, 16, 16, 1)>>>
            (tensor2view<float, 4>(tensor), tensor2view<float, 4>(result));
    else if (tensor.dtype() == at::kUInt16)
        avg_pool_downsample_integral_kernel<uint16_t><<<_LAUNCH_ARGS_3D(w, h, b, 16, 16, 1)>>>
            (tensor2view<uint16_t, 4>(tensor), tensor2view<uint16_t, 4>(result));
    else if (tensor.dtype() == at::kByte)
        avg_pool_downsample_integral_kernel<uint8_t><<<_LAUNCH_ARGS_3D(w, h, b, 16, 16, 1)>>>
            (tensor2view<uint8_t, 4>(tensor), tensor2view<uint8_t, 4>(result));
    else if (tensor.dtype() == at::kBool) {
        auto tensor_byte = tensor.view(at::kByte);
        auto result_byte = result.view(at::kByte);
        avg_pool_downsample_bool_kernel<<<_LAUNCH_ARGS_3D(w, h, b, 16, 16, 1)>>>
            (tensor2view<uint8_t, 4>(tensor_byte), tensor2view<uint8_t, 4>(result_byte));
    }
    else throw std::runtime_error("Unsupported dtype for avg_pool_downsample");
    
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return result;
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
std::tuple<
    at::Tensor,  // mean losses
    std::optional<at::Tensor>,  // loss map
    std::optional<std::tuple<
        std::optional<at::Tensor>, // render_rgb
        std::optional<at::Tensor>, // ref_rgb
        std::optional<at::Tensor>, // render_depth
        std::optional<at::Tensor>, // ref_depth
        std::optional<at::Tensor>, // render_normal
        std::optional<at::Tensor>, // depth_normal
        std::optional<at::Tensor>, // ref_normal
        std::optional<at::Tensor>, // render_Ts
        std::optional<at::Tensor>, // rgb_dist
        std::optional<at::Tensor>, // depth_dist
        std::optional<at::Tensor> // normal_dist
    >>,
    std::tuple<
        float,  // psnr value
        float  // ssim value
    >
> compute_multi_scale_per_pixel_losses_tensor(
    int num_loss_scales,
    std::optional<at::Tensor> render_rgb,
    std::optional<at::Tensor> ref_rgb,
    std::optional<at::Tensor> render_depth,
    std::optional<at::Tensor> ref_depth,
    std::optional<at::Tensor> render_normal,
    std::optional<at::Tensor> depth_normal,
    std::optional<at::Tensor> ref_normal,
    std::optional<at::Tensor> render_Ts,
    std::optional<at::Tensor> rgb_dist,
    std::optional<at::Tensor> depth_dist,
    std::optional<at::Tensor> normal_dist,
    std::optional<at::Tensor> ref_alpha,
    std::optional<at::Tensor> mask,
    std::optional<at::Tensor> depth_mask,
    std::optional<at::Tensor> normal_mask,
    std::optional<at::Tensor> alpha_mask,
    const std::array<float, (int)LossWeightIndex::length> loss_weights_0,
    const float w_ssim,
    std::optional<at::Tensor> v_losses,
    std::vector<bool> needs_input_grad,
    long num_train_images,
    std::optional<at::Tensor> camera_indices,
    bool return_loss_map
) {
    DEVICE_GUARD(render_rgb.value());

    long B = render_rgb.value().size(0);
    long H = render_rgb.value().size(1);
    long W = render_rgb.value().size(2);

    // Create scale tensors
    std::vector<std::optional<at::Tensor>> render_rgb_scales(num_loss_scales);
    std::vector<std::optional<at::Tensor>> ref_rgb_scales(num_loss_scales);
    std::vector<std::optional<at::Tensor>> render_depth_scales(num_loss_scales);
    std::vector<std::optional<at::Tensor>> ref_depth_scales(num_loss_scales);
    std::vector<std::optional<at::Tensor>> render_normal_scales(num_loss_scales);
    std::vector<std::optional<at::Tensor>> depth_normal_scales(num_loss_scales);
    std::vector<std::optional<at::Tensor>> ref_normal_scales(num_loss_scales);
    std::vector<std::optional<at::Tensor>> render_Ts_scales(num_loss_scales);
    std::vector<std::optional<at::Tensor>> rgb_dist_scales(num_loss_scales);
    std::vector<std::optional<at::Tensor>> depth_dist_scales(num_loss_scales);
    std::vector<std::optional<at::Tensor>> normal_dist_scales(num_loss_scales);
    std::vector<std::optional<at::Tensor>> ref_alpha_scales(num_loss_scales);
    std::vector<std::optional<at::Tensor>> mask_scales(num_loss_scales);
    std::vector<std::optional<at::Tensor>> depth_mask_scales(num_loss_scales);
    std::vector<std::optional<at::Tensor>> normal_mask_scales(num_loss_scales);
    std::vector<std::optional<at::Tensor>> alpha_mask_scales(num_loss_scales);

    render_rgb_scales[0] = render_rgb;
    ref_rgb_scales[0] = ref_rgb;
    render_depth_scales[0] = render_depth;
    ref_depth_scales[0] = ref_depth;
    render_normal_scales[0] = render_normal;
    depth_normal_scales[0] = depth_normal;
    ref_normal_scales[0] = ref_normal;
    render_Ts_scales[0] = render_Ts;
    rgb_dist_scales[0] = rgb_dist;
    depth_dist_scales[0] = depth_dist;
    normal_dist_scales[0] = normal_dist;
    ref_alpha_scales[0] = ref_alpha;
    mask_scales[0] = mask;
    depth_mask_scales[0] = depth_mask;
    normal_mask_scales[0] = normal_mask;
    alpha_mask_scales[0] = alpha_mask;

    // Downsample to create scales
    for (int scale = 1; scale < num_loss_scales; ++scale) {
        long H_prev = render_rgb_scales[scale-1].value().size(1);
        long W_prev = render_rgb_scales[scale-1].value().size(2);
        long H_new = H_prev / 2;
        long W_new = W_prev / 2;

        // Float tensors
        auto downsample_float = [&](std::optional<at::Tensor>& prev, std::optional<at::Tensor>& curr, int C) {
            if (prev.has_value()) {
                curr = at::empty({B, H_new, W_new, C}, prev.value().options());
                avg_pool_downsample_float_kernel<<<_LAUNCH_ARGS_3D(W_new, H_new, B, 16, 16, 1)>>>(
                    tensor2view<float, 4>(prev.value()),
                    tensor2view<float, 4>(curr.value())
                );
                CHECK_DEVICE_ERROR(cudaGetLastError());
            }
        };

        downsample_float(render_rgb_scales[scale-1], render_rgb_scales[scale], 3);
        downsample_float(ref_rgb_scales[scale-1], ref_rgb_scales[scale], 3);
        downsample_float(render_depth_scales[scale-1], render_depth_scales[scale], 1);
        downsample_float(ref_depth_scales[scale-1], ref_depth_scales[scale], 1);
        downsample_float(render_normal_scales[scale-1], render_normal_scales[scale], 3);
        downsample_float(depth_normal_scales[scale-1], depth_normal_scales[scale], 3);
        downsample_float(ref_normal_scales[scale-1], ref_normal_scales[scale], 3);
        downsample_float(render_Ts_scales[scale-1], render_Ts_scales[scale], 1);
        downsample_float(rgb_dist_scales[scale-1], rgb_dist_scales[scale], 3);
        downsample_float(depth_dist_scales[scale-1], depth_dist_scales[scale], 1);
        downsample_float(normal_dist_scales[scale-1], normal_dist_scales[scale], 3);

        // Bool tensors
        auto downsample_bool = [&](std::optional<at::Tensor>& prev, std::optional<at::Tensor>& curr) {
            if (prev.has_value()) {
                curr = at::empty({B, H_new, W_new, 1}, prev.value().options());
                auto prev_byte = prev.value().view(at::kByte);
                auto curr_byte = curr.value().view(at::kByte);
                avg_pool_downsample_bool_kernel<<<_LAUNCH_ARGS_3D(W_new, H_new, B, 16, 16, 1)>>>(
                    tensor2view<uint8_t, 4>(prev_byte),
                    tensor2view<uint8_t, 4>(curr_byte)
                );
                CHECK_DEVICE_ERROR(cudaGetLastError());
            }
        };

        downsample_bool(ref_alpha_scales[scale-1], ref_alpha_scales[scale]);
        downsample_bool(mask_scales[scale-1], mask_scales[scale]);
        downsample_bool(depth_mask_scales[scale-1], depth_mask_scales[scale]);
        downsample_bool(normal_mask_scales[scale-1], normal_mask_scales[scale]);
        downsample_bool(alpha_mask_scales[scale-1], alpha_mask_scales[scale]);
    }

    at::Tensor total_losses = at::zeros({(uint)LossIndex::length}, render_rgb.value().options());
    std::optional<at::Tensor> total_loss_map;
    float psnr_val = 0.0f;
    float ssim_val = 0.0f;

    std::optional<at::Tensor> v_render_rgb_acc = needs_input_grad[0] && render_rgb.has_value() ? (std::optional<at::Tensor>)at::empty({0}, render_rgb.value().options()) : std::nullopt;
    std::optional<at::Tensor> v_ref_rgb_acc = needs_input_grad[1] && ref_rgb.has_value() ? (std::optional<at::Tensor>)at::empty({0}, ref_rgb.value().options()) : std::nullopt;
    std::optional<at::Tensor> v_render_depth_acc = needs_input_grad[2] && render_depth.has_value() ? (std::optional<at::Tensor>)at::empty({0}, render_depth.value().options()) : std::nullopt;
    std::optional<at::Tensor> v_ref_depth_acc = needs_input_grad[3] && ref_depth.has_value() ? (std::optional<at::Tensor>)at::empty({0}, ref_depth.value().options()) : std::nullopt;
    std::optional<at::Tensor> v_render_normal_acc = needs_input_grad[4] && render_normal.has_value() ? (std::optional<at::Tensor>)at::empty({0}, render_normal.value().options()) : std::nullopt;
    std::optional<at::Tensor> v_depth_normal_acc = needs_input_grad[5] && depth_normal.has_value() ? (std::optional<at::Tensor>)at::empty({0}, depth_normal.value().options()) : std::nullopt;
    std::optional<at::Tensor> v_ref_normal_acc = needs_input_grad[6] && ref_normal.has_value() ? (std::optional<at::Tensor>)at::empty({0}, ref_normal.value().options()) : std::nullopt;
    std::optional<at::Tensor> v_render_Ts_acc = needs_input_grad[7] && render_Ts.has_value() ? (std::optional<at::Tensor>)at::empty({0}, render_Ts.value().options()) : std::nullopt;
    std::optional<at::Tensor> v_rgb_dist_acc = needs_input_grad[8] && rgb_dist.has_value() ? (std::optional<at::Tensor>)at::empty({0}, rgb_dist.value().options()) : std::nullopt;
    std::optional<at::Tensor> v_depth_dist_acc = needs_input_grad[9] && depth_dist.has_value() ? (std::optional<at::Tensor>)at::empty({0}, depth_dist.value().options()) : std::nullopt;
    std::optional<at::Tensor> v_normal_dist_acc = needs_input_grad[10] && normal_dist.has_value() ? (std::optional<at::Tensor>)at::empty({0}, normal_dist.value().options()) : std::nullopt;

    std::optional<std::tuple<
        std::optional<at::Tensor>, std::optional<at::Tensor>, std::optional<at::Tensor>, std::optional<at::Tensor>,
        std::optional<at::Tensor>, std::optional<at::Tensor>, std::optional<at::Tensor>, std::optional<at::Tensor>,
        std::optional<at::Tensor>, std::optional<at::Tensor>, std::optional<at::Tensor>
    >> grads;

    for (int scale = 0; scale < num_loss_scales; ++scale) {

        auto [losses, raw_losses, loss_map] = compute_per_pixel_losses_forward_tensor(
            render_rgb_scales[scale], ref_rgb_scales[scale], render_depth_scales[scale], ref_depth_scales[scale],
            render_normal_scales[scale], depth_normal_scales[scale], ref_normal_scales[scale], render_Ts_scales[scale],
            rgb_dist_scales[scale], depth_dist_scales[scale], normal_dist_scales[scale],
            ref_alpha_scales[scale], mask_scales[scale], depth_mask_scales[scale], normal_mask_scales[scale], alpha_mask_scales[scale],
            loss_weights_0, num_train_images, camera_indices, return_loss_map
        );
        auto grad_tuple = compute_per_pixel_losses_backward_tensor(
            render_rgb_scales[scale], ref_rgb_scales[scale], render_depth_scales[scale], ref_depth_scales[scale],
            render_normal_scales[scale], depth_normal_scales[scale], ref_normal_scales[scale], render_Ts_scales[scale],
            rgb_dist_scales[scale], depth_dist_scales[scale], normal_dist_scales[scale],
            ref_alpha_scales[scale], mask_scales[scale], depth_mask_scales[scale], normal_mask_scales[scale], alpha_mask_scales[scale],
            raw_losses, loss_weights_0, v_losses.value(), needs_input_grad, num_train_images, camera_indices
        );

        // Backward
        float ssim = fused_ssim_inplace(
            render_rgb_scales[scale].value(),
            ref_rgb_scales[scale].value(),
            mask_scales[scale],
            -w_ssim,
            std::get<0>(grad_tuple).value(),
            scale == 0,
            return_loss_map && loss_map.has_value() ?
                loss_map.value() :
                (std::optional<at::Tensor>)std::nullopt,
            w_ssim
        );
        if (scale == 0) {
            psnr_val = losses[1].item<float>();
            ssim_val = ssim;  // TODO: this isn't accurate in masked mode
        }

        total_losses += losses;
        if (return_loss_map && loss_map.has_value()) {
            if (scale == 0) {
                total_loss_map = loss_map;
            } else {
                long H_scale = loss_map.value().size(1);
                long W_scale = loss_map.value().size(2);
                avg_pool_upsample_float_kernel<<<_LAUNCH_ARGS_3D(W, H, B, 16, 16, 1)>>>(
                    tensor2view<float, 4>(total_loss_map.value()),
                    tensor2view<float, 4>(loss_map.value()),
                    1 << scale,
                    scale == 1 ? 1.0f / num_loss_scales : 1.0f,
                    1.0f / num_loss_scales
                );
                CHECK_DEVICE_ERROR(cudaGetLastError());
            }
        }

        auto upsample_grad = [&](std::optional<at::Tensor>& grad_scale, std::optional<at::Tensor>& grad_acc, int C) {
            if (grad_acc.has_value() && scale == 0) {
                grad_acc = grad_scale;
                return;
            }
            if (grad_scale.has_value() && grad_acc.has_value()) {
                long H_scale = grad_scale.value().size(1);
                long W_scale = grad_scale.value().size(2);
                float a = (scale == 1 ? 1.0f / num_loss_scales : 1.0f);
                float b = powf(0.25f, (float)scale) / num_loss_scales;
                avg_pool_upsample_float_kernel<<<_LAUNCH_ARGS_3D(W, H, B, 16, 16, 1)>>>(
                    tensor2view<float, 4>(grad_acc.value()),
                    tensor2view<float, 4>(grad_scale.value()),
                    1 << scale, a, b
                );
                CHECK_DEVICE_ERROR(cudaGetLastError());
            }
        };

        upsample_grad(std::get<0>(grad_tuple), v_render_rgb_acc, 3);
        upsample_grad(std::get<1>(grad_tuple), v_ref_rgb_acc, 3);
        upsample_grad(std::get<2>(grad_tuple), v_render_depth_acc, 1);
        upsample_grad(std::get<3>(grad_tuple), v_ref_depth_acc, 1);
        upsample_grad(std::get<4>(grad_tuple), v_render_normal_acc, 3);
        upsample_grad(std::get<5>(grad_tuple), v_depth_normal_acc, 3);
        upsample_grad(std::get<6>(grad_tuple), v_ref_normal_acc, 3);
        upsample_grad(std::get<7>(grad_tuple), v_render_Ts_acc, 1);
        upsample_grad(std::get<8>(grad_tuple), v_rgb_dist_acc, 3);
        upsample_grad(std::get<9>(grad_tuple), v_depth_dist_acc, 1);
        upsample_grad(std::get<10>(grad_tuple), v_normal_dist_acc, 3);

    }
    total_losses *= (1.0f / (float)num_loss_scales);

    if (v_losses.has_value())  // TODO: don't do backward when this is False
        grads = std::make_tuple(
            v_render_rgb_acc, v_ref_rgb_acc, v_render_depth_acc, v_ref_depth_acc,
            v_render_normal_acc, v_depth_normal_acc, v_ref_normal_acc, v_render_Ts_acc,
            v_rgb_dist_acc, v_depth_dist_acc, v_normal_dist_acc
        );

    return std::make_tuple(
        total_losses, total_loss_map, grads,
        std::make_tuple(psnr_val, ssim_val)
    );
}
