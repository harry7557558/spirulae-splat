#include "PerPixelLoss.cuh"

#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>
namespace cg = cooperative_groups;

#define TensorView _Slang_TensorView
#include "generated/slang_all.cu"
#undef TensorView

#include "common.cuh"


__global__ void per_pixel_losses_forward_kernel(
    const size_t num_pixels,
    const float3* __restrict__ render_rgb,
    const float3* __restrict__ ref_rgb,
    const float* __restrict__ render_depth,
    const float* __restrict__ ref_depth,
    const float3* __restrict__ render_normal,
    const float3* __restrict__ depth_normal,
    const float3* __restrict__ ref_normal,
    const float* __restrict__ render_alpha,
    const float3* __restrict__ rgb_dist,
    const float* __restrict__ depth_dist,
    const float3* __restrict__ normal_dist,
    const bool* __restrict__ ref_alpha,
    const bool* __restrict__ mask,
    const bool* __restrict__ depth_mask,
    const bool* __restrict__ normal_mask,
    float* __restrict__ out_losses
) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;

    FixedArray<float, (uint)RawLossIndex::length> losses;

    bool inside = idx < num_pixels;
    if (inside) {
        per_pixel_losses(
            render_rgb ? render_rgb[idx] : make_float3(0),
            ref_rgb ? ref_rgb[idx] : make_float3(0),
            render_depth ? render_depth[idx] : 1.f,
            ref_depth ? ref_depth[idx] : 1.f,
            render_normal ? render_normal[idx] : make_float3(0),
            depth_normal ? depth_normal[idx] : make_float3(0),
            ref_normal ? ref_normal[idx] : make_float3(0),
            render_alpha ? render_alpha[idx] : 0.f,
            rgb_dist ? rgb_dist[idx] : make_float3(0),
            depth_dist ? depth_dist[idx] : 0.f,
            normal_dist ? normal_dist[idx] : make_float3(0),
            ref_alpha ? ref_alpha[idx] : true,
            mask ? mask[idx] : true,
            depth_mask ? depth_mask[idx] : true,
            normal_mask ? normal_mask[idx] : true,
            &losses
        );
    }

    auto block = cg::this_thread_block();
    cg::thread_block_tile<WARP_SIZE> warp = cg::tiled_partition<WARP_SIZE>(block);
    uint warp_idx = block.thread_rank() / WARP_SIZE;

    __shared__ float atomic_reduce[WARP_SIZE];

    for (uint i = 0; i < (uint)RawLossIndex::length; i++) {
        float loss = inside ? losses[i] : 0.0f;
        float loss_reduced = cg::reduce(warp, loss, cg::plus<float>());
        if (warp.thread_rank() == 0)
            atomic_reduce[warp_idx] = loss_reduced;
        __syncthreads();
        loss = (warp_idx == 0) ? atomic_reduce[warp.thread_rank()] : 0.0f;
        if (__ballot_sync(~0u, loss != 0.0f) == 0)
            continue;
        loss_reduced = cg::reduce(warp, loss, cg::plus<float>());
        if (block.thread_rank() == 0 && loss_reduced != 0.0f)
            atomicAdd(out_losses+i, loss_reduced);
    }
}

__global__ void per_pixel_losses_backward_kernel(
    const size_t num_pixels,
    const float3* __restrict__ render_rgb,
    const float3* __restrict__ ref_rgb,
    const float* __restrict__ render_depth,
    const float* __restrict__ ref_depth,
    const float3* __restrict__ render_normal,
    const float3* __restrict__ depth_normal,
    const float3* __restrict__ ref_normal,
    const float* __restrict__ render_alpha,
    const float3* __restrict__ rgb_dist,
    const float* __restrict__ depth_dist,
    const float3* __restrict__ normal_dist,
    const bool* __restrict__ ref_alpha,
    const bool* __restrict__ mask,
    const bool* __restrict__ depth_mask,
    const bool* __restrict__ normal_mask,
    const float* __restrict__ v_out_losses,
    float3* __restrict__ v_render_rgb,
    float3* __restrict__ v_ref_rgb,
    float* __restrict__ v_render_depth,
    float* __restrict__ v_ref_depth,
    float3* __restrict__ v_render_normal,
    float3* __restrict__ v_depth_normal,
    float3* __restrict__ v_ref_normal,
    float* __restrict__ v_render_alpha,
    float3* __restrict__ v_rgb_dist,
    float* __restrict__ v_depth_dist,
    float3* __restrict__ v_normal_dist
) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;

    bool inside = idx < num_pixels;
    if (!inside) return;

    FixedArray<float, (uint)RawLossIndex::length> v_losses;
    for (uint i = 0; i < (uint)RawLossIndex::length; i++)
        v_losses[i] = v_out_losses[i];

    float3 temp_v_render_rgb;
    float3 temp_v_ref_rgb;
    float temp_v_render_depth;
    float temp_v_ref_depth;
    float3 temp_v_render_normal;
    float3 temp_v_depth_normal;
    float3 temp_v_ref_normal;
    float temp_v_render_alpha;
    float3 temp_v_rgb_dist;
    float temp_v_depth_dist;
    float3 temp_v_normal_dist;

    per_pixel_losses_bwd(
        render_rgb ? render_rgb[idx] : make_float3(0),
        ref_rgb ? ref_rgb[idx] : make_float3(0),
        render_depth ? render_depth[idx] : 1.f,
        ref_depth ? ref_depth[idx] : 1.f,
        render_normal ? render_normal[idx] : make_float3(0),
        depth_normal ? depth_normal[idx] : make_float3(0),
        ref_normal ? ref_normal[idx] : make_float3(0),
        render_alpha ? render_alpha[idx] : 0.f,
        rgb_dist ? rgb_dist[idx] : make_float3(0),
        depth_dist ? depth_dist[idx] : 0.f,
        normal_dist ? normal_dist[idx] : make_float3(0),
        ref_alpha ? ref_alpha[idx] : true,
        mask ? mask[idx] : true,
        depth_mask ? depth_mask[idx] : true,
        normal_mask ? normal_mask[idx] : true,
        &v_losses,
        &temp_v_render_rgb,
        &temp_v_ref_rgb,
        &temp_v_render_depth,
        &temp_v_ref_depth,
        &temp_v_render_normal,
        &temp_v_depth_normal,
        &temp_v_ref_normal,
        &temp_v_render_alpha,
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
    if (v_render_alpha) v_render_alpha[idx] = temp_v_render_alpha;
    if (v_rgb_dist) v_rgb_dist[idx] = temp_v_rgb_dist;
    if (v_depth_dist) v_depth_dist[idx] = temp_v_depth_dist;
    if (v_normal_dist) v_normal_dist[idx] = temp_v_normal_dist;
}

__global__ void per_pixel_losses_reduce_forward_kernel(
    FixedArray<float, (uint)RawLossIndex::length>* __restrict__ raw_losses,
    size_t num_pixels,
    FixedArray<float, (uint)LossWeightIndex::length> loss_weights,
    FixedArray<float, (uint)LossIndex::length>* __restrict__ losses
) {
    per_pixel_losses_reduce(
        raw_losses, num_pixels, &loss_weights,
        losses
    );
}

__global__ void per_pixel_losses_reduce_backward_kernel(
    FixedArray<float, (uint)RawLossIndex::length>* __restrict__ raw_losses,
    size_t num_pixels,
    FixedArray<float, (uint)LossWeightIndex::length> loss_weights,
    FixedArray<float, (uint)LossIndex::length>* __restrict__ v_losses,
    FixedArray<float, (uint)RawLossIndex::length>* __restrict__ v_raw_losses
) {
    per_pixel_losses_reduce_bwd(
        raw_losses, num_pixels, &loss_weights,
        v_losses, v_raw_losses
    );
}


std::tuple<torch::Tensor, torch::Tensor>
compute_per_pixel_losses_forward_tensor(
    std::optional<at::Tensor> render_rgb,
    std::optional<at::Tensor> ref_rgb,
    std::optional<at::Tensor> render_depth,
    std::optional<at::Tensor> ref_depth,
    std::optional<at::Tensor> render_normal,
    std::optional<at::Tensor> depth_normal,
    std::optional<at::Tensor> ref_normal,
    std::optional<at::Tensor> render_alpha,
    std::optional<at::Tensor> rgb_dist,
    std::optional<at::Tensor> depth_dist,
    std::optional<at::Tensor> normal_dist,
    std::optional<at::Tensor> ref_alpha,
    std::optional<at::Tensor> mask,
    std::optional<at::Tensor> depth_mask,
    std::optional<at::Tensor> normal_mask,
    const std::array<float, (uint)LossIndex::length> loss_weights_0
) {
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
    check_float("render_alpha", render_alpha, 1);
    check_float("rgb_dist", rgb_dist, 3);
    check_float("depth_dist", depth_dist, 1);
    check_float("normal_dist", normal_dist, 3);
    check_float("ref_alpha", ref_alpha, 1);
    check_float("mask", mask, 1);
    check_float("depth_mask", depth_mask, 1);
    check_float("normal_mask", normal_mask, 1);

    size_t num_pixels = render_rgb.value().numel() / 3;

    FixedArray<float, (uint)LossWeightIndex::length> loss_weights =
        *reinterpret_cast<const FixedArray<float, (uint)LossWeightIndex::length>*>(loss_weights_0.data());

    torch::Tensor raw_losses = torch::zeros({(uint)RawLossIndex::length}, render_rgb.value().options());
    torch::Tensor losses = torch::zeros({(uint)LossIndex::length}, render_rgb.value().options());

    per_pixel_losses_forward_kernel<<<_LAUNCH_ARGS_1D(num_pixels, WARP_SIZE*WARP_SIZE)>>>(
        num_pixels,
        render_rgb.has_value() ? (float3*)render_rgb.value().contiguous().data_ptr<float>() : nullptr,
        ref_rgb.has_value() ? (float3*)ref_rgb.value().contiguous().data_ptr<float>() : nullptr,
        render_depth.has_value() ? (float*)render_depth.value().contiguous().data_ptr<float>() : nullptr,
        ref_depth.has_value() ? (float*)ref_depth.value().contiguous().data_ptr<float>() : nullptr,
        render_normal.has_value() ? (float3*)render_normal.value().contiguous().data_ptr<float>() : nullptr,
        depth_normal.has_value() ? (float3*)depth_normal.value().contiguous().data_ptr<float>() : nullptr,
        ref_normal.has_value() ? (float3*)ref_normal.value().contiguous().data_ptr<float>() : nullptr,
        render_alpha.has_value() ? (float*)render_alpha.value().contiguous().data_ptr<float>() : nullptr,
        rgb_dist.has_value() ? (float3*)rgb_dist.value().contiguous().data_ptr<float>() : nullptr,
        depth_dist.has_value() ? (float*)depth_dist.value().contiguous().data_ptr<float>() : nullptr,
        normal_dist.has_value() ? (float3*)normal_dist.value().contiguous().data_ptr<float>() : nullptr,
        ref_alpha.has_value() ? (bool*)ref_alpha.value().contiguous().data_ptr<bool>() : nullptr,
        mask.has_value() ? (bool*)mask.value().contiguous().data_ptr<bool>() : nullptr,
        depth_mask.has_value() ? (bool*)depth_mask.value().contiguous().data_ptr<bool>() : nullptr,
        normal_mask.has_value() ? (bool*)normal_mask.value().contiguous().data_ptr<bool>() : nullptr,
        raw_losses.data_ptr<float>()
    );

    per_pixel_losses_reduce_forward_kernel<<<1, 1>>>(
        (FixedArray<float, (uint)RawLossIndex::length>*)raw_losses.data_ptr<float>(),
        num_pixels,
        loss_weights,
        (FixedArray<float, (uint)LossIndex::length>*)losses.data_ptr<float>()
    );

    return std::make_tuple(losses, raw_losses);
}


std::tuple<
    std::optional<at::Tensor>, // render_rgb
    std::optional<at::Tensor>, // ref_rgb
    std::optional<at::Tensor>, // render_depth
    std::optional<at::Tensor>, // ref_depth
    std::optional<at::Tensor>, // render_normal
    std::optional<at::Tensor>, // depth_normal
    std::optional<at::Tensor>, // ref_normal
    std::optional<at::Tensor>, // render_alpha
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
    std::optional<at::Tensor> render_alpha,
    std::optional<at::Tensor> rgb_dist,
    std::optional<at::Tensor> depth_dist,
    std::optional<at::Tensor> normal_dist,
    std::optional<at::Tensor> ref_alpha,
    std::optional<at::Tensor> mask,
    std::optional<at::Tensor> depth_mask,
    std::optional<at::Tensor> normal_mask,
    at::Tensor raw_losses,
    const std::array<float, (uint)LossIndex::length> loss_weights_0,
    at::Tensor v_losses
) {

    size_t num_pixels = render_rgb.value().numel() / 3;

    std::optional<at::Tensor> v_render_rgb = render_rgb.has_value() ? (std::optional<at::Tensor>)at::empty_like(render_rgb.value()) : std::nullopt;
    std::optional<at::Tensor> v_ref_rgb = ref_rgb.has_value() ? (std::optional<at::Tensor>)at::empty_like(ref_rgb.value()) : std::nullopt;
    std::optional<at::Tensor> v_render_depth = render_depth.has_value() ? (std::optional<at::Tensor>)at::empty_like(render_depth.value()) : std::nullopt;
    std::optional<at::Tensor> v_ref_depth = ref_depth.has_value() ? (std::optional<at::Tensor>)at::empty_like(ref_depth.value()) : std::nullopt;
    std::optional<at::Tensor> v_render_normal = render_normal.has_value() ? (std::optional<at::Tensor>)at::empty_like(render_normal.value()) : std::nullopt;
    std::optional<at::Tensor> v_depth_normal = depth_normal.has_value() ? (std::optional<at::Tensor>)at::empty_like(depth_normal.value()) : std::nullopt;
    std::optional<at::Tensor> v_ref_normal = ref_normal.has_value() ? (std::optional<at::Tensor>)at::empty_like(ref_normal.value()) : std::nullopt;
    std::optional<at::Tensor> v_render_alpha = render_alpha.has_value() ? (std::optional<at::Tensor>)at::empty_like(render_alpha.value()) : std::nullopt;
    std::optional<at::Tensor> v_rgb_dist = rgb_dist.has_value() ? (std::optional<at::Tensor>)at::empty_like(rgb_dist.value()) : std::nullopt;
    std::optional<at::Tensor> v_depth_dist = depth_dist.has_value() ? (std::optional<at::Tensor>)at::empty_like(depth_dist.value()) : std::nullopt;
    std::optional<at::Tensor> v_normal_dist = normal_dist.has_value() ? (std::optional<at::Tensor>)at::empty_like(normal_dist.value()) : std::nullopt;

    FixedArray<float, (uint)LossWeightIndex::length> loss_weights =
        *reinterpret_cast<const FixedArray<float, (uint)LossWeightIndex::length>*>(loss_weights_0.data());

    torch::Tensor v_raw_losses = torch::empty({(uint)RawLossIndex::length}, render_rgb.value().options());

    per_pixel_losses_reduce_backward_kernel<<<1, 1>>>(
        (FixedArray<float, (uint)RawLossIndex::length>*)raw_losses.data_ptr<float>(),
        num_pixels,
        loss_weights,
        (FixedArray<float, (uint)LossIndex::length>*)v_losses.data_ptr<float>(),
        (FixedArray<float, (uint)RawLossIndex::length>*)v_raw_losses.data_ptr<float>()
    );

    per_pixel_losses_backward_kernel<<<_LAUNCH_ARGS_1D(num_pixels, WARP_SIZE*WARP_SIZE)>>>(
        num_pixels,
        render_rgb.has_value() ? (float3*)render_rgb.value().contiguous().data_ptr<float>() : nullptr,
        ref_rgb.has_value() ? (float3*)ref_rgb.value().contiguous().data_ptr<float>() : nullptr,
        render_depth.has_value() ? (float*)render_depth.value().contiguous().data_ptr<float>() : nullptr,
        ref_depth.has_value() ? (float*)ref_depth.value().contiguous().data_ptr<float>() : nullptr,
        render_normal.has_value() ? (float3*)render_normal.value().contiguous().data_ptr<float>() : nullptr,
        depth_normal.has_value() ? (float3*)depth_normal.value().contiguous().data_ptr<float>() : nullptr,
        ref_normal.has_value() ? (float3*)ref_normal.value().contiguous().data_ptr<float>() : nullptr,
        render_alpha.has_value() ? (float*)render_alpha.value().contiguous().data_ptr<float>() : nullptr,
        rgb_dist.has_value() ? (float3*)rgb_dist.value().contiguous().data_ptr<float>() : nullptr,
        depth_dist.has_value() ? (float*)depth_dist.value().contiguous().data_ptr<float>() : nullptr,
        normal_dist.has_value() ? (float3*)normal_dist.value().contiguous().data_ptr<float>() : nullptr,
        ref_alpha.has_value() ? (bool*)ref_alpha.value().contiguous().data_ptr<bool>() : nullptr,
        mask.has_value() ? (bool*)mask.value().contiguous().data_ptr<bool>() : nullptr,
        depth_mask.has_value() ? (bool*)depth_mask.value().contiguous().data_ptr<bool>() : nullptr,
        normal_mask.has_value() ? (bool*)normal_mask.value().contiguous().data_ptr<bool>() : nullptr,
        v_raw_losses.data_ptr<float>(),
        v_render_rgb.has_value() ? (float3*)v_render_rgb.value().contiguous().data_ptr<float>() : nullptr,
        v_ref_rgb.has_value() ? (float3*)v_ref_rgb.value().contiguous().data_ptr<float>() : nullptr,
        v_render_depth.has_value() ? (float*)v_render_depth.value().contiguous().data_ptr<float>() : nullptr,
        v_ref_depth.has_value() ? (float*)v_ref_depth.value().contiguous().data_ptr<float>() : nullptr,
        v_render_normal.has_value() ? (float3*)v_render_normal.value().contiguous().data_ptr<float>() : nullptr,
        v_depth_normal.has_value() ? (float3*)v_depth_normal.value().contiguous().data_ptr<float>() : nullptr,
        v_ref_normal.has_value() ? (float3*)v_ref_normal.value().contiguous().data_ptr<float>() : nullptr,
        v_render_alpha.has_value() ? (float*)v_render_alpha.value().contiguous().data_ptr<float>() : nullptr,
        v_rgb_dist.has_value() ? (float3*)v_rgb_dist.value().contiguous().data_ptr<float>() : nullptr,
        v_depth_dist.has_value() ? (float*)v_depth_dist.value().contiguous().data_ptr<float>() : nullptr,
        v_normal_dist.has_value() ? (float3*)v_normal_dist.value().contiguous().data_ptr<float>() : nullptr
    );

    return std::make_tuple(
        v_render_rgb,
        v_ref_rgb,
        v_render_depth,
        v_ref_depth,
        v_render_normal,
        v_depth_normal,
        v_ref_normal,
        v_render_alpha,
        v_rgb_dist,
        v_depth_dist,
        v_normal_dist
    );
}


