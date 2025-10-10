#include "PerSplatLoss.cuh"

#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>
namespace cg = cooperative_groups;

#define TensorView _Slang_TensorView
#include "generated/slang_all.cu"
#undef TensorView

#include "common.cuh"

__global__ void per_splat_losses_forward_kernel(
    bool is_3dgs,
    const size_t num_points,
    const float* __restrict__ scales_buffer,
    const float* __restrict__ opacities_buffer,
    const float* __restrict__ quats_buffer,
    float* __restrict__ out_losses,
    float max_gauss_ratio,
    float scale_regularization_weight,
    float mcmc_opacity_reg_weight,
    float mcmc_scale_reg_weight,
    float erank_reg_weight,
    float erank_reg_weight_s3,
    float quat_norm_reg_weight
) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;

    FixedArray<float, kNumPerSplatLosses> losses;

    bool inside = idx < num_points;
    if (inside) {
        float3 scale;
        if (is_3dgs) scale = { scales_buffer[3*idx+0], scales_buffer[3*idx+1], scales_buffer[3*idx+2] };
        else scale = { scales_buffer[2*idx+0], scales_buffer[2*idx+1], 0.0f };
        float opacity = opacities_buffer[idx];
        float4 quat = { quats_buffer[4*idx+0], quats_buffer[4*idx+1], quats_buffer[4*idx+2], quats_buffer[4*idx+3] };
        per_splat_losses(
            is_3dgs,
            scale, opacity, quat,
            max_gauss_ratio,
            scale_regularization_weight,
            mcmc_opacity_reg_weight,
            mcmc_scale_reg_weight,
            erank_reg_weight,
            erank_reg_weight_s3,
            quat_norm_reg_weight,
            &losses
        );
    }

    auto block = cg::this_thread_block();
    cg::thread_block_tile<32> warp = cg::tiled_partition<32>(block);
    for (int i = 0; i < kNumPerSplatLosses; i++) {
        float loss = inside ? losses[i] : 0.0f;
        float loss_reduced = cg::reduce(warp, loss, cg::plus<float>());
        loss_reduced = loss_reduced / (float)num_points;
        if (warp.thread_rank() == 0) {
            atomicAdd(&out_losses[i], loss_reduced);
        }
    }
}


__global__ void per_splat_losses_backward_kernel(
    bool is_3dgs,
    const size_t num_points,
    const float* __restrict__ scales_buffer,
    const float* __restrict__ opacities_buffer,
    const float* __restrict__ quats_buffer,
    const float* __restrict__ v_out_losses,
    float* __restrict__ v_scales_buffer,
    float* __restrict__ v_opacities_buffer,
    float* __restrict__ v_quats_buffer,
    float max_gauss_ratio,
    float scale_regularization_weight,
    float mcmc_opacity_reg_weight,
    float mcmc_scale_reg_weight,
    float erank_reg_weight,
    float erank_reg_weight_s3,
    float quat_norm_reg_weight
) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;

    bool inside = idx < num_points;
    if (!inside) return;

    FixedArray<float, kNumPerSplatLosses> v_losses;
    for (int i = 0; i < kNumPerSplatLosses; i++)
        v_losses[i] = v_out_losses[i] / (float)num_points;

    float3 scale;
    if (is_3dgs) scale = { scales_buffer[3*idx+0], scales_buffer[3*idx+1], scales_buffer[3*idx+2] };
    else scale = { scales_buffer[2*idx+0], scales_buffer[2*idx+1], 0.0f };
    float opacity = opacities_buffer[idx];
    float4 quat = { quats_buffer[4*idx+0], quats_buffer[4*idx+1], quats_buffer[4*idx+2], quats_buffer[4*idx+3] };

    float3 v_scale;
    float v_opacity;
    float4 v_quat;

    per_splat_losses_bwd(
        is_3dgs,
        scale, opacity, quat,
        v_losses,
        &v_scale, &v_opacity, &v_quat,
        max_gauss_ratio,
        scale_regularization_weight,
        mcmc_opacity_reg_weight,
        mcmc_scale_reg_weight,
        erank_reg_weight,
        erank_reg_weight_s3,
        quat_norm_reg_weight
    );

    if (is_3dgs) v_scales_buffer[3*idx+0] = v_scale.x, v_scales_buffer[3*idx+1] = v_scale.y, v_scales_buffer[3*idx+2] = v_scale.z;
    else v_scales_buffer[2*idx+0] = v_scale.x, v_scales_buffer[2*idx+1] = v_scale.y;
    v_opacities_buffer[idx] = v_opacity;
    v_quats_buffer[4*idx+0] = v_quat.x, v_quats_buffer[4*idx+1] = v_quat.y, v_quats_buffer[4*idx+2] = v_quat.z, v_quats_buffer[4*idx+3] = v_quat.w;
}



torch::Tensor compute_per_splat_losses_forward_tensor(
    torch::Tensor &scales,  // [N, 3] or [N, 2]
    torch::Tensor &opacities,  // [N, 1]
    torch::Tensor &quats,  // [N, 4]
    float mcmc_opacity_reg_weight,
    float mcmc_scale_reg_weight,
    float max_gauss_ratio,
    float scale_regularization_weight,
    float erank_reg_weight,
    float erank_reg_weight_s3,
    float quat_norm_reg_weight
) {
    DEVICE_GUARD(scales);
    CHECK_INPUT(scales);
    CHECK_INPUT(opacities);
    CHECK_INPUT(quats);

    const size_t num_points = opacities.size(0);
    const bool is_3dgs = (scales.size(-1) == 3);

    if (scales.ndimension() != 2 || scales.size(0) != num_points ||
        (scales.size(1) != 3 && scales.size(1) != 2))
        AT_ERROR("scales shape must be (n, 2) or (n, 3)");
    if (opacities.ndimension() != 2 || opacities.size(0) != num_points || opacities.size(1) != 1)
        AT_ERROR("opacities shape must be (n, 1)");
    if (quats.ndimension() != 2 || quats.size(0) != num_points || quats.size(1) != 4)
        AT_ERROR("quats shape must be (n, 4)");

    torch::Tensor loss = torch::zeros({kNumPerSplatLosses}, opacities.options());

    per_splat_losses_forward_kernel<<<_LAUNGH_ARGS_1D(num_points)>>>(
        is_3dgs,
        num_points,
        scales.contiguous().data_ptr<float>(),
        opacities.contiguous().data_ptr<float>(),
        quats.contiguous().data_ptr<float>(),
        loss.data_ptr<float>(),
        mcmc_opacity_reg_weight,
        mcmc_scale_reg_weight,
        max_gauss_ratio,
        scale_regularization_weight,
        erank_reg_weight,
        erank_reg_weight_s3,
        quat_norm_reg_weight
    );

    return loss;
}


std::tuple<torch::Tensor, torch::Tensor, torch::Tensor>
compute_per_splat_losses_backward_tensor(
    torch::Tensor &scales,  // [N, 3] or [N, 2]
    torch::Tensor &opacities,  // [N, 1]
    torch::Tensor &quats,  // [N, 4]
    torch::Tensor &v_losses,  // [kNumPerSplatLosses]
    float mcmc_opacity_reg_weight,
    float mcmc_scale_reg_weight,
    float max_gauss_ratio,
    float scale_regularization_weight,
    float erank_reg_weight,
    float erank_reg_weight_s3,
    float quat_norm_reg_weight
) {
    DEVICE_GUARD(scales);
    CHECK_INPUT(scales);
    CHECK_INPUT(opacities);
    CHECK_INPUT(quats);
    CHECK_INPUT(v_losses);

    const size_t num_points = opacities.size(0);
    const bool is_3dgs = (scales.size(-1) == 3);

    if (scales.ndimension() != 2 || scales.size(0) != num_points ||
        (scales.size(1) != 3 && scales.size(1) != 2))
        AT_ERROR("scales shape must be (n, 2) or (n, 3)");
    if (opacities.ndimension() != 2 || opacities.size(0) != num_points || opacities.size(1) != 1)
        AT_ERROR("opacities shape must be (n, 1)");
    if (quats.ndimension() != 2 || quats.size(0) != num_points || quats.size(1) != 4)
        AT_ERROR("quats shape must be (n, 4)");
    if (v_losses.ndimension() != 1 || v_losses.size(0) != kNumPerSplatLosses)
        AT_ERROR("v_losses shape must be (kNumPerSplatLosses,)");

    torch::Tensor v_scales = torch::empty_like(scales);
    torch::Tensor v_opacities = torch::empty_like(opacities);
    torch::Tensor v_quats = torch::empty_like(quats);

    per_splat_losses_backward_kernel<<<_LAUNGH_ARGS_1D(num_points)>>>(
        is_3dgs,
        num_points,
        scales.contiguous().data_ptr<float>(),
        opacities.contiguous().data_ptr<float>(),
        quats.contiguous().data_ptr<float>(),
        v_losses.data_ptr<float>(),
        v_scales.contiguous().data_ptr<float>(),
        v_opacities.contiguous().data_ptr<float>(),
        v_quats.contiguous().data_ptr<float>(),
        mcmc_opacity_reg_weight,
        mcmc_scale_reg_weight,
        max_gauss_ratio,
        scale_regularization_weight,
        erank_reg_weight,
        erank_reg_weight_s3,
        quat_norm_reg_weight
    );

    return std::make_tuple(v_scales, v_opacities, v_quats);
}


