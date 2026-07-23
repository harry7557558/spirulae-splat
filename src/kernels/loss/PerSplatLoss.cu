#include "kernels/loss/PerSplatLoss.cuh"

#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>
namespace cg = cooperative_groups;

#include "generated/slang.cuh"
namespace SlangPerSplatLosses {
#include "generated/set_namespace.cuh"
#include "generated/per_splat_losses.cuh"
}

#if 0

#include "core/Common.cuh"

__global__ void per_splat_losses_forward_kernel(
    const size_t num_points,
    const float3* __restrict__ scales_buffer,
    const float* __restrict__ opacities_buffer,
    const float4* __restrict__ quats_buffer,
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
        float3 scale = scales_buffer[idx];
        float opacity = opacities_buffer[idx];
        float4 quat = quats_buffer[idx];
        SlangPerSplatLosses::per_splat_losses(
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
    cg::thread_block_tile<WARP_SIZE> warp = cg::tiled_partition<WARP_SIZE>(block);
    for (int i = 0; i < kNumPerSplatLosses; i++) {
        float loss = inside ? losses[i] : 0.0f;
        if (!isfinite(loss)) loss = 0.0f;
        float loss_reduced = cg::reduce(warp, loss, cg::plus<float>());
        loss_reduced = loss_reduced / (float)num_points;
        if (warp.thread_rank() == 0 && loss_reduced != 0.0f && isfinite(loss_reduced)) {
            atomicAdd(&out_losses[i], loss_reduced);
        }
    }
}


__global__ void per_splat_losses_backward_kernel(
    const size_t num_points,
    const float3* __restrict__ scales_buffer,
    const float* __restrict__ opacities_buffer,
    const float4* __restrict__ quats_buffer,
    const float* __restrict__ v_out_losses,
    float3* __restrict__ v_scales_buffer,
    float* __restrict__ v_opacities_buffer,
    float4* __restrict__ v_quats_buffer,
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

    float3 scale = scales_buffer[idx];
    float opacity = opacities_buffer[idx];
    float4 quat = quats_buffer[idx];

    float3 v_scale;
    float v_opacity;
    float4 v_quat;

    SlangPerSplatLosses::per_splat_losses_bwd(
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

    v_scales_buffer[idx] = isfinite(dot(v_scale, v_scale)) ? v_scale : make_float3(0.0f);
    v_opacities_buffer[idx] = isfinite(v_opacity) ? v_opacity : 0.0f;
    v_quats_buffer[idx] = isfinite(dot(v_quat, v_quat)) ? v_quat : make_float4(0.0f);
}


/*[AutoHeaderGeneratorExport]*/
at::Tensor compute_per_splat_losses_forward(
    at::Tensor &scales,  // [N, 3] or [N, 2]
    at::Tensor &opacities,  // [N, 1]
    at::Tensor &quats,  // [N, 4]
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

    if (scales.ndimension() != 2 || scales.size(0) != num_points || scales.size(1) != 3)
        AT_ERROR("scales shape must be (n, 3)");
    if (opacities.ndimension() != 2 || opacities.size(0) != num_points || opacities.size(1) != 1)
        AT_ERROR("opacities shape must be (n, 1)");
    if (quats.ndimension() != 2 || quats.size(0) != num_points || quats.size(1) != 4)
        AT_ERROR("quats shape must be (n, 4)");

    at::Tensor loss = at::zeros({kNumPerSplatLosses}, opacities.options());

    per_splat_losses_forward_kernel<<<_LAUNCH_ARGS_1D(num_points, 256)>>>(
        num_points,
        (float3*)scales.contiguous().data_ptr<float>(),
        opacities.contiguous().data_ptr<float>(),
        (float4*)quats.contiguous().data_ptr<float>(),
        loss.data_ptr<float>(),
        mcmc_opacity_reg_weight,
        mcmc_scale_reg_weight,
        max_gauss_ratio,
        scale_regularization_weight,
        erank_reg_weight,
        erank_reg_weight_s3,
        quat_norm_reg_weight
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return loss;
}


/*[AutoHeaderGeneratorExport]*/
std::tuple<at::Tensor, at::Tensor, at::Tensor>
compute_per_splat_losses_backward(
    at::Tensor &scales,  // [N, 3] or [N, 2]
    at::Tensor &opacities,  // [N, 1]
    at::Tensor &quats,  // [N, 4]
    at::Tensor &v_losses,  // [kNumPerSplatLosses]
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

    if (scales.ndimension() != 2 || scales.size(0) != num_points || scales.size(1) != 3)
        AT_ERROR("scales shape must be (n, 3)");
    if (opacities.ndimension() != 2 || opacities.size(0) != num_points || opacities.size(1) != 1)
        AT_ERROR("opacities shape must be (n, 1)");
    if (quats.ndimension() != 2 || quats.size(0) != num_points || quats.size(1) != 4)
        AT_ERROR("quats shape must be (n, 4)");
    if (v_losses.ndimension() != 1 || v_losses.size(0) != kNumPerSplatLosses)
        AT_ERROR("v_losses shape must be (kNumPerSplatLosses,)");

    at::Tensor v_scales = at::empty_like(scales);
    at::Tensor v_opacities = at::empty_like(opacities);
    at::Tensor v_quats = at::empty_like(quats);

    per_splat_losses_backward_kernel<<<_LAUNCH_ARGS_1D(num_points, 256)>>>(
        num_points,
        (float3*)scales.contiguous().data_ptr<float>(),
        opacities.contiguous().data_ptr<float>(),
        (float4*)quats.contiguous().data_ptr<float>(),
        v_losses.data_ptr<float>(),
        (float3*)v_scales.contiguous().data_ptr<float>(),
        v_opacities.contiguous().data_ptr<float>(),
        (float4*)v_quats.contiguous().data_ptr<float>(),
        mcmc_opacity_reg_weight,
        mcmc_scale_reg_weight,
        max_gauss_ratio,
        scale_regularization_weight,
        erank_reg_weight,
        erank_reg_weight_s3,
        quat_norm_reg_weight
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return std::make_tuple(v_scales, v_opacities, v_quats);
}

#endif
