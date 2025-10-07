#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>
namespace cg = cooperative_groups;

#define TensorView _Slang_TensorView
#include "generated/slang_all.cu"
#undef TensorView

#include "misc.cuh"


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


__global__ void blend_background_forward_kernel(
    const TensorView<float, 3> in_rgb,
    const TensorView<float, 3> in_alpha,
    const TensorView<float, 3> in_background,
    TensorView<float, 3> out_rgb
) {
    unsigned gid = blockIdx.x * blockDim.x + threadIdx.x;
    if (gid >= in_rgb.shape[0]*in_rgb.shape[1])
        return;
    unsigned y = gid / in_rgb.shape[1];
    unsigned x = gid % in_rgb.shape[1];

    float3 rgb = in_rgb.load3f(y, x);
    float alpha = in_alpha.load1f(y, x);
    float3 background = in_background.load3f(y, x);

    rgb = blend_background(rgb, alpha, background);

    out_rgb.store3f(y, x, rgb);
}


__global__ void blend_background_backward_kernel(
    const TensorView<float, 3> in_rgb,
    const TensorView<float, 3> in_alpha,
    const TensorView<float, 3> in_background,
    const TensorView<float, 3> v_out_rgb,
    TensorView<float, 3> v_in_rgb,
    TensorView<float, 3> v_in_alpha,
    TensorView<float, 3> v_in_background
) {
    unsigned gid = blockIdx.x * blockDim.x + threadIdx.x;
    if (gid >= in_rgb.shape[0]*in_rgb.shape[1])
        return;
    unsigned y = gid / in_rgb.shape[1];
    unsigned x = gid % in_rgb.shape[1];

    float3 rgb = in_rgb.load3f(y, x);
    float alpha = in_alpha.load1f(y, x);
    float3 background = in_background.load3f(y, x);

    float3 v_out = v_out_rgb.load3f(y, x);

    float3 v_rgb; float v_alpha; float3 v_background;
    blend_background_bwd(
        rgb, alpha, background,
        v_out,
        &v_rgb, &v_alpha, &v_background
    );

    v_in_rgb.store3f(y, x, v_rgb);
    v_in_alpha.store1f(y, x, v_alpha);
    v_in_background.store3f(y, x, v_background);

}
