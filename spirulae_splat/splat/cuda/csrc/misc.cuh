#include "utils.cuh"

static constexpr uint kNumPerSplatLosses = 5;


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



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
);


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
);


__global__ void blend_background_forward_kernel(
    const TensorView<float, 3> in_rgb,
    const TensorView<float, 3> in_alpha,
    const TensorView<float, 3> in_background,
    TensorView<float, 3> out_rgb
);


__global__ void blend_background_backward_kernel(
    const TensorView<float, 3> in_rgb,
    const TensorView<float, 3> in_alpha,
    const TensorView<float, 3> in_background,
    const TensorView<float, 3> v_out_rgb,
    TensorView<float, 3> v_in_rgb,
    TensorView<float, 3> v_in_alpha,
    TensorView<float, 3> v_in_background
);
