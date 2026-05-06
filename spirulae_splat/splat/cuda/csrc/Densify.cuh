#pragma once

#include <ATen/Tensor.h>

#include "types.cuh"


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



at::Tensor quantile_of_abs_of_finite_elements_tensor(
    at::Tensor inputs,
    float q,
    bool return_reciprocal
);


void inplace_index_tensor(
    at::Tensor indices,
    at::Tensor src,
    at::Tensor dst
);


void inplace_scatter_add_tensor(
    at::Tensor indices,
    at::Tensor src,
    at::Tensor dst
);


void inplace_scatter_max_tensor(
    at::Tensor indices,
    at::Tensor src,
    at::Tensor dst
);


at::Tensor weighted_sample_without_replacement_tensor(
    int64_t numel,
    at::Tensor weights,
    std::optional<at::Tensor> masks,
    uint32_t num_sample,
    uint32_t seed
);


at::Tensor cov_scale_init_tensor(
    at::Tensor points,  // [N, 3]
    at::Tensor is_fisheye,  // [C], bool
    at::Tensor sizes,  // [C, 2], int32
    at::Tensor intrins,  // [C, 4]
    at::Tensor viewmats,  // [C, 4, 4]
    CameraDistortionCoeffsTensor dist_coeffs // [C]
);


void densify_clip_scale_tensor(
    int64_t num_splats,
    at::Tensor radii,
    at::Tensor log_scales,
    std::optional<at::Tensor> logit_opacs,
    float max_scale2d,
    float clip_hardness,
    float max_scale3d
);


void densify_update_weight_tensor(
    int64_t num_splats,
    at::Tensor radii,
    at::Tensor opacs,
    at::Tensor accum_weight,
    at::Tensor accum_buffer,
    bool is_max_mode
);


void relocate_splats_with_long_axis_split_tensor(
    int64_t cur_num_splats,
    float min_opacity,
    at::Tensor means, at::Tensor quats, at::Tensor scales, at::Tensor opacs, at::Tensor features_dc, at::Tensor features_sh,
    at::Tensor g1_means, at::Tensor g1_quats, at::Tensor g1_scales, at::Tensor g1_opacs, at::Tensor g1_features_dc, at::Tensor g1_features_sh,
    at::Tensor g2_means, at::Tensor g2_quats, at::Tensor g2_scales, at::Tensor g2_opacs, at::Tensor g2_features_dc, at::Tensor g2_features_sh,
    at::Tensor densify_accum_buffer,
    uint32_t seed
);


void add_splats_with_long_axis_split_tensor(
    int64_t cur_num_splats,
    int64_t num_new_splats,
    at::Tensor means, at::Tensor quats, at::Tensor scales, at::Tensor opacs, at::Tensor features_dc, at::Tensor features_sh,
    at::Tensor g1_means, at::Tensor g1_quats, at::Tensor g1_scales, at::Tensor g1_opacs, at::Tensor g1_features_dc, at::Tensor g1_features_sh,
    at::Tensor g2_means, at::Tensor g2_quats, at::Tensor g2_scales, at::Tensor g2_opacs, at::Tensor g2_features_dc, at::Tensor g2_features_sh,
    at::Tensor densify_accum_buffer,
    uint32_t seed
);


void mcmc_add_noise_tensor(
    std::string primitive,
    float scaler, float min_opacity,
    at::Tensor &means,
    at::Tensor &log_scales,
    at::Tensor &quats,
    at::Tensor &opacs
);


std::tuple<at::Tensor, at::Tensor>
compute_relocation_tensor(
    // inputs
    at::Tensor opacities, // [N]
    at::Tensor scales,    // [N, 3]
    at::Tensor ratios,    // [N]
    at::Tensor binoms,    // [n_max, n_max]
    const int n_max
);


std::tuple<at::Tensor, at::Tensor, at::Tensor>
long_axis_split_tensor(
    std::string primitive,
    at::Tensor &log_scales,
    at::Tensor &logit_opacities,
    at::Tensor &quats
);


at::Tensor laplacian_edge_filter_tensor(
    at::Tensor &img_in
);


at::Tensor smoothed_laplacian_edge_filter_tensor(
    at::Tensor &img_in
);


at::Tensor canny_edge_filter_tensor(
    at::Tensor &img_in
);
