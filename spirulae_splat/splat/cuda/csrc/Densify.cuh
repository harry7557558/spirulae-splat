#pragma once


#include <Tensor.h>


// densify_update_weight score mode. Mirrors DensifyConfig::score_mode.
enum class DensifyScoreMode : int {
    Mean   = 0,
    Max    = 1,
    Median = 2,
    Geom   = 3,
};


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



void quantile_of_abs_of_finite_elements_tensor(
    DeviceVector<float> inputs,
    float q,
    bool return_reciprocal,
    DeviceVector<float> outputs
);


void normalize_by_median_inplace_tensor(
    DeviceVector<float> data
);


void inplace_index_tensor(
    DeviceVector<int32_t> indices,
    DeviceVector<float> src,
    DeviceVector<float> dst
);


void inplace_scatter_add_tensor(
    DeviceVector<int32_t> indices,
    DeviceVector<float> src,
    DeviceVector<float> dst
);


void inplace_scatter_max_tensor(
    DeviceVector<int32_t> indices,
    DeviceVector<float> src,
    DeviceVector<float> dst
);


void weighted_sample_without_replacement_tensor(
    int64_t numel,
    DeviceVector<float> weights,
    bool* masks_ptr,
    uint32_t num_sample,
    uint32_t seed,
    DeviceVector<int32_t> out_idx
);


void cov_scale_init_tensor(
    DeviceVector<float3> points,  // [N, 3]
    DeviceVector<bool> is_fisheye,  // [C], bool
    DeviceVector<int2> sizes,  // [C, 2], int32
    DeviceVector<float4> intrins,  // [C, 4]
    DeviceVector<float4> viewmats,  // [C, 4, 4] as 4*C float4 elements
    TorchTensorView dist_coeffs, // [C]
    DeviceVector<float> log_scales  // [N, 1] output
);


void densify_clip_scale_tensor(
    int64_t num_splats,
    DeviceVector<float> radii,
    DeviceVector<float3> log_scales,
    float* logit_opacs_ptr,
    float max_scale2d,
    float clip_hardness,
    float max_scale3d
);


void densify_update_weight(
    int64_t num_splats,
    DeviceVector<float> radii,
    float3* scales_ptr,
    float* opacs_ptr,
    DeviceVector<float> accum_weight,
    DeviceVector<float2> accum_buffer,
    int score_mode
);


void relocate_splats_with_long_axis_split_tensor(
    int64_t cur_num_splats,
    float min_opacity,
    DeviceVector<float3> means, DeviceVector<float4> quats, DeviceVector<float3> scales, DeviceVector<float> opacs, DeviceVector<float3> features_dc, DeviceVector<float3> features_sh,
    DeviceVector<float3> g1_means, DeviceVector<float4> g1_quats, DeviceVector<float3> g1_scales, DeviceVector<float> g1_opacs, DeviceVector<float3> g1_features_dc, DeviceVector<float3> g1_features_sh,
    DeviceVector<float3> g2_means, DeviceVector<float4> g2_quats, DeviceVector<float3> g2_scales, DeviceVector<float> g2_opacs, DeviceVector<float3> g2_features_dc, DeviceVector<float3> g2_features_sh,
    DeviceVector<float2> densify_accum_buffer,
    DeviceVector<int32_t> bias_correction_steps,
    int sh_optim_bits,
    int num_sh,
    // SH-quant bounds buffer + layout flag used to encode (g1=0, g2=0) into
    // the dst splats' packed bytes. Null when sh_optim_bits == 32 (no quant).
    DeviceVector<float4> sh_quant_bounds,
    bool sh_bounds_per_splat,
    // SH VALUE-quant codec copy params. When sh_value_bits != 32 we also do
    // a codec-aware src->dst copy of the SH coefficient bytes (decode src,
    // encode dst against current dst bounds; clipping is acceptable -- see
    // _copy_quant_sh_value_for_splat for the rationale).
    DeviceVector<uint8_t> sh_value_packed,
    DeviceVector<float2>  sh_value_bounds,
    int  sh_value_bits,             // 32 / 8 / 16
    bool sh_value_bounds_per_splat,
    int  num_sh_buffer,             // buffer stride for cell indexing
    uint32_t seed
);


void add_splats_with_long_axis_split_tensor(
    int64_t cur_num_splats,
    int64_t num_new_splats,
    DeviceVector<float3> means, DeviceVector<float4> quats, DeviceVector<float3> scales, DeviceVector<float> opacs, DeviceVector<float3> features_dc, DeviceVector<float3> features_sh,
    DeviceVector<float3> g1_means, DeviceVector<float4> g1_quats, DeviceVector<float3> g1_scales, DeviceVector<float> g1_opacs, DeviceVector<float3> g1_features_dc, DeviceVector<float3> g1_features_sh,
    DeviceVector<float3> g2_means, DeviceVector<float4> g2_quats, DeviceVector<float3> g2_scales, DeviceVector<float> g2_opacs, DeviceVector<float3> g2_features_dc, DeviceVector<float3> g2_features_sh,
    DeviceVector<float2> densify_accum_buffer,
    DeviceVector<int32_t> bias_correction_steps,
    int sh_optim_bits,
    int num_sh,
    DeviceVector<float4> sh_quant_bounds,
    bool sh_bounds_per_splat,
    DeviceVector<uint8_t> sh_value_packed,
    DeviceVector<float2>  sh_value_bounds,
    int  sh_value_bits,
    bool sh_value_bounds_per_splat,
    int  num_sh_buffer,
    uint32_t seed
);


void relocate_splats_mcmc_tensor(
    int64_t cur_num_splats,
    float min_opacity,
    DeviceVector<float3> means, DeviceVector<float4> quats, DeviceVector<float3> scales, DeviceVector<float> opacs, DeviceVector<float3> features_dc, DeviceVector<float3> features_sh,
    DeviceVector<float3> g1_means, DeviceVector<float4> g1_quats, DeviceVector<float3> g1_scales, DeviceVector<float> g1_opacs, DeviceVector<float3> g1_features_dc, DeviceVector<float3> g1_features_sh,
    DeviceVector<float3> g2_means, DeviceVector<float4> g2_quats, DeviceVector<float3> g2_scales, DeviceVector<float> g2_opacs, DeviceVector<float3> g2_features_dc, DeviceVector<float3> g2_features_sh,
    DeviceVector<int32_t> bias_correction_steps,
    int sh_optim_bits,
    int num_sh,
    DeviceVector<float4> sh_quant_bounds,
    bool sh_bounds_per_splat,
    DeviceVector<uint8_t> sh_value_packed,
    DeviceVector<float2>  sh_value_bounds,
    int  sh_value_bits,
    bool sh_value_bounds_per_splat,
    int  num_sh_buffer,
    uint32_t seed
);


void add_splats_mcmc_tensor(
    int64_t cur_num_splats,
    int64_t num_add,
    float min_opacity,
    DeviceVector<float3> means, DeviceVector<float4> quats, DeviceVector<float3> scales, DeviceVector<float> opacs, DeviceVector<float3> features_dc, DeviceVector<float3> features_sh,
    DeviceVector<float3> g1_means, DeviceVector<float4> g1_quats, DeviceVector<float3> g1_scales, DeviceVector<float> g1_opacs, DeviceVector<float3> g1_features_dc, DeviceVector<float3> g1_features_sh,
    DeviceVector<float3> g2_means, DeviceVector<float4> g2_quats, DeviceVector<float3> g2_scales, DeviceVector<float> g2_opacs, DeviceVector<float3> g2_features_dc, DeviceVector<float3> g2_features_sh,
    DeviceVector<int32_t> bias_correction_steps,
    int sh_optim_bits,
    int num_sh,
    DeviceVector<float4> sh_quant_bounds,
    bool sh_bounds_per_splat,
    DeviceVector<uint8_t> sh_value_packed,
    DeviceVector<float2>  sh_value_bounds,
    int  sh_value_bits,
    bool sh_value_bounds_per_splat,
    int  num_sh_buffer,
    uint32_t seed
);


void mcmc_add_noise_tensor(
    int64_t num_splats,
    float scaler,
    DeviceVector<float3> means,
    DeviceVector<float3> log_scales,
    DeviceVector<float4> quats,
    DeviceVector<float> opacs
);


void revised_add_noise_tensor(
    int64_t num_splats,
    float scaler,
    DeviceVector<float> radii,
    DeviceVector<float3> means,
    DeviceVector<float3> log_scales,
    DeviceVector<float4> quats,
    DeviceVector<float> opacs
);


void long_axis_split_tensor(
    std::string primitive,
    DeviceVector<float3> log_scales,
    DeviceVector<float> logit_opacities,
    DeviceVector<float4> quats,
    DeviceVector<float3> new_log_scales,
    DeviceVector<float> new_logit_opacities,
    DeviceVector<float3> mean_deltas
);


void laplacian_edge_filter_tensor(
    DeviceTensor3D<float3> img_in,
    DeviceTensor3D<float> img_out
);


void smoothed_laplacian_edge_filter_tensor(
    DeviceTensor3D<float3> img_in,
    DeviceTensor3D<float> img_out
);


void canny_edge_filter_tensor(
    DeviceTensor3D<float3> img_in,
    bool* mask_in_ptr,
    DeviceTensor3D<float> img_out
);
