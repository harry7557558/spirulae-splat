#pragma once

#include "Tensor.h"

#include "IntersectTile.cuh"
#include "BackgroundSphericalHarmonics.cuh"
#include "PerSplatLoss.cuh"
#include "PerPixelLoss.cuh"
#include "PixelWise.cuh"
#include "FusedSSIM.cuh"
#include "SplatTileIntersector.cuh"
#include "SVHash.cuh"
// #include "Projection.cuh"
#include "ProjectionFwd.cuh"
#include "ProjectionBwd.cuh"
#include "ProjectionPackedFwd.cuh"
#include "RasterizationFwd.cuh"
#include "RasterizationBwd.cuh"
#include "RasterizationEval3DFwd.cuh"
#include "RasterizationEval3DBwd.cuh"
#include "RasterizationSortedEval3DFwd.cuh"
#include "RasterizationSortedEval3DBwd.cuh"
#include "Optimizer.cuh"
#include "FusedProjectionBwdOptim.cuh"
#include "Densify.cuh"
#include "BilagridUtils.cuh"
#include "Visualizer.cuh"

#include <map>
#include <string>


// ============================================================
// Engine API — not auto-generated; manually maintained
// ============================================================

// --- Setup ---

void set_data_3dgs(
    int64_t num_splats,
    TorchTensorView means,
    TorchTensorView quats,
    TorchTensorView scales,
    TorchTensorView opacities,
    TorchTensorView features_dc,
    TorchTensorView features_sh
);

void set_camera_params(
    int width,
    int height,
    std::string camera_model,
    TorchTensorView viewmats,
    TorchTensorView intrins,
    TorchTensorView dist_coeffs
);

void set_training_data(
    TorchTensorView gt_rgb,          // [C, H, W, 3] float32 (sRGB 0-1)
    TorchTensorView gt_depth,        // [C, H, W, 1] float32, or null
    TorchTensorView gt_normal,       // [C, H, W, 3] float32, or null
    TorchTensorView gt_alpha,        // [C, H, W, 1] bool,    or null
    TorchTensorView gt_rgb_mask,     // [C, H, W, 1] bool,    or null
    TorchTensorView gt_depth_mask,   // [C, H, W, 1] bool,    or null
    TorchTensorView gt_normal_mask,  // [C, H, W, 1] bool,    or null
    TorchTensorView gt_alpha_mask    // [C, H, W, 1] bool,    or null
);

// --- Forward ---

void forward_3dgs(
    std::string primitive,
    int sh_degree,
    bool packed
);

// --- Loss + backward (combined) ---

std::map<std::string, float> engine_compute_loss_backward(
    int step,
    // loss weights
    std::array<float, (int)LossWeightIndex::length> loss_weights,
    float w_ssim,
    int num_loss_scales,
    bool compute_loss_map
);

// --- Optimizer step ---

void engine_optim_step(
    int step,
    // learning rates
    float lr_means,
    float lr_quats,
    float lr_scales,
    float lr_opacities,
    float lr_features_dc,
    float lr_features_sh,
    // regularization
    float max_gauss_ratio,
    float scale_regularization_weight,
    float mcmc_opacity_reg_weight,
    float mcmc_scale_reg_weight,
    float erank_reg_weight,
    float erank_reg_weight_s3,
    float quat_norm_reg_weight,
    float sh_reg_weight,
    bool use_scale_agnostic_mean,
    bool quantize_sh,
    bool use_per_splat_bias_correction
);

// --- Densification step ---

// Returns: number of splats added (0 if no densification this step)
int engine_densify_step(
    int step,
    int max_steps,
    // densification config
    int refine_start_iter,
    int refine_stop_num_iter,
    int refine_every,
    float growth_factor,
    float min_opacity,
    float max_screen_size,
    float max_screen_size_clip_hardness,
    float max_world_size,
    float noise_lr,
    float noise_lr_final,
    float relocate_heuristic_weight
);

// --- Fused training step (set_camera + set_training + forward + loss/bwd + optim + densify) ---

std::map<std::string, float> engine_train_step(
    int step, int max_steps,
    std::string primitive,
    int sh_degree,
    bool packed,
    int width, int height, std::string camera_model,
    TorchTensorView viewmats,
    TorchTensorView intrins,
    TorchTensorView dist_coeffs,
    TorchTensorView gt_rgb,
    TorchTensorView gt_depth,
    TorchTensorView gt_normal,
    TorchTensorView gt_alpha,
    TorchTensorView gt_rgb_mask,
    TorchTensorView gt_depth_mask,
    TorchTensorView gt_normal_mask,
    TorchTensorView gt_alpha_mask,
    std::array<float, (int)LossWeightIndex::length> loss_weights,
    float w_ssim,
    int num_loss_scales,
    bool compute_loss_map,
    float lr_means, float lr_quats, float lr_scales, float lr_opacities,
    float lr_features_dc, float lr_features_sh,
    float max_gauss_ratio, float scale_regularization_weight,
    float mcmc_opacity_reg_weight, float mcmc_scale_reg_weight,
    float erank_reg_weight, float erank_reg_weight_s3,
    float quat_norm_reg_weight, float sh_reg_weight,
    bool use_scale_agnostic_mean,
    bool quantize_sh,
    bool use_per_splat_bias_correction,
    int refine_start_iter, int refine_stop_num_iter, int refine_every,
    float growth_factor, float min_opacity,
    float max_screen_size, float max_screen_size_clip_hardness,
    float max_world_size,
    float noise_lr, float noise_lr_final,
    float relocate_heuristic_weight
);

// --- Debug rendering ---

void engine_debug_forward(
    TorchTensorView override_features_dc,  // [max_N, 3] custom DC color, or null to use original
    int override_sh_degree,                 // -1 to use original
    TorchTensorView out_rgb                 // [C, H, W, 3] output
);

// --- Query internal state ---

void engine_copy_accum_buffer(TorchTensorView dst);
int64_t engine_get_cur_num_splats();

void engine_copy_render_to_host(
    TorchTensorView out_rgb,    // [C, H, W, 3] float32, CPU
    TorchTensorView out_depth,  // [C, H, W, 1] float32, CPU
    TorchTensorView out_Ts      // [C, H, W, 1] float32, CPU
);

void engine_copy_splats_to_host(
    TorchTensorView means,
    TorchTensorView quats,
    TorchTensorView scales,
    TorchTensorView opacities,
    TorchTensorView features_dc,
    TorchTensorView features_sh
);

// Pool VRAM breakdown: [(key, used_bytes, cap_bytes), ...]
std::vector<std::tuple<std::string, size_t, size_t>> engine_get_pool_breakdown();
size_t engine_get_scratch_bytes();
