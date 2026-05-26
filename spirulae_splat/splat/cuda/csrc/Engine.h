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
    bool packed,
    TorchTensorView out_rgb,
    TorchTensorView out_depth,
    TorchTensorView out_Ts
);

// --- Loss + backward (combined) ---

std::map<std::string, float> engine_compute_loss_backward(
    int step,
    // loss weights
    std::array<float, (int)LossWeightIndex::length> loss_weights,
    float w_ssim,
    int num_loss_scales,
    bool compute_loss_map,
    // pre-allocated gradient output tensors (zero-initialized by Python)
    TorchTensorView v_means,
    TorchTensorView v_quats,
    TorchTensorView v_scales,
    TorchTensorView v_opacities,
    TorchTensorView v_features_dc,
    TorchTensorView v_features_sh
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
    // splat data (params + grads + optimizer state) as TorchTensorViews
    TorchTensorView means,    TorchTensorView v_means,    TorchTensorView g1_means,    TorchTensorView g2_means,
    TorchTensorView quats,    TorchTensorView v_quats,    TorchTensorView g1_quats,    TorchTensorView g2_quats,
    TorchTensorView scales,   TorchTensorView v_scales,   TorchTensorView g1_scales,   TorchTensorView g2_scales,
    TorchTensorView opacities,TorchTensorView v_opacities,TorchTensorView g1_opacities,TorchTensorView g2_opacities,
    TorchTensorView features_dc, TorchTensorView v_features_dc, TorchTensorView g1_features_dc, TorchTensorView g2_features_dc,
    TorchTensorView features_sh, TorchTensorView v_features_sh, TorchTensorView g1_features_sh, TorchTensorView g2_features_sh,
    TorchTensorView radii,
    TorchTensorView quant_bounds_sh  // null if !quantize_sh
);

// --- Backward only (existing, kept for compatibility) ---

void backward_3dgs(
    TorchTensorView v_rgb,
    TorchTensorView v_depth,
    TorchTensorView v_Ts,
    TorchTensorView v_means,
    TorchTensorView v_quats,
    TorchTensorView v_scales,
    TorchTensorView v_opacities,
    TorchTensorView v_features_dc,
    TorchTensorView v_features_sh
);
