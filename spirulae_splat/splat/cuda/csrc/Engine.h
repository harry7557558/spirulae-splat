#pragma once

#include "Tensor.h"

#include "IntersectTile.cuh"
#include "BackgroundSphericalHarmonics.cuh"
#include "PerSplatLoss.cuh"
#include "PerPixelLoss.cuh"
#include "PixelWise.cuh"
#include "FusedSSIM.cuh"
#include "SplatTileIntersector.cuh"
// #include "Projection.cuh"
#include "ProjectionFwd.cuh"
#include "ProjectionBwd.cuh"
#include "ProjectionPackedFwd.cuh"
#include "RasterizationFwd.cuh"
#include "RasterizationBwd.cuh"
#include "RasterizationEval3DFwd.cuh"
#include "RasterizationEval3DBwd.cuh"
#include "Optimizer.cuh"
#include "FusedProjectionBwdOptim.cuh"
#include "Densify.cuh"
#include "BilagridUtils.cuh"
#include "Visualizer.cuh"

#include "EngineConfig.h"

#include <map>
#include <string>
#include <vector>


// ============================================================
// Engine API — not auto-generated; manually maintained
// ============================================================

// --- Lifecycle ---
//
// Reset the EngineState singleton to a freshly-constructed state and free all
// device memory owned by the global pool + scratch buffer. Must be called
// between training runs that swap datasets (e.g. ss_benchmark looping over
// scenes) — without it the new run inherits the previous run's world splats,
// camera table, bilagrid / PPISP / background config, optimizer moments, and
// color-space matrices, which produces broken renders and wrong metrics.
void engine_reset();

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
    // RGB: float32 [0,1], uint8 [0,255], or uint16 [0,65535]. uint8/uint16
    // uploaded as raw bytes; converted to float on the GPU.
    TorchTensorView gt_rgb,          // [C, H, W, 3]
    // Depth: float32 (passthrough) or uint16 (cast on GPU, no scaling). Null OK.
    TorchTensorView gt_depth,        // [C, H, W, 1]
    // Normal: float32 (passthrough) or uint8 (x/127.5 - 1 on GPU). Null OK.
    TorchTensorView gt_normal,       // [C, H, W, 3]
    // External mask (per-pixel bool/uint8). Sets Buffers::has_mask. Null OK.
    // Drives the RGB mask + alpha-sup target inside per_pixel_losses.slang.
    TorchTensorView gt_alpha         // [C, H, W, 1]
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
    bool compute_loss_map,
    bool structure_only_loss_map,
    // Image-space overexposure regularization weight (see LossConfig).
    // Zero (default) disables the kernel launch entirely.
    float overexposure_reg_weight = 0.0f
);

// --- Optimizer step ---

void engine_optim_step(int step, const OptimConfig& cfg);

// --- Bilagrid (RGB / depth / normal) ---
// Allocates and identity-initializes a bilagrid per camera. Must be called
// after set_camera_params. RGB applies to the rendered prediction; depth and
// normal apply to GT (matching the Python flow in training_losses.py).

// quantize_optim: store Adam g1 / g2 as uint8 + per-256-cell float4 bounds
// (mirrors the SH-quantize path), trading a small numerical hit for ~3x
// smaller optimizer-state VRAM. Configurable per-type so RGB and geometry
// can be quantized independently.
void engine_init_bilagrid_rgb(int n_grids, std::string type, int L, int H, int W,
                              bool quantize_optim);
void engine_init_bilagrid_depth(int n_grids, int L, int H, int W,
                                bool quantize_optim);
void engine_init_bilagrid_normal(int n_grids, int L, int H, int W,
                                 bool quantize_optim);

// Apply forward bilagrid for the current batch. cam_indices is a [C_batch]
// int32 tensor of per-image camera-table indices; pass null/empty for identity
// (image i uses grid slot i). Must be called between forward_3dgs and
// engine_compute_loss_backward.
void engine_bilagrid_forward(TorchTensorView cam_indices);

// Adam step + optional TV-loss regularization for each enabled bilagrid type.
void engine_bilagrid_optim_step(int step, const BilagridStepConfig& cfg);

// --- PPISP (RGB only, applied AFTER bilagrid). ---
// Allocates a per-camera PPISP parameter table and seeds it with the type's
// default values ("original" -> 12 zeros, 12 zeros, 3x(a,a,b,0); "rqs" -> all 0).
void engine_init_ppisp(int n_grids, std::string param_type);

// Apply PPISP forward in place on the current rendered RGB; saves a pre-PPISP
// copy used by backward. cam_indices: [C_batch] int32, or null/empty for
// identity. Must be called after forward_3dgs (and after engine_bilagrid_forward
// when bilagrid is also enabled).
void engine_ppisp_forward(TorchTensorView cam_indices);

// Adam step over the PPISP parameter table. Also folds in the 6 PPISP
// regularization-loss gradients (using the provided per-loss weights) before
// the Adam update. Pass loss weights in PPISPRegLossIndex order.
void engine_ppisp_optim_step(int step, const PpispStepConfig& cfg);

// --- Background blending (applied BEFORE bilagrid/PPISP) ---
//
// Two modes — exactly one of these init calls activates background blending:
//   Noise: random per-pixel color, warmup-weighted via cfg.background.randomize_weight.
//   SH:    skybox = SH(world ray dir), DC (slot 0) seeds the "background_color"
//          plus higher-order bands; both updated by Adam when train_color=true.
//
// dc_color is the linear-space DC color used at SH init time (set slot 0).
void engine_init_background_noise(bool splat_color_is_linear);
void engine_init_background_sh(
    int sh_degree, bool splat_color_is_linear);

// Per-iter (seed, randomize_weight) for the next forward_3dgs background blend.
// Training calls this each step; the viewer/eval path can ignore it and reuse
// the last-stashed values (the engine defaults of 0/0 produce a uniform
// half-gray noise blend, avoiding per-frame flicker in noise-mode viewer renders).
void engine_set_background_step_params(uint32_t seed, float randomize_weight);

// Adam step over the SH coefficient table. No-op for Noise mode or when SH
// training was not enabled at init.
void engine_background_optim_step(int step, const BackgroundStepConfig& cfg);

// Copy the engine's background image for the current camera setup to a host
// (..., H, W, 3) float buffer. SH mode: returns the skybox rendered by the
// most recent forward_3dgs. Noise mode: returns a uniform mean-color image.
// Returns 1 on success, 0 if no engine background is active.
int engine_copy_background_to_host(TorchTensorView out_image);

// --- Linear / wide-gamut color space ---
//
// Configure linear / wide-gamut handling for both the splat (predicted) RGB
// and the GT image color space. When enabled, the splat-side conversion
// (linear or wide-gamut -> sRGB) runs at the end of every forward, and is
// inverted through the vjp on the loss-side gradient before raster bwd. The
// image-side conversion (same direction) runs once per upload in
// set_training_data. Pass row-major 3x3 source->Rec.709 matrices (matches
// get_color_transform_matrix in Python). Pass empty matrices when the
// corresponding side is not enabled.
void engine_init_color_space(
    bool splat_enabled,
    bool splat_is_linear,
    std::vector<float> splat_color_matrix,    // [9], row-major
    bool image_enabled,
    bool image_is_linear,
    std::vector<float> image_color_matrix     // [9], row-major
);

// --- Densification step ---

// Returns: number of splats added (0 if no densification this step)
int engine_densify_step(int step, int max_steps, const DensifyConfig& cfg);

// --- Fused training step (set_camera + set_training + forward + loss/bwd + optim + densify) ---
//
// Bilagrid + PPISP entries in `cfg` are no-ops when the corresponding
// engine_init_* was never called. `bilagrid_cam_indices` is a [C_batch] int32
// tensor for per-image grid-slot lookup (null/empty -> identity).
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
    TorchTensorView bilagrid_cam_indices,
    const EngineStepConfig& cfg
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

// `out_rgb_raw` is the pre-color-space-conversion render (linear / wide-gamut)
// stashed by the color-space forward hook. Null OK; when the engine has no
// color space configured, `out_rgb` already holds the same values so the
// caller can pass null. When the engine has a color space, `out_rgb` is the
// sRGB post-conversion render and `out_rgb_raw` is the working-space version.
void engine_copy_render_to_host(
    TorchTensorView out_rgb,      // [C, H, W, 3] float32, CPU
    TorchTensorView out_depth,    // [C, H, W, 1] float32, CPU
    TorchTensorView out_Ts,       // [C, H, W, 1] float32, CPU
    TorchTensorView out_rgb_raw   // [C, H, W, 3] float32, CPU, optional
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

// --- Checkpoint save ---
//
// Default mode writes into `output_dir`:
//   splat.ply              : cur_num_splats only, NaN+low-opacity-filtered
//   bilagrid_rgb.npy       : [N_cam, C, L, H, W] float (when bilagrid RGB enabled)
//   bilagrid_depth.npy     : [N_cam, 2, L, H, W] float (when enabled)
//   bilagrid_depth_scalars.npy : [N_cam] float
//   bilagrid_normal.npy    : [N_cam, 3, L, H, W] float (when enabled)
//   ppisp.npy              : [N_cam, P] float (when PPISP enabled)
//   meta.txt               : key=value training-state metadata
//
// `full_dump=true` additionally writes a `full/` subfolder containing every
// world-splat parameter (sized to max_num_splats), every bilagrid grid and
// PPISP table, and all Adam optimizer states (g1, g2; quantized variants when
// active). One .npy per buffer, for offline inspection.
void engine_save_checkpoint(
    std::string output_dir,
    bool full_dump,
    int step
);
