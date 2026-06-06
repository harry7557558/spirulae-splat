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

#include "DataManager.h"
#include "EngineConfig.h"

#include <map>
#include <string>
#include <vector>


// ============================================================
// Engine API -- not auto-generated; manually maintained
// ============================================================

// --- Lifecycle ---
//
// Reset the EngineState singleton to a freshly-constructed state and free all
// device memory owned by the global pool + scratch buffer. Must be called
// between training runs that swap datasets (e.g. ss_benchmark looping over
// scenes) -- without it the new run inherits the previous run's world splats,
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


// Warp-fused GT upload used by engine_train_step_managed when the current
// batch needs fisheye / equirectangular -> pinhole splitting. The GT RGB
// (and optionally mask) byte buffer is warped on GPU directly into a
// post-split float buffer; no full-res float intermediate is allocated.
// Throws when PPISP is enabled. Depth / normal are silently dropped on
// this path (the trainer must error before reaching here if they were
// requested).
void set_training_data_warped(
    std::string  input_model_name,    // "FISHEYE" / "EQUISOLID" / "EQUIRECTANGULAR"
    int B_in, int in_H, int in_W,
    int K, int out_H, int out_W,
    TorchTensorView gt_rgb,           // [B_in, in_H, in_W, 3] byte
    TorchTensorView gt_alpha,         // [B_in, mask_in_H, mask_in_W, 1] uint8 (nullable)
    int mask_in_H, int mask_in_W,
    TorchTensorView input_intrins,    // [B_in, 4] float
    TorchTensorView input_dist_coeffs,// [B_in, 10] float
    uint64_t axes_dev                 // device float ptr [K, 3, 3]
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

// quantize_optim: store optimizer state as uint8 + per-256-cell bounds
// (mirrors the SH-quantize path), trading a small numerical hit for ~3x
// smaller optimizer-state VRAM. Configurable per-type so RGB and geometry
// can be quantized independently.
//
// use_adagrad: if true, the bilagrid uses AdaGrad (lr_decay=0,
// weight_decay=0, initial_accumulator_value=0, eps=1e-15) instead of Adam.
// With quantize_optim=true, the per-cell accumulator is stored log-encoded
// via QuantizedTensorLog<8, 256> (1 byte/cell + 1 float2/256-cell block).
void engine_init_bilagrid_rgb(int n_grids, std::string type, int L, int H, int W,
                              bool quantize_optim, bool use_adagrad);
void engine_init_bilagrid_depth(int n_grids, int L, int H, int W,
                                bool quantize_optim, bool use_adagrad);
void engine_init_bilagrid_normal(int n_grids, int L, int H, int W,
                                 bool quantize_optim, bool use_adagrad);

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
// Two modes -- exactly one of these init calls activates background blending:
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


// Warp-fused fused train step. Used by engine_train_step_managed when a batch
// needs fisheye/equirectangular -> pinhole splitting. Camera table is set up
// at POST-split (B*K) PINHOLE; GT is byte-warped on GPU into a float
// post-split buffer with no intermediate full-res float image.
std::map<std::string, float> engine_train_step_warped(
    int step, int max_steps,
    std::string primitive, int sh_degree, bool packed,
    int out_W, int out_H,
    TorchTensorView post_viewmats,
    TorchTensorView post_intrins,
    TorchTensorView post_dist_coeffs,
    std::string input_camera_model,
    int B_in, int in_H, int in_W, int K,
    TorchTensorView input_intrins,
    TorchTensorView input_dist_coeffs,
    TorchTensorView gt_rgb_byte,
    TorchTensorView gt_alpha_byte,
    int mask_in_H, int mask_in_W,
    uint64_t axes_dev,
    TorchTensorView bilagrid_cam_indices,
    const EngineStepConfig& cfg
);


// --- DataManager-driven training step --------------------------------------
//
// Two-call setup + fused per-step pull. `engine_setup_data_manager` installs
// (or replaces) the engine's DataManager from a parsed dataset; subsequent
// `engine_train_step_managed(...)` calls pull the next training batch
// internally and run the same fused forward + loss/bwd + optim + densify
// pipeline as `engine_train_step`.
//
// This is the "engine handles data loading" entrypoint -- Python (or any
// future C++ trainer) no longer constructs per-step GT/camera tensors.

void engine_setup_data_manager(
    DataManagerConfig         cfg,
    // Per-camera model enum value (int; cast to CameraModelType internally).
    // Mixed pinhole + fisheye datasets just produce extra index groups.
    std::vector<int32_t>      camera_models,
    std::vector<std::string>  image_filenames,
    std::vector<std::string>  mask_filenames,
    std::vector<std::string>  depth_filenames,
    std::vector<std::string>  normal_filenames,
    std::vector<int32_t>      widths,
    std::vector<int32_t>      heights,
    // Per-input split factor + exclusive-prefix-sum offset into the
    // POST-split camera arrays below. Pass empty vectors to skip the warp
    // path entirely (K defaults to 1 per camera).
    std::vector<int32_t>      K_per_camera,
    std::vector<int32_t>      post_offsets,
    // POST-split camera params (length sum(K_per_camera) along the first dim).
    std::vector<float>        viewmats,
    std::vector<float>        intrins,
    std::vector<float>        dist_coeffs,
    // Per-INPUT intrins / dist_coeffs (length N). Required for the wide
    // warp kernel (fisheye / equisolid); pass empty when no fisheye-warp
    // camera is in the dataset.
    std::vector<float>        input_intrins,
    std::vector<float>        input_dist_coeffs,
    std::vector<int32_t>      train_indices,
    std::vector<int32_t>      val_indices);

std::map<std::string, float> engine_train_step_managed(
    int step, int max_steps,
    std::string primitive,
    int sh_degree,
    bool packed,
    const EngineStepConfig& cfg);

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
// Debug introspection: shape + D->H copy for engine().gt.rgb / gt.alpha and
// fwd.renders.rgb. Used by spirulae_splat.modules.debug_image to visualize
// what the engine actually sees per step, especially on the warp path.
std::tuple<int64_t, int64_t, int64_t, int64_t> engine_get_gt_rgb_shape();
std::tuple<int64_t, int64_t, int64_t, int64_t> engine_get_gt_alpha_shape();
std::tuple<int64_t, int64_t, int64_t, int64_t> engine_get_render_rgb_shape();
void engine_copy_gt_rgb_to_host(TorchTensorView out);   // float [B, H, W, 3]
void engine_copy_gt_alpha_to_host(TorchTensorView out); // bool  [B, H, W, 1]


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

// --- Viewer ---
//
// engine_viewer_init: lazy one-shot upload of the dataset-wide POST-split
// camera arrays (intrins / dist_coeffs / c2w / camera_models / W / H) and
// allocation of the device-side per-camera thumbnail cache + done-mask.
// After this, every set_training_data / set_training_data_warped also runs
// the thumbnail-fill kernel for any new post-camera ids in the batch (cheap
// kernel; gated on host counter once all cams have a thumbnail). Idempotent;
// safe to call multiple times (re-uploads + resets thumbnail state).
//
// camera_models / intrins / dist_coeffs / camera_to_worlds / widths /
// heights are at POST-split layout (length N_post == sum over input
// cameras of K), matching the DataManager's post-split arrays.
//
// camera_size: the visualization-frustum scale (knn-style), in world units.
void engine_viewer_init(
    TorchTensorView camera_models,   // [N_post]   int32
    TorchTensorView intrins,         // [N_post,4] float32
    TorchTensorView dist_coeffs,     // [N_post,10] float32
    TorchTensorView camera_to_worlds,// [N_post,3,4] float32 (y/z-flipped form)
    TorchTensorView widths,          // [N_post]   int32
    TorchTensorView heights,         // [N_post]   int32
    float camera_size);

// engine_blit_view: fused viewer-render annotation. Caches the BVH on the
// first show_training_cameras=true call; reads thumbnails out of the
// engine's device-side cache (no thumbnails arg). buffer_key selects which
// channel of the rendered output drives the colormap path:
//   "rgb"          render_rgbs is [H,W,3] float -- passed through.
//   "depth"        render_rgbs is [H,W,1] float -- turbo-colormapped.
//   "alpha"        render_rgbs is [H,W,1] float -- viridis-colormapped.
//   "normal"       [H,W,3] float in [0,1] -- passed through.
//   "depth_normal" [H,W,3] float in [0,1] -- passed through.
// depth + alpha are always [H,W,1] float and provide the occlusion test.
// view_* fields describe the visualization camera; out_rgb is [H,W,3] uint8
// (pre-allocated CUDA buffer).
void engine_blit_view(
    std::string     buffer_key,
    TorchTensorView render_buffer,   // [H,W,1|3] float (channel matches buffer_key)
    TorchTensorView render_depth,    // [H,W,1] float, may be null when buffer_key=="alpha"
    TorchTensorView render_alpha,    // [H,W,1] float
    int             view_camera_model,
    TorchTensorView view_intrins,    // [1,4]   float
    TorchTensorView view_viewmat,    // [4,4]   float
    TorchTensorView view_dist_coeffs,
    bool            show_training_cameras,
    TorchTensorView out_rgb          // [H,W,3] uint8, pre-allocated CUDA
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
