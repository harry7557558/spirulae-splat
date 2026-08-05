#pragma once

// The training config: every flag `ssplat train` accepts, in one table.
//
// This file is the single source of truth. Adding a row to
// SSPLAT_CONFIG_FIELDS makes the field appear in the CLI parser, in `--help`,
// in the GUI's "All Options" editor, in the run's config.json and in
// TrainerCore -- with no other edit. The struct is expanded from the same
// table, so a field cannot exist in one and not the other.
//
// Row: X(type, member, default, group, choices, help)
//
//   type     one of the scalar types below, or SsplatVec3i / SsplatVec3f.
//            std::array<T, N> cannot appear here: its comma would split the
//            macro argument.
//   member   the struct member AND, stringified, the CLI flag. The parser
//            treats '-' and '_' alike, so --sh-degree sets sh_degree.
//   default  a constant expression. Vector defaults go through ssplat_v3i() /
//            ssplat_v3f() for the same comma reason.
//   group    the section the flag is listed and nested under. Rows must stay
//            contiguous per group: --help and the GUI stream group headers as
//            they walk the table rather than sorting first.
//   choices  '|'-separated list for string fields; "" is free-form; "none"
//            means the empty string is allowed and displays as `none`.
//   help     one line, shown by --help and as the GUI tooltip. Keep the first
//            sentence self-contained -- --help truncates at the first ". ".

#include <array>
#include <limits>
#include <optional>
#include <string>
#include <string_view>

// Vector field types and their makers, so the table stays comma-free.
using SsplatVec3i = std::array<int, 3>;
using SsplatVec3f = std::array<float, 3>;
constexpr SsplatVec3i ssplat_v3i(int a, int b, int c) { return {a, b, c}; }
constexpr SsplatVec3f ssplat_v3f(float a, float b, float c) { return {a, b, c}; }

inline constexpr float kSsplatInf = std::numeric_limits<float>::infinity();

// The run's config.json keys the fields by flag name, with one exception:
// datamanager's split_batch would collide with model.split_batch, so its flag
// is dm_split_batch while the on-disk key stays split_batch. config.json is a
// read-back format (`ssplat mesh` and --resume parse it), so this mapping is
// compatibility, not policy -- a new field never needs an entry here.
constexpr const char* ssplat_json_key(const char* flag) {
    return std::string_view(flag) == "dm_split_batch" ? "split_batch" : flag;
}


// ===========================================================================
// The field table
// ===========================================================================

#define SSPLAT_CONFIG_FIELDS(X) \
                                                                              \
    /* ==== trainer -- run control: output, checkpoints, viewer ==== */       \
    X(std::string, data, {}, "trainer", "",                                   \
      "Path to dataset. Can be a Nerfstudio or a COLMAP dataset.")            \
    X(std::string, resume, "", "trainer", "none",                             \
      "Resume training from a checkpoint. Pass a run output dir (latest step-*.ckpt is used) or a specific step-*.ckpt dir. The run's config.json supplies the architecture/model/data config (run-control flags like num_iterations / save cadence / viewer are still taken from the CLI); requires the checkpoint to have been written with save_full_checkpoint.") \
    X(std::string, output_dir_prefix, "outputs", "trainer", "",               \
      "Prefix to output directory")                                           \
    X(std::string, output_dir_name, "", "trainer", "none",                    \
      "Output directory name relative to output_dir_prefix. If not specified, will set a generic combining current timestamp and dataset name.") \
    X(int, steps_per_save, 2000, "trainer", "",                               \
      "Save checkpoint every this number of steps. If -1, save only at the end. If zero, never save (used in benchmark).") \
    X(bool, save_only_latest_checkpoint, true, "trainer", "",                 \
      "Whether to save only last checkpoint")                                 \
    X(bool, save_full_checkpoint, false, "trainer", "",                       \
      "If True, each checkpoint's `state.tar` additionally includes the Resume slots -- world raw parameters (at max_num_splats) plus all optimizer state -- making the checkpoint sufficient to resume training, not just for inference. If False, only the Always slots (appearance/inference params) are saved alongside `splat.ply`.") \
    X(bool, save_eval_images, false, "trainer", "",                           \
      "Whether to save eval images at end of training")                       \
    X(int, num_iterations, 30000, "trainer", "",                              \
      "Number of training iterations")                                        \
    X(int, viewer_port, 7007, "trainer", "",                                  \
      "Port used by the web viewer")                                          \
    X(bool, disable_viewer, false, "trainer", "",                             \
      "If True, ss_trainer skips starting the viewer thread. Used by ss_benchmark so each scene runs without competing for the viewer port.") \
    X(bool, keep_viewer_alive, true, "trainer", "",                           \
      "If True, ss_trainer keeps the process (and thus the viewer) running after training + eval finish, so the result can still be inspected in the browser. Press Ctrl-C to exit. Ignored when disable_viewer=True.") \
                                                                              \
    /* ==== dataparser -- where the dataset is and how poses are normalized ==== */ \
    X(std::string, data_format, "", "dataparser", "colmap|nerfstudio|metashape|none", \
      "Data format, leave None to auto detect")                               \
    X(std::string, colmap_recon_dir, "", "dataparser", "none",                \
      "Path to COLMAP reconstruction relative to dataset directory (e.g. sparse/0). Will auto detect if not specified (picking the model with the most registered images when several exist).") \
    X(std::string, image_dir, "images", "dataparser", "",                     \
      "Path to images relative to dataset directory, used for COLMAP and Metashape datasets") \
    X(std::string, mask_dir, "masks", "dataparser", "",                       \
      "Path to masks relative to dataset directory, used for COLMAP and Metashape datasets") \
    X(std::string, depth_dir, "depths", "dataparser", "",                     \
      "Path to depth maps relative to dataset directory, used for COLMAP and Metashape datasets") \
    X(std::string, normal_dir, "normals", "dataparser", "",                   \
      "Path to normal maps relative to dataset directory, used for COLMAP and Metashape datasets") \
    X(std::string, metashape_xml, "", "dataparser", "none",                   \
      "Path to the Metashape xml file. Will automatically detect if not specified.") \
    X(std::string, metashape_ply, "", "dataparser", "none",                   \
      "Path to the Metashape point export ply file. Will automatically detect if not specified.") \
    X(std::string, metashape_psx, "", "dataparser", "none",                   \
      "Path to Metashape PSX file, used to resolve file name ambiguity when there are multiple images with the same file name") \
    X(float, rescale_camera_to_fit, 0.0f, "dataparser", "",                   \
      "Whether to check if image resolution match camera resolution and scale camera intrinsics accordingly if not. Set this to a number to divide intrinsics by that number, e.g. Mip-NeRF 360 and Zip-NeRF with images_(2|4) Set this to True to detect resolution, e.g. tankt_db [CLI: 0 = off, -1 = auto-detect from image resolution, > 0 = divide intrinsics by this]") \
    X(std::string, downscale_rounding_mode, "floor", "dataparser", "floor|ceil|round", \
      "Rounding mode applied to camera width/height when dividing by `rescale_camera_to_fit`. Use `round` to match the convention used by most image downscalers (e.g. Mip-NeRF 360 images_(2|4|8)).") \
    X(float, scene_scale, 1.0f, "dataparser", "",                             \
      "How much to scale the region of interest by.")                         \
    X(std::string, orientation_method, "up", "dataparser", "pca|up|vertical|none|gsplat", \
      "The method to use for orientation.")                                   \
    X(std::string, center_method, "poses", "dataparser", "poses|focus|none|gsplat", \
      "The method to use to center the poses.")                               \
    X(bool, auto_scale_poses, true, "dataparser", "",                         \
      "Whether to automatically scale the poses to fit in +/- 1 bounding box.") \
    X(float, outlier_threshold, kSsplatInf, "dataparser", "",                 \
      "Threshold to reject outlier camera poses.")                            \
    X(std::string, train_frame, "points", "dataparser", "normalized|camera|points", \
      "Coordinate frame in which splats are trained.")                        \
    X(std::string, eval_mode, "all", "dataparser", "fraction|filename|interval|all", \
      "The method to use for splitting the dataset into train and eval. Fraction splits based on a percentage for train and the remaining for eval. Filename splits based on filenames containing train/eval. Interval uses every nth frame for eval. All uses all the images for any split.") \
    X(float, train_split_fraction, 0.9f, "dataparser", "",                    \
      "The percentage of the dataset to use for training. Only used when eval_mode is train-split-fraction.") \
    X(int, eval_interval, 8, "dataparser", "",                                \
      "The interval between frames to use for eval. Only used when eval_mode is eval-interval.") \
    X(float, depth_unit_scale_factor, 0.001f, "dataparser", "",               \
      "Scales the depth values to meters. Default value is 0.001 for a millimeter to meter conversion.") \
    X(float, validation_fraction, 0.0f, "dataparser", "",                     \
      "Use this fraction of training images for validation. Stop training when performance on validation images start to drop.") \
                                                                              \
    /* ==== datamanager -- image caching, masks, warping ==== */              \
    X(int, max_batch_per_epoch, 800, "datamanager", "",                       \
      "Maximum number of batches per epoch, used for configuring batch size") \
    X(bool, dm_split_batch, false, "datamanager", "",                         \
      "Whether to one large batch into many small batches to avoid OOM, at cost of slower training") \
    X(std::string, cache_images, "disk", "datamanager", "cpu-pageable|cpu|gpu|disk", \
      "Image cache location. If \"cpu\", caches on cpu. If \"gpu\", caches on device. If \"cpu-pageable\", cache on cpu pageable memory (saves RAM but may cause error if spill to swap memory). If \"disk\", cache on disk (limited support).") \
    X(bool, load_depths, true, "datamanager", "",                             \
      "Whether to load depth maps, if exist")                                 \
    X(bool, load_normals, true, "datamanager", "",                            \
      "Whether to load normal maps, if exist")                                \
    X(float, mask_boundary_offset, 0.0f, "datamanager", "",                   \
      "Signed boundary offset applied to binarized masks at decode time, as a fraction of sqrt(W*H) of the decoded mask. Positive dilates (grows) foreground, negative erodes (shrinks). Runs on CPU during data loading via separable Felzenszwalb-Huttenlocher squared-Euclidean DT (exact, O(N) per row + col).") \
    X(bool, warp_to_pinhole, false, "datamanager", "",                        \
      "Whether to split a fisheye image into 5 undistorted pinhole images. Can sometimes give better quality and compatibility for dataset captured by fisheye/360 cameras.") \
    X(bool, warp_spherical_to_pinhole, true, "datamanager", "",               \
      "Whether to split an equirectangular (spherical panorama) image into 6 pinhole cubemap faces for training. When True (default), equirectangular images are split into 6 undistorted pinhole sub-images (the historical behavior). When False, train directly on the equirectangular image using an equirectangular projection/unprojection plugged into the linear/UT projection pipeline (supports 3dgs/mip and 3dgut primitives). Direct equirectangular training does not support depth/normal supervision.") \
    X(bool, deblur_training_images, false, "datamanager", "",                 \
      "Whether to use a custom trained deep learning model to deblur images before training") \
                                                                              \
    /* ==== model -- the splat model, losses, densification, regularizers ==== */ \
    X(std::string, primitive, "3dgs", "model", "3dgs|mip|3dgut",              \
      "Splat primitive to use")                                               \
    X(int, sh_degree, 3, "model", "",                                         \
      "Maximum degree of spherical harmonics to use.")                        \
    X(int, sh_degree_warmup_every, 1000, "model", "",                         \
      "Increase SH degree every this number of iterations")                   \
    X(std::string, background_mode, "black", "model", "black|noise|sh",       \
      "Background mode, black per convention, noise to discourage transparency, sh for skybox.") \
    X(int, background_noise_warmup, 2000, "model", "",                        \
      "Number of steps to warmup background noise. This applies when background_mode is noise") \
    X(float, background_noise_pre_warmup, 0.25f, "model", "",                 \
      "Weight of background noise at start of training (0 to 1). Higher value reduce the chance of washing away splat opacities near the beginning of training.") \
    X(int, background_sh_degree, 4, "model", "",                              \
      "SH degree for background color, only used when background_mode is sh.") \
    X(std::optional<float>, relative_scale, std::nullopt, "model", "",        \
      "Manually set scale when a scene is poorly scaled, i.e. increase this for large datasets. If not set, will use a scale agnostic optimizer. To prevent this, set it to 1.0.") \
    X(float, l1_weight, 1.0f, "model", "",                                    \
      "Weight of L1 loss, default 1.0")                                       \
    X(float, l2_weight, 0.0f, "model", "",                                    \
      "Weight of L2 loss, default 0.0")                                       \
    X(float, ssim_lambda, 0.2f, "model", "",                                  \
      "Weight of ssim loss; 0.2 for academic baseline, higher for potentially more high-frequency details, lower for less blurry background in outdoor scenes") \
    X(float, l1_weight_y, 0.0f, "model", "",                                  \
      "Weight of per-pixel BT.601 luma (Y) L1 loss.")                         \
    X(float, l2_weight_y, 0.0f, "model", "",                                  \
      "Weight of per-pixel BT.601 luma (Y) L2 loss.")                         \
    X(float, l2_weight_u, 0.0f, "model", "",                                  \
      "Weight of per-pixel BT.601 chroma U L2 loss.")                         \
    X(float, l2_weight_v, 0.0f, "model", "",                                  \
      "Weight of per-pixel BT.601 chroma V L2 loss.")                         \
    X(int, num_loss_scales, 0, "model", "",                                   \
      "Number of scales for image loss. For multi-scale loss, image is downscaled by 2 this number of times, and losses are averaged across scales. Improves convergence for high-resolution images.") \
    X(int, loss_scale_min_pixels, 1920, "model", "",                          \
      "If positive, overrides num_loss_scales per image based on resolution, in units of pixels. num_loss_scales is chosen so the smallest image dimension is halved down toward (but not below) this many pixels. e.g. with 2000: min dim 1999 -> num_loss_scales=0, 2000 -> 1, 4000 -> 2, 8000 -> 3, etc. Adapts per training step, so datasets with mixed image resolutions get the right count per image automatically.") \
    X(bool, use_camera_optimizer, false, "model", "",                         \
      "Whether to use camera optimizer Note: this only works well in patch batching mode") \
    X(bool, packed, true, "model", "",                                        \
      "Pack projection outputs, reduce VRAM usage at large batch size but can be slightly slower") \
    X(bool, use_bvh, false, "model", "",                                      \
      "Use BVH for splat-patch intersection test, may be faster when batching large number of small patches") \
    X(bool, use_fused_proj_bwd_optim, true, "model", "",                      \
      "Whether to use fused projection backward and optimizer. More memory efficient for large number of Gaussians, with slight performance hit.") \
    X(bool, split_batch, true, "model", "",                                   \
      "Split the camera batch into one-camera sub-batches inside the C++ train step. Per-splat grads accumulate via atomicAdd across sub-batches; a single optim+densify pass at the end consumes the accumulator with grad_scale = 1/B. Drops peak VRAM for the immediate projection / rasterization buffers by roughly 1/B. Per-image grad magnitude vs regularization weight stays batch-size invariant. Not compatible with use_fused_proj_bwd_optim or with the warped train-step path.") \
    X(int, quantization_level, 1, "model", "",                                \
      "SH quantization level: a single int that selects one of two (param bits, optim bits) configurations. 0 = off : 32-bit param, fp32 optim 1 = light : 16-bit param, 8-bit packed optim (2 B / cell) Collapsing the prior independent param+optim bit controls into a single level minimizes the FPBO kernel instantiations.") \
    X(std::string, optimizer_offload, "", "model", "sh|all|none",             \
      "Whether to offload optimizer momentum to CPU to save VRAM. This is only supported for Adam optimizer.") \
    X(int, resolution_schedule, 3000, "model", "",                            \
      "training starts at 1/d resolution, every n steps this is doubled")     \
    X(int, num_downscales, 0, "model", "",                                    \
      "At the beginning, resolution is 1/2^d, where d is this number")        \
    X(bool, use_mcmc, true, "model", "",                                      \
      "Must be True for 3DGS methods.")                                       \
    X(bool, preallocate_splat_tensors, true, "model", "",                     \
      "Whether to pre-allocate Gaussian attribute tensors to cap_max to avoid OOM during densification") \
    X(int, cap_max, 1000000, "model", "",                                     \
      "maximum number of splats, dataset-specific tuning required")           \
    X(float, min_init_fraction, 0.0f, "model", "",                            \
      "minimum fraction of splats out of cap_max at initialization")          \
    X(int, refine_every, 100, "model", "",                                    \
      "Densify every this number of steps")                                   \
    X(int, refine_start_iter, 500, "model", "",                               \
      "Start densification at this number of steps")                          \
    X(int, refine_stop_num_iter, 5000, "model", "",                           \
      "Stop densification at this number of steps before maximum number of training iterations") \
    X(int, refine_stop_iter, 25000, "model", "",                              \
      "Densification runs until max(this, num_iterations - refine_stop_num_iter). Without this floor, runs shorter than refine_stop_num_iter would never densify at all (num_iterations - refine_stop_num_iter goes negative), which confuses users.") \
    X(float, noise_lr, 80.0f, "model", "",                                    \
      "MCMC-like noise injection magnitude at start of training")             \
    X(float, noise_lr_final, 0.8f, "model", "",                               \
      "MCMC-like noise injection magnitude at end of training")               \
    X(float, min_opacity, 0.005f, "model", "",                                \
      "Minimum Gaussian opacity before relocation")                           \
    X(float, growth_factor, 1.05f, "model", "",                               \
      "Multiply number of splats by this number at each densification step")  \
    X(bool, use_revised_densification, true, "model", "",                     \
      "Whether to use revised densification instead of original MCMC.")       \
    X(std::string, densify_score_mode, "mean", "model", "mean|max|median|geom", \
      "How to accumulate per-splat scores across iterations for densification. \"mean\": running mean of |w|. \"max\": running max of |w|. \"median\": running median of |w| (approximation). \"geom\": running geometric mean of |w|.") \
    X(float, densify_score_blend_world_grad, 0.0f, "model", "",               \
      "Blend weight `w` in [0, 1] between the image-space loss score and the world-space gradient score for densification. The per-step score is (image-space accum_weight)^(1-w) * (||dL/dmean_world|| * max post-exp world scale)^w. The world-grad term favors world-space-large splats (e.g. distant background in unbounded outdoor scenes) that the image-space score under-weights; the geometric blend is invariant to each score's global scale so no cross-normalization is needed. 0 (default) = image-space score only, identical cost and behavior to before. 1 = world-grad score only; the per-pixel densification loss map (densify_loss_map_mode) is skipped entirely. In between, both scores are computed (one extra float per splat of VRAM).") \
    X(std::string, densify_loss_map_mode, "ssim_structure", "model", "none|loss_full|ssim_full|ssim_cs|ssim_structure|edge_aware|robust_edge_aware", \
      "What gets accumulated into the per-pixel densification loss map. The loss map is read by raster bwd to weight the per-splat accum_weight. Only active when use_revised_densification. Modes: \"none\": no loss map (uniform alpha*T accumulation). \"loss_full\": per-pixel L1/L2 + auxiliary supervisory terms + full SSIM (luminance*contrast*structure). \"ssim_full\": full SSIM only. \"ssim_cs\": contrast*structure SSIM (no luminance). \"ssim_structure\": structure-only SSIM, biases toward pattern/edge mismatches and ignores brightness/contrast errors. \"edge_aware\": canny edge magnitude of GT rgb (Plenoxels-style, https://arxiv.org/abs/2603.08661). Biases densification toward GT edges directly, regardless of how well the splats already reconstruct them. \"robust_edge_aware\": RobustNeRF-style Tukey biweight on the BT.601 luma of |render - GT|, capped at the per-image `densify_robust_edge_aware_quantile`, then canny. Near-zero where the render already matches GT, zeroed past the quantile cutoff so distractor pixels (people/cars/operator) don't pull splats toward them, and luminance-shift tolerant since a global DC residual has no spatial gradient. For num_loss_scales > 0 the map is computed per scale and the per-scale results are averaged (matches the multi-scale loss accumulation). Affects loss_map only; training gradients and scalar losses unchanged.") \
    X(float, densify_robust_edge_aware_quantile, 0.9f, "model", "",           \
      "Per-image quantile of the luma residual used as the Tukey biweight cutoff in `robust_edge_aware` mode. Pixels whose residual exceeds this quantile get zero densification weight. Lower values are more aggressive outlier rejection (good for distractor-heavy datasets); higher are more permissive (good for clean datasets where real edges may produce large residuals). Ignored unless `densify_loss_map_mode == \"robust_edge_aware\"`.") \
    X(bool, use_long_axis_split, true, "model", "",                           \
      "whether to use long-axis split described in https://arxiv.org/abs/2508.12313 for relocation and sample add. When combined with use_revised_densification, this can give less blurry background details for unbounded outdoor scenes.") \
    X(SsplatVec3f, long_axis_split_opacity_k, ssplat_v3f(0.6f, 0.6f, 4500.0f), "model", "", \
      "opacity split factor `k` for long-axis split, as (initial, final, warmup_steps). Each split child keeps opacity `logit^-1(k / (1 + exp(-logit_opacity) - k))`; `k` is linearly scheduled from `initial` to `final` over the first `warmup_steps` training steps, then held at `final`. Larger `k` preserves more opacity per child (denser, sharper); smaller `k` fades children faster.") \
    X(float, relocate_screen_size, kSsplatInf, "model", "",                   \
      "if a gaussian is more than this fraction of screen space, relocate it Useful for fisheye with 3DGUT, may drop PSNR for conventional cameras For likely better quality, use max_screen_size instead") \
    X(float, max_screen_size, 0.3f, "model", "",                              \
      "if a gaussian is more than this fraction of screen space, clip scale and increase opacity Intended to be an MCMC-friendly alternative of relocate_screen_size") \
    X(float, max_screen_size_clip_hardness, 1.5f, "model", "",                \
      "clip hardness for Gaussians with large screen space size, between 1 and infinity, larger is harder") \
    X(float, max_world_size, kSsplatInf, "model", "",                         \
      "if a gaussian is more than this of world space, clip scale Useful if you see huge floaters at a distance in large indoor space") \
    X(int, reset_alpha_every, 30, "model", "",                                \
      "Every this many refinement steps, reset the alpha. Only applies for opaque triangle splatting.") \
    X(bool, use_bilateral_grid, true, "model", "",                            \
      "If True, use bilateral grid to handle the ISP changes in the image space. This technique was introduced in the paper 'Bilateral Guided Radiance Field Processing' (https://bilarfpro.github.io/).") \
    X(SsplatVec3i, bilagrid_shape, ssplat_v3i(16, 16, 8), "model", "",        \
      "Shape of the bilateral grid, typically `16 16 8`, or `8 8 4` for scenes with low-texture surfaces.") \
    X(std::string, bilagrid_type, "ppisp", "model", "affine|ppisp|loglinear", \
      "What the bilateral grid predicts. affine: 4x3 matrix per original bilateral grid. ppisp: PPISP exposure and color parameters, generally gives less color shift but can be less numerically stable. loglinear: 3x3 linear transformation matrix with log-encoded diagonals, balances color shift and numerical stability.") \
    X(bool, use_bilateral_grid_for_geometry, true, "model", "",               \
      "If True, use bilateral grid for depth and normal (e.g. AI generated biased ones)") \
    X(SsplatVec3i, bilagrid_shape_geometry, ssplat_v3i(8, 8, 4), "model", "", \
      "Shape of the bilateral grid for depth and normal (X, Y, W)")           \
    X(bool, use_adagrad_bilagrid_optim, true, "model", "",                    \
      "Use AdaGrad (lr_decay=0, weight_decay=0, initial_accumulator_value=0, eps=1e-15) instead of Adam for all bilateral-grid parameters (RGB + depth + normal). When True, the bilagrid LR fields read from OptimizerConfig switch to ``bilagrid_adagrad_*_lr``. Bilagrid bit depths are coupled to `quantization_level`: level 0 = fp32; level 1 = 16-bit value + 8x2-bit optimizer state across all three bilagrid types.") \
    X(float, bilagrid_tv_loss_weight, 10.0f, "model", "",                     \
      "Total variation loss weight for bilateral grid used for radiance")     \
    X(float, color_shift_reg_weight, 0.0f, "model", "",                       \
      "Color-shift regularizer for the combined bilagrid + PPISP color transform. Penalizes the dataset-wide mean of sign(post - pre) per channel, where `pre` is the splat-side rendered RGB (input to whichever of bilagrid / PPISP runs first) and `post` is the final image fed to the photometric loss: R = w * ||EMA[mean_p sign(post - pre)]||^2. The gradient is injected on the POST-transforms v_render_rgb buffer and flows through each transform's vjp into its parameters (and as a small leak into the splats). Active when at least one of bilagrid_rgb / PPISP is enabled; 0 disables. Typical values: 0.01--1.0.") \
    X(int, color_shift_reg_ema_period, 750, "model", "",                      \
      "EMA decay length for the color-shift regularizer, in BATCHES. beta = max(0, 1 - 1/period). Should be roughly the number of batches per epoch so the EMA estimates the dataset-wide mean. Ignored when color_shift_reg_weight = 0.") \
    X(float, bilagrid_tv_loss_weight_geometry, 10.0f, "model", "",            \
      "Total variation loss weight for bilateral grid used for geometry")     \
    X(bool, use_ppisp, true, "model", "",                                     \
      "If True, use the PPISP model (https://research.nvidia.com/labs/sil/projects/ppisp/) to handle per-pixel color distortions.") \
    X(std::string, ppisp_param_type, "no_crf", "model", "original|rqs|no_crf", \
      "Parameterization for PPISP. \"original\" implements the original paper, \"rqs\" uses a parameterization that is more friendly to optimization and can produce better results in darker areas, \"no_crf\" omits the tone-curve (CRF) stage entirely and just clamps colors to [0,1] after exposure/vignetting/color correction.") \
    X(bool, use_adagrad_ppisp_optim, true, "model", "",                       \
      "Use unscheduled AdaGrad (lr_decay=0, weight_decay=0, initial_accumulator_value=0, eps=1e-15) instead of Adam for the PPISP parameter table. When True, the PPISP LR reads from ``OptimizerConfig.ppisp_adagrad_lr`` (constant) instead of the scheduled ``ppisp_lr``. No quantization path either way.") \
    X(bool, apply_ppisp_before_bilagrid, true, "model", "",                   \
      "When True, the PPISP forward runs BEFORE the RGB bilagrid (and PPISP backward runs AFTER bilagrid backward), i.e. render -> PPISP -> bilagrid -> loss. Otherwise: render -> bilagrid -> PPISP -> loss. Only meaningful when both ``use_bilateral_grid`` and ``use_ppisp`` are enabled.") \
    X(float, ppisp_reg_exposure_mean, 1.0f, "model", "",                      \
      "Encourage exposure mean ~ 0 to resolve SH <-> exposure ambiguity in PPISP.") \
    X(float, ppisp_reg_vig_center, 0.02f, "model", "",                        \
      "Encourage vignetting optical center near image center in PPISP.")      \
    X(float, ppisp_reg_vig_non_pos, 0.01f, "model", "",                       \
      "Penalize positive vignetting alpha coefficients in PPISP (should be <= 0).") \
    X(float, ppisp_reg_vig_channel_var, 0.1f, "model", "",                    \
      "Encourage similar vignetting across RGB channels in PPISP.")           \
    X(float, ppisp_reg_color_mean, 1.0f, "model", "",                         \
      "Encourage color correction mean ~ 0 across frames in PPISP.")          \
    X(float, ppisp_reg_crf_channel_var, 0.1f, "model", "",                    \
      "Encourage similar CRF parameters across RGB channels in PPISP.")       \
    X(bool, image_color_is_linear, false, "model", "",                        \
      "Whether to assume training images are in linear color space.")         \
    X(std::string, image_color_gamut, "", "model", "ACES2065-1|ACEScg|Rec.2020|AdobeRGB|DCI-P3|none", \
      "Color gamut of input images. If None, Rec.709 will be used. Note that tonemap is not applied.") \
    X(std::optional<bool>, splat_color_is_linear, std::nullopt, "model", "",  \
      "Whether to train splats in linear color space. If None, will use same as images.") \
    X(std::string, splat_color_gamut, "", "model", "Rec.709|ACES2065-1|ACEScg|Rec.2020|AdobeRGB|DCI-P3|none", \
      "Color gamut of trained splats. If None, will use same as images. Note that tonemap is not applied.") \
    X(std::optional<bool>, convert_initial_point_cloud_color, std::nullopt, "model", "", \
      "If True, this will assume color in initial point cloud is sRGB, and convert if images are in a linear or wide-gamut color space.") \
    X(std::optional<float>, scale_init, std::nullopt, "model", "",            \
      "Initial scale. If not set, auto decide")                               \
    X(std::optional<float>, opacity_init, std::nullopt, "model", "",          \
      "Initial opacity. If not set, auto decide")                             \
    X(bool, suppress_initial_scales, false, "model", "",                      \
      "Whether to suppress scales during initialization to discourage large floaters in vacant areas") \
    X(float, scale_regularization_weight, 0.0f, "model", "",                  \
      "If enabled, a scale regularization introduced in PhysGauss (https://xpandora.github.io/PhysGaussian/) is used for reducing huge spikey gaussians.") \
    X(float, max_gauss_ratio, 10.0f, "model", "",                             \
      "Threshold of ratio of gaussian max to min scale before applying regularization loss from the PhysGaussian paper") \
    X(float, depth_distortion_reg, 0.0f, "model", "",                         \
      "Weight for depth distortion regularizer")                              \
    X(float, normal_distortion_reg, 0.0f, "model", "",                        \
      "Weight for normal distortion regularizer")                             \
    X(float, rgb_distortion_reg, 0.0f, "model", "",                           \
      "Weight for rgb distortion regularizer")                                \
    X(int, distortion_reg_warmup, 6000, "model", "",                          \
      "warmup steps for depth regularizer, regularization weight ramps up")   \
    X(float, normal_reg_weight, 0.04f, "model", "",                           \
      "Weight for normal regularizer")                                        \
    X(int, normal_reg_warmup, 6000, "model", "",                              \
      "warmup steps for normal regularizer, regularization weight ramps up")  \
    X(float, alpha_reg_weight, 0.0f, "model", "",                             \
      "Weight for alpha regularizer (encourage alpha to go to either 0 or 1)") \
    X(int, alpha_reg_warmup, 12000, "model", "",                              \
      "warmup steps for alpha regularizer, regularization weight ramps up")   \
    X(int, reg_warmup_length, 0, "model", "",                                 \
      "Warmup steps for depth, normal, and alpha regularizers. only apply regularizers after this many steps.") \
    X(bool, apply_loss_for_mask, false, "model", "",                          \
      "Set this to False to use masks to ignore distractors (e.g. people and cars, area outside fisheye circle, over exposure) Set this to True to remove background (e.g. sky, background outside centered object)") \
    X(bool, enable_sky_masking, true, "model", "",                            \
      "If enabled, sky from depth map will be used for masking. Alpha loss will always be applied for sky.") \
    X(float, alpha_loss_weight, 0.01f, "model", "",                           \
      "Loss weight for alpha, applies when rendered alpha is above reference alpha") \
    X(float, alpha_loss_weight_under, 0.0f, "model", "",                      \
      "Loss weight for alpha, applies when rendered alpha is below reference alpha") \
    X(float, opacity_reg, 0.01f, "model", "",                                 \
      "Encourage low opacity to aid densification, per MCMC.")                \
    X(float, scale_reg, 0.01f, "model", "",                                   \
      "Encourage low scale, per MCMC.")                                       \
    X(float, opacity_decay, 0.0f, "model", "",                                \
      "Decay opacity to aid densification, per MRNF.")                        \
    X(float, scale_decay, 0.0f, "model", "",                                  \
      "Decay scale to aid densification, per MRNF.")                          \
    X(float, erank_reg, 0.0f, "model", "",                                    \
      "erank regularization weight, for 3DGS only - see *Effective Rank Analysis and Regularization for Enhanced 3D Gaussian Splatting, Hyung et al.*") \
    X(float, erank_reg_s3, 0.0f, "model", "",                                 \
      "erank regularization weight for smallest dimension, for 3DGS only")    \
    X(float, quat_norm_reg, 0.01f, "model", "",                               \
      "Weight to regularize quaternion norm to identity")                     \
    X(float, sh_reg, 0.001f, "model", "",                                     \
      "Regularize SH magnitude to find a balance between bilagrid/PPISP and improve generalizability.") \
    X(float, overexposure_reg, 0.0f, "model", "",                             \
      "Image-space L2 penalty on rendered RGB outside [0, 1]: L = w * mean(max(-x, x-1, 0)^2) over all pixels and channels of the raw rendered image (pre-bilagrid / pre-PPISP / pre-color-space). When non-zero a dedicated CUDA kernel adds dL/dx directly into the rendered-RGB gradient; the scalar loss value is never materialized.") \
    X(int, supervision_warmup, 0, "model", "",                                \
      "Start using foundation model depth at this number of steps")           \
    X(float, depth_supervision_weight, 0.0f, "model", "",                     \
      "Weight for depth supervision by comparing rendered depth with depth predicted by a foundation model Warn that this can reduce quality if AI generated depth is heavily biased") \
    X(bool, input_depth_is_ray_depth, false, "model", "",                     \
      "Whether the input/supervision depth maps store ray depth (Euclidean distance along the camera ray) rather than linear (z) depth. The rasterizer renders ray depth, so when this is False (the common case, e.g. most foundation-model depths) the GT depth is converted from linear to ray depth in place on the GPU before the depth bilateral grid / loss. Set True for depth maps already in ray depth, e.g. >180deg fisheye captures where linear depth is ill-defined.") \
    X(float, normal_supervision_weight, 0.01f, "model", "",                   \
      "Weight for normal supervision by comparing normal from rendered depth with normal from depth predicted by a foundation model") \
    X(float, mean_median_depth_weight, 0.0f, "model", "",                     \
      "L1 between the mean (expected) depth and the median depth, where both are nonzero.") \
    X(float, median_depth_normal_reg_weight, 0.0f, "model", "",               \
      "normal_loss between the normal from the median depth and the normal from the mean depth.") \
    X(float, median_normal_supervision_weight, 0.0f, "model", "",             \
      "normal_loss between the normal from the median depth and the reference (foundation-model) normal.") \
    X(float, median_render_normal_reg_weight, 0.0f, "model", "",              \
      "normal_loss between the normal from the median depth and the rendered normal (placeholder until render_normal exists).") \
    X(int, median_warmup, 6000, "model", "",                                  \
      "Linear warmup length (steps) shared by all four median-depth loss weights: each ramps 0->full over this many steps.") \
    X(std::string, overfit_score_aggregation_mode, "min", "model", "max|min|mean", \
      "Mode to aggregate multiple overfitting objectives. Use max for more aggressive early stopping, min for more conservative early stopping, and mean for something in between.") \
    X(int, validation_loss_average_window, 500, "model", "",                  \
      "Window to calculate moving average validation loss for early stop")    \
    X(int, early_stop_patience, 1000, "model", "",                            \
      "Stop training if overfitting score remains positive for this number of iterations") \
    X(int, early_stop_warmup, 12000, "model", "",                             \
      "Warmup steps for early stop, will not early stop before this number of steps Recommend setting this number no less than regularization warmups") \
                                                                              \
    /* ==== optimizer -- learning rates and schedules ==== */                 \
    X(std::optional<int>, max_steps, std::nullopt, "optimizer", "",           \
      "")                                                                     \
    X(bool, use_scale_agnostic_mean, true, "optimizer", "",                   \
      "")                                                                     \
    X(bool, use_per_splat_bias_correction, true, "optimizer", "",             \
      "")                                                                     \
    X(float, means_lr, 0.000128f, "optimizer", "",                            \
      "")                                                                     \
    X(std::optional<float>, means_lr_final, 1.6e-06f, "optimizer", "",        \
      "")                                                                     \
    X(float, scales_lr, 0.02f, "optimizer", "",                               \
      "")                                                                     \
    X(std::optional<float>, scales_lr_final, 0.005f, "optimizer", "",         \
      "")                                                                     \
    X(float, quats_lr, 0.0015f, "optimizer", "",                              \
      "")                                                                     \
    X(float, opacities_lr, 0.025f, "optimizer", "",                           \
      "")                                                                     \
    X(float, features_dc_lr, 0.005f, "optimizer", "",                         \
      "")                                                                     \
    X(float, features_sh_lr, 0.00025f, "optimizer", "",                       \
      "")                                                                     \
    X(float, background_dc_lr, 0.0025f, "optimizer", "",                      \
      "")                                                                     \
    X(float, background_sh_lr, 0.0005f, "optimizer", "",                      \
      "")                                                                     \
    X(float, bilagrid_lr, 0.002f, "optimizer", "",                            \
      "")                                                                     \
    X(std::optional<float>, bilagrid_lr_final, 0.0001f, "optimizer", "",      \
      "")                                                                     \
    X(int, bilagrid_lr_warmup, 1000, "optimizer", "",                         \
      "")                                                                     \
    X(float, bilagrid_depth_lr, 0.002f, "optimizer", "",                      \
      "")                                                                     \
    X(std::optional<float>, bilagrid_depth_lr_final, 0.0001f, "optimizer", "", \
      "")                                                                     \
    X(int, bilagrid_depth_lr_warmup, 2000, "optimizer", "",                   \
      "")                                                                     \
    X(float, bilagrid_normal_lr, 0.0005f, "optimizer", "",                    \
      "")                                                                     \
    X(std::optional<float>, bilagrid_normal_lr_final, 4e-05f, "optimizer", "", \
      "")                                                                     \
    X(int, bilagrid_normal_lr_warmup, 2000, "optimizer", "",                  \
      "")                                                                     \
    X(float, bilagrid_adagrad_lr, 0.04f, "optimizer", "",                     \
      "")                                                                     \
    X(float, bilagrid_adagrad_depth_lr, 0.04f, "optimizer", "",               \
      "")                                                                     \
    X(float, bilagrid_adagrad_normal_lr, 0.01f, "optimizer", "",              \
      "")                                                                     \
    X(float, ppisp_lr, 0.002f, "optimizer", "",                               \
      "")                                                                     \
    X(std::optional<float>, ppisp_lr_final, 2e-05f, "optimizer", "",          \
      "")                                                                     \
    X(int, ppisp_lr_warmup, 500, "optimizer", "",                             \
      "")                                                                     \
    X(float, ppisp_adagrad_lr, 0.1f, "optimizer", "",                         \
      "")                                                                     \
    X(float, camera_opt_lr, 0.0001f, "optimizer", "",                         \
      "")                                                                     \
    X(std::optional<float>, camera_opt_lr_final, 5e-07f, "optimizer", "",     \
      "")                                                                     \
    X(int, camera_opt_lr_warmup, 1000, "optimizer", "",                       \
      "")                                                                     \
    /* end */


// ===========================================================================
// The config struct, expanded from the table above
// ===========================================================================

struct SsplatConfig {
#define SSPLAT_DECLARE_FIELD(type, member, default_, group, choices, help) \
    type member = default_;
    SSPLAT_CONFIG_FIELDS(SSPLAT_DECLARE_FIELD)
#undef SSPLAT_DECLARE_FIELD
};

// Fields whose default ({}) is not a usable value. Checked after flag
// parsing, so --help still works without them.
#define SSPLAT_CONFIG_REQUIRED_FIELDS(X) \
    X(data) \
    /* end */


// ===========================================================================
// Presets -- named bundles of default overrides, selected as
// `ssplat train <preset>`. "3dgs" is the base config and applies nothing.
// ===========================================================================

struct SsplatPresetInfo { const char* name; const char* help; };
inline constexpr SsplatPresetInfo kSsplatPresets[] = {
    {"3dgs", "Generic method that works well for most datasets."},
    {"360-camera", "Preset for training on original distorted images captured by 360 cameras (e.g. Insta360, DJI Osmo). Recommended if your dataset contains fisheye images with a circle visible."},
    {"in-the-wild", "Preset for datasets consisting of internet images, with extreme lighting variation, with un-masked outliers, and/or shot with long focal lengths."},
    {"linear-color", "Preset for training splats in linear color spaces (e.g. ACEScg)."},
    {"synthetic", "Preset for training splats on synthetic datasets rendered with constant exposure."},
    {"meshing", "Preset for training splats for meshing. Use `spirulae-meshing` to convert trained splats to mesh."},
    {"academic-baseline", "Preset that replicates 3DGS MCMC as faithful as possible."},
};

// Returns false for an unknown preset name.
inline bool ssplat_apply_preset(SsplatConfig& c, const std::string& name) {
    if (name == "3dgs") {
        return true;
    }
    if (name == "360-camera") {
        c.warp_to_pinhole = true;
        c.mask_boundary_offset = -0.025f;
        c.primitive = "mip";
        c.long_axis_split_opacity_k = {0.5f, 0.6f, 15000.0f};
        c.input_depth_is_ray_depth = true;
        return true;
    }
    if (name == "in-the-wild") {
        c.center_method = "focus";
        c.outlier_threshold = 10.0f;
        c.load_depths = true;
        c.load_normals = true;
        c.mask_boundary_offset = -0.025f;
        c.densify_score_mode = "median";
        c.densify_loss_map_mode = "robust_edge_aware";
        c.densify_robust_edge_aware_quantile = 0.75f;
        c.ssim_lambda = 0.1f;
        c.rgb_distortion_reg = 0.1f;
        c.depth_distortion_reg = 0.01f;
        c.sh_degree_warmup_every = 0;
        c.long_axis_split_opacity_k = {0.5f, 0.6f, 30000.0f};
        c.noise_lr = 10.0f;
        c.noise_lr_final = 0.1f;
        c.erank_reg = 0.1f;
        c.means_lr = 5e-05f;
        c.means_lr_final = 1e-07f;
        return true;
    }
    if (name == "linear-color") {
        c.splat_color_gamut = "ACEScg";
        c.splat_color_is_linear = true;
        c.image_color_gamut = "Rec.2020";
        c.image_color_is_linear = false;
        c.background_mode = "noise";
        return true;
    }
    if (name == "synthetic") {
        c.min_init_fraction = 0.1f;
        c.use_bilateral_grid = false;
        c.use_ppisp = false;
        c.use_bilateral_grid_for_geometry = false;
        c.long_axis_split_opacity_k = {0.5f, 0.6f, 25000.0f};
        return true;
    }
    if (name == "meshing") {
        c.primitive = "3dgut";
        c.sh_degree = 0;
        c.sh_reg = 10.0f;
        c.overexposure_reg = 10.0f;
        c.background_mode = "noise";
        c.depth_distortion_reg = 0.01f;
        c.normal_distortion_reg = 0.01f;
        c.mean_median_depth_weight = 0.01f;
        c.median_depth_normal_reg_weight = 0.01f;
        c.normal_supervision_weight = 0.01f;
        c.median_normal_supervision_weight = 0.01f;
        c.median_render_normal_reg_weight = 0.01f;
        c.erank_reg = 0.01f;
        c.erank_reg_s3 = 0.01f;
        return true;
    }
    if (name == "academic-baseline") {
        c.eval_mode = "interval";
        c.eval_interval = 8;
        c.center_method = "gsplat";
        c.orientation_method = "gsplat";
        c.max_batch_per_epoch = 387420489;
        c.load_depths = false;
        c.load_normals = false;
        c.primitive = "3dgs";
        c.relative_scale = 1.0f;
        c.use_bilateral_grid = false;
        c.use_bilateral_grid_for_geometry = false;
        c.use_ppisp = false;
        c.use_revised_densification = false;
        c.densify_loss_map_mode = "none";
        c.use_long_axis_split = false;
        c.use_fused_proj_bwd_optim = false;
        c.quantization_level = 0;
        c.max_screen_size = kSsplatInf;
        c.max_world_size = kSsplatInf;
        c.suppress_initial_scales = false;
        c.scale_init = 0.1f;
        c.opacity_init = 0.5f;
        c.depth_distortion_reg = 0.0f;
        c.normal_reg_weight = 0.0f;
        c.alpha_reg_weight = 0.0f;
        c.alpha_loss_weight = 0.0f;
        c.alpha_loss_weight_under = 0.0f;
        c.erank_reg = 0.0f;
        c.erank_reg_s3 = 0.0f;
        c.quat_norm_reg = 0.0f;
        c.sh_reg = 0.0f;
        c.normal_supervision_weight = 0.0f;
        c.opacity_reg = 0.01f;
        c.scale_reg = 0.01f;
        c.max_steps = 30000;
        c.use_scale_agnostic_mean = false;
        c.use_per_splat_bias_correction = false;
        c.means_lr = 0.00016f;
        c.means_lr_final = 1.6e-06f;
        c.scales_lr = 0.005f;
        c.scales_lr_final = std::nullopt;
        c.quats_lr = 0.001f;
        c.opacities_lr = 0.05f;
        c.features_dc_lr = 0.0025f;
        c.features_sh_lr = 0.000125f;
        return true;
    }
    (void)c;
    return false;
}
