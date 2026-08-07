#pragma once

// The training config: every flag `spirula train` accepts, in one table.
//
// This file is the single source of truth. Adding a row to
// SS_CONFIG_FIELDS makes the field appear in the CLI parser, in `--help`,
// in the GUI's "All Options" editor, in the run's config.json and in
// TrainerCore -- with no other edit. The struct is expanded from the same
// table, so a field cannot exist in one and not the other.
//
// Row: X(type, member, default, group, choices, help)
//
//   type     one of the scalar types below, or TrainVec3i / TrainVec3f.
//            std::array<T, N> cannot appear here: its comma would split the
//            macro argument.
//   member   the struct member AND, stringified, the CLI flag. The parser
//            treats '-' and '_' alike, so --sh-degree sets sh_degree.
//   default  a constant expression. Vector defaults go through train_v3i() /
//            train_v3f() for the same comma reason.
//   group    the section the flag is listed and nested under. Rows must stay
//            contiguous per group: --help and the GUI stream group headers as
//            they walk the table rather than sorting first.
//   choices  '|'-separated list for string fields; "" is free-form; "none"
//            means the empty string is allowed and displays as `none`.
//   help     one line, shown by --help and as the GUI tooltip. Written for
//            someone training a scene, not for someone reading the kernels:
//            say what the flag does to the result and which way to move it,
//            not how it is implemented. No paper citations, no formulas, no
//            "True"/"False" (the GUI draws a checkbox). Keep the first
//            sentence self-contained and under ~105 characters -- --help
//            prints only up to the first ". " and elides past 110. That also
//            means no "e.g." in the first sentence: it ends it early.

#include <array>
#include <limits>
#include <optional>
#include <string>
#include <string_view>

// Vector field types and their makers, so the table stays comma-free.
using TrainVec3i = std::array<int, 3>;
using TrainVec3f = std::array<float, 3>;
constexpr TrainVec3i train_v3i(int a, int b, int c) { return {a, b, c}; }
constexpr TrainVec3f train_v3f(float a, float b, float c) { return {a, b, c}; }

inline constexpr float kTrainInf = std::numeric_limits<float>::infinity();

// The run's config.json keys the fields by flag name, with one exception:
// datamanager's split_batch would collide with model.split_batch, so its flag
// is dm_split_batch while the on-disk key stays split_batch. config.json is a
// read-back format (`spirula mesh` and --resume parse it), so this mapping is
// compatibility, not policy -- a new field never needs an entry here.
constexpr const char* train_json_key(const char* flag) {
    return std::string_view(flag) == "dm_split_batch" ? "split_batch" : flag;
}


// ===========================================================================
// The field table
// ===========================================================================

#define SS_CONFIG_FIELDS(X) \
                                                                              \
    /* ==== trainer -- run control: output, checkpoints, viewer ==== */       \
    X(std::string, data, {}, "trainer", "",                                   \
      "Folder holding the dataset to train on. COLMAP, Nerfstudio and Metashape layouts are all recognized.") \
    X(std::string, resume, "", "trainer", "none",                             \
      "Continue a previous run instead of starting from scratch. Point this at a run's output folder to pick up its newest checkpoint, or at one specific step-*.ckpt folder. The model and dataset settings come from that run, while run length, save cadence and viewer settings still come from the command line. Only works if the earlier run was saved with save_full_checkpoint turned on.") \
    X(std::string, output_dir_prefix, "outputs", "trainer", "",               \
      "Folder that each run's output directory is created inside.")           \
    X(std::string, output_dir_name, "", "trainer", "none",                    \
      "Name of this run's folder inside the output prefix. Leave empty to get a name built from the dataset name and the current time.") \
    X(int, steps_per_save, 2000, "trainer", "",                               \
      "How often a checkpoint is written, in steps. Use -1 to save only when training finishes, or 0 to never save. Frequent saves cost disk space and a little time.") \
    X(bool, save_only_latest_checkpoint, true, "trainer", "",                 \
      "Keep only the newest checkpoint and delete older ones as training goes. Turn off to keep the whole history, which uses considerably more disk space.") \
    X(bool, save_full_checkpoint, false, "trainer", "",                       \
      "Also store everything needed to resume training later, not just the finished splats. Checkpoints get much larger because they carry every splat slot and the optimizer state. Leave off if you only want the exported splat file.") \
    X(bool, save_eval_images, false, "trainer", "",                           \
      "Write rendered and reference images for the held-out views when training finishes. Useful for judging quality, at the cost of a little disk space.") \
    X(int, num_iterations, 30000, "trainer", "",                              \
      "How long to train, in steps. More steps often give better visual results with diminishing returns, and a proportionally longer wait.") \
    X(int, viewer_port, 7007, "trainer", "",                                  \
      "Network port the built-in web viewer listens on. Change it if that port is already taken.") \
    X(bool, disable_viewer, false, "trainer", "",                             \
      "Do not start the web viewer. Frees a little memory and avoids port conflicts when several trainings run at once.") \
    X(bool, keep_viewer_alive, true, "trainer", "",                           \
      "Keep the program running after training finishes so the result stays open in the viewer. Press Ctrl-C to exit. Has no effect when the viewer is disabled.") \
                                                                              \
    /* ==== dataparser -- where the dataset is and how poses are normalized ==== */ \
    X(std::string, data_format, "", "dataparser", "colmap|nerfstudio|metashape|none", \
      "Which dataset layout to read. Leave empty to detect it from the folder contents.") \
    X(std::string, colmap_recon_dir, "", "dataparser", "none",                \
      "Which COLMAP reconstruction to read, relative to the dataset folder, such as sparse/0. Leave empty to pick automatically: the reconstruction with the most registered images wins.") \
    X(std::string, image_dir, "images", "dataparser", "",                     \
      "Subfolder holding the training images, for COLMAP and Metashape datasets.") \
    X(std::string, mask_dir, "masks", "dataparser", "",                       \
      "Subfolder holding the image masks, for COLMAP and Metashape datasets. What a mask means is set by apply_loss_for_mask.") \
    X(std::string, depth_dir, "depths", "dataparser", "",                     \
      "Subfolder holding depth maps, for COLMAP and Metashape datasets. Only read when load_depths is on.") \
    X(std::string, normal_dir, "normals", "dataparser", "",                   \
      "Subfolder holding normal maps, for COLMAP and Metashape datasets. Only read when load_normals is on.") \
    X(std::string, metashape_xml, "", "dataparser", "none",                   \
      "Metashape camera export to read. Leave empty to find it automatically inside the dataset folder.") \
    X(std::string, metashape_ply, "", "dataparser", "none",                   \
      "Metashape point cloud export used to seed the splats. Leave empty to find it automatically.") \
    X(std::string, metashape_psx, "", "dataparser", "none",                   \
      "Metashape project file, used to resolve ambiguity when several images in the project share a file name.") \
    X(float, rescale_camera_to_fit, 0.0f, "dataparser", "",                   \
      "Fix a mismatch between image size and the camera parameters stored in the dataset. Set it to the factor the images were shrunk by, such as 2 when training on images_2, or 0 to leave the cameras alone. Auto-detection (-1) is not supported yet.") \
    X(std::string, downscale_rounding_mode, "floor", "dataparser", "floor|ceil|round", \
      "How image size is rounded when divided by rescale_camera_to_fit. Most image downscalers round, so switch to `round` if a pre-shrunk dataset comes out a pixel off and the render looks slightly shifted.") \
    X(std::string, orientation_method, "up", "dataparser", "pca|up|vertical|none|gsplat", \
      "How the scene is rotated to stand upright. This only affects how the result is framed for viewing, not the splats themselves. Anything other than `up` is approximated for now.") \
    X(std::string, center_method, "poses", "dataparser", "poses|focus|none|gsplat", \
      "How the scene's origin is chosen. Like the orientation setting, this affects framing rather than the splats themselves. Anything other than `poses` is approximated for now.") \
    X(bool, auto_scale_poses, true, "dataparser", "",                         \
      "Normalize the scene so the cameras fit in a unit-sized box. This keeps learning rates and regularizers meaningful across scenes of very different physical size; turn it off only if the dataset is already scaled the way you want.") \
    X(float, outlier_threshold, kTrainInf, "dataparser", "",                 \
      "Discard cameras that sit far outside the rest of the capture. Lowering this rejects more, which helps when a few badly estimated poses stretch the scene and throw its scale off. Leave at infinity to keep every camera.") \
    X(std::string, train_frame, "points", "dataparser", "normalized|camera|points", \
      "Coordinate frame the splats are trained in. Only `points`, the dataset's own frame, is supported.") \
    X(std::string, eval_mode, "all", "dataparser", "fraction|filename|interval|all", \
      "How images are split between training and held-out evaluation. `all` trains on every image and reports numbers on views it has already seen. `interval` holds out every Nth image, `fraction` holds out a share of them, and `filename` looks for train or eval in the file names. Holding images out gives honest quality numbers at the cost of a few training views.") \
    X(float, train_split_fraction, 0.9f, "dataparser", "",                    \
      "Share of the images used for training when eval_mode is `fraction`; the rest are held out.") \
    X(int, eval_interval, 8, "dataparser", "",                                \
      "Hold out every Nth image when eval_mode is `interval`. Larger values keep more images for training.") \
    X(float, depth_unit_scale_factor, 0.001f, "dataparser", "",               \
      "Multiplier that turns the stored depth values into scene units. The default reads them as millimeters.") \
    X(float, validation_fraction, 0.0f, "dataparser", "",                     \
      "Share of the training images set aside to watch for overfitting. Training can then stop before quality starts to drop. Early stopping is not active yet: the images are held out but the run always goes to the end.") \
                                                                              \
    /* ==== datamanager -- image caching, masks, warping ==== */              \
    X(int, max_batch_per_epoch, 800, "datamanager", "",                       \
      "Target number of steps per pass over the dataset, which decides how many images each step uses. Raising it makes each step lighter and cheaper; lowering it groups more images into a step, which is steadier but slower and needs more memory. Datasets smaller than this simply use one image per step.") \
    X(std::string, cache_images, "disk", "datamanager", "cpu|gpu|disk", \
      "Where decoded training images are kept between steps. `disk` re-reads them and uses the least memory, `cpu` keeps them in RAM for faster steps. `gpu` is not supported yet.") \
    X(bool, load_depths, true, "datamanager", "",                             \
      "Use the dataset's depth maps when they exist. They drive depth supervision; turn off to ignore them.") \
    X(bool, load_normals, true, "datamanager", "",                            \
      "Use the dataset's normal maps when they exist. They drive normal supervision; turn off to ignore them.") \
    X(float, mask_boundary_offset, 0.0f, "datamanager", "",                   \
      "Grow or shrink masks by this fraction of the image size. Negative values pull the mask edge inward, trimming halos and bad pixels along the boundary; positive values push it outward.") \
    X(bool, warp_to_pinhole, false, "datamanager", "",                        \
      "Split each fisheye image into five ordinary perspective views before training. Often gives better quality and wider compatibility for fisheye and 360 captures, at the cost of more images to process.") \
    X(bool, warp_spherical_to_pinhole, true, "datamanager", "",               \
      "Split each 360 panorama into six cube faces before training. Turn this off to train directly on the panorama, which keeps the original pixels but cannot use depth or normal supervision.") \
    X(bool, deblur_training_images, false, "datamanager", "",                 \
      "Sharpen blurry photos with a learned deblurring model before training. Not supported yet.") \
                                                                              \
    /* ==== model -- the splat model, losses, densification, regularizers ==== */ \
    X(std::string, primitive, "3dgs", "model", "3dgs|mip|3dgut",              \
      "Shape used for each splat. `3dgs` is the standard Gaussian most compatible with mainstream viewers, `mip` reduces aliasing when views differ a lot in distance or resolution, and `3dgut` handles wide-angle and distorted lenses more accurately and gives cleaner geometry for meshing.") \
    X(int, sh_degree, 3, "model", "",                                         \
      "How much the color of a splat may change with viewing angle. Higher shows sharper reflections and shading as the camera moves, at the cost of memory and file size; 0 gives flat, view-independent color. Values of 4 or higher have limited support in mainstream viewers.") \
    X(int, sh_degree_warmup_every, 1000, "model", "",                         \
      "How many steps between each step up in view-dependent color detail. Introducing it gradually keeps early training stable and stops reflections from being baked in too early. Small values reach full detail almost immediately.") \
    X(std::string, background_mode, "black", "model", "black|noise|sh",       \
      "What fills pixels no splat covers. `black` is the usual choice, `noise` discourages a semi-transparent haze from forming in empty space, and `sh` learns a skybox so distant background is represented instead of ignored.") \
    X(int, background_noise_warmup, 2000, "model", "",                        \
      "How many steps the background noise takes to reach full strength. Only used with the `noise` background.") \
    X(float, background_noise_pre_warmup, 0.25f, "model", "",                 \
      "How strong the background noise is at the very start, from 0 to 1. Higher values keep splats from being washed away in the first steps.") \
    X(int, background_sh_degree, 4, "model", "",                              \
      "How detailed the learned skybox may be. Higher captures finer sky and distant scenery. Only used with the `sh` background.") \
    X(std::optional<float>, relative_scale, std::nullopt, "model", "",        \
      "Multiply the whole scene by this factor before training. Raise it when a large capture comes out too small for detail to resolve. Leave unset to let the optimizer cope with scene scale on its own.") \
    X(float, l1_weight, 1.0f, "model", "",                                    \
      "Weight of plain per-pixel color error. This is the main term driving color accuracy.") \
    X(float, l2_weight, 0.0f, "model", "",                                    \
      "Weight of squared per-pixel color error. It punishes large mistakes harder than l1_weight, which makes color settle faster but also chases outliers such as moving objects.") \
    X(float, ssim_lambda, 0.2f, "model", "",                                  \
      "How much the loss cares about local structure instead of exact pixel color. Higher brings out fine texture and high-frequency detail; lower gives a smoother, less noisy background, which sometimes looks better in outdoor scenes.") \
    X(float, l1_weight_y, 0.0f, "model", "",                                  \
      "Extra weight on brightness error alone, ignoring hue. Raising it favors luminance detail over color accuracy.") \
    X(float, l2_weight_y, 0.0f, "model", "",                                  \
      "Extra weight on squared brightness error. Same idea as l1_weight_y, but large brightness mistakes count for much more.") \
    X(float, l2_weight_u, 0.0f, "model", "",                                  \
      "Extra weight on the blue-versus-yellow color error. Raising it tightens hue accuracy in that direction at the expense of detail elsewhere.") \
    X(float, l2_weight_v, 0.0f, "model", "",                                  \
      "Extra weight on the red-versus-cyan color error. Raising it tightens hue accuracy in that direction at the expense of detail elsewhere.") \
    X(int, num_loss_scales, 0, "model", "",                                   \
      "How many progressively smaller copies of each image the loss also compares. Looking at several sizes helps large smooth areas converge on high-resolution datasets instead of only fine detail. Normally left for loss_scale_min_pixels to decide.") \
    X(int, loss_scale_min_pixels, 1920, "model", "",                          \
      "Pick the number of loss scales automatically from image size. Images are halved until the shorter side is near this many pixels, so a dataset that mixes resolutions gets the right amount for each image. Set to 0 to use num_loss_scales instead.") \
    X(bool, use_camera_optimizer, false, "model", "",                         \
      "Let training nudge the camera poses to absorb small pose errors. Not supported yet.") \
    X(bool, packed, true, "model", "",                                        \
      "Store projection results compactly. Cuts GPU memory when many images are processed per step, sometimes at a small speed cost.") \
    X(bool, use_bvh, false, "model", "",                                      \
      "Use a spatial index for splat-tile intersection, which can help when batching many small patches. Not supported yet.") \
    X(bool, use_fused_proj_bwd_optim, true, "model", "",                      \
      "Merge the backward pass and the parameter update into one operation. Uses noticeably less memory at large splat counts, for a small speed cost.") \
    X(bool, split_batch, true, "model", "",                                   \
      "Process a step's images one at a time inside the training step. Cuts peak GPU memory roughly in proportion to how many images each step uses, and gives the same result as processing them together.") \
    X(int, quantization_level, 1, "model", "",                                \
      "How compactly splat colors are stored during training. 1 roughly halves the memory spent on view-dependent color with little visible difference; 0 keeps full precision.") \
    X(std::string, optimizer_offload, "", "model", "sh|all|none",             \
      "Move optimizer state to system memory to free up GPU memory. Not supported yet.") \
    X(bool, preallocate_splat_tensors, true, "model", "",                     \
      "Reserve memory for the maximum splat count up front. Avoids running out of GPU memory partway through as splats are added, at the cost of holding that memory from the start.") \
    X(int, cap_max, 1000000, "model", "",                                     \
      "Largest number of splats the scene may grow to. This is the main quality dial: raising it captures more detail and produces a bigger file that renders more slowly. Worth tuning per scene.") \
    X(float, min_init_fraction, 0.0f, "model", "",                            \
      "Smallest starting splat count, as a share of cap_max. Raise it when the initial point cloud is sparse, such as synthetic scenes, so training has enough to work with.") \
    X(int, refine_every, 100, "model", "",                                    \
      "How many steps between rounds of adding and relocating splats. Smaller reacts to missing detail sooner; larger is calmer and slightly cheaper.") \
    X(int, refine_start_iter, 500, "model", "",                               \
      "Step at which splats first start being added. Waiting a little lets the initial splats settle before the count starts growing.") \
    X(int, refine_stop_num_iter, 5000, "model", "",                           \
      "Stop adding splats this many steps before the end. The remaining steps polish what already exists instead of introducing new splats that never get refined.") \
    X(int, refine_stop_iter, 25000, "model", "",                              \
      "Earliest step at which splat growth may stop. Growth ends at whichever comes later, this step or refine_stop_num_iter before the end, so short runs still get to add splats at all.") \
    X(float, noise_lr, 80.0f, "model", "",                                    \
      "How much random jitter is applied to splat positions early in training. Jitter helps splats escape bad spots and spread into unfilled areas; too much of it blurs detail.") \
    X(float, noise_lr_final, 0.8f, "model", "",                               \
      "How much position jitter is left at the end of training. Lower lets detail settle and sharpen over the final steps.") \
    X(float, min_opacity, 0.005f, "model", "",                                \
      "Splats fainter than this get recycled into places that need them. Raising it prunes harder and keeps the splat budget on visible surfaces.") \
    X(float, growth_factor, 1.05f, "model", "",                               \
      "How fast the splat count grows at each round, as a multiplier. Higher reaches cap_max sooner; lower grows gradually, which tends to place splats more carefully.") \
    X(bool, use_revised_densification, true, "model", "",                     \
      "Use the improved rule for deciding where new splats go. It usually recovers missing detail faster; turn off to match the original method.") \
    X(std::string, densify_score_mode, "mean", "model", "mean|max|median|geom", \
      "How a splat's need-more-detail score builds up over time. `mean` is the balanced default, `max` reacts to a single bad view, `median` ignores the occasional odd view and helps when people or cars move through the scene, and `geom` sits between mean and median.") \
    X(float, densify_score_blend_world_grad, 0.0f, "model", "",               \
      "Balance between adding splats where the image looks wrong and where splats are physically large. Raise toward 1 to spend more splats on big distant structures that image-based scoring tends to starve; 0 uses image error alone.") \
    X(std::string, densify_loss_map_mode, "ssim_structure", "model", "none|loss_full|ssim_full|ssim_cs|ssim_structure|edge_aware|robust_edge_aware", \
      "What kind of error decides where new splats are added. `ssim_structure` targets mismatched patterns and edges while ignoring brightness differences. `ssim_full`, `ssim_cs` and `loss_full` fold in progressively more of the raw color error. `edge_aware` chases edges in the reference photos whether or not they are already reconstructed well. `robust_edge_aware` does the same but ignores the worst-matching pixels, so moving people and cars do not attract splats. `none` spreads new splats evenly.") \
    X(float, densify_robust_edge_aware_quantile, 0.9f, "model", "",           \
      "How much of the worst-matching image area is ignored when placing splats in `robust_edge_aware` mode. Lower ignores more, which suits captures full of moving distractions; higher keeps more, which suits clean captures where large errors are real detail.") \
    X(bool, use_long_axis_split, true, "model", "",                           \
      "Split stretched splats along their long axis when adding detail, rather than splitting them evenly. Gives less blurry distant background in large outdoor scenes.") \
    X(TrainVec3f, long_axis_split_opacity_k, train_v3f(0.5f, 0.6f, 8000.0f), "model", "", \
      "How much opacity each half keeps when a splat is split. Given as a starting value, a final value, and how many steps to move between them. Higher keeps the halves denser and sharper; lower encourages floaters to fade and relocate to where details are needed.") \
    X(float, max_screen_size, 0.3f, "model", "",                              \
      "Shrink splats that cover more than this share of the screen instead of letting them stay huge. Keeps big blobby splats from smearing across the image.") \
    X(float, max_screen_size_clip_hardness, 1.5f, "model", "",                \
      "How firmly the screen-size limit is enforced, from 1 upward. Higher clamps oversized splats decisively; lower eases them down.") \
    X(float, max_world_size, kTrainInf, "model", "",                         \
      "Shrink splats bigger than this in world units. Set it when huge floaters show up in the distance in large indoor spaces.") \
    X(bool, use_bilateral_grid, true, "model", "",                            \
      "Give each photo its own smooth color correction, absorbing exposure and white balance drift between shots. The splats then keep one consistent color instead of averaging every camera's quirks. Turn off for synthetic or already-consistent datasets.") \
    X(TrainVec3i, bilagrid_shape, train_v3i(16, 16, 8), "model", "",        \
      "How finely the per-photo color correction may vary, as width, height and brightness steps. Finer grids fix more localized shifts and use more VRAM; coarser is safer on flat, low-texture surfaces where a fine grid starts eating real detail.") \
    X(std::string, bilagrid_type, "ppisp", "model", "affine|ppisp|loglinear", \
      "What the per-photo color correction is allowed to do. `ppisp` adjusts exposure and color gain and shifts hue the least. `affine` is a full color matrix, the most flexible but the most prone to color drift. `loglinear` sits in between.") \
    X(bool, use_bilateral_grid_for_geometry, true, "model", "",               \
      "Apply the same per-photo correction to depth and normal maps. Biased AI-generated maps can then still be used without dragging the geometry off.") \
    X(TrainVec3i, bilagrid_shape_geometry, train_v3i(8, 8, 4), "model", "", \
      "How finely the depth and normal correction may vary. Same meaning as bilagrid_shape, for geometry rather than color.") \
    X(bool, use_adagrad_bilagrid_optim, true, "model", "",                    \
      "Use a steadier update rule for the per-photo color correction. Generally more stable and needs less tuning; turn off to use the scheduled learning rates instead.") \
    X(float, bilagrid_tv_loss_weight, 10.0f, "model", "",                     \
      "How smooth the per-photo color correction has to be. Higher keeps corrections gentle and global; lower lets them vary from place to place, which can start absorbing real image detail.") \
    X(float, color_shift_reg_weight, 0.0f, "model", "",                       \
      "Keep the per-photo corrections from tinting the result overall. Raise it if the finished splats come out consistently warmer, cooler, darker or brighter than the photos. 0 turns it off, and 0.01 to 1 is the useful range.") \
    X(int, color_shift_reg_ema_period, 750, "model", "",                      \
      "How many steps the color-shift check averages over. It should be roughly one pass over the dataset so the average reflects every photo. Ignored when color_shift_reg_weight is 0.") \
    X(float, bilagrid_tv_loss_weight_geometry, 10.0f, "model", "",            \
      "How smooth the depth and normal correction has to be. Higher keeps it gentle and global; lower lets it vary from place to place.") \
    X(bool, use_ppisp, true, "model", "",                                     \
      "Model per-pixel camera effects such as vignetting, exposure and lens color response. Keeps darkened corners and per-photo exposure shifts out of the splats themselves.") \
    X(std::string, ppisp_param_type, "no_crf", "model", "original|rqs|no_crf", \
      "Which camera effects get modeled. `no_crf` covers exposure, vignetting and color, then simply clips the result. `original` adds a tone curve on top. `rqs` uses a tone curve that behaves better in dark areas.") \
    X(bool, use_adagrad_ppisp_optim, true, "model", "",                       \
      "Use a steadier update rule for the camera-effect model. Generally more stable, needs less tuning, and leads to fewer floaters; turn off to use the scheduled learning rate instead.") \
    X(bool, apply_ppisp_before_bilagrid, true, "model", "",                   \
      "Run the camera-effect model before the per-photo color correction rather than after. Only matters when both are on, and decides which of the two absorbs a given color difference.") \
    X(float, ppisp_reg_exposure_mean, 1.0f, "model", "",                      \
      "Keep estimated exposures centered around neutral. Stops overall brightness from being counted twice between the splats and the camera model.") \
    X(float, ppisp_reg_vig_center, 0.02f, "model", "",                        \
      "Keep estimated vignetting centered near the middle of the image rather than drifting toward a corner.") \
    X(float, ppisp_reg_vig_non_pos, 0.01f, "model", "",                       \
      "Keep vignetting darkening the corners rather than brightening them, which is what real lenses do.") \
    X(float, ppisp_reg_vig_channel_var, 0.1f, "model", "",                    \
      "Keep vignetting similar across red, green and blue, so image corners do not pick up a color cast.") \
    X(float, ppisp_reg_color_mean, 1.0f, "model", "",                         \
      "Keep the per-photo color corrections centered, so no overall tint gets baked into the splats.") \
    X(float, ppisp_reg_crf_channel_var, 0.1f, "model", "",                    \
      "Keep the tone curve similar across red, green and blue, so brightness changes do not shift hue.") \
    X(bool, image_color_is_linear, false, "model", "",                        \
      "Treat the input images as linear light rather than ordinary display-encoded photos. Set this for renders or captures exported in linear.") \
    X(std::string, image_color_gamut, "", "model", "ACES2065-1|ACEScg|Rec.2020|AdobeRGB|DCI-P3|none", \
      "Color space the input images were captured in. Leave empty for ordinary sRGB or Rec.709 photos. No tone mapping is applied.") \
    X(std::optional<bool>, splat_color_is_linear, std::nullopt, "model", "",  \
      "Train splat colors in linear light. Leave unset to follow the input images. Linear color holds bright highlights better for HDR work.") \
    X(std::string, splat_color_gamut, "", "model", "Rec.709|ACES2065-1|ACEScg|Rec.2020|AdobeRGB|DCI-P3|none", \
      "Color space the trained splats are stored in. Leave unset to follow the input images. A wider gamut preserves saturated colors for later grading but needs a viewer that understands it. No tone mapping is applied.") \
    X(std::optional<bool>, convert_initial_point_cloud_color, std::nullopt, "model", "", \
      "Read the seed point cloud's colors as ordinary sRGB and convert them into the training color space. Turn on when starting colors look wrong in a linear or wide-gamut run.") \
    X(std::optional<float>, scale_init, std::nullopt, "model", "",            \
      "How big each splat starts out. Leave unset to derive it from the point cloud; larger fills space faster but starts blurrier.") \
    X(std::optional<float>, opacity_init, std::nullopt, "model", "",          \
      "How solid each splat starts out. Leave unset to choose automatically; lower makes early training more forgiving, higher locks geometry in sooner.") \
    X(bool, suppress_initial_scales, false, "model", "",                      \
      "Start splats small where the point cloud is sparse. Keeps them from blooming into large floaters over empty space.") \
    X(float, scale_regularization_weight, 0.0f, "model", "",                  \
      "Penalize splats that are far longer in one direction than another, which suppresses long spiky artifacts.") \
    X(float, max_gauss_ratio, 10.0f, "model", "",                             \
      "How stretched a splat may get before the spiky-splat penalty applies. Lower forces rounder splats.") \
    X(float, depth_distortion_reg, 0.0f, "model", "",                         \
      "Encourage each pixel's depth to come from one surface and discourage floaters. Gives crisper geometry and better meshes; too much flattens fine translucent detail.") \
    X(float, normal_distortion_reg, 0.0f, "model", "",                        \
      "The same idea as depth distortion, applied to surface direction. Each pixel settles on one consistent orientation instead of several.") \
    X(float, rgb_distortion_reg, 0.0f, "model", "",                           \
      "Encourage each pixel's color to come from one surface rather than blended layers. Helps discourage false transparency.") \
    X(int, distortion_reg_warmup, 6000, "model", "",                          \
      "How many steps the distortion penalties take to reach full strength. Ramping in lets coarse structure form before geometry is tightened.") \
    X(float, normal_reg_weight, 0.04f, "model", "",                           \
      "Encourage splats to lie flat along the surfaces they represent. Higher gives cleaner geometry and better meshes; too high flattens fine detail.") \
    X(int, normal_reg_warmup, 6000, "model", "",                              \
      "How many steps the surface-alignment penalty takes to reach full strength.") \
    X(float, alpha_reg_weight, 0.0f, "model", "",                             \
      "Push each pixel's coverage toward fully solid or fully empty instead of half transparent. Clears haze and gives cleaner background cutouts, at the risk of less stable training.") \
    X(int, alpha_reg_warmup, 12000, "model", "",                              \
      "How many steps the coverage penalty takes to reach full strength.")    \
    X(int, reg_warmup_length, 0, "model", "",                                 \
      "Hold the depth, normal and coverage penalties off for this many steps. Lets the scene take shape before geometry constraints start pulling on it.") \
    X(bool, apply_loss_for_mask, false, "model", "",                          \
      "Whether masked-out pixels are ignored or trained as empty space. Off ignores them, which is how you hide distractions such as people, cars, or the black area outside a fisheye circle. On trains them as empty, which removes the background and leaves just the subject.") \
    X(float, alpha_loss_weight, 0.01f, "model", "",                           \
      "How firmly to clear splats out of areas the mask says should be empty. Raise it if background creeps back in.") \
    X(float, alpha_loss_weight_under, 0.0f, "model", "",                      \
      "How firmly to fill in areas the mask says should be solid but the render leaves empty. Raise it if holes appear in the subject.") \
    X(float, opacity_reg, 0.01f, "model", "",                                 \
      "Gently push splat opacity down so weak splats get recycled where they are needed more. Higher recycles more aggressively.") \
    X(float, scale_reg, 0.01f, "model", "",                                   \
      "Gently push splat size down. Keeps splats compact so detail stays local; too high leaves large flat areas underfilled.") \
    X(float, opacity_decay, 0.0f, "model", "",                                \
      "Fade every splat's opacity slightly each step, freeing weak ones for reuse. An alternative to opacity_reg that acts on all splats equally.") \
    X(float, scale_decay, 0.0f, "model", "",                                  \
      "Shrink every splat slightly each step. An alternative to scale_reg that acts on all splats equally.") \
    X(float, erank_reg, 0.0f, "model", "",                                    \
      "Discourage needle-like splats in favor of rounder ones. Reduces spiky artifacts and flicker as the camera moves.") \
    X(float, erank_reg_s3, 0.0f, "model", "",                                 \
      "Discourage splats from collapsing into flat sheets, keeping some thickness in every direction.") \
    X(float, quat_norm_reg, 0.01f, "model", "",                               \
      "Keep splat rotations well formed. Guards against numerical drift and rarely needs changing.") \
    X(float, sh_reg, 0.001f, "model", "",                                     \
      "Hold view-dependent color in check so it does not absorb shifts the per-photo correction should handle. Higher gives more consistent color from every angle and better results on unseen views; too high flattens genuine reflections.") \
    X(float, overexposure_reg, 0.0f, "model", "",                             \
      "Penalize splat colors that fall outside the displayable range. Keeps blown-out highlights from hiding inside the splats, which matters when exporting or meshing.") \
    X(int, supervision_warmup, 0, "model", "",                                \
      "Step at which AI-predicted depth and normals start guiding training. Waiting lets photometric detail establish itself first.") \
    X(float, depth_supervision_weight, 0.0f, "model", "",                     \
      "How strongly AI-predicted depth guides the geometry. Helps in textureless areas, but a heavily biased prediction can pull quality down.") \
    X(bool, input_depth_is_ray_depth, false, "model", "",                     \
      "Whether the supplied depth maps measure distance along the camera ray instead of distance straight ahead. Most AI-predicted depth is the latter, so leave this off; turn it on for very wide fisheye captures where straight-ahead depth is meaningless.") \
    X(float, normal_supervision_weight, 0.01f, "model", "",                   \
      "How strongly AI-predicted surface direction guides the geometry. Helps flat surfaces come out flat.") \
    X(float, mean_median_depth_weight, 0.0f, "model", "",                     \
      "Encourage a pixel's average depth and its most-solid depth to agree, which pulls splats onto a single surface. Useful when the result is destined for a mesh.") \
    X(float, median_depth_normal_reg_weight, 0.0f, "model", "",               \
      "Encourage surface direction to agree between the average depth and the most-solid depth. Another crispness dial for meshing.") \
    X(float, median_normal_supervision_weight, 0.0f, "model", "",             \
      "Match the surface direction at the most-solid depth to AI-predicted normals.") \
    X(float, median_render_normal_reg_weight, 0.0f, "model", "",              \
      "Match the surface direction at the most-solid depth to the rendered normals.") \
    X(int, median_warmup, 6000, "model", "",                                  \
      "How many steps the four most-solid-depth penalties take to reach full strength.") \
                                                                              \
    /* ==== optimizer -- learning rates and schedules ==== */                 \
    X(std::optional<int>, max_steps, std::nullopt, "optimizer", "",           \
      "Length of the learning-rate schedule, in steps. Leave unset to match num_iterations. Set it to keep the rates on their usual curve when training longer or shorter than the schedule was designed for.") \
    X(bool, use_scale_agnostic_mean, true, "optimizer", "",                   \
      "Make how fast splats move independent of how large the scene is, so one setting works across datasets. Turn off to have it scale with the scene, matching the original 3DGS behavior.") \
    X(bool, use_per_splat_bias_correction, true, "optimizer", "",             \
      "Give newly created splats a fresh start in the optimizer. They then move at full speed instead of inheriting the momentum of the splat they came from, so new detail sharpens up faster.") \
    X(float, means_lr, 0.000128f, "optimizer", "",                            \
      "How fast splats move through space. Higher rearranges geometry quickly but can jitter and blur; lower stays closer to the starting point cloud.") \
    X(std::optional<float>, means_lr_final, 1.6e-06f, "optimizer", "",        \
      "How fast splats move by the end of training. Positions ease to a stop so detail can settle. Set to none to keep the rate constant.") \
    X(float, scales_lr, 0.02f, "optimizer", "",                               \
      "How fast splats change size. Higher adapts coverage quickly; lower keeps sizes near where they started.") \
    X(std::optional<float>, scales_lr_final, 0.005f, "optimizer", "",         \
      "How fast splats change size by the end of training. Set to none to keep the rate constant.") \
    X(float, quats_lr, 0.0015f, "optimizer", "",                              \
      "How fast splats rotate, which sets how quickly they align themselves to surfaces.") \
    X(float, opacities_lr, 0.025f, "optimizer", "",                           \
      "How fast splat transparency changes. Higher clears haze and prunes faint splats sooner; lower gives weak structure more time to prove itself.") \
    X(float, features_dc_lr, 0.005f, "optimizer", "",                         \
      "How fast the base color of a splat changes.")                          \
    X(float, features_sh_lr, 0.00025f, "optimizer", "",                       \
      "How fast view-dependent color changes. Kept well below the base color rate so reflections do not run away with the color.") \
    X(float, background_dc_lr, 0.0025f, "optimizer", "",                      \
      "How fast the skybox base color changes. Only used with the `sh` background.") \
    X(float, background_sh_lr, 0.0005f, "optimizer", "",                      \
      "How fast skybox detail changes. Only used with the `sh` background.")  \
    X(float, bilagrid_lr, 0.002f, "optimizer", "",                            \
      "How fast the per-photo color correction adapts. Higher tracks exposure changes sooner but can start absorbing real detail. Ignored when use_adagrad_bilagrid_optim is on.") \
    X(std::optional<float>, bilagrid_lr_final, 0.0001f, "optimizer", "",      \
      "How fast the per-photo color correction adapts by the end of training. Set to none to keep the rate constant.") \
    X(int, bilagrid_lr_warmup, 1000, "optimizer", "",                         \
      "How many steps the per-photo color correction takes to reach full adaptation speed. Ramping in stops it from claiming color before the splats have any.") \
    X(float, bilagrid_depth_lr, 0.002f, "optimizer", "",                      \
      "How fast the per-photo depth correction adapts. Ignored when use_adagrad_bilagrid_optim is on.") \
    X(std::optional<float>, bilagrid_depth_lr_final, 0.0001f, "optimizer", "", \
      "How fast the per-photo depth correction adapts by the end of training. Set to none to keep the rate constant.") \
    X(int, bilagrid_depth_lr_warmup, 2000, "optimizer", "",                   \
      "How many steps the per-photo depth correction takes to reach full adaptation speed.") \
    X(float, bilagrid_normal_lr, 0.0005f, "optimizer", "",                    \
      "How fast the per-photo normal correction adapts. Ignored when use_adagrad_bilagrid_optim is on.") \
    X(std::optional<float>, bilagrid_normal_lr_final, 4e-05f, "optimizer", "", \
      "How fast the per-photo normal correction adapts by the end of training. Set to none to keep the rate constant.") \
    X(int, bilagrid_normal_lr_warmup, 2000, "optimizer", "",                  \
      "How many steps the per-photo normal correction takes to reach full adaptation speed.") \
    X(float, bilagrid_adagrad_lr, 0.04f, "optimizer", "",                     \
      "How fast the per-photo color correction adapts when use_adagrad_bilagrid_optim is on. This rate is constant, with no schedule or warmup.") \
    X(float, bilagrid_adagrad_depth_lr, 0.04f, "optimizer", "",               \
      "How fast the per-photo depth correction adapts when use_adagrad_bilagrid_optim is on. This rate is constant, with no schedule or warmup.") \
    X(float, bilagrid_adagrad_normal_lr, 0.01f, "optimizer", "",              \
      "How fast the per-photo normal correction adapts when use_adagrad_bilagrid_optim is on. This rate is constant, with no schedule or warmup.") \
    X(float, ppisp_lr, 0.002f, "optimizer", "",                               \
      "How fast the camera-effect model adapts. Ignored when use_adagrad_ppisp_optim is on.") \
    X(std::optional<float>, ppisp_lr_final, 2e-05f, "optimizer", "",          \
      "How fast the camera-effect model adapts by the end of training. Set to none to keep the rate constant.") \
    X(int, ppisp_lr_warmup, 500, "optimizer", "",                             \
      "How many steps the camera-effect model takes to reach full adaptation speed. Ramping in stops it from claiming brightness before the splats have any.") \
    X(float, ppisp_adagrad_lr, 0.1f, "optimizer", "",                         \
      "How fast the camera-effect model adapts when use_adagrad_ppisp_optim is on. This rate is constant, with no schedule or warmup.") \
    /* end */


// ===========================================================================
// The config struct, expanded from the table above
// ===========================================================================

struct TrainConfig {
#define SS_DECLARE_FIELD(type, member, default_, group, choices, help)         \
    type member = default_;
    SS_CONFIG_FIELDS(SS_DECLARE_FIELD)
#undef SS_DECLARE_FIELD
};

// Fields whose default ({}) is not a usable value. Checked after flag
// parsing, so --help still works without them.
#define SS_CONFIG_REQUIRED_FIELDS(X) \
    X(data) \
    /* end */


// ===========================================================================
// Presets -- named bundles of default overrides, selected as
// `spirula train <preset>`. "3dgs" is the base config and applies nothing.
// ===========================================================================

// Just the names, which are what the command line and the checkpoints spell.
// The label and the sentence explaining each one are in
// i18n/catalog/Train.h: this header is included nearly everywhere and has no
// business dragging a translation catalog behind it. Both consumers
// (app/cli/main.cpp, app/gui/GuiApp.cpp) static_assert that the two lists are
// the same length.
struct TrainPresetInfo { const char* name; };
inline constexpr TrainPresetInfo kTrainPresets[] = {
    {"3dgs"},
    {"360-camera"},
    {"in-the-wild"},
    {"linear-color"},
    {"synthetic"},
    {"meshing"},
    {"academic-baseline"},
};

// Returns false for an unknown preset name.
inline bool train_apply_preset(TrainConfig& c, const std::string& name) {
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
        c.min_init_fraction = 0.02f;
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
        c.max_screen_size = kTrainInf;
        c.max_world_size = kTrainInf;
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
