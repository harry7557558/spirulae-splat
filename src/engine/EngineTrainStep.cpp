// Engine fused training step: orchestrates set_camera + set_training + forward +
// loss/bwd + bilagrid forward/optim + ppisp forward/optim + splat optim + densify.

#include "engine/Engine.h"
#include "engine/EngineInternal.h"
#include "engine/EngineState.h"

#include <map>
#include <stdexcept>
#include <string>


// Per-camera batch slice of a TorchTensorView. Returns a view at the k-th
// row of the leading dim with shape [1, rest...]. Null TVs pass through
// unchanged (signaling "absent" GT modality).
static TorchTensorView _slice_tv_first_dim(const TorchTensorView& tv, int64_t k) {
    uint64_t base = std::get<0>(tv);
    if (base == 0) return tv;
    uint32_t esize = std::get<1>(tv);
    const auto& shape = std::get<2>(tv);
    if (shape.empty()) return tv;
    int64_t inner = 1;
    for (size_t i = 1; i < shape.size(); ++i) inner *= shape[i];
    uint64_t offset_bytes = (uint64_t)k * (uint64_t)inner * (uint64_t)esize;
    std::vector<int64_t> new_shape = shape;
    new_shape[0] = 1;
    return TorchTensorView(base + offset_bytes, esize, std::move(new_shape));
}


// Slice [k_start : k_start + k_count] along the leading dim. Used by the
// warp split-batch path to carve out K consecutive post-split cameras per
// sub-batch (one slice for the per-input-image table, another for the
// K-per-input post-camera table). Null TVs pass through.
static TorchTensorView _slice_tv_range_first_dim(const TorchTensorView& tv,
                                                 int64_t k_start, int64_t k_count) {
    uint64_t base = std::get<0>(tv);
    if (base == 0) return tv;
    uint32_t esize = std::get<1>(tv);
    const auto& shape = std::get<2>(tv);
    if (shape.empty()) return tv;
    int64_t inner = 1;
    for (size_t i = 1; i < shape.size(); ++i) inner *= shape[i];
    uint64_t offset_bytes = (uint64_t)k_start * (uint64_t)inner * (uint64_t)esize;
    std::vector<int64_t> new_shape = shape;
    new_shape[0] = k_count;
    return TorchTensorView(base + offset_bytes, esize, std::move(new_shape));
}


// Forward + bilagrid/PPISP forward + loss + raster/proj backward only.
// No optimizer / densify. Used both as a per-sub-batch step in split mode
// and as the inner of the single-batch path. Returns the loss_dict from
// engine_compute_loss_backward; callers aggregate across sub-batches.
static std::map<std::string, float> _engine_step_fwd_bwd_only(
    int step,
    std::string primitive,
    int sh_degree,
    bool packed,
    TorchTensorView bilagrid_cam_indices,
    const EngineStepConfig& cfg
) {
    engine().optim.use_fused_proj_bwd_optim = cfg.optim.use_fused_proj_bwd_optim;
    // Block-wise quantized gradient accumulation: only on the non-FPBO path
    // (FPBO already avoids fp32 grad buffers) and only when quantization is on.
    // Preserves split_batch via the projection-bwd kernel's requantize-accumulate.
    engine().grad.quantize_grad =
        !cfg.optim.use_fused_proj_bwd_optim && (cfg.optim.quantization_level != 0);

    _set_cur_cam_indices(bilagrid_cam_indices);

    // Pick the distortion channel set from the per-channel weights, restricted
    // to the channels this primitive renders (normal is honored only once a
    // normal-rendering primitive exists). 3dgs/mip and 3dgut all consume the
    // distortion gradient in their backward.
    const DistortionType dist_type = engine_distortion_type(
        engine_primitive_pixel_type(primitive),
        cfg.loss.weights[(int)LossWeightIndex::RgbDistReg],
        cfg.loss.weights[(int)LossWeightIndex::DepthDistReg],
        cfg.loss.weights[(int)LossWeightIndex::NormalDistReg]);
    forward_3dgs(primitive, sh_degree, packed, /*output_median=*/false, (int)dist_type);

    engine().ppisp.cur_run_before_bilagrid = cfg.ppisp.run_before_bilagrid;
    const bool bg_enabled = engine().bilagrid_rgb.enabled ||
                            engine().bilagrid_depth.enabled ||
                            engine().bilagrid_normal.enabled;
    if (engine().ppisp.cur_run_before_bilagrid) {
        if (engine().ppisp.enabled) engine_ppisp_forward(bilagrid_cam_indices);
        if (bg_enabled)             engine_bilagrid_forward(bilagrid_cam_indices);
    } else {
        if (bg_enabled)             engine_bilagrid_forward(bilagrid_cam_indices);
        if (engine().ppisp.enabled) engine_ppisp_forward(bilagrid_cam_indices);
    }

    return engine_compute_loss_backward(
        step, cfg.loss.weights, cfg.loss.w_ssim,
        cfg.loss.num_loss_scales, cfg.loss.loss_scale_min_pixels,
        cfg.loss.compute_loss_map,
        cfg.loss.loss_map_mode,
        cfg.loss.robust_edge_aware_quantile,
        cfg.loss.nms_falloff,
        cfg.loss.loss_map_normalize,
        cfg.loss.loss_map_clip_quantile,
        cfg.loss.loss_map_power,
        cfg.loss.loss_map_accum_mode,
        cfg.loss.overexposure_reg_weight,
        cfg.loss.color_shift_reg_weight,
        cfg.loss.color_shift_reg_beta);
}


// Splat / background / bilagrid / PPISP optimizer steps + densify. Consumes
// whatever is in the per-splat grad accumulators (with engine().optim.grad_scale
// applied inside the splat Adam kernel) and the per-camera bilagrid/PPISP grad
// buffers (filled by their backward hooks during the sub-batch loop). Returns
// the loss_dict augmented with TV / PPISP reg readouts and density counts.
static void _engine_step_optim_and_densify(
    int step, int max_steps,
    const EngineStepConfig& cfg,
    std::map<std::string, float>& loss_dict
) {
    engine_optim_step(step, cfg.optim);

    if (engine().background.enabled) {
        engine_background_optim_step(step, cfg.background);
    }

    if (engine().bilagrid_rgb.enabled || engine().bilagrid_depth.enabled ||
        engine().bilagrid_normal.enabled) {
        engine_bilagrid_optim_step(step, cfg.bilagrid);

        float* tv_buf = DevicePool::global().acquire<float>(PoolSlot::EngBgTvReadout, 3);
        _engine_bilagrid_tv_into(tv_buf);
        static AsyncReadout<float> tv_readout(3);
        const float* h_tv = tv_readout.read_previous();
        if (h_tv) {
            if (engine().bilagrid_rgb.enabled)    loss_dict["bilagrid_tv"]        = h_tv[0];
            if (engine().bilagrid_depth.enabled)  loss_dict["bilagrid_depth_tv"]  = h_tv[1];
            if (engine().bilagrid_normal.enabled) loss_dict["bilagrid_normal_tv"] = h_tv[2];
        }
        tv_readout.issue(tv_buf);
    }

    if (engine().ppisp.enabled) {
        engine_ppisp_optim_step(step, cfg.ppisp);

        float* losses_buf = _engine_ppisp_reg_loss_into(cfg.ppisp.reg_weights,
                                                       /*compute_grad=*/false);
        static AsyncReadout<float> ppisp_reg_readout((int)PPISPRegLossIndex::length);
        const float* h_losses = ppisp_reg_readout.read_previous();
        if (h_losses) {
            loss_dict["ppisp_reg_exposure_mean"]   = h_losses[(int)PPISPRegLossIndex::ExposureMean];
            loss_dict["ppisp_reg_vig_center"]      = h_losses[(int)PPISPRegLossIndex::VignettingCenter];
            loss_dict["ppisp_reg_vig_non_pos"]     = h_losses[(int)PPISPRegLossIndex::VignettingNonPositivity];
            loss_dict["ppisp_reg_vig_channel_var"] = h_losses[(int)PPISPRegLossIndex::VignettingChannelVariance];
            loss_dict["ppisp_reg_color_mean"]      = h_losses[(int)PPISPRegLossIndex::ColorMean];
            loss_dict["ppisp_reg_crf_channel_var"] = h_losses[(int)PPISPRegLossIndex::CRFChannelVariance];
        }
        ppisp_reg_readout.issue(losses_buf);
    }

    int num_added = engine_densify_step(step, max_steps, cfg.densify);

    loss_dict["num_added"] = (float)num_added;
    loss_dict["cur_num_splats"] = (float)engine().cur_num_splats;
    loss_dict["max_num_splats"] = (float)engine().max_num_splats;
}


// Inner training step: everything from "data is now in engine().gt + camera
// table" through optim + densify. Both engine_train_step (caller stages GT
// via set_training_data) and engine_train_step_warped (caller stages GT via
// set_training_data_warped) call this.
std::map<std::string, float> _engine_train_step_after_setup(
    int step, int max_steps,
    std::string primitive,
    int sh_degree,
    bool packed,
    TorchTensorView bilagrid_cam_indices,
    const EngineStepConfig& cfg
) {
    // Viewer: capture thumbnails of the GT for any new post-camera ids in
    // this batch. No-op when the viewer hasn't been initialized or the
    // host counter is already drained (all cams cached). Must run AFTER
    // _set_cur_cam_indices uploads the device indices and AFTER
    // engine().gt.rgb has been populated by set_training_data{,_warped}
    // (called by the public engine_train_step{,_warped} before us). The
    // sub-batched fwd/bwd helper re-uploads cam_indices per sub-batch, so
    // we run this once here (with the FULL batch indices) before that
    // happens.
    _set_cur_cam_indices(bilagrid_cam_indices);
    engine_viewer_capture_thumbnails(bilagrid_cam_indices);

    if (engine().background.enabled) {
        engine_set_background_step_params(cfg.background.seed,
                                          cfg.background.randomize_weight);
    }

    auto loss_dict = _engine_step_fwd_bwd_only(
        step, std::move(primitive), sh_degree, packed,
        bilagrid_cam_indices, cfg);

    _engine_step_optim_and_densify(step, max_steps, cfg, loss_dict);
    return loss_dict;
}


// The mode this step's raster backward will use -- None when nothing asked
// for a loss map, so the accumulator stays one lane wide.
static DensifyAccumMode _step_accum_mode(const EngineStepConfig& cfg) {
    return cfg.loss.compute_loss_map ? (DensifyAccumMode)cfg.loss.loss_map_accum_mode
                                     : DensifyAccumMode::None;
}

// Fold one sub-batch's raster-bwd accum_weight into `dst` exactly as the
// unsplit kernel folds a batch's cameras: Max reduces with atomicMax, Sum and
// Avg with atomicAdd (Avg over both planar lanes, divided later).
static void _fold_accum_weight(DeviceVector<float>& dst, int64_t num_splats) {
    if (engine().fwd.accum_weight.data_ptr() == nullptr) return;
    const DensifyAccumMode mode = engine().fwd.accum_mode;
    const int64_t n = num_splats * accum_lanes(mode);
    if (mode == DensifyAccumMode::Max)
        float_max_into(dst, engine().fwd.accum_weight, n);
    else
        float_add_into(dst, engine().fwd.accum_weight, n);
}


// Sub-batched (one-camera-per-sub-batch) train step. Pre-conditions:
//   * !cfg.optim.use_fused_proj_bwd_optim   (FPBO path is incompatible)
//   * !cfg.optim.use_color_trust_region     (TR kernels not yet threaded)
//   * non-warped entry path                 (warped path handled separately)
//
// Stages GT + camera state for one camera per iteration, runs the full
// forward + loss + raster/proj backward chain, atomic-adds into the
// per-splat grad accumulators, then runs a SINGLE optimizer + densify pass
// at the end. The splat Adam kernel reads grad_scale = 1/B from engine
// state so per-image grad magnitude vs regularization weight stays
// batch-size invariant.

static std::map<std::string, float> _engine_train_step_split_one_per_camera(
    int step, int max_steps,
    std::string primitive,
    int sh_degree,
    bool packed,
    int width, int height, std::string camera_model, std::string distortion,
    TorchTensorView viewmats,      // [B, 4, 4]
    TorchTensorView intrins,       // [B, 4]
    TorchTensorView dist_coeffs,   // [B, K]
    TorchTensorView gt_rgb,        // [B, H, W, 3] or null
    TorchTensorView gt_depth,      // [B, Hd, Wd, 1] or null
    TorchTensorView gt_normal,     // [B, Hn, Wn, 3] or null
    TorchTensorView gt_alpha,      // [B, Hm, Wm, 1] or null (mask)
    TorchTensorView bilagrid_cam_indices, // [B] int32 or null
    const EngineStepConfig& cfg
) {
    if (cfg.optim.use_fused_proj_bwd_optim) {
        throw std::runtime_error(
            "split_batch is not compatible with "
            "use_fused_proj_bwd_optim (FPBO folds grad into the optim "
            "kernel; would need a separate non-fused-proj-bwd path).");
    }

    const auto& viewmats_shape = std::get<2>(viewmats);
    if (viewmats_shape.empty())
        throw std::runtime_error("split_batch: empty viewmats");
    int64_t B = viewmats_shape[0];
    if (B <= 0)
        throw std::runtime_error("split_batch: zero cameras in batch");

    // Single-camera batch -> no point splitting; fall through to the
    // standard path with grad_scale = 1.
    if (B == 1) {
        set_camera_params(width, height, camera_model, distortion,
                          viewmats, intrins, dist_coeffs);
        set_training_data(gt_rgb, gt_depth, gt_normal, gt_alpha,
                          cfg.loss.input_depth_is_ray_depth);
        return _engine_train_step_after_setup(
            step, max_steps, std::move(primitive), sh_degree, packed,
            bilagrid_cam_indices, cfg);
    }

    // Background per-iter params are set once for the full step (the
    // background blend reads them per-camera inside the forward kernel).
    if (engine().background.enabled) {
        engine_set_background_step_params(cfg.background.seed,
                                          cfg.background.randomize_weight);
    }

    // color_shift_reg EMA: the per-call decay is beta. In single-batch mode
    // that decay is applied ONCE per train step. In sub-batched mode the
    // hook is called B times per step, so to keep the effective per-step
    // decay identical, we apply beta^(1/B) inside each sub-batch. The
    // composition beta_sub^B = beta gives the same end-of-step EMA value
    // (modulo intra-step ordering which is FP-non-associative). This rebinds
    // through the local cfg copy below; the caller's cfg is untouched.
    EngineStepConfig cfg_sub = cfg;
    const float beta_orig = cfg.loss.color_shift_reg_beta;
    if (beta_orig > 0.0f && beta_orig < 1.0f) {
        cfg_sub.loss.color_shift_reg_beta = powf(beta_orig, 1.0f / (float)B);
    }

    // The raster_bwd pool buffer ("raster_bwd.accum_weight") is zeroed and
    // overwritten on every call, so engine().fwd.accum_weight would otherwise
    // carry only the LAST sub-batch. Folded here, then re-pointed for densify.
    int64_t N_splats = engine().cur_num_splats;
    DeviceVector<float> accum_weight_sum;
    accum_weight_sum.resize(PoolSlot::EngSubbatchAccumWeightSum,
                            N_splats * accum_lanes(_step_accum_mode(cfg)));
    accum_weight_sum.zero();
    // Drop the previous step's re-pointed view (set below for densify): when
    // this step's raster bwd doesn't produce accum_weight (e.g. score_blend_
    // world_grad == 1), the stale view would otherwise pass the null check in
    // the loop and dangle once the pool moves (densify growth, viewer render).
    engine().fwd.accum_weight = DeviceVector<float>();
    engine().fwd.accum_mode = DensifyAccumMode::None;

    std::map<std::string, float> agg;

    for (int64_t k = 0; k < B; ++k) {
        // First sub-batch: zero grad. Subsequent: skip the zero so atomic
        // adds accumulate across sub-batches.
        engine().optim.skip_grad_zero = (k > 0);

        TorchTensorView vmt_k = _slice_tv_first_dim(viewmats,    k);
        TorchTensorView itr_k = _slice_tv_first_dim(intrins,     k);
        TorchTensorView dst_k = _slice_tv_first_dim(dist_coeffs, k);
        TorchTensorView rgb_k = _slice_tv_first_dim(gt_rgb,      k);
        TorchTensorView dep_k = _slice_tv_first_dim(gt_depth,    k);
        TorchTensorView nrm_k = _slice_tv_first_dim(gt_normal,   k);
        TorchTensorView aph_k = _slice_tv_first_dim(gt_alpha,    k);
        TorchTensorView bgi_k = _slice_tv_first_dim(bilagrid_cam_indices, k);

        set_camera_params(width, height, camera_model, distortion, vmt_k, itr_k, dst_k);
        set_training_data(rgb_k, dep_k, nrm_k, aph_k,
                          cfg.loss.input_depth_is_ray_depth);

        // Thumbnails: capture per-sub-batch with the sliced cam index. The
        // viewer host counter naturally throttles re-captures, so doing this
        // B times per step is safe (and matches the unsplit behavior of
        // capturing once over all cams in the batch).
        _set_cur_cam_indices(bgi_k);
        engine_viewer_capture_thumbnails(bgi_k);

        auto ld = _engine_step_fwd_bwd_only(
            step, primitive, sh_degree, packed, bgi_k, cfg_sub);

        // Accumulate per-sub-batch accum_weight (per-splat densification
        // score from raster_bwd). engine().fwd.accum_weight points at the
        // raster_bwd pool buffer for this sub-batch; we sum it into our
        // persistent buffer before raster_bwd overwrites it next iter.
        if (engine().fwd.accum_weight.data_ptr() != nullptr) {
            _fold_accum_weight(accum_weight_sum, N_splats);
        }

        // Aggregate loss values across sub-batches as a mean.
        if (k == 0) {
            agg = std::move(ld);
        } else {
            for (auto& [key, val] : ld) {
                auto it = agg.find(key);
                if (it == agg.end()) agg.emplace(key, val);
                else                 it->second += val;
            }
        }
    }
    for (auto& [key, val] : agg) val /= (float)B;

    // Reduced across sub-batches the way the unsplit raster_bwd reduces the
    // cameras of one launch: atomicMax, not a sum.
    engine().fwd.accum_weight = accum_weight_sum;

    // Single end-of-step optim + densify pass. The splat Adam kernels read
    // grad_scale from engine state and divide the data gradient by B before
    // adding regularization terms; the fused zero_grad side effect leaves
    // the per-splat grad buffers at 0 for the next train step.
    engine().optim.skip_grad_zero      = false;
    engine().optim.grad_scale          = 1.0f / (float)B;
    engine().optim.zero_grad_in_optim  = true;

    _engine_step_optim_and_densify(step, max_steps, cfg, agg);

    // Reset back to single-batch defaults so any non-split train step that
    // follows behaves exactly as before. (engine_optim_step reads these
    // each call, so a stale 1/B here would silently corrupt a follow-up
    // single-batch step.)
    engine().optim.grad_scale         = 1.0f;
    engine().optim.zero_grad_in_optim = false;

    return agg;
}


// Heterogeneous (mixed-resolution) training step. Generalizes the one-camera-
// per-sub-batch split above to an OUTER loop over sub-batches of differing
// (W, H, model): each camera of each sub-batch runs the full fwd + loss +
// raster/proj backward and atomic-adds into the per-splat grad accumulators,
// then a SINGLE optim + densify pass applies grad_scale = 1 / (total cameras).
// This is how cross-group-packed steps (many small resolution groups) are
// consumed. Non-warp only; the DataManager never packs K > 1 groups together.
std::map<std::string, float> engine_train_step_hetero(
    int step, int max_steps,
    std::string primitive,
    int sh_degree,
    bool packed,
    const std::vector<HeteroSubBatch>& subs,
    const EngineStepConfig& cfg
) {
    if (cfg.optim.use_fused_proj_bwd_optim) {
        throw std::runtime_error(
            "engine_train_step_hetero requires the non-fused (split_batch) "
            "path; FPBO folds grads into the optim kernel and cannot span "
            "mixed-resolution sub-batches. Keep split_batch enabled.");
    }
    if (subs.empty())
        throw std::runtime_error("engine_train_step_hetero: empty step");

    // Total cameras across all sub-batches drives grad normalization + the
    // per-call color_shift EMA decay.
    int64_t total_cams = 0;
    for (const auto& s : subs) {
        const auto& vs = std::get<2>(s.viewmats);
        int64_t nb = vs.empty() ? 0 : vs[0];
        if (nb <= 0)
            throw std::runtime_error("engine_train_step_hetero: empty sub-batch");
        total_cams += nb;
    }

    if (engine().background.enabled) {
        engine_set_background_step_params(cfg.background.seed,
                                          cfg.background.randomize_weight);
    }

    // color_shift_reg EMA: the hook runs once per fwd/bwd, i.e. total_cams
    // times per step, so use beta^(1/total_cams) to preserve the end-of-step
    // EMA value (matching the single-batch and split-per-camera paths).
    EngineStepConfig cfg_sub = cfg;
    const float beta_orig = cfg.loss.color_shift_reg_beta;
    if (beta_orig > 0.0f && beta_orig < 1.0f) {
        cfg_sub.loss.color_shift_reg_beta = powf(beta_orig, 1.0f / (float)total_cams);
    }

    int64_t N_splats = engine().cur_num_splats;
    DeviceVector<float> accum_weight_sum;
    accum_weight_sum.resize(PoolSlot::EngSubbatchAccumWeightSum,
                            N_splats * accum_lanes(_step_accum_mode(cfg)));
    accum_weight_sum.zero();
    // Drop the previous step's re-pointed view (see per-sub-batch path).
    engine().fwd.accum_weight = DeviceVector<float>();
    engine().fwd.accum_mode = DensifyAccumMode::None;

    std::map<std::string, float> agg;
    int64_t kcam = 0;   // global camera index across the whole step

    for (const auto& s : subs) {
        const auto& vs = std::get<2>(s.viewmats);
        int64_t nb = vs[0];
        for (int64_t c = 0; c < nb; ++c) {
            // Zero grads + radii only on the VERY FIRST camera of the step;
            // every subsequent camera (any sub-batch) atomic-adds into them.
            engine().optim.skip_grad_zero = (kcam > 0);

            TorchTensorView vmt_c = _slice_tv_first_dim(s.viewmats,    c);
            TorchTensorView itr_c = _slice_tv_first_dim(s.intrins,     c);
            TorchTensorView dst_c = _slice_tv_first_dim(s.dist_coeffs, c);
            TorchTensorView rgb_c = _slice_tv_first_dim(s.gt_rgb,      c);
            TorchTensorView dep_c = _slice_tv_first_dim(s.gt_depth,    c);
            TorchTensorView nrm_c = _slice_tv_first_dim(s.gt_normal,   c);
            TorchTensorView aph_c = _slice_tv_first_dim(s.gt_alpha,    c);
            TorchTensorView bgi_c = _slice_tv_first_dim(s.bilagrid_cam_indices, c);

            set_camera_params(s.width, s.height, s.camera_model, s.distortion,
                              vmt_c, itr_c, dst_c);
            set_training_data(rgb_c, dep_c, nrm_c, aph_c,
                              cfg.loss.input_depth_is_ray_depth);

            _set_cur_cam_indices(bgi_c);
            engine_viewer_capture_thumbnails(bgi_c);

            auto ld = _engine_step_fwd_bwd_only(
                step, primitive, sh_degree, packed, bgi_c, cfg_sub);

            if (engine().fwd.accum_weight.data_ptr() != nullptr) {
                _fold_accum_weight(accum_weight_sum, N_splats);
            }

            if (kcam == 0) {
                agg = std::move(ld);
            } else {
                for (auto& [key, val] : ld) {
                    auto it = agg.find(key);
                    if (it == agg.end()) agg.emplace(key, val);
                    else                 it->second += val;
                }
            }
            ++kcam;
        }
    }
    for (auto& [key, val] : agg) val /= (float)total_cams;

    // Densify sees the per-splat score maxed across every camera of the step.
    engine().fwd.accum_weight = accum_weight_sum;

    engine().optim.skip_grad_zero      = false;
    engine().optim.grad_scale          = 1.0f / (float)total_cams;
    engine().optim.zero_grad_in_optim  = true;

    _engine_step_optim_and_densify(step, max_steps, cfg, agg);

    // Restore single-batch defaults for any follow-up non-split step.
    engine().optim.grad_scale         = 1.0f;
    engine().optim.zero_grad_in_optim = false;

    return agg;
}


// Standard (non-warp) public entry: stage camera + GT then run the inner step.
std::map<std::string, float> engine_train_step(
    int step, int max_steps,
    std::string primitive,
    int sh_degree,
    bool packed,
    int width, int height, std::string camera_model, std::string distortion,
    TorchTensorView viewmats,
    TorchTensorView intrins,
    TorchTensorView dist_coeffs,
    TorchTensorView gt_rgb,
    TorchTensorView gt_depth,
    TorchTensorView gt_normal,
    TorchTensorView gt_alpha,
    TorchTensorView bilagrid_cam_indices,
    const EngineStepConfig& cfg
) {
    if (cfg.optim.split_batch) {
        return _engine_train_step_split_one_per_camera(
            step, max_steps, std::move(primitive), sh_degree, packed,
            width, height, std::move(camera_model), std::move(distortion),
            viewmats, intrins, dist_coeffs,
            gt_rgb, gt_depth, gt_normal, gt_alpha,
            bilagrid_cam_indices, cfg);
    }
    set_camera_params(width, height, camera_model, distortion, viewmats, intrins, dist_coeffs);
    set_training_data(gt_rgb, gt_depth, gt_normal, gt_alpha,
                      cfg.loss.input_depth_is_ray_depth);
    return _engine_train_step_after_setup(
        step, max_steps, std::move(primitive), sh_degree, packed,
        bilagrid_cam_indices, cfg);
}


// Warped variant of the split-batch dispatcher. Splits per INPUT image
// (B_in chunks of K post-split cameras each). Numerically: each sub-batch's
// per-pixel loss kernel sees K cameras, normalizes by 1/(K*H*W). Accumulated
// across B_in sub-batches the per-splat grad ends up B_in x what an unsplit
// run on B_in*K post-cams would produce; grad_scale = 1/B_in in the optim
// step's splat Adam kernels brings it back to the unsplit-equivalent
// magnitude (mod atomicAdd reordering + the radii zero-on-first-sub-batch
// fix in EngineForward.cpp).
static std::map<std::string, float> _engine_train_step_split_warped(
    int step, int max_steps,
    std::string primitive, int sh_degree, bool packed,
    int out_W, int out_H,
    TorchTensorView post_viewmats,
    TorchTensorView post_intrins,
    TorchTensorView post_dist_coeffs,
    std::string input_camera_model,
    std::string input_distortion,
    int B_in, int in_H, int in_W,
    int K,
    TorchTensorView input_intrins,
    TorchTensorView input_dist_coeffs,
    TorchTensorView input_source_models,
    TorchTensorView input_source_params,
    TorchTensorView gt_rgb_byte,
    TorchTensorView gt_alpha_byte,
    int mask_in_H, int mask_in_W,
    TorchTensorView gt_depth_byte,
    int depth_in_H, int depth_in_W,
    TorchTensorView gt_normal_byte,
    int normal_in_H, int normal_in_W,
    uint64_t axes_dev,
    TorchTensorView bilagrid_cam_indices,
    const EngineStepConfig& cfg
) {
    if (cfg.optim.use_fused_proj_bwd_optim) {
        throw std::runtime_error(
            "split_batch is not compatible with use_fused_proj_bwd_optim "
            "(FPBO folds grad into the optim kernel).");
    }
    if (B_in <= 0)
        throw std::runtime_error("split_batch (warped): zero input images");
    if (K <= 0)
        throw std::runtime_error("split_batch (warped): K must be positive");

    // Cubemap faces are canonical pinhole. K == 1 is the re-distort path,
    // whose destination is the camera the parser fitted -- the input one.
    const std::string post_model = K > 1 ? "PINHOLE" : input_camera_model;
    const std::string post_dist  = K > 1 ? "NONE"    : input_distortion;

    // Single input image (B_in == 1) -> nothing to split; route through the
    // standard warped path so K cams are processed in one shot. grad_scale
    // stays at 1.0 (the normal post-batch normalization is correct).
    if (B_in == 1) {
        set_camera_params(out_W, out_H, post_model, post_dist,
                          post_viewmats, post_intrins, post_dist_coeffs);
        set_training_data_warped(input_camera_model, input_distortion,
                                 B_in, in_H, in_W, K, out_H, out_W,
                                 gt_rgb_byte, gt_alpha_byte,
                                 mask_in_H, mask_in_W,
                                 gt_depth_byte, depth_in_H, depth_in_W,
                                 gt_normal_byte, normal_in_H, normal_in_W,
                                 cfg.loss.input_depth_is_ray_depth,
                                 input_intrins, input_dist_coeffs,
                                 input_source_models, input_source_params,
                                 axes_dev);
        return _engine_train_step_after_setup(
            step, max_steps, std::move(primitive), sh_degree, packed,
            bilagrid_cam_indices, cfg);
    }

    if (engine().background.enabled) {
        engine_set_background_step_params(cfg.background.seed,
                                          cfg.background.randomize_weight);
    }

    EngineStepConfig cfg_sub = cfg;
    const float beta_orig = cfg.loss.color_shift_reg_beta;
    if (beta_orig > 0.0f && beta_orig < 1.0f) {
        cfg_sub.loss.color_shift_reg_beta = powf(beta_orig, 1.0f / (float)B_in);
    }

    int64_t N_splats = engine().cur_num_splats;
    DeviceVector<float> accum_weight_sum;
    accum_weight_sum.resize(PoolSlot::EngSubbatchAccumWeightSum,
                            N_splats * accum_lanes(_step_accum_mode(cfg)));
    accum_weight_sum.zero();
    // Drop the previous step's re-pointed view (see per-sub-batch path).
    engine().fwd.accum_weight = DeviceVector<float>();
    engine().fwd.accum_mode = DensifyAccumMode::None;

    std::map<std::string, float> agg;

    for (int64_t k = 0; k < B_in; ++k) {
        engine().optim.skip_grad_zero = (k > 0);

        // Per-input-image slices (B_in axis).
        TorchTensorView rgb_k   = _slice_tv_first_dim(gt_rgb_byte,        k);
        TorchTensorView mask_k  = _slice_tv_first_dim(gt_alpha_byte,      k);
        TorchTensorView dep_k   = _slice_tv_first_dim(gt_depth_byte,      k);
        TorchTensorView nrm_k   = _slice_tv_first_dim(gt_normal_byte,     k);
        TorchTensorView i_itr_k = _slice_tv_first_dim(input_intrins,      k);
        TorchTensorView i_dst_k = _slice_tv_first_dim(input_dist_coeffs,  k);
        TorchTensorView i_src_m_k = std::get<0>(input_source_models) == 0
            ? input_source_models : _slice_tv_first_dim(input_source_models, k);
        TorchTensorView i_src_p_k = std::get<0>(input_source_params) == 0
            ? input_source_params : _slice_tv_first_dim(input_source_params, k);

        // Per-input-image POST-split slices (B_post axis, K consecutive rows
        // per input).
        TorchTensorView p_vmt_k = _slice_tv_range_first_dim(post_viewmats,     k * K, K);
        TorchTensorView p_itr_k = _slice_tv_range_first_dim(post_intrins,      k * K, K);
        TorchTensorView p_dst_k = _slice_tv_range_first_dim(post_dist_coeffs,  k * K, K);
        TorchTensorView bgi_k   = _slice_tv_range_first_dim(bilagrid_cam_indices, k * K, K);

        set_camera_params(out_W, out_H, post_model, post_dist,
                          p_vmt_k, p_itr_k, p_dst_k);
        set_training_data_warped(input_camera_model, input_distortion,
                                 /*B_in=*/1, in_H, in_W, K, out_H, out_W,
                                 rgb_k, mask_k,
                                 mask_in_H, mask_in_W,
                                 dep_k, depth_in_H, depth_in_W,
                                 nrm_k, normal_in_H, normal_in_W,
                                 cfg.loss.input_depth_is_ray_depth,
                                 i_itr_k, i_dst_k, i_src_m_k, i_src_p_k,
                                 axes_dev);

        _set_cur_cam_indices(bgi_k);
        engine_viewer_capture_thumbnails(bgi_k);

        auto ld = _engine_step_fwd_bwd_only(
            step, primitive, sh_degree, packed, bgi_k, cfg_sub);

        if (engine().fwd.accum_weight.data_ptr() != nullptr) {
            _fold_accum_weight(accum_weight_sum, N_splats);
        }

        if (k == 0) {
            agg = std::move(ld);
        } else {
            for (auto& [key, val] : ld) {
                auto it = agg.find(key);
                if (it == agg.end()) agg.emplace(key, val);
                else                 it->second += val;
            }
        }
    }
    for (auto& [key, val] : agg) val /= (float)B_in;

    engine().fwd.accum_weight = accum_weight_sum;

    // grad_scale = 1/B_in (number of sub-batches), NOT 1/B_post: each
    // sub-batch's loss kernel already normalizes by 1/(K*H*W), so summing
    // B_in sub-batches over-counts by exactly B_in relative to the unsplit
    // B_post-camera batch (whose per-pixel grad is 1/(B_post*H*W)).
    engine().optim.skip_grad_zero      = false;
    engine().optim.grad_scale          = 1.0f / (float)B_in;
    engine().optim.zero_grad_in_optim  = true;

    _engine_step_optim_and_densify(step, max_steps, cfg, agg);

    engine().optim.grad_scale         = 1.0f;
    engine().optim.zero_grad_in_optim = false;

    return agg;
}


// Warp public entry: stage POST-split camera params at (out_W, out_H) and
// the byte GT through the fused warp kernel, then run the inner step.
std::map<std::string, float> engine_train_step_warped(
    int step, int max_steps,
    std::string primitive, int sh_degree, bool packed,
    // POST-split camera setup (B*K cameras, all PINHOLE).
    int out_W, int out_H,
    TorchTensorView post_viewmats,        // [B*K, 4, 4]
    TorchTensorView post_intrins,         // [B*K, 4]
    TorchTensorView post_dist_coeffs,     // [B*K, 8]
    // Warp inputs (operate on INPUT shape).
    std::string input_camera_model,       // "FISHEYE" / "EQUISOLID" / "EQUIRECTANGULAR"
    std::string input_distortion,
    int B_in, int in_H, int in_W,
    int K,
    TorchTensorView input_intrins,        // [B_in, 4]
    TorchTensorView input_dist_coeffs,    // [B_in, 8]
    TorchTensorView input_source_models,  // [B_in] int32 (nullable)
    TorchTensorView input_source_params,  // [B_in, 16]
    TorchTensorView gt_rgb_byte,          // [B_in, in_H, in_W, 3] u8/u16
    TorchTensorView gt_alpha_byte,        // [B_in, mask_in_H, mask_in_W, 1] u8 (nullable)
    int mask_in_H, int mask_in_W,
    TorchTensorView gt_depth_byte,        // [B_in, depth_in_H, depth_in_W, 1] u16/f32 (nullable)
    int depth_in_H, int depth_in_W,
    TorchTensorView gt_normal_byte,       // [B_in, normal_in_H, normal_in_W, 3] u8/f32 (nullable)
    int normal_in_H, int normal_in_W,
    uint64_t axes_dev,                    // device float ptr [K, 3, 3]
    // Bilagrid cam indices in the POST-split index space.
    TorchTensorView bilagrid_cam_indices,
    const EngineStepConfig& cfg
) {
    if (cfg.optim.split_batch) {
        return _engine_train_step_split_warped(
            step, max_steps, std::move(primitive), sh_degree, packed,
            out_W, out_H, post_viewmats, post_intrins, post_dist_coeffs,
            std::move(input_camera_model), std::move(input_distortion), B_in, in_H, in_W, K,
            input_intrins, input_dist_coeffs,
            input_source_models, input_source_params,
            gt_rgb_byte, gt_alpha_byte, mask_in_H, mask_in_W,
            gt_depth_byte, depth_in_H, depth_in_W,
            gt_normal_byte, normal_in_H, normal_in_W,
            axes_dev, bilagrid_cam_indices, cfg);
    }
    // Cubemap faces are canonical pinhole; K == 1 (re-distort) keeps the
    // camera the parser fitted, which is the input one.
    set_camera_params(out_W, out_H,
                      K > 1 ? "PINHOLE" : input_camera_model,
                      K > 1 ? "NONE"    : input_distortion,
                      post_viewmats, post_intrins, post_dist_coeffs);
    // GT is warped on the fly into a float [B*K, out_H, out_W, 3] buffer.
    set_training_data_warped(input_camera_model, input_distortion,
                             B_in, in_H, in_W, K, out_H, out_W,
                             gt_rgb_byte, gt_alpha_byte,
                             mask_in_H, mask_in_W,
                             gt_depth_byte, depth_in_H, depth_in_W,
                             gt_normal_byte, normal_in_H, normal_in_W,
                             cfg.loss.input_depth_is_ray_depth,
                             input_intrins, input_dist_coeffs,
                                 input_source_models, input_source_params,
                                 axes_dev);
    return _engine_train_step_after_setup(
        step, max_steps, std::move(primitive), sh_degree, packed,
        bilagrid_cam_indices, cfg);
}
