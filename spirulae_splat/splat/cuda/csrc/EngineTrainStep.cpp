// Engine fused training step: orchestrates set_camera + set_training + forward +
// loss/bwd + bilagrid forward/optim + ppisp forward/optim + splat optim + densify.

#include "Engine.h"
#include "EngineInternal.h"
#include "EngineState.h"

#include <map>
#include <string>


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
    // Propagate the fused projection-bwd+optim flag to engine state BEFORE
    // forward/loss so engine_compute_loss_backward can skip projection_*_backward
    // and stash v_splats_s for the fused optim call below.
    engine().optim.use_fused_proj_bwd_optim = cfg.optim.use_fused_proj_bwd_optim;

    // Populate the shared cam_indices buffer + background per-iter params
    // BEFORE the forward pass. The forward_3dgs path itself runs the
    // background blend (needs cam_indices for SH-mode rotation gather), so
    // both must be in place by the time we launch projection.
    _set_cur_cam_indices(bilagrid_cam_indices);
    if (engine().background.enabled) {
        engine_set_background_step_params(cfg.background.seed,
                                          cfg.background.randomize_weight);
    }

    // Forward (writes to pool; no D->H). Includes the background blend in
    // place on fwd.renders.rgb when an engine background mode is active.
    forward_3dgs(primitive, sh_degree, packed);

    // Bilagrid forward (between rendering and loss). No-op when disabled.
    if (engine().bilagrid_rgb.enabled || engine().bilagrid_depth.enabled ||
        engine().bilagrid_normal.enabled) {
        engine_bilagrid_forward(bilagrid_cam_indices);
    }

    // PPISP forward (in place on rendered RGB, AFTER bilagrid). No-op when
    // disabled.
    if (engine().ppisp.enabled) {
        engine_ppisp_forward(bilagrid_cam_indices);
    }

    // Stash the color-shift regularizer per-step args for the bilagrid bwd
    // hook to read. The hook runs inside engine_compute_loss_backward and
    // has no other access to BilagridStepConfig.
    if (engine().bilagrid_rgb.enabled) {
        engine().bilagrid_rgb.cur_shift_reg_weight = cfg.bilagrid.shift_reg_weight_rgb;
        engine().bilagrid_rgb.cur_shift_reg_beta   = cfg.bilagrid.shift_reg_beta_rgb;
    }

    // Loss + backward (reads pool, writes pool grads; D->H only for scalar loss values)
    std::map<std::string, float> loss_dict = engine_compute_loss_backward(
        step, cfg.loss.weights, cfg.loss.w_ssim,
        cfg.loss.num_loss_scales, cfg.loss.compute_loss_map,
        cfg.loss.structure_only_loss_map,
        cfg.loss.overexposure_reg_weight);

    // Optimizer (in-place on pool buffers, no copies)
    engine_optim_step(step, cfg.optim);

    // Background SH optim step (no-op for Noise mode or when train_color=false).
    if (engine().background.enabled) {
        engine_background_optim_step(step, cfg.background);
    }

    // Bilagrid Adam step + TV regularization (after splat optim, before densify).
    if (engine().bilagrid_rgb.enabled || engine().bilagrid_depth.enabled ||
        engine().bilagrid_normal.enabled) {
        engine_bilagrid_optim_step(step, cfg.bilagrid);

        // Read post-update TV losses for verbose display. Async pattern:
        // read previous iter's slot now (cheap event-sync), queue this iter's
        // D->H, never block the host on cudaMemcpy.
        float* tv_buf = DevicePool::global().acquire<float>("eng.bg.tv_readout", 3);
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

    // PPISP Adam step + regularization (after splat optim, before densify).
    // Also read post-update regularization losses for verbose display.
    if (engine().ppisp.enabled) {
        engine_ppisp_optim_step(step, cfg.ppisp);

        float* losses_buf = _engine_ppisp_reg_loss_into(cfg.ppisp.reg_weights,
                                                       /*compute_grad=*/false);
        // Async readout: previous iter's slot now, queue this iter's D->H.
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

    // Densify (in-place on pool buffers, no copies)
    int num_added = engine_densify_step(step, max_steps, cfg.densify);

    loss_dict["num_added"] = (float)num_added;
    loss_dict["cur_num_splats"] = (float)engine().cur_num_splats;
    loss_dict["max_num_splats"] = (float)engine().max_num_splats;
    return loss_dict;
}


// Standard (non-warp) public entry: stage camera + GT then run the inner step.
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
) {
    set_camera_params(width, height, camera_model, viewmats, intrins, dist_coeffs);
    set_training_data(gt_rgb, gt_depth, gt_normal, gt_alpha);
    return _engine_train_step_after_setup(
        step, max_steps, std::move(primitive), sh_degree, packed,
        bilagrid_cam_indices, cfg);
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
    TorchTensorView post_dist_coeffs,     // [B*K, 10]
    // Warp inputs (operate on INPUT shape).
    std::string input_camera_model,       // "FISHEYE" / "EQUISOLID" / "EQUIRECTANGULAR"
    int B_in, int in_H, int in_W,
    int K,
    TorchTensorView input_intrins,        // [B_in, 4]
    TorchTensorView input_dist_coeffs,    // [B_in, 10]
    TorchTensorView gt_rgb_byte,          // [B_in, in_H, in_W, 3] u8/u16
    TorchTensorView gt_alpha_byte,        // [B_in, mask_in_H, mask_in_W, 1] u8 (nullable)
    int mask_in_H, int mask_in_W,
    uint64_t axes_dev,                    // device float ptr [K, 3, 3]
    // Bilagrid cam indices in the POST-split index space.
    TorchTensorView bilagrid_cam_indices,
    const EngineStepConfig& cfg
) {
    // Camera table is set up at POST-split resolution + PINHOLE.
    set_camera_params(out_W, out_H, "PINHOLE",
                      post_viewmats, post_intrins, post_dist_coeffs);
    // GT is warped on the fly into a float [B*K, out_H, out_W, 3] buffer.
    set_training_data_warped(input_camera_model,
                             B_in, in_H, in_W, K, out_H, out_W,
                             gt_rgb_byte, gt_alpha_byte,
                             mask_in_H, mask_in_W,
                             input_intrins, input_dist_coeffs, axes_dev);
    return _engine_train_step_after_setup(
        step, max_steps, std::move(primitive), sh_degree, packed,
        bilagrid_cam_indices, cfg);
}
