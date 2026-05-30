// Engine loss + combined backward pass.

#include "Engine.h"
#include "EngineCommon.h"
#include "EngineInternal.h"
#include "EngineState.h"

#include "BilagridBindings.h"
#include "ProjectionBwd.cuh"
#include "RasterizationBwd.cuh"

#include <map>
#include <stdexcept>
#include <vector>


static void _alloc_grad_buffers() {
    int64_t N = engine().max_num_splats;
    int64_t K = engine().num_sh;
    engine().grad.means.resize("eng.v_means", N);
    engine().grad.quats.resize("eng.v_quats", N);
    engine().grad.scales.resize("eng.v_scales", N);
    engine().grad.opacities.resize("eng.v_opacities", N);
    engine().grad.features_dc.resize("eng.v_features_dc", N);
    engine().grad.features_sh.resize("eng.v_features_sh", N, K);
    engine().grad.means.zero();
    engine().grad.quats.zero();
    engine().grad.scales.zero();
    engine().grad.opacities.zero();
    engine().grad.features_dc.zero();
    engine().grad.features_sh.zero();
}


std::map<std::string, float> engine_compute_loss_backward(
    int step,
    std::array<float, (int)LossWeightIndex::length> loss_weights,
    float w_ssim,
    int num_loss_scales,
    bool compute_loss_map
) {
    // Validate that forward was run
    if (std::get<0>(engine().fwd.renders) .data_ptr() == nullptr)
        throw std::runtime_error("engine_compute_loss_backward: forward_3dgs must be called first");
    if (!engine().gt.has_gt)
        throw std::runtime_error("engine_compute_loss_backward: set_training_data must be called first");

    // Allocate and zero gradient buffers from pool
    _alloc_grad_buffers();

    int64_t C = engine().camera.num;
    int64_t H = engine().camera.height;
    int64_t W = engine().camera.width;

    // Pool-allocate intermediates for loss computation
    TorchTensorView loss_map_buf = compute_loss_map ?
        _pool_tv_zero("eng.loss_map", C, H, W, 1) : _tv_null();

    // v_losses: constant vector [1, 0, 1, 1, ...] (1 for all, 0 for psnr slot)
    // Initialized once; pool never shrinks so pointer is stable
    static bool v_losses_initialized = false;
    TorchTensorView v_losses_buf = _pool_tv_1d("eng.v_losses", (int)LossIndex::length);
    if (!v_losses_initialized) {
        float h_v[(int)LossIndex::length];
        for (int i = 0; i < (int)LossIndex::length; i++) h_v[i] = 1.0f;
        h_v[(int)LossIndex::RgbPSNR] = 0.0f;
        cudaMemcpy((void*)std::get<0>(v_losses_buf), h_v, sizeof(h_v), cudaMemcpyHostToDevice);
        v_losses_initialized = true;
    }

    // Render outputs from forward pass (pool-backed, already populated)
    TorchTensorView render_rgb = TorchTensorView(
        (uint64_t)std::get<0>(engine().fwd.renders).data_ptr(), 4, {C, H, W, 3});
    TorchTensorView render_depth = TorchTensorView(
        (uint64_t)std::get<1>(engine().fwd.renders).data_ptr(), 4, {C, H, W, 1});
    TorchTensorView render_Ts = TorchTensorView(
        (uint64_t)engine().fwd.render_Ts.data_ptr(), 4, {C, H, W, 1});

    // Render normal: not yet available from forward (needs distortion output)
    // TODO: pass distortion outputs from forward when enabled
    TorchTensorView render_normal = _tv_null();
    TorchTensorView depth_normal = _tv_null();
    TorchTensorView rgb_dist = _tv_null();
    TorchTensorView depth_dist = _tv_null();
    TorchTensorView normal_dist = _tv_null();

    // Depth -> normal: derive depth_normal from rendered depth when gt_normal is provided
    // (matches training_losses.py logic: pred_normal is None, pred_depth exists, gt_normal exists).
    bool compute_depth_normal = (engine().gt.normal.data_ptr() != nullptr);
    bool is_ray_depth = (engine().primitive != "3dgs" && engine().primitive != "mip");
    if (compute_depth_normal) {
        depth_normal = _pool_tv("eng.depth_normal", C, H, W, 3);
        depth_to_normal_forward(
            engine().camera.model_str,
            _dv_tv(engine().camera.intrins),
            _dt2d_tv(engine().camera.dist_coeffs),
            is_ray_depth,
            DeviceTensor3D<float>(render_depth),
            DeviceTensor3D<float3>(depth_normal)
        );
    }

    std::vector<bool> needs_input_grad = {
        true,                                  // pred_rgb
        false,                                 // gt_rgb
        true,                                  // pred_depth
        engine().bilagrid_depth.enabled,       // gt_depth (true when bilagrid depth)
        true,                                  // pred_normal
        true,                                  // pred_depth_normal
        engine().bilagrid_normal.enabled,      // gt_normal (true when bilagrid normal)
        true,                                  // pred_transmittance
        true, true, true,                      // distortion
    };

    PerPixelGrads pixel_grads = {};

    // Allocate per-pixel gradient outputs
    pixel_grads.v_render_rgb  = _pool_tv("eng.v_rgb",   C, H, W, 3);
    pixel_grads.v_render_depth = _pool_tv("eng.v_depth", C, H, W, 1);
    pixel_grads.v_render_Ts   = _pool_tv("eng.v_Ts",    C, H, W, 1);
    if (compute_depth_normal) {
        pixel_grads.v_depth_normal = _pool_tv("eng.v_depth_normal", C, H, W, 3);
    }
    if (engine().bilagrid_depth.enabled) {
        pixel_grads.v_ref_depth = _pool_tv("eng.v_ref_depth", C, H, W, 1);
    }
    if (engine().bilagrid_normal.enabled) {
        pixel_grads.v_ref_normal = _pool_tv("eng.v_ref_normal", C, H, W, 3);
    }
    // TODO: allocate normal/distortion grads when those features are enabled

    // --- Compute per-pixel losses + SSIM, get gradients ---
    LossValues lv = compute_multi_scale_per_pixel_losses(
        num_loss_scales,
        render_rgb,
        _dt3d_tv(engine().gt.rgb),
        render_depth,
        _dt3d_tv(engine().gt.depth),
        render_normal,
        depth_normal,
        _dt3d_tv(engine().gt.normal),
        render_Ts,
        rgb_dist,
        depth_dist,
        normal_dist,
        _dt3d_tv(engine().gt.alpha),
        engine().gt.has_mask,
        loss_weights,
        w_ssim,
        v_losses_buf,
        needs_input_grad,
        -1,  // num_train_images: -1 means use batch size
        _tv_null(),  // camera_indices: null means identity
        loss_map_buf,
        pixel_grads
    );

    // --- PPISP backward hook ---
    // Forward order is render -> bilagrid -> PPISP -> loss, so PPISP backward
    // runs FIRST: it rewrites v_render_rgb (post-PPISP -> pre-PPISP) and
    // accumulates the per-camera PPISP parameter gradient. The bilagrid hook
    // below then consumes the pre-PPISP v_render_rgb.
    if (engine().ppisp.enabled) {
        _ensure_ppisp_optim_state();
        _engine_ppisp_backward_hook(pixel_grads.v_render_rgb);
    }

    // --- Bilagrid backward hook ---
    // Transforms v_render_rgb (post-bilagrid -> pre-bilagrid) and accumulates
    // gradients into bilagrid_*_grads for any enabled bilagrid types. The
    // updated v_render_rgb then flows into rasterization backward as usual.
    // For depth / normal (gt-side), the loss-side gradient is consumed here
    // and only the bilagrid grid gradient is retained.
    if (engine().bilagrid_rgb.enabled || engine().bilagrid_depth.enabled ||
        engine().bilagrid_normal.enabled) {
        _ensure_bilagrid_optim_state();
        _engine_bilagrid_backward_hook(
            pixel_grads.v_render_rgb,
            pixel_grads.v_ref_depth,
            pixel_grads.v_ref_normal);
    }

    // --- Background blend backward hook ---
    // Forward order is render -> background -> bilagrid -> PPISP -> loss, so
    // background backward runs LAST (after PPISP+bilagrid hooks). It rewrites
    // v_render_rgb (post-blend -> pre-blend) and ADDS the blend's
    // transmittance gradient into v_render_Ts before raster bwd consumes it.
    if (engine().background.enabled) {
        _engine_background_backward_hook(
            pixel_grads.v_render_rgb,
            pixel_grads.v_render_Ts);
    }

    // TODO: color space conversion (rgb_to_srgb) forward/backward

    // Depth -> normal backward: propagate v_depth_normal grads into v_render_depth (in-place add)
    if (compute_depth_normal) {
        depth_to_normal_backward(
            engine().camera.model_str,
            _dv_tv(engine().camera.intrins),
            _dt2d_tv(engine().camera.dist_coeffs),
            is_ray_depth,
            DeviceTensor3D<float>(render_depth),
            DeviceTensor3D<float3>(pixel_grads.v_depth_normal),
            DeviceTensor3D<float>(pixel_grads.v_render_depth)
        );
    }

    // --- Rasterization + projection backward ---
    // Build v_render_outputs from pixel_grads
    RenderOutput::TensorTuple v_render_outputs = std::make_tuple(
        DeviceTensor3D<float3>(pixel_grads.v_render_rgb),
        DeviceTensor3D<float>(pixel_grads.v_render_depth),
        DeviceTensor3D<float3>()  // no normal gradient yet
    );
    DeviceTensor3D<float> v_render_Ts(pixel_grads.v_render_Ts);

    DeviceTensorFloatND v_fnd_opac;
    {
        TorchTensorView opac_tv((uint64_t)engine().grad.opacities.data_ptr(), 4,
            {engine().grad.opacities.size(), 1LL, 1LL});
        v_fnd_opac = DeviceTensorFloatND(opac_tv);  // ndim=2, shape=[N, 1]
    }
    std::vector<DeviceTensorFloatND> v_splats_w = {
        DeviceTensorFloatND(engine().grad.means),         // [N, 3]
        DeviceTensorFloatND(engine().grad.quats),         // [N, 4]
        DeviceTensorFloatND(engine().grad.scales),        // [N, 3]
        v_fnd_opac,                                       // [N, 1]
        DeviceTensorFloatND(engine().grad.features_dc),   // [N, 3]
        DeviceTensorFloatND(engine().grad.features_sh),   // [N, K, 3]
    };

    // Build accum_weight_map from loss_map (pixel-space -> per-splat mapping in raster bwd)
    DeviceTensor3D<float> accum_weight_map;
    if (compute_loss_map && _tv_valid(loss_map_buf)) {
        accum_weight_map = DeviceTensor3D<float>(loss_map_buf);
    }

    std::vector<DeviceTensorFloatND> v_splats_w_out, v_splats_s_out;

    if (engine().primitive == "3dgs" || engine().primitive == "mip") {
        auto [vw, vs, accum_weight] = rasterize_to_pixels_3dgs_bwd(
            engine().cur_num_splats,
            engine().fwd.splats_w,
            engine().fwd.splats_s,
            engine().fwd.gaussian_ids,
            (uint32_t)engine().camera.width,
            (uint32_t)engine().camera.height,
            engine().fwd.tile_offsets,
            engine().fwd.flatten_ids,
            engine().fwd.render_Ts,
            engine().fwd.last_ids,
            engine().fwd.renders,
            accum_weight_map,
            v_render_outputs,
            v_render_Ts,
            std::make_optional(v_splats_w),
            std::nullopt
        );
        v_splats_w_out = vw;
        v_splats_s_out = vs;
        if (accum_weight.data_ptr()) engine().fwd.accum_weight = accum_weight;
    } else if (engine().primitive == "3dgut") {
        auto [vw, vs, vviewmats, accum_weight] = rasterize_to_pixels_3dgut_bwd(
            engine().cur_num_splats,
            engine().fwd.splats_w,
            engine().fwd.splats_s,
            engine().fwd.gaussian_ids,
            _dt2d_tv(engine().camera.viewmats),
            _dv_tv(engine().camera.intrins),
            engine().camera.model_str,
            _dt2d_tv(engine().camera.dist_coeffs),
            (uint32_t)engine().camera.width,
            (uint32_t)engine().camera.height,
            engine().fwd.tile_offsets,
            engine().fwd.flatten_ids,
            engine().fwd.render_Ts,
            engine().fwd.last_ids,
            engine().fwd.renders,
            std::nullopt,  // render2_outputs
            DeviceTensor3D<float>(),  // loss_map
            accum_weight_map,
            v_render_outputs,
            v_render_Ts,
            std::nullopt,  // v_distortion_outputs
            std::make_optional(v_splats_w),
            std::nullopt,
            false
        );
        v_splats_w_out = vw;
        v_splats_s_out = vs;
        if (accum_weight.data_ptr()) engine().fwd.accum_weight = accum_weight;
    } else {
        throw std::runtime_error("engine_compute_loss_backward: unknown primitive: " + engine().primitive);
    }

    // --- Projection backward ---
    if (engine().primitive == "3dgs") {
        projection_3dgs_backward(
            engine().cur_num_splats, engine().sh_degree, engine().fwd.splats_w,
            _dt2d_tv(engine().camera.viewmats), _dv_tv(engine().camera.intrins),
            (uint32_t)engine().camera.width, (uint32_t)engine().camera.height,
            engine().camera.model_str, _dt2d_tv(engine().camera.dist_coeffs),
            engine().fwd.camera_ids, engine().fwd.gaussian_ids,
            engine().fwd.aabb, v_splats_s_out, v_splats_w_out, nullptr);
    } else if (engine().primitive == "mip") {
        projection_mip_backward(
            engine().cur_num_splats, engine().sh_degree, engine().fwd.splats_w,
            _dt2d_tv(engine().camera.viewmats), _dv_tv(engine().camera.intrins),
            (uint32_t)engine().camera.width, (uint32_t)engine().camera.height,
            engine().camera.model_str, _dt2d_tv(engine().camera.dist_coeffs),
            engine().fwd.camera_ids, engine().fwd.gaussian_ids,
            engine().fwd.aabb, v_splats_s_out, v_splats_w_out, nullptr);
    } else if (engine().primitive == "3dgut") {
        projection_3dgut_backward(
            engine().cur_num_splats, engine().sh_degree, engine().fwd.splats_w,
            _dt2d_tv(engine().camera.viewmats), _dv_tv(engine().camera.intrins),
            (uint32_t)engine().camera.width, (uint32_t)engine().camera.height,
            engine().camera.model_str, _dt2d_tv(engine().camera.dist_coeffs),
            engine().fwd.camera_ids, engine().fwd.gaussian_ids,
            engine().fwd.aabb, v_splats_s_out, v_splats_w_out, nullptr);
    }

    // --- Build loss dict for display ---
    auto sdiv = [](float x, float y) -> float { return y != 0.0f ? x / y : 0.0f; };

    std::map<std::string, float> loss_dict;
    loss_dict["rgb_loss"] = lv.rgb_loss + w_ssim * (1.0f - lv.ssim);
    loss_dict["psnr"] = lv.rgb_psnr;
    loss_dict["ssim"] = lv.ssim;
    loss_dict["depth_loss"] = sdiv(lv.depth_sup,
        loss_weights[(int)LossWeightIndex::DepthSup]);
    loss_dict["normal_loss"] = sdiv(lv.normal_sup,
        loss_weights[(int)LossWeightIndex::NormalSup]);
    loss_dict["alpha_loss"] = sdiv(lv.alpha_sup,
        loss_weights[(int)LossWeightIndex::AlphaSup] + loss_weights[(int)LossWeightIndex::AlphaSupUnder]);
    loss_dict["normal_reg"] = sdiv(lv.normal_reg,
        loss_weights[(int)LossWeightIndex::NormalReg]);
    loss_dict["alpha_reg"] = sdiv(lv.alpha_reg,
        loss_weights[(int)LossWeightIndex::AlphaReg]);
    loss_dict["rgb_dist_reg"] = sdiv(lv.rgb_dist_reg,
        loss_weights[(int)LossWeightIndex::RgbDistReg]);
    loss_dict["depth_dist_reg"] = sdiv(lv.depth_dist_reg,
        loss_weights[(int)LossWeightIndex::DepthDistReg]);
    loss_dict["normal_dist_reg"] = sdiv(lv.normal_dist_reg,
        loss_weights[(int)LossWeightIndex::NormalDistReg]);

    return loss_dict;
}
