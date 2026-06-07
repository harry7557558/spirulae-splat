// Engine Adam + per-splat regularization optimizer step.

#include "Engine.h"
#include "EngineCommon.h"
#include "EngineState.h"

#include <stdexcept>
#include <variant>
#include <vector>


// FPBO block granularity for SH-quant bounds: one float4 per kFpboBlock splats,
// covering all 3*K cells per splat. Mirrors the BLOCK_SIZE used inside the
// fused projection-bwd+optim kernel launcher.
static constexpr int64_t kFpboBlock = 256;


static void _ensure_optim_state(int sh_optim_bits, bool use_per_splat_bias_correction = false) {
    if (sh_optim_bits != 32 && sh_optim_bits != 8 && sh_optim_bits != 4) {
        throw std::runtime_error(
            "_ensure_optim_state: sh_optim_bits must be 32 (off), 4, or 8; got " +
            std::to_string(sh_optim_bits));
    }
    bool quantize_sh = (sh_optim_bits != 32);
    bool fused = engine().optim.use_fused_proj_bwd_optim;
    if (engine().optim.initialized && engine().optim.sh_optim_bits == sh_optim_bits
        && engine().optim.use_per_splat_bias_correction == use_per_splat_bias_correction
        && engine().optim.fused_state_active == fused)
        return;

    int64_t N = engine().max_num_splats;
    int64_t K = engine().num_sh;
    engine().optim.sh_optim_bits = sh_optim_bits;
    engine().optim.fused_state_active = fused;

    // g1 (exp_avg)
    engine().optim.g1_means.resize("eng.g1_means", N);
    engine().optim.g1_quats.resize("eng.g1_quats", N);
    engine().optim.g1_scales.resize("eng.g1_scales", N);
    engine().optim.g1_opacities.resize("eng.g1_opacities", N);
    engine().optim.g1_features_dc.resize("eng.g1_features_dc", N);
    if (!quantize_sh)
        engine().optim.g1_features_sh.resize("eng.g1_features_sh", N, K);
    else
        engine().optim.g1_features_sh = DeviceTensor2D<float3>();

    // g2 (exp_avg_sq)
    engine().optim.g2_means.resize("eng.g2_means", N);
    engine().optim.g2_quats.resize("eng.g2_quats", N);
    engine().optim.g2_scales.resize("eng.g2_scales", N);
    engine().optim.g2_opacities.resize("eng.g2_opacities", N);
    engine().optim.g2_features_dc.resize("eng.g2_features_dc", N);
    if (!quantize_sh)
        engine().optim.g2_features_sh.resize("eng.g2_features_sh", N, K);
    else
        engine().optim.g2_features_sh = DeviceTensor2D<float3>();

    // radii [max_N]
    engine().optim.radii.resize("eng.radii", N);

    // accum_buffer [max_N]
    engine().optim.accum_buffer.resize("eng.accum_buffer", N);

    // Joint (u, sqrt(g2)) Adam state for SH features.
    //   Non-fused path: Optimizer.cu's fused_adam_with_steps_8bit_kernel
    //   launches one thread per cell (= per (splat, coef, channel)), grouped
    //   256 threads per CUDA block -- one float4 bounds slot per 256 cells.
    //   Fused path: FPBO launches one thread per SPLAT, grouped 256 splats
    //   per block, and the kernel writes one bounds slot per block. So same
    //   packed-cell count, but 3*K x fewer bounds slots.
    int64_t sh_cells = (int64_t)N * K * 3;
    if (quantize_sh && !fused) {
        constexpr int64_t BLOCK_SIZE = QuantizedAdamState<8, 256>::kBlockSize;
        int64_t sh_bounds = (sh_cells + BLOCK_SIZE - 1) / BLOCK_SIZE;
        engine().optim.sh_quant_state.resize("eng.sh_quant", sh_cells, sh_bounds);
        engine().optim.sh_quant_state_fpbo = QuantizedAdamState<8, 256>();
    } else if (quantize_sh && fused) {
        int64_t sh_bounds = (N + kFpboBlock - 1) / kFpboBlock;
        engine().optim.sh_quant_state_fpbo.resize("eng.sh_quant_fpbo", sh_cells, sh_bounds);
        engine().optim.sh_quant_state = QuantizedAdamState<8, 256>();
    } else {
        engine().optim.sh_quant_state      = QuantizedAdamState<8, 256>();
        engine().optim.sh_quant_state_fpbo = QuantizedAdamState<8, 256>();
    }

    // bias_correction_steps
    engine().optim.use_per_splat_bias_correction = use_per_splat_bias_correction;
    if (use_per_splat_bias_correction) {
        engine().optim.bias_correction_steps.resize("eng.bias_correction_steps", N);
    } else {
        engine().optim.bias_correction_steps = DeviceVector<int32_t>();
    }

    // Zero everything on first init
    engine().optim.g1_means.zero();       engine().optim.g2_means.zero();
    engine().optim.g1_quats.zero();       engine().optim.g2_quats.zero();
    engine().optim.g1_scales.zero();      engine().optim.g2_scales.zero();
    engine().optim.g1_opacities.zero();   engine().optim.g2_opacities.zero();
    engine().optim.g1_features_dc.zero(); engine().optim.g2_features_dc.zero();
    if (quantize_sh && !fused) {
        engine().optim.sh_quant_state.zero();
    } else if (quantize_sh && fused) {
        engine().optim.sh_quant_state_fpbo.zero();
    } else {
        engine().optim.g1_features_sh.zero();
        engine().optim.g2_features_sh.zero();
    }
    engine().optim.accum_buffer.zero();
    engine().optim.bias_correction_steps.zero();

    engine().optim.initialized = true;
}


// Build a [N, P]-shaped DeviceTensorFloatND view over a raw float buffer used
// for the per-attribute Adam state (g1/g2). Mirrors the layout produced by
// `DeviceTensorFloatND(DeviceVector<floatK>)`.
static inline DeviceTensorFloatND _fnd_view(float* ptr, int64_t N, int64_t P) {
    TorchTensorView tv((uint64_t)ptr, 4, {N, P, 1LL});
    return DeviceTensorFloatND(tv);
}


// Fused projection-bwd + Adam-and-regularization optimizer step. Replaces
// the projection_*_backward + engine_optim_step pair when
// `cfg.use_fused_proj_bwd_optim` is set. Requires
// engine_compute_loss_backward to have stashed v_splats_s into
// engine().fwd.v_splats_s. World-space gradient buffers are not allocated
// in this path (raster_*_bwd's world atomicStore is null-pointer-guarded
// per channel, and the fused kernel's atomicLoad uses the same guard).
void engine_fused_proj_bwd_optim_step(int step, const OptimConfig& cfg) {
    _ensure_optim_state(cfg.sh_optim_bits,
                        cfg.use_per_splat_bias_correction);

    int64_t N = engine().cur_num_splats;
    if (engine().fwd.v_splats_s.empty())
        throw std::runtime_error(
            "engine_fused_proj_bwd_optim_step: v_splats_s not stashed; "
            "engine_compute_loss_backward must run first");

    // World-grad vector handed to the FPBO kernel. For 3dgs / mip the
    // per-channel slots stay null (raster_*_bwd never writes the world
    // buffer for those prims, so there is nothing to pick up). For 3dgut,
    // raster_*_bwd has already atomicAdded mean / quat / scale into the
    // engine grad buffers (_alloc_grad_buffers allocates them in that case);
    // pass those three slots so the kernel's atomicLoad picks them up and
    // sums them with the projection-bwd contribution. Other slots remain
    // null -- they were never written and reading 0 is correct.
    std::vector<DeviceTensorFloatND> v_splats_w;
    if (engine().grad.means.data_ptr() != nullptr) {
        v_splats_w = {
            DeviceTensorFloatND(engine().grad.means),         // [N, 3]
            DeviceTensorFloatND(engine().grad.quats),         // [N, 4]
            DeviceTensorFloatND(engine().grad.scales),        // [N, 3]
            DeviceTensorFloatND(),                            // opacities -- null
            DeviceTensorFloatND(),                            // features_dc -- null
            DeviceTensorFloatND(),                            // features_sh -- null
        };
    }

    // Adam state. For SH, fp32 g1/g2 OR quantized packed-bytes -- never both.
    int64_t M  = engine().max_num_splats;
    int64_t K3 = (int64_t)engine().num_sh * 3;

    // SH g1/g2 view: zero-shaped empty view when quantized SH is active so
    // the kernel takes its quant path (driven by sh_packed != nullptr).
    const bool quantize_sh = (cfg.sh_optim_bits != 32);
    DeviceTensorFloatND g1_sh = quantize_sh
        ? DeviceTensorFloatND()
        : _fnd_view((float*)engine().optim.g1_features_sh.data_ptr(), M, K3);
    DeviceTensorFloatND g2_sh = quantize_sh
        ? DeviceTensorFloatND()
        : _fnd_view((float*)engine().optim.g2_features_sh.data_ptr(), M, K3);

    std::vector<DeviceTensorFloatND> g1 = {
        _fnd_view((float*)engine().optim.g1_means.data_ptr(),       M, 3),
        _fnd_view((float*)engine().optim.g1_quats.data_ptr(),       M, 4),
        _fnd_view((float*)engine().optim.g1_scales.data_ptr(),      M, 3),
        _fnd_view((float*)engine().optim.g1_opacities.data_ptr(),   M, 1),
        _fnd_view((float*)engine().optim.g1_features_dc.data_ptr(), M, 3),
        g1_sh,
    };
    std::vector<DeviceTensorFloatND> g2 = {
        _fnd_view((float*)engine().optim.g2_means.data_ptr(),       M, 3),
        _fnd_view((float*)engine().optim.g2_quats.data_ptr(),       M, 4),
        _fnd_view((float*)engine().optim.g2_scales.data_ptr(),      M, 3),
        _fnd_view((float*)engine().optim.g2_opacities.data_ptr(),   M, 1),
        _fnd_view((float*)engine().optim.g2_features_dc.data_ptr(), M, 3),
        g2_sh,
    };

    // SH quantization: pass packed-byte + per-block bounds when active.
    std::optional<TorchTensorView> sh_packed_opt   = std::nullopt;
    std::optional<TorchTensorView> sh_bounds_opt   = std::nullopt;
    if (quantize_sh && engine().optim.sh_quant_state_fpbo.initialized()) {
        auto& qs = engine().optim.sh_quant_state_fpbo;
        sh_packed_opt = TorchTensorView(
            (uint64_t)qs.packed_ptr(), 1, {qs.packed_bytes()});
        sh_bounds_opt = TorchTensorView(
            (uint64_t)qs.bounds_ptr(), 4, {qs.n_bounds, 4LL});
    }

    // Per-splat bias correction: increment all steps by 1 before use.
    std::variant<int32_t, TorchTensorView> step_arg;
    if (cfg.use_per_splat_bias_correction) {
        increment_int32_inplace(engine().optim.bias_correction_steps,
                                engine().max_num_splats);
        step_arg = _dv_tv(engine().optim.bias_correction_steps);
    } else {
        step_arg = (int32_t)(step + 1);
    }

    // World-splat tensors come straight from fwd cache; they alias the
    // engine's mutable world params, which the fused kernel updates in place.
    auto& splats_w = engine().fwd.splats_w;

    auto aabb_nd = DeviceTensorFloatND(engine().fwd.aabb);

    // Dispatch on primitive. Pass engine().sh_degree (= sh_degree_to_use) so
    // the kernel-template SH degree matches the forward pass during warmup --
    // otherwise SH backward computes spurious gradients on unused bands which
    // contaminate v_mean through the viewdir chain.
    int max_sh_degree = engine().sh_degree;
    auto call_dispatch = [&](auto fn) {
        fn(
            N, max_sh_degree, splats_w,
            _dt2d_tv(engine().camera.viewmats),
            _dv_tv(engine().camera.intrins),
            (uint32_t)engine().camera.width,
            (uint32_t)engine().camera.height,
            engine().camera.model_str,
            _dt2d_tv(engine().camera.dist_coeffs),
            engine().fwd.camera_ids,
            engine().fwd.gaussian_ids,
            aabb_nd,
            v_splats_w, std::nullopt, std::nullopt,
            engine().fwd.v_splats_s, std::nullopt, std::nullopt,
            g1, g2,
            sh_packed_opt, sh_bounds_opt,
            engine().optim.radii,
            cfg.lr_means, cfg.lr_quats, cfg.lr_scales,
            cfg.lr_opacities, cfg.lr_features_dc, cfg.lr_features_sh,
            cfg.max_gauss_ratio, cfg.scale_regularization_weight,
            cfg.mcmc_opacity_reg_weight, cfg.mcmc_scale_reg_weight,
            cfg.erank_reg_weight, cfg.erank_reg_weight_s3,
            cfg.quat_norm_reg_weight, cfg.sh_reg_weight,
            cfg.use_scale_agnostic_mean,
            cfg.use_color_trust_region,
            cfg.color_is_linear,
            cfg.eps_tr,
            step_arg,
            cfg.sh_optim_bits
        );
    };

    if (engine().primitive == "3dgs") {
        call_dispatch(&fused_projection_bwd_optimizer_3dgs);
    } else if (engine().primitive == "mip") {
        call_dispatch(&fused_projection_bwd_optimizer_mip);
    } else if (engine().primitive == "3dgut") {
        call_dispatch(&fused_projection_bwd_optimizer_3dgut);
    } else {
        throw std::runtime_error(
            "engine_fused_proj_bwd_optim_step: unsupported primitive: " +
            engine().primitive);
    }
}


void engine_optim_step(int step, const OptimConfig& cfg) {
    if (engine().optim.use_fused_proj_bwd_optim) {
        engine_fused_proj_bwd_optim_step(step, cfg);
        return;
    }

    _ensure_optim_state(cfg.sh_optim_bits, cfg.use_per_splat_bias_correction);

    int64_t N = engine().cur_num_splats;

    // Per-splat bias correction: increment all steps by 1 before use
    DeviceVector<int32_t> per_splat_steps;
    if (cfg.use_per_splat_bias_correction) {
        increment_int32_inplace(engine().optim.bias_correction_steps, engine().max_num_splats);
        per_splat_steps = engine().optim.bias_correction_steps;
    }

    fused_optim_3dgs_geometry(
        N,
        engine().world.means,     engine().grad.means,     engine().optim.g1_means,     engine().optim.g2_means,
        engine().world.quats,     engine().grad.quats,     engine().optim.g1_quats,     engine().optim.g2_quats,
        engine().world.scales,    engine().grad.scales,    engine().optim.g1_scales,    engine().optim.g2_scales,
        engine().world.opacities, engine().grad.opacities, engine().optim.g1_opacities, engine().optim.g2_opacities,
        engine().optim.radii,
        cfg.lr_means, cfg.lr_quats, cfg.lr_scales, cfg.lr_opacities,
        cfg.max_gauss_ratio, cfg.scale_regularization_weight,
        cfg.mcmc_opacity_reg_weight, cfg.mcmc_scale_reg_weight,
        cfg.erank_reg_weight, cfg.erank_reg_weight_s3, cfg.quat_norm_reg_weight,
        cfg.use_scale_agnostic_mean,
        step + 1, per_splat_steps
    );

    // DC color: trust-region Adam when the splat works in a non-sRGB color
    // space, otherwise plain fused Adam. The TR variants live in
    // Optimizer.cu (fused_adamtr_(linear_)rgb_optim) and clip each step to
    // +/-kSh0*sqrt(4*eps_tr*c/opac) so the working-color-space update stays
    // inside the model's confidence radius.
    if (cfg.use_color_trust_region) {
        const int s1 = step + 1;
        const float beta1 = 0.9f, beta2 = 0.999f, eps = 1e-15f;
        (cfg.color_is_linear ? fused_adamtr_linear_rgb_optim : fused_adamtr_rgb_optim)(
            _dv_tv(engine().world.features_dc),
            _dv_tv(engine().grad.features_dc),
            _dv_tv(engine().optim.g1_features_dc),
            _dv_tv(engine().optim.g2_features_dc),
            _dv_tv(engine().world.opacities),
            cfg.lr_features_dc,
            beta1, beta2, eps, cfg.eps_tr, s1
        );
    } else {
        fused_adam_step(N,
            DeviceTensorFloatND(engine().world.features_dc),
            DeviceTensorFloatND(engine().grad.features_dc),
            DeviceTensorFloatND(engine().optim.g1_features_dc),
            DeviceTensorFloatND(engine().optim.g2_features_dc),
            cfg.lr_features_dc, step + 1, per_splat_steps,
            cfg.sh_reg_weight, 0.5f / 0.28209479177387814f);
    }

    if (cfg.sh_optim_bits != 32 && engine().optim.sh_quant_state.initialized()) {
        // SH-quant + TR is not supported by fused_adamtr_*_sh_optim (those
        // kernels need fp32 g1/g2). Fall back to the plain quantized kernel.
        // Document this as a limitation for callers; FPBO doesn't have this
        // restriction because the FPBO kernel handles both internally.
        fused_adam_step_quantized(N,
            DeviceTensorFloatND(engine().world.features_sh),
            DeviceTensorFloatND(engine().grad.features_sh),
            engine().optim.sh_quant_state.packed_ptr(),
            engine().optim.sh_quant_state.bounds_ptr(),
            cfg.lr_features_sh, step + 1, per_splat_steps,
            cfg.sh_reg_weight, 0.0f,
            cfg.sh_optim_bits);
    } else if (cfg.use_color_trust_region) {
        const int s1 = step + 1;
        const float beta1 = 0.9f, beta2 = 0.999f, eps = 1e-15f;
        (cfg.color_is_linear ? fused_adamtr_linear_rgb_sh_optim : fused_adamtr_rgb_sh_optim)(
            _dt2d_tv(engine().world.features_sh),
            _dt2d_tv(engine().grad.features_sh),
            _dt2d_tv(engine().optim.g1_features_sh),
            _dt2d_tv(engine().optim.g2_features_sh),
            _dv_tv(engine().world.features_dc),
            _dv_tv(engine().world.opacities),
            cfg.lr_features_sh,
            beta1, beta2, eps, cfg.eps_tr, s1
        );
    } else {
        fused_adam_step(N,
            DeviceTensorFloatND(engine().world.features_sh),
            DeviceTensorFloatND(engine().grad.features_sh),
            DeviceTensorFloatND(engine().optim.g1_features_sh),
            DeviceTensorFloatND(engine().optim.g2_features_sh),
            cfg.lr_features_sh, step + 1, per_splat_steps,
            cfg.sh_reg_weight, 0.0f);
    }
}
