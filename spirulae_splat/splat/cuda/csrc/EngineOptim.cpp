// Engine Adam + per-splat regularization optimizer step.

#include "Engine.h"
#include "EngineCommon.h"
#include "EngineState.h"


static void _ensure_optim_state(bool quantize_sh, bool use_per_splat_bias_correction = false) {
    if (engine().optim.initialized && engine().optim.quantize_sh == quantize_sh
        && engine().optim.use_per_splat_bias_correction == use_per_splat_bias_correction)
        return;

    int64_t N = engine().max_num_splats;
    int64_t K = engine().num_sh;
    engine().optim.quantize_sh = quantize_sh;

    // g1 (exp_avg)
    engine().optim.g1_means.resize("eng.g1_means", N);
    engine().optim.g1_quats.resize("eng.g1_quats", N);
    engine().optim.g1_scales.resize("eng.g1_scales", N);
    engine().optim.g1_opacities.resize("eng.g1_opacities", N);
    engine().optim.g1_features_dc.resize("eng.g1_features_dc", N);
    if (!quantize_sh)
        engine().optim.g1_features_sh.resize("eng.g1_features_sh", N, K);

    // g2 (exp_avg_sq)
    engine().optim.g2_means.resize("eng.g2_means", N);
    engine().optim.g2_quats.resize("eng.g2_quats", N);
    engine().optim.g2_scales.resize("eng.g2_scales", N);
    engine().optim.g2_opacities.resize("eng.g2_opacities", N);
    engine().optim.g2_features_dc.resize("eng.g2_features_dc", N);
    if (!quantize_sh)
        engine().optim.g2_features_sh.resize("eng.g2_features_sh", N, K);

    // radii [max_N]
    engine().optim.radii.resize("eng.radii", N);

    // accum_buffer [max_N]
    engine().optim.accum_buffer.resize("eng.accum_buffer", N);

    // Joint (u, sqrt(g2)) Adam state for SH features. The engine optim path
    // uses Optimizer.cu's fused_adam_with_steps_8bit_kernel which launches
    // one thread per cell (= per (splat, coef, channel)), grouped 256 threads
    // per CUDA block. So one float4 bounds slot covers 256 contiguous cells.
    if (quantize_sh) {
        constexpr int64_t BLOCK_SIZE = QuantizedAdamState<8, 256>::kBlockSize;
        int64_t sh_cells = (int64_t)N * K * 3;
        int64_t sh_bounds = (sh_cells + BLOCK_SIZE - 1) / BLOCK_SIZE;
        engine().optim.sh_quant_state.resize("eng.sh_quant", sh_cells, sh_bounds);
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
    if (quantize_sh) {
        engine().optim.sh_quant_state.zero();
    } else {
        engine().optim.g1_features_sh.zero();
        engine().optim.g2_features_sh.zero();
    }
    engine().optim.accum_buffer.zero();
    engine().optim.bias_correction_steps.zero();

    engine().optim.initialized = true;
}


void engine_optim_step(int step, const OptimConfig& cfg) {
    _ensure_optim_state(cfg.quantize_sh, cfg.use_per_splat_bias_correction);

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

    fused_adam_step(N,
        DeviceTensorFloatND(engine().world.features_dc),
        DeviceTensorFloatND(engine().grad.features_dc),
        DeviceTensorFloatND(engine().optim.g1_features_dc),
        DeviceTensorFloatND(engine().optim.g2_features_dc),
        cfg.lr_features_dc, step + 1, per_splat_steps,
        cfg.sh_reg_weight, 0.5f / 0.28209479177387814f);

    if (cfg.quantize_sh && engine().optim.sh_quant_state.initialized()) {
        fused_adam_step_8bit(N,
            DeviceTensorFloatND(engine().world.features_sh),
            DeviceTensorFloatND(engine().grad.features_sh),
            engine().optim.sh_quant_state.packed_ptr(),
            engine().optim.sh_quant_state.bounds_ptr(),
            cfg.lr_features_sh, step + 1, per_splat_steps,
            cfg.sh_reg_weight, 0.0f);
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
