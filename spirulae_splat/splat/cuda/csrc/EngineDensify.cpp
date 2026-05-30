// Engine densification step (MCMC + revised relocate/add splits).

#include "Engine.h"
#include "EngineCommon.h"
#include "EngineState.h"

#include <algorithm>
#include <cmath>


// Helper: build DeviceVector<T> from TorchTensorView (non-owning view).
template<typename T>
static inline DeviceVector<T> _dv(const TorchTensorView& tv) {
    return DeviceVector<T>(tv);
}
// For null-able TorchTensorViews, return default (null) DeviceVector.
template<typename T>
static inline DeviceVector<T> _dv_or_null(const TorchTensorView& tv) {
    if (std::get<0>(tv) == 0) return DeviceVector<T>();
    return DeviceVector<T>(tv);
}
// For [N, K, C] tensors like features_sh: flatten to [N*K, C] then build DeviceVector.
template<typename T>
static inline DeviceVector<T> _dv_flat(const TorchTensorView& tv) {
    if (std::get<0>(tv) == 0) return DeviceVector<T>();
    const auto& shape = std::get<2>(tv);
    int64_t flat_count = 1;
    for (size_t i = 0; i + 1 < shape.size(); i++) flat_count *= shape[i];
    TorchTensorView flat_tv(std::get<0>(tv), std::get<1>(tv),
        {flat_count, (int64_t)(sizeof(T) / std::get<1>(tv))});
    return DeviceVector<T>(flat_tv);
}


int engine_densify_step(int step, int max_steps, const DensifyConfig& cfg) {
    // _ensure_optim_state(engine().optim.quantize_sh);

    int64_t cur_num_splats = engine().cur_num_splats;
    int64_t max_num_splats = engine().max_num_splats;
    bool quantize_sh = engine().optim.quantize_sh;
    int num_sh = engine().num_sh;

    bool use_revised = (cfg.relocate_heuristic_weight >= 1.0f);
    bool densify_ongoing = (step < max_steps - cfg.refine_stop_num_iter);
    bool do_densify = densify_ongoing && (step > cfg.refine_start_iter && step % cfg.refine_every == 0);
    float progress = ((float)step + 0.5f) / (float)max_steps;

    // Use pool-backed DeviceVector/DeviceTensor from Buffers directly
    auto& dv_means = engine().world.means;
    auto& dv_quats = engine().world.quats;
    auto& dv_scales = engine().world.scales;
    auto& dv_opacs = engine().world.opacities;
    auto& dv_features_dc = engine().world.features_dc;
    // features_sh: DeviceTensor2D<float3> [N, K] -> flatten to DeviceVector<float3> [N*K]
    DeviceVector<float3> dv_features_sh;
    {
        auto& t = engine().world.features_sh;
        TorchTensorView tv((uint64_t)t.data_ptr(), 4,
            {t.size<0>() * t.size<1>(), 3LL});
        dv_features_sh = DeviceVector<float3>(tv);
    }
    auto& dv_g1_means = engine().optim.g1_means;
    auto& dv_g1_quats = engine().optim.g1_quats;
    auto& dv_g1_scales = engine().optim.g1_scales;
    auto& dv_g1_opacs = engine().optim.g1_opacities;
    auto& dv_g1_features_dc = engine().optim.g1_features_dc;
    auto& dv_g2_means = engine().optim.g2_means;
    auto& dv_g2_quats = engine().optim.g2_quats;
    auto& dv_g2_scales = engine().optim.g2_scales;
    auto& dv_g2_opacs = engine().optim.g2_opacities;
    auto& dv_g2_features_dc = engine().optim.g2_features_dc;
    // g1/g2 features_sh as DeviceVector<float3> with numel=N*K. The densify kernel
    // reinterprets the data as uchar3 when quantize_sh, so the element_size carrier
    // type is just for the constructor's shape check.
    DeviceVector<float3> dv_g1_features_sh, dv_g2_features_sh;
    int64_t sh_flat = (int64_t)engine().max_num_splats * (int64_t)engine().num_sh;
    if (quantize_sh) {
        // Joint (u, sqrt(g2)) AoS layout: bytes for both halves of every cell
        // sit in a single packed buffer. Densify only zeros bytes on dst
        // splats, so we pass the same packed pointer for both the g1 and g2
        // wrappers; the kernel template is instantiated with `short3` (6 bytes
        // per SH coef = 3 channels x 2 bytes), so each write zeros the joint
        // (u, sqrt_g2) record for one (splat, coef). The second write
        // (g2_features_sh[i] = 0) targets the same memory and is a no-op.
        uint8_t* packed = engine().optim.sh_quant_state.packed_ptr();
        TorchTensorView tv((uint64_t)packed, 1, {sh_flat, 12LL});
        dv_g1_features_sh = DeviceVector<float3>(tv);
        dv_g2_features_sh = DeviceVector<float3>(tv);
    } else {
        // Non-quantized: storage is float3 N*K elements (12 bytes each), shape [N*K, 3], element_size=4
        TorchTensorView tv1((uint64_t)engine().optim.g1_features_sh.data_ptr(), 4, {sh_flat, 3LL});
        TorchTensorView tv2((uint64_t)engine().optim.g2_features_sh.data_ptr(), 4, {sh_flat, 3LL});
        dv_g1_features_sh = DeviceVector<float3>(tv1);
        dv_g2_features_sh = DeviceVector<float3>(tv2);
    }
    auto& dv_radii = engine().optim.radii;
    auto& dv_accum_buf = engine().optim.accum_buffer;
    auto& dv_bias_steps = engine().optim.bias_correction_steps;

    // Clip large splats
    if (std::isfinite(cfg.max_screen_size) || std::isfinite(cfg.max_world_size)) {
        densify_clip_scale_tensor(
            cur_num_splats,
            dv_radii, dv_scales,
            nullptr,  // don't clip opacities
            cfg.max_screen_size,
            cfg.max_screen_size_clip_hardness,
            cfg.max_world_size
        );
    }

    // Update densification score using accum_weight from rasterization backward
    if (densify_ongoing && use_revised && dv_accum_buf.data_ptr() != nullptr) {
        DeviceVector<float> score;
        if (engine().fwd.accum_weight.data_ptr() != nullptr) {
            score = engine().fwd.accum_weight;
        } else {
            score = engine().grad.opacities.data_ptr() ? engine().grad.opacities : dv_opacs;
        }

        densify_update_weight_tensor(
            cur_num_splats,
            dv_radii,
            nullptr,
            (float*)dv_opacs.data_ptr(),
            score,
            dv_accum_buf,
            false
        );
    }

    int num_added = 0;

    if (do_densify && use_revised) {
        // Revised relocation (long axis split)
        relocate_splats_with_long_axis_split_tensor(
            cur_num_splats, cfg.min_opacity,
            dv_means, dv_quats, dv_scales, dv_opacs, dv_features_dc, dv_features_sh,
            dv_g1_means, dv_g1_quats, dv_g1_scales, dv_g1_opacs, dv_g1_features_dc, dv_g1_features_sh,
            dv_g2_means, dv_g2_quats, dv_g2_scales, dv_g2_opacs, dv_g2_features_dc, dv_g2_features_sh,
            dv_accum_buf, dv_bias_steps,
            quantize_sh, num_sh, 2 * step + 0
        );

        // Add more splats
        int64_t n_target = std::min(max_num_splats, (int64_t)(cfg.growth_factor * cur_num_splats));
        num_added = (int)std::max((int64_t)0, n_target - cur_num_splats);
        if (num_added > 0) {
            add_splats_with_long_axis_split_tensor(
                cur_num_splats, num_added,
                dv_means, dv_quats, dv_scales, dv_opacs, dv_features_dc, dv_features_sh,
                dv_g1_means, dv_g1_quats, dv_g1_scales, dv_g1_opacs, dv_g1_features_dc, dv_g1_features_sh,
                dv_g2_means, dv_g2_quats, dv_g2_scales, dv_g2_opacs, dv_g2_features_dc, dv_g2_features_sh,
                dv_accum_buf, dv_bias_steps,
                quantize_sh, num_sh, 2 * step + 1
            );
        }
    } else if (do_densify) {
        // MCMC relocation
        relocate_splats_mcmc_tensor(
            cur_num_splats, cfg.min_opacity,
            dv_means, dv_quats, dv_scales, dv_opacs, dv_features_dc, dv_features_sh,
            dv_g1_means, dv_g1_quats, dv_g1_scales, dv_g1_opacs, dv_g1_features_dc, dv_g1_features_sh,
            dv_g2_means, dv_g2_quats, dv_g2_scales, dv_g2_opacs, dv_g2_features_dc, dv_g2_features_sh,
            dv_bias_steps,
            quantize_sh, num_sh, 2 * step + 0
        );

        // MCMC sample add
        int64_t n_target = std::min(max_num_splats, (int64_t)(cfg.growth_factor * cur_num_splats));
        num_added = (int)std::max((int64_t)0, n_target - cur_num_splats);
        if (num_added > 0) {
            add_splats_mcmc_tensor(
                cur_num_splats, num_added, cfg.min_opacity,
                dv_means, dv_quats, dv_scales, dv_opacs, dv_features_dc, dv_features_sh,
                dv_g1_means, dv_g1_quats, dv_g1_scales, dv_g1_opacs, dv_g1_features_dc, dv_g1_features_sh,
                dv_g2_means, dv_g2_quats, dv_g2_scales, dv_g2_opacs, dv_g2_features_dc, dv_g2_features_sh,
                dv_bias_steps,
                quantize_sh, num_sh, 2 * step + 1
            );
        }
    }

    // Reset accum buffer after densification step
    if (do_densify && dv_accum_buf.data_ptr() != nullptr) {
        cudaMemset(dv_accum_buf.data_ptr(), 0, dv_accum_buf.size() * sizeof(float2));
    }

    // Add MCMC noise
    if (cfg.noise_lr > 0.0f && cfg.noise_lr_final > 0.0f) {
        float noise_scalar = cfg.noise_lr * powf(cfg.noise_lr_final / cfg.noise_lr, progress);
        if (use_revised) {
            revised_add_noise_tensor(
                cur_num_splats + num_added, noise_scalar,
                dv_radii, dv_means, dv_scales, dv_quats, dv_opacs);
        } else {
            mcmc_add_noise_tensor(
                cur_num_splats + num_added, noise_scalar,
                dv_means, dv_scales, dv_quats, dv_opacs);
        }
    }

    engine().cur_num_splats = cur_num_splats + num_added;
    return num_added;
}
