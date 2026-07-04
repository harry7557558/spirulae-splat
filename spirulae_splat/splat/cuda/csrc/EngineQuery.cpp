// Read-only / cheap query entrypoints (D->H copies for inspection).

#include "Engine.h"
#include "EngineCommon.h"
#include "EngineState.h"

#include <algorithm>


void engine_copy_accum_buffer(TorchTensorView dst) {
    if (engine().optim.accum_buffer.data_ptr() == nullptr || std::get<0>(dst) == 0)
        return;
    int64_t dst_n = std::get<2>(dst)[0];
    int64_t src_n = engine().optim.accum_buffer.size();
    int64_t n = std::min(dst_n, src_n);
    cudaMemcpy((void*)std::get<0>(dst),
               engine().optim.accum_buffer.data_ptr(),
               n * sizeof(float2), cudaMemcpyDeviceToHost);
}


int64_t engine_get_cur_num_splats() {
    return engine().cur_num_splats;
}


void engine_copy_render_to_host(
    TorchTensorView out_rgb,
    TorchTensorView out_depth,
    TorchTensorView out_Ts,
    TorchTensorView out_rgb_raw,
    TorchTensorView out_median
) {
    auto& renders = engine().fwd.renders;
    auto& rgb = std::get<0>(renders);
    auto& depth = std::get<1>(renders);
    if (rgb.data_ptr() && std::get<0>(out_rgb) != 0) {
        cudaMemcpy((void*)std::get<0>(out_rgb), rgb.data_ptr(),
                   rgb.numel() * sizeof(float3), cudaMemcpyDeviceToHost);
    }
    if (depth.data_ptr() && std::get<0>(out_depth) != 0) {
        cudaMemcpy((void*)std::get<0>(out_depth), depth.data_ptr(),
                   depth.numel() * sizeof(float), cudaMemcpyDeviceToHost);
    }
    if (engine().fwd.render_Ts.data_ptr() && std::get<0>(out_Ts) != 0) {
        cudaMemcpy((void*)std::get<0>(out_Ts), engine().fwd.render_Ts.data_ptr(),
                   engine().fwd.render_Ts.numel() * sizeof(float), cudaMemcpyDeviceToHost);
    }
    // Pre-conversion (linear / wide-gamut) render is stashed by the color
    // space forward hook into cs.fwd_pre. When the engine has no color
    // space configured, fwd_pre is empty; the caller mirrors out_rgb to
    // out_rgb_raw on the Python side to avoid a redundant D->H of the
    // identical buffer.
    auto& cs = engine().color_space;
    if (std::get<0>(out_rgb_raw) != 0 && cs.fwd_pre.data_ptr() != nullptr) {
        cudaMemcpy((void*)std::get<0>(out_rgb_raw), cs.fwd_pre.data_ptr(),
                   cs.fwd_pre.numel() * sizeof(float3), cudaMemcpyDeviceToHost);
    }
    if (engine().fwd.render_median.data_ptr() && std::get<0>(out_median) != 0) {
        cudaMemcpy((void*)std::get<0>(out_median), engine().fwd.render_median.data_ptr(),
                   engine().fwd.render_median.numel() * sizeof(float), cudaMemcpyDeviceToHost);
    }
}


// Per-pixel distortion image D = W*S - C^2 from the most recent forward
// (populated only when output_distortion was requested). RGB_D primitives
// emit rgb + depth only; normal is left empty. Mirrors
// engine_copy_render_to_host. Null out_* args are skipped.
void engine_copy_distortion_to_host(
    TorchTensorView out_rgb_dist,    // [C, H, W, 3] float32, CPU, optional
    TorchTensorView out_depth_dist   // [C, H, W, 1] float32, CPU, optional
) {
    auto& dist = engine().fwd.distortions;
    auto& rgb = std::get<0>(dist);
    auto& depth = std::get<1>(dist);
    if (rgb.data_ptr() && std::get<0>(out_rgb_dist) != 0) {
        cudaMemcpy((void*)std::get<0>(out_rgb_dist), rgb.data_ptr(),
                   rgb.numel() * sizeof(float3), cudaMemcpyDeviceToHost);
    }
    if (depth.data_ptr() && std::get<0>(out_depth_dist) != 0) {
        cudaMemcpy((void*)std::get<0>(out_depth_dist), depth.data_ptr(),
                   depth.numel() * sizeof(float), cudaMemcpyDeviceToHost);
    }
}


void engine_copy_splats_to_host(
    TorchTensorView means,
    TorchTensorView quats,
    TorchTensorView scales,
    TorchTensorView opacities,
    TorchTensorView features_dc,
    TorchTensorView features_sh
) {
    _dv_to_host(engine().world.means, means);
    _dv_to_host(engine().world.quats, quats);
    _dv_to_host(engine().world.scales, scales);
    _dv_to_host(engine().world.opacities, opacities);
    _dv_to_host(engine().world.features_dc, features_dc);
    // features_sh: DeviceTensor2D<float3>
    if (engine().world.features_sh.data_ptr() && std::get<0>(features_sh) != 0) {
        cudaMemcpy((void*)std::get<0>(features_sh), engine().world.features_sh.data_ptr(),
                   engine().world.features_sh.numel() * sizeof(float3), cudaMemcpyDeviceToHost);
    }
}


// Per-splat gradients accumulated by the most recent engine backward
// (engine_compute_loss_backward or engine_backward_from_render_grad). Buffers
// are sized max_num_splats and zeroed at the start of each backward. Mirrors
// engine_copy_splats_to_host but reads engine().grad.* instead of world.*.
void engine_copy_grads_to_host(
    TorchTensorView means,
    TorchTensorView quats,
    TorchTensorView scales,
    TorchTensorView opacities,
    TorchTensorView features_dc,
    TorchTensorView features_sh
) {
    _dv_to_host(engine().grad.means, means);
    _dv_to_host(engine().grad.quats, quats);
    _dv_to_host(engine().grad.scales, scales);
    _dv_to_host(engine().grad.opacities, opacities);
    _dv_to_host(engine().grad.features_dc, features_dc);
    // features_sh: DeviceTensor2D<float3>
    if (engine().grad.features_sh.data_ptr() && std::get<0>(features_sh) != 0) {
        cudaMemcpy((void*)std::get<0>(features_sh), engine().grad.features_sh.data_ptr(),
                   engine().grad.features_sh.numel() * sizeof(float3), cudaMemcpyDeviceToHost);
    }
}


// ---------------------------------------------------------------------------
// Debug introspection helpers. These pull the current engine training-data
// + rendered images back to host so a Python utility can save / display them
// for visual debugging (especially useful for the warp-to-pinhole path
// where the on-engine GT is generated on GPU and never seen by Python).
// ---------------------------------------------------------------------------

// Shape getters: callers use these to size the host buffer before calling
// the copy. Returns (B, H, W, C); zeros when the buffer is empty.
std::tuple<int64_t, int64_t, int64_t, int64_t> engine_get_gt_rgb_shape() {
    auto& t = engine().gt.rgb;
    if (t.data_ptr() == nullptr) return {0, 0, 0, 0};
    return {t.template size<0>(), t.template size<1>(), t.template size<2>(), 3LL};
}

std::tuple<int64_t, int64_t, int64_t, int64_t> engine_get_gt_alpha_shape() {
    auto& t = engine().gt.alpha;
    if (t.data_ptr() == nullptr) return {0, 0, 0, 0};
    return {t.template size<0>(), t.template size<1>(), t.template size<2>(), 1LL};
}

std::tuple<int64_t, int64_t, int64_t, int64_t> engine_get_render_rgb_shape() {
    auto& t = std::get<0>(engine().fwd.renders);
    if (t.data_ptr() == nullptr) return {0, 0, 0, 0};
    return {t.template size<0>(), t.template size<1>(), t.template size<2>(), 3LL};
}

void engine_copy_gt_rgb_to_host(TorchTensorView out) {
    auto& t = engine().gt.rgb;
    if (t.data_ptr() == nullptr || std::get<0>(out) == 0) return;
    cudaMemcpy((void*)std::get<0>(out), t.data_ptr(),
               t.numel() * sizeof(float3), cudaMemcpyDeviceToHost);
}

void engine_copy_gt_alpha_to_host(TorchTensorView out) {
    auto& t = engine().gt.alpha;
    if (t.data_ptr() == nullptr || std::get<0>(out) == 0) return;
    cudaMemcpy((void*)std::get<0>(out), t.data_ptr(),
               t.numel() * sizeof(bool), cudaMemcpyDeviceToHost);
}


std::vector<std::tuple<std::string, size_t, size_t>> engine_get_pool_breakdown() {
    return DevicePool::global().getBreakdown();
}

// Same buffers as engine_get_pool_breakdown(), but each row also carries its
// VRAM category as a string ("splat", "splat x img", "image", "appearance",
// "viewer", "other") sourced from the pool-slot metadata table -- so Python
// buckets by an authoritative tag instead of guessing from the key prefix.
std::vector<std::tuple<std::string, std::string, size_t, size_t>>
engine_get_pool_breakdown_categorized() {
    auto raw = DevicePool::global().getBreakdownCategorized();
    std::vector<std::tuple<std::string, std::string, size_t, size_t>> out;
    out.reserve(raw.size());
    for (auto& r : raw)
        out.emplace_back(std::get<0>(r), to_string((VramCategory)std::get<1>(r)),
                         std::get<2>(r), std::get<3>(r));
    return out;
}

size_t engine_get_scratch_bytes() {
    return DeviceScratch::global().capBytes();
}
