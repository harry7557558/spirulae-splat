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
    TorchTensorView out_Ts
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


std::vector<std::tuple<std::string, size_t, size_t>> engine_get_pool_breakdown() {
    return DevicePool::global().getBreakdown();
}

size_t engine_get_scratch_bytes() {
    return DeviceScratch::global().capBytes();
}
