// Engine setup entrypoints: copy splat/camera/GT inputs to device-resident pool.

#include "Engine.h"
#include "EngineCommon.h"
#include "EngineInternal.h"
#include "EngineState.h"

#include <stdexcept>


void set_data_3dgs(
    int64_t num_splats,
    TorchTensorView means,
    TorchTensorView quats,
    TorchTensorView scales,
    TorchTensorView opacities,
    TorchTensorView features_dc,
    TorchTensorView features_sh
) {
    int64_t max_num_splats = std::get<2>(means)[0];
    if (std::get<2>(quats)[0] != max_num_splats ||
        std::get<2>(scales)[0] != max_num_splats ||
        std::get<2>(opacities)[0] != max_num_splats ||
        std::get<2>(features_dc)[0] != max_num_splats ||
        std::get<2>(features_sh)[0] != max_num_splats)
        throw std::runtime_error("setData3DGS: max_num_splats mismatch across splat tensors");
    if (max_num_splats < num_splats)
        throw std::runtime_error("setData3DGS: tensor size smaller than num_splats");

    engine().cur_num_splats = num_splats;
    engine().max_num_splats = max_num_splats;

    auto sh_shape = std::get<2>(features_sh);
    engine().num_sh = (sh_shape.size() >= 2) ? (int)sh_shape[1] : 0;

    // Splats are set once and persist on device; subsequent calls are no-ops.
    if (engine().world.initialized) return;

    engine().world.means       = _hv_to_dv<float3>("world.means", means);
    engine().world.quats       = _hv_to_dv<float4>("world.quats", quats);
    engine().world.scales      = _hv_to_dv<float3>("world.scales", scales);
    engine().world.opacities   = _hv_to_dv<float>("world.opacities", opacities);
    engine().world.features_dc = _hv_to_dv<float3>("world.features_dc", features_dc);
    engine().world.features_sh = _hv_to_dt2d<float3>("world.features_sh", features_sh);

    engine().world.initialized = true;
}


void set_camera_params(
    int width,
    int height,
    std::string camera_model,
    TorchTensorView viewmats,
    TorchTensorView intrins,
    TorchTensorView dist_coeffs
) {
    engine().camera.width = width;
    engine().camera.height = height;
    engine().camera.model = cmt(camera_model);
    engine().camera.model_str = camera_model;
    engine().camera.num = std::get<2>(viewmats)[0];

    if (std::get<2>(intrins)[0] != engine().camera.num ||
        std::get<2>(dist_coeffs)[0] != engine().camera.num)
        throw std::runtime_error("setCameraParams: num_cameras mismatch");

    // viewmats: PyTorch shape [C, 4, 4] -> treat as DeviceTensor2D<float4> [C, 4]
    int64_t C = std::get<2>(viewmats)[0];
    {
        uint64_t src_ptr = std::get<0>(viewmats);
        float4* ptr;
        if (_is_device_ptr((void*)src_ptr)) {
            ptr = (float4*)src_ptr;
        } else {
            ptr = DevicePool::global().acquire<float4>("cam.viewmats", (size_t)(C * 4));
            cudaMemcpy(ptr, (void*)src_ptr, C * 4 * sizeof(float4), cudaMemcpyHostToDevice);
        }
        TorchTensorView dv_tv((uint64_t)ptr, 4, {C, 4LL, 4LL});
        engine().camera.viewmats = DeviceTensor2D<float4>(dv_tv);
    }
    engine().camera.intrins     = _hv_to_dv<float4>("cam.intrins", intrins);
    engine().camera.dist_coeffs = _hv_to_dt2d<float>("cam.dist_coeffs", dist_coeffs);
}


void set_training_data(
    TorchTensorView gt_rgb,
    TorchTensorView gt_depth,
    TorchTensorView gt_normal,
    TorchTensorView gt_alpha,
    bool input_depth_is_ray_depth
) {
    engine().gt.rgb    = _hv_to_dt3d_gt<float3>(gt_rgb,    "gt.rgb",    "rgb");
    engine().gt.depth  = _hv_to_dt3d_gt<float>(gt_depth,   "gt.depth",  "depth");
    engine().gt.normal = _hv_to_dt3d_gt<float3>(gt_normal, "gt.normal", "normal");
    // gt_alpha: small bool/uint8 buffer (the external mask). No conversion;
    // the slang kernel reads bool per pixel. Drives engine().gt.has_mask.
    engine().gt.alpha  = _hv_to_dt3d<bool>("gt.alpha", gt_alpha);
    engine().gt.has_gt    = (std::get<0>(gt_rgb) != 0);
    engine().gt.has_mask  = (std::get<0>(gt_alpha) != 0);

    // Linear (z) depth -> ray depth, in place on the GPU buffer, before the
    // depth bilateral grid / loss consume it. The rasterizer renders ray
    // depth, so linear supervision depth must be converted to match. Uses the
    // current camera intrinsics (set by set_camera_params just before this),
    // rescaled inside the kernel to the depth map resolution.
    if (!input_depth_is_ray_depth && engine().gt.depth.data_ptr() != nullptr) {
        linear_depth_to_ray_depth_inplace(
            engine().camera.model_str,
            _dv_tv(engine().camera.intrins),
            _dt2d_tv(engine().camera.dist_coeffs),
            engine().camera.width, engine().camera.height,
            engine().gt.depth
        );
    }

    // Apply image-side color space conversion (linear / wide-gamut -> sRGB)
    // in place on the freshly uploaded GT. No-op when no image color space
    // was configured.
    _engine_color_space_apply_to_gt();
}
