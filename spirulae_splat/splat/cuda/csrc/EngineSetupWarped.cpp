// Warp-fused GT upload path.
//
// Used by engine_train_step_managed when a batch's input camera model is
// FISHEYE / EQUISOLID (with warp_to_pinhole) or EQUIRECTANGULAR (always
// split). The byte-staged GT image is warped DIRECTLY into a float
// post-split output buffer; no intermediate full-resolution float image
// is ever materialized.
//
// Depth / normal GT are unsupported on this path (the warp kernels operate
// only on RGB + binary mask). PPISP IS supported: it runs render-side on
// the post-split pinhole sub-cameras, and N_grids must be sized to the
// post-split camera count (one parameter slot per sub-camera).

#include "Engine.h"
#include "EngineCommon.h"
#include "EngineState.h"
#include "Common.cuh"
#include "PixelWise.cuh"

#include <stdexcept>
#include <string>


// Stage a byte buffer to device. Reuses the same shared staging slot as
// _hv_to_dt3d_gt (key "gt.staging_u8" / "gt.staging_u16") so peak host->device
// memory stays the same in both paths.
static const void* _h2d_stage_byte(
    const TorchTensorView& src_tv, size_t numel, uint32_t elem_size,
    bool& out_is_u16)
{
    if (std::get<0>(src_tv) == 0) return nullptr;
    if (elem_size != 1 && elem_size != 2)
        throw std::runtime_error("set_training_data_warped: byte stage requires uint8 or uint16");
    out_is_u16 = (elem_size == 2);
    uint64_t src_ptr = std::get<0>(src_tv);
    if (_is_device_ptr((void*)src_ptr)) {
        return (const void*)src_ptr;
    }
    if (out_is_u16) {
        uint16_t* p = DevicePool::global().acquire<uint16_t>(PoolSlot::GtStagingU16, numel);
        backend::memcpy_sync(p, (void*)src_ptr, numel * sizeof(uint16_t), backend::MemcpyKind::HostToDevice);
        return (const void*)p;
    } else {
        uint8_t* p = DevicePool::global().acquire<uint8_t>(PoolSlot::GtStagingU8, numel);
        backend::memcpy_sync(p, (void*)src_ptr, numel * sizeof(uint8_t), backend::MemcpyKind::HostToDevice);
        return (const void*)p;
    }
}


// Stage a small float array (per-input intrins / dist_coeffs) to device. The
// wide warp kernels dereference these on-device as intrins[bid] /
// dist_coeffs.load(bid); when the managed C++ DataManager drives this path the
// views point at host std::vector memory (DecodedBatch hands out host
// TorchTensorViews), so they must be copied to device first -- otherwise the
// kernel faults on an illegal address. Zero-copy when already a device ptr.
static const float* _h2d_stage_floats(
    const TorchTensorView& src_tv, size_t numel, PoolSlot slot)
{
    uint64_t src_ptr = std::get<0>(src_tv);
    if (src_ptr == 0 || numel == 0) return nullptr;
    if (_is_device_ptr((void*)src_ptr)) {
        return (const float*)src_ptr;
    }
    float* p = DevicePool::global().acquire<float>(slot, numel);
    backend::memcpy_sync(p, (void*)src_ptr, numel * sizeof(float), backend::MemcpyKind::HostToDevice);
    return p;
}


void set_training_data_warped(
    // Input-side camera model name (FISHEYE / EQUISOLID / EQUIRECTANGULAR).
    // The kernel dispatcher uses this to pick wide vs equirectangular warp.
    std::string input_model_name,
    int B_in, int in_H, int in_W,
    int K, int out_H, int out_W,
    // GT RGB at INPUT shape, byte (uint8 or uint16). The kernel fuses
    // byte->float + warp and writes the result directly into engine().gt.rgb
    // at [B_in*K, out_H, out_W, 3] float.
    TorchTensorView gt_rgb,
    // Optional mask at (potentially different) input shape, uint8.
    TorchTensorView gt_alpha,
    int mask_in_H, int mask_in_W,
    // Optional GT depth (uint16 raw counts / float) and normal (uint8 /
    // float) at their own input shapes. Depth is warped to per-face ray
    // depth; normal is rotated into each face's camera frame.
    TorchTensorView gt_depth,
    int depth_in_H, int depth_in_W,
    TorchTensorView gt_normal,
    int normal_in_H, int normal_in_W,
    bool input_depth_is_ray_depth,
    // Per-INPUT-camera intrins + dist_coeffs (used by the wide warp; ignored
    // by the equirectangular kernel). Both are TorchTensorView so they can
    // be zero-copy when already on device.
    TorchTensorView input_intrins,
    TorchTensorView input_dist_coeffs,
    // Device pointer to the [K, 3, 3] cubemap axes table (lifetime managed
    // by DataManager).
    uint64_t axes_dev)
{
    // PPISP runs render-side on the POST-split pinhole sub-cameras (same
    // shape as the bilagrid RGB path), so it works unchanged here -- no
    // refusal needed.
    // Depth / normal are warped below; clear stale buffers first so a step
    // that supplies neither leaves the loss kernels seeing "no GT".
    engine().gt.depth  = DeviceTensor3D<float>();
    engine().gt.normal = DeviceTensor3D<float3>();

    if (std::get<0>(gt_rgb) == 0) {
        engine().gt.has_gt = false;
        engine().gt.has_mask = false;
        return;
    }

    // ---- RGB byte H2D staging ---------------------------------------------
    uint32_t rgb_elem = std::get<1>(gt_rgb);
    if (rgb_elem != 1 && rgb_elem != 2) {
        throw std::runtime_error(
            "set_training_data_warped: only uint8 / uint16 RGB are supported "
            "on the warp path");
    }
    bool rgb_u16 = false;
    size_t rgb_numel = (size_t)B_in * in_H * in_W * 3;
    const void* d_rgb_byte = _h2d_stage_byte(gt_rgb, rgb_numel, rgb_elem, rgb_u16);

    // ---- Allocate POST-split float output buffer + run warp ---------------
    int B_post = B_in * K;
    float* d_rgb_float = DevicePool::global().acquire<float>(
        PoolSlot::GtRgb,
        (size_t)B_post * out_H * out_W * 3);

    // Wide path uses input intrins / dist_coeffs; equirectangular path
    // doesn't (it ignores those args). Stage to device when host-resident --
    // the managed DataManager passes host views, and the wide warp kernels
    // dereference these pointers on-device.
    const float* d_intrins = _h2d_stage_floats(
        input_intrins, (size_t)B_in * 4, PoolSlot::WarpInputIntrins);
    const float* d_dist    = _h2d_stage_floats(
        input_dist_coeffs, (size_t)B_in * 10, PoolSlot::WarpInputDistCoeffs);

    CameraModelType cm = cmt(input_model_name);
    if (cm == CameraModelType::EQUIRECTANGULAR) {
        launch_warp_byte_to_float_equi(
            d_rgb_byte, rgb_u16, B_in, in_H, in_W, 3,
            d_rgb_float, K, out_H, out_W,
            (const float*)axes_dev);
    } else {
        // FISHEYE / EQUISOLID / PINHOLE -- the wide warp kernel does the
        // projection dispatch internally on `cm`. (PINHOLE doesn't normally
        // hit this path -- K==1 there -- but the kernel handles it anyway.)
        launch_warp_byte_to_float_wide(
            input_model_name, d_intrins, d_dist,
            d_rgb_byte, rgb_u16, B_in, in_H, in_W, 3,
            d_rgb_float, K, out_H, out_W,
            (const float*)axes_dev);
    }

    // Wrap as DeviceTensor3D<float3> [B_post, out_H, out_W, 1] (one float3
    // = 3 floats per pixel, matching the engine's existing layout).
    {
        TorchTensorView dv((uint64_t)d_rgb_float, 4,
                           {(int64_t)B_post, (int64_t)out_H, (int64_t)out_W, 3LL});
        engine().gt.rgb = DeviceTensor3D<float3>(dv);
    }

    // ---- Mask byte warp ---------------------------------------------------
    if (std::get<0>(gt_alpha) != 0) {
        size_t mask_numel = (size_t)B_in * mask_in_H * mask_in_W;
        bool   m_u16 = false;
        const void* d_mask_in = _h2d_stage_byte(gt_alpha, mask_numel,
                                                std::get<1>(gt_alpha), m_u16);
        if (m_u16) {
            throw std::runtime_error("set_training_data_warped: mask must be uint8");
        }
        uint8_t* d_mask_out = DevicePool::global().acquire<uint8_t>(
            PoolSlot::GtAlpha,
            (size_t)B_post * out_H * out_W);
        if (cm == CameraModelType::EQUIRECTANGULAR) {
            launch_warp_mask_equi(
                (const uint8_t*)d_mask_in, B_in, mask_in_H, mask_in_W,
                d_mask_out, K, out_H, out_W, (const float*)axes_dev);
        } else {
            launch_warp_mask_wide(
                input_model_name, d_intrins, d_dist,
                (const uint8_t*)d_mask_in, B_in, mask_in_H, mask_in_W,
                d_mask_out, K, out_H, out_W, (const float*)axes_dev);
        }
        TorchTensorView dv((uint64_t)d_mask_out, 1,
                           {(int64_t)B_post, (int64_t)out_H, (int64_t)out_W, 1LL});
        engine().gt.alpha = DeviceTensor3D<bool>(dv);
        engine().gt.has_mask = true;
    } else {
        engine().gt.alpha = DeviceTensor3D<bool>();
        engine().gt.has_mask = false;
    }

    // ---- Depth warp -------------------------------------------------------
    // Wide GT depth -> per-face ray depth float [B_post, out_H, out_W, 1].
    // Already ray depth on output, so the loss must NOT re-convert (the
    // caller sets the engine's ray-depth handling accordingly).
    if (std::get<0>(gt_depth) != 0) {
        uint32_t d_elem = std::get<1>(gt_depth);
        size_t d_numel  = (size_t)B_in * depth_in_H * depth_in_W;
        const void* d_depth_in;
        if (d_elem == 4) {
            d_depth_in = (const void*)_h2d_stage_floats(
                gt_depth, d_numel, PoolSlot::GtStagingU16);
        } else {
            bool d_u16 = false;
            d_depth_in = _h2d_stage_byte(gt_depth, d_numel, d_elem, d_u16);
        }
        float* d_depth_out = DevicePool::global().acquire<float>(
            PoolSlot::GtDepth, (size_t)B_post * out_H * out_W);
        if (cm == CameraModelType::EQUIRECTANGULAR) {
            launch_warp_depth_equi(
                d_depth_in, d_elem, B_in, depth_in_H, depth_in_W,
                d_depth_out, K, out_H, out_W,
                (const float*)axes_dev, input_depth_is_ray_depth);
        } else {
            launch_warp_depth_wide(
                input_model_name, d_intrins, d_dist,
                d_depth_in, d_elem, B_in, depth_in_H, depth_in_W,
                in_H, in_W, d_depth_out, K, out_H, out_W,
                (const float*)axes_dev, input_depth_is_ray_depth);
        }
        TorchTensorView dv((uint64_t)d_depth_out, 4,
                           {(int64_t)B_post, (int64_t)out_H, (int64_t)out_W, 1LL});
        engine().gt.depth = DeviceTensor3D<float>(dv);
    }

    // ---- Normal warp ------------------------------------------------------
    // Wide GT normal (INPUT camera frame) -> per-face camera-frame normal
    // float3 [B_post, out_H, out_W, 3], rotated by the face axes.
    if (std::get<0>(gt_normal) != 0) {
        uint32_t n_elem = std::get<1>(gt_normal);
        size_t n_numel  = (size_t)B_in * normal_in_H * normal_in_W * 3;
        const void* d_normal_in;
        if (n_elem == 4) {
            d_normal_in = (const void*)_h2d_stage_floats(
                gt_normal, n_numel, PoolSlot::GtStagingU8);
        } else {
            bool n_u16 = false;
            d_normal_in = _h2d_stage_byte(gt_normal, n_numel, n_elem, n_u16);
            if (n_u16)
                throw std::runtime_error(
                    "set_training_data_warped: normal must be uint8 or float32");
        }
        float* d_normal_out = DevicePool::global().acquire<float>(
            PoolSlot::GtNormal, (size_t)B_post * out_H * out_W * 3);
        if (cm == CameraModelType::EQUIRECTANGULAR) {
            launch_warp_normal_equi(
                d_normal_in, n_elem, B_in, normal_in_H, normal_in_W,
                d_normal_out, K, out_H, out_W, (const float*)axes_dev);
        } else {
            launch_warp_normal_wide(
                input_model_name, d_intrins, d_dist,
                d_normal_in, n_elem, B_in, normal_in_H, normal_in_W,
                in_H, in_W, d_normal_out, K, out_H, out_W,
                (const float*)axes_dev);
        }
        TorchTensorView dv((uint64_t)d_normal_out, 4,
                           {(int64_t)B_post, (int64_t)out_H, (int64_t)out_W, 3LL});
        engine().gt.normal = DeviceTensor3D<float3>(dv);
    }

    engine().gt.has_gt = true;

    // Apply image-side color space conversion AFTER warp (in place on the
    // post-warp float buffer). No-op when no color space is configured.
    _engine_color_space_apply_to_gt();
}
