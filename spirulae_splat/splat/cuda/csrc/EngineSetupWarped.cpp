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
        uint16_t* p = DevicePool::global().acquire<uint16_t>("gt.staging_u16", numel);
        cudaMemcpy(p, (void*)src_ptr, numel * sizeof(uint16_t), cudaMemcpyHostToDevice);
        return (const void*)p;
    } else {
        uint8_t* p = DevicePool::global().acquire<uint8_t>("gt.staging_u8", numel);
        cudaMemcpy(p, (void*)src_ptr, numel * sizeof(uint8_t), cudaMemcpyHostToDevice);
        return (const void*)p;
    }
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
    // Depth / normal: warp path doesn't supply them, but if the engine still
    // has stale ones around from a prior step, drop them so the loss kernels
    // see "no GT depth / normal" (rather than mismatched-shape garbage).
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
        "gt.rgb",
        (size_t)B_post * out_H * out_W * 3);

    // Wide path uses input intrins / dist_coeffs; equirectangular path
    // doesn't (it ignores those args).
    const float* d_intrins = (float*)std::get<0>(input_intrins);
    const float* d_dist    = (float*)std::get<0>(input_dist_coeffs);

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
            "gt.alpha",
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

    engine().gt.has_gt = true;

    // Apply image-side color space conversion AFTER warp (in place on the
    // post-warp float buffer). No-op when no color space is configured.
    _engine_color_space_apply_to_gt();
}
