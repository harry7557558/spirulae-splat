// Engine linear / wide-gamut color space conversion.
//
// When the splat (model) works in a non-sRGB color space, the rendered RGB
// is converted to sRGB before any downstream stage (bilagrid / PPISP /
// loss). When the input GT imagery is in a non-sRGB space, the GT is
// converted once at upload time and kept in sRGB thereafter. The 3x3
// transform matrix maps the source color space to Rec.709; identity when
// source is sRGB already and is_linear is the only active flag.

#include "engine/Engine.h"
#include "engine/EngineCommon.h"
#include "engine/EngineInternal.h"
#include "engine/EngineState.h"

#include "kernels/pixelwise/PixelWise.cuh"

#include <stdexcept>
#include <vector>


static void _alloc_and_upload_matrix(
    DeviceTensor2D<float3>& dst, const float* host9, PoolSlot key
) {
    dst.resize(key, 3, 3);
    backend::memcpy_sync(dst.data_ptr(), host9, 9 * sizeof(float), backend::MemcpyKind::HostToDevice);
}


/*[AutoHeaderGeneratorExport]*/
void engine_init_color_space(
    // Splat side
    bool splat_enabled,
    bool splat_is_linear,
    std::vector<float> splat_color_matrix,   // [9], row-major
    // Image side
    bool image_enabled,
    bool image_is_linear,
    std::vector<float> image_color_matrix    // [9], row-major
) {
    auto& cs = engine().color_space;
    cs.splat_enabled   = splat_enabled;
    cs.splat_is_linear = splat_is_linear;
    if (splat_enabled) {
        if (splat_color_matrix.size() != 9)
            throw std::runtime_error("engine_init_color_space: splat_color_matrix must have 9 elements");
        _alloc_and_upload_matrix(cs.splat_color_matrix, splat_color_matrix.data(),
                                 PoolSlot::ColorSpaceSplatMatrix);
    }

    cs.image_enabled   = image_enabled;
    cs.image_is_linear = image_is_linear;
    if (image_enabled) {
        if (image_color_matrix.size() != 9)
            throw std::runtime_error("engine_init_color_space: image_color_matrix must have 9 elements");
        _alloc_and_upload_matrix(cs.image_color_matrix, image_color_matrix.data(),
                                 PoolSlot::ColorSpaceImageMatrix);
    }
}


void _engine_color_space_forward() {
    auto& cs = engine().color_space;
    if (!cs.splat_enabled) return;

    auto& fwd_rgb_tensor = std::get<0>(engine().fwd.renders);
    if (fwd_rgb_tensor.data_ptr() == nullptr) return;

    // Mirror the background-blend pattern: alias the current rgb buffer as
    // "pre" (no D2D copy), allocate a fresh "post" buffer, run the
    // conversion out-of-place, then re-point fwd.renders.rgb at post. The
    // pre buffer (held only by cs.fwd_pre and the previous fwd_rgb_tensor
    // pool allocation) survives until the next forward.
    cs.fwd_pre = fwd_rgb_tensor;

    int64_t B = cs.fwd_pre.size<0>(), H = cs.fwd_pre.size<1>(), W = cs.fwd_pre.size<2>();
    DeviceTensor3D<float3> post_rgb;
    post_rgb.resize(PoolSlot::ColorSpaceFwdPost, B, H, W);
    rgb_to_srgb_forward(cs.splat_is_linear, cs.fwd_pre, cs.splat_color_matrix,
                        post_rgb);
    fwd_rgb_tensor = post_rgb;
}


void _engine_color_space_backward_hook(TorchTensorView v_render_rgb) {
    auto& cs = engine().color_space;
    if (!cs.splat_enabled) return;
    if (cs.fwd_pre.data_ptr() == nullptr) return;

    DeviceTensor3D<float3> v_rgb_t = DeviceTensor3D<float3>(v_render_rgb);
    if (v_rgb_t.data_ptr() == nullptr) return;

    // v_rgb_t holds dL/d(rgb_sRGB) from the downstream chain. After the vjp
    // it holds dL/d(rgb_linear), which is what raster bwd expects. Per-pixel
    // local, in-place on v_rgb_t is safe. The background backward hook reads
    // bg.fwd_pre_blend_rgb (its own alias to the pre-bg rgb), not
    // engine().fwd.renders.rgb, so we don't need to swap pointers back.
    rgb_to_srgb_backward(cs.splat_is_linear, cs.fwd_pre, cs.splat_color_matrix,
                         v_rgb_t, v_rgb_t);
}


void _engine_color_space_apply_to_gt() {
    auto& cs = engine().color_space;
    if (!cs.image_enabled) return;

    DeviceTensor3D<float3> rgb = engine().gt.rgb;
    if (rgb.data_ptr() == nullptr) return;

    rgb_to_srgb_forward(cs.image_is_linear, rgb, cs.image_color_matrix, rgb);
}
