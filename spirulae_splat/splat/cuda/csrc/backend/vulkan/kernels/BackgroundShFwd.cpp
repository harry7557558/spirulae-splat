// Vulkan implementation of the background spherical-harmonics launch API
// (csrc/BackgroundSphericalHarmonics.cuh). Mirrors the host dispatcher in
// BackgroundSphericalHarmonics.cu; the device work runs
// slang/vulkan/background_sh.slang with the SH degree as a specialization
// constant. The backward is training-phase work and throws for now.

#include <BackgroundSphericalHarmonics.cuh>
#include <Common.cuh>

#include "../VulkanInternal.h"
#include "../VulkanPipelines.h"

#include <stdexcept>

namespace {

// Mirrors BackgroundShFwdParams in slang/vulkan/background_sh.slang.
struct BackgroundShFwdParams {
    uint64_t viewmats, intrins, dist_coeffs, sh_coeffs, out_img;
    uint32_t width, height, B;
    int32_t camera_model;
};
static_assert(sizeof(BackgroundShFwdParams) == 5 * 8 + 4 * 4,
              "params layout must match the slang struct");

int64_t _batch_count(const TorchTensorView& t, int64_t per_item) {
    int64_t n = 1;
    for (auto s : std::get<2>(t)) n *= s;
    return n / per_item;
}

}  // namespace

/* API definitions matching csrc/BackgroundSphericalHarmonics.cuh */

void render_background_sh_forward(
    int w,
    int h,
    std::string camera_model,
    int sh_degree,
    TorchTensorView viewmats,
    TorchTensorView intrins,
    TorchTensorView dist_coeffs,
    TorchTensorView sh_coeffs,
    TorchTensorView out_color
) {
    if (sh_degree < 0 || sh_degree > 4)
        throw std::runtime_error(
            "render_background_sh_forward: sh_degree must be in [0, 4]");

    CameraModelType cam = cmt(camera_model);
    if ((int)cam < 0 || (int)cam > 3)
        throw std::runtime_error("Camera model " + camera_model +
                                 " is not supported for skybox");

    int64_t b = _batch_count(out_color, (int64_t)h * w * 3);
    if (b * h * w == 0) return;
    if (_batch_count(viewmats, 16) != b || _batch_count(intrins, 4) != b)
        throw std::runtime_error(
            "viewmats/intrins must have batch dim matching out_color");

    BackgroundShFwdParams p{};
    p.viewmats = std::get<0>(viewmats);
    p.intrins = std::get<0>(intrins);
    p.dist_coeffs = std::get<0>(dist_coeffs);
    p.sh_coeffs = std::get<0>(sh_coeffs);
    p.out_img = std::get<0>(out_color);
    p.width = (uint32_t)w;
    p.height = (uint32_t)h;
    p.B = (uint32_t)b;
    p.camera_model = (int32_t)cam;

    backend::vk::SpecList spec{(uint32_t)sh_degree};
    if (!backend::vk::dispatch(backend::kDefaultStream,
                               "background_sh.background_sh_fwd", spec,
                               (uint32_t)((w + 127) / 128), (uint32_t)h,
                               (uint32_t)b, &p, sizeof(p)))
        throw std::runtime_error(
            "Vulkan backend: background SH dispatch failed");
}

void render_background_sh_backward(
    int, int, std::string, int,
    TorchTensorView, TorchTensorView, TorchTensorView, TorchTensorView,
    TorchTensorView, TorchTensorView, TorchTensorView
) {
    throw std::runtime_error(
        "Vulkan backend: render_background_sh_backward is not implemented "
        "yet (training phase)");
}
