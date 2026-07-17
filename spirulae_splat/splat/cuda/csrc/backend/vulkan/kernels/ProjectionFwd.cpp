// Vulkan implementation of the projection-forward launch API
// (csrc/ProjectionFwd.cuh). Mirrors ProjectionFwd.cu's launcher logic; the
// device work runs slang/vulkan/projection_fwd.slang with camera model /
// SH degree / antialiased as specialization constants.

#include <ProjectionFwd.cuh>

#include "../VulkanInternal.h"
#include "../VulkanPipelines.h"

#include <cstring>
#include <stdexcept>

namespace {

// Mirrors ProjectionFwdParams in slang/vulkan/projection_fwd.slang.
struct ProjectionFwdParams {
    uint64_t means, quats, scales, opacities, features_dc, features_sh;
    uint64_t viewmats, intrins, dist_coeffs;
    uint64_t out_aabb, out_depths, out_radii;
    uint64_t s_xy, s_depth, s_conic, s_opac, s_rgb;
    uint64_t sh_packed, sh_bounds;
    int64_t sh_bounds_stride;
    uint32_t sh_stride;
    uint32_t C, N;
    uint32_t width, height;
    uint32_t wgs_per_row;
    uint32_t num_sh_buffer;
    uint32_t _pad0;
};
static_assert(sizeof(ProjectionFwdParams) == 20 * 8 + 8 * 4,
              "params layout must match the slang struct");

// Mirrors Projection3dgutParams in slang/vulkan/projection_fwd.slang.
struct Projection3dgutParams {
    uint64_t means, quats, scales, opacities, features_dc, features_sh;
    uint64_t viewmats, intrins, dist_coeffs;
    uint64_t out_aabb, out_depths, out_radii;
    uint64_t s_scale, s_opac, s_rgb;
    uint64_t sh_packed, sh_bounds;
    int64_t sh_bounds_stride;
    uint32_t sh_stride;
    uint32_t C, N;
    uint32_t width, height;
    uint32_t wgs_per_row;
    uint32_t num_sh_buffer;
    uint32_t _pad0;
};
static_assert(sizeof(Projection3dgutParams) == 18 * 8 + 8 * 4,
              "params layout must match the slang struct");

// Shared validation + resolution of the SH value-quant launch args. Returns
// the kShValueBits spec value (0 = fp32). Mirrors the CUDA kernel's stride
// default: 0 -> FPBO per-splat-block layout (256 * 3 * num_sh_buffer).
uint32_t resolve_sh_quant(
    const std::optional<TorchTensorView>& sh_value_packed,
    const std::optional<TorchTensorView>& sh_value_bounds,
    const uint32_t num_sh_buffer, const int sh_value_bits,
    const int64_t sh_bounds_stride,
    uint64_t* out_packed, uint64_t* out_bounds, int64_t* out_stride) {
    *out_packed = 0;
    *out_bounds = 0;
    *out_stride = 1;  // never 0 (the shader divides by it)
    if (sh_value_bits == 32)
        return 0;
    if (sh_value_bits != 8 && sh_value_bits != 16)
        throw std::runtime_error(
            "projection forward: sh_value_bits must be 8, 16 or 32");
    if (!sh_value_packed.has_value() || !sh_value_bounds.has_value())
        throw std::runtime_error(
            "projection forward: sh_value_bits != 32 requires "
            "sh_value_packed and sh_value_bounds");
    *out_packed = std::get<0>(sh_value_packed.value());
    *out_bounds = std::get<0>(sh_value_bounds.value());
    *out_stride = sh_bounds_stride > 0
                      ? sh_bounds_stride
                      : (int64_t)256 * 3 * (int64_t)num_sh_buffer;
    return (uint32_t)sh_value_bits;
}

std::tuple<DeviceTensor2D<float4>, DeviceTensor2D<float>,
           std::vector<DeviceTensorFloatND>>
launch_projection_fwd_vk(
    const bool antialiased,
    const int64_t N,
    const std::vector<DeviceTensorFloatND>& in_splats,
    const int max_sh_degree,
    TorchTensorView viewmats, TorchTensorView intrins,
    const uint32_t image_width, const uint32_t image_height,
    const std::string& camera_model, const TorchTensorView& dist_coeffs,
    DeviceVector<float>& radii,
    const std::optional<TorchTensorView>& sh_value_packed,
    const std::optional<TorchTensorView>& sh_value_bounds,
    const uint32_t num_sh_buffer,
    const int sh_value_bits,
    const int64_t sh_bounds_stride) {
    uint64_t q_packed, q_bounds;
    int64_t q_stride;
    const uint32_t spec_bits =
        resolve_sh_quant(sh_value_packed, sh_value_bounds, num_sh_buffer,
                         sh_value_bits, sh_bounds_stride, &q_packed,
                         &q_bounds, &q_stride);

    CameraModelType cam = cmt(camera_model);
    if ((int)cam < 0 || (int)cam > 3)
        throw std::runtime_error("Unsupported camera model");

    Vanilla3DGS<0>::WorldBuffer wb(in_splats);
    int sh_degree = wb.sh_degree();
    sh_degree = std::min(sh_degree, max_sh_degree);

    const uint32_t C = (uint32_t)std::get<2>(viewmats)[0];

    DeviceTensor2D<float4> aabb;
    aabb.resize(PoolSlot::ProjAabb, C, N);
    DeviceTensor2D<float> sorting_depths;
    sorting_depths.resize(PoolSlot::ProjDepths, C, N);
    std::vector<DeviceTensorFloatND> splats_screen =
        Vanilla3DGS<0>::ScreenBuffer::empty_pool(C * N,
                                                 PoolSlot::ProjScreen);
    Vanilla3DGS<0>::ScreenBuffer sb(splats_screen);

    const uint64_t total_threads = (uint64_t)C * (uint64_t)N;
    const uint32_t wgs = (uint32_t)((total_threads + 127) / 128);
    const uint32_t per_row = std::min(wgs, 65535u);
    const uint32_t rows = (wgs + per_row - 1) / per_row;

    ProjectionFwdParams p{};
    p.means = (uint64_t)wb.raw_data(0);
    p.quats = (uint64_t)wb.raw_data(1);
    p.scales = (uint64_t)wb.raw_data(2);
    p.opacities = (uint64_t)wb.raw_data(3);
    p.features_dc = (uint64_t)wb.raw_data(4);
    p.features_sh = (uint64_t)wb.raw_data(5);
    p.sh_stride = (uint32_t)wb.raw_stride(5);
    p.viewmats = std::get<0>(viewmats);
    p.intrins = std::get<0>(intrins);
    p.dist_coeffs = std::get<0>(dist_coeffs);
    p.out_aabb = (uint64_t)aabb.data_ptr();
    p.out_depths = (uint64_t)sorting_depths.data_ptr();
    p.out_radii = (uint64_t)radii.data_ptr();
    p.s_xy = (uint64_t)sb.raw_data(0);
    p.s_depth = (uint64_t)sb.raw_data(1);
    p.s_conic = (uint64_t)sb.raw_data(2);
    p.s_opac = (uint64_t)sb.raw_data(3);
    p.s_rgb = (uint64_t)sb.raw_data(4);
    p.sh_packed = q_packed;
    p.sh_bounds = q_bounds;
    p.sh_bounds_stride = q_stride;
    p.num_sh_buffer = num_sh_buffer;
    p.C = C;
    p.N = (uint32_t)N;
    p.width = image_width;
    p.height = image_height;
    p.wgs_per_row = per_row;

    uint64_t params_addr = 0;
    void* params_mapped = nullptr;
    if (!backend::vk::params_alloc(sizeof(p), &params_addr, &params_mapped))
        throw std::runtime_error("Vulkan backend: params ring failed");
    std::memcpy(params_mapped, &p, sizeof(p));

    backend::vk::SpecList spec{(uint32_t)cam, (uint32_t)sh_degree,
                               antialiased ? 1u : 0u, spec_bits};
    if (!backend::vk::dispatch(backend::kDefaultStream,
                               "projection_fwd.projection_fwd_3dgs", spec,
                               per_row, rows, 1, &params_addr,
                               sizeof(params_addr)))
        throw std::runtime_error("Vulkan backend: projection dispatch failed");

    return std::make_tuple(aabb, sorting_depths, splats_screen);
}

}  // namespace

/* API definitions matching csrc/ProjectionFwd.cuh */

std::tuple<
    DeviceTensor2D<float4>, DeviceTensor2D<float>, std::vector<DeviceTensorFloatND>
> projection_3dgs_forward(
    const int64_t num_splats, const int max_sh_degree,
    const std::vector<DeviceTensorFloatND> &in_splats,
    TorchTensorView viewmats, TorchTensorView intrins,
    const uint32_t image_width, const uint32_t image_height,
    const std::string camera_model, const TorchTensorView dist_coeffs,
    DeviceVector<float> radii,
    const std::optional<TorchTensorView> sh_value_packed,
    const std::optional<TorchTensorView> sh_value_bounds,
    const uint32_t num_sh_buffer,
    const int sh_value_bits,
    const int64_t sh_bounds_stride
) {
    return launch_projection_fwd_vk(
        /*antialiased=*/false, num_splats, in_splats, max_sh_degree,
        viewmats, intrins, image_width, image_height, camera_model,
        dist_coeffs, radii, sh_value_packed, sh_value_bounds,
        num_sh_buffer, sh_value_bits, sh_bounds_stride);
}

std::tuple<
    DeviceTensor2D<float4>, DeviceTensor2D<float>, std::vector<DeviceTensorFloatND>
> projection_mip_forward(
    const int64_t num_splats, const int max_sh_degree,
    const std::vector<DeviceTensorFloatND> &in_splats,
    TorchTensorView viewmats, TorchTensorView intrins,
    const uint32_t image_width, const uint32_t image_height,
    const std::string camera_model, const TorchTensorView dist_coeffs,
    DeviceVector<float> radii,
    const std::optional<TorchTensorView> sh_value_packed,
    const std::optional<TorchTensorView> sh_value_bounds,
    const uint32_t num_sh_buffer,
    const int sh_value_bits,
    const int64_t sh_bounds_stride
) {
    return launch_projection_fwd_vk(
        /*antialiased=*/true, num_splats, in_splats, max_sh_degree,
        viewmats, intrins, image_width, image_height, camera_model,
        dist_coeffs, radii, sh_value_packed, sh_value_bounds,
        num_sh_buffer, sh_value_bits, sh_bounds_stride);
}

std::tuple<
    DeviceTensor2D<float4>, DeviceTensor2D<float>, std::vector<DeviceTensorFloatND>
> projection_3dgut_forward(
    const int64_t num_splats, const int max_sh_degree,
    const std::vector<DeviceTensorFloatND> &in_splats,
    TorchTensorView viewmats, TorchTensorView intrins,
    const uint32_t image_width, const uint32_t image_height,
    const std::string camera_model, const TorchTensorView dist_coeffs,
    DeviceVector<float> radii,
    const std::optional<TorchTensorView> sh_value_packed,
    const std::optional<TorchTensorView> sh_value_bounds,
    const uint32_t num_sh_buffer,
    const int sh_value_bits,
    const int64_t sh_bounds_stride
) {
    uint64_t q_packed, q_bounds;
    int64_t q_stride;
    const uint32_t spec_bits =
        resolve_sh_quant(sh_value_packed, sh_value_bounds, num_sh_buffer,
                         sh_value_bits, sh_bounds_stride, &q_packed,
                         &q_bounds, &q_stride);

    CameraModelType cam = cmt(camera_model);
    if ((int)cam < 0 || (int)cam > 3)
        throw std::runtime_error("Unsupported camera model");

    const int64_t N = num_splats;
    Vanilla3DGUT<0>::WorldBuffer wb(in_splats);
    int sh_degree = wb.sh_degree();
    sh_degree = std::min(sh_degree, max_sh_degree);

    const uint32_t C = (uint32_t)std::get<2>(viewmats)[0];

    DeviceTensor2D<float4> aabb;
    aabb.resize(PoolSlot::ProjAabb, C, N);
    DeviceTensor2D<float> sorting_depths;
    sorting_depths.resize(PoolSlot::ProjDepths, C, N);
    std::vector<DeviceTensorFloatND> splats_screen =
        Vanilla3DGUT<0>::ScreenBuffer::empty_pool(C * N,
                                                  PoolSlot::ProjScreen);
    Vanilla3DGUT<0>::ScreenBuffer sb(splats_screen);

    const uint64_t total_threads = (uint64_t)C * (uint64_t)N;
    const uint32_t wgs = (uint32_t)((total_threads + 127) / 128);
    const uint32_t per_row = std::min(wgs, 65535u);
    const uint32_t rows = (wgs + per_row - 1) / per_row;

    Projection3dgutParams p{};
    p.means = (uint64_t)wb.raw_data(0);
    p.quats = (uint64_t)wb.raw_data(1);
    p.scales = (uint64_t)wb.raw_data(2);
    p.opacities = (uint64_t)wb.raw_data(3);
    p.features_dc = (uint64_t)wb.raw_data(4);
    p.features_sh = (uint64_t)wb.raw_data(5);
    p.sh_stride = (uint32_t)wb.raw_stride(5);
    p.viewmats = std::get<0>(viewmats);
    p.intrins = std::get<0>(intrins);
    p.dist_coeffs = std::get<0>(dist_coeffs);
    p.out_aabb = (uint64_t)aabb.data_ptr();
    p.out_depths = (uint64_t)sorting_depths.data_ptr();
    p.out_radii = (uint64_t)radii.data_ptr();
    p.s_scale = (uint64_t)sb.raw_data(0);
    p.s_opac = (uint64_t)sb.raw_data(1);
    p.s_rgb = (uint64_t)sb.raw_data(2);
    p.sh_packed = q_packed;
    p.sh_bounds = q_bounds;
    p.sh_bounds_stride = q_stride;
    p.num_sh_buffer = num_sh_buffer;
    p.C = C;
    p.N = (uint32_t)N;
    p.width = image_width;
    p.height = image_height;
    p.wgs_per_row = per_row;

    uint64_t params_addr = 0;
    void* params_mapped = nullptr;
    if (!backend::vk::params_alloc(sizeof(p), &params_addr, &params_mapped))
        throw std::runtime_error("Vulkan backend: params ring failed");
    std::memcpy(params_mapped, &p, sizeof(p));

    backend::vk::SpecList spec{(uint32_t)cam, (uint32_t)sh_degree, 0u,
                               spec_bits};
    if (!backend::vk::dispatch(backend::kDefaultStream,
                               "projection_fwd.projection_fwd_3dgut", spec,
                               per_row, rows, 1, &params_addr,
                               sizeof(params_addr)))
        throw std::runtime_error("Vulkan backend: projection dispatch failed");

    return std::make_tuple(aabb, sorting_depths, splats_screen);
}
