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
    uint32_t sh_stride;
    uint32_t C, N;
    uint32_t width, height;
    uint32_t wgs_per_row;
};
static_assert(sizeof(ProjectionFwdParams) == 17 * 8 + 6 * 4,
              "params layout must match the slang struct");

// Mirrors Projection3dgutParams in slang/vulkan/projection_fwd.slang.
struct Projection3dgutParams {
    uint64_t means, quats, scales, opacities, features_dc, features_sh;
    uint64_t viewmats, intrins, dist_coeffs;
    uint64_t out_aabb, out_depths, out_radii;
    uint64_t s_scale, s_opac, s_rgb;
    uint32_t sh_stride;
    uint32_t C, N;
    uint32_t width, height;
    uint32_t wgs_per_row;
};
static_assert(sizeof(Projection3dgutParams) == 15 * 8 + 6 * 4,
              "params layout must match the slang struct");

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
    const int sh_value_bits) {
    if (sh_value_bits != 32 || sh_value_packed.has_value())
        throw std::runtime_error(
            "Vulkan backend: SH value-quantization (sh_value_bits != 32) is "
            "not implemented yet");

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
                               antialiased ? 1u : 0u};
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
    (void)sh_value_bounds; (void)num_sh_buffer; (void)sh_bounds_stride;
    return launch_projection_fwd_vk(
        /*antialiased=*/false, num_splats, in_splats, max_sh_degree,
        viewmats, intrins, image_width, image_height, camera_model,
        dist_coeffs, radii, sh_value_packed, sh_value_bits);
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
    (void)sh_value_bounds; (void)num_sh_buffer; (void)sh_bounds_stride;
    return launch_projection_fwd_vk(
        /*antialiased=*/true, num_splats, in_splats, max_sh_degree,
        viewmats, intrins, image_width, image_height, camera_model,
        dist_coeffs, radii, sh_value_packed, sh_value_bits);
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
    (void)sh_value_bounds; (void)num_sh_buffer; (void)sh_bounds_stride;
    if (sh_value_bits != 32 || sh_value_packed.has_value())
        throw std::runtime_error(
            "Vulkan backend: SH value-quantization (sh_value_bits != 32) is "
            "not implemented yet");

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

    backend::vk::SpecList spec{(uint32_t)cam, (uint32_t)sh_degree, 0u};
    if (!backend::vk::dispatch(backend::kDefaultStream,
                               "projection_fwd.projection_fwd_3dgut", spec,
                               per_row, rows, 1, &params_addr,
                               sizeof(params_addr)))
        throw std::runtime_error("Vulkan backend: projection dispatch failed");

    return std::make_tuple(aabb, sorting_depths, splats_screen);
}
