// Vulkan implementation of the packed projection-forward launch API
// (csrc/ProjectionPackedFwd.cuh). Two-pass compaction mirroring
// ProjectionPackedFwd.cu, with two deliberate substitutions:
//  - the mask is int32 0/1 (bool would need 8-bit stores; see README rules)
//  - the prefix sum is backend::inclusive_sum<int32> (SortScan) instead of
//    CUB, with a 4-byte nnz readback.

#include <ProjectionPackedFwd.cuh>

#include "../../common/SortScan.h"
#include "../VulkanInternal.h"
#include "../VulkanPipelines.h"

#include <cstring>
#include <stdexcept>

namespace {

// Mirrors PackedMaskParams in slang/vulkan/projection_fwd.slang.
struct PackedMaskParams {
    uint64_t means, quats, scales, opacities;
    uint64_t viewmats, intrins, dist_coeffs;
    uint64_t out_mask;
    uint32_t C, N, width, height, wgs_per_row;
};
static_assert(sizeof(PackedMaskParams) == 8 * 8 + 5 * 4 + 4 /*pad*/,
              "params layout must match the slang struct");

// Mirrors PackedFwdParams in slang/vulkan/projection_fwd.slang.
struct PackedFwdParams {
    uint64_t means, quats, scales, opacities, features_dc, features_sh;
    uint64_t viewmats, intrins, dist_coeffs;
    uint64_t mask_scan;
    uint64_t out_camera_ids, out_gaussian_ids, out_aabb, out_depths,
        out_radii;
    uint64_t s0, s1, s2, s3, s4;
    uint32_t sh_stride;
    uint32_t C, N, width, height, wgs_per_row;
};
static_assert(sizeof(PackedFwdParams) == 20 * 8 + 6 * 4,
              "params layout must match the slang struct");

std::tuple<DeviceVector<int32_t>, DeviceVector<int32_t>, DeviceVector<float4>,
           DeviceVector<float>, std::vector<DeviceTensorFloatND>>
launch_projection_packed_vk(
    const bool eval3d, const bool antialiased,
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

    const uint64_t total_threads = (uint64_t)C * (uint64_t)N;
    const uint32_t wgs = (uint32_t)((total_threads + 127) / 128);
    const uint32_t per_row = std::min(std::max(wgs, 1u), 65535u);
    const uint32_t rows = (wgs + per_row - 1) / per_row;

    // --- pass 1: intersection mask (int32 0/1) + inclusive scan + nnz ---
    DeviceVector<int32_t> mask;
    mask.resize(PoolSlot::ProjMask, (int64_t)C * N);
    DeviceVector<int32_t> mask_scan;
    mask_scan.resize(PoolSlot::ProjScan, (int64_t)C * N);

    PackedMaskParams mp{};
    mp.means = (uint64_t)wb.raw_data(0);
    mp.quats = (uint64_t)wb.raw_data(1);
    mp.scales = (uint64_t)wb.raw_data(2);
    mp.opacities = (uint64_t)wb.raw_data(3);
    mp.viewmats = std::get<0>(viewmats);
    mp.intrins = std::get<0>(intrins);
    mp.dist_coeffs = std::get<0>(dist_coeffs);
    mp.out_mask = (uint64_t)mask.data_ptr();
    mp.C = C;
    mp.N = (uint32_t)N;
    mp.width = image_width;
    mp.height = image_height;
    mp.wgs_per_row = per_row;

    // Spec IDs: 0 = camera model, 1 = SH degree (unused by the mask),
    // 2 = antialiased, 3 = eval3d.
    backend::vk::SpecList spec{(uint32_t)cam, (uint32_t)sh_degree,
                               antialiased ? 1u : 0u, eval3d ? 1u : 0u};
    if (!backend::vk::dispatch(backend::kDefaultStream,
                               "projection_fwd.projection_packed_mask", spec,
                               per_row, rows, 1, &mp, sizeof(mp)))
        throw std::runtime_error("Vulkan backend: packed mask dispatch failed");

    backend::inclusive_sum<int32_t>(mask.data_ptr(), mask_scan.data_ptr(),
                                    (int64_t)C * N);
    int32_t nnz = 0;
    if ((int64_t)C * N > 0)
        backend::memcpy_sync(&nnz, mask_scan.data_ptr() + ((int64_t)C * N - 1),
                             sizeof(int32_t),
                             backend::MemcpyKind::DeviceToHost);

    // --- pass 2: compacted projection ---
    DeviceVector<int32_t> camera_ids;
    camera_ids.resize(PoolSlot::ProjCameraIds, nnz);
    DeviceVector<int32_t> gaussian_ids;
    gaussian_ids.resize(PoolSlot::ProjGaussianIds, nnz);
    DeviceVector<float4> aabb;
    aabb.resize(PoolSlot::ProjAabb, nnz);
    DeviceVector<float> sorting_depths;
    sorting_depths.resize(PoolSlot::ProjDepths, nnz);
    std::vector<DeviceTensorFloatND> splats_screen =
        eval3d ? Vanilla3DGUT<0>::ScreenBuffer::empty_pool(
                     nnz, PoolSlot::ProjScreen)
               : Vanilla3DGS<0>::ScreenBuffer::empty_pool(
                     nnz, PoolSlot::ProjScreen);

    if (nnz > 0) {
        PackedFwdParams p{};
        p.means = mp.means;
        p.quats = mp.quats;
        p.scales = mp.scales;
        p.opacities = mp.opacities;
        p.features_dc = (uint64_t)wb.raw_data(4);
        p.features_sh = (uint64_t)wb.raw_data(5);
        p.sh_stride = (uint32_t)wb.raw_stride(5);
        p.viewmats = mp.viewmats;
        p.intrins = mp.intrins;
        p.dist_coeffs = mp.dist_coeffs;
        p.mask_scan = (uint64_t)mask_scan.data_ptr();
        p.out_camera_ids = (uint64_t)camera_ids.data_ptr();
        p.out_gaussian_ids = (uint64_t)gaussian_ids.data_ptr();
        p.out_aabb = (uint64_t)aabb.data_ptr();
        p.out_depths = (uint64_t)sorting_depths.data_ptr();
        p.out_radii = (uint64_t)radii.data_ptr();
        if (eval3d) {
            Vanilla3DGUT<0>::ScreenBuffer sb(splats_screen);
            p.s0 = (uint64_t)sb.raw_data(0);
            p.s1 = (uint64_t)sb.raw_data(1);
            p.s2 = (uint64_t)sb.raw_data(2);
        } else {
            Vanilla3DGS<0>::ScreenBuffer sb(splats_screen);
            p.s0 = (uint64_t)sb.raw_data(0);
            p.s1 = (uint64_t)sb.raw_data(1);
            p.s2 = (uint64_t)sb.raw_data(2);
            p.s3 = (uint64_t)sb.raw_data(3);
            p.s4 = (uint64_t)sb.raw_data(4);
        }
        p.C = C;
        p.N = (uint32_t)N;
        p.width = image_width;
        p.height = image_height;
        p.wgs_per_row = per_row;

        uint64_t params_addr = 0;
        void* params_mapped = nullptr;
        if (!backend::vk::params_alloc(sizeof(p), &params_addr,
                                       &params_mapped))
            throw std::runtime_error("Vulkan backend: params ring failed");
        std::memcpy(params_mapped, &p, sizeof(p));
        if (!backend::vk::dispatch(backend::kDefaultStream,
                                   "projection_fwd.projection_packed_fwd",
                                   spec, per_row, rows, 1, &params_addr,
                                   sizeof(params_addr)))
            throw std::runtime_error(
                "Vulkan backend: packed fwd dispatch failed");
    }

    return std::make_tuple(camera_ids, gaussian_ids, aabb, sorting_depths,
                           splats_screen);
}

}  // namespace

/* API definitions matching csrc/ProjectionPackedFwd.cuh */

#define _SSPLAT_DEF_PACKED(name, eval3d, antialiased)                        \
std::tuple<                                                                  \
    DeviceVector<int32_t>, DeviceVector<int32_t>, DeviceVector<float4>,      \
    DeviceVector<float>, std::vector<DeviceTensorFloatND>                    \
> name(                                                                      \
    const int64_t num_splats,                                                \
    const int max_sh_degree,                                                 \
    const std::vector<DeviceTensorFloatND> &in_splats,                       \
    TorchTensorView viewmats,                                                \
    TorchTensorView intrins,                                                 \
    const uint32_t image_width,                                              \
    const uint32_t image_height,                                             \
    const std::string camera_model,                                          \
    const TorchTensorView dist_coeffs,                                       \
    DeviceVector<float> radii,                                               \
    const std::optional<TorchTensorView> sh_value_packed,                    \
    const std::optional<TorchTensorView> sh_value_bounds,                    \
    const uint32_t num_sh_buffer,                                            \
    const int sh_value_bits,                                                 \
    const int64_t sh_bounds_stride                                           \
) {                                                                          \
    (void)sh_value_bounds; (void)num_sh_buffer; (void)sh_bounds_stride;      \
    return launch_projection_packed_vk(                                      \
        eval3d, antialiased, num_splats, in_splats, max_sh_degree,           \
        viewmats, intrins, image_width, image_height, camera_model,          \
        dist_coeffs, radii, sh_value_packed, sh_value_bits);                 \
}

_SSPLAT_DEF_PACKED(projection_3dgs_packed_forward, false, false)
_SSPLAT_DEF_PACKED(projection_mip_packed_forward, false, true)
_SSPLAT_DEF_PACKED(projection_3dgut_packed_forward, true, false)

#undef _SSPLAT_DEF_PACKED
