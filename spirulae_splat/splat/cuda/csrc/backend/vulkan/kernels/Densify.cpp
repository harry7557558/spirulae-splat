// Vulkan implementation of the densification launch APIs (csrc/Densify.cuh).
// Currently densify_update_weight only (called from the optimizer step);
// splitting/relocation/MCMC land with the densify/MCMC port. Device work:
// slang/vulkan/densify.slang.

#include <Densify.cuh>

#include "KernelCommon.h"

namespace {

// Mirrors DensifyUpdateWeightParams in slang/vulkan/densify.slang.
struct DensifyUpdateWeightParams {
    uint64_t radii, opacs, accum_weight, accum_weight2, accum_buffer;
    float blend_w;
    int32_t score_mode;
    uint32_t use_opacs, use_w2, num_splats, wgs_per_row;
};
static_assert(sizeof(DensifyUpdateWeightParams) == 5 * 8 + 6 * 4,
              "params layout must match the slang struct");

}  // namespace

/* API definitions matching csrc/Densify.cuh (engine-referenced subset) */

void densify_update_weight(
    int64_t num_splats,
    DeviceVector<float> radii,
    float3* scales_ptr,
    float* opacs_ptr,
    DeviceVector<float> accum_weight,
    DeviceVector<float> accum_weight2,
    float blend_w,
    DeviceVector<float2> accum_buffer,
    int score_mode
) {
    (void)scales_ptr;  // unused by the kernel, as in CUDA
    if (num_splats <= 0) return;
    DensifyUpdateWeightParams p{};
    p.radii = (uint64_t)radii.data_ptr();
    p.opacs = vkk::or_fallback(opacs_ptr);
    p.accum_weight = (uint64_t)accum_weight.data_ptr();
    p.accum_weight2 = vkk::or_fallback(accum_weight2.data_ptr());
    p.accum_buffer = (uint64_t)accum_buffer.data_ptr();
    p.blend_w = blend_w;
    p.score_mode = score_mode;
    p.use_opacs = opacs_ptr ? 1u : 0u;
    p.use_w2 = accum_weight2.data_ptr() ? 1u : 0u;
    p.num_splats = (uint32_t)num_splats;
    vkk::dispatch_flat("densify.densify_update_weight",
                       backend::vk::SpecList{}, num_splats, 256, &p,
                       sizeof(p), &p.wgs_per_row);
}
