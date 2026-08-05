#pragma once
// Shared helpers for the Vulkan kernel launchers in backend/vulkan/kernels/.
// Each launcher TU mirrors one CUDA launch API (src/kernels/**/*.cuh); these helpers
// cover the idioms every TU repeats: 1D grid folding under the 65535
// per-dimension limit, dispatch that throws (with optional per-dispatch
// debug sync), the params-ring indirection for structs above the
// push-constant floor, and SH value-quant argument resolution shared by the
// projection launchers.

#include "backend/vulkan/VulkanInternal.h"
#include "backend/vulkan/VulkanPipelines.h"
#include "backend/common/SortScan.h"

#include <core/Tensor.h>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <optional>
#include <stdexcept>
#include <string>
#include "core/Env.h"

namespace vkk {

// A 1D range of `threads` threads folded into a (per_row, rows) workgroup
// grid: grid dimensions cap at 65535, so long ranges fold into rows and the
// shader recovers the flat index as
// ((gy * wgs_per_row + gx) * block + tid).
struct Fold {
    uint32_t per_row, rows;
};

inline Fold fold_1d(int64_t threads, uint32_t block) {
    uint32_t wgs = (uint32_t)((threads + block - 1) / block);
    uint32_t per_row = std::min(std::max(wgs, 1u), 65535u);
    return {per_row, (wgs + per_row - 1) / per_row};
}

// SS_VK_DEBUG_SYNC=1: device-sync + report after every dispatch, to
// bisect device-lost errors to a kernel. SS_VIS_DEBUG_SYNC is honored
// as a legacy alias.
inline bool debug_sync_enabled() {
    static const bool v = spirula::env("VK_DEBUG_SYNC") != nullptr ||
                          spirula::env("VIS_DEBUG_SYNC") != nullptr;
    return v;
}

// backend::vk::dispatch that throws (with the backend's sticky error detail)
// instead of returning false. `params` goes straight into push constants —
// use dispatch_ring for structs above the device push floor.
inline void dispatch(const char* entry, const backend::vk::SpecList& spec,
                     uint32_t gx, uint32_t gy, uint32_t gz,
                     const void* params, uint32_t size) {
    if (debug_sync_enabled())
        std::fprintf(stderr, "[vk-sync] %s (%u,%u,%u)...\n", entry, gx, gy,
                     gz);
    if (!backend::vk::dispatch(backend::kDefaultStream, entry, spec, gx, gy,
                               gz, params, size)) {
        const char* detail = backend::last_error();
        throw std::runtime_error(std::string("Vulkan backend: dispatch of ") +
                                 entry + " failed" +
                                 (detail ? std::string(": ") + detail : ""));
    }
    if (debug_sync_enabled()) {
        backend::device_synchronize();
        const char* err = backend::last_error();
        std::fprintf(stderr, "[vk-sync] %s: %s\n", entry, err ? err : "ok");
    }
}

// Flat-1D dispatch: folds `total` threads into `block`-wide groups. If
// `wgs_per_row_field` is given it is set inside the params struct before
// dispatch (the shader needs the fold to reconstruct the flat index).
// No-op when total <= 0.
inline void dispatch_flat(const char* entry,
                          const backend::vk::SpecList& spec, int64_t total,
                          uint32_t block, void* params, uint32_t size,
                          uint32_t* wgs_per_row_field = nullptr) {
    if (total <= 0) return;
    Fold f = fold_1d(total, block);
    if (wgs_per_row_field) *wgs_per_row_field = f.per_row;
    dispatch(entry, spec, f.per_row, f.rows, 1, params, size);
}

// Writes `params` into the kernel-params ring and dispatches with only the
// 8-byte device address pushed — for param structs above the push-constant
// floor. The shader side declares `uniform Params* p`.
inline void dispatch_ring(const char* entry,
                          const backend::vk::SpecList& spec, uint32_t gx,
                          uint32_t gy, uint32_t gz, const void* params,
                          uint32_t size) {
    uint64_t addr = 0;
    void* mapped = nullptr;
    if (!backend::vk::params_alloc(size, &addr, &mapped))
        throw std::runtime_error("Vulkan backend: params ring failed");
    std::memcpy(mapped, params, size);
    dispatch(entry, spec, gx, gy, gz, &addr, sizeof(addr));
}

// Device address of a small zeroed allocation that stands in for null
// optional pointers in kernel params. Rationale (learned on llvmpipe): the
// CPU JIT can SPECULATE loads that sit behind a pointer-null (or flag)
// branch, so a never-taken branch over a null device address still
// segfaults. Launchers therefore never put a null address in a params
// pointer field a shader might load through: "null means zeros" pointers
// (e.g. dist_coeffs) read the fallback with identical results, and "null
// means skip" semantics move to explicit flag fields or specialization
// constants. 4 KB covers every constant-offset load a JIT could hoist.
inline uint64_t null_fallback() {
    static const uint64_t addr = [] {
        constexpr size_t kBytes = 4096;
        void* p = backend::device_malloc(kBytes);
        if (!p)
            throw std::runtime_error(
                "Vulkan backend: null-fallback allocation failed");
        backend::memset_sync(p, 0, kBytes);
        return (uint64_t)p;
    }();
    return addr;
}

// addr if non-null, else the zeroed fallback.
inline uint64_t or_fallback(uint64_t addr) {
    return addr ? addr : null_fallback();
}
inline uint64_t or_fallback(const void* p) {
    return or_fallback((uint64_t)p);
}

// Validation + resolution of the SH value-quant launch args shared by the
// projection launchers (fwd, packed fwd; bwd reuses it in the training
// phase). Returns the kShValueBits spec value (0 = fp32). Mirrors the CUDA
// kernels' stride default: 0 -> FPBO per-splat-block layout
// (256 * 3 * num_sh_buffer). out_stride is never 0 (the shaders divide by
// it even in the fp32 specialization).
inline uint32_t resolve_sh_quant(
    const std::optional<TorchTensorView>& sh_value_packed,
    const std::optional<TorchTensorView>& sh_value_bounds,
    const uint32_t num_sh_buffer, const int sh_value_bits,
    const int64_t sh_bounds_stride, uint64_t* out_packed,
    uint64_t* out_bounds, int64_t* out_stride) {
    *out_packed = 0;
    *out_bounds = 0;
    *out_stride = 1;
    if (sh_value_bits == 32)
        return 0;
    if (sh_value_bits != 8 && sh_value_bits != 16)
        throw std::runtime_error(
            "projection: sh_value_bits must be 8, 16 or 32");
    if (!sh_value_packed.has_value() || !sh_value_bounds.has_value())
        throw std::runtime_error(
            "projection: sh_value_bits != 32 requires "
            "sh_value_packed and sh_value_bounds");
    *out_packed = std::get<0>(sh_value_packed.value());
    *out_bounds = std::get<0>(sh_value_bounds.value());
    *out_stride = sh_bounds_stride > 0
                      ? sh_bounds_stride
                      : (int64_t)256 * 3 * (int64_t)num_sh_buffer;
    return (uint32_t)sh_value_bits;
}

// Packed-mode preprocessing shared by the grad-quant and FPBO projection
// backward launchers: identity permutation sorted by gaussian id (so the
// splat-parallel kernels can loop each splat's intersections) plus the
// per-splat [N+1] camera ranges. Mirrors the CUDA launchers' CUB steps; the
// small kernels live in projection_qgrad.slang. NOTE: like CUB, the radix
// sort may clobber the input gaussian_ids buffer (ping-pong).
struct PackedCameraRanges {
    DeviceVector<int32_t> camera_id_bounds;
    DeviceVector<int32_t> gauss_sorted, perm, perm_sorted;
    const int32_t* sorted_perm = nullptr;
};

inline PackedCameraRanges build_packed_camera_ranges(
    const DeviceVector<int32_t>& gaussian_ids, int64_t N) {
    struct IotaParams {
        uint64_t buf;
        uint32_t n;
        uint32_t wgs_per_row;
    };
    struct CamBoundsParams {
        uint64_t gaussian_ids_sorted;
        uint64_t camera_id_bounds;
        uint32_t nnz;
        uint32_t N;
        uint32_t wgs_per_row;
        uint32_t _pad0;
    };

    PackedCameraRanges out;
    const int64_t nnz = gaussian_ids.size();
    out.gauss_sorted.resize(PoolSlot::FusedProjBwdGaussSorted, nnz);
    out.perm.resize(PoolSlot::FusedProjBwdPerm, nnz);
    out.perm_sorted.resize(PoolSlot::FusedProjBwdPermSorted, nnz);

    IotaParams ip{};
    ip.buf = (uint64_t)out.perm.data_ptr();
    ip.n = (uint32_t)nnz;
    dispatch_flat("projection_qgrad.qgrad_iota", {}, nnz, 256, &ip,
                  sizeof(ip), &ip.wgs_per_row);

    backend::DoubleBuffer<int32_t> d_keys(
        const_cast<int32_t*>(gaussian_ids.data_ptr()),
        out.gauss_sorted.data_ptr());
    backend::DoubleBuffer<int32_t> d_values(out.perm.data_ptr(),
                                            out.perm_sorted.data_ptr());
    int n_bits = 0;
    while ((1u << n_bits) <= (uint64_t)N)
        ++n_bits;
    backend::sort_pairs(d_keys, d_values, nnz, 0, n_bits);
    out.sorted_perm = d_values.current();

    out.camera_id_bounds.resize(PoolSlot::FusedProjBwdCamBounds, N + 1);
    CamBoundsParams bp{};
    bp.gaussian_ids_sorted = (uint64_t)d_keys.current();
    bp.camera_id_bounds = (uint64_t)out.camera_id_bounds.data_ptr();
    bp.nnz = (uint32_t)nnz;
    bp.N = (uint32_t)N;
    dispatch_flat("projection_qgrad.qgrad_camera_id_bounds", {}, nnz + 1, 256,
                  &bp, sizeof(bp), &bp.wgs_per_row);
    return out;
}

}  // namespace vkk
