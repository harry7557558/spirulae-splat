#pragma once

#include "backend/api/BackendTypes.h"
#include <cstdint>
#include <cmath>


// Per-cell decode wrapper for bilagrid value (16-bit linear quantization).
// Bilagrid sampler kernels take this by value and read cells via operator[],
// so the existing call sites `bilagrid[idx]` stay unchanged.
//
// Storage selector: when `fp32` is non-null the reader passes through;
// otherwise it decodes from the 16-bit packed buffer using per-block bounds.
// Block size is hard-coded to 256 (matches QuantizedTensor<16, 256>).
//
// The branch is uniform across the warp (every thread sees the same
// pointers), so the decode path adds ~5 ops per read but no warp divergence.
struct BilagridReader {
    const float*    __restrict__ fp32   = nullptr;
    const uint16_t* __restrict__ q16    = nullptr;
    const float2*   __restrict__ bounds = nullptr;

    // Implicit ctor from a raw fp32 pointer (and default ctor) so existing
    // launcher call sites that pass `grids.data_ptr()` keep working without
    // any change. The 16-bit path is opt-in via the static factory.
    BilagridReader() = default;
    __host__ __device__ BilagridReader(const float* p)
        : fp32(p), q16(nullptr), bounds(nullptr) {}

    __device__ __forceinline__ float operator[](int64_t idx) const {
        if (fp32 != nullptr) return fp32[idx];
        float2 mm = bounds[idx >> 8];          // 256-cell blocks
        float v   = (float)q16[idx];
        return mm.x + (mm.y - mm.x) * (v * (1.0f / 65535.0f));
    }

    // Pointer view for kernels that need raw addresses (TV neighbor reads
    // in BilagridFusedAdam). Unused on the quant path -- callers branch on
    // `value_quantized()` first.
    __device__ __forceinline__ const float* raw_fp32() const { return fp32; }
    __device__ __forceinline__ bool value_quantized() const { return fp32 == nullptr; }

    static __host__ __device__ __forceinline__ BilagridReader from_q16(
        const uint16_t* q, const float2* b
    ) {
        BilagridReader r; r.q16 = q; r.bounds = b; return r;
    }
};


// 16-bit linear-codec encode of a single cell value into the packed buffer.
// Endpoint-exact within [mm.x, mm.y]; out-of-range values clamp.
__device__ __forceinline__ void bilagrid_encode_q16(
    uint16_t* __restrict__ q16, int64_t idx, float val, float2 mm
) {
    float range = fmaxf(mm.y - mm.x, 1e-30f);
    float qf    = roundf(65535.0f * (val - mm.x) / range);
    q16[idx]    = (uint16_t)fminf(fmaxf(qf, 0.0f), 65535.0f);
}
