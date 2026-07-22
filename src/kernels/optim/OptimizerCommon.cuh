#pragma once

// Shared preamble for the optimizer kernels.
//
// Optimizer.cu was one file (~2.7k lines) holding every optimizer variant.
// It is split by function; the parts all declare into the single
// Optimizer.cuh via HEADER_SOURCES['Optimizer'] in generate_headers.py.
// See docs/codegen.md.
//
//   TensorSetZero.cu           set_zero_tensor
//   AdamOptim.cu               Adam, multi-tensor Adam, Riemannian quat Adam,
//                              8-bit / stepped variants
//   NewtonOptim.cu             Newton
//   ScaleAgnosticMeanOptim.cu  scale-agnostic mean Adam
//   FusedGeometryOptim.cu      fused 3DGS geometry optimizer
//   FusedAppearanceOptim.cu    fused appearance optimizer (quantized SH)
//   TrustRegion3DGS2Optim.cu   3DGS^2-TR mean/scale/color/opacity/quaternion
//   ColorOptim.cu              linear-RGB Adam, trust-region RGB / SH Adam

#include "kernels/optim/Optimizer.cuh"

#include "generated/slang.cuh"
namespace SlangProjectionUtils {
#include "generated/set_namespace.cuh"
#include "generated/projection_utils.cuh"
}
namespace SlangPixelWise {
#include "generated/set_namespace.cuh"
#include "generated/pixel_wise.cuh"
}
namespace SlangPerSplatLosses {
#include "generated/set_namespace.cuh"
#include "generated/per_splat_losses.cuh"
}

#include <core/Common.cuh>
#include <core/Tensor.h>
#include "core/NonShQuantState.h"
#include "core/GradQuant.cuh"
#include "kernels/projection/ProjectionBwdQuantGrad.cuh"   // GradQuantBuffers
#include <stdexcept>
#include <string>

#include <cooperative_groups.h>
namespace cg = cooperative_groups;

#if defined(__INTEL_COMPILER) || defined(__GNUC__) || defined(__clang__) || defined(_MSC_VER)
#if defined(__AVX2__) || defined(__SSE2__)
#include <immintrin.h>
#endif
#endif




inline constexpr float kSh0 = 0.28209479177387814f;

// Helper: compute total number of elements from a TorchTensorView
static inline int64_t _tv_numel(const TorchTensorView& tv) {
    int64_t n = 1;
    for (auto d : std::get<2>(tv)) n *= d;
    return n;
}
// Helper: get float* from TorchTensorView
static inline float* _tv_f(const TorchTensorView& tv) { return (float*)std::get<0>(tv); }

__forceinline__ __device__ float3 fmaxf(float3 v, float k) {
    return {
        fmaxf(v.x, k),
        fmaxf(v.y, k),
        fmaxf(v.z, k)
    };
}
__forceinline__ __device__ float3 sqrtf(float3 v) {
    return {
        sqrtf(fmaxf(v.x, 0.0f)),
        sqrtf(fmaxf(v.y, 0.0f)),
        sqrtf(fmaxf(v.z, 0.0f))
    };
}

__forceinline__ __device__ float4 sqrtf(float4 v) {
    return {
        sqrtf(fmaxf(v.x, 0.0f)),
        sqrtf(fmaxf(v.y, 0.0f)),
        sqrtf(fmaxf(v.z, 0.0f)),
        sqrtf(fmaxf(v.w, 0.0f))
    };
}
