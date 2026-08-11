#pragma once

// Ray -> input-image pixel, for the two ways a GT image can be resampled.
//
// `kFromSource == false` is the ordinary path: the input image belongs to a
// camera the engine can project directly. `kFromSource == true` is a camera
// whose lens model no distortion tier represents (COLMAP FOV / DIVISION /
// EUCM / RAD_TAN_THIN_PRISM_FISHEYE); the parser fitted it onto a tier, and the
// pixels have to come from the TRUE source projection, not from the fit.
//
// Both warps take this as a template argument so an ordinary dataset pays
// neither the branch nor the 16 registers the source parameters occupy.

#include "core/Common.cuh"
#include "core/CameraDistortion.cuh"

#ifdef __CUDACC__

template<CameraDistortionType D, bool kFromSource>
struct RayToPixel {
    CameraModelType camera_model;
    float4 intrin;
    typename SlangDistortion<D>::Coeffs coeffs;
    int    source_model  = -1;
    const float* __restrict__ source_params = nullptr;   // [16], this image's row

    __device__ __forceinline__ bool operator()(float3 raydir, float2* uv) const {
        if constexpr (kFromSource) {
            FixedArray<float, 16> p;
            #pragma unroll
            for (int i = 0; i < 16; i++) p[i] = source_params[i];
            return SlangCameraSource::source_project(source_model, p, raydir, uv);
        } else {
            return camera_proj_nav<D>(camera_model, raydir, intrin, coeffs, uv);
        }
    }
};

// Build one for image `bid`. `source_models` / `source_params` are null on the
// ordinary path.
template<CameraDistortionType D, bool kFromSource>
__device__ __forceinline__ RayToPixel<D, kFromSource> make_ray_to_pixel(
    long bid,
    CameraModelType camera_model,
    const float4* __restrict__ intrins,
    const CameraDistortionCoeffsBuffer& dist_coeffs_buffer,
    const int* __restrict__ source_models,
    const float* __restrict__ source_params
) {
    RayToPixel<D, kFromSource> m;
    m.camera_model = camera_model;
    if constexpr (kFromSource) {
        m.source_model  = source_models[bid];
        m.source_params = source_params + 16 * bid;
    } else {
        m.intrin = intrins[bid];
        m.coeffs = dist_coeffs_buffer.load<D>(bid);
    }
    return m;
}

// _SS_DISPATCH_DISTORTION with the source axis folded in; BODY takes
// (tier, from_source). Defined here rather than in PixelWiseCommon.cuh so a TU
// that never re-distorts does not see it.
#define _SS_DISPATCH_SOURCE_CASE(BODY, TIER, HAS)                              \
    case CameraDistortionType::TIER:                                           \
        /* braced: BODY is itself an if/else in the depth and normal warps */  \
        if (HAS) { BODY(CameraDistortionType::TIER, true); }                   \
        else     { BODY(CameraDistortionType::TIER, false); }                  \
        break;

#define _SS_DISPATCH_SOURCE(name, has_source, BODY)                            \
    do { const bool _ss_has = (has_source);                                    \
         switch (cdt(name)) {                                                  \
        _SS_DISPATCH_SOURCE_CASE(BODY, None, _ss_has)                          \
        _SS_DISPATCH_SOURCE_CASE(BODY, OpenCV, _ss_has)                        \
        _SS_DISPATCH_SOURCE_CASE(BODY, ThinPrism, _ss_has)                     \
        _SS_DISPATCH_SOURCE_CASE(BODY, Rational, _ss_has)                      \
        default: throw std::runtime_error(                                     \
            "Unknown camera distortion: " + std::string(name));                \
    } } while (0)

#endif  // __CUDACC__
