#pragma once

// CUDA-free camera-model enum + name helpers.
//
// CameraModelType and camera_model_from_name were historically defined in
// Common.cuh / Camera.h, which under the CUDA backend transitively pull in
// <cuda_runtime.h> (via backend/api/BackendTypes.h). The
// standalone dataset parsers (data/parsers/) only need these two trivial
// symbols, and are also compiled for WebAssembly (viewer/), where CUDA headers
// are unavailable. They include THIS header instead of Camera.h.
//
// The CUDA engine build is unaffected: it never includes this file (it keeps
// using Common.cuh's CameraModelType and Camera.h's helpers), and no single
// translation unit includes both this header and Common.cuh/Camera.h.
//
// The enum values MUST stay in sync with Common.cuh:130 and
// projection_utils.slang.

#include <string>

enum class CameraModelType {
    PINHOLE = 0,
    FISHEYE = 1,
    EQUISOLID = 2,
    EQUIRECTANGULAR = 3,
};

// Lens distortion, orthogonal to CameraModelType and chosen per camera group at
// load time: the parsers pick the cheapest tier that represents the source
// camera exactly, and fit + re-distort when none does. Values must match
// CameraDistortionType in shaders/projection_utils.slang and Common.cuh.
//
// Coefficient order is each tier's own -- a slot index does NOT mean the same
// thing across tiers:
//   OpenCV      k1 k2 p1 p2                 (COLMAP OPENCV)
//   ThinPrism   k1 k2 k3 k4 p1 p2 sx1 sy1   (COLMAP THIN_PRISM_FISHEYE)
//   Rational    k1 k2 k3 k4 k5 k6 p1 p2     (COLMAP FULL_OPENCV; k4..k6 divide)
enum class CameraDistortionType {
    None = 0,
    OpenCV = 1,
    ThinPrism = 2,
    Rational = 3,
};

// Storage width of the per-camera coefficient row. Every tier reads a prefix of
// it; the tail is zero. MUST match kCameraDistortionParams in Common.cuh.
inline constexpr int kCameraDistortionParams = 8;

inline constexpr int camera_distortion_num_params(CameraDistortionType d) {
    return d == CameraDistortionType::None      ? 0 :
           d == CameraDistortionType::OpenCV    ? 4 : 8;
}

inline const char* camera_distortion_to_string(CameraDistortionType d) {
    switch (d) {
        case CameraDistortionType::None:      return "NONE";
        case CameraDistortionType::OpenCV:    return "OPENCV";
        case CameraDistortionType::ThinPrism: return "THIN_PRISM";
        case CameraDistortionType::Rational:  return "RATIONAL";
        default:                              return "UNKNOWN";
    }
}

inline CameraDistortionType camera_distortion_from_name(const std::string& name) {
    if (name == "NONE")       return CameraDistortionType::None;
    if (name == "OPENCV")     return CameraDistortionType::OpenCV;
    if (name == "THIN_PRISM") return CameraDistortionType::ThinPrism;
    if (name == "RATIONAL")   return CameraDistortionType::Rational;
    return (CameraDistortionType)-1;
}

// Cheapest tier that still represents `coeffs` exactly, given the tier it was
// written in. Zero-extension OpenCV -> ThinPrism is exact, so a ThinPrism
// camera whose k3/k4/sx1/sy1 vanish demotes; a Rational camera whose
// denominator vanishes is a plain polynomial and demotes as well.
inline CameraDistortionType camera_distortion_demote(
    CameraDistortionType tier, const float* coeffs, float* out)
{
    auto zero = [](float v) { return v == 0.0f; };
    if (tier == CameraDistortionType::Rational) {
        // 1/(1 + 0) == 1, so an all-zero denominator is the polynomial form
        // with k1..k3 and no k4 term.
        if (zero(coeffs[3]) && zero(coeffs[4]) && zero(coeffs[5])) {
            float k1 = coeffs[0], k2 = coeffs[1], k3 = coeffs[2];
            float p1 = coeffs[6], p2 = coeffs[7];
            if (zero(k3)) {
                out[0] = k1; out[1] = k2; out[2] = p1; out[3] = p2;
                for (int i = 4; i < kCameraDistortionParams; i++) out[i] = 0.0f;
                return CameraDistortionType::OpenCV;
            }
            out[0] = k1; out[1] = k2; out[2] = k3; out[3] = 0.0f;
            out[4] = p1; out[5] = p2; out[6] = 0.0f; out[7] = 0.0f;
            return CameraDistortionType::ThinPrism;
        }
    } else if (tier == CameraDistortionType::ThinPrism) {
        if (zero(coeffs[2]) && zero(coeffs[3]) && zero(coeffs[6]) && zero(coeffs[7])) {
            float k1 = coeffs[0], k2 = coeffs[1], p1 = coeffs[4], p2 = coeffs[5];
            out[0] = k1; out[1] = k2; out[2] = p1; out[3] = p2;
            for (int i = 4; i < kCameraDistortionParams; i++) out[i] = 0.0f;
            return camera_distortion_demote(CameraDistortionType::OpenCV, out, out);
        }
    } else if (tier == CameraDistortionType::OpenCV) {
        if (zero(coeffs[0]) && zero(coeffs[1]) && zero(coeffs[2]) && zero(coeffs[3])) {
            for (int i = 0; i < kCameraDistortionParams; i++) out[i] = 0.0f;
            return CameraDistortionType::None;
        }
    } else {
        for (int i = 0; i < kCameraDistortionParams; i++) out[i] = 0.0f;
        return CameraDistortionType::None;
    }
    if (out != coeffs)
        for (int i = 0; i < kCameraDistortionParams; i++) out[i] = coeffs[i];
    return tier;
}

// Rewrite `coeffs` from `from` into `to`, which must be a tier that contains it
// (see camera_distortion_kernel_tier). Returns false when it does not.
inline bool camera_distortion_promote(CameraDistortionType from,
                                      CameraDistortionType to,
                                      const float* coeffs, float* out) {
    if (from == to) {
        if (out != coeffs)
            for (int i = 0; i < kCameraDistortionParams; i++) out[i] = coeffs[i];
        return true;
    }
    float in[kCameraDistortionParams];
    for (int i = 0; i < kCameraDistortionParams; i++) in[i] = coeffs[i];
    for (int i = 0; i < kCameraDistortionParams; i++) out[i] = 0.0f;
    if (from == CameraDistortionType::None)
        return true;   // all-zero coefficients are the identity in every tier
    if (from == CameraDistortionType::OpenCV && to == CameraDistortionType::ThinPrism) {
        out[0] = in[0]; out[1] = in[1];          // k1 k2
        out[4] = in[2]; out[5] = in[3];          // p1 p2
        return true;
    }
    if (from == CameraDistortionType::OpenCV && to == CameraDistortionType::Rational) {
        out[0] = in[0]; out[1] = in[1];          // k1 k2 (numerator)
        out[6] = in[2]; out[7] = in[3];          // p1 p2
        return true;
    }
    // ThinPrism -> Rational would need a 4th numerator term and drop sx1/sy1;
    // Rational -> ThinPrism cannot express the denominator at all.
    return false;
}

// Which (camera model, distortion tier) pairs the CUDA kernels are instantiated
// for. Every tier is available on Vulkan (a specialization constant costs no
// SPIR-V), so this bounds the CUDA binary only, and it must match
// kCameraVariants in tools/codegen/generate_kernel_instantiation.py and the
// export lists in shaders/primitive_3dgs.slang.
//
// The gaps are not arbitrary: no COLMAP fisheye model is rational, and
// EQUIRECTANGULAR has no lens distortion.
inline bool camera_distortion_is_compiled(CameraModelType m, CameraDistortionType d) {
    switch (m) {
        case CameraModelType::PINHOLE:         return true;
        case CameraModelType::FISHEYE:
        case CameraModelType::EQUISOLID:       return d != CameraDistortionType::Rational;
        case CameraModelType::EQUIRECTANGULAR: return d == CameraDistortionType::None;
        default:                               return false;
    }
}

// COLMAP / NerfStudio camera-model string -> CameraModelType.
// Returns CameraModelType(-1) for unknown / unsupported models; callers
// should validate and raise. Mirrors Camera.h:90-105.
inline CameraModelType camera_model_from_name(const std::string& name) {
    if (name == "PINHOLE" ||
        name == "SIMPLE_PINHOLE" ||
        name == "SIMPLE_RADIAL" ||
        name == "RADIAL" ||
        name == "OPENCV")             return CameraModelType::PINHOLE;
    if (name == "FISHEYE" ||
        name == "SIMPLE_FISHEYE" ||
        name == "SIMPLE_RADIAL_FISHEYE" ||
        name == "RADIAL_FISHEYE" ||
        name == "OPENCV_FISHEYE" ||
        name == "THIN_PRISM_FISHEYE") return CameraModelType::FISHEYE;
    if (name == "EQUISOLID")          return CameraModelType::EQUISOLID;
    if (name == "EQUIRECTANGULAR")    return CameraModelType::EQUIRECTANGULAR;
    return (CameraModelType)-1;
}

inline const char* camera_model_to_string(CameraModelType m) {
    switch (m) {
        case CameraModelType::PINHOLE:         return "PINHOLE";
        case CameraModelType::FISHEYE:         return "FISHEYE";
        case CameraModelType::EQUISOLID:       return "EQUISOLID";
        case CameraModelType::EQUIRECTANGULAR: return "EQUIRECTANGULAR";
        default:                               return "UNKNOWN";
    }
}
