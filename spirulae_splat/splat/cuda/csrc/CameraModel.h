#pragma once

// CUDA-free camera-model enum + name helpers.
//
// CameraModelType and camera_model_from_name were historically defined in
// Common.cuh / Camera.h, which under the CUDA backend transitively pull in
// <cuda_runtime.h> (via backend/api/BackendTypes.h). The
// standalone dataset parsers (app/ColmapParser.cpp, app/NerfstudioParser.cpp,
// app/MetashapeParser.cpp, app/DatasetCommon.cpp) only need these two trivial
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
