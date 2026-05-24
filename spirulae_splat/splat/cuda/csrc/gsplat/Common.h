#pragma once

#include <algorithm>
#include <cstdint>
#include <glm/gtc/type_ptr.hpp>

#include <Tensor.h>

namespace ssplat {

// Routes CUB temporary storage through DeviceScratch (monotonically growing,
// never fragmented) instead of the PyTorch caching allocator.
#define CUB_WRAPPER(func, ...)                                                 \
    do {                                                                       \
        size_t _cub_temp_bytes = 0;                                            \
        func(nullptr, _cub_temp_bytes, __VA_ARGS__);                           \
        func(DeviceScratch::global().acquire(_cub_temp_bytes),                 \
             _cub_temp_bytes, __VA_ARGS__);                                    \
    } while (false)

//
// Convenience typedefs for CUDA types
//
using vec2 = glm::vec<2, float>;
using vec3 = glm::vec<3, float>;
using vec4 = glm::vec<4, float>;
using mat2 = glm::mat<2, 2, float>;
using mat3 = glm::mat<3, 3, float>;
using mat4 = glm::mat<4, 4, float>;
using mat3x2 = glm::mat<3, 2, float>;

//
// Legacy Camera Types
//
// This must match projection_utils.slang
enum class CameraModelType {
    PINHOLE = 0,
    FISHEYE = 1,
    EQUISOLID = 2,
    EQUIRECTANGULAR = 3,
};

// #define N_THREADS_PACKED 256
// #define ALPHA_THRESHOLD (1.f / 255.f)

} // namespace ssplat

#include <string>

inline ssplat::CameraModelType cmt(const std::string &s) {
    return (s == "PINHOLE") ? ssplat::CameraModelType::PINHOLE :
        (s == "FISHEYE") ? ssplat::CameraModelType::FISHEYE :
        (s == "EQUISOLID") ? ssplat::CameraModelType::EQUISOLID :
        (ssplat::CameraModelType)-1;
}
