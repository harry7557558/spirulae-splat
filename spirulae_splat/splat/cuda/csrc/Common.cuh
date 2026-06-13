#pragma once

#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_fp16.h>
#include <vector_types.h>
#include <type_traits>

#include <cstdint>
#include <string>
#include <stdexcept>
#include <tuple>
#include <variant>
#include <array>
#include <cmath>

// MSVC's <cmath> doesn't expose M_PI unless _USE_MATH_DEFINES is set before
// the first include. We can't rely on that ordering in TUs, so just define
// the constants we use ourselves.
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifdef _MSC_VER
/* Old compatibility names for C types.  */
typedef unsigned long int ulong;
typedef unsigned short int ushort;
typedef unsigned int uint;
#endif

using std::max;
using std::min;

inline constexpr int WARP_SIZE = 32;

inline constexpr int TILE_SIZE = 16;
inline constexpr float ALPHA_THRESHOLD = (1.f/255.f);


#ifdef __CUDACC__

#include "generated/slang.cuh"

#else  // #ifdef __CUDACC__

template<typename T, size_t SIZE>
struct FixedArray
{
    T m_data[SIZE];
};

#endif  // #ifdef __CUDACC__


#ifdef __CUDACC__
// k1 k2 k3 k4 p1 p2 sx1 sy1 b1 b2
typedef FixedArray<float, 10> CameraDistortionCoeffs;
#endif

#include <Tensor.h>

struct CameraDistortionCoeffsBuffer {
    float* __restrict__ coeffs;

    CameraDistortionCoeffsBuffer(const TorchTensorView &tensor) {
        if (std::get<2>(tensor).size() != 2 || std::get<2>(tensor)[1] != 10)
            throw std::runtime_error("dist coeffs must have shape (C, 10)");
        if (std::get<1>(tensor) != sizeof(float))
            throw std::runtime_error("dist coeffs must be float");
        coeffs = (float*)std::get<0>(tensor);
    }

    CameraDistortionCoeffsBuffer(float* ptr) : coeffs(ptr) {}

    #ifdef __CUDACC__
    __device__ CameraDistortionCoeffs load(long idx) const {
        CameraDistortionCoeffs res;
        if (coeffs == nullptr) {
            #pragma unroll
            for (int i = 0; i < 10; i++)
                res[i] = 0.0f;
        } else {
            float* c = coeffs + 10 * idx;
            #pragma unroll
            for (int i = 0; i < 10; i++)
                res[i] = c[i];
        }
        return res;
    }
    #endif
};

// Camera Types
// This must match projection_utils.slang
enum class CameraModelType {
    PINHOLE = 0,
    FISHEYE = 1,
    EQUISOLID = 2,
    EQUIRECTANGULAR = 3,
};

inline CameraModelType cmt(const std::string &s) {
    return (s == "PINHOLE") ? CameraModelType::PINHOLE :
        (s == "FISHEYE") ? CameraModelType::FISHEYE :
        (s == "EQUISOLID") ? CameraModelType::EQUISOLID :
        (s == "EQUIRECTANGULAR") ? CameraModelType::EQUIRECTANGULAR :
        (CameraModelType)-1;
}


// Kernel launcher

#define CHECK_INPUT(x)
#define DEVICE_GUARD(_ten)

#define CHECK_DEVICE_ERROR(call)                                    \
do {                                                                \
    cudaError_t err = call;                                         \
    if (err != cudaSuccess) {                                       \
        fprintf(stderr, "\033[41mCUDA Error at %s:%d: %s\033[m\n",  \
                __FILE__, __LINE__, cudaGetErrorString(err));       \
        exit(EXIT_FAILURE);                                         \
    }                                                               \
} while (0)

#define _CEIL_DIV(n,m) (((n)+(m)-1)/(m))

#define _LAUNCH_ARGS_1D(n,b) ((n)==0?1:_CEIL_DIV(n,b)),b,0,(cudaStream_t)0
#define _LAUNCH_ARGS_2D(nx,ny,bx,by) dim3(_CEIL_DIV(nx,bx),_CEIL_DIV(ny,by),1),dim3(bx,by),0,(cudaStream_t)0
#define _LAUNCH_ARGS_3D(nx,ny,nz,bx,by,bz) dim3(_CEIL_DIV(nx,bx),_CEIL_DIV(ny,by),_CEIL_DIV(nz,bz)),dim3(bx,by,bz),0,(cudaStream_t)0

//--------------
#define CUDA_CALL(x)                                                           \
    do {                                                                       \
        if ((x) != cudaSuccess) {                                              \
            printf(                                                            \
                "Error at %s:%d - %s\n",                                       \
                __FILE__,                                                      \
                __LINE__,                                                      \
                cudaGetErrorString(cudaGetLastError())                         \
            );                                                                 \
            exit(EXIT_FAILURE);                                                \
        }                                                                      \
    } while (0)


// Routes CUB temporary storage through DeviceScratch (monotonically growing,
// never fragmented) instead of the PyTorch caching allocator.
#define CUB_WRAPPER(func, ...)                                                 \
    do {                                                                       \
        size_t _cub_temp_bytes = 0;                                            \
        func(nullptr, _cub_temp_bytes, __VA_ARGS__);                           \
        func(DeviceScratch::global().acquire(_cub_temp_bytes),                 \
             _cub_temp_bytes, __VA_ARGS__);                                    \
    } while (false)


#ifdef __CUDACC__

// Vector types

#define _DEF_GENERIC_VEC_FUNCTIONAL_UNARY_OP(dtype, o) \
    __host__ __device__ __forceinline__ dtype##2 operator o(dtype##2 v) \
        { return { o v.x, o v.y }; } \
    __host__ __device__ __forceinline__ dtype##3 operator o(dtype##3 v) \
        { return { o v.x, o v.y, o v.z }; } \
    __host__ __device__ __forceinline__ dtype##4 operator o(dtype##4 v) \
        { return { o v.x, o v.y, o v.z, o v.w }; } \

#define _DEF_GENERIC_VEC_VEC_FUNCTIONAL_BINARY_OP(dtype, o) \
    __host__ __device__ __forceinline__ dtype##2 operator o(dtype##2 a, dtype##2 b) \
        { return { a.x o b.x, a.y o b.y }; } \
    __host__ __device__ __forceinline__ dtype##3 operator o(dtype##3 a, dtype##3 b) \
        { return { a.x o b.x, a.y o b.y, a.z o b.z }; } \
    __host__ __device__ __forceinline__ dtype##4 operator o(dtype##4 a, dtype##4 b) \
        { return { a.x o b.x, a.y o b.y, a.z o b.z, a.w o b.w }; } \

#define _DEF_GENERIC_VEC_SCALAR_FUNCTIONAL_BINARY_OP(dtype, o) \
    __host__ __device__ __forceinline__ dtype##2 operator o(dtype##2 a, dtype b) \
        { return { a.x o b, a.y o b }; } \
    __host__ __device__ __forceinline__ dtype##3 operator o(dtype##3 a, dtype b) \
        { return { a.x o b, a.y o b, a.z o b }; } \
    __host__ __device__ __forceinline__ dtype##4 operator o(dtype##4 a, dtype b) \
        { return { a.x o b, a.y o b, a.z o b, a.w o b }; } \
    __host__ __device__ __forceinline__ dtype##2 operator o(dtype a, dtype##2 b) \
        { return { a o b.x, a o b.y }; } \
    __host__ __device__ __forceinline__ dtype##3 operator o(dtype a, dtype##3 b) \
        { return { a o b.x, a o b.y, a o b.z }; } \
    __host__ __device__ __forceinline__ dtype##4 operator o(dtype a, dtype##4 b) \
        { return { a o b.x, a o b.y, a o b.z, a o b.w }; } \

#define _DEF_GENERIC_VEC_INPLACE_BINARY_OP(dtype, o) \
    __host__ __device__ __forceinline__ dtype##2 operator o##=(dtype##2 &a, dtype##2 b) \
        { return { a.x o##= b.x, a.y o##= b.y }; } \
    __host__ __device__ __forceinline__ dtype##3 operator o##=(dtype##3 &a, dtype##3 b) \
        { return { a.x o##= b.x, a.y o##= b.y, a.z o##= b.z }; } \
    __host__ __device__ __forceinline__ dtype##4 operator o##=(dtype##4 &a, dtype##4 b) \
        { return { a.x o##= b.x, a.y o##= b.y, a.z o##= b.z, a.w o##= b.w }; } \

#define _DEF_BOOLEAN_VEC_FUNCTIONAL_BINARY_OP(dtype, oc, ox) \
    __device__ __forceinline__ bool operator oc(dtype##2 a, dtype##2 b) \
        { return { (a.x oc b.x) ox (a.y oc b.y) }; } \
    __device__ __forceinline__ bool operator oc(dtype##3 a, dtype##3 b) \
        { return { (a.x oc b.x) ox (a.y oc b.y) ox (a.z oc b.z) }; } \
    __device__ __forceinline__ bool operator oc(dtype##4 a, dtype##4 b) \
        { return { (a.x oc b.x) ox (a.y oc b.y) ox (a.z oc b.z) ox (a.w oc b.w) }; } \

// TODO: this may break in different Slang versions
#ifdef SLANG_PRELUDE_EXPORT

#define _DEF_GENERIC_VEC_UNARY_OP(dtype, o)

#define _DEF_GENERIC_VEC_BINARY_OP(dtype, o) \
    _DEF_GENERIC_VEC_SCALAR_FUNCTIONAL_BINARY_OP(dtype, o) \
    _DEF_GENERIC_VEC_INPLACE_BINARY_OP(dtype, o) \

#define _DEF_BOOLEAN_VEC_BINARY_OP(dtype, oc, ox)

#else  // #ifdef SLANG_PRELUDE_EXPORT

#define _DEF_GENERIC_VEC_UNARY_OP(dtype, o) \
    _DEF_GENERIC_VEC_FUNCTIONAL_UNARY_OP(dtype, o)

#define _DEF_GENERIC_VEC_BINARY_OP(dtype, o) \
    _DEF_GENERIC_VEC_VEC_FUNCTIONAL_BINARY_OP(dtype, o) \
    _DEF_GENERIC_VEC_SCALAR_FUNCTIONAL_BINARY_OP(dtype, o) \
    _DEF_GENERIC_VEC_INPLACE_BINARY_OP(dtype, o) \

#define _DEF_BOOLEAN_VEC_BINARY_OP(dtype, oc, ox) \
    _DEF_BOOLEAN_VEC_FUNCTIONAL_BINARY_OP(dtype, oc, ox)

#endif  // #ifdef SLANG_PRELUDE_EXPORT

#define _DEF_GENERIC_VEC_OP(dtype) \
    _DEF_GENERIC_VEC_UNARY_OP(dtype, -) \
    _DEF_GENERIC_VEC_BINARY_OP(dtype, +) \
    _DEF_GENERIC_VEC_BINARY_OP(dtype, -) \
    _DEF_GENERIC_VEC_BINARY_OP(dtype, *) \
    _DEF_GENERIC_VEC_BINARY_OP(dtype, /) \
    _DEF_BOOLEAN_VEC_BINARY_OP(dtype, ==, &&) \
    _DEF_BOOLEAN_VEC_BINARY_OP(dtype, !=, ||) \

#define _DEF_FLOAT_VEC_OP(dtype) \
    _DEF_GENERIC_VEC_OP(dtype) \

#define _DEF_INTEGRAL_VEC_OP(dtype) \
    _DEF_GENERIC_VEC_OP(dtype) \
    _DEF_GENERIC_VEC_BINARY_OP(dtype, %) \
    _DEF_GENERIC_VEC_BINARY_OP(dtype, &) \
    _DEF_GENERIC_VEC_BINARY_OP(dtype, |) \
    _DEF_GENERIC_VEC_BINARY_OP(dtype, ^) \
    _DEF_GENERIC_VEC_BINARY_OP(dtype, <<) \
    _DEF_GENERIC_VEC_BINARY_OP(dtype, >>) \

_DEF_FLOAT_VEC_OP(float)
_DEF_INTEGRAL_VEC_OP(int)
_DEF_INTEGRAL_VEC_OP(uint)

#ifndef SLANG_PRELUDE_EXPORT

#define _DEF_MAKE_VEC_FROM_SCALAR(dtype) \
    __host__ __device__ __forceinline__ dtype##2 make_##dtype##2(dtype v) \
        { return { v, v }; } \
    __host__ __device__ __forceinline__ dtype##3 make_##dtype##3(dtype v) \
        { return { v, v, v }; } \
    __host__ __device__ __forceinline__ dtype##4 make_##dtype##4(dtype v) \
        { return { v, v, v, v }; } \

_DEF_MAKE_VEC_FROM_SCALAR(float)
_DEF_MAKE_VEC_FROM_SCALAR(int)
_DEF_MAKE_VEC_FROM_SCALAR(uint)

#endif  // #ifndef SLANG_PRELUDE_EXPORT

#define _DEF_MAKE_VEC_FROM_VEC(dtype1, dtype2) \
    __host__ __device__ __forceinline__ dtype1##2 make_##dtype1##2(dtype2##2 v) \
        { return { (dtype1)v.x, (dtype1)v.y }; } \
    __host__ __device__ __forceinline__ dtype1##3 make_##dtype1##3(dtype2##3 v) \
        { return { (dtype1)v.x, (dtype1)v.y, (dtype1)v.z }; } \
    __host__ __device__ __forceinline__ dtype1##4 make_##dtype1##4(dtype2##4 v) \
        { return { (dtype1)v.x, (dtype1)v.y, (dtype1)v.z, (dtype1)v.w }; } \

_DEF_MAKE_VEC_FROM_VEC(float, int)
_DEF_MAKE_VEC_FROM_VEC(float, uint)
_DEF_MAKE_VEC_FROM_VEC(int, float)
_DEF_MAKE_VEC_FROM_VEC(uint, float)

#define _DEF_GENERIC_VEC_UNARY_FUN(dtype, fun, efun) \
    __device__ __forceinline__ dtype##2 fun(dtype##2 v) \
        { return { efun(v.x), efun(v.y) }; } \
    __device__ __forceinline__ dtype##3 fun(dtype##3 v) \
        { return { efun(v.x), efun(v.y), efun(v.z) }; } \
    __device__ __forceinline__ dtype##4 fun(dtype##4 v) \
        { return { efun(v.x), efun(v.y), efun(v.z), efun(v.w) }; } \

#define _DEF_GENERIC_VEC_BINARY_FUN(dtype, fun, efun) \
    __device__ __forceinline__ dtype##2 fun(dtype##2 a, dtype##2 b) \
        { return make_##dtype##2(efun(a.x, b.x), efun(a.y, b.y)); } \
    __device__ __forceinline__ dtype##3 fun(dtype##3 a, dtype##3 b) \
        { return make_##dtype##3(efun(a.x, b.x), efun(a.y, b.y), efun(a.z, b.z)); } \
    __device__ __forceinline__ dtype##4 fun(dtype##4 a, dtype##4 b) \
        { return make_##dtype##4(efun(a.x, b.x), efun(a.y, b.y), efun(a.z, b.z), efun(a.w, b.w)); } \
    __device__ __forceinline__ dtype##2 fun(dtype##2 a, dtype b) \
        { return make_##dtype##2(efun(a.x, b), efun(a.y, b)); } \
    __device__ __forceinline__ dtype##3 fun(dtype##3 a, dtype b) \
        { return make_##dtype##3(efun(a.x, b), efun(a.y, b), efun(a.z, b)); } \
    __device__ __forceinline__ dtype##4 fun(dtype##4 a, dtype b) \
        { return make_##dtype##4(efun(a.x, b), efun(a.y, b), efun(a.z, b), efun(a.w, b)); } \

#define _DEF_FLOAT_VEC_TO_SCALAR_FUN(dtype) \
    __device__ __forceinline__ dtype length(dtype##2 v) \
        { return { sqrtf(v.x*v.x + v.y*v.y) }; } \
    __device__ __forceinline__ dtype length(dtype##3 v) \
        { return { sqrtf(v.x*v.x + v.y*v.y + v.z*v.z) }; } \
    __device__ __forceinline__ dtype length(dtype##4 v) \
        { return { sqrtf(v.x*v.x + v.y*v.y + v.z*v.z + v.w*v.w) }; } \

#define _DEF_FLOAT_VEC_TO_VEC_FUN(dtype) \
    __device__ __forceinline__ dtype##2 normalize(dtype##2 v) \
        { return { v * (1.f/sqrtf(v.x*v.x + v.y*v.y)) }; } \
    __device__ __forceinline__ dtype##3 normalize(dtype##3 v) \
        { return { v * (1.f/sqrtf(v.x*v.x + v.y*v.y + v.z*v.z)) }; } \
    __device__ __forceinline__ dtype##4 normalize(dtype##4 v) \
        { return { v * (1.f/sqrtf(v.x*v.x + v.y*v.y + v.z*v.z + v.w*v.w)) }; } \

#define _DEF_GENERIC_VEC_VEC_TO_SCALAR_FUN(dtype) \
    __device__ __forceinline__ dtype dot(dtype##2 a, dtype##2 b) \
        { return { a.x * b.x + a.y * b.y }; } \
    __device__ __forceinline__ dtype dot(dtype##3 a, dtype##3 b) \
        { return { a.x * b.x + a.y * b.y + a.z * b.z }; } \
    __device__ __forceinline__ dtype dot(dtype##4 a, dtype##4 b) \
        { return { a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w }; } \

#define _DEF_GENERIC_VEC_FUN(dtype) \
    _DEF_GENERIC_VEC_VEC_TO_SCALAR_FUN(dtype) \

#define _DEF_FLOAT_VEC_FUN(dtype) \
    _DEF_GENERIC_VEC_FUN(dtype) \
    _DEF_GENERIC_VEC_UNARY_FUN(dtype, fabs, fabsf) \
    _DEF_GENERIC_VEC_BINARY_FUN(dtype, fmin, fminf) \
    _DEF_GENERIC_VEC_BINARY_FUN(dtype, fmax, fmaxf) \
    _DEF_FLOAT_VEC_TO_SCALAR_FUN(dtype) \
    _DEF_FLOAT_VEC_TO_VEC_FUN(dtype) \

#ifdef SLANG_PRELUDE_EXPORT

#define _DEF_INTEGRAL_VEC_FUN(dtype) \
    _DEF_GENERIC_VEC_FUN(dtype)

#else  // #ifdef SLANG_PRELUDE_EXPORT

#define _DEF_INTEGRAL_VEC_FUN(dtype) \
    _DEF_GENERIC_VEC_FUN(dtype) \
    _DEF_GENERIC_VEC_BINARY_FUN(dtype, min, min) \
    _DEF_GENERIC_VEC_BINARY_FUN(dtype, max, max) \

#endif  // #ifdef SLANG_PRELUDE_EXPORT

_DEF_FLOAT_VEC_FUN(float)
_DEF_INTEGRAL_VEC_FUN(int)
_DEF_INTEGRAL_VEC_FUN(uint)

__host__ __device__ __forceinline__ float3 cross(float3 a, float3 b) {
    return {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
}

#ifdef SLANG_PRELUDE_EXPORT
typedef Matrix<float, 2, 2> float2x2;
typedef Matrix<float, 3, 3> float3x3;
typedef Matrix<float, 4, 4> float4x4;
#endif


///////////////////////////////
// Reduce / Atomic
///////////////////////////////

__forceinline__ __device__ uint32_t _warpIdx() {
    // uint32_t warpId;
    // asm("mov.u32 %0, %%warpid;" : "=r"(warpId));
    // return warpId;
    return (threadIdx.x + blockDim.x * (threadIdx.y + blockDim.y * threadIdx.z)) / WARP_SIZE;
}
__forceinline__ __device__ uint32_t _laneIdx() {
    uint32_t laneId;
    asm("mov.u32 %0, %%laneid;" : "=r"(laneId));
    return laneId;
}

///////////////////////////////
// reduce from gsplat, updates values for all threads
///////////////////////////////

#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>

namespace cg = cooperative_groups;

template <uint32_t DIM, typename WarpT>
inline __device__ void warpSum(float *val, WarpT &warp) {
    #pragma unroll
    for (uint32_t i = 0; i < DIM; i++) {
        val[i] = cg::reduce(warp, val[i], cg::plus<float>());
    }
}

template <typename WarpT> inline __device__ void warpSum(float &val, WarpT &warp) {
    val = cg::reduce(warp, val, cg::plus<float>());
}

template <typename WarpT> inline __device__ void warpSum(float4 &val, WarpT &warp) {
    val.x = cg::reduce(warp, val.x, cg::plus<float>());
    val.y = cg::reduce(warp, val.y, cg::plus<float>());
    val.z = cg::reduce(warp, val.z, cg::plus<float>());
    val.w = cg::reduce(warp, val.w, cg::plus<float>());
}

template <typename WarpT> inline __device__ void warpSum(float3 &val, WarpT &warp) {
    val.x = cg::reduce(warp, val.x, cg::plus<float>());
    val.y = cg::reduce(warp, val.y, cg::plus<float>());
    val.z = cg::reduce(warp, val.z, cg::plus<float>());
}

template <typename WarpT> inline __device__ void warpSum(float2 &val, WarpT &warp) {
    val.x = cg::reduce(warp, val.x, cg::plus<float>());
    val.y = cg::reduce(warp, val.y, cg::plus<float>());
}

template <typename WarpT> inline __device__ void warpMax(float &val, WarpT &warp) {
    val = cg::reduce(warp, val, cg::greater<float>());
}

template <typename WarpT> inline __device__ void warpSum(int &val, WarpT &warp) {
    val = cg::reduce(warp, val, cg::plus<int>());
}


///////////////////////////////
// warp reduce, update value to first lane
///////////////////////////////

inline __device__ void warpSum(float &val) {
    #pragma unroll
    for (int offset = WARP_SIZE >> 1; offset > 0; offset >>= 1) {
        val += __shfl_down_sync(0xFFFFFFFF, val, offset);
    }
}

inline __device__ void warpSum(float4 &val) {
    warpSum(val.x); warpSum(val.y); warpSum(val.z); warpSum(val.w);
}

inline __device__ void warpSum(float3 &val) {
    warpSum(val.x); warpSum(val.y); warpSum(val.z);
}

inline __device__ void warpSum(float2 &val) {
    warpSum(val.x); warpSum(val.y);
}

inline __device__ void warpAtomicAdd(float* addr, float val) {
    #pragma unroll
    for (int offset = WARP_SIZE >> 1; offset > 0; offset >>= 1) {
        val += __shfl_down_sync(~0u, val, offset);
    }
    if (_laneIdx() == 0 && val != 0.0f)
        atomicAdd(addr, val);
}

inline __device__ void warpAtomicAdd(int* addr, int val) {
    #pragma unroll
    for (int offset = WARP_SIZE >> 1; offset > 0; offset >>= 1) {
        val += __shfl_down_sync(~0u, val, offset);
    }
    if (_laneIdx() == 0 && val != 0)
        atomicAdd(addr, val);
}

inline __device__ void warpAtomicAdd(float4* addr, float4 val) {
#if 0
    warpAtomicAdd(&addr->x, val.x);
    warpAtomicAdd(&addr->y, val.y);
    warpAtomicAdd(&addr->z, val.z);
    warpAtomicAdd(&addr->w, val.w);
#else
    #pragma unroll
    for (int offset = WARP_SIZE >> 1; offset > 0; offset >>= 1) {
        val.x += __shfl_down_sync(~0u, val.x, offset);
        val.y += __shfl_down_sync(~0u, val.y, offset);
        val.z += __shfl_down_sync(~0u, val.z, offset);
        val.w += __shfl_down_sync(~0u, val.w, offset);
    }
    if (_laneIdx() == 0) {
        if (val.x) atomicAdd(&addr->x, val.x);
        if (val.y) atomicAdd(&addr->y, val.y);
        if (val.z) atomicAdd(&addr->z, val.z);
        if (val.w) atomicAdd(&addr->w, val.w);
    }
#endif
}

inline __device__ void warpAtomicAdd(float3* addr, float3 val) {
#if 0
    warpAtomicAdd(&addr->x, val.x);
    warpAtomicAdd(&addr->y, val.y);
    warpAtomicAdd(&addr->z, val.z);
#else
    #pragma unroll
    for (int offset = WARP_SIZE >> 1; offset > 0; offset >>= 1) {
        val.x += __shfl_down_sync(~0u, val.x, offset);
        val.y += __shfl_down_sync(~0u, val.y, offset);
        val.z += __shfl_down_sync(~0u, val.z, offset);
    }
    if (_laneIdx() == 0) {
        if (val.x) atomicAdd(&addr->x, val.x);
        if (val.y) atomicAdd(&addr->y, val.y);
        if (val.z) atomicAdd(&addr->z, val.z);
    }
#endif
}

inline __device__ void warpAtomicAdd(float2* addr, float2 &val) {
#if 0
    warpAtomicAdd(&addr->x, val.x);
    warpAtomicAdd(&addr->y, val.y);
#else
    #pragma unroll
    for (int offset = WARP_SIZE >> 1; offset > 0; offset >>= 1) {
        val.x += __shfl_down_sync(~0u, val.x, offset);
        val.y += __shfl_down_sync(~0u, val.y, offset);
    }
    if (_laneIdx() == 0) {
        if (val.x) atomicAdd(&addr->x, val.x);
        if (val.y) atomicAdd(&addr->y, val.y);
    }
#endif
}


///////////////////////////////
// block reduce / atomic
///////////////////////////////

template<int BLOCK_SIZE>
inline __device__ void blockAtomicAdd(float* addr, float val) {
    static_assert(BLOCK_SIZE > WARP_SIZE && BLOCK_SIZE <= WARP_SIZE * WARP_SIZE);
    static_assert(BLOCK_SIZE % WARP_SIZE == 0);

    static __shared__ float sharedSums[BLOCK_SIZE / WARP_SIZE];

    uint laneId = _laneIdx();
    uint warpId = _warpIdx();

    #pragma unroll
    for (int stride = WARP_SIZE >> 1; stride > 0; stride >>= 1) {
        val += __shfl_down_sync(0xFFFFFFFF, val, stride);
    }

    if (laneId == 0)
        sharedSums[warpId] = val;
    __syncthreads();

    if (warpId == 0) {
        val = laneId < BLOCK_SIZE / WARP_SIZE ?
            sharedSums[laneId] : 0.0f;
        #pragma unroll
        for (int stride = (BLOCK_SIZE / WARP_SIZE) >> 1; stride > 0; stride >>= 1)
            val += __shfl_down_sync(0xFFFFFFFF, val, stride);
        if (laneId == 0 && val != 0.0f)
            atomicAdd(addr, val);
    }
}

template<int BLOCK_SIZE>
inline __device__ void blockAtomicAdd(float4* addr, float4 val) {
    static_assert(BLOCK_SIZE >= 4 * WARP_SIZE && BLOCK_SIZE <= WARP_SIZE * WARP_SIZE);
    static_assert(BLOCK_SIZE % WARP_SIZE == 0);

    constexpr int NUM_WARPS = BLOCK_SIZE / WARP_SIZE;
    static __shared__ float sharedSums[4 * NUM_WARPS];

    uint laneId = _laneIdx();
    uint warpId = _warpIdx();

    #pragma unroll
    for (int stride = WARP_SIZE >> 1; stride > 0; stride >>= 1) {
        val.x += __shfl_down_sync(0xFFFFFFFF, val.x, stride);
        val.y += __shfl_down_sync(0xFFFFFFFF, val.y, stride);
        val.z += __shfl_down_sync(0xFFFFFFFF, val.z, stride);
        val.w += __shfl_down_sync(0xFFFFFFFF, val.w, stride);
    }

    if (laneId == 0) {
        sharedSums[warpId + 0 * NUM_WARPS] = val.x;
        sharedSums[warpId + 1 * NUM_WARPS] = val.y;
        sharedSums[warpId + 2 * NUM_WARPS] = val.z;
        sharedSums[warpId + 3 * NUM_WARPS] = val.w;
    }
    __syncthreads();

    if (warpId < 4) {
        float val = laneId < NUM_WARPS ?
            sharedSums[laneId + warpId * NUM_WARPS] : 0.0f;
        #pragma unroll
        for (int stride = NUM_WARPS >> 1; stride > 0; stride >>= 1)
            val += __shfl_down_sync(0xFFFFFFFF, val, stride);
        if (laneId == 0 && val != 0.0f)
            atomicAdd((float*)addr + warpId, val);
    }
}

template<int BLOCK_SIZE>
inline __device__ void blockAtomicAdd(float3* addr, float3 val) {
    static_assert(BLOCK_SIZE >= 3 * WARP_SIZE && BLOCK_SIZE <= WARP_SIZE * WARP_SIZE);
    static_assert(BLOCK_SIZE % WARP_SIZE == 0);

    constexpr int NUM_WARPS = BLOCK_SIZE / WARP_SIZE;
    static __shared__ float sharedSums[3 * NUM_WARPS];

    uint laneId = _laneIdx();
    uint warpId = _warpIdx();

    #pragma unroll
    for (int stride = WARP_SIZE >> 1; stride > 0; stride >>= 1) {
        val.x += __shfl_down_sync(0xFFFFFFFF, val.x, stride);
        val.y += __shfl_down_sync(0xFFFFFFFF, val.y, stride);
        val.z += __shfl_down_sync(0xFFFFFFFF, val.z, stride);
    }

    if (laneId == 0) {
        sharedSums[warpId + 0 * NUM_WARPS] = val.x;
        sharedSums[warpId + 1 * NUM_WARPS] = val.y;
        sharedSums[warpId + 2 * NUM_WARPS] = val.z;
    }
    __syncthreads();

    if (warpId < 3) {
        float val = laneId < NUM_WARPS ?
            sharedSums[laneId + warpId * NUM_WARPS] : 0.0f;
        #pragma unroll
        for (int stride = NUM_WARPS >> 1; stride > 0; stride >>= 1)
            val += __shfl_down_sync(0xFFFFFFFF, val, stride);
        if (laneId == 0 && val != 0.0f)
            atomicAdd((float*)addr + warpId, val);
    }
}

template<int BLOCK_SIZE>
inline __device__ void blockAtomicAdd(float2* addr, float2 val) {
    static_assert(BLOCK_SIZE >= 2 * WARP_SIZE && BLOCK_SIZE <= WARP_SIZE * WARP_SIZE);
    static_assert(BLOCK_SIZE % WARP_SIZE == 0);

    constexpr int NUM_WARPS = BLOCK_SIZE / WARP_SIZE;
    static __shared__ float sharedSums[2 * NUM_WARPS];

    uint laneId = _laneIdx();
    uint warpId = _warpIdx();

    #pragma unroll
    for (int stride = WARP_SIZE >> 1; stride > 0; stride >>= 1) {
        val.x += __shfl_down_sync(0xFFFFFFFF, val.x, stride);
        val.y += __shfl_down_sync(0xFFFFFFFF, val.y, stride);
    }

    if (laneId == 0) {
        sharedSums[warpId + 0 * NUM_WARPS] = val.x;
        sharedSums[warpId + 1 * NUM_WARPS] = val.y;
    }
    __syncthreads();

    if (warpId < 2) {
        float val = laneId < NUM_WARPS ?
            sharedSums[laneId + warpId * NUM_WARPS] : 0.0f;
        #pragma unroll
        for (int stride = NUM_WARPS >> 1; stride > 0; stride >>= 1)
            val += __shfl_down_sync(0xFFFFFFFF, val, stride);
        if (laneId == 0 && val != 0.0f)
            atomicAdd((float*)addr + warpId, val);
    }
}


///////////////////////////////
// templated reduce / atomic
///////////////////////////////

inline __device__ float atomicMin(float* p, float v) {
    return (__float_as_int(v) >= 0) ?
        __int_as_float(atomicMin((int*)p, __float_as_int(v))) :
        __uint_as_float(atomicMax((unsigned*)p, __float_as_uint(v)));
}
inline __device__ float atomicMax(float* p, float v) {
    return (__float_as_int(v) >= 0) ?
        __int_as_float(atomicMax((int*)p, __float_as_int(v))) :
        __uint_as_float(atomicMin((unsigned*)p, __float_as_uint(v)));
}

template<int reduce = 1>
inline __device__ void atomicAddFVec(float* p, float v) {
    static_assert(reduce == 1 || reduce % WARP_SIZE == 0);
    if (p == nullptr) return;
    v = isfinite(v) ? v : 0.0f;
    if (reduce == 1)
        { if (v != 0.0f) atomicAdd(p, v); }
    else if (reduce == WARP_SIZE)
        warpAtomicAdd(p, v);
    else
        blockAtomicAdd<reduce == 1 || reduce == WARP_SIZE ? 2 * WARP_SIZE : reduce>(p, v);
}

template<int reduce = 1>
inline __device__ void atomicAddFVec(float2* p, float2 v) {
    static_assert(reduce == 1 || reduce == WARP_SIZE || reduce >= 2 * WARP_SIZE);
    if (p == nullptr) return;
    v.x = isfinite(v.x) ? v.x : 0.0f;
    v.y = isfinite(v.y) ? v.y : 0.0f;
    if (reduce == 1) {
        if (v.x != 0.0f) atomicAdd(&p->x, v.x);
        if (v.y != 0.0f) atomicAdd(&p->y, v.y);
    }
    else if (reduce == WARP_SIZE)
        warpAtomicAdd(p, v);
    else
        blockAtomicAdd<reduce == 1 || reduce == WARP_SIZE ? 2 * WARP_SIZE : reduce>(p, v);
}

template<int reduce = 1>
inline __device__ void atomicAddFVec(float3* p, float3 v) {
    static_assert(reduce == 1 || reduce % WARP_SIZE == 0);
    if (p == nullptr) return;
    v.x = isfinite(v.x) ? v.x : 0.0f;
    v.y = isfinite(v.y) ? v.y : 0.0f;
    v.z = isfinite(v.z) ? v.z : 0.0f;
    if (reduce == 1) {
        if (v.x != 0.0f) atomicAdd(&p->x, v.x);
        if (v.y != 0.0f) atomicAdd(&p->y, v.y);
        if (v.z != 0.0f) atomicAdd(&p->z, v.z);
    }
    else if (reduce == WARP_SIZE)
        warpAtomicAdd(p, v);
    else
        blockAtomicAdd<reduce == 1 || reduce == WARP_SIZE ? 3 * WARP_SIZE : reduce>(p, v);
}

template<int reduce = 1>
inline __device__ void atomicAddFVec(float4* p, float4 v) {
    static_assert(reduce == 1 || reduce % WARP_SIZE == 0);
    if (p == nullptr) return;
    v.x = isfinite(v.x) ? v.x : 0.0f;
    v.y = isfinite(v.y) ? v.y : 0.0f;
    v.z = isfinite(v.z) ? v.z : 0.0f;
    v.w = isfinite(v.w) ? v.w : 0.0f;
    if (reduce == 1) {
        if (v.x != 0.0f) atomicAdd(&p->x, v.x);
        if (v.y != 0.0f) atomicAdd(&p->y, v.y);
        if (v.z != 0.0f) atomicAdd(&p->z, v.z);
        if (v.w != 0.0f) atomicAdd(&p->w, v.w);
    }
    else if (reduce == WARP_SIZE)
        warpAtomicAdd(p, v);
    else
        blockAtomicAdd<reduce == 1 || reduce == WARP_SIZE ? 4 * WARP_SIZE : reduce>(p, v);
}


// compute a*x*a where a is tiny number and x is large number, without underflowing floating points
__forceinline__ __device__ float fmul_axa(float a, float x) {
    return __fmul_rn(__fmul_rn(a, x), a);
}
__forceinline__ __device__ float2 fmul_axa(float2 a, float x) {
    return {
        fmul_axa(a.x, x),
        fmul_axa(a.y, x)
    };
}
__forceinline__ __device__ float3 fmul_axa(float3 a, float x) {
    return {
        fmul_axa(a.x, x),
        fmul_axa(a.y, x),
        fmul_axa(a.z, x)
    };
}
__forceinline__ __device__ float4 fmul_axa(float4 a, float x) {
    return {
        fmul_axa(a.x, x),
        fmul_axa(a.y, x),
        fmul_axa(a.z, x),
        fmul_axa(a.w, x)
    };
}

#endif  // #ifdef __CUDACC__



///////////////////////////////
// Non-Contiguous Tensor
///////////////////////////////


template<typename T, int ndim>
struct TensorView {
    T* __restrict__ data;
    long shape[ndim];
    long strides[ndim];

    using vec2_t =
        typename std::conditional<std::is_same<T, float>::value, float2,
        typename std::conditional<std::is_same<T, int>::value, int2,
        typename std::conditional<std::is_same<T, uint>::value, uint2,
        typename std::conditional<std::is_same<T, uint8_t>::value, uchar2,
        typename std::conditional<std::is_same<T, uint16_t>::value, ushort2,
        void>::type>::type>::type>::type>::type;
    using vec3_t =
        typename std::conditional<std::is_same<T, float>::value, float3,
        typename std::conditional<std::is_same<T, int>::value, int3,
        typename std::conditional<std::is_same<T, uint>::value, uint3,
        typename std::conditional<std::is_same<T, uint8_t>::value, uchar3,
        typename std::conditional<std::is_same<T, uint16_t>::value, ushort3,
        void>::type>::type>::type>::type>::type;
    using vec4_t =
        typename std::conditional<std::is_same<T, float>::value, float4,
        typename std::conditional<std::is_same<T, int>::value, int4,
        typename std::conditional<std::is_same<T, uint>::value, uint4,
        typename std::conditional<std::is_same<T, uint8_t>::value, uchar4,
        typename std::conditional<std::is_same<T, uint16_t>::value, ushort4,
        void>::type>::type>::type>::type>::type;

    __device__ T at(long i) const
        { static_assert(ndim == 1); return data[i * strides[0]]; }
    __device__ T& at(long i)
        { static_assert(ndim == 1); return data[i * strides[0]]; }

    __device__ T at(long i0, long i1) const
        { static_assert(ndim == 2); return data[i0 * strides[0] + i1 * strides[1]]; }
    __device__ T& at(long i0, long i1)
        { static_assert(ndim == 2); return data[i0 * strides[0] + i1 * strides[1]]; }

    __device__ T at(long i0, long i1, long i2) const
        { static_assert(ndim == 3); return data[i0 * strides[0] + i1 * strides[1] + i2 * strides[2]]; }
    __device__ T& at(long i0, long i1, long i2)
        { static_assert(ndim == 3); return data[i0 * strides[0] + i1 * strides[1] + i2 * strides[2]]; }

    __device__ T at(long i0, long i1, long i2, long i3) const
        { static_assert(ndim == 4); return data[i0 * strides[0] + i1 * strides[1] + i2 * strides[2] + i3 * strides[3]]; }
    __device__ T& at(long i0, long i1, long i2, long i3)
        { static_assert(ndim == 4); return data[i0 * strides[0] + i1 * strides[1] + i2 * strides[2] + i3 * strides[3]]; }

    __device__ T at(long i0, long i1, long i2, long i3, long i4) const
        { static_assert(ndim == 5); return data[i0 * strides[0] + i1 * strides[1] + i2 * strides[2] + i3 * strides[3] + i4 * strides[4]]; }
    __device__ T& at(long i0, long i1, long i2, long i3, long i4)
        { static_assert(ndim == 5); return data[i0 * strides[0] + i1 * strides[1] + i2 * strides[2] + i3 * strides[3] + i4 * strides[4]]; }

    __device__ T load1(long i) const
        { static_assert(ndim == 2); return at(i, 0); }
    __device__ T load1(long i0, long i1) const
        { static_assert(ndim == 3); return at(i0, i1, 0); }
    __device__ T load1(long i0, long i1, long i2) const
        { static_assert(ndim == 4); return at(i0, i1, i2, 0); }

    __device__ void store1(long i, T v)
        { static_assert(ndim == 2); at(i, 0) = v; }
    __device__ void store1(long i0, long i1, T v)
        { static_assert(ndim == 3); at(i0, i1, 0) = v; }
    __device__ void store1(long i0, long i1, long i2, T v)
        { static_assert(ndim == 4); at(i0, i1, i2, 0) = v; }

    __device__ void atomicStore1(long i, T v)
        { static_assert(ndim == 2); if (v) atomicAdd(&at(i, 0), v); }
    __device__ void atomicStore1(long i0, long i1, T v)
        { static_assert(ndim == 3); if (v) atomicAdd(&at(i0, i1, 0), v); }
    __device__ void atomicStore1(long i0, long i1, long i2, T v)
        { static_assert(ndim == 4); if (v) atomicAdd(&at(i0, i1, i2, 0), v); }

    __device__ vec2_t load2(long i) const
        { static_assert(ndim == 2); return {at(i, 0), at(i, 1)}; }
    __device__ vec2_t load2(long i0, long i1) const
        { static_assert(ndim == 3); return {at(i0, i1, 0), at(i0, i1, 1)}; }
    __device__ vec2_t load2(long i0, long i1, long i2) const
        { static_assert(ndim == 4); return {at(i0, i1, i2, 0), at(i0, i1, i2, 1), at(i0, i1, i2, 2)}; }

    __device__ vec3_t load3(long i) const
        { static_assert(ndim == 2); return {at(i, 0), at(i, 1), at(i, 2)}; }
    __device__ vec3_t load3(long i0, long i1) const
        { static_assert(ndim == 3); return {at(i0, i1, 0), at(i0, i1, 1), at(i0, i1, 2)}; }
    __device__ vec3_t load3(long i0, long i1, long i2) const
        { static_assert(ndim == 4); return {at(i0, i1, i2, 0), at(i0, i1, i2, 1), at(i0, i1, i2, 2)}; }

    __device__ vec4_t load4(long i) const
        { static_assert(ndim == 2); return {at(i, 0), at(i, 1), at(i, 2), at(i, 3)}; }
    __device__ vec4_t load4(long i0, long i1) const
        { static_assert(ndim == 3); return {at(i0, i1, 0), at(i0, i1, 1), at(i0, i1, 2), at(i0, i1, 3)}; }
    __device__ vec4_t load4(long i0, long i1, long i2) const
        { static_assert(ndim == 4); return {at(i0, i1, i2, 0), at(i0, i1, i2, 1), at(i0, i1, i2, 2), at(i0, i1, i2, 3)}; }

    __device__ void store2(long i, vec2_t v)
        { static_assert(ndim == 2); at(i, 0) = v.x; at(i, 1) = v.y; }
    __device__ void store2(long i0, long i1, vec2_t v)
        { static_assert(ndim == 3); at(i0, i1, 0) = v.x; at(i0, i1, 1) = v.y; }
    __device__ void store2(long i0, long i1, long i2, vec2_t v)
        { static_assert(ndim == 4); at(i0, i1, i2, 0) = v.x; at(i0, i1, i2, 1) = v.y; }

    __device__ void store3(long i, vec3_t v)
        { static_assert(ndim == 2); at(i, 0) = v.x; at(i, 1) = v.y; at(i, 2) = v.z; }
    __device__ void store3(long i0, long i1, vec3_t v)
        { static_assert(ndim == 3); at(i0, i1, 0) = v.x; at(i0, i1, 1) = v.y; at(i0, i1, 2) = v.z; }
    __device__ void store3(long i0, long i1, long i2, vec3_t v)
        { static_assert(ndim == 4); at(i0, i1, i2, 0) = v.x; at(i0, i1, i2, 1) = v.y; at(i0, i1, i2, 2) = v.z; }

    __device__ void store4(long i, vec4_t v)
        { static_assert(ndim == 2); at(i, 0) = v.x; at(i, 1) = v.y; at(i, 2) = v.z; at(i, 4) = v.w; }
    __device__ void store4(long i0, long i1, vec4_t v)
        { static_assert(ndim == 3); at(i0, i1, 0) = v.x; at(i0, i1, 1) = v.y; at(i0, i1, 2) = v.z; at(i0, i1, 3) = v.w; }
    __device__ void store4(long i0, long i1, long i2, vec4_t v)
        { static_assert(ndim == 4); at(i0, i1, i2, 0) = v.x; at(i0, i1, i2, 1) = v.y; at(i0, i1, i2, 2) = v.z; at(i0, i1, i2, 3) = v.w; }

};


#ifdef __CUDACC__

__forceinline__ __device__ uint32_t hash_uint3(uint32_t a, uint32_t b, uint32_t c) {
    uint32_t hash = a;

    hash *= 0x01000193;
    hash = (hash << 16) | (hash >> 16);

    hash ^= b;
    hash *= 0x01000193;
    hash = (hash << 16) | (hash >> 16);

    hash ^= c;
    hash *= 0x01000193;

    hash ^= hash >> 15;
    hash *= 0x85ebca6b;
    hash ^= hash >> 13;
    hash *= 0xc2b2ae35;
    hash ^= hash >> 16;

    return hash;
}

__forceinline__ __device__ float sigmoid(float x) {
    return 1.0f / (1.0f + __expf(-x));
}

__forceinline__ __device__ float logit(float x) {
    return __logf(x / (1.0f - x));
}

#endif  // #ifdef __CUDACC__
