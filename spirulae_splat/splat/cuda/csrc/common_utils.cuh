#pragma once

#include <glm/glm.hpp>

#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_fp16.h>
#include <vector_types.h>
#include <type_traits>


#ifndef SLANG_PRELUDE_EXPORT

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

#define _DEF_GENERIC_VEC_UNARY_OP(dtype, o) \
    __host__ __device__ __forceinline__ dtype##2 operator o(dtype##2 v) \
        { return { o v.x, o v.y }; } \
    __host__ __device__ __forceinline__ dtype##3 operator o(dtype##3 v) \
        { return { o v.x, o v.y, o v.z }; } \
    __host__ __device__ __forceinline__ dtype##4 operator o(dtype##4 v) \
        { return { o v.x, o v.y, o v.z, o v.w }; } \

#define _DEF_GENERIC_VEC_BINARY_OP(dtype, o) \
    __host__ __device__ __forceinline__ dtype##2 operator o(dtype##2 a, dtype##2 b) \
        { return { a.x o b.x, a.y o b.y }; } \
    __host__ __device__ __forceinline__ dtype##3 operator o(dtype##3 a, dtype##3 b) \
        { return { a.x o b.x, a.y o b.y, a.z o b.z }; } \
    __host__ __device__ __forceinline__ dtype##4 operator o(dtype##4 a, dtype##4 b) \
        { return { a.x o b.x, a.y o b.y, a.z o b.z, a.w o b.w }; } \
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

#define _DEF_BOOLEAN_VEC_BINARY_OP(dtype, oc, ox) \
    __device__ __forceinline__ bool operator oc(dtype##2 a, dtype##2 b) \
        { return { (a.x oc b.x) ox (a.y oc b.y) }; } \
    __device__ __forceinline__ bool operator oc(dtype##3 a, dtype##3 b) \
        { return { (a.x oc b.x) ox (a.y oc b.y) ox (a.z oc b.z) }; } \
    __device__ __forceinline__ bool operator oc(dtype##4 a, dtype##4 b) \
        { return { (a.x oc b.x) ox (a.y oc b.y) ox (a.z oc b.z) ox (a.w oc b.w) }; } \

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

#define _DEF_FLOAT_VEC_TO_SCALAR_FUN(dtype) \
    __device__ __forceinline__ dtype length(dtype##2 v) \
        { return { sqrtf(v.x*v.x + v.y*v.y) }; } \
    __device__ __forceinline__ dtype length(dtype##3 v) \
        { return { sqrtf(v.x*v.x + v.y*v.y + v.z*v.z) }; } \
    __device__ __forceinline__ dtype length(dtype##4 v) \
        { return { sqrtf(v.x*v.x + v.y*v.y + v.z*v.z + v.w*v.w) }; } \

#define _DEF_FLOAT_VEC_TO_VEC_FUN(dtype) \
    __device__ __forceinline__ dtype##2 normalize(dtype##2 v) \
        { return { v * rsqrtf(v.x*v.x + v.y*v.y) }; } \
    __device__ __forceinline__ dtype##3 normalize(dtype##3 v) \
        { return { v * rsqrtf(v.x*v.x + v.y*v.y + v.z*v.z) }; } \
    __device__ __forceinline__ dtype##4 normalize(dtype##4 v) \
        { return { v * rsqrtf(v.x*v.x + v.y*v.y + v.z*v.z + v.w*v.w) }; } \

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

#define _DEF_INTEGRAL_VEC_FUN(dtype) \
    _DEF_GENERIC_VEC_FUN(dtype) \
    _DEF_GENERIC_VEC_BINARY_FUN(dtype, min, min) \
    _DEF_GENERIC_VEC_BINARY_FUN(dtype, max, max) \

_DEF_FLOAT_VEC_FUN(float)
_DEF_INTEGRAL_VEC_FUN(int)
_DEF_INTEGRAL_VEC_FUN(uint)

#endif

template<typename T, int ndim>
struct TensorView {
    T* __restrict__ data;
    long shape[ndim];
    long strides[ndim];

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

    __device__ float load1f(long i) const
        { static_assert(ndim == 2 && std::is_same<T, float>::value); return at(i, 0); }
    __device__ float load1f(long i0, long i1) const
        { static_assert(ndim == 3 && std::is_same<T, float>::value); return at(i0, i1, 0); }

    __device__ float2 load2f(long i) const
        { static_assert(ndim == 2 && std::is_same<T, float>::value); return {at(i, 0), at(i, 1)}; }
    __device__ float2 load2f(long i0, long i1) const
        { static_assert(ndim == 3 && std::is_same<T, float>::value); return {at(i0, i1, 0), at(i0, i1, 1)}; }

    __device__ float3 load3f(long i) const
        { static_assert(ndim == 2 && std::is_same<T, float>::value); return {at(i, 0), at(i, 1), at(i, 2)}; }
    __device__ float3 load3f(long i0, long i1) const
        { static_assert(ndim == 3 && std::is_same<T, float>::value); return {at(i0, i1, 0), at(i0, i1, 1), at(i0, i1, 2)}; }

    __device__ float4 load4f(long i) const
        { static_assert(ndim == 2 && std::is_same<T, float>::value); return {at(i, 0), at(i, 1), at(i, 2), at(i, 3)}; }
    __device__ float4 load4f(long i0, long i1) const
        { static_assert(ndim == 3 && std::is_same<T, float>::value); return {at(i0, i1, 0), at(i0, i1, 1), at(i0, i1, 2), at(i0, i1, 3)}; }

    __device__ void store1f(long i, float v)
        { static_assert(ndim == 2 && std::is_same<T, float>::value); at(i, 0) = v; }
    __device__ void store1f(long i0, long i1, float v)
        { static_assert(ndim == 3 && std::is_same<T, float>::value); at(i0, i1, 0) = v; }

    __device__ void store2f(long i, float2 v)
        { static_assert(ndim == 2 && std::is_same<T, float>::value); at(i, 0) = v.x; at(i, 1) = v.y; }
    __device__ void store2f(long i0, long i1, float2 v)
        { static_assert(ndim == 3 && std::is_same<T, float>::value); at(i0, i1, 0) = v.x; at(i0, i1, 1) = v.y; }

    __device__ void store3f(long i, float3 v)
        { static_assert(ndim == 2 && std::is_same<T, float>::value); at(i, 0) = v.x; at(i, 1) = v.y; at(i, 2) = v.z; }
    __device__ void store3f(long i0, long i1, float3 v)
        { static_assert(ndim == 3 && std::is_same<T, float>::value); at(i0, i1, 0) = v.x; at(i0, i1, 1) = v.y; at(i0, i1, 2) = v.z; }

    __device__ void store4f(long i, float4 v)
        { static_assert(ndim == 2 && std::is_same<T, float>::value); at(i, 0) = v.x; at(i, 1) = v.y; at(i, 2) = v.z; at(i, 4) = v.w; }
    __device__ void store4f(long i0, long i1, float4 v)
        { static_assert(ndim == 3 && std::is_same<T, float>::value); at(i0, i1, 0) = v.x; at(i0, i1, 1) = v.y; at(i0, i1, 2) = v.z; at(i0, i1, 3) = v.w; }

};

