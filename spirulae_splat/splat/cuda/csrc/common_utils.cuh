#pragma once

#include <glm/glm.hpp>

#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_fp16.h>
#include <vector_types.h>
#include <type_traits>


#ifdef __CUDACC__

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

inline __device__ void atomicAddFVec(float2* p, float2 v) {
    if (v.x != 0.0f) atomicAdd(&p->x, v.x);
    if (v.y != 0.0f) atomicAdd(&p->y, v.y);
}
inline __device__ void atomicAddFVec(float3* p, float3 v) {
    if (v.x != 0.0f) atomicAdd(&p->x, v.x);
    if (v.y != 0.0f) atomicAdd(&p->y, v.y);
    if (v.z != 0.0f) atomicAdd(&p->z, v.z);
}
inline __device__ void atomicAddFVec(float4* p, float4 v) {
    if (v.x != 0.0f) atomicAdd(&p->x, v.x);
    if (v.y != 0.0f) atomicAdd(&p->y, v.y);
    if (v.z != 0.0f) atomicAdd(&p->z, v.z);
    if (v.w != 0.0f) atomicAdd(&p->w, v.w);
}

#endif  // #ifdef __CUDACC__


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

#endif  // #ifndef SLANG_PRELUDE_EXPORT

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


#ifdef __CUDACC__
#ifdef SLANG_PRELUDE_EXPORT
typedef Matrix<float, 2, 2> float2x2;
typedef Matrix<float, 3, 3> float3x3;
typedef Matrix<float, 4, 4> float4x4;
#endif
#endif


///////////////////////////////
// Reduce (from gsplat)
///////////////////////////////

#ifndef SSPLAT_HOST_ONLY

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

template <typename WarpT> inline __device__ void warpSum(glm::vec4 &val, WarpT &warp) {
    val.x = cg::reduce(warp, val.x, cg::plus<float>());
    val.y = cg::reduce(warp, val.y, cg::plus<float>());
    val.z = cg::reduce(warp, val.z, cg::plus<float>());
    val.w = cg::reduce(warp, val.w, cg::plus<float>());
}

template <typename WarpT> inline __device__ void warpSum(glm::vec3 &val, WarpT &warp) {
    val.x = cg::reduce(warp, val.x, cg::plus<float>());
    val.y = cg::reduce(warp, val.y, cg::plus<float>());
    val.z = cg::reduce(warp, val.z, cg::plus<float>());
}

template <typename WarpT> inline __device__ void warpSum(glm::vec2 &val, WarpT &warp) {
    val.x = cg::reduce(warp, val.x, cg::plus<float>());
    val.y = cg::reduce(warp, val.y, cg::plus<float>());
}

template <typename WarpT> inline __device__ void warpSum(glm::mat4 &val, WarpT &warp) {
    warpSum(val[0], warp);
    warpSum(val[1], warp);
    warpSum(val[2], warp);
    warpSum(val[3], warp);
}

template <typename WarpT> inline __device__ void warpSum(glm::mat3 &val, WarpT &warp) {
    warpSum(val[0], warp);
    warpSum(val[1], warp);
    warpSum(val[2], warp);
}

template <typename WarpT> inline __device__ void warpSum(glm::mat2 &val, WarpT &warp) {
    warpSum(val[0], warp);
    warpSum(val[1], warp);
}

template <typename WarpT> inline __device__ void warpMax(float &val, WarpT &warp) {
    val = cg::reduce(warp, val, cg::greater<float>());
}

template <typename WarpT> inline __device__ void warpSum(int &val, WarpT &warp) {
    val = cg::reduce(warp, val, cg::plus<int>());
}

#endif  // #ifndef SSPLAT_HOST_ONLY



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
        void>::type>::type>::type;
    using vec3_t =
        typename std::conditional<std::is_same<T, float>::value, float3,
        typename std::conditional<std::is_same<T, int>::value, int3,
        typename std::conditional<std::is_same<T, uint>::value, uint3,
        void>::type>::type>::type;
    using vec4_t =
        typename std::conditional<std::is_same<T, float>::value, float4,
        typename std::conditional<std::is_same<T, int>::value, int4,
        typename std::conditional<std::is_same<T, uint>::value, uint4,
        void>::type>::type>::type;

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

