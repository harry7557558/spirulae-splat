#include <vector_types.h>
#include <type_traits>

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

