#pragma once

#ifdef __CUDACC__
#include "generated/slang.cuh"
namespace SlangProjectionUtils {
#include "generated/set_namespace.cuh"
#include "generated/projection_utils.cuh"
}
#endif

#include "Primitive.cuh"


template<int _sh_degree>
struct _BasePrimitive3DGS {
    class WorldBuffer;
    class ScreenBuffer;

    #ifdef __CUDACC__
    struct World {
        World& operator=(const World&) = default;

        float3 mean;
        float4 quat;
        float3 scale;  // log
        float opacity;  // logit
        float3 features_dc;
        float* __restrict__ features_sh;  // don't load everything into register

        static constexpr int num_sh() {
            return _sh_degree * (_sh_degree + 2);
        }

        static __device__ __forceinline__ World zero() {
            World w;
            w.mean = make_float3(0.f);
            w.quat = make_float4(0.f);
            w.scale = make_float3(0.f);
            w.opacity = 0.f;
            w.features_dc = make_float3(0.f);
            w.features_sh = nullptr;
            return w;
        }

        static __device__ __forceinline__ World zero(WorldBuffer &buffer, int64_t i) {
            World w;
            w.mean = make_float3(0.f);
            w.quat = make_float4(0.f);
            w.scale = make_float3(0.f);
            w.opacity = 0.f;
            w.features_dc = make_float3(0.f);
            w.features_sh = buffer.features_sh(i);
            return w;
        }

        __device__ __forceinline__ void load(WorldBuffer &buffer, int64_t i) {
            mean = buffer.means(i);
            quat = buffer.quats(i);
            scale = buffer.scales(i);
            opacity = buffer.opacities(i);
            features_dc = buffer.features_dc(i);
            features_sh = buffer.features_sh(i);
        }

        __device__ __forceinline__ void atomicLoad(WorldBuffer &buffer, int64_t i) {
            if (&buffer.means(0)) mean += buffer.means(i);
            if (&buffer.quats(0)) quat += buffer.quats(i);
            if (&buffer.scales(0)) scale += buffer.scales(i);
            if (&buffer.opacities(0)) opacity += buffer.opacities(i);
            if (&buffer.features_dc(0)) features_dc += buffer.features_dc(i);
            // This function is used in fused projection backwards and optimizer
            // features_sh must be handled by calling code manually
        }

        __device__ __forceinline__ void atomicStore(WorldBuffer &buffer, int64_t i) const {
            if (&buffer.means(0)) atomicAddFVec(&buffer.means(i), mean);
            if (&buffer.quats(0)) atomicAddFVec(&buffer.quats(i), quat);
            if (&buffer.scales(0)) atomicAddFVec(&buffer.scales(i), scale);
            if (&buffer.opacities(0)) atomicAddFVec(&buffer.opacities(i), opacity);
            if (&buffer.features_dc(0)) atomicAddFVec(&buffer.features_dc(i), features_dc);
            // This function is used in projection backward
            // features_sh is properly handled by projection function in general
        }
    };
    #endif

    class WorldBuffer : public TensorArray<6> {
    public:
        using TensorArray<6>::TensorArray;
    #ifdef __CUDACC__
        __forceinline__ __device__ float3& means(int64_t i)
            { return *reinterpret_cast<float3*>(&_data[0][3*i]); }
        __forceinline__ __device__ float4& quats(int64_t i)
            { return *reinterpret_cast<float4*>(&_data[1][4*i]); }
        __forceinline__ __device__ float3& scales(int64_t i)
            { return *reinterpret_cast<float3*>(&_data[2][3*i]); }
        __forceinline__ __device__ float& opacities(int64_t i)
            { return *reinterpret_cast<float*>(&_data[3][i]); }
        __forceinline__ __device__ float3& features_dc(int64_t i)
            { return *reinterpret_cast<float3*>(&_data[4][3*i]); }
        __forceinline__ __device__ float* features_sh(int64_t i)
            { return &_data[5][_strides[5]*i]; }

        __forceinline__ __device__ float3 means(int64_t i) const
            { return *reinterpret_cast<float3*>(&_data[0][3*i]); }
        __forceinline__ __device__ float4 quats(int64_t i) const
            { return *reinterpret_cast<float4*>(&_data[1][4*i]); }
        __forceinline__ __device__ float3 scales(int64_t i) const
            { return *reinterpret_cast<float3*>(&_data[2][3*i]); }
        __forceinline__ __device__ float opacities(int64_t i) const
            { return *reinterpret_cast<float*>(&_data[3][i]); }
        __forceinline__ __device__ float3 features_dc(int64_t i) const
            { return *reinterpret_cast<float3*>(&_data[4][3*i]); }
        __forceinline__ __device__ const float* features_sh(int64_t i) const
            { return &_data[5][_strides[5]*i]; }

    #endif  // #ifdef __CUDACC__

        int sh_degree() {
            int num_sh = _strides[5] / 3;
            return (
                num_sh == 3 ? 1 :
                num_sh == 8 ? 2 :
                num_sh == 15 ? 3 :
                num_sh == 24 ? 4 :
                0
            );
        }
        static std::vector<DeviceTensorFloatND> empty_pool(int64_t size, int num_sh, PoolSlot key_prefix) {
            return TensorArray<6>::empty_pool(size, {3, 4, 3, 1, 3, (int32_t)(3*num_sh)}, key_prefix);
        }
        static std::vector<DeviceTensorFloatND> zeros_pool(int64_t size, int num_sh, PoolSlot key_prefix) {
            return TensorArray<6>::zeros_pool(size, {3, 4, 3, 1, 3, (int32_t)(3*num_sh)}, key_prefix);
        }
        static std::vector<DeviceTensorFloatND> zeros_pool(
            const std::vector<DeviceTensorFloatND>& tmpl, PoolSlot key_prefix) {
            return TensorArray<6>::zeros_pool(tmpl, key_prefix);
        }
    };

    #ifdef __CUDACC__
    struct Screen {
        Screen& operator=(const Screen&) = default;

        float2 xy;
        float depth;
        float3 conic;
        float opac;
        float3 rgb;

        static __device__ __forceinline__ Screen zero() {
            Screen s;
            s.xy = make_float2(0.f);
            s.depth = 0.f;
            s.conic = make_float3(0.f);
            s.opac = 0.f;
            s.rgb = make_float3(0.f);
            return s;
        }

        __device__ __forceinline__ void load(const ScreenBuffer &buffer, int64_t i) {
            xy = buffer.xy(i);
            depth = buffer.depth(i);
            conic = buffer.conic(i);
            opac = buffer.opac(i);
            rgb = buffer.rgb(i);
        }

        __device__ __forceinline__ void store(ScreenBuffer &buffer, int64_t i) const {
            if (&buffer.xy(0)) buffer.xy(i) = xy;
            if (&buffer.depth(0)) buffer.depth(i) = depth;
            if (&buffer.conic(0)) buffer.conic(i) = conic;
            if (&buffer.opac(0)) buffer.opac(i) = opac;
            if (&buffer.rgb(0)) buffer.rgb(i) = rgb;
        }

        __device__ __forceinline__ void atomicStore(ScreenBuffer &buffer, int64_t i) const {
            if (&buffer.xy(0)) atomicAddFVec(&buffer.xy(i), xy);
            if (&buffer.depth(0)) atomicAddFVec(&buffer.depth(i), depth);
            if (&buffer.conic(0)) atomicAddFVec(&buffer.conic(i), conic);
            if (&buffer.opac(0)) atomicAddFVec(&buffer.opac(i), opac);
            if (&buffer.rgb(0)) atomicAddFVec(&buffer.rgb(i), rgb);
        }
    };
    #endif

    class ScreenBuffer : public TensorArray<5> {
    public:
        using TensorArray<5>::TensorArray;
    #ifdef __CUDACC__
        __forceinline__ __device__ float2& xy(int64_t i)
            { return *reinterpret_cast<float2*>(&_data[0][2*i]); }
        __forceinline__ __device__ float& depth(int64_t i)
            { return *reinterpret_cast<float*>(&_data[1][i]); }
        __forceinline__ __device__ float3& conic(int64_t i)
            { return *reinterpret_cast<float3*>(&_data[2][3*i]); }
        __forceinline__ __device__ float& opac(int64_t i)
            { return *reinterpret_cast<float*>(&_data[3][i]); }
        __forceinline__ __device__ float3& rgb(int64_t i)
            { return *reinterpret_cast<float3*>(&_data[4][3*i]); }
    
        __forceinline__ __device__ float2 xy(int64_t i) const
            { return *reinterpret_cast<float2*>(&_data[0][2*i]); }
        __forceinline__ __device__ float depth(int64_t i) const
            { return *reinterpret_cast<float*>(&_data[1][i]); }
        __forceinline__ __device__ float3 conic(int64_t i) const    
            { return *reinterpret_cast<float3*>(&_data[2][3*i]); }
        __forceinline__ __device__ float opac(int64_t i) const
            { return *reinterpret_cast<float*>(&_data[3][i]); }
        __forceinline__ __device__ float3 rgb(int64_t i) const
            { return *reinterpret_cast<float3*>(&_data[4][3*i]); }

    #endif  // #ifdef __CUDACC__

        static std::vector<DeviceTensorFloatND> empty_pool(int64_t size, PoolSlot key_prefix) {
            return TensorArray<5>::empty_pool(size, {2, 1, 3, 1, 3}, key_prefix);
        }
        static std::vector<DeviceTensorFloatND> zeros_pool(int64_t size, PoolSlot key_prefix) {
            return TensorArray<5>::zeros_pool(size, {2, 1, 3, 1, 3}, key_prefix);
        }
        static std::vector<DeviceTensorFloatND> zeros_pool(
            const std::vector<DeviceTensorFloatND>& tmpl, PoolSlot key_prefix) {
            return TensorArray<5>::zeros_pool(tmpl, key_prefix);
        }
    };

    #ifdef __CUDACC__
    struct FragmentFwd {
        FragmentFwd& operator=(const FragmentFwd&) = default;

        float2 xy;
        float depth;
        float3 conic;
        float opac;
        float3 rgb;

        static __device__ __forceinline__ Screen zero() {
            Screen s;
            s.xy = make_float2(0.f);
            s.depth = 0.f;
            s.conic = make_float3(0.f);
            s.opac = 0.f;
            s.rgb = make_float3(0.f);
            return s;
        }
    };
    using FragmentBwd = FragmentFwd;
    #endif

};



template<int _sh_degree>
struct _BasePrimitive3DGUT : _BasePrimitive3DGS<_sh_degree> {

    using WorldBuffer = typename _BasePrimitive3DGS<_sh_degree>::WorldBuffer;

    #ifdef __CUDACC__
    struct World : public _BasePrimitive3DGS<_sh_degree>::World {
        World() = default;
        __device__ World(const typename _BasePrimitive3DGS<_sh_degree>::World& other) {
            _BasePrimitive3DGS<_sh_degree>::World::operator=(other);
        }
        __device__ World& operator=(const typename _BasePrimitive3DGS<_sh_degree>::World& other) {
            _BasePrimitive3DGS<_sh_degree>::World::operator=(other);
            return *this;
        }
    };
    #endif

    class ScreenBuffer;

    #ifdef __CUDACC__
    struct Screen {
        float3 scale;  // post exp
        float opacity;  // post sigmoid
        float3 rgb;

        static __device__ __forceinline__ Screen zero() {
            Screen s;
            s.scale = make_float3(0.f);
            s.opacity = 0.f;
            s.rgb = make_float3(0.f);
            return s;
        }

        __device__ __forceinline__ void load(const ScreenBuffer &buffer, int64_t i) {
            scale = buffer.scales(i);
            opacity = buffer.opacities(i);
            rgb = buffer.colors(i);
        }

        __device__ __forceinline__ void store(ScreenBuffer &buffer, int64_t i) const {
            if (&buffer.scales(0)) buffer.scales(i) = scale;
            if (&buffer.opacities(0)) buffer.opacities(i) = opacity;
            if (&buffer.colors(0)) buffer.colors(i) = rgb;
        }

        __device__ __forceinline__ void atomicStore(ScreenBuffer &buffer, int64_t i) const {
            if (&buffer.scales(0)) atomicAddFVec(&buffer.scales(i), scale);
            if (&buffer.opacities(0)) atomicAddFVec(&buffer.opacities(i), opacity);
            if (&buffer.colors(0)) atomicAddFVec(&buffer.colors(i), rgb);
        }
    };
    #endif

    class ScreenBuffer : public TensorArray<3> {
    public:
        using TensorArray<3>::TensorArray;
    #ifdef __CUDACC__
        __forceinline__ __device__ float3& scales(int64_t i)
            { return *reinterpret_cast<float3*>(&_data[0][3*i]); }
        __forceinline__ __device__ float& opacities(int64_t i)
            { return *reinterpret_cast<float*>(&_data[1][i]); }
        __forceinline__ __device__ float3& colors(int64_t i)
            { return *reinterpret_cast<float3*>(&_data[2][3*i]); }
    
        __forceinline__ __device__ float3 scales(int64_t i) const
            { return *reinterpret_cast<float3*>(&_data[0][3*i]); }
        __forceinline__ __device__ float opacities(int64_t i) const
            { return *reinterpret_cast<float*>(&_data[1][i]); }
        __forceinline__ __device__ float3 colors(int64_t i) const
            { return *reinterpret_cast<float3*>(&_data[2][3*i]); }

    #endif  // #ifdef __CUDACC__

        static std::vector<DeviceTensorFloatND> empty_pool(int64_t size, PoolSlot key_prefix) {
            return TensorArray<3>::empty_pool(size, {3, 1, 3}, key_prefix);
        }
        static std::vector<DeviceTensorFloatND> zeros_pool(int64_t size, PoolSlot key_prefix) {
            return TensorArray<3>::zeros_pool(size, {-1, 1, 3}, key_prefix);
        }
        static std::vector<DeviceTensorFloatND> zeros_pool(
            const std::vector<DeviceTensorFloatND>& tmpl, PoolSlot key_prefix) {
            return TensorArray<3>::zeros_pool(tmpl, key_prefix);
        }
    };

};
