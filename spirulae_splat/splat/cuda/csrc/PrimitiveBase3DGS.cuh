#pragma once

#ifdef __CUDACC__
#include "generated/slang.cuh"
namespace SlangProjectionUtils {
#include "generated/set_namespace.cuh"
#include "generated/projection_utils.cuh"
}
#endif

#include "Primitive.cuh"


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
        // TODO: register pressure
        FixedArray<float3, 16> sh_coeffs;

        static __device__ __forceinline__ World zero() {
            World w;
            w.mean = make_float3(0.f);
            w.quat = make_float4(0.f);
            w.scale = make_float3(0.f);
            w.opacity = 0.f;
            for (int i = 0; i < 16; ++i)
                w.sh_coeffs[i] = make_float3(0.f);
            return w;
        }

        __device__ __forceinline__ void load(const WorldBuffer &buffer, int64_t i) {
            mean = buffer.means(i);
            quat = buffer.quats(i);
            scale = buffer.scales(i);
            opacity = buffer.opacities(i);
            sh_coeffs[0] = buffer.features_dc(i);
            for (int j = 0; j < 15; ++j)
                sh_coeffs[j+1] = j < buffer.num_sh() ? buffer.features_sh(i, j) : make_float3(0.f);
        }

        __device__ __forceinline__ void store(WorldBuffer &buffer, int64_t i) const {
            if (&buffer.means(0)) buffer.means(i) = mean;
            if (&buffer.quats(0)) buffer.quats(i) = quat;
            if (&buffer.scales(0)) buffer.scales(i) = scale;
            if (&buffer.opacities(0)) buffer.opacities(i) = opacity;
            if (&buffer.features_dc(0)) buffer.features_dc(i) = sh_coeffs[0];
            for (int j = 0; j < buffer.num_sh(); ++j)
                if (&buffer.features_sh(0, 0)) buffer.features_sh(i, j) = sh_coeffs[j+1];
        }

        __device__ __forceinline__ void atomicStore(WorldBuffer &buffer, int64_t i) const {
            if (&buffer.means(0)) atomicAddFVec(&buffer.means(i), mean);
            if (&buffer.quats(0)) atomicAddFVec(&buffer.quats(i), quat);
            if (&buffer.scales(0)) atomicAddFVec(&buffer.scales(i), scale);
            if (&buffer.opacities(0)) atomicAddFVec(&buffer.opacities(i), opacity);
            if (&buffer.features_dc(0)) atomicAddFVec(&buffer.features_dc(i), sh_coeffs[0]);
            for (int j = 0; j < buffer.num_sh(); ++j)
                if (&buffer.features_sh(0, 0)) atomicAddFVec(&buffer.features_sh(i, j), sh_coeffs[j+1]);
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
        __forceinline__ __device__ int32_t num_sh() const
            { return _strides[5] / 3; }
        __forceinline__ __device__ float3& features_sh(int64_t i, int64_t j)
            { return *reinterpret_cast<float3*>(&_data[5][_strides[5]*i+3*j]); }

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
        __forceinline__ __device__ float3 features_sh(int64_t i, int64_t j) const
            { return *reinterpret_cast<float3*>(&_data[5][_strides[5]*i+3*j]); }

    #endif  // #ifdef __CUDACC__

    #ifndef NO_TORCH
        static TensorList empty(int64_t size, int num_sh) {
            return TensorArray<6>::empty(size, {3, 4, 3, 1, 3, (int32_t)(3*num_sh)});
        }
        static TensorList zeros(int64_t size, int num_sh) {
            return TensorArray<6>::zeros(size, {3, 4, 3, 1, 3, (int32_t)(3*num_sh)});
        }
    #endif  // #ifndef NO_TORCH
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

    #ifndef NO_TORCH
        static TensorList empty(int64_t size) {
            return TensorArray<5>::empty(size, {2, 1, 3, 1, 3});
        }
        static TensorList zeros(int64_t size) {
            return TensorArray<5>::zeros(size, {2, 1, 3, 1, 3});
        }
    #endif  // #ifndef NO_TORCH
    };

    #ifdef __CUDACC__
    struct Fragment {
        Fragment& operator=(const Fragment&) = default;

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
    #endif

};



struct _BasePrimitive3DGUT : _BasePrimitive3DGS {

    class WorldBuffer : public _BasePrimitive3DGS::WorldBuffer {
        using _BasePrimitive3DGS::WorldBuffer::WorldBuffer;
    };

    #ifdef __CUDACC__
    struct World : public _BasePrimitive3DGS::World {
        __device__ World() = default;
        __device__ World(const _BasePrimitive3DGS::World& other) {
            _BasePrimitive3DGS::World::operator=(other);
        }
        __device__ World& operator=(const _BasePrimitive3DGS::World& other) {
            _BasePrimitive3DGS::World::operator=(other);
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

    #ifndef NO_TORCH
        static TensorList empty(int64_t size) {
            return TensorArray<3>::empty(size, {3, 1, 3});
        }
        static TensorList zeros(int64_t size) {
            return TensorArray<3>::zeros(size, {3, 1, 3});
        }
    #endif  // #ifndef NO_TORCH
    };

};
