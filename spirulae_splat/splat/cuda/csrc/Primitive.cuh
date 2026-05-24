#pragma once

#include "types.cuh"
#include "common.cuh"

#include <array>
#include <tuple>
#include <optional>
#include <vector>

#include <Tensor.h>


enum class RenderOutputType {
    RGB, RGB_D, RGB_DN
};

class RenderOutput {
    static constexpr bool _has_depth(RenderOutputType type)
        { return type == RenderOutputType::RGB_D || type == RenderOutputType::RGB_DN; }
    static constexpr bool _has_normal(RenderOutputType type)
        { return type == RenderOutputType::RGB_DN; }
    static constexpr float _default_depth = 0.0f;
    static constexpr float3 _default_normal = {0.0f, 0.0f, 0.0f};

public:

    float3 rgb;
    float depth;
    float3 normal;

    struct Buffer;

    typedef std::tuple<
        DeviceTensor3D<float3>,
        DeviceTensor3D<float>,
        DeviceTensor3D<float3>
    > TensorTuple;

    template<RenderOutputType type>
    static void resize(
        TensorTuple& tensors,
        int64_t batch, int64_t height, int64_t width,
        const std::string& key
    ) {
        std::get<0>(tensors).resize(key + ".rgb", batch, height, width);
        if (_has_depth(type))
            std::get<1>(tensors).resize(key + ".depth", batch, height, width);
        if (_has_normal(type))
            std::get<2>(tensors).resize(key + ".normal", batch, height, width);
    }

    struct Tensor {
        DeviceTensor3D<float3> rgbs;
        DeviceTensor3D<float> depths;
        DeviceTensor3D<float3> normals;

        Tensor() {}

        Tensor(const TensorTuple& images) {
            rgbs = std::get<0>(images);
            depths = std::get<1>(images);
            normals = std::get<2>(images);
        }

        TensorTuple tuple() const {
            return std::make_tuple(rgbs, depths, normals);
        }

        bool has_value() const {
            return rgbs.data_ptr() != nullptr;
        }

        long width() const {
            return rgbs.size<2>();
        }
        long height() const {
            return rgbs.size<1>();
        }
        long batchSize() const {
            return rgbs.size<0>();
        }

        Buffer buffer() { return Buffer(*this); }
    };

    struct Buffer {
        float3* __restrict__ rgbs;
        float* __restrict__ depths;
        float3* __restrict__ normals;

        Buffer() : rgbs(nullptr), depths(nullptr), normals(nullptr) {}

        Buffer(const Tensor& tensors) {
            rgbs = tensors.rgbs.data_ptr();
            depths = tensors.depths.data_ptr();
            normals = tensors.normals.data_ptr();
        }

#ifdef __CUDACC__

    template<RenderOutputType type>
    __device__ RenderOutput load(long idx) const {
        return {
            rgbs[idx],
            _has_depth(type) ? depths[idx] : _default_depth,
            _has_normal(type) ? normals[idx] : _default_normal,
        };
    }

#endif

    };

#ifdef __CUDACC__

    static __device__ __forceinline__ RenderOutput zero() {
        return {
            {0.f, 0.f, 0.f},
            0.f,
            {0.f, 0.f, 0.f},
        };
    }

    __device__ __forceinline__ void operator+=(const RenderOutput &other) {
        rgb += other.rgb;
        depth += other.depth;
        normal += other.normal;
    }

    __device__ __forceinline__ RenderOutput operator*(float k) const {
        return {
            rgb * k,
            depth * k,
            normal * k,
        };
    }

    __device__ __forceinline__ RenderOutput operator+(const RenderOutput &other) const {
        return {
            rgb + other.rgb,
            depth + other.depth,
            normal + other.normal,
        };
    }

    __device__ __forceinline__ RenderOutput operator*(const RenderOutput &other) const {
        return {
            rgb * other.rgb,
            depth * other.depth,
            normal * other.normal,
        };
    }

    __device__ __forceinline__ float dot(const RenderOutput &other) const {
        float val = 0.0f;
        val += rgb.x * other.rgb.x + rgb.y * other.rgb.y + rgb.z * other.rgb.z;
        val += depth * other.depth;
        val += (normal.x * other.normal.x + normal.y * other.normal.y + normal.z * other.normal.z);
        return val;
    }

    template<RenderOutputType type>
    __device__ void saveParamsToBuffer(Buffer &buffer, long idx) {
        buffer.rgbs[idx] = rgb;
        if (_has_depth(type)) buffer.depths[idx] = depth;
        if (_has_normal(type)) buffer.normals[idx] = normal;
    }

#endif  // #ifdef __CUDACC__

};



#ifdef __CUDACC__

struct ProjCamera {
    float3x3 R;
    float3 t;
    float fx, fy, cx, cy;
    uint width, height;
    CameraDistortionCoeffs dist_coeffs;
};

#endif



template<int N>
class TensorArray {
protected:
    int64_t _size;
    float* _data[N];
    int32_t _strides[N];

    static_assert(N >= 1);

public:
    TensorArray() {
        _size = 0;
        for (int i = 0; i < N; ++i) {
            _data[i] = nullptr;
            _strides[i] = 0;
        }
    }

    TensorArray(const TensorArray& other) {
        _size = other._size;
        for (int i = 0; i < N; ++i) {
            _data[i] = other._data[i];
            _strides[i] = other._strides[i];
        }
    }

    int64_t size() const {
        return _size;
    }

    TensorArray(std::vector<DeviceTensorFloatND> tensors) {
        if (tensors.size() == 0) {
            _size = 0;
            for (int i = 0; i < N; ++i) {
                _data[i] = nullptr;
                _strides[i] = 0;
            }
            return;
        }
        if (tensors.size() != N)
            throw std::runtime_error("Number of tensors mismatch: Expect "
                + std::to_string(N) + ", got " + std::to_string(tensors.size()));
        _size = 0;
        for (int i = 0; i < N; ++i)
            if (tensors[i].data_ptr() != nullptr) { _size = tensors[i].size(0); break; }
        for (int i = 0; i < N; ++i) {
            _data[i] = tensors[i].data_ptr();
            if (_data[i] != nullptr) {
                if (tensors[i].size(0) != _size)
                    throw std::runtime_error("Tensor size mismatch");
                _strides[i] = _size == 0 ? 0 : (int32_t)(tensors[i].numel() / _size);
            } else {
                _strides[i] = 0;
            }
        }
    }

    static std::vector<DeviceTensorFloatND> empty_pool(
        int64_t size, std::array<int32_t, N> strides, const std::string& key_prefix
    ) {
        std::vector<DeviceTensorFloatND> res;
        for (int i = 0; i < N; ++i) {
            if (strides[i] <= 0) { res.push_back(DeviceTensorFloatND()); continue; }
            std::string key = key_prefix + "." + std::to_string(i);
            switch (strides[i]) {
                case 1: { DeviceVector<float>  b; b.resize(key, size); res.push_back(DeviceTensorFloatND(b)); break; }
                case 2: { DeviceVector<float2> b; b.resize(key, size); res.push_back(DeviceTensorFloatND(b)); break; }
                case 3: { DeviceVector<float3> b; b.resize(key, size); res.push_back(DeviceTensorFloatND(b)); break; }
                case 4: { DeviceVector<float4> b; b.resize(key, size); res.push_back(DeviceTensorFloatND(b)); break; }
                default: { DeviceTensor2D<float> b; b.resize(key, size, (int64_t)strides[i]); res.push_back(DeviceTensorFloatND(b)); break; }
            }
        }
        return res;
    }

    static std::vector<DeviceTensorFloatND> zeros_pool(
        int64_t size, std::array<int32_t, N> strides, const std::string& key_prefix
    ) {
        auto res = empty_pool(size, strides, key_prefix);
        for (auto& t : res)
            if (t.data_ptr() != nullptr)
                cudaMemset(t.data_ptr(), 0, (size_t)t.numel() * sizeof(float));
        return res;
    }

    static std::vector<DeviceTensorFloatND> zeros_pool(
        const std::vector<DeviceTensorFloatND>& tmpl_vec, const std::string& key_prefix
    ) {
        TensorArray<N> tmpl(tmpl_vec);
        std::array<int32_t, N> strides_arr;
        for (int i = 0; i < N; i++) strides_arr[i] = tmpl._strides[i];
        return zeros_pool(tmpl._size, strides_arr, key_prefix);
    }

#if 0

    static std::vector<DeviceTensorFloatND> empty_like(const TensorArray& other) {
        std::vector<DeviceTensorFloatND> res;
        for (int i = 0; i < N; ++i)
            res.push_back(
                at::empty({other._size, other._strides[i]}, kTensorOptionF32())
            );
        return res;
    }

    static std::vector<DeviceTensorFloatND> zeros_like(const TensorArray& other) {
        std::vector<DeviceTensorFloatND> res;
        for (int i = 0; i < N; ++i) {
            res.push_back(
                at::empty({other._size, other._strides[i]}, kTensorOptionF32())
            );
            set_zero_tensor(res.back().value());
        }
        return res;
    }

    static std::vector<DeviceTensorFloatND> empty(int64_t size, std::array<int32_t, N> strides) {
        std::vector<DeviceTensorFloatND> res;
        for (int i = 0; i < N; ++i)
            res.push_back(
                at::empty({size, strides[i]}, kTensorOptionF32())
            );
        return res;
    }

    static std::vector<DeviceTensorFloatND> zeros(int64_t size, std::array<int32_t, N> strides) {
        std::vector<DeviceTensorFloatND> res;
        for (int i = 0; i < N; ++i) {
            if (strides[i] < 0) {
                res.push_back(std::nullopt);
                continue;
            }
            res.push_back(
                at::empty({size, strides[i]}, kTensorOptionF32())
            );
            set_zero_tensor(res.back().value());
        }
        return res;
    }
#endif

};


namespace GlobalDeviceTensors {
    inline RenderOutput::TensorTuple renders;
    inline RenderOutput::TensorTuple renders2;
    inline RenderOutput::TensorTuple distortions;
    inline DeviceTensor3D<float> render_Ts;
    inline DeviceTensor3D<int32_t> render_last_ids;
}
