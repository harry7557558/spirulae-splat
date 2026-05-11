#pragma once

#include "types.cuh"
#include "common.cuh"

#include <array>
#include <tuple>
#include <optional>
#include <vector>

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

#ifndef NO_TORCH

    typedef std::tuple<
        at::Tensor,
        std::optional<at::Tensor>,
        std::optional<at::Tensor>
    > TensorTuple;

    struct Tensor {
        at::Tensor rgbs;
        std::optional<at::Tensor> depths = std::nullopt;
        std::optional<at::Tensor> normals = std::nullopt;

        Tensor() {}

        Tensor(const TensorTuple& images) {
            rgbs = std::get<0>(images);
            depths = std::get<1>(images);
            normals = std::get<2>(images);
        }

        TensorTuple tuple() const {
            return std::make_tuple(rgbs, depths, normals);
        }

        template<RenderOutputType type>
        static Tensor empty(at::DimVector dims) {
            at::DimVector rgbs_dims(dims); rgbs_dims.append({3});
            at::DimVector depths_dims(dims); depths_dims.append({1});
            at::DimVector normals_dims(dims); normals_dims.append({3});
            return Tensor(std::make_tuple(
                at::empty(rgbs_dims, kTensorOptionF32()),
                _has_depth(type) ? at::empty(depths_dims, kTensorOptionF32()) :
                    (std::optional<at::Tensor>)std::nullopt,
                _has_normal(type) ? at::empty(normals_dims, kTensorOptionF32()) :
                    (std::optional<at::Tensor>)std::nullopt
            ));
        }

        auto options() const {
            return rgbs.options();
        }
        long width() const {
            return rgbs.size(-2);
        }
        long height() const {
            return rgbs.size(-3);
        }
        long batchSize() const {
            return rgbs.numel() / (3*width()*height());
        }

        Buffer buffer() { return Buffer(*this); }
    };

#endif  // #ifndef NO_TORCH

    struct Buffer {
        float3* __restrict__ rgbs;
        float* __restrict__ depths;
        float3* __restrict__ normals;

        Buffer() : rgbs(nullptr), depths(nullptr), normals(nullptr) {}

        #ifndef NO_TORCH
        Buffer(const Tensor& tensors) {
            DEVICE_GUARD(tensors.rgbs);
            CHECK_INPUT(tensors.rgbs);
            rgbs = (float3*)tensors.rgbs.data_ptr<float>();
            if (tensors.depths) {
                CHECK_INPUT(tensors.depths.value());
                depths = tensors.depths.value().data_ptr<float>();
            } else depths = nullptr;
            if (tensors.normals) {
                CHECK_INPUT(tensors.normals.value());
                normals = (float3*)tensors.normals.value().data_ptr<float>();
            }
        }
        #endif

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


#ifndef NO_TORCH
typedef std::vector<std::optional<at::Tensor>> TensorList;
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

#ifndef NO_TORCH
    TensorArray(std::vector<at::Tensor> tensors) {
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
        DEVICE_GUARD(tensors[0]);
        for (int i = 0; i < N; ++i) {
            CHECK_INPUT(tensors[i]);
            if (i == 0)
                _size = tensors[i].size(0);
            else if (_size != tensors[i].size(0))
                throw std::runtime_error("Tensor size mismatch");
            _data[i] = tensors[i].data_ptr<float>();
            if (_data[i] != nullptr)
                _strides[i] = tensors[i].numel() / _size;
        }
    }

    TensorArray(TensorList tensors) {
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
        bool saw_first = false;
        for (int i = 0; i < N; ++i) {
            if (tensors[i].has_value()) {
                if (!saw_first) {
                    DEVICE_GUARD(tensors[i].value());
                    _size = tensors[i].value().size(0);
                    saw_first = true;
                } else if (_size != tensors[i].value().size(0)) {
                    throw std::runtime_error("Tensor size mismatch");
                }
                CHECK_INPUT(tensors[i].value());
                _data[i] = tensors[i].value().data_ptr<float>();
                _strides[i] = tensors[i].value().numel() / _size;
            } else {
                _data[i] = nullptr;
                _strides[i] = 0;
            }
        }
        if (!saw_first)
            _size = 0;
    }

    static TensorList empty_like(const TensorArray& other) {
        TensorList res;
        for (int i = 0; i < N; ++i)
            res.push_back(
                at::empty({other._size, other._strides[i]}, kTensorOptionF32())
            );
        return res;
    }

    static TensorList zeros_like(const TensorArray& other) {
        TensorList res;
        for (int i = 0; i < N; ++i) {
            res.push_back(
                at::empty({other._size, other._strides[i]}, kTensorOptionF32())
            );
            set_zero_tensor(res.back().value());
        }
        return res;
    }

    static TensorList empty(int64_t size, std::array<int32_t, N> strides) {
        TensorList res;
        for (int i = 0; i < N; ++i)
            res.push_back(
                at::empty({size, strides[i]}, kTensorOptionF32())
            );
        return res;
    }

    static TensorList zeros(int64_t size, std::array<int32_t, N> strides) {
        TensorList res;
        for (int i = 0; i < N; ++i) {
            res.push_back(
                at::empty({size, strides[i]}, kTensorOptionF32())
            );
            set_zero_tensor(res.back().value());
        }
        return res;
    }
#endif

};
