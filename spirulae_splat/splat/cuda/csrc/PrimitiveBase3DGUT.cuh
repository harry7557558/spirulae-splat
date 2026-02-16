#pragma once

#ifdef __CUDACC__
#include "generated/slang.cuh"
namespace SlangProjectionUtils {
#include "generated/set_namespace.cuh"
#include "generated/projection_utils.cuh"
}
#endif

#include "types.cuh"
#include "common.cuh"

#include <tuple>


struct Base3DGUT {
    struct World;
    struct Screen;
    struct RenderOutput;

#ifdef __CUDACC__

    struct FwdProjCamera {
        float3x3 R;
        float3 t;
        float fx, fy, cx, cy;
        uint width, height;
        float near_plane, far_plane;
        CameraDistortionCoeffs dist_coeffs;
    };

    struct BwdProjCamera {
        float3x3 R;
        float3 t;
        float fx, fy, cx, cy;
        uint width, height;
        CameraDistortionCoeffs dist_coeffs;
    };

#endif  // #ifdef __CUDACC__

};

struct Base3DGUT::World {
    float3 mean;
    float4 quat;
    float3 scale;
    float opacity;
};


struct Base3DGUT::RenderOutput {

    float3 rgb;
    float depth;

    typedef std::tuple<at::Tensor, at::Tensor> TensorTuple;

    struct Buffer;

    struct Tensor {
        at::Tensor rgbs;
        at::Tensor depths;

        Tensor() {}

        Tensor(const TensorTuple& images) {
            rgbs = std::get<0>(images);
            depths = std::get<1>(images);
        }

        TensorTuple tuple() const {
            return std::make_tuple(rgbs, depths);
        }

        static Tensor empty(at::DimVector dims, at::TensorOptions opt) {
            at::DimVector rgbs_dims(dims); rgbs_dims.append({3});
            at::DimVector depths_dims(dims); depths_dims.append({1});
            return Tensor(std::make_tuple(
                at::empty(rgbs_dims, opt),
                at::empty(depths_dims, opt)
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

    struct Buffer {
        float3* __restrict__ rgbs;
        float* __restrict__ depths;

        Buffer() : rgbs(nullptr), depths(nullptr) {}

        Buffer(const Tensor& tensors) {
            DEVICE_GUARD(tensors.rgbs);
            CHECK_INPUT(tensors.rgbs);
            CHECK_INPUT(tensors.depths);
            rgbs = (float3*)tensors.rgbs.data_ptr<float>();
            depths = tensors.depths.data_ptr<float>();
        }
    };

#ifdef __CUDACC__


    static __device__ RenderOutput load(const Buffer &buffer, long idx) {
        return {
            buffer.rgbs[idx],
            buffer.depths[idx]
        };
    }

    static __device__ __forceinline__ RenderOutput zero() {
        return {
            {0.f, 0.f, 0.f},
            0.f
        };
    }

    __device__ __forceinline__ void operator+=(const RenderOutput &other) {
        rgb += other.rgb;
        depth += other.depth;
    }

    __device__ __forceinline__ RenderOutput operator*(float k) const {
        return {rgb * k, depth * k};
    }

    __device__ __forceinline__ RenderOutput operator+(const RenderOutput &other) const {
        return {rgb + other.rgb, depth + other.depth};
    }

    __device__ __forceinline__ RenderOutput operator*(const RenderOutput &other) const {
        return {rgb * other.rgb, depth * other.depth};
    }

    __device__ __forceinline__ float dot(const RenderOutput &other) const {
        return (rgb.x * other.rgb.x + rgb.y * other.rgb.y + rgb.z * other.rgb.z) + depth * other.depth;
    }

    __device__ void saveParamsToBuffer(Buffer &buffer, long idx) {
        buffer.rgbs[idx] = rgb;
        buffer.depths[idx] = depth;
    }

#endif  // #ifdef __CUDACC__
};


struct Base3DGUT::Screen {

    // from world
    float3 mean;
    float4 quat;
    // from screen
    float depth;
    float3 scale;
    float opacity;
    float3 rgb;
    // precompute
#ifdef __CUDACC__
    float3x3 iscl_rot;
#endif

    typedef std::tuple<at::Tensor, at::Tensor, at::Tensor, at::Tensor> TensorTupleProj;
    typedef std::tuple<at::Tensor, at::Tensor, at::Tensor, at::Tensor, at::Tensor, at::Tensor> TensorTuple;

    struct Buffer;

    struct Tensor {
        bool hasWorld;
        at::Tensor means;
        at::Tensor quats;
        at::Tensor depths;
        at::Tensor scales;
        at::Tensor opacities;
        at::Tensor rgbs;

        Tensor() {}

        Tensor(const TensorTuple& splats) : hasWorld(true) {
            means = std::get<0>(splats);
            quats = std::get<1>(splats);
            depths = std::get<2>(splats);
            scales = std::get<3>(splats);
            opacities = std::get<4>(splats);
            rgbs = std::get<5>(splats);
        }

        Tensor(const TensorTupleProj& splats) : hasWorld(false) {
            depths = std::get<0>(splats);
            scales = std::get<1>(splats);
            opacities = std::get<2>(splats);
            rgbs = std::get<3>(splats);
        }

        TensorTupleProj tupleProjFwd() const {
            return std::make_tuple(depths, scales, opacities, rgbs);
        }

        TensorTuple tupleProjFwdPacked() const {
            return std::make_tuple(means, quats, depths, scales, opacities, rgbs);
        }

        TensorTuple tupleRasterBwd() const {
            return std::make_tuple(means, quats, depths, scales, opacities, rgbs);
        }

        static TensorTupleProj allocProjFwd(long C, long N, c10::TensorOptions opt) {
            return std::make_tuple(
                at::empty({C, N}, opt),
                at::empty({C, N, 3}, opt),
                at::empty({C, N}, opt),
                at::empty({C, N, 3}, opt)
            );
        }

        static TensorTuple allocProjFwdPacked(long nnz, c10::TensorOptions opt) {
            return std::make_tuple(
                at::empty({nnz, 3}, opt),
                at::empty({nnz, 4}, opt),
                at::empty({nnz}, opt),
                at::empty({nnz, 3}, opt),
                at::empty({nnz}, opt),
                at::empty({nnz, 3}, opt)
            );
        }

        Tensor allocRasterBwd() const {
            if (!hasWorld)
                throw std::runtime_error("!hasWorld");
            Tensor result = Tensor(std::make_tuple(
                at::zeros_like(means),
                at::zeros_like(quats),
                at::zeros_like(depths),
                at::zeros_like(scales),
                at::zeros_like(opacities),
                at::zeros_like(rgbs)
            ));
            return result;
        }

        auto options() const {
            return rgbs.options();
        }
        bool isPacked() const {
            return rgbs.dim() == 2;
        }
        long size() const {
            return rgbs.size(-2);
        }

        Buffer buffer() { return Buffer(*this); }
    };

    struct Buffer {
        float3* __restrict__ means;
        float4* __restrict__ quats;
        float* __restrict__ depths;
        float3* __restrict__ scales;
        float* __restrict__ opacities;
        float3* __restrict__ rgbs;
        long size;

        Buffer() {}  // uninitialized

        Buffer(const Tensor& tensors) {
            DEVICE_GUARD(tensors.depths);
            if (tensors.hasWorld) {
                CHECK_INPUT(tensors.means);
                CHECK_INPUT(tensors.quats);
                means = (float3*)tensors.means.data_ptr<float>();
                quats = (float4*)tensors.quats.data_ptr<float>();
            } else { means = nullptr; quats = nullptr; };
            CHECK_INPUT(tensors.depths);
            CHECK_INPUT(tensors.scales);
            CHECK_INPUT(tensors.opacities);
            CHECK_INPUT(tensors.rgbs);
            depths = tensors.depths.data_ptr<float>();
            scales = (float3*)tensors.scales.data_ptr<float>();
            opacities = tensors.opacities.data_ptr<float>();
            rgbs = (float3*)tensors.rgbs.data_ptr<float>();
            size = tensors.hasWorld ?
                tensors.quats.numel() / 4 : tensors.depths.numel();
        }
    };

#ifdef __CUDACC__

    static __device__ Screen load(const Buffer &buffer, long idx) {
        return {
            buffer.means ? buffer.means[idx % buffer.size] : make_float3(0.f),
            buffer.quats ? buffer.quats[idx % buffer.size] : make_float4(0.f),
            buffer.depths[idx],
            buffer.scales[idx],
            buffer.opacities[idx],
            buffer.rgbs[idx],
        };
    }

    static __device__ Screen loadWithPrecompute(const Buffer &buffer, long idx) {
        Screen splat = Screen::load(buffer, idx);
        splat.iscl_rot = SlangProjectionUtils::compute_3dgut_iscl_rot(splat.quat, splat.scale);
        return splat;
    }

    static __device__ __forceinline__ Screen zero() {
        return {
            {0.f, 0.f, 0.f},
            {0.f, 0.f, 0.f, 0.f},
            0.0f,
            {0.f, 0.f, 0.f},
            0.0f,
            {0.f, 0.f, 0.f},
            {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f},
        };
    }

    __device__ __forceinline__ void addGradient(const Screen &grad, float weight=1.0f) {
        mean += grad.mean * weight;
        quat += grad.quat * weight;
        depth += grad.depth * weight;
        scale += grad.scale * weight;
        opacity += grad.opacity * weight;
        rgb += grad.rgb * weight;
        #pragma unroll
        for (int i = 0; i < 3; ++i)
            iscl_rot[i] = iscl_rot[i] + grad.iscl_rot[i] * weight;
    }

    __device__ __forceinline__ void addGaussNewtonHessianDiagonal(const Screen &grad, float weight=1.0f) {
        float4 v_quat; float3 v_scale;
        SlangProjectionUtils::compute_3dgut_iscl_rot_vjp(quat, scale, grad.iscl_rot, &v_quat, &v_scale);
        mean += grad.mean * grad.mean * weight;
        quat += (grad.quat + v_quat) * (grad.quat + v_quat) * weight;
        depth += grad.depth * grad.depth * weight;
        scale += (grad.scale + v_scale) * (grad.scale + v_scale) * weight;
        opacity += grad.opacity * grad.opacity * weight;
        rgb += grad.rgb * grad.rgb * weight;
    }

    __device__ void saveParamsToBuffer(Buffer &buffer, long idx) {
        if (buffer.means) buffer.means[idx % buffer.size] = mean;
        if (buffer.quats) buffer.quats[idx % buffer.size] = quat;
        buffer.depths[idx] = depth;
        buffer.scales[idx] = scale;
        buffer.opacities[idx] = opacity;
        buffer.rgbs[idx] = rgb;
        // iscl_rot is not saved
    }

    __device__ void atomicAddGradientToBuffer(const Screen &grad, Buffer &buffer, long idx) const {
        float4 v_quat; float3 v_scale;
        SlangProjectionUtils::compute_3dgut_iscl_rot_vjp(quat, scale, grad.iscl_rot, &v_quat, &v_scale);
        atomicAddFVec(buffer.means + idx % buffer.size, grad.mean);
        atomicAddFVec(buffer.quats + idx % buffer.size, grad.quat + v_quat);
        atomicAddFVec(buffer.depths + idx, grad.depth);
        atomicAddFVec(buffer.scales + idx, grad.scale + v_scale);
        atomicAddFVec(buffer.opacities + idx, grad.opacity);
        atomicAddFVec(buffer.rgbs + idx, grad.rgb);
    }

    __device__ void atomicAddGaussNewtonHessianDiagonalToBuffer(const Screen &grad, Buffer &buffer, long idx, float weight=1.0f) const {
        float4 v_quat; float3 v_scale;
        SlangProjectionUtils::compute_3dgut_iscl_rot_vjp(quat, scale, grad.iscl_rot, &v_quat, &v_scale);
        atomicAddFVec(buffer.means + idx % buffer.size, grad.mean * grad.mean * weight);
        atomicAddFVec(buffer.quats + idx % buffer.size, (grad.quat + v_quat) * (grad.quat + v_quat) * weight);
        atomicAddFVec(buffer.depths + idx, grad.depth * grad.depth * weight);
        atomicAddFVec(buffer.scales + idx, (grad.scale + v_scale) * (grad.scale + v_scale) * weight);
        atomicAddFVec(buffer.opacities + idx, grad.opacity * grad.opacity * weight);
        atomicAddFVec(buffer.rgbs + idx, grad.rgb * grad.rgb * weight);
    }

    __device__ __forceinline__ float evaluate_alpha(float3 ray_o, float3 ray_d) {
        if (dot(mean-ray_o, ray_d) <= 0.0f)
            return 0.0;
        return SlangProjectionUtils::evaluate_alpha_3dgs(
            mean, iscl_rot, opacity,
            ray_o, ray_d
        );
    }

    __device__ __forceinline__ Screen evaluate_alpha_vjp(
        float3 ray_o, float3 ray_d, float v_alpha,
        float3 &v_ray_o, float3 &v_ray_d
    ) {
        Screen v_splat = Screen::zero();
        if (dot(mean-ray_o, ray_d) <= 0.0f) {
            v_ray_o = v_ray_d = make_float3(0.f);
            return v_splat;
        }
        SlangProjectionUtils::evaluate_alpha_3dgs_vjp(
            mean, iscl_rot, opacity,
            ray_o, ray_d, v_alpha,
            &v_splat.mean, &v_splat.iscl_rot, &v_splat.opacity,
            &v_ray_o, &v_ray_d
        );
        return v_splat;
    }

    __device__ __forceinline__ Base3DGUT::RenderOutput evaluate_color(float3 ray_o, float3 ray_d) {
        // return {rgb, depth};
        float3 out_rgb; float out_depth;
        SlangProjectionUtils::evaluate_color_3dgs(
            mean, iscl_rot, opacity, rgb,
            ray_o, ray_d, &out_rgb, &out_depth
        );
        return {out_rgb, out_depth};
    }

    __device__ __forceinline__ Screen evaluate_color_vjp(
        float3 ray_o, float3 ray_d, Base3DGUT::RenderOutput v_render,
        float3 &v_ray_o, float3 &v_ray_d
    ) {
        Screen v_splat = Screen::zero();
    #if 0
        v_splat.rgb = v_render.rgb;
        v_splat.depth = v_render.depth;
        v_ray_o = make_float3(0.f);
        v_ray_d = make_float3(0.f);
    #else
        SlangProjectionUtils::evaluate_color_3dgs_vjp(
            mean, iscl_rot, opacity, rgb,
            ray_o, ray_d,
            v_render.rgb, v_render.depth,
            &v_splat.mean, &v_splat.iscl_rot, &v_splat.opacity, &v_splat.rgb,
            &v_ray_o, &v_ray_d
        );
    #endif
        return v_splat;
    }

#endif  // #ifdef __CUDACC__

};
