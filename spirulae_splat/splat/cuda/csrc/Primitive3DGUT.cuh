#pragma once

#ifdef __CUDACC__
#define TensorView _Slang_TensorView
#include "generated/primitive.cu"
#undef TensorView
#endif

#include "types.cuh"
#include "common.cuh"

#include <tuple>


struct Vanilla3DGUT {
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

    inline static __device__ void project_persp(
        World world, FwdProjCamera cam,
        Screen& proj, int4& aabb
    );

    inline static __device__ void project_fisheye(
        World world, FwdProjCamera cam,
        Screen& proj, int4& aabb
    );

    struct BwdProjCamera {
        float3x3 R;
        float3 t;
        float fx, fy, cx, cy;
        uint width, height;
        CameraDistortionCoeffs dist_coeffs;
    };

    inline static __device__ void project_persp_vjp(
        World world, BwdProjCamera cam,
        Screen v_proj,
        World& v_world, float3x3 &v_R, float3 &v_t
    );

    inline static __device__ void project_fisheye_vjp(
        World world, BwdProjCamera cam,
        Screen v_proj,
        World& v_world, float3x3 &v_R, float3 &v_t
    );

#endif  // #ifdef __CUDACC__

    static constexpr float eps2d = 0.3f;
};

struct Vanilla3DGUT::World {

    float3 mean;
    float4 quat;
    float3 scale;
    float opacity;
    FixedArray<float3, 16> sh_coeffs;

    typedef std::tuple<at::Tensor, at::Tensor, at::Tensor, at::Tensor, at::Tensor, at::Tensor> TensorTuple;

    struct Buffer;

    struct Tensor {
        at::Tensor means;
        at::Tensor quats;
        at::Tensor scales;
        at::Tensor opacities;
        at::Tensor features_dc;
        at::Tensor features_sh;

        Tensor(const TensorTuple& splats) {
            means = std::get<0>(splats);
            quats = std::get<1>(splats);
            scales = std::get<2>(splats);
            opacities = std::get<3>(splats);
            features_dc = std::get<4>(splats);
            features_sh = std::get<5>(splats);
        }

        TensorTuple tuple() const {
            return std::make_tuple(means, quats, scales, opacities, features_dc, features_sh);
        }

        Tensor zeros_like() const {
            return Tensor(std::make_tuple(
                at::zeros_like(means),
                at::zeros_like(quats),
                at::zeros_like(scales),
                at::zeros_like(opacities),
                at::zeros_like(features_dc),
                at::zeros_like(features_sh)
            ));
        }

        auto options() const {
            return means.options();
        }
        long size() const {
            return quats.size(-2);
        }
        long batchSize() const {
            return quats.numel() / (4*size());
        }

        Buffer buffer() { return Buffer(*this); }
    };

    struct Buffer {
        float3* __restrict__ means;
        float4* __restrict__ quats;
        float3* __restrict__ scales;
        float* __restrict__ opacities;
        float3* __restrict__ features_dc;
        float3* __restrict__ features_sh;
        uint num_sh;

        Buffer(const Tensor& tensors) {
            DEVICE_GUARD(tensors.means);
            CHECK_INPUT(tensors.means);
            CHECK_INPUT(tensors.quats);
            CHECK_INPUT(tensors.scales);
            CHECK_INPUT(tensors.opacities);
            CHECK_INPUT(tensors.features_dc);
            CHECK_INPUT(tensors.features_sh);
            means = (float3*)tensors.means.data_ptr<float>();
            quats = (float4*)tensors.quats.data_ptr<float>();
            scales = (float3*)tensors.scales.data_ptr<float>();
            opacities = tensors.opacities.data_ptr<float>();
            features_dc = (float3*)tensors.features_dc.data_ptr<float>();
            features_sh = (float3*)tensors.features_sh.data_ptr<float>();
            num_sh = tensors.features_sh.size(-2);
        }
    };

#ifdef __CUDACC__

    static __device__ World load(const Buffer &buffer, long idx) {
        World world = {
            buffer.means[idx],
            buffer.quats[idx],
            buffer.scales[idx],
            buffer.opacities[idx]
        };
        world.sh_coeffs[0] = buffer.features_dc[idx];
        for (int i = 0; i < 15; i++)
            world.sh_coeffs[i+1] = i < buffer.num_sh ?
                buffer.features_sh[idx*buffer.num_sh+i] : make_float3(0);
        return world;
    }

    static __device__ __forceinline__ World zero() {
        World world = {
            {0.f, 0.f, 0.f},
            {0.f, 0.f, 0.f, 0.f},
            {0.f, 0.f, 0.f},
            0.f,
        };
        for (int i = 0; i < 16; i++)
            world.sh_coeffs[i] = make_float3(0);
        return world;
    }

    template<typename Partition>
    __device__ void reduce(Partition& partition) {
        warpSum(mean, partition);
        warpSum(quat, partition);
        warpSum(scale, partition);
        warpSum(opacity, partition);
        for (int i = 0; i < 16; i++)
            warpSum(sh_coeffs[i], partition);
    }
    
    __device__ void saveParamsToBuffer(Buffer &buffer, long idx) {
        buffer.means[idx] = mean;
        buffer.quats[idx] = quat;
        buffer.scales[idx] = scale;
        buffer.opacities[idx] = opacity;
        buffer.features_dc[idx] = sh_coeffs[0];
        for (int i = 0; i < buffer.num_sh; i++)
            buffer.features_sh[idx*buffer.num_sh + i] = sh_coeffs[i+1];
    }

    __device__ void atomicAddGradientToBuffer(Buffer &buffer, long idx) {
        atomicAddFVec(buffer.means + idx, mean);
        atomicAddFVec(buffer.quats + idx, quat);
        atomicAddFVec(buffer.scales + idx, scale);
        atomicAddFVec(buffer.opacities + idx, opacity);
        atomicAddFVec(buffer.features_dc + idx, sh_coeffs[0]);
        for (int i = 0; i < buffer.num_sh; i++)
            atomicAddFVec(buffer.features_sh + idx*buffer.num_sh + i, sh_coeffs[i+1]);
    }

#endif  // #ifdef __CUDACC__
};


struct Vanilla3DGUT::RenderOutput {

    float3 rgb;
    float depth;

    typedef std::tuple<at::Tensor, at::Tensor> TensorTuple;

    struct Buffer;

    struct Tensor {
        at::Tensor rgbs;
        at::Tensor depths;

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


struct Vanilla3DGUT::Screen {

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

        TensorTuple tupleRasterBwd() const {
            return std::make_tuple(means, quats, depths, scales, opacities, rgbs);
        }

        static Tensor allocProjFwd(long C, long N, c10::TensorOptions opt) {
            return std::make_tuple(
                at::empty({N, 3}, opt),
                at::empty({N, 4}, opt),
                C == -1 ? at::empty({N}, opt) : at::empty({C, N}, opt),
                C == -1 ? at::empty({N, 3}, opt) : at::empty({C, N, 3}, opt),
                C == -1 ? at::empty({N}, opt) : at::empty({C, N}, opt),
                C == -1 ? at::empty({N, 3}, opt) : at::empty({C, N, 3}, opt)
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
            return rgbs.numel() / 3;
        }

        Buffer buffer() { return Buffer(*this); }
    };

    struct Buffer {
        float3* __restrict__ means = nullptr;
        float4* __restrict__ quats = nullptr;
        float* __restrict__ depths;
        float3* __restrict__ scales;
        float* __restrict__ opacities;
        float3* __restrict__ rgbs;
        long size;

        Buffer(const Tensor& tensors) {
            DEVICE_GUARD(tensors.depths);
            if (tensors.hasWorld) {
                CHECK_INPUT(tensors.means);
                CHECK_INPUT(tensors.quats);
                means = (float3*)tensors.means.data_ptr<float>();
                quats = (float4*)tensors.quats.data_ptr<float>();
            }
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
        splat.iscl_rot = compute_3dgut_iscl_rot(splat.quat, splat.scale);
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

    __device__ __forceinline__ void operator+=(const Screen &other) {
        mean += other.mean;
        quat += other.quat;
        depth += other.depth;
        scale += other.scale;
        opacity += other.opacity;
        rgb += other.rgb;
        iscl_rot = iscl_rot + other.iscl_rot;
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
        compute_3dgut_iscl_rot_vjp(quat, scale, grad.iscl_rot, &v_quat, &v_scale);
        atomicAddFVec(buffer.means + idx % buffer.size, grad.mean);
        atomicAddFVec(buffer.quats + idx % buffer.size, grad.quat + v_quat);
        atomicAddFVec(buffer.depths + idx, grad.depth);
        atomicAddFVec(buffer.scales + idx, grad.scale + v_scale);
        atomicAddFVec(buffer.opacities + idx, grad.opacity);
        atomicAddFVec(buffer.rgbs + idx, grad.rgb);
    }

    __device__ __forceinline__ float evaluate_alpha(float3 ray_o, float3 ray_d) {
        if (dot(mean-ray_o, ray_d) <= 0.0f)
            return 0.0;
        return evaluate_alpha_3dgs(
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
        evaluate_alpha_3dgs_vjp(
            mean, iscl_rot, opacity,
            ray_o, ray_d, v_alpha,
            &v_splat.mean, &v_splat.iscl_rot, &v_splat.opacity,
            &v_ray_o, &v_ray_d
        );
        return v_splat;
    }

    __device__ __forceinline__ Vanilla3DGUT::RenderOutput evaluate_color(float3 ray_o, float3 ray_d) {
        return {rgb, depth};
    }

    __device__ __forceinline__ Screen evaluate_color_vjp(
        float3 ray_o, float3 ray_d, Vanilla3DGUT::RenderOutput v_render,
        float3 &v_ray_o, float3 &v_ray_d
    ) {
        Screen v_splat = Screen::zero();
        v_splat.rgb = v_render.rgb;
        v_splat.depth = v_render.depth;
        return v_splat;
    }

#endif  // #ifdef __CUDACC__

};


#ifdef __CUDACC__

inline __device__ void Vanilla3DGUT::project_persp(
    Vanilla3DGUT::World world, Vanilla3DGUT::FwdProjCamera cam,
    Vanilla3DGUT::Screen& proj, int4& aabb
) {
    float2 xy;
    projection_3dgs_eval3d_persp(
        false,
        world.mean, world.quat, world.scale, world.opacity, &world.sh_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, &cam.dist_coeffs,
        cam.width, cam.height, cam.near_plane, cam.far_plane,
        &aabb, &xy, &proj.depth, &proj.scale, &proj.opacity, &proj.rgb
    );
}

inline __device__ void Vanilla3DGUT::project_fisheye(
    Vanilla3DGUT::World world, Vanilla3DGUT::FwdProjCamera cam,
    Vanilla3DGUT::Screen& proj, int4& aabb
) {
    float2 xy;
    projection_3dgs_eval3d_fisheye(
        false,
        world.mean, world.quat, world.scale, world.opacity, &world.sh_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, &cam.dist_coeffs,
        cam.width, cam.height, cam.near_plane, cam.far_plane,
        &aabb, &xy, &proj.depth, &proj.scale, &proj.opacity, &proj.rgb
    );
}

inline __device__ void Vanilla3DGUT::project_persp_vjp(
    Vanilla3DGUT::World world, Vanilla3DGUT::BwdProjCamera cam,
    Vanilla3DGUT::Screen v_proj,
    Vanilla3DGUT::World& v_world, float3x3 &v_R, float3 &v_t
) {
    projection_3dgs_eval3d_persp_vjp(
        false,
        world.mean, world.quat, world.scale, world.opacity, &world.sh_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, &cam.dist_coeffs,
        cam.width, cam.height,
        make_float2(0), v_proj.depth, v_proj.scale, v_proj.opacity, v_proj.rgb,
        &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity, &v_world.sh_coeffs,
        &v_R, &v_t
    );
}

inline __device__ void Vanilla3DGUT::project_fisheye_vjp(
    Vanilla3DGUT::World world, Vanilla3DGUT::BwdProjCamera cam,
    Vanilla3DGUT::Screen v_proj,
    Vanilla3DGUT::World& v_world, float3x3 &v_R, float3 &v_t
) {
    projection_3dgs_eval3d_fisheye_vjp(
        false,
        world.mean, world.quat, world.scale, world.opacity, &world.sh_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, &cam.dist_coeffs,
        cam.width, cam.height,
        make_float2(0), v_proj.depth, v_proj.scale, v_proj.opacity, v_proj.rgb,
        &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity, &v_world.sh_coeffs,
        &v_R, &v_t
    );
}

#endif  // #ifdef __CUDACC__
