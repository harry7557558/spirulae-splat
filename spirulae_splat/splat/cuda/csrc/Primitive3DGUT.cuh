#pragma once

#ifdef __CUDACC__
#include "generated/slang.cuh"
namespace Slang3DGS {
#include "generated/set_namespace.cuh"
#include "generated/primitive_3dgs.cuh"
}
#endif

#include "PrimitiveBase3DGUT.cuh"


struct Vanilla3DGUT : public Base3DGUT {
    struct World;

#ifdef __CUDACC__

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

    inline static __device__ void project_persp_vjp(
        World world, BwdProjCamera cam,
        Screen v_proj, Screen vr_proj, Screen h_proj,
        World& v_world, float3x3 &v_R, float3 &v_t,
        float3 &vr_world_pos, float3 &h_world_pos
    );

    inline static __device__ void project_persp_vjp(
        World world, BwdProjCamera cam,
        Screen v_proj, Screen vr_proj, Screen h_proj,
        World& v_world, float3x3 &v_R, float3 &v_t,
        World& vr_world, World& h_world
    );

    inline static __device__ void project_fisheye_vjp(
        World world, BwdProjCamera cam,
        Screen v_proj,
        World& v_world, float3x3 &v_R, float3 &v_t
    );

    inline static __device__ void project_fisheye_vjp(
        World world, BwdProjCamera cam,
        Screen v_proj, Screen vr_proj, Screen h_proj,
        World& v_world, float3x3 &v_R, float3 &v_t,
        float3 &vr_world_pos, float3 &h_world_pos
    );

    inline static __device__ void project_fisheye_vjp(
        World world, BwdProjCamera cam,
        Screen v_proj, Screen vr_proj, Screen h_proj,
        World& v_world, float3x3 &v_R, float3 &v_t,
        World& vr_world, World& h_world
    );

#endif  // #ifdef __CUDACC__

};


struct Vanilla3DGUT::World : public Base3DGUT::World {

#ifdef __CUDACC__
    FixedArray<float3, 16> sh_coeffs;
#endif

    typedef std::tuple<at::Tensor, at::Tensor, at::Tensor, at::Tensor, at::Tensor, std::optional<at::Tensor>> TensorTuple;

    struct Buffer;

    struct Tensor {
        at::Tensor means;
        at::Tensor quats;
        at::Tensor scales;
        at::Tensor opacities;
        at::Tensor features_dc;
        std::optional<at::Tensor> features_sh;

        Tensor() {}

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

        Tensor allocProjBwd(bool is_hess_diag) const {
            return Tensor(std::make_tuple(
                at::zeros_like(means),
                at::zeros_like(quats),
                at::zeros_like(scales),
                at::zeros_like(opacities),
                at::zeros_like(features_dc),
                features_sh.has_value() && !is_hess_diag ?
                    (std::optional<at::Tensor>)at::zeros_like(features_sh.value()) :
                    (std::optional<at::Tensor>)std::nullopt
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

        Buffer() {}

        Buffer(const Tensor& tensors) {
            DEVICE_GUARD(tensors.means);
            CHECK_INPUT(tensors.means);
            CHECK_INPUT(tensors.quats);
            CHECK_INPUT(tensors.scales);
            CHECK_INPUT(tensors.opacities);
            CHECK_INPUT(tensors.features_dc);
            if (tensors.features_sh.has_value())
                CHECK_INPUT(tensors.features_sh.value());
            means = (float3*)tensors.means.data_ptr<float>();
            quats = (float4*)tensors.quats.data_ptr<float>();
            scales = (float3*)tensors.scales.data_ptr<float>();
            opacities = tensors.opacities.data_ptr<float>();
            features_dc = (float3*)tensors.features_dc.data_ptr<float>();
            features_sh = tensors.features_sh.has_value() ?
                (float3*)tensors.features_sh.value().template data_ptr<float>() : nullptr;
            num_sh = tensors.features_sh.has_value() ?
                tensors.features_sh.value().size(-2) : 0;
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
        World world;
        world.mean = {0.f, 0.f, 0.f};
        world.quat = {0.f, 0.f, 0.f, 0.f};
        world.scale = {0.f, 0.f, 0.f};
        world.opacity = 0.f;
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


#ifdef __CUDACC__

inline __device__ void Vanilla3DGUT::project_persp(
    Vanilla3DGUT::World world, Vanilla3DGUT::FwdProjCamera cam,
    Vanilla3DGUT::Screen& proj, int4& aabb
) {
    float2 xy;
    Slang3DGS::projection_3dgut_persp(
        false,
        world.mean, world.quat, world.scale, world.opacity, world.sh_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
        cam.width, cam.height, cam.near_plane, cam.far_plane,
        &aabb, &xy, &proj.depth, &proj.scale, &proj.opacity, &proj.rgb
    );
    proj.mean = world.mean, proj.quat = world.quat;
}

inline __device__ void Vanilla3DGUT::project_fisheye(
    Vanilla3DGUT::World world, Vanilla3DGUT::FwdProjCamera cam,
    Vanilla3DGUT::Screen& proj, int4& aabb
) {
    float2 xy;
    Slang3DGS::projection_3dgut_fisheye(
        false,
        world.mean, world.quat, world.scale, world.opacity, world.sh_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
        cam.width, cam.height, cam.near_plane, cam.far_plane,
        &aabb, &xy, &proj.depth, &proj.scale, &proj.opacity, &proj.rgb
    );
    proj.mean = world.mean, proj.quat = world.quat;
}

inline __device__ void Vanilla3DGUT::project_persp_vjp(
    Vanilla3DGUT::World world, Vanilla3DGUT::BwdProjCamera cam,
    Vanilla3DGUT::Screen v_proj,
    Vanilla3DGUT::World& v_world, float3x3 &v_R, float3 &v_t
) {
    Slang3DGS::projection_3dgut_persp_vjp(
        false,
        world.mean, world.quat, world.scale, world.opacity, world.sh_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
        cam.width, cam.height,
        make_float2(0), v_proj.depth, v_proj.scale, v_proj.opacity, v_proj.rgb,
        &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity, &v_world.sh_coeffs,
        &v_R, &v_t
    );
    v_world.mean = v_proj.mean, v_world.quat = v_proj.quat;
}

inline __device__ void Vanilla3DGUT::project_persp_vjp(
    Vanilla3DGUT::World world, Vanilla3DGUT::BwdProjCamera cam,
    Vanilla3DGUT::Screen v_proj, Vanilla3DGUT::Screen vr_proj, Vanilla3DGUT::Screen h_proj,
    Vanilla3DGUT::World& v_world, float3x3 &v_R, float3 &v_t,
    float3& vr_world_pos, float3& h_world_pos
) {
    Slang3DGS::projection_3dgut_persp_vjp(
        false,
        world.mean, world.quat, world.scale, world.opacity, world.sh_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
        cam.width, cam.height,
        make_float2(0), v_proj.depth, v_proj.scale, v_proj.opacity, v_proj.rgb,
        make_float2(0), vr_proj.depth, vr_proj.scale, vr_proj.opacity, vr_proj.rgb,
        make_float2(0), h_proj.depth, h_proj.scale, h_proj.opacity, h_proj.rgb,
        &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity, &v_world.sh_coeffs,
        &v_R, &v_t, &vr_world_pos, &h_world_pos
    );
    v_world.mean = v_proj.mean, v_world.quat = v_proj.quat;
}

inline __device__ void Vanilla3DGUT::project_persp_vjp(
    Vanilla3DGUT::World world, Vanilla3DGUT::BwdProjCamera cam,
    Vanilla3DGUT::Screen v_proj, Vanilla3DGUT::Screen vr_proj, Vanilla3DGUT::Screen h_proj,
    Vanilla3DGUT::World& v_world, float3x3 &v_R, float3 &v_t,
    Vanilla3DGUT::World& vr_world, Vanilla3DGUT::World& h_world
) {
    Slang3DGS::projection_3dgut_persp_vjp(
        false,
        world.mean, world.quat, world.scale, world.opacity, world.sh_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
        cam.width, cam.height,
        make_float2(0), v_proj.depth, v_proj.scale, v_proj.opacity, v_proj.rgb,
        make_float2(0), vr_proj.depth, vr_proj.scale, vr_proj.opacity, vr_proj.rgb,
        make_float2(0), h_proj.depth, h_proj.scale, h_proj.opacity, h_proj.rgb,
        &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity, &v_world.sh_coeffs,
        &v_R, &v_t,
        &vr_world.mean, &vr_world.quat, &vr_world.scale, &vr_world.opacity, &vr_world.sh_coeffs[0],
        &h_world.mean, &h_world.quat, &h_world.scale, &h_world.opacity, &h_world.sh_coeffs[0]
    );
}

inline __device__ void Vanilla3DGUT::project_fisheye_vjp(
    Vanilla3DGUT::World world, Vanilla3DGUT::BwdProjCamera cam,
    Vanilla3DGUT::Screen v_proj,
    Vanilla3DGUT::World& v_world, float3x3 &v_R, float3 &v_t
) {
    Slang3DGS::projection_3dgut_fisheye_vjp(
        false,
        world.mean, world.quat, world.scale, world.opacity, world.sh_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
        cam.width, cam.height,
        make_float2(0), v_proj.depth, v_proj.scale, v_proj.opacity, v_proj.rgb,
        &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity, &v_world.sh_coeffs,
        &v_R, &v_t
    );
    v_world.mean = v_proj.mean, v_world.quat = v_proj.quat;
}

inline __device__ void Vanilla3DGUT::project_fisheye_vjp(
    Vanilla3DGUT::World world, Vanilla3DGUT::BwdProjCamera cam,
    Vanilla3DGUT::Screen v_proj, Vanilla3DGUT::Screen vr_proj, Vanilla3DGUT::Screen h_proj,
    Vanilla3DGUT::World& v_world, float3x3 &v_R, float3 &v_t,
    float3& vr_world_pos, float3& h_world_pos
) {
    Slang3DGS::projection_3dgut_fisheye_vjp(
        false,
        world.mean, world.quat, world.scale, world.opacity, world.sh_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
        cam.width, cam.height,
        make_float2(0), v_proj.depth, v_proj.scale, v_proj.opacity, v_proj.rgb,
        make_float2(0), vr_proj.depth, vr_proj.scale, vr_proj.opacity, vr_proj.rgb,
        make_float2(0), h_proj.depth, h_proj.scale, h_proj.opacity, h_proj.rgb,
        &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity, &v_world.sh_coeffs,
        &v_R, &v_t, &vr_world_pos, &h_world_pos
    );
    v_world.mean = v_proj.mean, v_world.quat = v_proj.quat;
}

inline __device__ void Vanilla3DGUT::project_fisheye_vjp(
    Vanilla3DGUT::World world, Vanilla3DGUT::BwdProjCamera cam,
    Vanilla3DGUT::Screen v_proj, Vanilla3DGUT::Screen vr_proj, Vanilla3DGUT::Screen h_proj,
    Vanilla3DGUT::World& v_world, float3x3 &v_R, float3 &v_t,
    Vanilla3DGUT::World& vr_world, Vanilla3DGUT::World& h_world
) {
    Slang3DGS::projection_3dgs_fisheye_vjp(
        false,
        world.mean, world.quat, world.scale, world.opacity, world.sh_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
        cam.width, cam.height,
        make_float2(0), v_proj.depth, v_proj.scale, v_proj.opacity, v_proj.rgb,
        make_float2(0), vr_proj.depth, vr_proj.scale, vr_proj.opacity, vr_proj.rgb,
        make_float2(0), h_proj.depth, h_proj.scale, h_proj.opacity, h_proj.rgb,
        &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity, &v_world.sh_coeffs,
        &v_R, &v_t,
        &vr_world.mean, &vr_world.quat, &vr_world.scale, &vr_world.opacity, &vr_world.sh_coeffs[0],
        &h_world.mean, &h_world.quat, &h_world.scale, &h_world.opacity, &h_world.sh_coeffs[0]
    );
}

#endif  // #ifdef __CUDACC__
