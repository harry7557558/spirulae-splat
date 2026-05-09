#pragma once

#ifdef __CUDACC__
#include "generated/slang.cuh"
namespace Slang3DGS {
#include "generated/set_namespace.cuh"
#include "generated/primitive_3dgs.cuh"
}
#endif

#include "PrimitiveBase3DGS.cuh"

struct Vanilla3DGUT : public _BasePrimitive3DGS {
    static constexpr RenderOutputType pixelType = RenderOutputType::RGB_D;
};

#if 0
struct Vanilla3DGUT : public Base3DGUT {
    struct World;

#ifdef __CUDACC__

    inline static __device__ void project_persp(
        World world, ProjCamera cam,
        Screen& proj, float& sorting_depth, float4& aabb
    );

    inline static __device__ void project_fisheye(
        World world, ProjCamera cam,
        Screen& proj, float& sorting_depth, float4& aabb
    );

    inline static __device__ void project_persp_vjp(
        World world, ProjCamera cam,
        Screen v_proj,
        World& v_world, float3x3 &v_R, float3 &v_t
    );

    inline static __device__ void project_persp_vjp(
        World world, ProjCamera cam,
        Screen v_proj, Screen vr_proj, Screen h_proj,
        World& v_world, float3x3 &v_R, float3 &v_t,
        float3 &vr_world_pos, float3 &h_world_pos
    );

    inline static __device__ void project_persp_vjp(
        World world, ProjCamera cam,
        Screen v_proj, Screen vr_proj, Screen h_proj,
        World& v_world, float3x3 &v_R, float3 &v_t,
        World& vr_world, World& h_world
    );

    inline static __device__ void project_fisheye_vjp(
        World world, ProjCamera cam,
        Screen v_proj,
        World& v_world, float3x3 &v_R, float3 &v_t
    );

    inline static __device__ void project_fisheye_vjp(
        World world, ProjCamera cam,
        Screen v_proj, Screen vr_proj, Screen h_proj,
        World& v_world, float3x3 &v_R, float3 &v_t,
        float3 &vr_world_pos, float3 &h_world_pos
    );

    inline static __device__ void project_fisheye_vjp(
        World world, ProjCamera cam,
        Screen v_proj, Screen vr_proj, Screen h_proj,
        World& v_world, float3x3 &v_R, float3 &v_t,
        World& vr_world, World& h_world
    );

#endif  // #ifdef __CUDACC__

};


struct Vanilla3DGUT::World : public Base3DGUT::World {

#ifdef __CUDACC__

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

    __device__ void saveParamsToBuffer(Buffer &buffer, long idx) {
        buffer.means(idx) = mean;
        if (buffer.hasQuats()) buffer.quats(idx) = quat;
        buffer.scales(idx) = scale;
        buffer.opacities(idx) = opacity;
        buffer.features_dc(idx) = sh_coeffs[0];
        for (int i = 0; i < buffer.num_sh; i++)
            buffer.features_sh(idx, i) = sh_coeffs[i+1];
    }

    __device__ void atomicAddGradientToBuffer(Buffer &buffer, long idx) {
        atomicAddFVec(buffer.means(idx), mean);
        if (buffer.hasQuats()) atomicAddFVec(buffer.quats(idx), quat);
        atomicAddFVec(buffer.scales(idx), scale);
        atomicAddFVec(buffer.opacities(idx), opacity);
        atomicAddFVec(buffer.features_dc(idx), sh_coeffs[0]);
        for (int i = 0; i < buffer.num_sh; i++)
            atomicAddFVec(buffer.features_sh(idx, i), sh_coeffs[i+1]);
    }

#endif  // #ifdef __CUDACC__
};
#endif


#ifdef __CUDACC__

inline __device__ void project_persp(
    Vanilla3DGUT::World world, ProjCamera cam,
    Vanilla3DGUT::Screen& proj, float& sorting_depth, float4& aabb
) {
    float2 xy;
    Slang3DGS::projection_3dgut_persp(
        false,
        world.mean, world.quat, world.scale, world.opacity, world.sh_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
        cam.width, cam.height,
        &aabb, &xy, &sorting_depth, &proj.scale, &proj.opacity, &proj.rgb
    );
}

inline __device__ void project_fisheye(
    Vanilla3DGUT::World world, ProjCamera cam,
    Vanilla3DGUT::Screen& proj, float& sorting_depth, float4& aabb
) {
    float2 xy;
    Slang3DGS::projection_3dgut_fisheye(
        false,
        world.mean, world.quat, world.scale, world.opacity, world.sh_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
        cam.width, cam.height,
        &aabb, &xy, &sorting_depth, &proj.scale, &proj.opacity, &proj.rgb
    );
}

inline __device__ void project_persp_vjp(
    Vanilla3DGUT::World world, ProjCamera cam,
    Vanilla3DGUT::Screen v_proj,
    Vanilla3DGUT::World& v_world, float3x3 &v_R, float3 &v_t
) {
    Slang3DGS::projection_3dgut_persp_vjp(
        false,
        world.mean, world.quat, world.scale, world.opacity, world.sh_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
        cam.width, cam.height,
        make_float2(0), 0.0f, v_proj.scale, v_proj.opacity, v_proj.rgb,
        &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity, &v_world.sh_coeffs,
        &v_R, &v_t
    );
}

inline __device__ void project_persp_vjp(
    Vanilla3DGUT::World world, ProjCamera cam,
    Vanilla3DGUT::Screen v_proj, Vanilla3DGUT::Screen vr_proj, Vanilla3DGUT::Screen h_proj,
    Vanilla3DGUT::World& v_world, float3x3 &v_R, float3 &v_t,
    float3& vr_world_pos, float3& h_world_pos
) {
    Slang3DGS::projection_3dgut_persp_vjp(
        false,
        world.mean, world.quat, world.scale, world.opacity, world.sh_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
        cam.width, cam.height,
        make_float2(0), 0.0f, v_proj.scale, v_proj.opacity, v_proj.rgb,
        make_float2(0), 0.0f, vr_proj.scale, vr_proj.opacity, vr_proj.rgb,
        make_float2(0), 0.0f, h_proj.scale, h_proj.opacity, h_proj.rgb,
        &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity, &v_world.sh_coeffs,
        &v_R, &v_t, &vr_world_pos, &h_world_pos
    );
}

inline __device__ void project_persp_vjp(
    Vanilla3DGUT::World world, ProjCamera cam,
    Vanilla3DGUT::Screen v_proj, Vanilla3DGUT::Screen vr_proj, Vanilla3DGUT::Screen h_proj,
    Vanilla3DGUT::World& v_world, float3x3 &v_R, float3 &v_t,
    Vanilla3DGUT::World& vr_world, Vanilla3DGUT::World& h_world
) {
    Slang3DGS::projection_3dgut_persp_vjp(
        false,
        world.mean, world.quat, world.scale, world.opacity, world.sh_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
        cam.width, cam.height,
        make_float2(0), 0.0f, v_proj.scale, v_proj.opacity, v_proj.rgb,
        make_float2(0), 0.0f, vr_proj.scale, vr_proj.opacity, vr_proj.rgb,
        make_float2(0), 0.0f, h_proj.scale, h_proj.opacity, h_proj.rgb,
        &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity, &v_world.sh_coeffs,
        &v_R, &v_t,
        &vr_world.mean, &vr_world.quat, &vr_world.scale, &vr_world.opacity, &vr_world.sh_coeffs[0],
        &h_world.mean, &h_world.quat, &h_world.scale, &h_world.opacity, &h_world.sh_coeffs[0]
    );
}

inline __device__ void project_fisheye_vjp(
    Vanilla3DGUT::World world, ProjCamera cam,
    Vanilla3DGUT::Screen v_proj,
    Vanilla3DGUT::World& v_world, float3x3 &v_R, float3 &v_t
) {
    Slang3DGS::projection_3dgut_fisheye_vjp(
        false,
        world.mean, world.quat, world.scale, world.opacity, world.sh_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
        cam.width, cam.height,
        make_float2(0), 0.0f, v_proj.scale, v_proj.opacity, v_proj.rgb,
        &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity, &v_world.sh_coeffs,
        &v_R, &v_t
    );
}

inline __device__ void project_fisheye_vjp(
    Vanilla3DGUT::World world, ProjCamera cam,
    Vanilla3DGUT::Screen v_proj, Vanilla3DGUT::Screen vr_proj, Vanilla3DGUT::Screen h_proj,
    Vanilla3DGUT::World& v_world, float3x3 &v_R, float3 &v_t,
    float3& vr_world_pos, float3& h_world_pos
) {
    Slang3DGS::projection_3dgut_fisheye_vjp(
        false,
        world.mean, world.quat, world.scale, world.opacity, world.sh_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
        cam.width, cam.height,
        make_float2(0), 0.0f, v_proj.scale, v_proj.opacity, v_proj.rgb,
        make_float2(0), 0.0f, vr_proj.scale, vr_proj.opacity, vr_proj.rgb,
        make_float2(0), 0.0f, h_proj.scale, h_proj.opacity, h_proj.rgb,
        &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity, &v_world.sh_coeffs,
        &v_R, &v_t, &vr_world_pos, &h_world_pos
    );
}

inline __device__ void project_fisheye_vjp(
    Vanilla3DGUT::World world, ProjCamera cam,
    Vanilla3DGUT::Screen v_proj, Vanilla3DGUT::Screen vr_proj, Vanilla3DGUT::Screen h_proj,
    Vanilla3DGUT::World& v_world, float3x3 &v_R, float3 &v_t,
    Vanilla3DGUT::World& vr_world, Vanilla3DGUT::World& h_world
) {
    Slang3DGS::projection_3dgs_fisheye_vjp(
        false,
        world.mean, world.quat, world.scale, world.opacity, world.sh_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
        cam.width, cam.height,
        make_float2(0), 0.0f, v_proj.scale, v_proj.opacity, v_proj.rgb,
        make_float2(0), 0.0f, vr_proj.scale, vr_proj.opacity, vr_proj.rgb,
        make_float2(0), 0.0f, h_proj.scale, h_proj.opacity, h_proj.rgb,
        &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity, &v_world.sh_coeffs,
        &v_R, &v_t,
        &vr_world.mean, &vr_world.quat, &vr_world.scale, &vr_world.opacity, &vr_world.sh_coeffs[0],
        &h_world.mean, &h_world.quat, &h_world.scale, &h_world.opacity, &h_world.sh_coeffs[0]
    );
}

#endif  // #ifdef __CUDACC__
