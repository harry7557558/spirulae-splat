#pragma once

#ifdef __CUDACC__
#include "generated/slang.cuh"
namespace Slang3DGS {
#include "generated/set_namespace.cuh"
#include "generated/primitive_3dgs.cuh"
}
#endif

#include "PrimitiveBase3DGS.cuh"

struct Vanilla3DGUT : public _BasePrimitive3DGUT {
    static constexpr RenderOutputType pixelType = RenderOutputType::RGB_D;

    class WorldBuffer : public _BasePrimitive3DGUT::WorldBuffer {
        using _BasePrimitive3DGUT::WorldBuffer::WorldBuffer;
    };
    class ScreenBuffer : public _BasePrimitive3DGUT::ScreenBuffer {
        using _BasePrimitive3DGUT::ScreenBuffer::ScreenBuffer;
    };

#ifdef __CUDACC__

    struct Screen;

    struct World : public _BasePrimitive3DGUT::World{
        __device__ World() = default;
        __device__ World(const _BasePrimitive3DGS::World& other) {
            _BasePrimitive3DGS::World::operator=(other);
        }
        __device__ World& operator=(const _BasePrimitive3DGS::World& other) {
            _BasePrimitive3DGS::World::operator=(other);
            return *this;
        }
        __device__ World(const _BasePrimitive3DGUT::World& other) {
            _BasePrimitive3DGUT::World::operator=(other);
        }
        __device__ World& operator=(const _BasePrimitive3DGUT::World& other) {
            _BasePrimitive3DGUT::World::operator=(other);
            return *this;
        }

        template<ssplat::CameraModelType camera_model>
        inline __device__ void project(
            ProjCamera cam,
            Vanilla3DGUT::Screen& proj, float4& aabb, float& sorting_depth
        ) const {
            float2 xy;
            if constexpr (camera_model == ssplat::CameraModelType::PINHOLE)
                Slang3DGS::projection_3dgut_persp(
                    false,
                    mean, quat, scale, opacity, sh_coeffs,
                    cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
                    cam.width, cam.height,
                    &aabb, &xy, &sorting_depth, &proj.scale, &proj.opacity, &proj.rgb
                );
            else if constexpr (camera_model == ssplat::CameraModelType::FISHEYE)
                Slang3DGS::projection_3dgut_fisheye(
                    false,
                    mean, quat, scale, opacity, sh_coeffs,
                    cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
                    cam.width, cam.height,
                    &aabb, &xy, &sorting_depth, &proj.scale, &proj.opacity, &proj.rgb
                );
        }

        template<ssplat::CameraModelType camera_model>
        inline __device__ void project_vjp(
            ProjCamera cam,
            Vanilla3DGUT::Screen v_proj,
            Vanilla3DGUT::World& v_world, float3x3 &v_R, float3 &v_t
        ) const {
            if constexpr (camera_model == ssplat::CameraModelType::PINHOLE)
                Slang3DGS::projection_3dgut_persp_vjp(
                    false,
                    mean, quat, scale, opacity, sh_coeffs,
                    cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
                    cam.width, cam.height,
                    make_float2(0), 0.0f, v_proj.scale, v_proj.opacity, v_proj.rgb,
                    &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity, &v_world.sh_coeffs,
                    &v_R, &v_t
                );
            else if constexpr (camera_model == ssplat::CameraModelType::FISHEYE)
                Slang3DGS::projection_3dgut_fisheye_vjp(
                    false,
                    mean, quat, scale, opacity, sh_coeffs,
                    cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
                    cam.width, cam.height,
                    make_float2(0), 0.0f, v_proj.scale, v_proj.opacity, v_proj.rgb,
                    &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity, &v_world.sh_coeffs,
                    &v_R, &v_t
                );
        }

        template<ssplat::CameraModelType camera_model>
        inline __device__ void project_vjp_h_pos(
            ProjCamera cam,
            Vanilla3DGUT::Screen v_proj, Vanilla3DGUT::Screen vr_proj, Vanilla3DGUT::Screen h_proj,
            Vanilla3DGUT::World& v_world, float3x3 &v_R, float3 &v_t,
            float3& vr_world_pos, float3& h_world_pos
        ) const {
            if constexpr (camera_model == ssplat::CameraModelType::PINHOLE)
                Slang3DGS::projection_3dgut_persp_vjp(
                    false,
                    mean, quat, scale, opacity, sh_coeffs,
                    cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
                    cam.width, cam.height,
                    make_float2(0), 0.0f, v_proj.scale, v_proj.opacity, v_proj.rgb,
                    make_float2(0), 0.0f, vr_proj.scale, vr_proj.opacity, vr_proj.rgb,
                    make_float2(0), 0.0f, h_proj.scale, h_proj.opacity, h_proj.rgb,
                    &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity, &v_world.sh_coeffs,
                    &v_R, &v_t, &vr_world_pos, &h_world_pos
                );
            else if constexpr (camera_model == ssplat::CameraModelType::FISHEYE)
                Slang3DGS::projection_3dgut_fisheye_vjp(
                    false,
                    mean, quat, scale, opacity, sh_coeffs,
                    cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
                    cam.width, cam.height,
                    make_float2(0), 0.0f, v_proj.scale, v_proj.opacity, v_proj.rgb,
                    make_float2(0), 0.0f, vr_proj.scale, vr_proj.opacity, vr_proj.rgb,
                    make_float2(0), 0.0f, h_proj.scale, h_proj.opacity, h_proj.rgb,
                    &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity, &v_world.sh_coeffs,
                    &v_R, &v_t, &vr_world_pos, &h_world_pos
                );
        }

        template<ssplat::CameraModelType camera_model>
        inline __device__ void project_vjp_h_all(
            ProjCamera cam,
            Vanilla3DGUT::Screen v_proj, Vanilla3DGUT::Screen vr_proj, Vanilla3DGUT::Screen h_proj,
            Vanilla3DGUT::World& v_world, float3x3 &v_R, float3 &v_t,
            Vanilla3DGUT::World& vr_world, Vanilla3DGUT::World& h_world
        ) const {
            if constexpr (camera_model == ssplat::CameraModelType::PINHOLE)
                Slang3DGS::projection_3dgut_persp_vjp(
                    false,
                    mean, quat, scale, opacity, sh_coeffs,
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
            else if constexpr (camera_model == ssplat::CameraModelType::FISHEYE)
                Slang3DGS::projection_3dgut_fisheye_vjp(
                    false,
                    mean, quat, scale, opacity, sh_coeffs,
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

    };

    struct Screen : public _BasePrimitive3DGUT::Screen{
        __device__ Screen() = default;
        __device__ Screen(const _BasePrimitive3DGUT::Screen& other) {
            _BasePrimitive3DGUT::Screen::operator=(other);
        }
        __device__ Screen& operator=(const _BasePrimitive3DGUT::Screen& other) {
            _BasePrimitive3DGUT::Screen::operator=(other);
            return *this;
        }
    };

    struct Fragment {
        float3 mean;
        // float3x3 iscl_rot;
        float4 quat;
        float3 scale;
        float opacity;
        float3 rgb;

        static __device__ __forceinline__ Fragment zero() {
            Fragment f;
            f.mean = make_float3(0.f);
            // f.iscl_rot = float3x3{0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
            f.quat = make_float4(0.f);
            f.scale = make_float3(0.f);
            f.opacity = 0.f;
            f.rgb = make_float3(0.f);
            return f;
        }

        __device__ __forceinline__ void load(
            const WorldBuffer &wbuffer,
            const ScreenBuffer &sbuffer,
            int64_t wi, int64_t si
        ) {
            mean = wbuffer.means(wi);
            // iscl_rot = SlangProjectionUtils::compute_3dgut_iscl_rot(wbuffer.quats(wi), sbuffer.scales(si));
            quat = wbuffer.quats(wi);
            scale = sbuffer.scales(si);
            opacity = sbuffer.opacities(si);
            rgb = sbuffer.colors(si);
        }

        __device__ __forceinline__ void store(
            WorldBuffer &wbuffer,
            ScreenBuffer &sbuffer,
            int64_t wi, int64_t si
        ) const {
            if (&wbuffer.means(0)) wbuffer.means(wi) = mean;
            if (&wbuffer.quats(0)) wbuffer.quats(wi) = quat;
            if (&sbuffer.scales(0)) sbuffer.scales(si) = scale;
            if (&sbuffer.opacities(0)) sbuffer.opacities(si) = opacity;
            if (&sbuffer.colors(0)) sbuffer.colors(si) = rgb;
        }

        __device__ __forceinline__ void atomicStore(
            WorldBuffer &wbuffer,
            ScreenBuffer &sbuffer,
            int64_t wi, int64_t si
        ) const {
            if (&wbuffer.means(0)) atomicAddFVec(&wbuffer.means(wi), mean);
            if (&wbuffer.quats(0)) atomicAddFVec(&wbuffer.quats(wi), quat);
            if (&sbuffer.scales(0)) atomicAddFVec(&sbuffer.scales(si), scale);
            if (&sbuffer.opacities(0)) atomicAddFVec(&sbuffer.opacities(si), opacity);
            if (&sbuffer.colors(0)) atomicAddFVec(&sbuffer.colors(si), rgb);
        }

        __device__ __forceinline__ float evaluate_alpha(
            float3 ray_o, float3 ray_d
        ) const {
            // if (dot(mean-ray_o, ray_d) <= 0.0f)
            //     return 0.0;
            float3x3 iscl_rot = SlangProjectionUtils::compute_3dgut_iscl_rot(quat, scale);
            return SlangProjectionUtils::evaluate_alpha_3dgs(
                mean, iscl_rot, opacity,
                ray_o, ray_d
            );
        }

        __device__ __forceinline__ void evaluate_alpha_vjp(
            float3 ray_o, float3 ray_d, float v_alpha,
            Fragment &v_frag, float3 &v_ray_o, float3 &v_ray_d
        ) const {
            // if (dot(mean-ray_o, ray_d) <= 0.0f) {
            //     v_ray_o = v_ray_d = make_float3(0.f);
            //     return;
            // }
            float3 v_mean; float3x3 v_iscl_rot; float v_opacity;
            float3 v_ray_o_t, v_ray_d_t;
            float3x3 iscl_rot = SlangProjectionUtils::compute_3dgut_iscl_rot(quat, scale);
            SlangProjectionUtils::evaluate_alpha_3dgs_vjp(
                mean, iscl_rot, opacity,
                ray_o, ray_d, v_alpha,
                &v_mean, &v_iscl_rot, &v_opacity,
                &v_ray_o_t, &v_ray_d_t
            );
            float4 v_quat; float3 v_scale;
            SlangProjectionUtils::compute_3dgut_iscl_rot_vjp(quat, scale, v_iscl_rot, &v_quat, &v_scale);
            v_frag.mean += v_mean;
            // v_frag.iscl_rot = v_frag.iscl_rot + v_iscl_rot;
            v_frag.quat += v_quat;
            v_frag.scale += v_scale;
            v_frag.opacity += v_opacity;
            v_ray_o += v_ray_o_t;
            v_ray_d += v_ray_d_t;
        }

        __device__ __forceinline__ RenderOutput evaluate_color(
            float3 ray_o, float3 ray_d
        ) const {
            float3 out_rgb; float out_depth;
            float3x3 iscl_rot = SlangProjectionUtils::compute_3dgut_iscl_rot(quat, scale);
            SlangProjectionUtils::evaluate_color_3dgs(
                mean, iscl_rot, opacity, rgb,
                ray_o, ray_d, &out_rgb, &out_depth
            );
            return { out_rgb, out_depth, make_float3(0.0f) };
        }

        __device__ __forceinline__ void evaluate_color_vjp(
            float3 ray_o, float3 ray_d, RenderOutput v_render,
            Fragment &v_frag, float3 &v_ray_o, float3 &v_ray_d
        ) const {
            float3 v_mean; float3x3 v_iscl_rot; float v_opacity; float3 v_rgb;
            float3 v_ray_o_t, v_ray_d_t;
            float3x3 iscl_rot = SlangProjectionUtils::compute_3dgut_iscl_rot(quat, scale);
            SlangProjectionUtils::evaluate_color_3dgs_vjp(
                mean, iscl_rot, opacity, rgb,
                ray_o, ray_d,
                v_render.rgb, v_render.depth,
                &v_mean, &v_iscl_rot, &v_opacity, &v_rgb,
                &v_ray_o_t, &v_ray_d_t
            );
            float4 v_quat; float3 v_scale;
            SlangProjectionUtils::compute_3dgut_iscl_rot_vjp(quat, scale, v_iscl_rot, &v_quat, &v_scale);
            v_frag.mean += v_mean;
            // v_frag.iscl_rot = v_frag.iscl_rot + v_iscl_rot;
            v_frag.quat += v_quat;
            v_frag.scale += v_scale;
            v_frag.opacity += v_opacity;
            v_frag.rgb += v_rgb;
            v_ray_o += v_ray_o_t;
            v_ray_d += v_ray_d_t;
        }

    };

#endif  // #ifdef __CUDACC__

};
