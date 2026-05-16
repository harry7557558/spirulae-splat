#pragma once

#ifdef __CUDACC__
#include "generated/slang.cuh"
namespace Slang3DGS {
#include "generated/set_namespace.cuh"
#include "generated/primitive_3dgs.cuh"
}
#endif

#include "Primitive.cuh"
#include "PrimitiveBase3DGS.cuh"


template<bool antialiased>
struct _Base3DGS : public _BasePrimitive3DGS {
    static constexpr RenderOutputType pixelType = RenderOutputType::RGB_D;

    class WorldBuffer : public _BasePrimitive3DGS::WorldBuffer {
        using _BasePrimitive3DGS::WorldBuffer::WorldBuffer;
    };
    class ScreenBuffer : public _BasePrimitive3DGS::ScreenBuffer {
        using _BasePrimitive3DGS::ScreenBuffer::ScreenBuffer;
    };

#ifdef __CUDACC__

    struct World : public _BasePrimitive3DGS::World{
        __device__ World() = default;
        __device__ World(const _BasePrimitive3DGS::World& other) {
            _BasePrimitive3DGS::World::operator=(other);
        }
        __device__ World& operator=(const _BasePrimitive3DGS::World& other) {
            _BasePrimitive3DGS::World::operator=(other);
            return *this;
        }

        template<ssplat::CameraModelType camera_model>
        inline __device__ void project(
            ProjCamera cam,
            typename _Base3DGS<antialiased>::Screen& screen, float4& aabb, float& sorting_depth, float& radius
        ) const {
            if constexpr (camera_model == ssplat::CameraModelType::PINHOLE)
                Slang3DGS::projection_3dgs_persp(
                    antialiased,
                    mean, quat, scale, opacity, sh_coeffs,
                    cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy,
                    cam.dist_coeffs,
                    cam.width, cam.height,
                    &aabb, &sorting_depth, &radius, &screen.xy, &screen.depth, &screen.conic, &screen.opac, &screen.rgb
                );
            else if constexpr (camera_model == ssplat::CameraModelType::FISHEYE)
                Slang3DGS::projection_3dgs_fisheye(
                    antialiased,
                    mean, quat, scale, opacity, sh_coeffs,
                    cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
                    cam.width, cam.height,
                    &aabb, &sorting_depth, &radius, &screen.xy, &screen.depth, &screen.conic, &screen.opac, &screen.rgb
                );
            else if constexpr (camera_model == ssplat::CameraModelType::EQUISOLID)
                Slang3DGS::projection_3dgs_equisolid(
                    antialiased,
                    mean, quat, scale, opacity, sh_coeffs,
                    cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
                    cam.width, cam.height,
                    &aabb, &sorting_depth, &radius, &screen.xy, &screen.depth, &screen.conic, &screen.opac, &screen.rgb
                );
        }

        template<ssplat::CameraModelType camera_model>
        inline __device__ void project_vjp(
            ProjCamera cam,
            typename _Base3DGS<antialiased>::Screen v_screen,
            typename _Base3DGS<antialiased>::World& v_world, float3x3 &v_R, float3 &v_t
        ) const {
            if constexpr (camera_model == ssplat::CameraModelType::PINHOLE)
                Slang3DGS::projection_3dgs_persp_vjp(
                    antialiased,
                    mean, quat, scale, opacity, sh_coeffs,
                    cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
                    cam.width, cam.height,
                    v_screen.xy, v_screen.depth, v_screen.conic, v_screen.opac, v_screen.rgb,
                    &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity, &v_world.sh_coeffs,
                    &v_R, &v_t
                );
            else if constexpr (camera_model == ssplat::CameraModelType::FISHEYE)
                Slang3DGS::projection_3dgs_fisheye_vjp(
                    antialiased,
                    mean, quat, scale, opacity, sh_coeffs,
                    cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
                    cam.width, cam.height,
                    v_screen.xy, v_screen.depth, v_screen.conic, v_screen.opac, v_screen.rgb,
                    &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity, &v_world.sh_coeffs,
                    &v_R, &v_t
                );
            else if constexpr (camera_model == ssplat::CameraModelType::EQUISOLID)
                Slang3DGS::projection_3dgs_equisolid_vjp(
                    antialiased,
                    mean, quat, scale, opacity, sh_coeffs,
                    cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
                    cam.width, cam.height,
                    v_screen.xy, v_screen.depth, v_screen.conic, v_screen.opac, v_screen.rgb,
                    &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity, &v_world.sh_coeffs,
                    &v_R, &v_t
                );
        }

        template<ssplat::CameraModelType camera_model>
        inline __device__ void project_vjp_h_pos(
            ProjCamera cam,
            typename _Base3DGS<antialiased>::Screen v_proj, typename _Base3DGS<antialiased>::Screen vr_proj, typename _Base3DGS<antialiased>::Screen h_proj,
            typename _Base3DGS<antialiased>::World& v_world, float3x3 &v_R, float3 &v_t,
            float3 &vr_world_pos, float3 &h_world_pos
        ) const {
            if constexpr (camera_model == ssplat::CameraModelType::PINHOLE)
                Slang3DGS::projection_3dgs_persp_vjp(
                    antialiased,
                    mean, quat, scale, opacity, sh_coeffs,
                    cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
                    cam.width, cam.height,
                    v_proj.xy, v_proj.depth, v_proj.conic, v_proj.opac, v_proj.rgb,
                    vr_proj.xy, vr_proj.depth, vr_proj.conic, vr_proj.opac, vr_proj.rgb,
                    h_proj.xy, h_proj.depth, h_proj.conic, h_proj.opac, h_proj.rgb,
                    &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity, &v_world.sh_coeffs,
                    &v_R, &v_t, &vr_world_pos, &h_world_pos
                );
            else if constexpr (camera_model == ssplat::CameraModelType::FISHEYE)
                Slang3DGS::projection_3dgs_fisheye_vjp(
                    antialiased,
                    mean, quat, scale, opacity, sh_coeffs,
                    cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
                    cam.width, cam.height,
                    v_proj.xy, v_proj.depth, v_proj.conic, v_proj.opac, v_proj.rgb,
                    vr_proj.xy, vr_proj.depth, vr_proj.conic, vr_proj.opac, vr_proj.rgb,
                    h_proj.xy, h_proj.depth, h_proj.conic, h_proj.opac, h_proj.rgb,
                    &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity, &v_world.sh_coeffs,
                    &v_R, &v_t, &vr_world_pos, &h_world_pos
                );
            else if constexpr (camera_model == ssplat::CameraModelType::EQUISOLID)
                Slang3DGS::projection_3dgs_equisolid_vjp(
                    antialiased,
                    mean, quat, scale, opacity, sh_coeffs,
                    cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
                    cam.width, cam.height,
                    v_proj.xy, v_proj.depth, v_proj.conic, v_proj.opac, v_proj.rgb,
                    vr_proj.xy, vr_proj.depth, vr_proj.conic, vr_proj.opac, vr_proj.rgb,
                    h_proj.xy, h_proj.depth, h_proj.conic, h_proj.opac, h_proj.rgb,
                    &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity, &v_world.sh_coeffs,
                    &v_R, &v_t, &vr_world_pos, &h_world_pos
                );
        }

        template<ssplat::CameraModelType camera_model>
        inline __device__ void project_vjp_h_all(
            ProjCamera cam,
            typename _Base3DGS<antialiased>::Screen v_proj, typename _Base3DGS<antialiased>::Screen vr_proj, typename _Base3DGS<antialiased>::Screen h_proj,
            typename _Base3DGS<antialiased>::World& v_world, float3x3 &v_R, float3 &v_t,
            typename _Base3DGS<antialiased>::World& vr_world, typename _Base3DGS<antialiased>::World& h_world
        ) const {
            if constexpr (camera_model == ssplat::CameraModelType::PINHOLE)
                Slang3DGS::projection_3dgs_persp_vjp(
                    antialiased,
                    mean, quat, scale, opacity, sh_coeffs,
                    cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
                    cam.width, cam.height,
                    v_proj.xy, v_proj.depth, v_proj.conic, v_proj.opac, v_proj.rgb,
                    vr_proj.xy, vr_proj.depth, vr_proj.conic, vr_proj.opac, vr_proj.rgb,
                    h_proj.xy, h_proj.depth, h_proj.conic, h_proj.opac, h_proj.rgb,
                    &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity, &v_world.sh_coeffs,
                    &v_R, &v_t,
                    &vr_world.mean, &vr_world.quat, &vr_world.scale, &vr_world.opacity, &vr_world.sh_coeffs[0],
                    &h_world.mean, &h_world.quat, &h_world.scale, &h_world.opacity, &h_world.sh_coeffs[0]
                );
            else if constexpr (camera_model == ssplat::CameraModelType::FISHEYE)
                Slang3DGS::projection_3dgs_fisheye_vjp(
                    antialiased,
                    mean, quat, scale, opacity, sh_coeffs,
                    cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
                    cam.width, cam.height,
                    v_proj.xy, v_proj.depth, v_proj.conic, v_proj.opac, v_proj.rgb,
                    vr_proj.xy, vr_proj.depth, vr_proj.conic, vr_proj.opac, vr_proj.rgb,
                    h_proj.xy, h_proj.depth, h_proj.conic, h_proj.opac, h_proj.rgb,
                    &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity, &v_world.sh_coeffs,
                    &v_R, &v_t,
                    &vr_world.mean, &vr_world.quat, &vr_world.scale, &vr_world.opacity, &vr_world.sh_coeffs[0],
                    &h_world.mean, &h_world.quat, &h_world.scale, &h_world.opacity, &h_world.sh_coeffs[0]
                );
            else if constexpr (camera_model == ssplat::CameraModelType::EQUISOLID)
                Slang3DGS::projection_3dgs_equisolid_vjp(
                    antialiased,
                    mean, quat, scale, opacity, sh_coeffs,
                    cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
                    cam.width, cam.height,
                    v_proj.xy, v_proj.depth, v_proj.conic, v_proj.opac, v_proj.rgb,
                    vr_proj.xy, vr_proj.depth, vr_proj.conic, vr_proj.opac, vr_proj.rgb,
                    h_proj.xy, h_proj.depth, h_proj.conic, h_proj.opac, h_proj.rgb,
                    &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity, &v_world.sh_coeffs,
                    &v_R, &v_t,
                    &vr_world.mean, &vr_world.quat, &vr_world.scale, &vr_world.opacity, &vr_world.sh_coeffs[0],
                    &h_world.mean, &h_world.quat, &h_world.scale, &h_world.opacity, &h_world.sh_coeffs[0]
                );
        }

    };

    struct Screen : public _BasePrimitive3DGS::Screen {
        __device__ Screen() = default;
        __device__ Screen(const _BasePrimitive3DGS::Screen& other) {
            _BasePrimitive3DGS::Screen::operator=(other);
        }
        __device__ Screen& operator=(const _BasePrimitive3DGS::Screen& other) {
            _BasePrimitive3DGS::Screen::operator=(other);
            return *this;
        }
    };

    struct FragmentFwd : public _BasePrimitive3DGS::Screen {
        __device__ FragmentFwd() = default;
        __device__ FragmentFwd(const _BasePrimitive3DGS::Screen& other) {
            _BasePrimitive3DGS::Screen::operator=(other);
        }
        __device__ FragmentFwd& operator=(const _BasePrimitive3DGS::Screen& other) {
            _BasePrimitive3DGS::Screen::operator=(other);
            return *this;
        }

        __device__ __forceinline__ void load(
            const WorldBuffer &wbuffer,
            const ScreenBuffer &sbuffer,
            int64_t wi, int64_t si
        ) {
            _BasePrimitive3DGS::Screen::load(sbuffer, si);
        }

        __device__ __forceinline__ void store(
            WorldBuffer &wbuffer,
            ScreenBuffer &sbuffer,
            int64_t wi, int64_t si
        ) const {
            _BasePrimitive3DGS::Screen::store(sbuffer, si);
        }

        __device__ __forceinline__ void atomicStore(
            WorldBuffer &wbuffer,
            ScreenBuffer &sbuffer,
            int64_t wi, int64_t si
        ) const {
            _BasePrimitive3DGS::Screen::atomicStore(sbuffer, si);
        }

        __device__ __forceinline__ float evaluate_alpha(
            float px, float py
        ) const {
            float2 delta = {xy.x - px, xy.y - py};
            float sigma = 0.5f * (conic.x * delta.x * delta.x +
                                    conic.z * delta.y * delta.y) +
                            conic.y * delta.x * delta.y;
            float vis = __expf(-sigma);
            float alpha = min(0.999f, opac * vis);
            return sigma < 0.f ? 0.f : alpha;
        }

        __device__ __forceinline__ void evaluate_alpha_vjp(
            float px, float py, float v_alpha,
            FragmentFwd& v_frag
        ) const {
            float2 delta = {xy.x - px, xy.y - py};
            float sigma = 0.5f * (conic.x * delta.x * delta.x +
                                    conic.z * delta.y * delta.y) +
                            conic.y * delta.x * delta.y;
            float vis = __expf(-sigma);
            float alpha = min(0.999f, opac * vis);

            if (sigma >= 0.f && opac * vis <= 0.999f) {
                const float v_sigma = -opac * vis * v_alpha;
                v_frag.conic += float3{
                    0.5f * v_sigma * delta.x * delta.x,
                    v_sigma * delta.x * delta.y,
                    0.5f * v_sigma * delta.y * delta.y
                };
                v_frag.xy += float2{
                    v_sigma * (conic.x * delta.x + conic.y * delta.y),
                    v_sigma * (conic.y * delta.x + conic.z * delta.y)
                };
                v_frag.opac += vis * v_alpha;
            }
        }

        __device__ __forceinline__ RenderOutput evaluate_color(
            float px, float py
        ) const {
            return { rgb, depth, make_float3(0.0f) };
        }

        __device__ __forceinline__ void evaluate_color_vjp(
            float px, float py, RenderOutput v_render,
            FragmentFwd& v_frag
        ) const {
            v_frag.rgb += v_render.rgb;
            v_frag.depth += v_render.depth;
        }

    };

    struct FragmentBwd : public FragmentFwd {
        __device__ FragmentBwd() = default;
        __device__ FragmentBwd(const _BasePrimitive3DGS::Screen& other) {
            _BasePrimitive3DGS::Screen::operator=(other);
        }
        __device__ FragmentBwd& operator=(const _BasePrimitive3DGS::Screen& other) {
            _BasePrimitive3DGS::Screen::operator=(other);
            return *this;
        }

        static __device__ __forceinline__ FragmentBwd zero(const FragmentBwd& bwd) {
            return _BasePrimitive3DGS::Screen::zero();
        }
    };

#endif  // #ifdef __CUDACC__

};

typedef _Base3DGS<false> Vanilla3DGS;
typedef _Base3DGS<true> MipSplatting;


