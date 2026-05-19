#pragma once

#ifdef __CUDACC__
#include "generated/slang.cuh"
namespace Slang3DGS {
#include "generated/set_namespace.cuh"
#include "generated/primitive_3dgs.cuh"
}
namespace SlangHarmonics {
#include "generated/set_namespace.cuh"
#include "generated/harmonics.cuh"
}
#endif

#include "Primitive.cuh"
#include "PrimitiveBase3DGS.cuh"


template<int sh_degree, bool antialiased>
struct _Base3DGS : public _BasePrimitive3DGS<sh_degree> {
    static constexpr RenderOutputType pixelType = RenderOutputType::RGB_D;

    class WorldBuffer : public _BasePrimitive3DGS<sh_degree>::WorldBuffer {
        using _BasePrimitive3DGS<sh_degree>::WorldBuffer::WorldBuffer;
    };
    class ScreenBuffer : public _BasePrimitive3DGS<sh_degree>::ScreenBuffer {
        using _BasePrimitive3DGS<sh_degree>::ScreenBuffer::ScreenBuffer;
    };

#ifdef __CUDACC__

    struct World : public _BasePrimitive3DGS<sh_degree>::World{
        __device__ World() = default;
        __device__ World(const typename _BasePrimitive3DGS<sh_degree>::World& other) {
            _BasePrimitive3DGS<sh_degree>::World::operator=(other);
        }
        __device__ World& operator=(const typename _BasePrimitive3DGS<sh_degree>::World& other) {
            _BasePrimitive3DGS<sh_degree>::World::operator=(other);
            return *this;
        }

        template<ssplat::CameraModelType camera_model>
        inline __device__ void project(
            ProjCamera cam,
            typename _Base3DGS<sh_degree, antialiased>::Screen& screen, float4& aabb, float& sorting_depth, float& radius
        ) {
            if constexpr (camera_model == ssplat::CameraModelType::PINHOLE)
                Slang3DGS::projection_3dgs_persp(
                    antialiased,
                    this->mean, this->quat, this->scale, this->opacity,
                    cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy,
                    cam.dist_coeffs,
                    cam.width, cam.height,
                    &aabb, &sorting_depth, &radius, &screen.xy, &screen.depth, &screen.conic, &screen.opac
                );
            else if constexpr (camera_model == ssplat::CameraModelType::FISHEYE)
                Slang3DGS::projection_3dgs_fisheye(
                    antialiased,
                    this->mean, this->quat, this->scale, this->opacity,
                    cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
                    cam.width, cam.height,
                    &aabb, &sorting_depth, &radius, &screen.xy, &screen.depth, &screen.conic, &screen.opac
                );
            else if constexpr (camera_model == ssplat::CameraModelType::EQUISOLID)
                Slang3DGS::projection_3dgs_equisolid(
                    antialiased,
                    this->mean, this->quat, this->scale, this->opacity,
                    cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
                    cam.width, cam.height,
                    &aabb, &sorting_depth, &radius, &screen.xy, &screen.depth, &screen.conic, &screen.opac
                );
            if (aabb.z > aabb.x && aabb.w > aabb.y) {
                // screen.rgb = SlangHarmonics::sh_coeffs_to_color(
                //     sh_degree, this->mean, cam.R, cam.t, this->features_dc, this->features_sh
                // );
                if constexpr (sh_degree == 0) screen.rgb = SlangHarmonics::sh0_to_color
                    (this->mean, cam.R, cam.t, this->features_dc, (float3*)this->features_sh);
                else if constexpr (sh_degree == 1) screen.rgb = SlangHarmonics::sh1_to_color
                    (this->mean, cam.R, cam.t, this->features_dc, (float3*)this->features_sh);
                else if constexpr (sh_degree == 2) screen.rgb = SlangHarmonics::sh2_to_color
                    (this->mean, cam.R, cam.t, this->features_dc, (float3*)this->features_sh);
                else if constexpr (sh_degree == 3) screen.rgb = SlangHarmonics::sh3_to_color
                    (this->mean, cam.R, cam.t, this->features_dc, (float3*)this->features_sh);
                else if constexpr (sh_degree == 4) screen.rgb = SlangHarmonics::sh4_to_color
                    (this->mean, cam.R, cam.t, this->features_dc, (float3*)this->features_sh);
            }
        }

        template<ssplat::CameraModelType camera_model, bool atomic=true>
        inline __device__ void project_vjp(
            ProjCamera cam,
            typename _Base3DGS<sh_degree, antialiased>::Screen v_screen,
            typename _Base3DGS<sh_degree, antialiased>::World& v_world, float3x3 &v_R, float3 &v_t
        ) {
            if constexpr (camera_model == ssplat::CameraModelType::PINHOLE)
                Slang3DGS::projection_3dgs_persp_vjp(
                    antialiased,
                    this->mean, this->quat, this->scale, this->opacity,
                    cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
                    cam.width, cam.height,
                    v_screen.xy, v_screen.depth, v_screen.conic, v_screen.opac,
                    &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity,
                    &v_R, &v_t
                );
            else if constexpr (camera_model == ssplat::CameraModelType::FISHEYE)
                Slang3DGS::projection_3dgs_fisheye_vjp(
                    antialiased,
                    this->mean, this->quat, this->scale, this->opacity,
                    cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
                    cam.width, cam.height,
                    v_screen.xy, v_screen.depth, v_screen.conic, v_screen.opac,
                    &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity,
                    &v_R, &v_t
                );
            else if constexpr (camera_model == ssplat::CameraModelType::EQUISOLID)
                Slang3DGS::projection_3dgs_equisolid_vjp(
                    antialiased,
                    this->mean, this->quat, this->scale, this->opacity,
                    cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
                    cam.width, cam.height,
                    v_screen.xy, v_screen.depth, v_screen.conic, v_screen.opac,
                    &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity,
                    &v_R, &v_t
                );
            // SH: atomic for global memory, add for local/shared memory
            #define _ARGS ( \
                this->mean, cam.R, cam.t, this->features_dc, (float3*)this->features_sh, \
                v_screen.rgb, &v_world.features_dc, (float3*)v_world.features_sh, \
                &v_world.mean, &v_R, &v_t \
            );
            if constexpr (atomic) {
                if constexpr (sh_degree == 0) SlangHarmonics::sh0_to_color_vjp_atomic _ARGS
                else if constexpr (sh_degree == 1) SlangHarmonics::sh1_to_color_vjp_atomic _ARGS
                else if constexpr (sh_degree == 2) SlangHarmonics::sh2_to_color_vjp_atomic _ARGS
                else if constexpr (sh_degree == 3) SlangHarmonics::sh3_to_color_vjp_atomic _ARGS
                else if constexpr (sh_degree == 4) SlangHarmonics::sh4_to_color_vjp_atomic _ARGS
            }
            else {
                if constexpr (sh_degree == 0) SlangHarmonics::sh0_to_color_vjp_inplace _ARGS
                else if constexpr (sh_degree == 1) SlangHarmonics::sh1_to_color_vjp_inplace _ARGS
                else if constexpr (sh_degree == 2) SlangHarmonics::sh2_to_color_vjp_inplace _ARGS
                else if constexpr (sh_degree == 3) SlangHarmonics::sh3_to_color_vjp_inplace _ARGS
                else if constexpr (sh_degree == 4) SlangHarmonics::sh4_to_color_vjp_inplace _ARGS
            }
            #undef _ARGS
        }

        template<ssplat::CameraModelType camera_model, bool atomic=true>
        inline __device__ void project_vjp_h_pos(
            ProjCamera cam,
            typename _Base3DGS<sh_degree, antialiased>::Screen v_proj,
            typename _Base3DGS<sh_degree, antialiased>::Screen vr_proj,
            typename _Base3DGS<sh_degree, antialiased>::Screen h_proj,
            typename _Base3DGS<sh_degree, antialiased>::World& v_world, float3x3 &v_R, float3 &v_t,
            float3 &vr_world_pos, float3 &h_world_pos
        ) {
            if constexpr (camera_model == ssplat::CameraModelType::PINHOLE)
                Slang3DGS::projection_3dgs_persp_vjp(
                    antialiased,
                    this->mean, this->quat, this->scale, this->opacity,
                    cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
                    cam.width, cam.height,
                    v_proj.xy, v_proj.depth, v_proj.conic, v_proj.opac,
                    vr_proj.xy, vr_proj.depth, vr_proj.conic, vr_proj.opac,
                    h_proj.xy, h_proj.depth, h_proj.conic, h_proj.opac,
                    &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity,
                    &v_R, &v_t, &vr_world_pos, &h_world_pos
                );
            else if constexpr (camera_model == ssplat::CameraModelType::FISHEYE)
                Slang3DGS::projection_3dgs_fisheye_vjp(
                    antialiased,
                    this->mean, this->quat, this->scale, this->opacity,
                    cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
                    cam.width, cam.height,
                    v_proj.xy, v_proj.depth, v_proj.conic, v_proj.opac,
                    vr_proj.xy, vr_proj.depth, vr_proj.conic, vr_proj.opac,
                    h_proj.xy, h_proj.depth, h_proj.conic, h_proj.opac,
                    &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity,
                    &v_R, &v_t, &vr_world_pos, &h_world_pos
                );
            else if constexpr (camera_model == ssplat::CameraModelType::EQUISOLID)
                Slang3DGS::projection_3dgs_equisolid_vjp(
                    antialiased,
                    this->mean, this->quat, this->scale, this->opacity,
                    cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
                    cam.width, cam.height,
                    v_proj.xy, v_proj.depth, v_proj.conic, v_proj.opac,
                    vr_proj.xy, vr_proj.depth, vr_proj.conic, vr_proj.opac,
                    h_proj.xy, h_proj.depth, h_proj.conic, h_proj.opac,
                    &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity,
                    &v_R, &v_t, &vr_world_pos, &h_world_pos
                );
            // TODO: features_dc
        }

        template<ssplat::CameraModelType camera_model, bool atomic=true>
        inline __device__ void project_vjp_h_all(
            ProjCamera cam,
            typename _Base3DGS<sh_degree, antialiased>::Screen v_proj,
            typename _Base3DGS<sh_degree, antialiased>::Screen vr_proj,
            typename _Base3DGS<sh_degree, antialiased>::Screen h_proj,
            typename _Base3DGS<sh_degree, antialiased>::World& v_world, float3x3 &v_R, float3 &v_t,
            typename _Base3DGS<sh_degree, antialiased>::World& vr_world,
            typename _Base3DGS<sh_degree, antialiased>::World& h_world
        ) {
            if constexpr (camera_model == ssplat::CameraModelType::PINHOLE)
                Slang3DGS::projection_3dgs_persp_vjp(
                    antialiased,
                    this->mean, this->quat, this->scale, this->opacity,
                    cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
                    cam.width, cam.height,
                    v_proj.xy, v_proj.depth, v_proj.conic, v_proj.opac,
                    vr_proj.xy, vr_proj.depth, vr_proj.conic, vr_proj.opac,
                    h_proj.xy, h_proj.depth, h_proj.conic, h_proj.opac,
                    &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity,
                    &v_R, &v_t,
                    &vr_world.mean, &vr_world.quat, &vr_world.scale, &vr_world.opacity,
                    &h_world.mean, &h_world.quat, &h_world.scale, &h_world.opacity
                );
            else if constexpr (camera_model == ssplat::CameraModelType::FISHEYE)
                Slang3DGS::projection_3dgs_fisheye_vjp(
                    antialiased,
                    this->mean, this->quat, this->scale, this->opacity,
                    cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
                    cam.width, cam.height,
                    v_proj.xy, v_proj.depth, v_proj.conic, v_proj.opac,
                    vr_proj.xy, vr_proj.depth, vr_proj.conic, vr_proj.opac,
                    h_proj.xy, h_proj.depth, h_proj.conic, h_proj.opac,
                    &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity,
                    &v_R, &v_t,
                    &vr_world.mean, &vr_world.quat, &vr_world.scale, &vr_world.opacity,
                    &h_world.mean, &h_world.quat, &h_world.scale, &h_world.opacity
                );
            else if constexpr (camera_model == ssplat::CameraModelType::EQUISOLID)
                Slang3DGS::projection_3dgs_equisolid_vjp(
                    antialiased,
                    this->mean, this->quat, this->scale, this->opacity,
                    cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
                    cam.width, cam.height,
                    v_proj.xy, v_proj.depth, v_proj.conic, v_proj.opac,
                    vr_proj.xy, vr_proj.depth, vr_proj.conic, vr_proj.opac,
                    h_proj.xy, h_proj.depth, h_proj.conic, h_proj.opac,
                    &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity,
                    &v_R, &v_t,
                    &vr_world.mean, &vr_world.quat, &vr_world.scale, &vr_world.opacity,
                    &h_world.mean, &h_world.quat, &h_world.scale, &h_world.opacity
                );
            // TODO: features_dc
        }

    };

    struct Screen : public _BasePrimitive3DGS<sh_degree>::Screen {
        __device__ Screen() = default;
        __device__ Screen(const typename _BasePrimitive3DGS<sh_degree>::Screen& other) {
            _BasePrimitive3DGS<sh_degree>::Screen::operator=(other);
        }
        __device__ Screen& operator=(const typename _BasePrimitive3DGS<sh_degree>::Screen& other) {
            _BasePrimitive3DGS<sh_degree>::Screen::operator=(other);
            return *this;
        }
    };

    struct FragmentFwd : public _BasePrimitive3DGS<sh_degree>::Screen {
        __device__ FragmentFwd() = default;
        __device__ FragmentFwd(const typename _BasePrimitive3DGS<sh_degree>::Screen& other) {
            _BasePrimitive3DGS<sh_degree>::Screen::operator=(other);
        }
        __device__ FragmentFwd& operator=(const typename _BasePrimitive3DGS<sh_degree>::Screen& other) {
            _BasePrimitive3DGS<sh_degree>::Screen::operator=(other);
            return *this;
        }

        __device__ __forceinline__ void load(
            const WorldBuffer &wbuffer,
            const ScreenBuffer &sbuffer,
            int64_t wi, int64_t si
        ) {
            _BasePrimitive3DGS<sh_degree>::Screen::load(sbuffer, si);
        }

        __device__ __forceinline__ void store(
            WorldBuffer &wbuffer,
            ScreenBuffer &sbuffer,
            int64_t wi, int64_t si
        ) const {
            _BasePrimitive3DGS<sh_degree>::Screen::store(sbuffer, si);
        }

        __device__ __forceinline__ void atomicStore(
            WorldBuffer &wbuffer,
            ScreenBuffer &sbuffer,
            int64_t wi, int64_t si
        ) const {
            _BasePrimitive3DGS<sh_degree>::Screen::atomicStore(sbuffer, si);
        }

        __device__ __forceinline__ float evaluate_alpha(
            float px, float py
        ) const {
            float2 delta = {this->xy.x - px, this->xy.y - py};
            float sigma = 0.5f * (this->conic.x * delta.x * delta.x +
                                    this->conic.z * delta.y * delta.y) +
                            this->conic.y * delta.x * delta.y;
            float vis = __expf(-sigma);
            float alpha = min(0.999f, this->opac * vis);
            return sigma < 0.f ? 0.f : alpha;
        }

        __device__ __forceinline__ void evaluate_alpha_vjp(
            float px, float py, float v_alpha,
            FragmentFwd& v_frag
        ) const {
            float2 delta = {this->xy.x - px, this->xy.y - py};
            float sigma = 0.5f * (this->conic.x * delta.x * delta.x +
                                    this->conic.z * delta.y * delta.y) +
                            this->conic.y * delta.x * delta.y;
            float vis = __expf(-sigma);
            float alpha = min(0.999f, this->opac * vis);

            if (sigma >= 0.f && this->opac * vis <= 0.999f) {
                const float v_sigma = -this->opac * vis * v_alpha;
                v_frag.conic += float3{
                    0.5f * v_sigma * delta.x * delta.x,
                    v_sigma * delta.x * delta.y,
                    0.5f * v_sigma * delta.y * delta.y
                };
                v_frag.xy += float2{
                    v_sigma * (this->conic.x * delta.x + this->conic.y * delta.y),
                    v_sigma * (this->conic.y * delta.x + this->conic.z * delta.y)
                };
                v_frag.opac += vis * v_alpha;
            }
        }

        __device__ __forceinline__ RenderOutput evaluate_color(
            float px, float py
        ) const {
            return { this->rgb, this->depth, make_float3(0.0f) };
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
        __device__ FragmentBwd(const typename _BasePrimitive3DGS<sh_degree>::Screen& other) {
            _BasePrimitive3DGS<sh_degree>::Screen::operator=(other);
        }
        __device__ FragmentBwd& operator=(const typename _BasePrimitive3DGS<sh_degree>::Screen& other) {
            _BasePrimitive3DGS<sh_degree>::Screen::operator=(other);
            return *this;
        }

        static __device__ __forceinline__ FragmentBwd zero(const FragmentBwd& bwd) {
            return _BasePrimitive3DGS<sh_degree>::Screen::zero();
        }
    };

#endif  // #ifdef __CUDACC__

};

template<int sh_degree>
using Vanilla3DGS = _Base3DGS<sh_degree, false>;

template<int sh_degree>
using MipSplatting = _Base3DGS<sh_degree, true>;
