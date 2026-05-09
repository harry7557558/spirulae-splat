#pragma once

#ifdef __CUDACC__
#include "generated/slang.cuh"
namespace Slang3DGSSV {
#include "generated/set_namespace.cuh"
#include "generated/primitive_3dgs_sv.cuh"
}
#endif

#include "PrimitiveBase3DGS.cuh"

template<int num_sv>
struct SphericalVoronoi3DGUT : public _BasePrimitive3DGS {
    static constexpr RenderOutputType pixelType = RenderOutputType::RGB_D;
};


#if 0
template<int num_sv>
struct SphericalVoronoi3DGUT : public Base3DGUT {
    struct World;

#ifdef __CUDACC__

    inline static __device__ void project_persp(
        World world, ProjCamera cam,
        Screen& proj, float4& aabb
    );

    inline static __device__ void project_fisheye(
        World world, ProjCamera cam,
        Screen& proj, float4& aabb
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


template<int num_sv>
struct SphericalVoronoi3DGUT<num_sv>::World : public Base3DGUT::World {

#ifdef __CUDACC__
    FixedArray<float3, num_sv> sv_sites;
    FixedArray<float3, num_sv> sv_colors;
#endif

    #ifndef NO_TORCH
    typedef std::tuple<at::Tensor, at::Tensor, at::Tensor, at::Tensor, at::Tensor, at::Tensor> TensorTuple;
    #endif

    struct Buffer;

    #ifndef NO_TORCH
    struct Tensor {
        at::Tensor means;
        at::Tensor quats;
        at::Tensor scales;
        at::Tensor opacities;
        at::Tensor sv_sites;
        at::Tensor sv_colors;

        Tensor() {}

        Tensor(const TensorTuple& splats) {
            means = std::get<0>(splats);
            quats = std::get<1>(splats);
            scales = std::get<2>(splats);
            opacities = std::get<3>(splats);
            sv_sites = std::get<4>(splats);
            sv_colors = std::get<5>(splats);
        }

        TensorTuple tuple() const {
            return std::make_tuple(means, quats, scales, opacities, sv_sites, sv_colors);
        }

        Tensor allocProjBwd(bool is_hess_diag) const {
            return Tensor(std::make_tuple(
                zeros_like_tensor(means),
                zeros_like_tensor(quats),
                zeros_like_tensor(scales),
                zeros_like_tensor(opacities),
                zeros_like_tensor(sv_sites),
                zeros_like_tensor(sv_colors)
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
    #endif

    struct Buffer : public TensorArray<6> {
        using TensorArray<6>::TensorArray;

        Buffer() {}

        #ifndef NO_TORCH
        Buffer(const Tensor& tensors)
            : TensorArray<6>(std::vector<std::optional<at::Tensor>>{
                tensors.means, tensors.quats, tensors.scales,
                tensors.opacities, tensors.sv_sites, tensors.sv_colors
            }) {}
        #endif

        #ifdef __CUDACC__
        __forceinline__ __device__ float3& means(int64_t i)
            { return *reinterpret_cast<float3*>(&_data[0][3*i]); }
        __forceinline__ __device__ float4& quats(int64_t i)
            { return *reinterpret_cast<float4*>(&_data[1][4*i]); }
        __forceinline__ __device__ float3& scales(int64_t i)
            { return *reinterpret_cast<float3*>(&_data[2][3*i]); }
        __forceinline__ __device__ float& opacities(int64_t i)
            { return *reinterpret_cast<float*>(&_data[3][i]); }
        __forceinline__ __device__ float3& sv_sites(int64_t i)
            { return *reinterpret_cast<float3*>(&_data[4][3*i]); }
        __forceinline__ __device__ float3& sv_colors(int64_t i)
            { return *reinterpret_cast<float3*>(&_data[5][3*i]); }
        #endif
    };

#ifdef __CUDACC__

    static __device__ World load(const Buffer &buffer, long idx) {
        World world = {
            buffer.means(idx),
            buffer.quats(idx),
            buffer.scales(idx),
            buffer.opacities(idx)
        };
        for (int i = 0; i < num_sv; i++) {
            world.sv_sites[i] = buffer.sv_sites(idx*num_sv+i);
            world.sv_colors[i] = buffer.sv_colors(idx*num_sv+i);
        }
        return world;
    }

    static __device__ __forceinline__ World zero() {
        World world;
        world.mean = {0.f, 0.f, 0.f};
        world.quat = {0.f, 0.f, 0.f, 0.f};
        world.scale = {0.f, 0.f, 0.f};
        world.opacity = 0.f;
        for (int i = 0; i < num_sv; i++) {
            world.sv_sites[i] = make_float3(0);
            world.sv_colors[i] = make_float3(0);
        }
        return world;
    }

    __device__ void saveParamsToBuffer(Buffer &buffer, long idx) {
        buffer.means(idx) = mean;
        buffer.quats(idx) = quat;
        buffer.scales(idx) = scale;
        buffer.opacities(idx) = opacity;
        for (int i = 0; i < num_sv; i++) {
            buffer.sv_sites(idx*num_sv+i) = sv_sites[i];
            buffer.sv_colors(idx*num_sv+i) = sv_colors[i];
        }
    }

    __device__ void atomicAddGradientToBuffer(Buffer &buffer, long idx) {
        atomicAddFVec(buffer.means(idx), mean);
        atomicAddFVec(buffer.quats(idx), quat);
        atomicAddFVec(buffer.scales(idx), scale);
        atomicAddFVec(buffer.opacities(idx), opacity);
        for (int i = 0; i < num_sv; i++) {
            atomicAddFVec(buffer.sv_sites(idx*num_sv+i), sv_sites[i]);
            atomicAddFVec(buffer.sv_colors(idx*num_sv+i), sv_colors[i]);
        }
    }

#endif  // #ifdef __CUDACC__
};
#endif


#ifdef __CUDACC__

template<int num_sv>
inline __device__ void project_persp(
    typename SphericalVoronoi3DGUT<num_sv>::World world, ProjCamera cam,
    typename SphericalVoronoi3DGUT<num_sv>::Screen& proj, float4& aabb
) {
    float2 xy;
    Slang3DGSSV::projection_3dgut_sv_persp(
        false,
        world.mean, world.quat, world.scale, world.opacity, world.sv_sites, world.sv_colors,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
        cam.width, cam.height,
        &aabb, &xy, &proj.depth, &proj.scale, &proj.opacity, &proj.rgb
    );
    proj.mean = world.mean, proj.quat = world.quat;
}

template<int num_sv>
inline __device__ void project_fisheye(
    typename SphericalVoronoi3DGUT<num_sv>::World world, ProjCamera cam,
    typename SphericalVoronoi3DGUT<num_sv>::Screen& proj, float4& aabb
) {
    float2 xy;
    Slang3DGSSV::projection_3dgut_sv_fisheye(
        false,
        world.mean, world.quat, world.scale, world.opacity, world.sv_sites, world.sv_colors,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
        cam.width, cam.height,
        &aabb, &xy, &proj.depth, &proj.scale, &proj.opacity, &proj.rgb
    );
    proj.mean = world.mean, proj.quat = world.quat;
}

template<int num_sv>
inline __device__ void project_persp_vjp(
    typename SphericalVoronoi3DGUT<num_sv>::World world, ProjCamera cam,
    typename SphericalVoronoi3DGUT<num_sv>::Screen v_proj,
    typename SphericalVoronoi3DGUT<num_sv>::World& v_world, float3x3 &v_R, float3 &v_t
) {
    Slang3DGSSV::projection_3dgut_sv_persp_vjp(
        false,
        world.mean, world.quat, world.scale, world.opacity, world.sv_sites, world.sv_colors,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
        cam.width, cam.height,
        make_float2(0), v_proj.depth, v_proj.scale, v_proj.opacity, v_proj.rgb,
        &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity, &v_world.sv_sites, &v_world.sv_colors,
        &v_R, &v_t
    );
    v_world.mean = v_proj.mean, v_world.quat = v_proj.quat;
}

template<int num_sv>
inline __device__ void project_persp_vjp(
    typename SphericalVoronoi3DGUT<num_sv>::World world, ProjCamera cam,
    typename SphericalVoronoi3DGUT<num_sv>::Screen v_proj, typename SphericalVoronoi3DGUT<num_sv>::Screen vr_proj, typename SphericalVoronoi3DGUT<num_sv>::Screen h_proj,
    typename SphericalVoronoi3DGUT<num_sv>::World& v_world, float3x3 &v_R, float3 &v_t,
    float3 &vr_world_pos, float3 &h_world_pos
) {}  // TODO

template<int num_sv>
inline __device__ void project_persp_vjp(
    typename SphericalVoronoi3DGUT<num_sv>::World world, ProjCamera cam,
    typename SphericalVoronoi3DGUT<num_sv>::Screen v_proj, typename SphericalVoronoi3DGUT<num_sv>::Screen vr_proj, typename SphericalVoronoi3DGUT<num_sv>::Screen h_proj,
    typename SphericalVoronoi3DGUT<num_sv>::World& v_world, float3x3 &v_R, float3 &v_t,
    typename SphericalVoronoi3DGUT<num_sv>::World& vr_world, typename SphericalVoronoi3DGUT<num_sv>::World& h_world
) {}  // TODO

template<int num_sv>
inline __device__ void project_fisheye_vjp(
    typename SphericalVoronoi3DGUT<num_sv>::World world, ProjCamera cam,
    typename SphericalVoronoi3DGUT<num_sv>::Screen v_proj,
    typename SphericalVoronoi3DGUT<num_sv>::World& v_world, float3x3 &v_R, float3 &v_t
) {
    Slang3DGSSV::projection_3dgut_sv_fisheye_vjp(
        false,
        world.mean, world.quat, world.scale, world.opacity, world.sv_sites, world.sv_colors,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
        cam.width, cam.height,
        make_float2(0), v_proj.depth, v_proj.scale, v_proj.opacity, v_proj.rgb,
        &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity, &v_world.sv_sites, &v_world.sv_colors,
        &v_R, &v_t
    );
    v_world.mean = v_proj.mean, v_world.quat = v_proj.quat;
}

template<int num_sv>
inline __device__ void project_fisheye_vjp(
    typename SphericalVoronoi3DGUT<num_sv>::World world, ProjCamera cam,
    typename SphericalVoronoi3DGUT<num_sv>::Screen v_proj, typename SphericalVoronoi3DGUT<num_sv>::Screen vr_proj, typename SphericalVoronoi3DGUT<num_sv>::Screen h_proj,
    typename SphericalVoronoi3DGUT<num_sv>::World& v_world, float3x3 &v_R, float3 &v_t,
    float3 &vr_world_pos, float3 &h_world_pos
) {}  // TODO

template<int num_sv>
inline __device__ void project_fisheye_vjp(
    typename SphericalVoronoi3DGUT<num_sv>::World world, ProjCamera cam,
    typename SphericalVoronoi3DGUT<num_sv>::Screen v_proj, typename SphericalVoronoi3DGUT<num_sv>::Screen vr_proj, typename SphericalVoronoi3DGUT<num_sv>::Screen h_proj,
    typename SphericalVoronoi3DGUT<num_sv>::World& v_world, float3x3 &v_R, float3 &v_t,
    typename SphericalVoronoi3DGUT<num_sv>::World& vr_world, typename SphericalVoronoi3DGUT<num_sv>::World& h_world
) {}  // TODO

#endif  // #ifdef __CUDACC__

typedef SphericalVoronoi3DGUT<2> SphericalVoronoi3DGUT_Default;
