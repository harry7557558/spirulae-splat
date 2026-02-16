#pragma once

#ifdef __CUDACC__
#include "generated/slang.cuh"
namespace Slang3DGSSV {
#include "generated/set_namespace.cuh"
#include "generated/primitive_3dgs_sv.cuh"
}
#endif

#include "PrimitiveBase3DGUT.cuh"


template<int num_sv>
struct SphericalVoronoi3DGUT : public Base3DGUT {
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


template<int num_sv>
struct SphericalVoronoi3DGUT<num_sv>::World : public Base3DGUT::World {

#ifdef __CUDACC__
    FixedArray<float3, num_sv> sv_sites;
    FixedArray<float3, num_sv> sv_colors;
#endif

    typedef std::tuple<at::Tensor, at::Tensor, at::Tensor, at::Tensor, at::Tensor, at::Tensor> TensorTuple;

    struct Buffer;

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
                at::zeros_like(means),
                at::zeros_like(quats),
                at::zeros_like(scales),
                at::zeros_like(opacities),
                at::zeros_like(sv_sites),
                at::zeros_like(sv_colors)
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
        float3* __restrict__ sv_sites;
        float3* __restrict__ sv_colors;

        Buffer() {}

        Buffer(const Tensor& tensors) {
            DEVICE_GUARD(tensors.means);
            CHECK_INPUT(tensors.means);
            CHECK_INPUT(tensors.quats);
            CHECK_INPUT(tensors.scales);
            CHECK_INPUT(tensors.opacities);
            CHECK_INPUT(tensors.sv_sites);
            CHECK_INPUT(tensors.sv_colors);
            means = (float3*)tensors.means.template data_ptr<float>();
            quats = (float4*)tensors.quats.template data_ptr<float>();
            scales = (float3*)tensors.scales.template data_ptr<float>();
            opacities = tensors.opacities.template data_ptr<float>();
            sv_sites = (float3*)tensors.sv_sites.template data_ptr<float>();
            sv_colors = (float3*)tensors.sv_colors.template data_ptr<float>();
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
        for (int i = 0; i < num_sv; i++) {
            world.sv_sites[i] = buffer.sv_sites[idx*num_sv+i];
            world.sv_colors[i] = buffer.sv_colors[idx*num_sv+i];
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

    template<typename Partition>
    __device__ void reduce(Partition& partition) {
        warpSum(mean, partition);
        warpSum(quat, partition);
        warpSum(scale, partition);
        warpSum(opacity, partition);
        for (int i = 0; i < num_sv; i++) {
            warpSum(sv_sites[i], partition);
            warpSum(sv_colors[i], partition);
        }
    }

    __device__ void saveParamsToBuffer(Buffer &buffer, long idx) {
        buffer.means[idx] = mean;
        buffer.quats[idx] = quat;
        buffer.scales[idx] = scale;
        buffer.opacities[idx] = opacity;
        for (int i = 0; i < num_sv; i++) {
            buffer.sv_sites[idx*num_sv+i] = sv_sites[i];
            buffer.sv_colors[idx*num_sv+i] = sv_colors[i];
        }
    }

    __device__ void atomicAddGradientToBuffer(Buffer &buffer, long idx) {
        atomicAddFVec(buffer.means + idx, mean);
        atomicAddFVec(buffer.quats + idx, quat);
        atomicAddFVec(buffer.scales + idx, scale);
        atomicAddFVec(buffer.opacities + idx, opacity);
        for (int i = 0; i < num_sv; i++) {
            atomicAddFVec(buffer.sv_sites + idx*num_sv + i, sv_sites[i]);
            atomicAddFVec(buffer.sv_colors + idx*num_sv + i, sv_colors[i]);
        }
    }

#endif  // #ifdef __CUDACC__
};


#ifdef __CUDACC__

template<int num_sv>
inline __device__ void SphericalVoronoi3DGUT<num_sv>::project_persp(
    SphericalVoronoi3DGUT<num_sv>::World world, SphericalVoronoi3DGUT<num_sv>::FwdProjCamera cam,
    SphericalVoronoi3DGUT<num_sv>::Screen& proj, int4& aabb
) {
    float2 xy;
    Slang3DGSSV::projection_3dgut_sv_persp(
        false,
        world.mean, world.quat, world.scale, world.opacity, world.sv_sites, world.sv_colors,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
        cam.width, cam.height, cam.near_plane, cam.far_plane,
        &aabb, &xy, &proj.depth, &proj.scale, &proj.opacity, &proj.rgb
    );
    proj.mean = world.mean, proj.quat = world.quat;
}

template<int num_sv>
inline __device__ void SphericalVoronoi3DGUT<num_sv>::project_fisheye(
    SphericalVoronoi3DGUT<num_sv>::World world, SphericalVoronoi3DGUT<num_sv>::FwdProjCamera cam,
    SphericalVoronoi3DGUT<num_sv>::Screen& proj, int4& aabb
) {
    float2 xy;
    Slang3DGSSV::projection_3dgut_sv_fisheye(
        false,
        world.mean, world.quat, world.scale, world.opacity, world.sv_sites, world.sv_colors,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
        cam.width, cam.height, cam.near_plane, cam.far_plane,
        &aabb, &xy, &proj.depth, &proj.scale, &proj.opacity, &proj.rgb
    );
    proj.mean = world.mean, proj.quat = world.quat;
}

template<int num_sv>
inline __device__ void SphericalVoronoi3DGUT<num_sv>::project_persp_vjp(
    SphericalVoronoi3DGUT<num_sv>::World world, SphericalVoronoi3DGUT<num_sv>::BwdProjCamera cam,
    SphericalVoronoi3DGUT<num_sv>::Screen v_proj,
    SphericalVoronoi3DGUT<num_sv>::World& v_world, float3x3 &v_R, float3 &v_t
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
inline __device__ void SphericalVoronoi3DGUT<num_sv>::project_persp_vjp(
    SphericalVoronoi3DGUT<num_sv>::World world, SphericalVoronoi3DGUT<num_sv>::BwdProjCamera cam,
    SphericalVoronoi3DGUT<num_sv>::Screen v_proj, SphericalVoronoi3DGUT<num_sv>::Screen vr_proj, SphericalVoronoi3DGUT<num_sv>::Screen h_proj,
    SphericalVoronoi3DGUT<num_sv>::World& v_world, float3x3 &v_R, float3 &v_t,
    float3 &vr_world_pos, float3 &h_world_pos
) {}  // TODO

template<int num_sv>
inline __device__ void SphericalVoronoi3DGUT<num_sv>::project_persp_vjp(
    SphericalVoronoi3DGUT<num_sv>::World world, SphericalVoronoi3DGUT<num_sv>::BwdProjCamera cam,
    SphericalVoronoi3DGUT<num_sv>::Screen v_proj, SphericalVoronoi3DGUT<num_sv>::Screen vr_proj, SphericalVoronoi3DGUT<num_sv>::Screen h_proj,
    SphericalVoronoi3DGUT<num_sv>::World& v_world, float3x3 &v_R, float3 &v_t,
    SphericalVoronoi3DGUT<num_sv>::World& vr_world, SphericalVoronoi3DGUT<num_sv>::World& h_world
) {}  // TODO

template<int num_sv>
inline __device__ void SphericalVoronoi3DGUT<num_sv>::project_fisheye_vjp(
    SphericalVoronoi3DGUT<num_sv>::World world, SphericalVoronoi3DGUT<num_sv>::BwdProjCamera cam,
    SphericalVoronoi3DGUT<num_sv>::Screen v_proj,
    SphericalVoronoi3DGUT<num_sv>::World& v_world, float3x3 &v_R, float3 &v_t
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
inline __device__ void SphericalVoronoi3DGUT<num_sv>::project_fisheye_vjp(
    SphericalVoronoi3DGUT<num_sv>::World world, SphericalVoronoi3DGUT<num_sv>::BwdProjCamera cam,
    SphericalVoronoi3DGUT<num_sv>::Screen v_proj, SphericalVoronoi3DGUT<num_sv>::Screen vr_proj, SphericalVoronoi3DGUT<num_sv>::Screen h_proj,
    SphericalVoronoi3DGUT<num_sv>::World& v_world, float3x3 &v_R, float3 &v_t,
    float3 &vr_world_pos, float3 &h_world_pos
) {}  // TODO

template<int num_sv>
inline __device__ void SphericalVoronoi3DGUT<num_sv>::project_fisheye_vjp(
    SphericalVoronoi3DGUT<num_sv>::World world, SphericalVoronoi3DGUT<num_sv>::BwdProjCamera cam,
    SphericalVoronoi3DGUT<num_sv>::Screen v_proj, SphericalVoronoi3DGUT<num_sv>::Screen vr_proj, SphericalVoronoi3DGUT<num_sv>::Screen h_proj,
    SphericalVoronoi3DGUT<num_sv>::World& v_world, float3x3 &v_R, float3 &v_t,
    SphericalVoronoi3DGUT<num_sv>::World& vr_world, SphericalVoronoi3DGUT<num_sv>::World& h_world
) {}  // TODO

#endif  // #ifdef __CUDACC__

typedef SphericalVoronoi3DGUT<2> SphericalVoronoi3DGUT_Default;
