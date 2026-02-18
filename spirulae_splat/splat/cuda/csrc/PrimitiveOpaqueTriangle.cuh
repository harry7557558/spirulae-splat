#pragma once

#ifdef __CUDACC__
#include "generated/slang.cuh"
namespace SlangOpaqueTriangle {
#include "generated/set_namespace.cuh"
#include "generated/primitive_opaque_triangle_eval3d.cuh"
}
#endif

#include "types.cuh"
#include "common.cuh"

#include <tuple>



struct OpaqueTriangle {
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

struct OpaqueTriangle::World {

    float3 mean;
    float4 quat;
    float3 scale;
    // float3 vert0, vert1, vert2;
    float2 hardness;
    FixedArray<float3, 16> sh_coeffs;
    FixedArray<float3, 2> ch_coeffs;

    typedef std::tuple<at::Tensor, at::Tensor, at::Tensor, at::Tensor, at::Tensor, at::Tensor, at::Tensor> TensorTuple;
    // typedef std::tuple<at::Tensor, at::Tensor, at::Tensor, at::Tensor, at::Tensor> TensorTuple;

    struct Buffer;

    struct Tensor {
        at::Tensor means;
        at::Tensor quats;
        at::Tensor scales;
        // at::Tensor verts;
        at::Tensor hardness;
        at::Tensor features_dc;
        at::Tensor features_sh;
        at::Tensor features_ch;

        Tensor() {}

        Tensor(const TensorTuple& splats) {
            means = std::get<0>(splats);
            quats = std::get<1>(splats);
            scales = std::get<2>(splats);
            hardness = std::get<3>(splats);
            // verts = std::get<0>(splats);
            // hardness = std::get<1>(splats);
            features_dc = std::get<4>(splats);
            features_sh = std::get<5>(splats);
            features_ch = std::get<6>(splats);
        }

        TensorTuple tuple() const {
            return std::make_tuple(means, quats, scales, hardness, features_dc, features_sh, features_ch);
            // return std::make_tuple(verts, hardness, features_dc, features_sh, features_ch);
        }

        Tensor allocProjBwd(bool is_hess_diag) const {
            return Tensor(std::make_tuple(
                at::zeros_like(means),
                at::zeros_like(quats),
                at::zeros_like(scales),
                // at::zeros_like(verts),
                at::zeros_like(hardness),
                at::zeros_like(features_dc),
                at::zeros_like(features_sh),  // TODO: exclude when is_hess_diag
                at::zeros_like(features_ch)
            ));
        }

        auto options() const {
            return means.options();
            // return verts.options();
        }
        long size() const {
            return quats.size(-2);
            // return verts.size(-3);
        }
        long batchSize() const {
            return quats.numel() / (4*size());
            // return verts.numel() / (3*3*size());
        }

        Buffer buffer() { return Buffer(*this); }
    };

    struct Buffer {
        float3* __restrict__ means;
        float4* __restrict__ quats;
        float3* __restrict__ scales;
        // float3* __restrict__ verts;
        float2* __restrict__ hardness;
        float3* __restrict__ features_dc;
        float3* __restrict__ features_sh;
        float3* __restrict__ features_ch;
        uint num_sh;

        Buffer() {}

        Buffer(const Tensor& tensors) {
            DEVICE_GUARD(tensors.means);
            CHECK_INPUT(tensors.means);
            CHECK_INPUT(tensors.quats);
            CHECK_INPUT(tensors.scales);
            // DEVICE_GUARD(tensors.verts);
            // CHECK_INPUT(tensors.verts);
            CHECK_INPUT(tensors.hardness);
            CHECK_INPUT(tensors.features_dc);
            CHECK_INPUT(tensors.features_sh);
            CHECK_INPUT(tensors.features_ch);
            means = (float3*)tensors.means.data_ptr<float>();
            quats = (float4*)tensors.quats.data_ptr<float>();
            scales = (float3*)tensors.scales.data_ptr<float>();
            // verts = (float3*)tensors.verts.data_ptr<float>();
            hardness = (float2*)tensors.hardness.data_ptr<float>();
            features_dc = (float3*)tensors.features_dc.data_ptr<float>();
            features_sh = (float3*)tensors.features_sh.data_ptr<float>();
            num_sh = tensors.features_sh.size(-2);
            features_ch = (float3*)tensors.features_ch.data_ptr<float>();
        }
    };

#ifdef __CUDACC__

    static __device__ World load(const Buffer &buffer, long idx) {
        World world = {
            buffer.means[idx],
            buffer.quats[idx],
            buffer.scales[idx],
            // buffer.verts[3*idx+0],
            // buffer.verts[3*idx+1],
            // buffer.verts[3*idx+2],
            buffer.hardness[idx]
        };
        world.sh_coeffs[0] = buffer.features_dc[idx];
        for (int i = 0; i < 15; i++)
            world.sh_coeffs[i+1] = i < buffer.num_sh ?
                buffer.features_sh[idx*buffer.num_sh+i] : make_float3(0);
        world.ch_coeffs[0] = buffer.features_ch[2*idx+0];
        world.ch_coeffs[1] = buffer.features_ch[2*idx+1];
        return world;
    }

    static __device__ __forceinline__ World zero() {
        World world = {
            {0.f, 0.f, 0.f},
            {0.f, 0.f, 0.f, 0.f},
            // {0.f, 0.f, 0.f},
            {0.f, 0.f, 0.f},
            0.f
        };
        for (int i = 0; i < 16; i++)
            world.sh_coeffs[i] = make_float3(0);
        world.ch_coeffs[0] = world.ch_coeffs[1] = make_float3(0);
        return world;
    }

    template<typename Partition>
    __device__ void reduce(Partition& partition) {
        warpSum(mean, partition);
        warpSum(quat, partition);
        warpSum(scale, partition);
        // warpSum(vert0, partition);
        // warpSum(vert1, partition);
        // warpSum(vert2, partition);
        warpSum(hardness, partition);
        for (int i = 0; i < 16; i++)
            warpSum(sh_coeffs[i], partition);
        warpSum(ch_coeffs[0], partition);
        warpSum(ch_coeffs[1], partition);
    }
    
    __device__ void saveParamsToBuffer(Buffer &buffer, long idx) {
        buffer.means[idx] = mean;
        buffer.quats[idx] = quat;
        buffer.scales[idx] = scale;
        // buffer.verts[3*idx+0] = vert0;
        // buffer.verts[3*idx+1] = vert1;
        // buffer.verts[3*idx+2] = vert2;
        buffer.hardness[idx] = hardness;
        buffer.features_dc[idx] = sh_coeffs[0];
        for (int i = 0; i < buffer.num_sh; i++)
            buffer.features_sh[idx*buffer.num_sh + i] = sh_coeffs[i+1];
        buffer.features_ch[2*idx+0] = ch_coeffs[0];
        buffer.features_ch[2*idx+1] = ch_coeffs[1];
    }

    __device__ void atomicAddGradientToBuffer(Buffer &buffer, long idx) {
        atomicAddFVec(buffer.means + idx, mean);
        atomicAddFVec(buffer.quats + idx, quat);
        atomicAddFVec(buffer.scales + idx, scale);
        // atomicAddFVec(buffer.verts + 3*idx+0, vert0);
        // atomicAddFVec(buffer.verts + 3*idx+1, vert1);
        // atomicAddFVec(buffer.verts + 3*idx+2, vert2);
        atomicAddFVec(buffer.hardness + idx, hardness);
        atomicAddFVec(buffer.features_dc + idx, sh_coeffs[0]);
        for (int i = 0; i < buffer.num_sh; i++)
            atomicAddFVec(buffer.features_sh + idx*buffer.num_sh + i, sh_coeffs[i+1]);
        atomicAddFVec(buffer.features_ch + 2*idx+0, ch_coeffs[0]);
        atomicAddFVec(buffer.features_ch + 2*idx+1, ch_coeffs[1]);
    }

#endif  // #ifdef __CUDACC__
};


struct OpaqueTriangle::RenderOutput {

    float3 rgb;
    float depth;
    float3 normal;

    typedef std::tuple<at::Tensor, at::Tensor, at::Tensor> TensorTuple;

    struct Buffer;

    struct Tensor {
        at::Tensor rgbs;
        at::Tensor depths;
        at::Tensor normals;

        Tensor() {}

        Tensor(const TensorTuple& images) {
            rgbs = std::get<0>(images);
            depths = std::get<1>(images);
            normals = std::get<2>(images);
        }

        TensorTuple tuple() const {
            return std::make_tuple(rgbs, depths, normals);
        }

        static Tensor empty(at::DimVector dims, at::TensorOptions opt) {
            at::DimVector rgbs_dims(dims); rgbs_dims.append({3});
            at::DimVector depths_dims(dims); depths_dims.append({1});
            return Tensor(std::make_tuple(
                at::empty(rgbs_dims, opt),
                at::empty(depths_dims, opt),
                at::empty(rgbs_dims, opt)
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
        float3* __restrict__ normals;

        Buffer() : rgbs(nullptr), depths(nullptr), normals(nullptr) {}

        Buffer(const Tensor& tensors) {
            DEVICE_GUARD(tensors.rgbs);
            CHECK_INPUT(tensors.rgbs);
            CHECK_INPUT(tensors.depths);
            CHECK_INPUT(tensors.normals);
            rgbs = (float3*)tensors.rgbs.data_ptr<float>();
            depths = tensors.depths.data_ptr<float>();
            normals = (float3*)tensors.normals.data_ptr<float>();
        }
    };

#ifdef __CUDACC__


    static __device__ RenderOutput load(const Buffer &buffer, long idx) {
        return {
            buffer.rgbs[idx],
            buffer.depths[idx],
            buffer.normals[idx]
        };
    }

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
        return {rgb * k, depth * k, normal * k};
    }

    __device__ __forceinline__ RenderOutput operator+(const RenderOutput &other) const {
        return {rgb + other.rgb, depth + other.depth, normal + other.normal};
    }

    __device__ __forceinline__ RenderOutput operator*(const RenderOutput &other) const {
        return {rgb * other.rgb, depth * other.depth, normal * other.normal};
    }

    __device__ __forceinline__ float dot(const RenderOutput &other) const {
        return (rgb.x * other.rgb.x + rgb.y * other.rgb.y + rgb.z * other.rgb.z)
            + depth * other.depth
            + (normal.x * other.normal.x + normal.y * other.normal.y + normal.z * other.normal.z);
    }

    __device__ void saveParamsToBuffer(Buffer &buffer, long idx) {
        buffer.rgbs[idx] = rgb;
        buffer.depths[idx] = depth;
        buffer.normals[idx] = normal;
    }

#endif  // #ifdef __CUDACC__
};


struct OpaqueTriangle::Screen {

    // from world
    float2 hardness;
    // from screen
    float depth;  // for sorting only
#ifdef __CUDACC__
    FixedArray<float3, 3> verts;
    FixedArray<float3, 3> rgbs;
#endif
    float3 normal;

    typedef std::tuple<at::Tensor, at::Tensor, at::Tensor, at::Tensor> TensorTupleProj;
    typedef std::tuple<at::Tensor, at::Tensor, at::Tensor, at::Tensor, at::Tensor> TensorTuple;

    struct Buffer;

    struct Tensor {
        bool hasWorld;
        at::Tensor hardness;
        at::Tensor depths;
        at::Tensor verts;
        at::Tensor rgbs;
        at::Tensor normals;

        Tensor() {}

        Tensor(const TensorTuple& splats) : hasWorld(true) {
            hardness = std::get<0>(splats);
            depths = std::get<1>(splats);
            verts = std::get<2>(splats);
            rgbs = std::get<3>(splats);
            normals = std::get<4>(splats);
        }

        Tensor(const TensorTupleProj& splats) : hasWorld(false) {
            depths = std::get<0>(splats);
            verts = std::get<1>(splats);
            rgbs = std::get<2>(splats);
            normals = std::get<3>(splats);
        }

        TensorTupleProj tupleProjFwd() const {
            return std::make_tuple(depths, verts, rgbs, normals);
        }

        TensorTuple tupleProjFwdPacked() const {
            return std::make_tuple(hardness, depths, verts, rgbs, normals);
        }

        TensorTuple tupleRasterBwd() const {
            return std::make_tuple(hardness, depths, verts, rgbs, normals);
        }

        static TensorTupleProj allocProjFwd(long C, long N, c10::TensorOptions opt) {
            return std::make_tuple(
                at::empty({C, N}, opt),
                at::empty({C, N, 3, 3}, opt),
                at::empty({C, N, 3, 3}, opt),
                at::empty({C, N, 3}, opt)
            );
        }

        static TensorTuple allocProjFwdPacked(long N, c10::TensorOptions opt) {
            return std::make_tuple(
                at::empty({N, 2}, opt),
                at::empty({N}, opt),
                at::empty({N, 3, 3}, opt),
                at::empty({N, 3, 3}, opt),
                at::empty({N, 3}, opt)
            );
        }

        Tensor allocRasterBwd() const {
            if (!hasWorld)
                throw std::runtime_error("!hasWorld");
            Tensor result = Tensor(std::make_tuple(
                at::zeros_like(hardness),
                at::zeros_like(depths),
                at::zeros_like(verts),
                at::zeros_like(rgbs),
                at::zeros_like(normals)
            ));
            return result;
        }

        auto options() const {
            return rgbs.options();
        }
        bool isPacked() const {
            return rgbs.dim() == 3;
        }
        long size() const {
            return rgbs.size(-3);
        }
        long batchSize() const {
            return rgbs.numel() / (9*size());
        }

        Buffer buffer() { return Buffer(*this); }
    };

    struct Buffer {
        float2* __restrict__ hardness = nullptr;
        float* __restrict__ depths;
        float3* __restrict__ verts;
        float3* __restrict__ rgbs;
        float3* __restrict__ normals;
        long size;

        Buffer() {}

        Buffer(const Tensor& tensors) {
            DEVICE_GUARD(tensors.verts);
            if (tensors.hasWorld) {
                CHECK_INPUT(tensors.hardness);
                hardness = (float2*)tensors.hardness.data_ptr<float>();
            }
            CHECK_INPUT(tensors.depths);
            CHECK_INPUT(tensors.verts);
            CHECK_INPUT(tensors.rgbs);
            depths = tensors.depths.data_ptr<float>();
            verts = (float3*)tensors.verts.data_ptr<float>();
            rgbs = (float3*)tensors.rgbs.data_ptr<float>();
            normals = (float3*)tensors.normals.data_ptr<float>();
            size = tensors.hasWorld ?
                tensors.hardness.numel() / 2 : tensors.verts.numel() / 9;
        }
    };

#ifdef __CUDACC__

    static __device__ Screen load(const Buffer &buffer, long idx) {
        return {
            buffer.hardness ? buffer.hardness[idx % buffer.size] : make_float2(0.f),
            buffer.depths[idx],
            { buffer.verts[3*idx+0], buffer.verts[3*idx+1], buffer.verts[3*idx+2] },
            { buffer.rgbs[3*idx+0], buffer.rgbs[3*idx+1], buffer.rgbs[3*idx+2] },
            buffer.normals[idx]
        };
    }

    static __device__ __forceinline__ Screen loadWithPrecompute(const Buffer &buffer, long idx) {
        return Screen::load(buffer, idx);
    }

    static __device__ __forceinline__ Screen zero() {
        return {
            {0.f, 0.f},
            0.f,
            {make_float3(0.f), make_float3(0.f), make_float3(0.f)},
            {make_float3(0.f), make_float3(0.f), make_float3(0.f)},
            {0.f, 0.f, 0.f},
        };
    }

    __device__ __forceinline__ void addGradient(const Screen &grad, float weight=1.0f) {
        hardness += grad.hardness * weight;
        depth += grad.depth * weight;
        #pragma unroll
        for (int i = 0; i < 3; i++) {
            verts[i] += grad.verts[i] * weight;
            rgbs[i] += grad.rgbs[i] * weight;
        }
        normal += grad.normal * weight;
    }

    __device__ __forceinline__ void addGaussNewtonHessianDiagonal(Screen &result, const Screen &grad, float weight=1.0f) const {
        result.hardness += grad.hardness * grad.hardness * weight;
        result.depth += grad.depth * grad.depth * weight;
        #pragma unroll
        for (int i = 0; i < 3; i++) {
            result.verts[i] += grad.verts[i] * grad.verts[i] * weight;
            result.rgbs[i] += grad.rgbs[i] * grad.rgbs[i] * weight;
        }
        result.normal += grad.normal * grad.normal * weight;
    }

    __device__ void saveParamsToBuffer(Buffer &buffer, long idx) {
        if (buffer.hardness != nullptr && idx < buffer.size)
            buffer.hardness[idx] = hardness;
        buffer.depths[idx] = depth;
        #pragma unroll
        for (int i = 0; i < 3; i++) {
            buffer.verts[3*idx+i] = verts[i];
            buffer.rgbs[3*idx+i] = rgbs[i];
        }
        buffer.normals[idx] = normal;
    }

    static __device__ void atomicAddGradientToBuffer(const Screen &grad, Buffer &buffer, long idx) {
        if (buffer.hardness != nullptr)
            atomicAddFVec(buffer.hardness + idx % buffer.size, grad.hardness);
        atomicAddFVec(buffer.depths + idx, grad.depth);
        #pragma unroll
        for (int i = 0; i < 3; i++) {
            atomicAddFVec(buffer.verts + 3*idx+i, grad.verts[i]);
            atomicAddFVec(buffer.rgbs + 3*idx+i, grad.rgbs[i]);
        }
        atomicAddFVec(buffer.normals + idx, grad.normal);
    }

    static __device__ void atomicAddAccumulatedGradientToBuffer(const Screen &grad, Buffer &buffer, long idx) {
        atomicAddGradientToBuffer(grad, buffer, idx);
    }

    static __device__ __forceinline__ void atomicAddGaussNewtonHessianDiagonalToBuffer(const Screen &grad, Buffer &buffer, long idx, float weight=1.0f) {
        if (buffer.hardness != nullptr)
            atomicAddFVec(buffer.hardness + idx % buffer.size, grad.hardness * grad.hardness * weight);
        atomicAddFVec(buffer.depths + idx, grad.depth * grad.depth * weight);
        #pragma unroll
        for (int i = 0; i < 3; i++) {
            atomicAddFVec(buffer.verts + 3*idx+i, grad.verts[i] * grad.verts[i] * weight);
            atomicAddFVec(buffer.rgbs + 3*idx+i, grad.rgbs[i] * grad.rgbs[i] * weight);
        }
        atomicAddFVec(buffer.normals + idx, grad.normal * grad.normal * weight);
    }

    __device__ __forceinline__ float evaluate_alpha(float3 ray_o, float3 ray_d) {
        return SlangOpaqueTriangle::evaluate_alpha_opaque_triangle(verts, hardness, ray_o, ray_d);
    }

    __device__ __forceinline__ Screen evaluate_alpha_vjp(
        float3 ray_o, float3 ray_d, float v_alpha,
        float3 &v_ray_o, float3 &v_ray_d
    ) {
        Screen v_splat = Screen::zero();
        SlangOpaqueTriangle::evaluate_alpha_opaque_triangle_vjp(
            verts, hardness,
            ray_o, ray_d, v_alpha,
            &v_splat.verts, &v_splat.hardness,
            &v_ray_o, &v_ray_d
        );
        return v_splat;
    }

    __device__ __forceinline__ float evaluate_sorting_depth(float3 ray_o, float3 ray_d) {
        return SlangOpaqueTriangle::evaluate_sorting_depth_opaque_triangle(
            verts, rgbs, ray_o, ray_d
        );
    }

    __device__ __forceinline__ OpaqueTriangle::RenderOutput evaluate_color(float3 ray_o, float3 ray_d) {
        float3 out_rgb; float out_depth;
        SlangOpaqueTriangle::evaluate_color_opaque_triangle(
            verts, rgbs, ray_o, ray_d,
            &out_rgb, &out_depth
        );
        return {out_rgb, out_depth, normal};
    }

    __device__ __forceinline__ Screen evaluate_color_vjp(
        float3 ray_o, float3 ray_d, OpaqueTriangle::RenderOutput v_render,
        float3 &v_ray_o, float3 &v_ray_d
    ) {
        Screen v_splat = Screen::zero();
        SlangOpaqueTriangle::evaluate_color_opaque_triangle_vjp(
            verts, rgbs, ray_o, ray_d,
            v_render.rgb, v_render.depth,
            &v_splat.verts, &v_splat.rgbs,
            &v_ray_o, &v_ray_d
        );
        v_splat.normal = v_render.normal;
        return v_splat;
    }

#endif  // #ifdef __CUDACC__

};


#ifdef __CUDACC__

inline __device__ void OpaqueTriangle::project_persp(
    OpaqueTriangle::World world, OpaqueTriangle::FwdProjCamera cam,
    OpaqueTriangle::Screen& proj, int4& aabb
) {
    SlangOpaqueTriangle::projection_opaque_triangle_eval3d_persp(
        world.mean, world.quat, world.scale, world.hardness, world.sh_coeffs, world.ch_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
        cam.width, cam.height, cam.near_plane, cam.far_plane,
        &aabb, &proj.depth, &proj.verts, &proj.rgbs, &proj.normal
    );
    proj.hardness = world.hardness;
}

inline __device__ void OpaqueTriangle::project_fisheye(
    OpaqueTriangle::World world, OpaqueTriangle::FwdProjCamera cam,
    OpaqueTriangle::Screen& proj, int4& aabb
) {
    SlangOpaqueTriangle::projection_opaque_triangle_eval3d_fisheye(
        world.mean, world.quat, world.scale, world.hardness, world.sh_coeffs, world.ch_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
        cam.width, cam.height, cam.near_plane, cam.far_plane,
        &aabb, &proj.depth, &proj.verts, &proj.rgbs, &proj.normal
    );
    proj.hardness = world.hardness;
}

inline __device__ void OpaqueTriangle::project_persp_vjp(
    OpaqueTriangle::World world, OpaqueTriangle::BwdProjCamera cam,
    OpaqueTriangle::Screen v_proj,
    OpaqueTriangle::World& v_world, float3x3 &v_R, float3 &v_t
) {
    SlangOpaqueTriangle::projection_opaque_triangle_eval3d_persp_vjp(
        world.mean, world.quat, world.scale, world.hardness, world.sh_coeffs, world.ch_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
        cam.width, cam.height,
        v_proj.depth, v_proj.verts, v_proj.rgbs, v_proj.normal,
        &v_world.mean, &v_world.quat, &v_world.scale, &v_world.hardness, &v_world.sh_coeffs, &v_world.ch_coeffs,
        &v_R, &v_t
    );
    v_world.hardness = v_proj.hardness;
}

inline __device__ void OpaqueTriangle::project_persp_vjp(
    OpaqueTriangle::World world, OpaqueTriangle::BwdProjCamera cam,
    OpaqueTriangle::Screen v_proj, OpaqueTriangle::Screen vr_proj, OpaqueTriangle::Screen h_proj,
    OpaqueTriangle::World& v_world, float3x3 &v_R, float3 &v_t,
    float3 &vr_world_pos, float3 &h_world_pos
) {}  // TODO

inline __device__ void OpaqueTriangle::project_persp_vjp(
    OpaqueTriangle::World world, OpaqueTriangle::BwdProjCamera cam,
    OpaqueTriangle::Screen v_proj, OpaqueTriangle::Screen vr_proj, OpaqueTriangle::Screen h_proj,
    OpaqueTriangle::World& v_world, float3x3 &v_R, float3 &v_t,
    OpaqueTriangle::World& vr_world, OpaqueTriangle::World& h_world
) {}  // TODO

inline __device__ void OpaqueTriangle::project_fisheye_vjp(
    OpaqueTriangle::World world, OpaqueTriangle::BwdProjCamera cam,
    OpaqueTriangle::Screen v_proj,
    OpaqueTriangle::World& v_world, float3x3 &v_R, float3 &v_t
) {
    SlangOpaqueTriangle::projection_opaque_triangle_eval3d_fisheye_vjp(
        world.mean, world.quat, world.scale, world.hardness, world.sh_coeffs, world.ch_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
        cam.width, cam.height,
        v_proj.depth, v_proj.verts, v_proj.rgbs, v_proj.normal,
        &v_world.mean, &v_world.quat, &v_world.scale, &v_world.hardness, &v_world.sh_coeffs, &v_world.ch_coeffs,
        &v_R, &v_t
    );
    v_world.hardness = v_proj.hardness;
}

inline __device__ void OpaqueTriangle::project_fisheye_vjp(
    OpaqueTriangle::World world, OpaqueTriangle::BwdProjCamera cam,
    OpaqueTriangle::Screen v_proj, OpaqueTriangle::Screen vr_proj, OpaqueTriangle::Screen h_proj,
    OpaqueTriangle::World& v_world, float3x3 &v_R, float3 &v_t,
    float3 &vr_world_pos, float3 &h_world_pos
) {}  // TODO

inline __device__ void OpaqueTriangle::project_fisheye_vjp(
    OpaqueTriangle::World world, OpaqueTriangle::BwdProjCamera cam,
    OpaqueTriangle::Screen v_proj, OpaqueTriangle::Screen vr_proj, OpaqueTriangle::Screen h_proj,
    OpaqueTriangle::World& v_world, float3x3 &v_R, float3 &v_t,
    OpaqueTriangle::World& vr_world, OpaqueTriangle::World& h_world
) {}  // TODO

#endif  // #ifdef __CUDACC__
