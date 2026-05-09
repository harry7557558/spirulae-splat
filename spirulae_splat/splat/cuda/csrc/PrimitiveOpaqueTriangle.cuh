#pragma once

#ifdef __CUDACC__
#include "generated/slang.cuh"
namespace SlangOpaqueTriangle {
#include "generated/set_namespace.cuh"
#include "generated/primitive_opaque_triangle_eval3d.cuh"
}
#endif

#include "Primitive.cuh"

#if 0
struct OpaqueTriangle {
    struct World;
    struct Screen;
    static constexpr RenderOutputType pixelType = RenderOutputType::RGB_DN;

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

struct OpaqueTriangle::World {

    float3 mean;
    float4 quat;
    float3 scale;
    // float3 vert0, vert1, vert2;
    float2 hardness;
    FixedArray<float3, 16> sh_coeffs;
    FixedArray<float3, 2> ch_coeffs;

    #ifndef NO_TORCH
    typedef std::tuple<at::Tensor, at::Tensor, at::Tensor, at::Tensor, at::Tensor, at::Tensor, at::Tensor> TensorTuple;
    // typedef std::tuple<at::Tensor, at::Tensor, at::Tensor, at::Tensor, at::Tensor> TensorTuple;
    #endif

    struct Buffer;

    struct Buffer : public TensorArray<7> {
        using TensorArray<7>::TensorArray;
        uint num_sh;

        Buffer() : num_sh(0) {}

        #ifndef NO_TORCH
        Buffer(const Tensor& tensors)
            : TensorArray<7>(std::vector<std::optional<at::Tensor>>{
                tensors.means, tensors.quats, tensors.scales,
                tensors.hardness, tensors.features_dc, tensors.features_sh, tensors.features_ch
            }) {
            num_sh = tensors.features_sh.has_value() ?
                tensors.features_sh.value().size(-2) : 0;
        }
        #endif

        #ifdef __CUDACC__
        __forceinline__ __device__ float2& hardness(int64_t i)
            { return *reinterpret_cast<float2*>(&_data[3][2*i]); }
        __forceinline__ __device__ float3& features_ch(int64_t i)
            { return *reinterpret_cast<float3*>(&_data[6][3*i]); }
        __forceinline__ __device__ float3& features_sh(int64_t i, int64_t j)
            { return *reinterpret_cast<float3*>(&_data[5][_strides[5]*i+j]); }
        #endif
    };

#ifdef __CUDACC__

    static __device__ World load(const Buffer &buffer, long idx) {
        World world = {
            buffer.means(idx),
            buffer.quats(idx),
            buffer.scales(idx),
            buffer.hardness(idx)
        };
        world.sh_coeffs[0] = buffer.features_dc(idx);
        for (int i = 0; i < 15; i++)
            world.sh_coeffs[i+1] = i < buffer.num_sh ?
                buffer.features_sh(idx, i) : make_float3(0);
        world.ch_coeffs[0] = buffer.features_ch(2*idx+0);
        world.ch_coeffs[1] = buffer.features_ch(2*idx+1);
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

    __device__ void saveParamsToBuffer(Buffer &buffer, long idx) {
        buffer.means(idx) = mean;
        buffer.quats(idx) = quat;
        buffer.scales(idx) = scale;
        buffer.hardness(idx) = hardness;
        buffer.features_dc(idx) = sh_coeffs[0];
        for (int i = 0; i < buffer.num_sh; i++)
            buffer.features_sh(idx, i) = sh_coeffs[i+1];
        buffer.features_ch(2*idx+0) = ch_coeffs[0];
        buffer.features_ch(2*idx+1) = ch_coeffs[1];
    }

    __device__ void atomicAddGradientToBuffer(Buffer &buffer, long idx) {
        atomicAddFVec(buffer.means(idx), mean);
        atomicAddFVec(buffer.quats(idx), quat);
        atomicAddFVec(buffer.scales(idx), scale);
        atomicAddFVec(buffer.hardness(idx), hardness);
        atomicAddFVec(buffer.features_dc(idx), sh_coeffs[0]);
        for (int i = 0; i < buffer.num_sh; i++)
            atomicAddFVec(buffer.features_sh(idx, i), sh_coeffs[i+1]);
        atomicAddFVec(buffer.features_ch(2*idx+0), ch_coeffs[0]);
        atomicAddFVec(buffer.features_ch(2*idx+1), ch_coeffs[1]);
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

    #ifndef NO_TORCH
    typedef std::tuple<at::Tensor, at::Tensor, at::Tensor, at::Tensor> TensorTupleProj;
    typedef std::tuple<at::Tensor, at::Tensor, at::Tensor, at::Tensor, at::Tensor> TensorTuple;
    #endif

    struct Buffer;

    #ifndef NO_TORCH
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
                zeros_like_tensor(hardness),
                zeros_like_tensor(depths),
                zeros_like_tensor(verts),
                zeros_like_tensor(rgbs),
                zeros_like_tensor(normals)
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
    #endif

    struct Buffer : public TensorArray<5> {
        long size;

        using TensorArray<5>::TensorArray;

        Buffer() : size(0) {}

        #ifndef NO_TORCH
        Buffer(const Tensor& tensors)
            : TensorArray<5>(std::vector<std::optional<at::Tensor>>{
                tensors.hardness, tensors.depths, tensors.verts, tensors.rgbs, tensors.normals
            }) {
            size = tensors.hasWorld ?
                tensors.hardness.value().numel() / 2 : tensors.verts.numel() / 9;
        }
        #endif

        #ifdef __CUDACC__
        __forceinline__ __device__ bool hasHardness() const {
            return _data[0] != nullptr;
        }
        __forceinline__ __device__ float2& hardness(int64_t i)
            { return *reinterpret_cast<float2*>(&_data[0][2*i]); }
        __forceinline__ __device__ float& depths(int64_t i)
            { return *reinterpret_cast<float*>(&_data[1][i]); }
        __forceinline__ __device__ float3& verts(int64_t i)
            { return *reinterpret_cast<float3*>(&_data[2][3*i]); }
        __forceinline__ __device__ float3& rgbs(int64_t i)
            { return *reinterpret_cast<float3*>(&_data[3][3*i]); }
        __forceinline__ __device__ float3& normals(int64_t i)
            { return *reinterpret_cast<float3*>(&_data[4][3*i]); }
        #endif
    };

#ifdef __CUDACC__

    static __device__ Screen load(const Buffer &buffer, long idx, const uint32_t* gaussian_ids) {
        uint32_t idx0 = gaussian_ids ? gaussian_ids[idx] : idx % buffer.size;
        return {
            buffer.hasHardness() ? buffer.hardness(idx0) : make_float2(0.f),
            buffer.depths(idx),
            { buffer.verts(3*idx+0), buffer.verts(3*idx+1), buffer.verts(3*idx+2) },
            { buffer.rgbs(3*idx+0), buffer.rgbs(3*idx+1), buffer.rgbs(3*idx+2) },
            buffer.normals(idx)
        };
    }

    static __device__ __forceinline__ Screen loadWithPrecompute(const Buffer &buffer, long idx, const uint32_t* gaussian_ids) {
        return Screen::load(buffer, idx, gaussian_ids);
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

    __device__ __forceinline__ void precomputeBackward(Screen& grad) const {
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

    __device__ __forceinline__ void addGaussNewtonHessianDiagonal(const Screen &grad, float weight=1.0f) {
        hardness += grad.hardness * grad.hardness * weight;
        depth += grad.depth * grad.depth * weight;
        #pragma unroll
        for (int i = 0; i < 3; i++) {
            verts[i] += grad.verts[i] * grad.verts[i] * weight;
            rgbs[i] += grad.rgbs[i] * grad.rgbs[i] * weight;
        }
        normal += grad.normal * grad.normal * weight;
    }

    __device__ void saveParamsToBuffer(Buffer &buffer, long idx, const uint32_t* gaussian_ids) {
        if (buffer.hasHardness() && idx < buffer.size)
            buffer.hardness(idx) = hardness;
        buffer.depths(idx) = depth;
        #pragma unroll
        for (int i = 0; i < 3; i++) {
            buffer.verts(3*idx+i) = verts[i];
            buffer.rgbs(3*idx+i) = rgbs[i];
        }
        buffer.normals(idx) = normal;
    }

    __device__ void atomicAddToBuffer(Buffer &buffer, long idx, const uint32_t* gaussian_ids) const {
        uint32_t idx0 = gaussian_ids ? gaussian_ids[idx] : idx % buffer.size;
        if (buffer.hasHardness())
            atomicAddFVec(buffer.hardness(idx0), hardness);
        atomicAddFVec(buffer.depths(idx), depth);
        #pragma unroll
        for (int i = 0; i < 3; i++) {
            atomicAddFVec(buffer.verts(3*idx+i), verts[i]);
            atomicAddFVec(buffer.rgbs(3*idx+i), rgbs[i]);
        }
        atomicAddFVec(buffer.normals(idx), normal);
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

    __device__ __forceinline__ RenderOutput evaluate_color(float3 ray_o, float3 ray_d) {
        float3 out_rgb; float out_depth;
        SlangOpaqueTriangle::evaluate_color_opaque_triangle(
            verts, rgbs, ray_o, ray_d,
            &out_rgb, &out_depth
        );
        return {out_rgb, out_depth, normal};
    }

    __device__ __forceinline__ Screen evaluate_color_vjp(
        float3 ray_o, float3 ray_d, RenderOutput v_render,
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
    OpaqueTriangle::World world, ProjCamera cam,
    OpaqueTriangle::Screen& proj, float4& aabb
) {
    SlangOpaqueTriangle::projection_opaque_triangle_eval3d_persp(
        world.mean, world.quat, world.scale, world.hardness, world.sh_coeffs, world.ch_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
        cam.width, cam.height,
        &aabb, &proj.depth, &proj.verts, &proj.rgbs, &proj.normal
    );
    proj.hardness = world.hardness;
}

inline __device__ void OpaqueTriangle::project_fisheye(
    OpaqueTriangle::World world, ProjCamera cam,
    OpaqueTriangle::Screen& proj, float4& aabb
) {
    SlangOpaqueTriangle::projection_opaque_triangle_eval3d_fisheye(
        world.mean, world.quat, world.scale, world.hardness, world.sh_coeffs, world.ch_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
        cam.width, cam.height,
        &aabb, &proj.depth, &proj.verts, &proj.rgbs, &proj.normal
    );
    proj.hardness = world.hardness;
}

inline __device__ void OpaqueTriangle::project_persp_vjp(
    OpaqueTriangle::World world, ProjCamera cam,
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
    OpaqueTriangle::World world, ProjCamera cam,
    OpaqueTriangle::Screen v_proj, OpaqueTriangle::Screen vr_proj, OpaqueTriangle::Screen h_proj,
    OpaqueTriangle::World& v_world, float3x3 &v_R, float3 &v_t,
    float3 &vr_world_pos, float3 &h_world_pos
) {}  // TODO

inline __device__ void OpaqueTriangle::project_persp_vjp(
    OpaqueTriangle::World world, ProjCamera cam,
    OpaqueTriangle::Screen v_proj, OpaqueTriangle::Screen vr_proj, OpaqueTriangle::Screen h_proj,
    OpaqueTriangle::World& v_world, float3x3 &v_R, float3 &v_t,
    OpaqueTriangle::World& vr_world, OpaqueTriangle::World& h_world
) {}  // TODO

inline __device__ void OpaqueTriangle::project_fisheye_vjp(
    OpaqueTriangle::World world, ProjCamera cam,
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
    OpaqueTriangle::World world, ProjCamera cam,
    OpaqueTriangle::Screen v_proj, OpaqueTriangle::Screen vr_proj, OpaqueTriangle::Screen h_proj,
    OpaqueTriangle::World& v_world, float3x3 &v_R, float3 &v_t,
    float3 &vr_world_pos, float3 &h_world_pos
) {}  // TODO

inline __device__ void OpaqueTriangle::project_fisheye_vjp(
    OpaqueTriangle::World world, ProjCamera cam,
    OpaqueTriangle::Screen v_proj, OpaqueTriangle::Screen vr_proj, OpaqueTriangle::Screen h_proj,
    OpaqueTriangle::World& v_world, float3x3 &v_R, float3 &v_t,
    OpaqueTriangle::World& vr_world, OpaqueTriangle::World& h_world
) {}  // TODO

#endif  // #ifdef __CUDACC__
#endif