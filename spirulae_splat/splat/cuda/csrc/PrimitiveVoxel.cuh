#pragma once

#ifdef __CUDACC__
#include "generated/slang.cuh"
namespace SlangVoxel {
#include "generated/set_namespace.cuh"
#include "generated/primitive_voxel.cuh"
}
#endif

#include "types.cuh"
#include "common.cuh"

#include <tuple>



struct VoxelPrimitive {
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

struct VoxelPrimitive::World {

    float3 pos;
    float size;
    FixedArray<float, 8> densities;
    FixedArray<float3, 16> sh_coeffs;

    typedef std::tuple<std::optional<at::Tensor>, std::optional<at::Tensor>, at::Tensor, at::Tensor> TensorTuple;

    struct Buffer;

    struct Tensor {
        std::optional<at::Tensor> pos_size;
        std::optional<at::Tensor> densities;
        at::Tensor features_dc;
        at::Tensor features_sh;

        Tensor() {}

        Tensor(const TensorTuple& splats) {
            pos_size = std::get<0>(splats);
            densities = std::get<1>(splats);
            features_dc = std::get<2>(splats);
            features_sh = std::get<3>(splats);
        }

        TensorTuple tuple() const {
            return std::make_tuple(pos_size, densities, features_dc, features_sh);
        }

        Tensor allocProjBwd(bool is_hess_diag) const {
            return Tensor(std::make_tuple(
                (std::optional<at::Tensor>)at::nullopt,
                (std::optional<at::Tensor>)at::nullopt,
                at::zeros_like(features_dc),
                at::zeros_like(features_sh)  // TODO: exclude SH when is_hess_diag
            ));
        }

        auto options() const {
            return features_dc.options();
        }
        long size() const {
            return features_dc.size(-2);
        }
        long batchSize() const {
            return features_dc.numel() / (3*size());
        }

        Buffer buffer() { return Buffer(*this); }
    };

    struct Buffer {
        float4* __restrict__ pos_size;
        float4* __restrict__ densities;  // 8 densities = 2x float4
        float3* __restrict__ features_dc;
        float3* __restrict__ features_sh;
        uint num_sh;

        Buffer() {}

        Buffer(const Tensor& tensors) {
            DEVICE_GUARD(tensors.features_dc);
            if (tensors.pos_size.has_value())
                CHECK_INPUT(tensors.pos_size.value());
            if (tensors.pos_size.has_value())
                CHECK_INPUT(tensors.densities.value());
            CHECK_INPUT(tensors.features_dc);
            CHECK_INPUT(tensors.features_sh);
            pos_size = tensors.pos_size.has_value() ?
                (float4*)tensors.pos_size.value().data_ptr<float>() : nullptr;
            densities = tensors.densities.has_value() ?
                (float4*)tensors.densities.value().data_ptr<float>() : nullptr;
            features_dc = (float3*)tensors.features_dc.data_ptr<float>();
            features_sh = (float3*)tensors.features_sh.data_ptr<float>();
            num_sh = tensors.features_sh.size(-2);
        }
    };

#ifdef __CUDACC__

    static __device__ World load(const Buffer &buffer, long idx) {
        float4 pos_size = buffer.pos_size[idx];
        float4 d0 = buffer.densities[2*idx],
            d1 = buffer.densities[2*idx+1];
        World world = {
            {pos_size.x, pos_size.y, pos_size.z},
            pos_size.w,
            {d0.x, d0.y, d0.z, d0.w, d1.x, d1.y, d1.z, d1.w}
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
            0.f,
            {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f},
        };
        for (int i = 0; i < 16; i++)
            world.sh_coeffs[i] = make_float3(0);
        return world;
    }

    template<typename Partition>
    __device__ void reduce(Partition& partition) {
        // no need to do this for pos/size in gradient calculation
        for (int i = 0; i < 8; i++)
            warpSum(densities[i], partition);
        for (int i = 0; i < 16; i++)
            warpSum(sh_coeffs[i], partition);
    }

    __device__ void saveParamsToBuffer(Buffer &buffer, long idx) {
        if (buffer.pos_size != nullptr)
            buffer.pos_size[idx] = {pos.x, pos.y, pos.z, size};
        buffer.densities[2*idx] = {densities[0], densities[1], densities[2], densities[3]};
        buffer.densities[2*idx+1] = {densities[4], densities[5], densities[6], densities[7]};
        buffer.features_dc[idx] = sh_coeffs[0];
        for (int i = 0; i < buffer.num_sh; i++)
            buffer.features_sh[idx*buffer.num_sh + i] = sh_coeffs[i+1];
    }

    __device__ void atomicAddGradientToBuffer(Buffer &buffer, long idx) {
        if (buffer.pos_size != nullptr)
            atomicAddFVec(buffer.pos_size + idx, {pos.x, pos.y, pos.z, size});
        atomicAddFVec(buffer.densities + 2*idx, {densities[0], densities[1], densities[2], densities[3]});
        atomicAddFVec(buffer.densities + 2*idx+1, {densities[4], densities[5], densities[6], densities[7]});
        atomicAddFVec(buffer.features_dc + idx, sh_coeffs[0]);
        for (int i = 0; i < buffer.num_sh; i++)
            atomicAddFVec(buffer.features_sh + idx*buffer.num_sh + i, sh_coeffs[i+1]);
    }

#endif  // #ifdef __CUDACC__
};


struct VoxelPrimitive::RenderOutput {

    float3 rgb;
    float depth;

    typedef std::tuple<at::Tensor, at::Tensor> TensorTuple;

    struct Buffer;

    struct Tensor {
        at::Tensor rgbs;
        at::Tensor depths;

        Tensor() {}

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
        return (rgb.x * other.rgb.x + rgb.y * other.rgb.y + rgb.z * other.rgb.z)
            + depth * other.depth;
    }

    __device__ void saveParamsToBuffer(Buffer &buffer, long idx) {
        buffer.rgbs[idx] = rgb;
        buffer.depths[idx] = depth;
    }

#endif  // #ifdef __CUDACC__
};



struct VoxelPrimitive::Screen {

    float3 pos;
    float size;
    float depth;  // non differentiable
    FixedArray<float, 8> densities;
    float3 rgb;
    float density_abs;

    typedef std::tuple<std::optional<at::Tensor>, at::Tensor> TensorTupleProj;
    typedef std::tuple<std::optional<at::Tensor>, std::optional<at::Tensor>, std::optional<at::Tensor>, at::Tensor> TensorTuple;

    struct Buffer;

    struct Tensor {
        bool hasWorld;
        std::optional<at::Tensor> pos_size;
        std::optional<at::Tensor> depths;
        std::optional<at::Tensor> densities;
        at::Tensor rgbs;
        std::optional<at::Tensor> absgrad;

        Tensor() {}

        Tensor(const TensorTuple& splats) : hasWorld(true) {
            pos_size = std::get<0>(splats);
            depths = std::get<1>(splats);
            densities = std::get<2>(splats);
            rgbs = std::get<3>(splats);
        }

        Tensor(const TensorTupleProj& splats) : hasWorld(false) {
            depths = std::get<0>(splats);
            rgbs = std::get<1>(splats);
        }

        TensorTupleProj tupleProjFwd() const {
            return std::make_tuple(depths, rgbs);
        }

        TensorTuple tupleProjFwdPacked() const {
            return std::make_tuple(pos_size, depths, densities, rgbs);
        }

        TensorTuple tupleRasterBwd() const {
            return std::make_tuple(pos_size, depths, densities, rgbs);
        }

        static TensorTupleProj allocProjFwd(long C, long N, c10::TensorOptions opt) {
            return std::make_tuple(
                (std::optional<at::Tensor>)at::empty({C, N}, opt),
                at::empty({C, N, 3}, opt)
            );
        }

        static TensorTuple allocProjFwdPacked(long N, c10::TensorOptions opt) {
            return std::make_tuple(
                (std::optional<at::Tensor>)at::empty({N, 4}, opt),
                (std::optional<at::Tensor>)at::empty({N}, opt),
                (std::optional<at::Tensor>)at::empty({N, 8}, opt),
                at::empty({N, 3}, opt)
            );
        }

        Tensor allocRasterBwd() const {
            if (!hasWorld)
                throw std::runtime_error("!hasWorld");
            Tensor result = Tensor(std::make_tuple(
                (std::optional<at::Tensor>)at::nullopt,
                (std::optional<at::Tensor>)at::nullopt,
                (std::optional<at::Tensor>)at::zeros_like(densities.value()),
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
            return rgbs.size(-2);
        }
        long batchSize() const {
            return rgbs.numel() / (3*size());
        }

        Buffer buffer() { return Buffer(*this); }
    };

    struct Buffer {
        float4* __restrict__ pos_size = nullptr;
        float* __restrict__ depths = nullptr;
        float4* __restrict__ densities = nullptr;
        float3* __restrict__ rgbs;
        long size;

        Buffer() {}

        Buffer(const Tensor& tensors) {
            DEVICE_GUARD(tensors.rgbs);
            if (tensors.pos_size.has_value()) {
                CHECK_INPUT(tensors.pos_size.value());
                pos_size = (float4*)tensors.pos_size.value().data_ptr<float>();
            }
            if (tensors.depths.has_value()) {
                CHECK_INPUT(tensors.depths.value());
                depths = tensors.depths.value().data_ptr<float>();
            }
            if (tensors.densities.has_value())
                CHECK_INPUT(tensors.densities.value());
            CHECK_INPUT(tensors.rgbs);
            densities = tensors.densities.has_value() ?
                (float4*)tensors.densities.value().data_ptr<float>() : nullptr;
            rgbs = (float3*)tensors.rgbs.data_ptr<float>();
            size = tensors.rgbs.size(-2);
        }
    };

#ifdef __CUDACC__

    static __device__ Screen load(const Buffer &buffer, long idx) {
        long idx0 = idx % buffer.size;
        float4 pos_size = buffer.pos_size ? buffer.pos_size[idx0] : make_float4(0, 0, 0, 0);
        float4 d0 = buffer.densities ? buffer.densities[2*idx0] : make_float4(0, 0, 0, 0),
            d1 = buffer.densities ? buffer.densities[2*idx0+1] : make_float4(0, 0, 0, 0);
        return {
            {pos_size.x, pos_size.y, pos_size.z},
            pos_size.w,
            buffer.depths ? buffer.depths[idx] : 0.0f,
            {d0.x, d0.y, d0.z, d0.w, d1.x, d1.y, d1.z, d1.w},
            buffer.rgbs[idx]
        };
    }

    static __device__ __forceinline__ Screen loadWithPrecompute(const Buffer &buffer, long idx) {
        return Screen::load(buffer, idx);
    }

    static __device__ __forceinline__ Screen zero() {
        return {
            {0.f, 0.f, 0.f},
            0.f,
            0.f,
            {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f},
            {0.f, 0.f, 0.f},
        };
    }

    __device__ __forceinline__ void addGradient(const Screen &grad, float weight=1.0f) {
        pos += grad.pos * weight;
        size += grad.size * weight;
        depth += grad.depth * weight;
        #pragma unroll
        for (int i = 0; i < 8; i++)
            densities[i] += grad.densities[i] * weight;
        rgb += grad.rgb * weight;
    }

    __device__ __forceinline__ void addGaussNewtonHessianDiagonal(Screen &result, const Screen &grad, float weight=1.0f) const {
        result.pos += grad.pos * grad.pos * weight;
        result.size += grad.size * grad.size * weight;
        result.depth += grad.depth * grad.depth * weight;
        #pragma unroll
        for (int i = 0; i < 8; i++)
            result.densities[i] += grad.densities[i] * grad.densities[i] * weight;
        result.rgb += grad.rgb * grad.rgb * weight;
    }

    __device__ void saveParamsToBuffer(Buffer &buffer, long idx) {
        if (buffer.pos_size != nullptr && idx < buffer.size)
            buffer.pos_size[idx] = {pos.x, pos.y, pos.z, size};
        if (buffer.depths != nullptr)
            buffer.depths[idx] = depth;
        if (buffer.densities != nullptr && idx < buffer.size) {
            buffer.densities[2*idx] = {densities[0], densities[1], densities[2], densities[3]};
            buffer.densities[2*idx+1] = {densities[4], densities[5], densities[6], densities[7]};
        }
        buffer.rgbs[idx] = rgb;
    }

    static __device__ void atomicAddGradientToBuffer(const Screen &grad, Buffer &buffer, long idx) {
        long idx0 = idx % buffer.size;
        if (buffer.pos_size != nullptr)
            atomicAddFVec(buffer.pos_size + idx0, {grad.pos.x, grad.pos.y, grad.pos.z, grad.size});
        if (buffer.depths != nullptr)
            atomicAddFVec(buffer.depths + idx, grad.depth);
        // if (buffer.densities != nullptr) {
        if (true) {  // should not be nullptr, spot bug if crash
            atomicAddFVec(buffer.densities + 2*idx0, {grad.densities[0], grad.densities[1], grad.densities[2], grad.densities[3]});
            atomicAddFVec(buffer.densities + 2*idx0+1, {grad.densities[4], grad.densities[5], grad.densities[6], grad.densities[7]});
        }
        atomicAddFVec(buffer.rgbs + idx, grad.rgb);
    }

    static __device__ void atomicAddAccumulatedGradientToBuffer(const Screen &grad, Buffer &buffer, long idx) {
        atomicAddGradientToBuffer(grad, buffer, idx);
    }

    static __device__ __forceinline__ void atomicAddGaussNewtonHessianDiagonalToBuffer(const Screen &grad, Buffer &buffer, long idx, float weight=1.0f) {
        long idx0 = idx % buffer.size;
        if (buffer.pos_size != nullptr) {
            float4 temp = {grad.pos.x, grad.pos.y, grad.pos.z, grad.size};
            atomicAddFVec(buffer.pos_size + idx0, temp * temp * weight);
        }
        if (buffer.depths != nullptr)
            atomicAddFVec(buffer.depths + idx, grad.depth * grad.depth * weight);
        // if (buffer.densities != nullptr) {
        if (true) {  // should not be nullptr, spot bug if crash
            float4 temp = {grad.densities[0], grad.densities[1], grad.densities[2], grad.densities[3]};
            atomicAddFVec(buffer.densities + 2*idx0, temp * temp * weight);
            temp = {grad.densities[4], grad.densities[5], grad.densities[6], grad.densities[7]};
            atomicAddFVec(buffer.densities + 2*idx0+1, temp * temp * weight);
        }
        atomicAddFVec(buffer.rgbs + idx, grad.rgb * grad.rgb * weight);
    }

    __device__ __forceinline__ float evaluate_alpha(float3 ray_o, float3 ray_d) {
        return SlangVoxel::evaluate_alpha_voxel(pos, size, densities, ray_o, ray_d);
    }

    __device__ __forceinline__ Screen evaluate_alpha_vjp(
        float3 ray_o, float3 ray_d, float v_alpha,
        float3 &v_ray_o, float3 &v_ray_d
    ) {
        Screen v_splat = Screen::zero();
        SlangVoxel::evaluate_alpha_voxel_vjp(
            pos, size, densities,
            ray_o, ray_d, v_alpha,
            &v_splat.densities,
            &v_ray_o, &v_ray_d
        );
        return v_splat;
    }

    __device__ __forceinline__ VoxelPrimitive::RenderOutput evaluate_color(float3 ray_o, float3 ray_d) {
        float3 out_rgb; float out_depth;
        SlangVoxel::evaluate_color_voxel(
            pos, size, densities, rgb, ray_o, ray_d,
            &out_rgb, &out_depth
        );
        return {out_rgb, out_depth};
    }

    __device__ __forceinline__ Screen evaluate_color_vjp(
        float3 ray_o, float3 ray_d, VoxelPrimitive::RenderOutput v_render,
        float3 &v_ray_o, float3 &v_ray_d
    ) {
        Screen v_splat = Screen::zero();
        SlangVoxel::evaluate_color_voxel_vjp(
            pos, size, densities, rgb, ray_o, ray_d,
            v_render.rgb, v_render.depth,
            &v_splat.densities, &v_splat.rgb,
            &v_ray_o, &v_ray_d
        );
        return v_splat;
    }

#endif  // #ifdef __CUDACC__

};


#ifdef __CUDACC__

inline __device__ void VoxelPrimitive::project_persp(
    VoxelPrimitive::World world, VoxelPrimitive::FwdProjCamera cam,
    VoxelPrimitive::Screen& proj, int4& aabb
) {
    SlangVoxel::projection_voxel_eval3d_persp(
        world.pos, world.size, world.densities, world.sh_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
        cam.width, cam.height, cam.near_plane, cam.far_plane,
        &aabb, &proj.depth, &proj.rgb
    );
    proj.pos = world.pos, proj.size = world.size, proj.densities = world.densities;
}

inline __device__ void VoxelPrimitive::project_fisheye(
    VoxelPrimitive::World world, VoxelPrimitive::FwdProjCamera cam,
    VoxelPrimitive::Screen& proj, int4& aabb
) {
    SlangVoxel::projection_voxel_eval3d_fisheye(
        world.pos, world.size, world.densities, world.sh_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
        cam.width, cam.height, cam.near_plane, cam.far_plane,
        &aabb, &proj.depth, &proj.rgb
    );
    proj.pos = world.pos, proj.size = world.size, proj.densities = world.densities;
}

inline __device__ void VoxelPrimitive::project_persp_vjp(
    VoxelPrimitive::World world, VoxelPrimitive::BwdProjCamera cam,
    VoxelPrimitive::Screen v_proj,
    VoxelPrimitive::World& v_world, float3x3 &v_R, float3 &v_t
) {
    SlangVoxel::projection_voxel_eval3d_persp_vjp(
        world.pos, world.size, world.densities, world.sh_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
        cam.width, cam.height,
        v_proj.rgb,
        &v_world.densities, &v_world.sh_coeffs,
        &v_R, &v_t
    );
    v_world.pos = v_proj.pos, v_world.size = v_proj.size, v_world.densities = v_proj.densities;
}

inline __device__ void VoxelPrimitive::project_persp_vjp(
    VoxelPrimitive::World world, VoxelPrimitive::BwdProjCamera cam,
    VoxelPrimitive::Screen v_proj, VoxelPrimitive::Screen vr_proj, VoxelPrimitive::Screen h_proj,
    VoxelPrimitive::World& v_world, float3x3 &v_R, float3 &v_t,
    float3 &vr_world_pos, float3 &h_world_pos
) {}  // TODO

inline __device__ void VoxelPrimitive::project_persp_vjp(
    VoxelPrimitive::World world, VoxelPrimitive::BwdProjCamera cam,
    VoxelPrimitive::Screen v_proj, VoxelPrimitive::Screen vr_proj, VoxelPrimitive::Screen h_proj,
    VoxelPrimitive::World& v_world, float3x3 &v_R, float3 &v_t,
    VoxelPrimitive::World& vr_world, VoxelPrimitive::World& h_world
) {}  // TODO

inline __device__ void VoxelPrimitive::project_fisheye_vjp(
    VoxelPrimitive::World world, VoxelPrimitive::BwdProjCamera cam,
    VoxelPrimitive::Screen v_proj,
    VoxelPrimitive::World& v_world, float3x3 &v_R, float3 &v_t
) {
    SlangVoxel::projection_voxel_eval3d_fisheye_vjp(
        world.pos, world.size, world.densities, world.sh_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
        cam.width, cam.height,
        v_proj.rgb,
        &v_world.densities, &v_world.sh_coeffs,
        &v_R, &v_t
    );
    v_world.pos = v_proj.pos, v_world.size = v_proj.size, v_world.densities = v_proj.densities;
}

inline __device__ void VoxelPrimitive::project_fisheye_vjp(
    VoxelPrimitive::World world, VoxelPrimitive::BwdProjCamera cam,
    VoxelPrimitive::Screen v_proj, VoxelPrimitive::Screen vr_proj, VoxelPrimitive::Screen h_proj,
    VoxelPrimitive::World& v_world, float3x3 &v_R, float3 &v_t,
    float3 &vr_world_pos, float3 &h_world_pos
) {}  // TODO

inline __device__ void VoxelPrimitive::project_fisheye_vjp(
    VoxelPrimitive::World world, VoxelPrimitive::BwdProjCamera cam,
    VoxelPrimitive::Screen v_proj, VoxelPrimitive::Screen vr_proj, VoxelPrimitive::Screen h_proj,
    VoxelPrimitive::World& v_world, float3x3 &v_R, float3 &v_t,
    VoxelPrimitive::World& vr_world, VoxelPrimitive::World& h_world
) {}  // TODO

#endif  // #ifdef __CUDACC__
