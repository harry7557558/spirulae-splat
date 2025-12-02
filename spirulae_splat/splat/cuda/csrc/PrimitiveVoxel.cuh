#pragma once

#ifdef __CUDACC__
#define TensorView _Slang_TensorView
#include "generated/primitive.cu"
#undef TensorView
#endif

#include "types.cuh"
#include "common.cuh"

#include <tuple>



struct VoxelPrimitive {
    struct World;
    struct WorldEval3D;
    struct RenderOutput;

#ifdef __CUDACC__

    struct FwdProjCamera {
        float3x3 R;
        float3 t;
        float fx, fy, cx, cy;
        uint width, height, antialiased;
        float near_plane, far_plane;
        CameraDistortionCoeffs dist_coeffs;
    };

    inline static __device__ void project_persp_eval3d(
        World world, FwdProjCamera cam,
        WorldEval3D& proj, int4& aabb
    );

    inline static __device__ void project_fisheye_eval3d(
        World world, FwdProjCamera cam,
        WorldEval3D& proj, int4& aabb
    );

    struct BwdProjCamera {
        float3x3 R;
        float3 t;
        float fx, fy, cx, cy;
        uint width, height, antialiased;
        CameraDistortionCoeffs dist_coeffs;
    };

    inline static __device__ void project_persp_eval3d_vjp(
        World world, BwdProjCamera cam,
        WorldEval3D v_proj,
        World& v_world, float3x3 &v_R, float3 &v_t
    );

    inline static __device__ void project_fisheye_eval3d_vjp(
        World world, BwdProjCamera cam,
        WorldEval3D v_proj,
        World& v_world, float3x3 &v_R, float3 &v_t
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

        Tensor(const TensorTuple& splats) {
            pos_size = std::get<0>(splats);
            densities = std::get<1>(splats);
            features_dc = std::get<2>(splats);
            features_sh = std::get<3>(splats);
        }

        TensorTuple tuple() const {
            return std::make_tuple(pos_size, densities, features_dc, features_sh);
        }

        Tensor zeros_like() const {
            return Tensor(std::make_tuple(
                (std::optional<at::Tensor>)at::nullopt,
                (std::optional<at::Tensor>)at::nullopt,
                at::zeros_like(features_dc),
                at::zeros_like(features_sh)
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



struct VoxelPrimitive::WorldEval3D {

    float3 pos;
    float size;
    float depth;  // non differentiable
    FixedArray<float, 8> densities;
    float3 rgb;

    typedef std::tuple<std::optional<at::Tensor>, at::Tensor> TensorTupleProj;
    typedef std::tuple<std::optional<at::Tensor>, std::optional<at::Tensor>, std::optional<at::Tensor>, at::Tensor> TensorTuple;

    struct Buffer;

    struct Tensor {
        bool hasWorld;
        std::optional<at::Tensor> pos_size;
        std::optional<at::Tensor> depths;
        std::optional<at::Tensor> densities;
        at::Tensor rgbs;

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

        TensorTuple tupleAll() const {
            return std::make_tuple(pos_size, depths, densities, rgbs);
        }

        TensorTupleProj tupleProj() const {
            return std::make_tuple(depths, rgbs);
        }

        Tensor zeros_like() const {
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

        static Tensor empty(long C, long N, c10::TensorOptions opt) {
            return std::make_tuple(
                (std::optional<at::Tensor>)at::empty({N, 4}, opt),
                (std::optional<at::Tensor>)(C == -1 ? at::empty({N}, opt) : at::empty({C, N}, opt)),
                (std::optional<at::Tensor>)at::empty({N, 8}, opt),
                C == -1 ? at::empty({N, 3}, opt) : at::empty({C, N, 3}, opt)
            );
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

        Buffer(const Tensor& tensors) {
            DEVICE_GUARD(tensors.rgbs);
            if (tensors.pos_size.has_value()) {
                CHECK_INPUT(tensors.pos_size.value());
                pos_size = (float4*)tensors.pos_size.value().data_ptr<float>();
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

    static __device__ WorldEval3D load(const Buffer &buffer, long idx) {
        long idx0 = idx % buffer.size;
        float4 pos_size = buffer.pos_size ? buffer.pos_size[idx0] : make_float4(0, 0, 0, 0);
        float4 d0 = buffer.densities[2*idx0],
            d1 = buffer.densities[2*idx0+1];
        return {
            {pos_size.x, pos_size.y, pos_size.z},
            pos_size.w,
            buffer.depths ? buffer.depths[idx] : 0.0f,
            {d0.x, d0.y, d0.z, d0.w, d1.x, d1.y, d1.z, d1.w},
            buffer.rgbs[idx]
        };
    }

    static __device__ __forceinline__ WorldEval3D loadWithPrecompute(const Buffer &buffer, long idx) {
        return WorldEval3D::load(buffer, idx);
    }

    static __device__ __forceinline__ WorldEval3D zero() {
        return {
            {0.f, 0.f, 0.f},
            0.f,
            0.f,
            {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f},
            {0.f, 0.f, 0.f},
        };
    }

    __device__ __forceinline__ void operator+=(const WorldEval3D &other) {
        pos += other.pos;
        size += other.size;
        depth += other.depth;
        #pragma unroll
        for (int i = 0; i < 8; i++)
            densities[i] += other.densities[i];
        rgb += other.rgb;
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

    static __device__ void atomicAddGradientToBuffer(const WorldEval3D &grad, Buffer &buffer, long idx) {
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

    __device__ __forceinline__ float evaluate_alpha(float3 ray_o, float3 ray_d) {
        return evaluate_alpha_voxel(pos, size, &densities, ray_o, ray_d);
    }

    __device__ __forceinline__ WorldEval3D evaluate_alpha_vjp(
        float3 ray_o, float3 ray_d, float v_alpha,
        float3 &v_ray_o, float3 &v_ray_d
    ) {
        WorldEval3D v_splat = WorldEval3D::zero();
        evaluate_alpha_voxel_vjp(
            pos, size, &densities,
            ray_o, ray_d, v_alpha,
            &v_splat.densities,
            &v_ray_o, &v_ray_d
        );
        return v_splat;
    }

    __device__ __forceinline__ VoxelPrimitive::RenderOutput evaluate_color(float3 ray_o, float3 ray_d) {
        float3 out_rgb; float out_depth;
        evaluate_color_voxel(
            pos, size, rgb, ray_o, ray_d,
            &out_rgb, &out_depth
        );
        return {out_rgb, out_depth};
    }

    __device__ __forceinline__ WorldEval3D evaluate_color_vjp(
        float3 ray_o, float3 ray_d, VoxelPrimitive::RenderOutput v_render,
        float3 &v_ray_o, float3 &v_ray_d
    ) {
        WorldEval3D v_splat = WorldEval3D::zero();
        evaluate_color_voxel_vjp(
            pos, size, rgb, ray_o, ray_d,
            v_render.rgb, v_render.depth,
            &v_splat.rgb,
            &v_ray_o, &v_ray_d
        );
        return v_splat;
    }

#endif  // #ifdef __CUDACC__

};


#ifdef __CUDACC__

inline __device__ void VoxelPrimitive::project_persp_eval3d(
    VoxelPrimitive::World world, VoxelPrimitive::FwdProjCamera cam,
    VoxelPrimitive::WorldEval3D& proj, int4& aabb
) {
    projection_voxel_eval3d_persp(
        world.pos, world.size, &world.densities, &world.sh_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, &cam.dist_coeffs,
        cam.width, cam.height, cam.near_plane, cam.far_plane,
        &aabb, &proj.depth, &proj.rgb
    );
}

inline __device__ void VoxelPrimitive::project_fisheye_eval3d(
    VoxelPrimitive::World world, VoxelPrimitive::FwdProjCamera cam,
    VoxelPrimitive::WorldEval3D& proj, int4& aabb
) {
    projection_voxel_eval3d_fisheye(
        world.pos, world.size, &world.densities, &world.sh_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, &cam.dist_coeffs,
        cam.width, cam.height, cam.near_plane, cam.far_plane,
        &aabb, &proj.depth, &proj.rgb
    );
}

inline __device__ void VoxelPrimitive::project_persp_eval3d_vjp(
    VoxelPrimitive::World world, VoxelPrimitive::BwdProjCamera cam,
    VoxelPrimitive::WorldEval3D v_proj,
    VoxelPrimitive::World& v_world, float3x3 &v_R, float3 &v_t
) {
    projection_voxel_eval3d_persp_vjp(
        world.pos, world.size, &world.densities, &world.sh_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, &cam.dist_coeffs,
        cam.width, cam.height,
        v_proj.rgb,
        &v_world.densities, &v_world.sh_coeffs,
        &v_R, &v_t
    );
}

inline __device__ void VoxelPrimitive::project_fisheye_eval3d_vjp(
    VoxelPrimitive::World world, VoxelPrimitive::BwdProjCamera cam,
    VoxelPrimitive::WorldEval3D v_proj,
    VoxelPrimitive::World& v_world, float3x3 &v_R, float3 &v_t
) {
    projection_voxel_eval3d_fisheye_vjp(
        world.pos, world.size, &world.densities, &world.sh_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, &cam.dist_coeffs,
        cam.width, cam.height,
        v_proj.rgb,
        &v_world.densities, &v_world.sh_coeffs,
        &v_R, &v_t
    );
}

#endif  // #ifdef __CUDACC__
