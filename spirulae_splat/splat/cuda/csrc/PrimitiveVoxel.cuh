#pragma once

#ifdef __CUDACC__
#include "generated/slang.cuh"
namespace SlangVoxel {
#include "generated/set_namespace.cuh"
#include "generated/primitive_voxel.cuh"
}
#endif

#include "Primitive.cuh"

#if 0
struct VoxelPrimitive {
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

struct VoxelPrimitive::World {

    float3 pos;
    float size;
    FixedArray<float, 8> densities;
    FixedArray<float3, 16> sh_coeffs;

    #ifndef NO_TORCH
    typedef std::tuple<std::optional<at::Tensor>, std::optional<at::Tensor>, at::Tensor, at::Tensor> TensorTuple;
    #endif

    struct Buffer;

    #ifndef NO_TORCH
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
                zeros_like_tensor(features_dc),
                zeros_like_tensor(features_sh)  // TODO: exclude SH when is_hess_diag
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
    #endif

    struct Buffer : public TensorArray<4> {
        uint num_sh;

        using TensorArray<4>::TensorArray;

        Buffer() : num_sh(0) {}

        #ifndef NO_TORCH
        Buffer(const Tensor& tensors)
            : TensorArray<4>(std::vector<std::optional<at::Tensor>>{
                tensors.pos_size, tensors.densities, tensors.features_dc, tensors.features_sh
            }) {
            num_sh = tensors.features_sh.has_value() ? tensors.features_sh.value().size(-2) : 0;
        }
        #endif

        #ifdef __CUDACC__
        __forceinline__ __device__ bool hasPosSize() const { return _data[0] != nullptr; }
        __forceinline__ __device__ bool hasDensities() const { return _data[1] != nullptr; }
        __forceinline__ __device__ float4& pos_size(int64_t i)
            { return *reinterpret_cast<float4*>(&_data[0][4*i]); }
        __forceinline__ __device__ float4& densities(int64_t i)
            { return *reinterpret_cast<float4*>(&_data[1][4*i]); }
        __forceinline__ __device__ float3& features_dc(int64_t i)
            { return *reinterpret_cast<float3*>(&_data[2][3*i]); }
        __forceinline__ __device__ float3& features_sh(int64_t i, int64_t j)
            { return *reinterpret_cast<float3*>(&_data[3][_strides[3]*i+j]); }
        #endif
    };

#ifdef __CUDACC__

    static __device__ World load(const Buffer &buffer, long idx) {
        float4 pos_size = buffer.hasPosSize() ? buffer.pos_size(idx) : make_float4(0,0,0,0);
        float4 d0 = buffer.hasDensities() ? buffer.densities(2*idx) : make_float4(0,0,0,0);
        float4 d1 = buffer.hasDensities() ? buffer.densities(2*idx+1) : make_float4(0,0,0,0);
        World world = {
            {pos_size.x, pos_size.y, pos_size.z},
            pos_size.w,
            {d0.x, d0.y, d0.z, d0.w, d1.x, d1.y, d1.z, d1.w}
        };
        world.sh_coeffs[0] = buffer.features_dc(idx);
        for (int i = 0; i < 15; i++)
            world.sh_coeffs[i+1] = i < buffer.num_sh ?
                buffer.features_sh(idx, i) : make_float3(0);
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

    __device__ void saveParamsToBuffer(Buffer &buffer, long idx) {
        if (buffer.hasPosSize())
            buffer.pos_size(idx) = {pos.x, pos.y, pos.z, size};
        if (buffer.hasDensities()) {
            buffer.densities(2*idx) = {densities[0], densities[1], densities[2], densities[3]};
            buffer.densities(2*idx+1) = {densities[4], densities[5], densities[6], densities[7]};
        }
        buffer.features_dc(idx) = sh_coeffs[0];
        for (int i = 0; i < buffer.num_sh; i++)
            buffer.features_sh(idx, i) = sh_coeffs[i+1];
    }

    __device__ void atomicAddGradientToBuffer(Buffer &buffer, long idx) {
        if (buffer.hasPosSize())
            atomicAddFVec(&buffer.pos_size(idx), {pos.x, pos.y, pos.z, size});
        if (buffer.hasDensities()) {
            atomicAddFVec(&buffer.densities(2*idx), {densities[0], densities[1], densities[2], densities[3]});
            atomicAddFVec(&buffer.densities(2*idx+1), {densities[4], densities[5], densities[6], densities[7]});
        }
        atomicAddFVec(&buffer.features_dc(idx), sh_coeffs[0]);
        for (int i = 0; i < buffer.num_sh; i++)
            atomicAddFVec(&buffer.features_sh(idx, i), sh_coeffs[i+1]);
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

    #ifndef NO_TORCH
    typedef std::tuple<std::optional<at::Tensor>, at::Tensor> TensorTupleProj;
    typedef std::tuple<std::optional<at::Tensor>, std::optional<at::Tensor>, std::optional<at::Tensor>, at::Tensor> TensorTuple;
    #endif

    struct Buffer;

    #ifndef NO_TORCH
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
                (std::optional<at::Tensor>)zeros_like_tensor(densities.value()),
                zeros_like_tensor(rgbs)
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
    #endif

    struct Buffer : public TensorArray<4> {
        long size;

        using TensorArray<4>::TensorArray;

        Buffer() : size(0) {}

        #ifndef NO_TORCH
        Buffer(const Tensor& tensors)
            : TensorArray<4>(std::vector<std::optional<at::Tensor>>{
                tensors.pos_size, tensors.depths, tensors.densities, tensors.rgbs
            }) {
            size = tensors.rgbs.size(-2);
        }
        #endif

        #ifdef __CUDACC__
        __forceinline__ __device__ bool hasPosSize() const { return _data[0] != nullptr; }
        __forceinline__ __device__ bool hasDepths() const { return _data[1] != nullptr; }
        __forceinline__ __device__ bool hasDensities() const { return _data[2] != nullptr; }
        __forceinline__ __device__ float4& pos_size(int64_t i)
            { return *reinterpret_cast<float4*>(&_data[0][4*i]); }
        __forceinline__ __device__ float& depths(int64_t i)
            { return *reinterpret_cast<float*>(&_data[1][i]); }
        __forceinline__ __device__ float4& densities(int64_t i)
            { return *reinterpret_cast<float4*>(&_data[2][4*i]); }
        __forceinline__ __device__ float3& rgbs(int64_t i)
            { return *reinterpret_cast<float3*>(&_data[3][3*i]); }
        #endif
    };

#ifdef __CUDACC__

    static __device__ Screen load(const Buffer &buffer, long idx, const uint32_t* gaussian_ids) {
        uint32_t idx0 = gaussian_ids ? gaussian_ids[idx] : idx % buffer.size;
        float4 pos_size = buffer.hasPosSize() ? buffer.pos_size(idx0) : make_float4(0, 0, 0, 0);
        float4 d0 = buffer.hasDensities() ? buffer.densities(2*idx0) : make_float4(0, 0, 0, 0);
        float4 d1 = buffer.hasDensities() ? buffer.densities(2*idx0+1) : make_float4(0, 0, 0, 0);
        return {
            {pos_size.x, pos_size.y, pos_size.z},
            pos_size.w,
            buffer.hasDepths() ? buffer.depths(idx) : 0.0f,
            {d0.x, d0.y, d0.z, d0.w, d1.x, d1.y, d1.z, d1.w},
            buffer.rgbs(idx)
        };
    }

    static __device__ __forceinline__ Screen loadWithPrecompute(const Buffer &buffer, long idx, const uint32_t* gaussian_ids) {
        return Screen::load(buffer, idx, gaussian_ids);
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

    __device__ __forceinline__ void precomputeBackward(Screen& grad) const {
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

    __device__ __forceinline__ void addGaussNewtonHessianDiagonal(const Screen &grad, float weight=1.0f) {
        pos += grad.pos * grad.pos * weight;
        size += grad.size * grad.size * weight;
        depth += grad.depth * grad.depth * weight;
        #pragma unroll
        for (int i = 0; i < 8; i++)
            densities[i] += grad.densities[i] * grad.densities[i] * weight;
        rgb += grad.rgb * grad.rgb * weight;
    }

    __device__ void saveParamsToBuffer(Buffer &buffer, long idx, const uint32_t* gaussian_ids) {
        if (buffer.hasPosSize() && idx < buffer.size)
            buffer.pos_size(idx) = {pos.x, pos.y, pos.z, size};
        if (buffer.hasDepths())
            buffer.depths(idx) = depth;
        if (buffer.hasDensities() && idx < buffer.size) {
            buffer.densities(2*idx) = {densities[0], densities[1], densities[2], densities[3]};
            buffer.densities(2*idx+1) = {densities[4], densities[5], densities[6], densities[7]};
        }
        buffer.rgbs(idx) = rgb;
    }

    __device__ void atomicAddToBuffer(Buffer &buffer, long idx, const uint32_t* gaussian_ids) const {
        uint32_t idx0 = gaussian_ids ? gaussian_ids[idx] : idx % buffer.size;
        if (buffer.hasPosSize())
            atomicAddFVec(&buffer.pos_size(idx0), {pos.x, pos.y, pos.z, size});
        if (buffer.hasDepths())
            atomicAddFVec(&buffer.depths(idx), depth);
        if (buffer.hasDensities()) {
            atomicAddFVec(&buffer.densities(2*idx0), {densities[0], densities[1], densities[2], densities[3]});
            atomicAddFVec(&buffer.densities(2*idx0+1), {densities[4], densities[5], densities[6], densities[7]});
        }
        atomicAddFVec(&buffer.rgbs(idx), rgb);
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

    __device__ __forceinline__ RenderOutput evaluate_color(float3 ray_o, float3 ray_d) {
        float3 out_rgb; float out_depth;
        SlangVoxel::evaluate_color_voxel(
            pos, size, densities, rgb, ray_o, ray_d,
            &out_rgb, &out_depth
        );
        return {out_rgb, out_depth, make_float3(0.0f)};
    }

    __device__ __forceinline__ Screen evaluate_color_vjp(
        float3 ray_o, float3 ray_d, RenderOutput v_render,
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
    VoxelPrimitive::World world, ProjCamera cam,
    VoxelPrimitive::Screen& proj, float4& aabb
) {
    SlangVoxel::projection_voxel_eval3d_persp(
        world.pos, world.size, world.densities, world.sh_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
        cam.width, cam.height,
        &aabb, &proj.depth, &proj.rgb
    );
    proj.pos = world.pos, proj.size = world.size, proj.densities = world.densities;
}

inline __device__ void VoxelPrimitive::project_fisheye(
    VoxelPrimitive::World world, ProjCamera cam,
    VoxelPrimitive::Screen& proj, float4& aabb
) {
    SlangVoxel::projection_voxel_eval3d_fisheye(
        world.pos, world.size, world.densities, world.sh_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
        cam.width, cam.height,
        &aabb, &proj.depth, &proj.rgb
    );
    proj.pos = world.pos, proj.size = world.size, proj.densities = world.densities;
}

inline __device__ void VoxelPrimitive::project_persp_vjp(
    VoxelPrimitive::World world, ProjCamera cam,
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
    VoxelPrimitive::World world, ProjCamera cam,
    VoxelPrimitive::Screen v_proj, VoxelPrimitive::Screen vr_proj, VoxelPrimitive::Screen h_proj,
    VoxelPrimitive::World& v_world, float3x3 &v_R, float3 &v_t,
    float3 &vr_world_pos, float3 &h_world_pos
) {}  // TODO

inline __device__ void VoxelPrimitive::project_persp_vjp(
    VoxelPrimitive::World world, ProjCamera cam,
    VoxelPrimitive::Screen v_proj, VoxelPrimitive::Screen vr_proj, VoxelPrimitive::Screen h_proj,
    VoxelPrimitive::World& v_world, float3x3 &v_R, float3 &v_t,
    VoxelPrimitive::World& vr_world, VoxelPrimitive::World& h_world
) {}  // TODO

inline __device__ void VoxelPrimitive::project_fisheye_vjp(
    VoxelPrimitive::World world, ProjCamera cam,
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
    VoxelPrimitive::World world, ProjCamera cam,
    VoxelPrimitive::Screen v_proj, VoxelPrimitive::Screen vr_proj, VoxelPrimitive::Screen h_proj,
    VoxelPrimitive::World& v_world, float3x3 &v_R, float3 &v_t,
    float3 &vr_world_pos, float3 &h_world_pos
) {}  // TODO

inline __device__ void VoxelPrimitive::project_fisheye_vjp(
    VoxelPrimitive::World world, ProjCamera cam,
    VoxelPrimitive::Screen v_proj, VoxelPrimitive::Screen vr_proj, VoxelPrimitive::Screen h_proj,
    VoxelPrimitive::World& v_world, float3x3 &v_R, float3 &v_t,
    VoxelPrimitive::World& vr_world, VoxelPrimitive::World& h_world
) {}  // TODO

#endif  // #ifdef __CUDACC__
#endif
