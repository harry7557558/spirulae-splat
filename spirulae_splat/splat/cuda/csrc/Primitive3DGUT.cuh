#pragma once

#ifdef __CUDACC__
#define TensorView _Slang_TensorView
#include "generated/primitive.cuh"
#undef TensorView
#endif

#include "types.cuh"
#include "common.cuh"

#include <tuple>


struct Base3DGUT {
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

    struct BwdProjCamera {
        float3x3 R;
        float3 t;
        float fx, fy, cx, cy;
        uint width, height;
        CameraDistortionCoeffs dist_coeffs;
    };

#endif  // #ifdef __CUDACC__

};

struct Base3DGUT::World {
    float3 mean;
    float4 quat;
    float3 scale;
    float opacity;
};


struct Base3DGUT::RenderOutput {

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
        return (rgb.x * other.rgb.x + rgb.y * other.rgb.y + rgb.z * other.rgb.z) + depth * other.depth;
    }

    __device__ void saveParamsToBuffer(Buffer &buffer, long idx) {
        buffer.rgbs[idx] = rgb;
        buffer.depths[idx] = depth;
    }

#endif  // #ifdef __CUDACC__
};


struct Base3DGUT::Screen {

    // from world
    float3 mean;
    float4 quat;
    // from screen
    float depth;
    float3 scale;
    float opacity;
    float3 rgb;
    // precompute
#ifdef __CUDACC__
    float3x3 iscl_rot;
#endif

    typedef std::tuple<at::Tensor, at::Tensor, at::Tensor, at::Tensor> TensorTupleProj;
    typedef std::tuple<at::Tensor, at::Tensor, at::Tensor, at::Tensor, at::Tensor, at::Tensor> TensorTuple;

    struct Buffer;

    struct Tensor {
        bool hasWorld;
        at::Tensor means;
        at::Tensor quats;
        at::Tensor depths;
        at::Tensor scales;
        at::Tensor opacities;
        at::Tensor rgbs;

        Tensor() {}

        Tensor(const TensorTuple& splats) : hasWorld(true) {
            means = std::get<0>(splats);
            quats = std::get<1>(splats);
            depths = std::get<2>(splats);
            scales = std::get<3>(splats);
            opacities = std::get<4>(splats);
            rgbs = std::get<5>(splats);
        }

        Tensor(const TensorTupleProj& splats) : hasWorld(false) {
            depths = std::get<0>(splats);
            scales = std::get<1>(splats);
            opacities = std::get<2>(splats);
            rgbs = std::get<3>(splats);
        }

        TensorTupleProj tupleProjFwd() const {
            return std::make_tuple(depths, scales, opacities, rgbs);
        }

        TensorTuple tupleProjFwdPacked() const {
            return std::make_tuple(means, quats, depths, scales, opacities, rgbs);
        }

        TensorTuple tupleRasterBwd() const {
            return std::make_tuple(means, quats, depths, scales, opacities, rgbs);
        }

        static TensorTupleProj allocProjFwd(long C, long N, c10::TensorOptions opt) {
            return std::make_tuple(
                at::empty({C, N}, opt),
                at::empty({C, N, 3}, opt),
                at::empty({C, N}, opt),
                at::empty({C, N, 3}, opt)
            );
        }

        static TensorTuple allocProjFwdPacked(long nnz, c10::TensorOptions opt) {
            return std::make_tuple(
                at::empty({nnz, 3}, opt),
                at::empty({nnz, 4}, opt),
                at::empty({nnz}, opt),
                at::empty({nnz, 3}, opt),
                at::empty({nnz}, opt),
                at::empty({nnz, 3}, opt)
            );
        }

        Tensor allocRasterBwd() const {
            if (!hasWorld)
                throw std::runtime_error("!hasWorld");
            Tensor result = Tensor(std::make_tuple(
                at::zeros_like(means),
                at::zeros_like(quats),
                at::zeros_like(depths),
                at::zeros_like(scales),
                at::zeros_like(opacities),
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

        Buffer buffer() { return Buffer(*this); }
    };

    struct Buffer {
        float3* __restrict__ means;
        float4* __restrict__ quats;
        float* __restrict__ depths;
        float3* __restrict__ scales;
        float* __restrict__ opacities;
        float3* __restrict__ rgbs;
        long size;

        Buffer() {}  // uninitialized

        Buffer(const Tensor& tensors) {
            DEVICE_GUARD(tensors.depths);
            if (tensors.hasWorld) {
                CHECK_INPUT(tensors.means);
                CHECK_INPUT(tensors.quats);
                means = (float3*)tensors.means.data_ptr<float>();
                quats = (float4*)tensors.quats.data_ptr<float>();
            } else { means = nullptr; quats = nullptr; };
            CHECK_INPUT(tensors.depths);
            CHECK_INPUT(tensors.scales);
            CHECK_INPUT(tensors.opacities);
            CHECK_INPUT(tensors.rgbs);
            depths = tensors.depths.data_ptr<float>();
            scales = (float3*)tensors.scales.data_ptr<float>();
            opacities = tensors.opacities.data_ptr<float>();
            rgbs = (float3*)tensors.rgbs.data_ptr<float>();
            size = tensors.hasWorld ?
                tensors.quats.numel() / 4 : tensors.depths.numel();
        }
    };

#ifdef __CUDACC__

    static __device__ Screen load(const Buffer &buffer, long idx) {
        return {
            buffer.means ? buffer.means[idx % buffer.size] : make_float3(0.f),
            buffer.quats ? buffer.quats[idx % buffer.size] : make_float4(0.f),
            buffer.depths[idx],
            buffer.scales[idx],
            buffer.opacities[idx],
            buffer.rgbs[idx],
        };
    }

    static __device__ Screen loadWithPrecompute(const Buffer &buffer, long idx) {
        Screen splat = Screen::load(buffer, idx);
        splat.iscl_rot = compute_3dgut_iscl_rot(splat.quat, splat.scale);
        return splat;
    }

    static __device__ __forceinline__ Screen zero() {
        return {
            {0.f, 0.f, 0.f},
            {0.f, 0.f, 0.f, 0.f},
            0.0f,
            {0.f, 0.f, 0.f},
            0.0f,
            {0.f, 0.f, 0.f},
            {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f},
        };
    }

    __device__ __forceinline__ void addGradient(const Screen &grad, float weight=1.0f) {
        mean += grad.mean * weight;
        quat += grad.quat * weight;
        depth += grad.depth * weight;
        scale += grad.scale * weight;
        opacity += grad.opacity * weight;
        rgb += grad.rgb * weight;
        #pragma unroll
        for (int i = 0; i < 3; ++i)
            iscl_rot[i] = iscl_rot[i] + grad.iscl_rot[i] * weight;
    }

    __device__ __forceinline__ void addGaussNewtonHessianDiagonal(const Screen &grad, float weight=1.0f) {
        float4 v_quat; float3 v_scale;
        compute_3dgut_iscl_rot_vjp(quat, scale, grad.iscl_rot, &v_quat, &v_scale);
        mean += grad.mean * grad.mean * weight;
        quat += (grad.quat + v_quat) * (grad.quat + v_quat) * weight;
        depth += grad.depth * grad.depth * weight;
        scale += (grad.scale + v_scale) * (grad.scale + v_scale) * weight;
        opacity += grad.opacity * grad.opacity * weight;
        rgb += grad.rgb * grad.rgb * weight;
    }

    __device__ void saveParamsToBuffer(Buffer &buffer, long idx) {
        if (buffer.means) buffer.means[idx % buffer.size] = mean;
        if (buffer.quats) buffer.quats[idx % buffer.size] = quat;
        buffer.depths[idx] = depth;
        buffer.scales[idx] = scale;
        buffer.opacities[idx] = opacity;
        buffer.rgbs[idx] = rgb;
        // iscl_rot is not saved
    }

    __device__ void atomicAddGradientToBuffer(const Screen &grad, Buffer &buffer, long idx) const {
        float4 v_quat; float3 v_scale;
        compute_3dgut_iscl_rot_vjp(quat, scale, grad.iscl_rot, &v_quat, &v_scale);
        atomicAddFVec(buffer.means + idx % buffer.size, grad.mean);
        atomicAddFVec(buffer.quats + idx % buffer.size, grad.quat + v_quat);
        atomicAddFVec(buffer.depths + idx, grad.depth);
        atomicAddFVec(buffer.scales + idx, grad.scale + v_scale);
        atomicAddFVec(buffer.opacities + idx, grad.opacity);
        atomicAddFVec(buffer.rgbs + idx, grad.rgb);
    }

    __device__ void atomicAddGaussNewtonHessianDiagonalToBuffer(const Screen &grad, Buffer &buffer, long idx, float weight=1.0f) const {
        float4 v_quat; float3 v_scale;
        compute_3dgut_iscl_rot_vjp(quat, scale, grad.iscl_rot, &v_quat, &v_scale);
        atomicAddFVec(buffer.means + idx % buffer.size, grad.mean * grad.mean * weight);
        atomicAddFVec(buffer.quats + idx % buffer.size, (grad.quat + v_quat) * (grad.quat + v_quat) * weight);
        atomicAddFVec(buffer.depths + idx, grad.depth * grad.depth * weight);
        atomicAddFVec(buffer.scales + idx, (grad.scale + v_scale) * (grad.scale + v_scale) * weight);
        atomicAddFVec(buffer.opacities + idx, grad.opacity * grad.opacity * weight);
        atomicAddFVec(buffer.rgbs + idx, grad.rgb * grad.rgb * weight);
    }

    __device__ __forceinline__ float evaluate_alpha(float3 ray_o, float3 ray_d) {
        if (dot(mean-ray_o, ray_d) <= 0.0f)
            return 0.0;
        return evaluate_alpha_3dgs(
            mean, iscl_rot, opacity,
            ray_o, ray_d
        );
    }

    __device__ __forceinline__ Screen evaluate_alpha_vjp(
        float3 ray_o, float3 ray_d, float v_alpha,
        float3 &v_ray_o, float3 &v_ray_d
    ) {
        Screen v_splat = Screen::zero();
        if (dot(mean-ray_o, ray_d) <= 0.0f) {
            v_ray_o = v_ray_d = make_float3(0.f);
            return v_splat;
        }
        evaluate_alpha_3dgs_vjp(
            mean, iscl_rot, opacity,
            ray_o, ray_d, v_alpha,
            &v_splat.mean, &v_splat.iscl_rot, &v_splat.opacity,
            &v_ray_o, &v_ray_d
        );
        return v_splat;
    }

    __device__ __forceinline__ Base3DGUT::RenderOutput evaluate_color(float3 ray_o, float3 ray_d) {
        // return {rgb, depth};
        float3 out_rgb; float out_depth;
        evaluate_color_3dgs(
            mean, iscl_rot, opacity, rgb,
            ray_o, ray_d, &out_rgb, &out_depth
        );
        return {out_rgb, out_depth};
    }

    __device__ __forceinline__ Screen evaluate_color_vjp(
        float3 ray_o, float3 ray_d, Base3DGUT::RenderOutput v_render,
        float3 &v_ray_o, float3 &v_ray_d
    ) {
        Screen v_splat = Screen::zero();
    #if 0
        v_splat.rgb = v_render.rgb;
        v_splat.depth = v_render.depth;
        v_ray_o = make_float3(0.f);
        v_ray_d = make_float3(0.f);
    #else
        evaluate_color_3dgs_vjp(
            mean, iscl_rot, opacity, rgb,
            ray_o, ray_d,
            v_render.rgb, v_render.depth,
            &v_splat.mean, &v_splat.iscl_rot, &v_splat.opacity, &v_splat.rgb,
            &v_ray_o, &v_ray_d
        );
    #endif
        return v_splat;
    }

#endif  // #ifdef __CUDACC__

};



struct Vanilla3DGUT : public Base3DGUT {
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
        Screen v_proj, Screen h_proj,
        World& v_world, float3x3 &v_R, float3 &v_t,
        float3 &h_world_pos
    );

    inline static __device__ void project_fisheye_vjp(
        World world, BwdProjCamera cam,
        Screen v_proj,
        World& v_world, float3x3 &v_R, float3 &v_t
    );

    inline static __device__ void project_fisheye_vjp(
        World world, BwdProjCamera cam,
        Screen v_proj, Screen h_proj,
        World& v_world, float3x3 &v_R, float3 &v_t,
        float3 &h_world_pos
    );

#endif  // #ifdef __CUDACC__

};


struct Vanilla3DGUT::World : public Base3DGUT::World {

#ifdef __CUDACC__
    FixedArray<float3, 16> sh_coeffs;
#endif

    typedef std::tuple<at::Tensor, at::Tensor, at::Tensor, at::Tensor, at::Tensor, at::Tensor> TensorTuple;

    struct Buffer;

    struct Tensor {
        at::Tensor means;
        at::Tensor quats;
        at::Tensor scales;
        at::Tensor opacities;
        at::Tensor features_dc;
        at::Tensor features_sh;

        Tensor() {}

        Tensor(const TensorTuple& splats) {
            means = std::get<0>(splats);
            quats = std::get<1>(splats);
            scales = std::get<2>(splats);
            opacities = std::get<3>(splats);
            features_dc = std::get<4>(splats);
            features_sh = std::get<5>(splats);
        }

        TensorTuple tuple() const {
            return std::make_tuple(means, quats, scales, opacities, features_dc, features_sh);
        }

        Tensor zeros_like() const {
            return Tensor(std::make_tuple(
                at::zeros_like(means),
                at::zeros_like(quats),
                at::zeros_like(scales),
                at::zeros_like(opacities),
                at::zeros_like(features_dc),
                at::zeros_like(features_sh)
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
        float3* __restrict__ features_dc;
        float3* __restrict__ features_sh;
        uint num_sh;

        Buffer() {}

        Buffer(const Tensor& tensors) {
            DEVICE_GUARD(tensors.means);
            CHECK_INPUT(tensors.means);
            CHECK_INPUT(tensors.quats);
            CHECK_INPUT(tensors.scales);
            CHECK_INPUT(tensors.opacities);
            CHECK_INPUT(tensors.features_dc);
            CHECK_INPUT(tensors.features_sh);
            means = (float3*)tensors.means.data_ptr<float>();
            quats = (float4*)tensors.quats.data_ptr<float>();
            scales = (float3*)tensors.scales.data_ptr<float>();
            opacities = tensors.opacities.data_ptr<float>();
            features_dc = (float3*)tensors.features_dc.data_ptr<float>();
            features_sh = (float3*)tensors.features_sh.data_ptr<float>();
            num_sh = tensors.features_sh.size(-2);
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
        world.sh_coeffs[0] = buffer.features_dc[idx];
        for (int i = 0; i < 15; i++)
            world.sh_coeffs[i+1] = i < buffer.num_sh ?
                buffer.features_sh[idx*buffer.num_sh+i] : make_float3(0);
        return world;
    }

    static __device__ __forceinline__ World zero() {
        World world;
        world.mean = {0.f, 0.f, 0.f};
        world.quat = {0.f, 0.f, 0.f, 0.f};
        world.scale = {0.f, 0.f, 0.f};
        world.opacity = 0.f;
        for (int i = 0; i < 16; i++)
            world.sh_coeffs[i] = make_float3(0);
        return world;
    }

    template<typename Partition>
    __device__ void reduce(Partition& partition) {
        warpSum(mean, partition);
        warpSum(quat, partition);
        warpSum(scale, partition);
        warpSum(opacity, partition);
        for (int i = 0; i < 16; i++)
            warpSum(sh_coeffs[i], partition);
    }
    
    __device__ void saveParamsToBuffer(Buffer &buffer, long idx) {
        buffer.means[idx] = mean;
        buffer.quats[idx] = quat;
        buffer.scales[idx] = scale;
        buffer.opacities[idx] = opacity;
        buffer.features_dc[idx] = sh_coeffs[0];
        for (int i = 0; i < buffer.num_sh; i++)
            buffer.features_sh[idx*buffer.num_sh + i] = sh_coeffs[i+1];
    }

    __device__ void atomicAddGradientToBuffer(Buffer &buffer, long idx) {
        atomicAddFVec(buffer.means + idx, mean);
        atomicAddFVec(buffer.quats + idx, quat);
        atomicAddFVec(buffer.scales + idx, scale);
        atomicAddFVec(buffer.opacities + idx, opacity);
        atomicAddFVec(buffer.features_dc + idx, sh_coeffs[0]);
        for (int i = 0; i < buffer.num_sh; i++)
            atomicAddFVec(buffer.features_sh + idx*buffer.num_sh + i, sh_coeffs[i+1]);
    }

#endif  // #ifdef __CUDACC__
};


#ifdef __CUDACC__

inline __device__ void Vanilla3DGUT::project_persp(
    Vanilla3DGUT::World world, Vanilla3DGUT::FwdProjCamera cam,
    Vanilla3DGUT::Screen& proj, int4& aabb
) {
    float2 xy;
    projection_3dgut_persp(
        false,
        world.mean, world.quat, world.scale, world.opacity, &world.sh_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, &cam.dist_coeffs,
        cam.width, cam.height, cam.near_plane, cam.far_plane,
        &aabb, &xy, &proj.depth, &proj.scale, &proj.opacity, &proj.rgb
    );
    proj.mean = world.mean, proj.quat = world.quat;
}

inline __device__ void Vanilla3DGUT::project_fisheye(
    Vanilla3DGUT::World world, Vanilla3DGUT::FwdProjCamera cam,
    Vanilla3DGUT::Screen& proj, int4& aabb
) {
    float2 xy;
    projection_3dgut_fisheye(
        false,
        world.mean, world.quat, world.scale, world.opacity, &world.sh_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, &cam.dist_coeffs,
        cam.width, cam.height, cam.near_plane, cam.far_plane,
        &aabb, &xy, &proj.depth, &proj.scale, &proj.opacity, &proj.rgb
    );
    proj.mean = world.mean, proj.quat = world.quat;
}

inline __device__ void Vanilla3DGUT::project_persp_vjp(
    Vanilla3DGUT::World world, Vanilla3DGUT::BwdProjCamera cam,
    Vanilla3DGUT::Screen v_proj,
    Vanilla3DGUT::World& v_world, float3x3 &v_R, float3 &v_t
) {
    projection_3dgut_persp_vjp(
        false,
        world.mean, world.quat, world.scale, world.opacity, &world.sh_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, &cam.dist_coeffs,
        cam.width, cam.height,
        make_float2(0), v_proj.depth, v_proj.scale, v_proj.opacity, v_proj.rgb,
        &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity, &v_world.sh_coeffs,
        &v_R, &v_t
    );
    v_world.mean = v_proj.mean, v_world.quat = v_proj.quat;
}

inline __device__ void Vanilla3DGUT::project_persp_vjp(
    Vanilla3DGUT::World world, Vanilla3DGUT::BwdProjCamera cam,
    Vanilla3DGUT::Screen v_proj, Vanilla3DGUT::Screen h_proj,
    Vanilla3DGUT::World& v_world, float3x3 &v_R, float3 &v_t,
    float3& h_world_pos
) {
    projection_3dgut_persp_vjp(
        false,
        world.mean, world.quat, world.scale, world.opacity, &world.sh_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, &cam.dist_coeffs,
        cam.width, cam.height,
        make_float2(0), v_proj.depth, v_proj.scale, v_proj.opacity, v_proj.rgb,
        make_float2(0), h_proj.depth, h_proj.scale, h_proj.opacity, h_proj.rgb,
        &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity, &v_world.sh_coeffs,
        &v_R, &v_t, &h_world_pos
    );
    v_world.mean = v_proj.mean, v_world.quat = v_proj.quat;
}

inline __device__ void Vanilla3DGUT::project_fisheye_vjp(
    Vanilla3DGUT::World world, Vanilla3DGUT::BwdProjCamera cam,
    Vanilla3DGUT::Screen v_proj,
    Vanilla3DGUT::World& v_world, float3x3 &v_R, float3 &v_t
) {
    projection_3dgut_fisheye_vjp(
        false,
        world.mean, world.quat, world.scale, world.opacity, &world.sh_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, &cam.dist_coeffs,
        cam.width, cam.height,
        make_float2(0), v_proj.depth, v_proj.scale, v_proj.opacity, v_proj.rgb,
        &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity, &v_world.sh_coeffs,
        &v_R, &v_t
    );
    v_world.mean = v_proj.mean, v_world.quat = v_proj.quat;
}

inline __device__ void Vanilla3DGUT::project_fisheye_vjp(
    Vanilla3DGUT::World world, Vanilla3DGUT::BwdProjCamera cam,
    Vanilla3DGUT::Screen v_proj, Vanilla3DGUT::Screen h_proj,
    Vanilla3DGUT::World& v_world, float3x3 &v_R, float3 &v_t,
    float3& h_world_pos
) {
    projection_3dgut_fisheye_vjp(
        false,
        world.mean, world.quat, world.scale, world.opacity, &world.sh_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, &cam.dist_coeffs,
        cam.width, cam.height,
        make_float2(0), v_proj.depth, v_proj.scale, v_proj.opacity, v_proj.rgb,
        make_float2(0), h_proj.depth, h_proj.scale, h_proj.opacity, h_proj.rgb,
        &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity, &v_world.sh_coeffs,
        &v_R, &v_t, &h_world_pos
    );
    v_world.mean = v_proj.mean, v_world.quat = v_proj.quat;
}

#endif  // #ifdef __CUDACC__




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
        Screen v_proj, Screen h_proj,
        World& v_world, float3x3 &v_R, float3 &v_t,
        float3 &h_world_pos
    );

    inline static __device__ void project_fisheye_vjp(
        World world, BwdProjCamera cam,
        Screen v_proj,
        World& v_world, float3x3 &v_R, float3 &v_t
    );

    inline static __device__ void project_fisheye_vjp(
        World world, BwdProjCamera cam,
        Screen v_proj, Screen h_proj,
        World& v_world, float3x3 &v_R, float3 &v_t,
        float3 &h_world_pos
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

        Tensor zeros_like() const {
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
    projection_3dgut_sv_persp(
        false,
        world.mean, world.quat, world.scale, world.opacity, &world.sv_sites, &world.sv_colors,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, &cam.dist_coeffs,
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
    projection_3dgut_sv_fisheye(
        false,
        world.mean, world.quat, world.scale, world.opacity, &world.sv_sites, &world.sv_colors,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, &cam.dist_coeffs,
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
    projection_3dgut_sv_persp_vjp(
        false,
        world.mean, world.quat, world.scale, world.opacity, &world.sv_sites, &world.sv_colors,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, &cam.dist_coeffs,
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
    SphericalVoronoi3DGUT<num_sv>::Screen v_proj, SphericalVoronoi3DGUT<num_sv>::Screen h_proj,
    SphericalVoronoi3DGUT<num_sv>::World& v_world, float3x3 &v_R, float3 &v_t,
    float3 &h_world_pos
) {}  // TODO

template<int num_sv>
inline __device__ void SphericalVoronoi3DGUT<num_sv>::project_fisheye_vjp(
    SphericalVoronoi3DGUT<num_sv>::World world, SphericalVoronoi3DGUT<num_sv>::BwdProjCamera cam,
    SphericalVoronoi3DGUT<num_sv>::Screen v_proj,
    SphericalVoronoi3DGUT<num_sv>::World& v_world, float3x3 &v_R, float3 &v_t
) {
    projection_3dgut_sv_fisheye_vjp(
        false,
        world.mean, world.quat, world.scale, world.opacity, &world.sv_sites, &world.sv_colors,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, &cam.dist_coeffs,
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
    SphericalVoronoi3DGUT<num_sv>::Screen v_proj, SphericalVoronoi3DGUT<num_sv>::Screen h_proj,
    SphericalVoronoi3DGUT<num_sv>::World& v_world, float3x3 &v_R, float3 &v_t,
    float3 &h_world_pos
) {}  // TODO

#endif  // #ifdef __CUDACC__

typedef SphericalVoronoi3DGUT<2> SphericalVoronoi3DGUT_Default;

