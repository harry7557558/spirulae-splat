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

// template<bool antialiased>
// struct _Base3DGS  {
//     struct World;
//     struct Screen;
//     static constexpr RenderOutputType pixelType = RenderOutputType::RGB_D;

// #ifdef __CUDACC__

// #endif  // #ifdef __CUDACC__

// };


template<bool antialiased>
struct _Base3DGS : public _BasePrimitive3DGS {
    static constexpr RenderOutputType pixelType = RenderOutputType::RGB_D;
};

typedef _Base3DGS<false> Vanilla3DGS;
typedef _Base3DGS<true> MipSplatting;


#if 0
struct _Base3DGS<antialiased>::World {

    float3 mean;
    float4 quat;
    float3 scale;
    float opacity;
    FixedArray<float3, 16> sh_coeffs;

    #ifndef NO_TORCH
    typedef std::tuple<at::Tensor, at::Tensor, at::Tensor, at::Tensor, at::Tensor, std::optional<at::Tensor>> TensorTuple;
    #endif

    struct Buffer;

    #ifndef NO_TORCH
    struct Tensor {
        at::Tensor means;
        at::Tensor quats;
        at::Tensor scales;
        at::Tensor opacities;
        at::Tensor features_dc;
        std::optional<at::Tensor> features_sh;

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

        Tensor allocProjBwd(bool is_hess_diag) const {
            return Tensor(std::make_tuple(
                zeros_like_tensor(means),
                zeros_like_tensor(quats),
                zeros_like_tensor(scales),
                zeros_like_tensor(opacities),
                zeros_like_tensor(features_dc),
                features_sh.has_value() && !is_hess_diag ?
                    (std::optional<at::Tensor>)zeros_like_tensor(features_sh.value()) :
                    (std::optional<at::Tensor>)std::nullopt
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

    struct Buffer : public TensorArray<5> {
        float3* __restrict__ features_sh;
        uint num_sh;

        using TensorArray<5>::TensorArray;

        Buffer() : features_sh(nullptr), num_sh(0) {}

        #ifndef NO_TORCH
        Buffer(const Tensor& tensors)
            : TensorArray<5>(std::vector<std::optional<at::Tensor>>{
                tensors.means, tensors.quats, tensors.scales, tensors.opacities, tensors.features_dc
            }) {
            features_sh = tensors.features_sh.has_value() ?
                (float3*)tensors.features_sh.value().template data_ptr<float>() : nullptr;
            num_sh = tensors.features_sh.has_value() ?
                tensors.features_sh.value().size(-2) : 0;
        }
        #endif
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
        World world = {
            {0.f, 0.f, 0.f},
            {0.f, 0.f, 0.f, 0.f},
            {0.f, 0.f, 0.f},
            0.f,
        };
        for (int i = 0; i < 16; i++)
            world.sh_coeffs[i] = make_float3(0);
        return world;
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

template<bool antialiased>
struct _Base3DGS<antialiased>::Screen {

    float2 xy;
    float depth;
    float3 conic;
    float opac;
    float3 rgb;

    #ifndef NO_TORCH
    typedef std::tuple<at::Tensor, at::Tensor, at::Tensor, at::Tensor, at::Tensor> TensorTupleProj;
    typedef std::tuple<at::Tensor, at::Tensor, at::Tensor, at::Tensor, at::Tensor> TensorTuple;
    #endif

    struct Buffer;

    #ifndef NO_TORCH
    struct Tensor {
        at::Tensor means2d;
        at::Tensor depths;
        at::Tensor conics;
        at::Tensor opacities;
        at::Tensor rgbs;

        Tensor() {}

        Tensor(const TensorTuple& splats) {
            means2d = std::get<0>(splats);
            depths = std::get<1>(splats);
            conics = std::get<2>(splats);
            opacities = std::get<3>(splats);
            rgbs = std::get<4>(splats);
        }

        TensorTupleProj tupleProjFwd() const {
            return std::make_tuple(means2d, depths, conics, opacities, rgbs);
        }

        TensorTuple tupleProjFwdPacked() const {
            return std::make_tuple(
                means2d, depths, conics, opacities, rgbs
            );
        }

        TensorTuple tupleRasterBwd() const {
            return std::make_tuple(means2d, depths, conics, opacities, rgbs);
        }

        static TensorTupleProj allocProjFwd(long C, long N, c10::TensorOptions opt) {
            return std::make_tuple(
                at::empty({C, N, 2}, opt),  // means2d
                at::empty({C, N}, opt),  // depths
                at::empty({C, N, 3}, opt),  // conics
                at::empty({C, N}, opt),  // opacities
                at::empty({C, N, 3}, opt)  // rgbs
            );
        }

        static TensorTuple allocProjFwdPacked(long nnz, c10::TensorOptions opt) {
            return std::make_tuple(
                at::empty({nnz, 2}, opt),
                at::empty({nnz}, opt),
                at::empty({nnz, 3}, opt),
                at::empty({nnz}, opt),
                at::empty({nnz, 3}, opt)
            );
        }

        Tensor allocRasterBwd() const {
            return std::make_tuple(
                zeros_like_tensor(means2d),
                zeros_like_tensor(depths),
                zeros_like_tensor(conics),
                zeros_like_tensor(opacities),
                zeros_like_tensor(rgbs)
            );
        }

        auto options() const {
            return means2d.options();
        }
        bool isPacked() const {
            return means2d.dim() == 2;
        }
        long size() const {
            return means2d.size(-2);
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
                tensors.means2d, tensors.depths, tensors.conics, tensors.opacities, tensors.rgbs
            }) {
            size = tensors.opacities.numel();
        }
        #endif
    };

#ifdef __CUDACC__

    static __device__ Screen load(const Buffer &buffer, long idx, const uint32_t* gaussian_ids) {
        return {
            buffer.means2d[idx],
            buffer.depths[idx],
            buffer.conics[idx],
            buffer.opacities[idx],
            buffer.rgbs[idx],
        };
    }

    static __device__ Screen loadWithPrecompute(const Buffer &buffer, long idx, const uint32_t* gaussian_ids) {
        return load(buffer, idx, gaussian_ids);
    }

    static __device__ __forceinline__ Screen zero() {
        return {
            {0.f, 0.f},
            0.0f,
            {0.f, 0.f, 0.f},
            0.0f,
            {0.f, 0.f, 0.f},
        };
    }

    __device__ __forceinline__ void precomputeBackward(Screen& grad) const {
    }

    __device__ __forceinline__ void addGradient(const Screen &grad, float weight=1.0f) {
        xy += grad.xy * weight;
        depth += grad.depth * weight;
        conic += grad.conic * weight;
        opac += grad.opac * weight;
        rgb += grad.rgb * weight;
    }

    __device__ __forceinline__ void addGaussNewtonHessianDiagonal(const Screen &grad, float weight=1.0f) {
        xy += fmul_axa(grad.xy, weight);
        depth += fmul_axa(grad.depth, weight);
        conic += fmul_axa(grad.conic, weight);
        opac += fmul_axa(grad.opac, weight);
        rgb += fmul_axa(grad.rgb, weight);
    }

    __device__ void saveParamsToBuffer(Buffer &buffer, long idx, const uint32_t* gaussian_ids) {
        buffer.means2d[idx] = xy;
        buffer.depths[idx] = depth;
        buffer.conics[idx] = conic;
        buffer.opacities[idx] = opac;
        buffer.rgbs[idx] = rgb;
    }

    __device__ void atomicAddToBuffer(Buffer &buffer, long idx, const uint32_t* gaussian_ids) const {
        atomicAddFVec(buffer.means2d + idx, xy);
        atomicAddFVec(buffer.depths + idx, depth);
        atomicAddFVec(buffer.conics + idx, conic);
        atomicAddFVec(buffer.opacities + idx, opac);
        atomicAddFVec(buffer.rgbs + idx, rgb);
    }

    __device__ __forceinline__ float evaluate_alpha(float px, float py) {
        float2 delta = {xy.x - px, xy.y - py};
        float sigma = 0.5f * (conic.x * delta.x * delta.x +
                                conic.z * delta.y * delta.y) +
                        conic.y * delta.x * delta.y;
        float vis = __expf(-sigma);
        float alpha = min(0.999f, opac * vis);
        return sigma < 0.f ? 0.f : alpha;
    }

    __device__ __forceinline__ Screen evaluate_alpha_vjp(float px, float py, float v_alpha) {
        float2 delta = {xy.x - px, xy.y - py};
        float sigma = 0.5f * (conic.x * delta.x * delta.x +
                                conic.z * delta.y * delta.y) +
                        conic.y * delta.x * delta.y;
        float vis = __expf(-sigma);
        float alpha = min(0.999f, opac * vis);

        Screen v_splat = Screen::zero();
        if (sigma >= 0.f && opac * vis <= 0.999f) {
            const float v_sigma = -opac * vis * v_alpha;
            v_splat.conic = {
                0.5f * v_sigma * delta.x * delta.x,
                v_sigma * delta.x * delta.y,
                0.5f * v_sigma * delta.y * delta.y
            };
            v_splat.xy = {
                v_sigma * (conic.x * delta.x + conic.y * delta.y),
                v_sigma * (conic.y * delta.x + conic.z * delta.y)
            };
            v_splat.opac = vis * v_alpha;
        }
        return v_splat;
    }

    __device__ __forceinline__ RenderOutput evaluate_color(float px, float py) {
        return { rgb, depth };
    }

    __device__ __forceinline__ Screen evaluate_color_vjp(float px, float py, RenderOutput v_render) {
        Screen v_splat = Screen::zero();
        v_splat.rgb = v_render.rgb;
        v_splat.depth = v_render.depth;
        return v_splat;
    }

#endif  // #ifdef __CUDACC__

};
#endif


#ifdef __CUDACC__

template<bool antialiased>
inline __device__ void project_persp(
    typename _Base3DGS<antialiased>::World world, ProjCamera cam,
    typename _Base3DGS<antialiased>::Screen& screen, float4& aabb
) {
    Slang3DGS::projection_3dgs_persp(
        antialiased,
        world.mean, world.quat, world.scale, world.opacity, world.sh_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy,
        cam.dist_coeffs,
        cam.width, cam.height,
        &aabb, &screen.xy, &screen.depth, &screen.conic, &screen.opac, &screen.rgb
    );
}

template<bool antialiased>
inline __device__ void project_ortho(
    typename _Base3DGS<antialiased>::World world, ProjCamera cam,
    typename _Base3DGS<antialiased>::Screen& screen, float4& aabb
) {
    Slang3DGS::projection_3dgs_ortho(
        antialiased,
        world.mean, world.quat, world.scale, world.opacity, world.sh_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
        cam.width, cam.height,
        &aabb, &screen.xy, &screen.depth, &screen.conic, &screen.opac, &screen.rgb
    );
}

template<bool antialiased>
inline __device__ void project_fisheye(
    typename _Base3DGS<antialiased>::World world, ProjCamera cam,
    typename _Base3DGS<antialiased>::Screen& screen, float4& aabb
) {
    Slang3DGS::projection_3dgs_fisheye(
        antialiased,
        world.mean, world.quat, world.scale, world.opacity, world.sh_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
        cam.width, cam.height,
        &aabb, &screen.xy, &screen.depth, &screen.conic, &screen.opac, &screen.rgb
    );
}


template<bool antialiased>
inline __device__ void project_persp_vjp(
    typename _Base3DGS<antialiased>::World world, ProjCamera cam,
    typename _Base3DGS<antialiased>::Screen v_screen,
    typename _Base3DGS<antialiased>::World& v_world, float3x3 &v_R, float3 &v_t
) {
    Slang3DGS::projection_3dgs_persp_vjp(
        antialiased,
        world.mean, world.quat, world.scale, world.opacity, world.sh_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
        cam.width, cam.height,
        v_screen.xy, v_screen.depth, v_screen.conic, v_screen.opac, v_screen.rgb,
        &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity, &v_world.sh_coeffs,
        &v_R, &v_t
    );
}

template<bool antialiased>
inline __device__ void project_persp_vjp(
    typename _Base3DGS<antialiased>::World world, ProjCamera cam,
    typename _Base3DGS<antialiased>::Screen v_proj, typename _Base3DGS<antialiased>::Screen vr_proj, typename _Base3DGS<antialiased>::Screen h_proj,
    typename _Base3DGS<antialiased>::World& v_world, float3x3 &v_R, float3 &v_t,
    float3 &vr_world_pos, float3 &h_world_pos
) {
    Slang3DGS::projection_3dgs_persp_vjp(
        antialiased,
        world.mean, world.quat, world.scale, world.opacity, world.sh_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
        cam.width, cam.height,
        v_proj.xy, v_proj.depth, v_proj.conic, v_proj.opac, v_proj.rgb,
        vr_proj.xy, vr_proj.depth, vr_proj.conic, vr_proj.opac, vr_proj.rgb,
        h_proj.xy, h_proj.depth, h_proj.conic, h_proj.opac, h_proj.rgb,
        &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity, &v_world.sh_coeffs,
        &v_R, &v_t, &vr_world_pos, &h_world_pos
    );
}

template<bool antialiased>
inline __device__ void project_persp_vjp(
    typename _Base3DGS<antialiased>::World world, ProjCamera cam,
    typename _Base3DGS<antialiased>::Screen v_proj, typename _Base3DGS<antialiased>::Screen vr_proj, typename _Base3DGS<antialiased>::Screen h_proj,
    typename _Base3DGS<antialiased>::World& v_world, float3x3 &v_R, float3 &v_t,
    typename _Base3DGS<antialiased>::World& vr_world, typename _Base3DGS<antialiased>::World& h_world
) {
    Slang3DGS::projection_3dgs_persp_vjp(
        antialiased,
        world.mean, world.quat, world.scale, world.opacity, world.sh_coeffs,
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

template<bool antialiased>
inline __device__ void project_fisheye_vjp(
    typename _Base3DGS<antialiased>::World world, ProjCamera cam,
    typename _Base3DGS<antialiased>::Screen v_screen,
    typename _Base3DGS<antialiased>::World& v_world, float3x3 &v_R, float3 &v_t
) {
    Slang3DGS::projection_3dgs_fisheye_vjp(
        antialiased,
        world.mean, world.quat, world.scale, world.opacity, world.sh_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
        cam.width, cam.height,
        v_screen.xy, v_screen.depth, v_screen.conic, v_screen.opac, v_screen.rgb,
        &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity, &v_world.sh_coeffs,
        &v_R, &v_t
    );
}

template<bool antialiased>
inline __device__ void project_fisheye_vjp(
    typename _Base3DGS<antialiased>::World world, ProjCamera cam,
    typename _Base3DGS<antialiased>::Screen v_proj, typename _Base3DGS<antialiased>::Screen vr_proj, typename _Base3DGS<antialiased>::Screen h_proj,
    typename _Base3DGS<antialiased>::World& v_world, float3x3 &v_R, float3 &v_t,
    float3 &vr_world_pos, float3 &h_world_pos
) {
    Slang3DGS::projection_3dgs_fisheye_vjp(
        antialiased,
        world.mean, world.quat, world.scale, world.opacity, world.sh_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy, cam.dist_coeffs,
        cam.width, cam.height,
        v_proj.xy, v_proj.depth, v_proj.conic, v_proj.opac, v_proj.rgb,
        vr_proj.xy, vr_proj.depth, vr_proj.conic, vr_proj.opac, vr_proj.rgb,
        h_proj.xy, h_proj.depth, h_proj.conic, h_proj.opac, h_proj.rgb,
        &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity, &v_world.sh_coeffs,
        &v_R, &v_t, &vr_world_pos, &h_world_pos
    );
}

template<bool antialiased>
inline __device__ void project_fisheye_vjp(
    typename _Base3DGS<antialiased>::World world, ProjCamera cam,
    typename _Base3DGS<antialiased>::Screen v_proj, typename _Base3DGS<antialiased>::Screen vr_proj, typename _Base3DGS<antialiased>::Screen h_proj,
    typename _Base3DGS<antialiased>::World& v_world, float3x3 &v_R, float3 &v_t,
    typename _Base3DGS<antialiased>::World& vr_world, typename _Base3DGS<antialiased>::World& h_world
) {
    Slang3DGS::projection_3dgs_fisheye_vjp(
        antialiased,
        world.mean, world.quat, world.scale, world.opacity, world.sh_coeffs,
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

#endif  // #ifdef __CUDACC__

