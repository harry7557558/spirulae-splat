#pragma once

#ifdef __CUDACC__
#define TensorView _Slang_TensorView
#include "generated/projection.cu"
#undef TensorView
#endif

#include "common.cuh"

#include <tuple>


struct Vanilla3DGS {
    struct World;
    struct Screen;

#ifdef __CUDACC__

    struct FwdProjCamera {
        float3x3 R;
        float3 t;
        float fx, fy, cx, cy;
        uint width, height, antialiased;
        float near_plane, far_plane;
        float4 radial_coeffs = {0, 0, 0, 0};
        float2 tangential_coeffs = {0, 0};
        float2 thin_prism_coeffs = {0, 0};
    };

    inline static __device__ void project_persp(
        World world, FwdProjCamera cam,
        Screen& screen, int4& aabb, float& depth, float3& normal
    );

    inline static __device__ void project_ortho(
        World world, FwdProjCamera cam,
        Screen& screen, int4& aabb, float& depth, float3& normal
    );

    inline static __device__ void project_fisheye(
        World world, FwdProjCamera cam,
        Screen& screen, int4& aabb, float& depth, float3& normal
    );

    struct BwdProjCamera {
        float3x3 R;
        float3 t;
        float fx, fy, cx, cy;
        uint width, height, antialiased;
        float4 radial_coeffs = {0, 0, 0, 0};
        float2 tangential_coeffs = {0, 0};
        float2 thin_prism_coeffs = {0, 0};
    };

    inline static __device__ void project_persp_vjp(
        World world, BwdProjCamera cam,
        Screen v_screen, float v_depth, float3 v_normal,
        World& v_world, float3x3 &v_R, float3 &v_t
    );

    inline static __device__ void project_ortho_vjp(
        World world, BwdProjCamera cam,
        Screen v_screen, float v_depth, float3 v_normal,
        World& v_world, float3x3 &v_R, float3 &v_t
    );

    inline static __device__ void project_fisheye_vjp(
        World world, BwdProjCamera cam,
        Screen v_screen, float v_depth, float3 v_normal,
        World& v_world, float3x3 &v_R, float3 &v_t
    );

#endif  // #ifdef __CUDACC__

    static constexpr float eps2d = 0.3f;
};

struct Vanilla3DGS::World {

    typedef std::tuple<at::Tensor, at::Tensor, at::Tensor, at::Tensor> TensorTuple;

    struct Buffer;

    struct Tensor {
        at::Tensor means;
        at::Tensor quats;
        at::Tensor scales;
        at::Tensor opacities;

        Tensor(const TensorTuple& splats) {
            means = std::get<0>(splats);
            quats = std::get<1>(splats);
            scales = std::get<2>(splats);
            opacities = std::get<3>(splats);
        }

        TensorTuple tuple() const {
            return std::make_tuple(means, quats, scales, opacities);
        }

        Tensor zeros_like() const {
            return Tensor(std::make_tuple(
                at::zeros_like(means),
                at::zeros_like(quats),
                at::zeros_like(scales),
                at::zeros_like(opacities)
            ));
        }

        auto options() {
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

        Buffer(const Tensor& tensors) {
            DEVICE_GUARD(tensors.means);
            CHECK_INPUT(tensors.means);
            CHECK_INPUT(tensors.quats);
            CHECK_INPUT(tensors.scales);
            CHECK_INPUT(tensors.opacities);
            means = (float3*)tensors.means.data_ptr<float>();
            quats = (float4*)tensors.quats.data_ptr<float>();
            scales = (float3*)tensors.scales.data_ptr<float>();
            opacities = tensors.opacities.data_ptr<float>();
        }
    };

    float3 mean;
    float4 quat;
    float3 scale;
    float opacity;

#ifdef __CUDACC__

    static __device__ World load(const Buffer &buffer, long idx) {
        return {
            buffer.means[idx],
            buffer.quats[idx],
            buffer.scales[idx],
            buffer.opacities[idx]
        };
    }

    static __device__ __forceinline__ World zero() {
        return {
            {0.f, 0.f, 0.f},
            {0.f, 0.f, 0.f, 0.f},
            {0.f, 0.f, 0.f},
            0.f
        };
    }

    template<typename Partition>
    __device__ void reduce(Partition& partition) {
        warpSum(mean, partition);
        warpSum(quat, partition);
        warpSum(scale, partition);
        warpSum(opacity, partition);
    }
    
    __device__ void atomicAddBuffer(Buffer &buffer, long idx) {
        atomicAddFVec(buffer.means + idx, mean);
        atomicAddFVec(buffer.quats + idx, quat);
        atomicAddFVec(buffer.scales + idx, scale);
        atomicAdd(buffer.opacities + idx, opacity);
    }

#endif  // #ifdef __CUDACC__
};

struct Vanilla3DGS::Screen {

    typedef std::tuple<at::Tensor, at::Tensor, at::Tensor> TensorTuple;

    struct Buffer;

    struct Tensor {
        at::Tensor means2d;
        at::Tensor conics;
        at::Tensor opacities;
        std::optional<at::Tensor> absgrad;

        Tensor(const TensorTuple& splats) {
            means2d = std::get<0>(splats);
            conics = std::get<1>(splats);
            opacities = std::get<2>(splats);
        }

        TensorTuple tuple() const {
            return std::make_tuple(means2d, conics, opacities);
        }

        Tensor zeros_like(bool absgrad) const {
            Tensor result = Tensor(std::make_tuple(
                at::zeros_like(means2d),
                at::zeros_like(conics),
                at::zeros_like(opacities)
            ));
            if (absgrad)
                result.absgrad = at::zeros_like(means2d);
            return result;
        }

        static Tensor empty(long C, long N, c10::TensorOptions opt) {
            return std::make_tuple(
                at::empty({C, N, 2}, opt),
                at::empty({C, N, 3}, opt),
                at::empty({C, N}, opt)
            );
        }

        auto options() {
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

    struct Buffer {
        float2* __restrict__ means2d;  // [I, N, 2] or [nnz, 2]
        float3* __restrict__ conics;  // [I, N, 3] or [nnz, 3]
        float* __restrict__ opacities;  // [I, N] or [nnz]
        float2* __restrict__ absgrad;

        Buffer(const Tensor& tensors) {
            DEVICE_GUARD(tensors.means2d);
            CHECK_INPUT(tensors.means2d);
            CHECK_INPUT(tensors.conics);
            CHECK_INPUT(tensors.opacities);
            means2d = (float2*)tensors.means2d.data_ptr<float>();
            conics = (float3*)tensors.conics.data_ptr<float>();
            opacities = tensors.opacities.data_ptr<float>();
            absgrad = tensors.absgrad.has_value() ?
                (float2*)tensors.absgrad.value().data_ptr<float>()
                : nullptr;
        }
    };

    float2 xy;
    float3 conic;
    float opac;
    float2 xy_abs;

#ifdef __CUDACC__

    static __device__ Screen load(const Buffer &buffer, long idx) {
        return {
            buffer.means2d[idx],
            buffer.conics[idx],
            buffer.opacities[idx],
            // xy_abs is undefined
        };
    }

    static __device__ __forceinline__ Screen zero() {
        return {
            {0.f, 0.f},
            {0.f, 0.f, 0.f},
            0.0f,
            {0.f, 0.f}
        };
    }

    __device__ __forceinline__ void operator+=(const Screen &other) {
        xy += other.xy;
        conic += other.conic;
        opac += other.opac;
        xy_abs += fabs(other.xy);
    }

    __device__ void saveBuffer(Buffer &buffer, long idx) {
        buffer.means2d[idx] = xy;
        buffer.conics[idx] = conic;
        buffer.opacities[idx] = opac;
        if (buffer.absgrad != nullptr)
            buffer.absgrad[idx] = xy_abs;
    }

    __device__ void atomicAddBuffer(Buffer &buffer, long idx) {
        atomicAddFVec(buffer.means2d + idx, xy);
        atomicAddFVec(buffer.conics + idx, conic);
        atomicAdd(buffer.opacities + idx, opac);
        if (buffer.absgrad != nullptr)
            atomicAddFVec(buffer.absgrad + idx, xy_abs);
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

#endif  // #ifdef __CUDACC__

};


#ifdef __CUDACC__

inline __device__ void Vanilla3DGS::project_persp(
    Vanilla3DGS::World world, Vanilla3DGS::FwdProjCamera cam,
    Vanilla3DGS::Screen& screen, int4& aabb, float& depth, float3& normal
) {
    projection_3dgs_persp(
        bool(cam.antialiased),
        world.mean, world.quat, world.scale, world.opacity,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy,
        cam.radial_coeffs, cam.tangential_coeffs, cam.thin_prism_coeffs,
        cam.width, cam.height, cam.near_plane, cam.far_plane,
        &aabb, &depth, &screen.xy, &screen.conic, &screen.opac
    );
}

inline __device__ void Vanilla3DGS::project_ortho(
    Vanilla3DGS::World world, Vanilla3DGS::FwdProjCamera cam,
    Vanilla3DGS::Screen& screen, int4& aabb, float& depth, float3& normal
) {
    projection_3dgs_ortho(
        bool(cam.antialiased),
        world.mean, world.quat, world.scale, world.opacity,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy,
        cam.radial_coeffs, cam.tangential_coeffs, cam.thin_prism_coeffs,
        cam.width, cam.height, cam.near_plane, cam.far_plane,
        &aabb, &depth, &screen.xy, &screen.conic, &screen.opac
    );
}

inline __device__ void Vanilla3DGS::project_fisheye(
    Vanilla3DGS::World world, Vanilla3DGS::FwdProjCamera cam,
    Vanilla3DGS::Screen& screen, int4& aabb, float& depth, float3& normal
) {
    projection_3dgs_fisheye(
        bool(cam.antialiased),
        world.mean, world.quat, world.scale, world.opacity,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy,
        cam.radial_coeffs, cam.tangential_coeffs, cam.thin_prism_coeffs,
        cam.width, cam.height, cam.near_plane, cam.far_plane,
        &aabb, &depth, &screen.xy, &screen.conic, &screen.opac
    );
}

inline __device__ void Vanilla3DGS::project_persp_vjp(
    Vanilla3DGS::World world, Vanilla3DGS::BwdProjCamera cam,
    Vanilla3DGS::Screen v_screen, float v_depth, float3 v_normal,
    Vanilla3DGS::World& v_world, float3x3 &v_R, float3 &v_t
) {
    projection_3dgs_persp_vjp(
        bool(cam.antialiased),
        world.mean, world.quat, world.scale, world.opacity,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy,
        cam.radial_coeffs, cam.tangential_coeffs, cam.thin_prism_coeffs,
        cam.width, cam.height,
        v_depth, v_screen.xy, v_screen.conic, v_screen.opac,
        &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity,
        &v_R, &v_t
    );
}

inline __device__ void Vanilla3DGS::project_ortho_vjp(
    Vanilla3DGS::World world, Vanilla3DGS::BwdProjCamera cam,
    Vanilla3DGS::Screen v_screen, float v_depth, float3 v_normal,
    Vanilla3DGS::World& v_world, float3x3 &v_R, float3 &v_t
) {
    projection_3dgs_ortho_vjp(
        bool(cam.antialiased),
        world.mean, world.quat, world.scale, world.opacity,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy,
        cam.radial_coeffs, cam.tangential_coeffs, cam.thin_prism_coeffs,
        cam.width, cam.height,
        v_depth, v_screen.xy, v_screen.conic, v_screen.opac,
        &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity,
        &v_R, &v_t
    );
}

inline __device__ void Vanilla3DGS::project_fisheye_vjp(
    Vanilla3DGS::World world, Vanilla3DGS::BwdProjCamera cam,
    Vanilla3DGS::Screen v_screen, float v_depth, float3 v_normal,
    Vanilla3DGS::World& v_world, float3x3 &v_R, float3 &v_t
) {
    projection_3dgs_fisheye_vjp(
        bool(cam.antialiased),
        world.mean, world.quat, world.scale, world.opacity,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy,
        cam.radial_coeffs, cam.tangential_coeffs, cam.thin_prism_coeffs,
        cam.width, cam.height,
        v_depth, v_screen.xy, v_screen.conic, v_screen.opac,
        &v_world.mean, &v_world.quat, &v_world.scale, &v_world.opacity,
        &v_R, &v_t
    );
}

#endif  // #ifdef __CUDACC__
