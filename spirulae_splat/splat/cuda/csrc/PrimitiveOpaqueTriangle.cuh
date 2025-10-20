#pragma once

#ifdef __CUDACC__
#define TensorView _Slang_TensorView
#include "generated/projection.cu"
#undef TensorView
#endif

#include "common.cuh"

#include <tuple>



struct OpaqueTriangle {
    struct World;
    struct Screen;

#ifdef __CUDACC__

    struct FwdProjCamera {
        float3x3 R;
        float3 t;
        float fx, fy, cx, cy;
        uint width, height, antialiased;
        float near_plane, far_plane;
    };

    inline static __device__ void project_persp(
        World world, FwdProjCamera cam,
        Screen& screen, int4& aabb, float& depth, float3& normal
    );

    struct BwdProjCamera {
        float3x3 R;
        float3 t;
        float fx, fy, cx, cy;
        uint width, height, antialiased;
    };

    inline static __device__ void project_persp_vjp(
        World world, BwdProjCamera cam,
        Screen v_screen, float v_depth, float3 v_normal,
        World& v_world, float3x3 &v_R, float3 &v_t
    );

#endif  // #ifdef __CUDACC__

};

struct OpaqueTriangle::World {

    typedef std::tuple<at::Tensor, at::Tensor, at::Tensor, at::Tensor> TensorTuple;
    // typedef std::tuple<at::Tensor, at::Tensor> TensorTuple;

    struct Buffer;

    struct Tensor {
        at::Tensor means;
        at::Tensor quats;
        at::Tensor scales;
        // at::Tensor verts;
        at::Tensor hardness;

        Tensor(const TensorTuple& splats) {
            means = std::get<0>(splats);
            quats = std::get<1>(splats);
            scales = std::get<2>(splats);
            hardness = std::get<3>(splats);
            // verts = std::get<0>(splats);
            // hardness = std::get<1>(splats);
        }

        TensorTuple tuple() const {
            return std::make_tuple(means, quats, scales, hardness);
            // return std::make_tuple(verts, hardness);
        }

        Tensor zeros_like() const {
            return Tensor(std::make_tuple(
                at::zeros_like(means),
                at::zeros_like(quats),
                at::zeros_like(scales),
                // at::zeros_like(verts),
                at::zeros_like(hardness)
            ));
        }

        auto options() {
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

        Buffer(const Tensor& tensors) {
            DEVICE_GUARD(tensors.means);
            CHECK_INPUT(tensors.means);
            CHECK_INPUT(tensors.quats);
            CHECK_INPUT(tensors.scales);
            // DEVICE_GUARD(tensors.verts);
            // CHECK_INPUT(tensors.verts);
            CHECK_INPUT(tensors.hardness);
            means = (float3*)tensors.means.data_ptr<float>();
            quats = (float4*)tensors.quats.data_ptr<float>();
            scales = (float3*)tensors.scales.data_ptr<float>();
            // verts = (float3*)tensors.verts.data_ptr<float>();
            hardness = (float2*)tensors.hardness.data_ptr<float>();
        }
    };

    float3 mean;
    float4 quat;
    float3 scale;
    // float3 vert0, vert1, vert2;
    float2 hardness;

#ifdef __CUDACC__

    static __device__ World load(const Buffer &buffer, long idx) {
        return {
            buffer.means[idx],
            buffer.quats[idx],
            buffer.scales[idx],
            // buffer.verts[3*idx+0],
            // buffer.verts[3*idx+1],
            // buffer.verts[3*idx+2],
            buffer.hardness[idx]
        };
    }

    static __device__ __forceinline__ World zero() {
        return {
            {0.f, 0.f, 0.f},
            {0.f, 0.f, 0.f, 0.f},
            // {0.f, 0.f, 0.f},
            {0.f, 0.f, 0.f},
            0.f
        };
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
    }
    
    __device__ void atomicAddBuffer(Buffer &buffer, long idx) {
        atomicAddFVec(buffer.means + idx, mean);
        atomicAddFVec(buffer.quats + idx, quat);
        atomicAddFVec(buffer.scales + idx, scale);
        // atomicAddFVec(buffer.verts + 3*idx+0, vert0);
        // atomicAddFVec(buffer.verts + 3*idx+1, vert1);
        // atomicAddFVec(buffer.verts + 3*idx+2, vert2);
        atomicAddFVec(buffer.hardness + idx, hardness);
    }

#endif  // #ifdef __CUDACC__
};

struct OpaqueTriangle::Screen {

    typedef std::tuple<at::Tensor, at::Tensor> TensorTuple;

    struct Buffer;

    struct Tensor {
        at::Tensor vertices;  // [..., 3, 2]
        at::Tensor hardness;  // [..., 2]
        std::optional<at::Tensor> absgrad;

        Tensor(const TensorTuple& splats) {
            vertices = std::get<0>(splats);
            hardness = std::get<1>(splats);
        }

        TensorTuple tuple() const {
            return std::make_tuple(vertices, hardness);
        }

        Tensor zeros_like(bool absgrad) const {
            Tensor result = Tensor(std::make_tuple(
                at::zeros_like(vertices),
                at::zeros_like(hardness)
            ));
            if (absgrad) {
                // TODO: batched tensor
                result.absgrad = at::zeros({size(), 2}, vertices.options());
            }
            return result;
        }

        static Tensor empty(long C, long N, c10::TensorOptions opt) {
            return std::make_tuple(
                at::empty({C, N, 3, 2}, opt),
                at::empty({C, N, 2}, opt)
            );
        }

        auto options() {
            return vertices.options();
        }
        bool isPacked() const {
            return vertices.dim() == 3;
        }
        long size() const {
            return vertices.size(-3);
        }

        Buffer buffer() { return Buffer(*this); }
    };

    struct Buffer {
        float2* __restrict__ vertices;  // [I, N, 3, 2] or [nnz, 3, 2]
        float2* __restrict__ hardness;  // [I, N, 2] or [nnz, 2]
        float2* __restrict__ absgrad;  // [I, N, 2] or [nnz, 2]

        Buffer(const Tensor& tensors) {
            DEVICE_GUARD(tensors.vertices);
            CHECK_INPUT(tensors.vertices);
            CHECK_INPUT(tensors.hardness);
            vertices = (float2*)tensors.vertices.data_ptr<float>();
            hardness = (float2*)tensors.hardness.data_ptr<float>();
            absgrad = tensors.absgrad.has_value() ?
                (float2*)tensors.absgrad.value().data_ptr<float>()
                : nullptr;
        }
    };

    float2 vert0;
    float2 vert1;
    float2 vert2;
    float2 hardness;
    float2 absgrad;

#ifdef __CUDACC__

    static __device__ Screen load(const Buffer &buffer, long idx) {
        return {
            buffer.vertices[3*idx+0],
            buffer.vertices[3*idx+1],
            buffer.vertices[3*idx+2],
            buffer.hardness[idx]
            // absgrad is undefined
        };
    }

    static __device__ __forceinline__ Screen zero() {
        return {
            {0.f, 0.f},
            {0.f, 0.f},
            {0.f, 0.f},
            {0.f, 0.f},
            {0.f, 0.f}
        };
    }

    __device__ __forceinline__ void operator+=(const Screen &other) {
        vert0 += other.vert0;
        vert1 += other.vert1;
        vert2 += other.vert2;
        hardness += other.hardness;
        absgrad += fabs(other.vert0 + other.vert1 + other.vert2) / 3.0f;
    }

    __device__ void saveBuffer(Buffer &buffer, long idx) {
        buffer.vertices[3*idx+0] = vert0;
        buffer.vertices[3*idx+1] = vert2;
        buffer.vertices[3*idx+2] = vert1;
        buffer.hardness[idx] = hardness;
        if (buffer.absgrad != nullptr)
            buffer.absgrad[idx] = absgrad;
    }

    __device__ void atomicAddBuffer(Buffer &buffer, long idx) {
        atomicAddFVec(buffer.vertices + 3*idx+0, vert0);
        atomicAddFVec(buffer.vertices + 3*idx+1, vert1);
        atomicAddFVec(buffer.vertices + 3*idx+2, vert2);
        atomicAddFVec(buffer.hardness + idx, hardness);
        if (buffer.absgrad != nullptr)
            atomicAddFVec(buffer.absgrad + idx, absgrad);
    }

    __device__ __forceinline__ float evaluate_alpha(float px, float py) {
        // return evaluate_alpha_opaque_triangle_fast(
        return evaluate_alpha_opaque_triangle_precise(
            vert0, vert1, vert2, hardness,
            make_float2(px, py)
        );
    }

    __device__ __forceinline__ Screen evaluate_alpha_vjp(float px, float py, float v_alpha) {
        Screen v_splat = Screen::zero();
        // evaluate_alpha_opaque_triangle_fast_vjp(
        evaluate_alpha_opaque_triangle_precise_vjp(
            vert0, vert1, vert2, hardness,
            make_float2(px, py), v_alpha,
            &v_splat.vert0, &v_splat.vert1, &v_splat.vert2, &v_splat.hardness
        );
        return v_splat;
    }

#endif  // #ifdef __CUDACC__

};


#ifdef __CUDACC__

inline __device__ void OpaqueTriangle::project_persp(
    OpaqueTriangle::World world, OpaqueTriangle::FwdProjCamera cam,
    OpaqueTriangle::Screen& screen, int4& aabb, float& depth, float3& normal
) {
    projection_opaque_triangle_persp(
        world.mean, world.quat, world.scale, world.hardness,
        // world.vert0, world.vert1, world.vert2, world.hardness,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy,
        cam.width, cam.height, cam.near_plane, cam.far_plane,
        &aabb, &depth, &normal, &screen.vert0, &screen.vert1, &screen.vert2, &screen.hardness
    );
}

inline __device__ void OpaqueTriangle::project_persp_vjp(
    OpaqueTriangle::World world, OpaqueTriangle::BwdProjCamera cam,
    OpaqueTriangle::Screen v_screen, float v_depth, float3 v_normal,
    OpaqueTriangle::World& v_world, float3x3 &v_R, float3 &v_t
) {
    projection_opaque_triangle_persp_vjp(
        world.mean, world.quat, world.scale, world.hardness,
        // world.vert0, world.vert1, world.vert2, world.hardness,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy,
        cam.width, cam.height,
        v_depth, v_normal, v_screen.vert0, v_screen.vert1, v_screen.vert2, v_screen.hardness,
        &v_world.mean, &v_world.quat, &v_world.scale, &v_world.hardness,
        // &v_world.vert0, &v_world.vert1, &v_world.vert2, &v_world.hardness,
        &v_R, &v_t
    );
}

#endif  // #ifdef __CUDACC__
