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
    struct RenderOutput;

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
        Screen& screen, int4& aabb
    );

    inline static __device__ void project_fisheye(
        World world, FwdProjCamera cam,
        Screen& screen, int4& aabb
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
        Screen v_screen,
        World& v_world, float3x3 &v_R, float3 &v_t
    );

    inline static __device__ void project_fisheye_vjp(
        World world, BwdProjCamera cam,
        Screen v_screen,
        World& v_world, float3x3 &v_R, float3 &v_t
    );

#endif  // #ifdef __CUDACC__

};

struct OpaqueTriangle::World {

    float3 mean;
    float4 quat;
    float3 scale;
    // float3 vert0, vert1, vert2;
    float2 hardness;
#ifdef __CUDACC__
    FixedArray<float3, 16> sh_coeffs;
#endif

    typedef std::tuple<at::Tensor, at::Tensor, at::Tensor, at::Tensor, at::Tensor, at::Tensor> TensorTuple;
    // typedef std::tuple<at::Tensor, at::Tensor, at::Tensor, at::Tensor> TensorTuple;

    struct Buffer;

    struct Tensor {
        at::Tensor means;
        at::Tensor quats;
        at::Tensor scales;
        // at::Tensor verts;
        at::Tensor hardness;
        at::Tensor features_dc;
        at::Tensor features_sh;

        Tensor(const TensorTuple& splats) {
            means = std::get<0>(splats);
            quats = std::get<1>(splats);
            scales = std::get<2>(splats);
            hardness = std::get<3>(splats);
            // verts = std::get<0>(splats);
            // hardness = std::get<1>(splats);
            features_dc = std::get<4>(splats);
            features_sh = std::get<5>(splats);
        }

        TensorTuple tuple() const {
            return std::make_tuple(means, quats, scales, hardness, features_dc, features_sh);
            // return std::make_tuple(verts, hardness, features_dc, features_sh);
        }

        Tensor zeros_like() const {
            return Tensor(std::make_tuple(
                at::zeros_like(means),
                at::zeros_like(quats),
                at::zeros_like(scales),
                // at::zeros_like(verts),
                at::zeros_like(hardness),
                at::zeros_like(features_dc),
                at::zeros_like(features_sh)
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
        uint num_sh;

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
            means = (float3*)tensors.means.data_ptr<float>();
            quats = (float4*)tensors.quats.data_ptr<float>();
            scales = (float3*)tensors.scales.data_ptr<float>();
            // verts = (float3*)tensors.verts.data_ptr<float>();
            hardness = (float2*)tensors.hardness.data_ptr<float>();
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
            // buffer.verts[3*idx+0],
            // buffer.verts[3*idx+1],
            // buffer.verts[3*idx+2],
            buffer.hardness[idx]
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
            // {0.f, 0.f, 0.f},
            {0.f, 0.f, 0.f},
            0.f
        };
        for (int i = 0; i < 16; i++)
            world.sh_coeffs[i] = make_float3(0);
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
    }
    
    __device__ void atomicAddBuffer(Buffer &buffer, long idx) {
        atomicAddFVec(buffer.means + idx, mean);
        atomicAddFVec(buffer.quats + idx, quat);
        atomicAddFVec(buffer.scales + idx, scale);
        // atomicAddFVec(buffer.verts + 3*idx+0, vert0);
        // atomicAddFVec(buffer.verts + 3*idx+1, vert1);
        // atomicAddFVec(buffer.verts + 3*idx+2, vert2);
        atomicAddFVec(buffer.hardness + idx, hardness);
        atomicAdd(buffer.features_dc + idx, sh_coeffs[0]);
        for (int i = 0; i < buffer.num_sh; i++)
            atomicAddFVec(buffer.features_sh + idx*buffer.num_sh + i, sh_coeffs[i+1]);
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

    __device__ __forceinline__ float dot(const RenderOutput &other) const {
        return (rgb.x * other.rgb.x + rgb.y * other.rgb.y + rgb.z * other.rgb.z)
            + depth * other.depth
            + (normal.x * other.normal.x + normal.y * other.normal.y + normal.z * other.normal.z);
    }

    __device__ void saveBuffer(Buffer &buffer, long idx) {
        buffer.rgbs[idx] = rgb;
        buffer.depths[idx] = depth;
        buffer.normals[idx] = normal;
    }

#endif  // #ifdef __CUDACC__
};


struct OpaqueTriangle::Screen {

    float2 vert0;
    float2 vert1;
    float2 vert2;
    float depth;
    float2 hardness;
    float3 rgb;
    float3 normal;
    float2 absgrad;

    typedef std::tuple<at::Tensor, at::Tensor, at::Tensor, at::Tensor, at::Tensor> TensorTuple;

    struct Buffer;

    struct Tensor {
        at::Tensor vertices;  // [..., 3, 2]
        at::Tensor depths;  // [..., 2]
        at::Tensor hardness;  // [..., 2]
        at::Tensor rgbs;  // [..., 3]
        at::Tensor normals;  // [..., 3]
        std::optional<at::Tensor> absgrad;

        Tensor(const TensorTuple& splats) {
            vertices = std::get<0>(splats);
            depths = std::get<1>(splats);
            hardness = std::get<2>(splats);
            rgbs = std::get<3>(splats);
            normals = std::get<4>(splats);
        }

        TensorTuple tuple() const {
            return std::make_tuple(vertices, depths, hardness, rgbs, normals);
        }

        Tensor zeros_like(bool absgrad) const {
            Tensor result = Tensor(std::make_tuple(
                at::zeros_like(vertices),
                at::zeros_like(depths),
                at::zeros_like(hardness),
                at::zeros_like(rgbs),
                at::zeros_like(normals)
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
                at::empty({C, N}, opt),
                at::empty({C, N, 2}, opt),
                at::empty({C, N, 3}, opt),
                at::empty({C, N, 3}, opt)
            );
        }

        auto options() const {
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
        float* __restrict__ depths;  // [I, N] or [nnz]
        float2* __restrict__ hardness;  // [I, N, 2] or [nnz, 2]
        float3* __restrict__ rgbs;  // [I, N, 3] or [nnz, 3]
        float3* __restrict__ normals;  // [I, N, 3] or [nnz, 3]
        float2* __restrict__ absgrad;  // [I, N, 2] or [nnz, 2]

        Buffer(const Tensor& tensors) {
            DEVICE_GUARD(tensors.vertices);
            CHECK_INPUT(tensors.vertices);
            CHECK_INPUT(tensors.depths);
            CHECK_INPUT(tensors.hardness);
            CHECK_INPUT(tensors.rgbs);
            CHECK_INPUT(tensors.normals);
            vertices = (float2*)tensors.vertices.data_ptr<float>();
            depths = tensors.depths.data_ptr<float>();
            hardness = (float2*)tensors.hardness.data_ptr<float>();
            rgbs = (float3*)tensors.rgbs.data_ptr<float>();
            normals = (float3*)tensors.normals.data_ptr<float>();
            absgrad = tensors.absgrad.has_value() ?
                (float2*)tensors.absgrad.value().data_ptr<float>()
                : nullptr;
        }
    };

#ifdef __CUDACC__

    static __device__ Screen load(const Buffer &buffer, long idx) {
        return {
            buffer.vertices[3*idx+0],
            buffer.vertices[3*idx+1],
            buffer.vertices[3*idx+2],
            buffer.depths[idx],
            buffer.hardness[idx],
            buffer.rgbs[idx],
            buffer.normals[idx]
            // absgrad is undefined
        };
    }

    static __device__ __forceinline__ Screen zero() {
        return {
            {0.f, 0.f},
            {0.f, 0.f},
            {0.f, 0.f},
            0.0f,
            {0.f, 0.f},
            {0.f, 0.f, 0.f},
            {0.f, 0.f, 0.f},
            {0.f, 0.f}
        };
    }

    __device__ __forceinline__ void operator+=(const Screen &other) {
        vert0 += other.vert0;
        vert1 += other.vert1;
        vert2 += other.vert2;
        depth += other.depth;
        hardness += other.hardness;
        rgb += other.rgb;
        normal += other.normal;
        absgrad += fabs(other.vert0 + other.vert1 + other.vert2) / 3.0f;
    }

    __device__ void saveBuffer(Buffer &buffer, long idx) {
        buffer.vertices[3*idx+0] = vert0;
        buffer.vertices[3*idx+1] = vert2;
        buffer.vertices[3*idx+2] = vert1;
        buffer.depths[idx] = depth;
        buffer.hardness[idx] = hardness;
        buffer.rgbs[idx] = rgb;
        buffer.normals[idx] = normal;
        if (buffer.absgrad != nullptr)
            buffer.absgrad[idx] = absgrad;
    }

    __device__ void atomicAddBuffer(Buffer &buffer, long idx) {
        atomicAddFVec(buffer.vertices + 3*idx+0, vert0);
        atomicAddFVec(buffer.vertices + 3*idx+1, vert1);
        atomicAddFVec(buffer.vertices + 3*idx+2, vert2);
        atomicAdd(buffer.depths + idx, depth);
        atomicAddFVec(buffer.hardness + idx, hardness);
        atomicAddFVec(buffer.rgbs + idx, rgb);
        atomicAddFVec(buffer.normals + idx, normal);
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

    __device__ __forceinline__ OpaqueTriangle::RenderOutput evaluate_color(float px, float py) {
        return { rgb, depth, normal };
    }

    __device__ __forceinline__ Screen evaluate_color_vjp(float px, float py, OpaqueTriangle::RenderOutput v_render) {
        Screen v_splat = Screen::zero();
        v_splat.rgb = v_render.rgb;
        v_splat.depth = v_render.depth;
        v_splat.normal = v_render.normal;
        return v_splat;
    }

#endif  // #ifdef __CUDACC__

};


#ifdef __CUDACC__

inline __device__ void OpaqueTriangle::project_persp(
    OpaqueTriangle::World world, OpaqueTriangle::FwdProjCamera cam,
    OpaqueTriangle::Screen& screen, int4& aabb
) {
    projection_opaque_triangle_persp(
        world.mean, world.quat, world.scale, world.hardness, &world.sh_coeffs,
        // world.vert0, world.vert1, world.vert2, world.hardness, &world.sh_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy,
        cam.radial_coeffs, cam.tangential_coeffs, cam.thin_prism_coeffs,
        cam.width, cam.height, cam.near_plane, cam.far_plane,
        &aabb, &screen.vert0, &screen.vert1, &screen.vert2, &screen.depth, &screen.hardness, &screen.rgb, &screen.normal
    );
}

inline __device__ void OpaqueTriangle::project_fisheye(
    OpaqueTriangle::World world, OpaqueTriangle::FwdProjCamera cam,
    OpaqueTriangle::Screen& screen, int4& aabb
) {
    projection_opaque_triangle_fisheye(
        world.mean, world.quat, world.scale, world.hardness, &world.sh_coeffs,
        // world.vert0, world.vert1, world.vert2, world.hardness, &world.sh_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy,
        cam.radial_coeffs, cam.tangential_coeffs, cam.thin_prism_coeffs,
        cam.width, cam.height, cam.near_plane, cam.far_plane,
        &aabb, &screen.vert0, &screen.vert1, &screen.vert2, &screen.depth, &screen.hardness, &screen.rgb, &screen.normal
    );
}

inline __device__ void OpaqueTriangle::project_persp_vjp(
    OpaqueTriangle::World world, OpaqueTriangle::BwdProjCamera cam,
    OpaqueTriangle::Screen v_screen,
    OpaqueTriangle::World& v_world, float3x3 &v_R, float3 &v_t
) {
    projection_opaque_triangle_persp_vjp(
        world.mean, world.quat, world.scale, world.hardness, &world.sh_coeffs,
        // world.vert0, world.vert1, world.vert2, world.hardness, &world.sh_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy,
        cam.radial_coeffs, cam.tangential_coeffs, cam.thin_prism_coeffs,
        cam.width, cam.height,
        v_screen.vert0, v_screen.vert1, v_screen.vert2, v_screen.depth, v_screen.hardness, v_screen.rgb, v_screen.normal,
        &v_world.mean, &v_world.quat, &v_world.scale, &v_world.hardness, &v_world.sh_coeffs,
        // &v_world.vert0, &v_world.vert1, &v_world.vert2, &v_world.hardness, &v_world.sh_coeffs,
        &v_R, &v_t
    );
}

inline __device__ void OpaqueTriangle::project_fisheye_vjp(
    OpaqueTriangle::World world, OpaqueTriangle::BwdProjCamera cam,
    OpaqueTriangle::Screen v_screen,
    OpaqueTriangle::World& v_world, float3x3 &v_R, float3 &v_t
) {
    projection_opaque_triangle_fisheye_vjp(
        world.mean, world.quat, world.scale, world.hardness, &world.sh_coeffs,
        // world.vert0, world.vert1, world.vert2, world.hardness, &world.sh_coeffs,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy,
        cam.radial_coeffs, cam.tangential_coeffs, cam.thin_prism_coeffs,
        cam.width, cam.height,
        v_screen.vert0, v_screen.vert1, v_screen.vert2, v_screen.depth, v_screen.hardness, v_screen.rgb, v_screen.normal,
        &v_world.mean, &v_world.quat, &v_world.scale, &v_world.hardness, &v_world.sh_coeffs,
        // &v_world.vert0, &v_world.vert1, &v_world.vert2, &v_world.hardness, &v_world.sh_coeffs,
        &v_R, &v_t
    );
}

#endif  // #ifdef __CUDACC__
