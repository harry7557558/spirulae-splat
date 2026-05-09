#pragma once

#ifdef __CUDACC__
#include "generated/slang.cuh"
namespace SlangProjectionUtils {
#include "generated/set_namespace.cuh"
#include "generated/projection_utils.cuh"
}
#endif

#include "Primitive.cuh"


struct _BasePrimitive3DGS {
    class WorldBuffer;
    class ScreenBuffer;

    #ifdef __CUDACC__
    struct World {
        float3 mean;
        float4 quat;
        float3 scale;  // log
        float opacity;  // logit
        // TODO: register pressure
        FixedArray<float3, 16> sh_coeffs;
    };
    #endif

    class WorldBuffer : public TensorArray<6> {
    public:
        using TensorArray<6>::TensorArray;
    #ifdef __CUDACC__
        __forceinline__ __device__ float3& means(int64_t i)
            { return *reinterpret_cast<float3*>(&_data[0][3*i]); }
        __forceinline__ __device__ float4& quats(int64_t i)
            { return *reinterpret_cast<float4*>(&_data[1][4*i]); }
        __forceinline__ __device__ float3& scales(int64_t i)
            { return *reinterpret_cast<float3*>(&_data[2][3*i]); }
        __forceinline__ __device__ float& opacities(int64_t i)
            { return *reinterpret_cast<float*>(&_data[3][i]); }
        __forceinline__ __device__ float3& features_dc(int64_t i)
            { return *reinterpret_cast<float3*>(&_data[4][3*i]); }
        __forceinline__ __device__ int32_t num_sh() const
            { return _strides[5] / 3; }
        __forceinline__ __device__ float3& features_sh(int64_t i, int64_t j)
            { return *reinterpret_cast<float3*>(&_data[5][_strides[5]*i+j]); }

        __forceinline__ __device__ World load(int64_t i) {
            World splat {
                means(i),
                quats(i),
                scales(i),
                opacities(i),
            };
            splat.sh_coeffs[0] = features_dc(i);
            for (int j = 0; j < 15; ++j)
                splat.sh_coeffs[j+1] = j < num_sh() ?
                    features_sh(i, j) : make_float3(0);
            return splat;
        }

        __forceinline__ __device__ void store(int64_t i, const World &world) {
            if (_data[0]) means(i) = world.mean;
            if (_data[1]) quats(i) = world.quat;
            if (_data[2]) scales(i) = world.scale;
            if (_data[3]) opacities(i) = world.opacity;
            if (_data[4]) features_dc(i) = world.sh_coeffs[0];
            if (_data[5]) for (int j = 0; j < num_sh(); ++j)
                features_sh(i, j) = world.sh_coeffs[j+1];
        }

        __forceinline__ __device__ void atomicStore(int64_t i, const World &world) {
            if (_data[0]) atomicAddFVec(&means(i), world.mean);
            if (_data[1]) atomicAddFVec(&quats(i), world.quat);
            if (_data[2]) atomicAddFVec(&scales(i), world.scale);
            if (_data[3]) atomicAddFVec(&opacities(i), world.opacity);
            if (_data[4]) atomicAddFVec(&features_dc(i), world.sh_coeffs[0]);
            if (_data[5]) for (int j = 0; j < num_sh(); ++j)
                atomicAddFVec(&features_sh(i, j), world.sh_coeffs[j+1]);
        }
        
    #endif  // #ifdef __CUDACC__

    #ifndef NO_TORCH
        static TensorList empty(int64_t size, int num_sh) {
            return TensorArray<6>::empty(size, {3, 4, 3, 1, 3, (int32_t)(3*num_sh)});
        }
        static TensorList zeros(int64_t size, int num_sh) {
            return TensorArray<6>::zeros(size, {3, 4, 3, 1, 3, (int32_t)(3*num_sh)});
        }
    #endif  // #ifndef NO_TORCH
    };

    #ifdef __CUDACC__
    struct Screen {
        float3 scale;  // post exp
        float opacity;  // post sigmoid
        float3 rgb;
    };
    #endif

    class ScreenBuffer : public TensorArray<3> {
    public:
        using TensorArray<3>::TensorArray;
    #ifdef __CUDACC__
        __forceinline__ __device__ float3& scales(int64_t i)
            { return *reinterpret_cast<float3*>(&_data[0][3*i]); }
        __forceinline__ __device__ float& opacities(int64_t i)
            { return *reinterpret_cast<float*>(&_data[1][i]); }
        __forceinline__ __device__ float3& colors(int64_t i)
            { return *reinterpret_cast<float3*>(&_data[2][3*i]); }

        __forceinline__ __device__ Screen load(int64_t i) {
            Screen screen;
            screen.scale = scales(i);
            screen.opacity = opacities(i);
            screen.rgb = colors(i);
            return screen;
        }

        __forceinline__ __device__ void store(int64_t i, const Screen &screen) {
            if (_data[0]) scales(i) = screen.scale;
            if (_data[1]) opacities(i) = screen.opacity;
            if (_data[2]) colors(i) = screen.rgb;
        }

        __forceinline__ __device__ void atomicStore(int64_t i, const Screen &screen) {
            if (_data[0]) atomicAddFVec(&scales(i), screen.scale);
            if (_data[1]) atomicAddFVec(&opacities(i), screen.opacity);
            if (_data[2]) atomicAddFVec(&colors(i), screen.rgb);
        }
    #endif  // #ifdef __CUDACC__

    #ifndef NO_TORCH
        static TensorList empty(int64_t size) {
            return TensorArray<3>::empty(size, {3, 1, 3});
        }
        static TensorList zeros(int64_t size) {
            return TensorArray<3>::zeros(size, {3, 1, 3});
        }
    #endif  // #ifndef NO_TORCH
    };

};


// struct _Base3DGS : public _BasePrimitive3DGS {

//     static constexpr RenderOutputType pixelType = RenderOutputType::RGB_D;


//     inline static __device__ void project_persp(
//         World world, ProjCamera cam,
//         Screen& screen, float4& aabb
//     );

//     inline static __device__ void project_ortho(
//         World world, ProjCamera cam,
//         Screen& screen, float4& aabb
//     );

//     inline static __device__ void project_fisheye(
//         World world, ProjCamera cam,
//         Screen& screen, float4& aabb
//     );

//     inline static __device__ void project_persp_vjp(
//         World world, ProjCamera cam,
//         Screen v_screen,
//         World& v_world, float3x3 &v_R, float3 &v_t
//     );

//     inline static __device__ void project_persp_vjp(
//         World world, ProjCamera cam,
//         Screen v_proj, Screen vr_proj, Screen h_proj,
//         World& v_world, float3x3 &v_R, float3 &v_t,
//         float3 &vr_world, float3 &h_world
//     );

//     inline static __device__ void project_persp_vjp(
//         World world, ProjCamera cam,
//         Screen v_proj, Screen vr_proj, Screen h_proj,
//         World& v_world, float3x3 &v_R, float3 &v_t,
//         World& vr_world, World& h_world
//     );

//     inline static __device__ void project_ortho_vjp(
//         World world, ProjCamera cam,
//         Screen v_screen,
//         World& v_world, float3x3 &v_R, float3 &v_t
//     );

//     inline static __device__ void project_fisheye_vjp(
//         World world, ProjCamera cam,
//         Screen v_screen,
//         World& v_world, float3x3 &v_R, float3 &v_t
//     );

//     inline static __device__ void project_fisheye_vjp(
//         World world, ProjCamera cam,
//         Screen v_proj, Screen vr_proj, Screen h_proj,
//         World& v_world, float3x3 &v_R, float3 &v_t,
//         float3 &vr_world_pos, float3 &h_world_pos
//     );

//     inline static __device__ void project_fisheye_vjp(
//         World world, ProjCamera cam,
//         Screen v_proj, Screen vr_proj, Screen h_proj,
//         World& v_world, float3x3 &v_R, float3 &v_t,
//         World& vr_world, World& h_world
//     );

// };


// struct _Base3DGUT : public _BasePrimitive3DGS {

// }

#if 0

struct Base3DGUT::Fragment {

    // from world
    float3 mean;
    float4 quat;
    // from screen
    float3 scale;  // post exp
    float opacity;  // post sigmoid
    float3 rgb;
    // precomputed
#ifdef __CUDACC__
    float3x3 iscl_rot;
#endif

    struct Buffer {
        float3* __restrict__ scales;
        float* __restrict__ opacities;
        float3* __restrict__ rgbs;
        long size;

        Buffer() {}  // uninitialized

    #ifndef NO_TORCH
        Buffer(const Tensor& tensors) {
            DEVICE_GUARD(tensors.scales);
            CHECK_INPUT(tensors.scales);
            CHECK_INPUT(tensors.opacities);
            CHECK_INPUT(tensors.rgbs);
            scales = (float3*)tensors.scales.data_ptr<float>();
            opacities = tensors.opacities.data_ptr<float>();
            rgbs = (float3*)tensors.rgbs.data_ptr<float>();
            size = tensors.opacities.numel();
        }
    #endif
    };

#ifdef __CUDACC__

    static __device__ Screen load(const Buffer &buffer, long idx, const uint32_t* gaussian_ids) {
        uint32_t idx0 = gaussian_ids ? gaussian_ids[idx] : idx % buffer.size;
        return {
            buffer.means ? buffer.means[idx0] : make_float3(0.f),
            buffer.quats ? buffer.quats[idx0] : make_float4(0.f),
            buffer.depths ? buffer.depths[idx] : 0.0f,
            buffer.scales[idx],
            buffer.opacities[idx],
            // buffer.scales[idx % buffer.size],
            // buffer.opacities[idx % buffer.size],
            buffer.rgbs[idx],
        };
    }

    static __device__ Screen loadWithPrecompute(const Buffer &buffer, long idx, const uint32_t* gaussian_ids) {
        Screen splat = Screen::load(buffer, idx, gaussian_ids);
        splat.iscl_rot = SlangProjectionUtils::compute_3dgut_iscl_rot(splat.quat, splat.scale);
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

    __device__ __forceinline__ void precomputeBackward(Screen& grad) const {
        float4 v_quat; float3 v_scale;
        SlangProjectionUtils::compute_3dgut_iscl_rot_vjp(quat, scale, grad.iscl_rot, &v_quat, &v_scale);
        grad.quat += v_quat;
        grad.scale += v_scale;
    }

    __device__ __forceinline__ void addGradient(const Screen &grad, float weight=1.0f) {
        // optional precomputeBackward
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
        // precomputeBackward needed
        mean += fmul_axa(grad.mean, weight);
        quat += fmul_axa(grad.quat, weight);
        depth += fmul_axa(grad.depth, weight);
        scale += fmul_axa(grad.scale, weight);
        opacity += fmul_axa(grad.opacity, weight);
        rgb += fmul_axa(grad.rgb, weight);
    }

    __device__ void saveParamsToBuffer(Buffer &buffer, long idx, const uint32_t* gaussian_ids) {
        uint32_t idx0 = gaussian_ids ? gaussian_ids[idx] : idx % buffer.size;
        if (buffer.means) buffer.means[idx0] = mean;
        if (buffer.quats) buffer.quats[idx0] = quat;
        if (buffer.depths) buffer.depths[idx] = depth;
        buffer.scales[idx] = scale;
        buffer.opacities[idx] = opacity;
        // buffer.scales[idx % buffer.size] = scale;
        // buffer.opacities[idx % buffer.size] = opacity;
        buffer.rgbs[idx] = rgb;
        // iscl_rot is not saved
    }

    template<int reduce = 1>
    __device__ void atomicAddToBuffer(Buffer &buffer, long idx, const uint32_t* gaussian_ids) const {
        uint32_t idx0 = gaussian_ids ? gaussian_ids[idx] : idx % buffer.size;
        // precomputeBackward needed
        atomicAddFVec<reduce>(buffer.means + idx0, mean);
        atomicAddFVec<reduce>(buffer.quats + idx0, quat);
        // atomicAddFVec<reduce>(buffer.depths + idx, depth);
        atomicAddFVec<reduce>(buffer.scales + idx, scale);
        atomicAddFVec<reduce>(buffer.opacities + idx, opacity);
        // atomicAddFVec<reduce>(buffer.scales + idx % buffer.size, scale);
        // atomicAddFVec<reduce>(buffer.opacities + idx % buffer.size, opacity);
        atomicAddFVec<reduce>(buffer.rgbs + idx, rgb);
    }

    __device__ __forceinline__ float evaluate_alpha(float3 ray_o, float3 ray_d) {
        // if (dot(mean-ray_o, ray_d) <= 0.0f)
        //     return 0.0;
        return SlangProjectionUtils::evaluate_alpha_3dgs(
            mean, iscl_rot, opacity,
            ray_o, ray_d
        );
    }

    __device__ __forceinline__ Screen evaluate_alpha_vjp(
        float3 ray_o, float3 ray_d, float v_alpha,
        float3 &v_ray_o, float3 &v_ray_d
    ) {
        Screen v_splat = Screen::zero();
        // if (dot(mean-ray_o, ray_d) <= 0.0f) {
        //     v_ray_o = v_ray_d = make_float3(0.f);
        //     return v_splat;
        // }
        SlangProjectionUtils::evaluate_alpha_3dgs_vjp(
            mean, iscl_rot, opacity,
            ray_o, ray_d, v_alpha,
            &v_splat.mean, &v_splat.iscl_rot, &v_splat.opacity,
            &v_ray_o, &v_ray_d
        );
        return v_splat;
    }

    __device__ __forceinline__ RenderOutput evaluate_color(float3 ray_o, float3 ray_d) {
        // return {rgb, depth};
        float3 out_rgb; float out_depth;
        SlangProjectionUtils::evaluate_color_3dgs(
            mean, iscl_rot, opacity, rgb,
            ray_o, ray_d, &out_rgb, &out_depth
        );
        return {out_rgb, out_depth, make_float3(0.0f)};
    }

    __device__ __forceinline__ Screen evaluate_color_vjp(
        float3 ray_o, float3 ray_d, RenderOutput v_render,
        float3 &v_ray_o, float3 &v_ray_d
    ) {
        Screen v_splat = Screen::zero();
    #if 0
        v_splat.rgb = v_render.rgb;
        v_splat.depth = v_render.depth;
        v_ray_o = make_float3(0.f);
        v_ray_d = make_float3(0.f);
    #else
        SlangProjectionUtils::evaluate_color_3dgs_vjp(
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

#endif
