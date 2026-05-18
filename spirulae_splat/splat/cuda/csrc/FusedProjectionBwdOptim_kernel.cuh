#include <cuda_runtime.h>
#include <cstdint>

#include <gsplat/Common.h>
#include <gsplat/Utils.cuh>

#ifndef NO_TORCH
#define NO_TORCH
#endif

#include "generated/slang.cuh"
namespace SlangPerSplatLosses {
#include "generated/set_namespace.cuh"
#include "generated/per_splat_losses.cuh"
}
namespace SlangDensify {
#include "generated/set_namespace.cuh"
#include "generated/densify.cuh"
}

#include "types.cuh"

#include <cooperative_groups.h>
namespace cg = cooperative_groups;


inline constexpr float kSh0 = 0.28209479177387814f;

__forceinline__ __device__ float3 fmaxf(float3 v, float k) {
    return {
        fmaxf(v.x, k),
        fmaxf(v.y, k),
        fmaxf(v.z, k)
    };
}
__forceinline__ __device__ float3 sqrtf(float3 v) {
    return {
        sqrtf(fmaxf(v.x, 0.0f)),
        sqrtf(fmaxf(v.y, 0.0f)),
        sqrtf(fmaxf(v.z, 0.0f))
    };
}

__forceinline__ __device__ float4 sqrtf(float4 v) {
    return {
        sqrtf(fmaxf(v.x, 0.0f)),
        sqrtf(fmaxf(v.y, 0.0f)),
        sqrtf(fmaxf(v.z, 0.0f)),
        sqrtf(fmaxf(v.w, 0.0f))
    };
}


template<
    typename SplatPrimitive,
    ssplat::CameraModelType camera_model,
    HessianDiagonalOutputMode hessian_diagonal_output_mode,
    bool use_scale_agnostic_mean
>
#if 0
__global__ void __launch_bounds__(512) fused_projection_bwd_optimizer_3dgs_kernel
#else
__global__ void fused_projection_bwd_optimizer_3dgs_kernel
#endif
(
    // fwd inputs
    const uint32_t C,
    const uint32_t N,
    typename SplatPrimitive::WorldBuffer splats_world,
    const float *__restrict__ viewmats, // [C, 4, 4]
    const float4 *__restrict__ intrins,  // [C, 4], fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const uint32_t image_width,
    const uint32_t image_height,
    // fwd outputs
    const int32_t *__restrict__ camera_id_bounds,   // [N]
    const int32_t *__restrict__ camera_ids,   // [nnz]
    const float4 *__restrict__ aabb,   // [C, N, 4] or [nnz, 4]
    // grad outputs from rasterization
    typename SplatPrimitive::WorldBuffer v_splats_world,
    typename SplatPrimitive::WorldBuffer vr_splats_world,
    typename SplatPrimitive::WorldBuffer h_splats_world,
    typename SplatPrimitive::ScreenBuffer v_splats_screen,
    typename SplatPrimitive::ScreenBuffer vr_splats_screen,
    typename SplatPrimitive::ScreenBuffer h_splats_screen,
    // optimizer states
    typename SplatPrimitive::WorldBuffer g1_splats_world,
    typename SplatPrimitive::WorldBuffer g2_splats_world,
    // float *__restrict__ v_viewmats // [C, 4, 4] optional
    // optimizer params
    const float* __restrict__ radii,
    const float lr_means,
    const float lr_quats,
    const float lr_scales,
    const float lr_opacs,
    const float lr_features_dc,
    const float lr_features_sh,
    const float mcmc_noise_scalar,
    const float min_opacity,
    const float max_gauss_ratio,
    const float scale_regularization_weight,
    const float mcmc_opacity_reg_weight,
    const float mcmc_scale_reg_weight,
    const float erank_reg_weight,
    const float erank_reg_weight_s3,
    const float quat_norm_reg_weight,
    const float mrnf_opacity_decay_factor,
    const float mrnf_scale_decay_factor,
    const int32_t scalar_step,
    const int32_t* __restrict__ steps
) {
    constexpr uint32_t BLOCK_SIZE = WARP_SIZE;

    uint32_t gid = blockIdx.x * BLOCK_SIZE + threadIdx.x;
    if (gid >= N) return;

    bool packed = (camera_id_bounds != nullptr && camera_ids != nullptr);
    int cid_0 = 0, cid_1 = C;
    if (packed) {
        cid_0 = (gid == 0 ? 0 : camera_id_bounds[gid-1]);
        cid_1 = camera_id_bounds[gid];
    }

    // Load splat
    typename SplatPrimitive::World splat_world;
    splat_world.load(splats_world, gid);

    // TODO: put SH degree as a template argument
    constexpr int num_sh = SplatPrimitive::World::num_sh();
    // __shared__ float3 v_sh_coeffs[num_sh*BLOCK_SIZE];
    float3 v_sh_coeffs[num_sh == 0 ? 1 : num_sh];

    typename SplatPrimitive::World v_splat_world = SplatPrimitive::World::zero();
    v_splat_world.atomicLoad(v_splats_world, gid);
    // v_splat_world.sh_degree = splat_world.sh_degree;
    // v_splat_world.features_sh = &v_sh_coeffs[num_sh*threadIdx.x];
    v_splat_world.features_sh = &v_sh_coeffs[0];
    #pragma unroll
    for (int i = 0; i < num_sh; ++i)
        v_splat_world.features_sh[i] = make_float3(0.0f);

    float3x3 v_R = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
    float3 v_t = {0.f, 0.f, 0.f};
    float3 vr_world_pos = {0.f, 0.f, 0.f};
    float3 h_world_pos = {0.f, 0.f, 0.f};
    if constexpr (hessian_diagonal_output_mode == HessianDiagonalOutputMode::Position) {
        vr_world_pos = vr_splats_world.means(gid);
        h_world_pos = h_splats_world.means(gid);
    }
    typename SplatPrimitive::World vr_splat_world = SplatPrimitive::World::zero();
    typename SplatPrimitive::World h_splat_world = SplatPrimitive::World::zero();
    if constexpr (hessian_diagonal_output_mode == HessianDiagonalOutputMode::AllReasonable) {
        vr_splat_world.atomicLoad(vr_splat_world, gid);
        h_splat_world.atomicLoad(h_splat_world, gid);
    }

    // Loop over intersections
    for (int cid_t = cid_0; cid_t < cid_1; ++cid_t) {
        int cid = (packed ? camera_ids[cid_t] : cid_t);
        int idx = packed ? cid_t : cid_t * N + gid;
        if (aabb[idx].z <= aabb[idx].x || aabb[idx].w <= aabb[idx].y)
            continue;

        // Load camera
        float4 intrin = intrins[cid];
        float3x3 R = {
            viewmats[cid*16+0], viewmats[cid*16+1], viewmats[cid*16+2],  // 1st row
            viewmats[cid*16+4], viewmats[cid*16+5], viewmats[cid*16+6],  // 2nd row
            viewmats[cid*16+8], viewmats[cid*16+9], viewmats[cid*16+10],  // 3rd row
        };
        float3 t = { viewmats[cid*16+3], viewmats[cid*16+7], viewmats[cid*16+11] };
        float fx = intrin.x, fy = intrin.y, cx = intrin.z, cy = intrin.w;
        ProjCamera cam = {
            R, t, fx, fy, cx, cy,
            image_width, image_height,
        };
        cam.dist_coeffs = dist_coeffs_buffer.load(cid);

        // Load splat gradient
        typename SplatPrimitive::Screen v_splat_screen;
        typename SplatPrimitive::Screen vr_splat_screen;
        typename SplatPrimitive::Screen h_splat_screen;
        v_splat_screen.load(v_splats_screen, idx);
        if (hessian_diagonal_output_mode != HessianDiagonalOutputMode::None) {
            vr_splat_screen.load(vr_splats_screen, idx);
            h_splat_screen.load(h_splats_screen, idx);
        }

        // Accumulate gradient
        if constexpr (hessian_diagonal_output_mode == HessianDiagonalOutputMode::None) {
            splat_world.template project_vjp<camera_model>(cam, v_splat_screen, v_splat_world, v_R, v_t);
        } else if constexpr (hessian_diagonal_output_mode == HessianDiagonalOutputMode::Position) {
            splat_world.template project_vjp_h_pos<camera_model>(cam, v_splat_screen,
                vr_splat_screen, h_splat_screen, v_splat_world, v_R, v_t, vr_world_pos, h_world_pos);
        } else if constexpr (hessian_diagonal_output_mode == HessianDiagonalOutputMode::AllReasonable) {
            splat_world.template project_vjp_h_all<camera_model>(cam, v_splat_screen,
                vr_splat_screen, h_splat_screen, v_splat_world, v_R, v_t, vr_splat_world, h_splat_world);
        }

        // TODO: affects performance, refactor this as a template argument?
    #if 0
        // gradient to viewmats
        if (v_viewmats != nullptr) {
            auto warp = cg::tiled_partition<WARP_SIZE>(cg::this_thread_block());
            auto warp_group_c = cg::labeled_partition(warp, cid);
            warpSum(v_R[0], warp_group_c);
            warpSum(v_R[1], warp_group_c);
            warpSum(v_R[2], warp_group_c);
            warpSum(v_t, warp_group_c);
            if (warp_group_c.thread_rank() == 0) {
                v_viewmats += cid * 16;
                #pragma unroll
                for (uint32_t i = 0; i < 3; i++) { // rows
                    atomicAdd(v_viewmats + i * 4 + 0, v_R[i].x);
                    atomicAdd(v_viewmats + i * 4 + 1, v_R[i].y);
                    atomicAdd(v_viewmats + i * 4 + 2, v_R[i].z);
                }
                atomicAdd(v_viewmats + 0 * 4 + 3, v_t.x);
                atomicAdd(v_viewmats + 1 * 4 + 3, v_t.y);
                atomicAdd(v_viewmats + 2 * 4 + 3, v_t.z);
            }
        }
    #endif

    }

    // TODO: second-order optimizer

    // add regularization to gradient
    static constexpr int kNumPerSplatLosses = 5;
    FixedArray<float, kNumPerSplatLosses> v_losses = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
    float3 v_scale_t = make_float3(0.0f);
    float4 v_quat_t = make_float4(0.0f);
    float v_opac_t = 0.0f;
    SlangPerSplatLosses::per_splat_losses_bwd(
        splat_world.scale, splat_world.opacity, splat_world.quat,
        v_losses,
        &v_scale_t, &v_opac_t, &v_quat_t,
        mcmc_opacity_reg_weight,
        mcmc_scale_reg_weight,
        max_gauss_ratio,
        scale_regularization_weight,
        erank_reg_weight,
        erank_reg_weight_s3,
        quat_norm_reg_weight
    );
    v_splat_world.scale += v_scale_t;
    v_splat_world.quat += v_quat_t;
    v_splat_world.opacity += v_opac_t;

    // optimizer for mean/quat/scale

    static constexpr float eps = 1e-15f;
    static constexpr float beta1 = 0.9f;
    static constexpr float beta2 = 0.999f;
    float step = (float)(steps != nullptr ? steps[gid] : scalar_step);
    float inv_bias_correction1 = 1.0f / (1.0f - powf(beta1, step));
    float inv_bias_correction2 = 1.0f / (1.0f - powf(beta2, step));

    // update scales
    float3 g1_scale = beta1 * g1_splats_world.scales(gid) + (1.f - beta1) * v_splat_world.scale;
    float3 g2_scale = beta2 * g2_splats_world.scales(gid) + (1.f - beta2) * v_splat_world.scale*v_splat_world.scale;
    float3 new_scale = splat_world.scale - lr_scales * inv_bias_correction1 * g1_scale / (sqrtf(g2_scale * inv_bias_correction2) + eps);
    if (mrnf_scale_decay_factor != 1.0f)
        new_scale += make_float3(__logf(mrnf_scale_decay_factor));
    splats_world.scales(gid) = new_scale;
    g1_splats_world.scales(gid) = g1_scale;
    g2_splats_world.scales(gid) = g2_scale;

    // update quats (Riemannian)
    float4 new_quat = normalize(splat_world.quat);
    v_splat_world.quat -= dot(new_quat, v_splat_world.quat) * new_quat;
    float4 g1_quat = beta1 * g1_splats_world.quats(gid) + (1.f - beta1) * v_splat_world.quat;
    float4 g2_quat = beta2 * g2_splats_world.quats(gid) + (1.f - beta2) * v_splat_world.quat*v_splat_world.quat;
    new_quat -= lr_quats * inv_bias_correction1 * g1_quat / (sqrtf(g2_quat * inv_bias_correction2) + eps);
    splats_world.quats(gid) = normalize(new_quat);
    g1_splats_world.quats(gid) = g1_quat;
    g2_splats_world.quats(gid) = g2_quat;

    // update opacs
    float g1_opac = beta1 * g1_splats_world.opacities(gid) + (1.f - beta1) * v_splat_world.opacity;
    float g2_opac = beta2 * g2_splats_world.opacities(gid) + (1.f - beta2) * v_splat_world.opacity*v_splat_world.opacity;
    float new_opac = splat_world.opacity - lr_opacs * inv_bias_correction1 * g1_opac / (sqrtf(g2_opac * inv_bias_correction2) + eps);
    if (mrnf_opacity_decay_factor != 0.0f) {
        new_opac = sigmoid(new_opac);
        new_opac = fmaxf(new_opac + mrnf_opacity_decay_factor, 1e-12f);
        new_opac = logit(new_opac);
    }
    splats_world.opacities(gid) = new_opac;
    g1_splats_world.opacities(gid) = g1_opac;
    g2_splats_world.opacities(gid) = g2_opac;

    // update means (scale agnostic)
    float3 g1_mean = g1_splats_world.means(gid);
    float3 g2_mean = g2_splats_world.means(gid);
    float noise_lr_scalar = 1.0f;
    if constexpr (use_scale_agnostic_mean) {
        float3 v_mean_scaled_num = SlangProjectionUtils::apply_covar_to_vec(
            splat_world.quat,
            {expf(0.5f*splat_world.scale.x), expf(0.5f*splat_world.scale.y), expf(0.5f*splat_world.scale.z)},  // unit: L^0.5
            v_splat_world.mean  // unit: L^-1
        ) * sqrtf(2.0f * __logf(fmaxf(255.0f * sigmoid(splat_world.opacity), 1.00001f)));  // unit: dimensionless
        float v_mean_scaled_den = radii[gid] * 0.6f;  // unit: dimensionless
        g1_mean = beta1 * g1_mean + (1.f - beta1) * v_mean_scaled_num;  // unit: dimensionless
        g2_mean = beta2 * g2_mean + (1.f - beta2) * v_splat_world.mean*v_splat_world.mean * v_mean_scaled_den*v_mean_scaled_den;  // unit: L^-2
        // TODO: probably better use a globally consistent one; (MRNF chooses based on median scale)
        noise_lr_scalar = length(v_mean_scaled_num) / (length(v_splat_world.mean) * v_mean_scaled_den + eps);
    } else {
        g1_mean = beta1 * g1_mean + (1.f - beta1) * v_splat_world.mean;  // unit: L^-1
        g2_mean = beta2 * g2_mean + (1.f - beta2) * v_splat_world.mean*v_splat_world.mean;  // unit: L^-2
    }
    float3 new_mean = splat_world.mean;
    SlangDensify::mcmc_add_noise_3dgs(
        mcmc_noise_scalar * noise_lr_scalar, min_opacity,
        &new_mean, splat_world.scale, splat_world.quat, sigmoid(splat_world.opacity)
        // official MCMC use scale/quat/opac after optimizer step/densification, shouldn't matter in practice
    );
    splats_world.means(gid) = new_mean - lr_means * inv_bias_correction1 * g1_mean /
        (sqrtf(g2_mean * inv_bias_correction2) + eps); // unit: L or dimensionless for dimensionless lr_means
    g1_splats_world.means(gid) = g1_mean;
    g2_splats_world.means(gid) = g2_mean;

    // update features_dc
    float3 g1_feature_dc = beta1 * g1_splats_world.features_dc(gid) + (1.f - beta1) * v_splat_world.features_dc;
    float3 g2_feature_dc = beta2 * g2_splats_world.features_dc(gid) + (1.f - beta2) * v_splat_world.features_dc*v_splat_world.features_dc;
    splats_world.features_dc(gid) = splat_world.features_dc - lr_features_dc * inv_bias_correction1
        * g1_feature_dc / (sqrtf(g2_feature_dc * inv_bias_correction2) + eps);
    g1_splats_world.features_dc(gid) = g1_feature_dc;
    g2_splats_world.features_dc(gid) = g2_feature_dc;

    // update features_sh
    // TODO: more cache-friendly memory access pattern
    #pragma unroll
    for (int i = 0; i < num_sh; ++i) {
        float3 v_sh_coeff = v_splat_world.features_sh[i];
        float3 g1_feature_sh = beta1 * g1_splats_world.features_sh(gid, i) + (1.f - beta1) * v_sh_coeff;
        float3 g2_feature_sh = beta2 * g2_splats_world.features_sh(gid, i) + (1.f - beta2) * v_sh_coeff*v_sh_coeff;
        splats_world.features_sh(gid, i) = splat_world.features_sh[i] - lr_features_sh * inv_bias_correction1
            * g1_feature_sh / (sqrtf(g2_feature_sh * inv_bias_correction2) + eps);
        g1_splats_world.features_sh(gid, i) = g1_feature_sh;
        g2_splats_world.features_sh(gid, i) = g2_feature_sh;
    }
}

template<
    typename SplatPrimitive,
    ssplat::CameraModelType camera_model,
    HessianDiagonalOutputMode hessian_diagonal_output_mode,
    const bool use_scale_agnostic_mean
>
void fused_projection_bwd_optimizer_3dgs_kernel_wrapper(
    cudaStream_t stream,
    // fwd inputs
    const uint32_t C,
    const uint32_t N,
    typename SplatPrimitive::WorldBuffer splats_world,
    const float *__restrict__ viewmats, // [C, 4, 4]
    const float4 *__restrict__ intrins,  // [C, 4], fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const uint32_t image_width,
    const uint32_t image_height,
    // fwd outputs
    const int32_t *__restrict__ camera_id_bounds,   // [N]
    const int32_t *__restrict__ camera_ids,   // [nnz]
    const float4 *__restrict__ aabb,   // [C, N, 4] or [nnz, 4]
    // grad outputs from rasterization
    typename SplatPrimitive::WorldBuffer v_splats_world,
    typename SplatPrimitive::WorldBuffer vr_splats_world,
    typename SplatPrimitive::WorldBuffer h_splats_world,
    typename SplatPrimitive::ScreenBuffer v_splats_screen,
    typename SplatPrimitive::ScreenBuffer vr_splats_screen,
    typename SplatPrimitive::ScreenBuffer h_splats_screen,
    // optimizer states
    typename SplatPrimitive::WorldBuffer g1_splats_world,
    typename SplatPrimitive::WorldBuffer g2_splats_world,
    // float *__restrict__ v_viewmats // [C, 4, 4] optional
    // optimizer params
    const float* __restrict__ radii,
    const float lr_means,
    const float lr_quats,
    const float lr_scales,
    const float lr_opacs,
    const float lr_features_dc,
    const float lr_features_sh,
    const float mcmc_noise_scalar,
    const float min_opacity,
    const float max_gauss_ratio,
    const float scale_regularization_weight,
    const float mcmc_opacity_reg_weight,
    const float mcmc_scale_reg_weight,
    const float erank_reg_weight,
    const float erank_reg_weight_s3,
    const float quat_norm_reg_weight,
    const float mrnf_opacity_decay_factor,
    const float mrnf_scale_decay_factor,
    const int32_t scalar_step,
    const int32_t* __restrict__ steps
) {
    fused_projection_bwd_optimizer_3dgs_kernel<
        SplatPrimitive, camera_model, hessian_diagonal_output_mode, use_scale_agnostic_mean
    >
    <<<_CEIL_DIV(N, WARP_SIZE), WARP_SIZE, 0, stream>>>(
        C, N,
        splats_world, viewmats, intrins, dist_coeffs_buffer, image_width, image_height,
        camera_id_bounds, camera_ids, aabb,
        v_splats_world, vr_splats_world, h_splats_world, v_splats_screen, vr_splats_screen, h_splats_screen,
        g1_splats_world, g2_splats_world, //v_viewmats,
        radii, lr_means, lr_quats, lr_scales, lr_opacs, lr_features_dc, lr_features_sh,
        mcmc_noise_scalar * lr_means,
        min_opacity,
        max_gauss_ratio,
        scale_regularization_weight / (float)N,  // TODO: use actual N when max gaussian count hasn't reached
        mcmc_opacity_reg_weight / (float)N,
        mcmc_scale_reg_weight / (float)N,
        erank_reg_weight / (float)N,
        erank_reg_weight_s3 / (float)N,
        quat_norm_reg_weight / (float)N,
        mrnf_opacity_decay_factor,
        mrnf_scale_decay_factor,
        scalar_step, steps
    );
}
