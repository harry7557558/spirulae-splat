#include <cuda_runtime.h>
#include <cstdint>

#include <Common.cuh>

#ifndef NO_TORCH
#define NO_TORCH
#endif

#include "generated/slang.cuh"
namespace SlangPerSplatLosses {
#include "generated/set_namespace.cuh"
#include "generated/per_splat_losses.cuh"
}
namespace SlangPixelWise {
#include "generated/set_namespace.cuh"
#include "generated/pixel_wise.cuh"
}

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
__forceinline__ __device__ float3 fminf(float3 v, float k) {
    return {
        fminf(v.x, k),
        fminf(v.y, k),
        fminf(v.z, k)
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
    CameraModelType camera_model,
    HessianDiagonalOutputMode hessian_diagonal_output_mode,
    bool use_scale_agnostic_mean,
    // Merged flag for the two color-space variants. In Python both
    // `OptimConfig::use_color_trust_region` and `OptimConfig::color_is_linear`
    // are wired to the same `splat_color_is_linear` flag (core.py), so
    // collapsing the bool2 -> bool here halves the FPBO instantiation count
    // on the color-space axis. Gates BOTH the trust-region clip on the
    // per-channel Adam delta AND the linear-to-sRGB Jacobian inversion on
    // the per-channel gradient.
    bool color_trust_linear,
    int BLOCK_SIZE,
    int QUANT_BITS = 8,    // SH-Adam quant bit depth: 4 or 8 (ignored when BLOCK_SIZE == 0)
    int VALUE_BITS = 32    // SH-VALUE quant bit depth: 32 (off), 8, or 16. Per-splat-block
                            // bounds layout (one float2 per 256 splats) when != 32; only
                            // wired through the BLOCK_SIZE != 0 path for now -- combining
                            // VALUE_BITS != 32 with BLOCK_SIZE == 0 (fp32 optim state)
                            // is deferred to a follow-up; the kernel falls back to
                            // fp32 SH reads/writes in that combination.
>
#if 1
__global__ void __launch_bounds__(512) fused_projection_bwd_optimizer_3dgs_kernel
#else
__global__ void fused_projection_bwd_optimizer_3dgs_kernel
#endif
(
    // fwd inputs
    const uint32_t C,
    const uint32_t N,
    const uint32_t num_sh_buffer,  // stride into sh_packed; may exceed template num_sh during SH warmup
    typename SplatPrimitive::WorldBuffer splats_world,
    const float *__restrict__ viewmats, // [C, 4, 4]
    const float4 *__restrict__ intrins,  // [C, 4], fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const uint32_t image_width,
    const uint32_t image_height,
    // fwd outputs
    const int32_t *__restrict__ camera_id_bounds,   // [N+1]
    const int32_t *__restrict__ camera_ids,   // [nnz] -- ORIGINAL (unsorted) order
    const int32_t *__restrict__ perm,         // [nnz] -- sorted_pos -> original_pos
    const float4 *__restrict__ aabb,   // [C, N, 4] or [nnz, 4] (original order)
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
    const uint8_t* __restrict__ sh_packed,      // AoS (u, sqrt_g2) packed SH state
    float4* __restrict__ sh_quant_bounds,
    // SH VALUE quant (per-splat-block layout: 1 float2 per 256 splats covering
    // all 3*K cells per splat). Active when VALUE_BITS != 32. Read+written
    // in place by the SH section below. The byte stride is the BUFFER's full
    // SH stride (3 * num_sh_buffer cells per splat), same as sh_packed.
    const uint8_t* __restrict__ sh_value_packed,
    float2* __restrict__ sh_value_bounds,
    // float *__restrict__ v_viewmats // [C, 4, 4] optional
    // optimizer params
    const float* __restrict__ radii,
    const float lr_means,
    const float lr_quats,
    const float lr_scales,
    const float lr_opacs,
    const float lr_features_dc,
    const float lr_features_sh,
    const float max_gauss_ratio,
    const float scale_regularization_weight,
    const float mcmc_opacity_reg_weight,
    const float mcmc_scale_reg_weight,
    const float erank_reg_weight,
    const float erank_reg_weight_s3,
    const float quat_norm_reg_weight,
    const float sh_reg_weight,
    const float eps_tr,
    const int32_t scalar_step,
    const int32_t* __restrict__ steps
) {
    uint32_t gid = blockIdx.x * (BLOCK_SIZE == 0 ? blockDim.x : BLOCK_SIZE) + threadIdx.x;
    bool inside = (gid < N);

    bool packed = (camera_id_bounds != nullptr && camera_ids != nullptr);
    int cid_0 = 0, cid_1 = C;
    if (packed) {
        cid_0 = inside ? camera_id_bounds[gid] : 0;
        cid_1 = inside ? camera_id_bounds[gid+1] : 0;
    }

    // Load splat
    typename SplatPrimitive::World splat_world;
    if (inside) splat_world.load(splats_world, gid);

    constexpr int num_sh = SplatPrimitive::World::num_sh();
    float v_sh_coeffs[num_sh == 0 ? 1 : 3*num_sh];

    typename SplatPrimitive::World v_splat_world = SplatPrimitive::World::zero();
    if (inside) v_splat_world.atomicLoad(v_splats_world, gid);
    v_splat_world.features_sh = &v_sh_coeffs[0];
    #pragma unroll
    for (int i = 0; i < 3*num_sh; ++i)
        v_splat_world.features_sh[i] = 0.0f;

    float3x3 v_R = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
    float3 v_t = {0.f, 0.f, 0.f};
    float3 vr_world_pos = {0.f, 0.f, 0.f};
    float3 h_world_pos = {0.f, 0.f, 0.f};
    if constexpr (hessian_diagonal_output_mode == HessianDiagonalOutputMode::Position) {
        if (inside) vr_world_pos = vr_splats_world.means(gid);
        if (inside) h_world_pos = h_splats_world.means(gid);
    }
    typename SplatPrimitive::World vr_splat_world = SplatPrimitive::World::zero();
    typename SplatPrimitive::World h_splat_world = SplatPrimitive::World::zero();
    if constexpr (hessian_diagonal_output_mode == HessianDiagonalOutputMode::AllReasonable) {
        if (inside) vr_splat_world.atomicLoad(vr_splat_world, gid);
        if (inside) h_splat_world.atomicLoad(h_splat_world, gid);
    }

    // Loop over intersections.
    // In packed mode `cid_t` is a SORTED position (sorted by gaussian_id),
    // but `aabb` and `v_splats_screen` are stored in the ORIGINAL out_idx
    // order from the forward's intersection-mask scan -- and `camera_ids`
    // is also the ORIGINAL (unsorted) buffer here. `perm[cid_t]` recovers
    // the original out_idx so all three reads land on the right entry.
    for (int cid_t = cid_0; cid_t < cid_1; ++cid_t) {
        int idx = packed ? perm[cid_t] : cid_t * N + gid;
        int cid = (packed ? camera_ids[idx] : cid_t);
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

        // Accumulate gradient. When VALUE_BITS != 32 the canonical SH storage
        // is the packed codec buffer, NOT the fp32 features_sh array; pass the
        // codec args through so project_vjp evaluates v_dir against the actual
        // (decoded) SH coefficients. Without this, v_dir reads the stale fp32
        // features_sh (left untouched by value-quant FPBO writeback) which is
        // typically all zero -> biased v_means / v_R / v_t every step.
        if constexpr (hessian_diagonal_output_mode == HessianDiagonalOutputMode::None) {
            if constexpr (VALUE_BITS == 32) {
                splat_world.template project_vjp<camera_model, false, 32>(cam, v_splat_screen, v_splat_world, v_R, v_t);
            } else {
                const int64_t sh_base_vjp = (int64_t)3 * (int64_t)num_sh_buffer * (int64_t)gid;
                const int64_t sh_bounds_stride_vjp = (int64_t)256 * 3 * (int64_t)num_sh_buffer;
                splat_world.template project_vjp<camera_model, false, VALUE_BITS>(
                    cam, v_splat_screen, v_splat_world, v_R, v_t,
                    const_cast<uint8_t*>(sh_value_packed), sh_value_bounds,
                    sh_base_vjp, sh_bounds_stride_vjp);
            }
        } else if constexpr (hessian_diagonal_output_mode == HessianDiagonalOutputMode::Position) {
            splat_world.template project_vjp_h_pos<camera_model, false>(cam, v_splat_screen,
                vr_splat_screen, h_splat_screen, v_splat_world, v_R, v_t, vr_world_pos, h_world_pos);
        } else if constexpr (hessian_diagonal_output_mode == HessianDiagonalOutputMode::AllReasonable) {
            splat_world.template project_vjp_h_all<camera_model, false>(cam, v_splat_screen,
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
    v_splat_world.features_dc += sh_reg_weight * (
        fmaxf(splat_world.features_dc - make_float3(0.5f / 0.28209479177387814f), 0.0f) +
        fminf(splat_world.features_dc + make_float3(0.5f / 0.28209479177387814f), 0.0f)
    );

    // optimizer for mean/quat/scale

    static constexpr float eps = 1e-15f;
    static constexpr float beta1 = 0.9f;
    static constexpr float beta2 = 0.999f;
    float step = inside ? (float)(steps != nullptr ? steps[gid] : scalar_step) : 0;
    float inv_bias_correction1 = 1.0f / (1.0f - powf(beta1, step));
    float inv_bias_correction2 = 1.0f / (1.0f - powf(beta2, step));

    // update scales
    if (inside) {
        float3 g1_scale = beta1 * g1_splats_world.scales(gid) + (1.f - beta1) * v_splat_world.scale;
        float3 g2_scale = beta2 * g2_splats_world.scales(gid) + (1.f - beta2) * v_splat_world.scale*v_splat_world.scale;
        float3 new_scale = splat_world.scale - lr_scales * inv_bias_correction1 * g1_scale / (sqrtf(g2_scale * inv_bias_correction2) + eps);
        splats_world.scales(gid) = new_scale;
        g1_splats_world.scales(gid) = g1_scale;
        g2_splats_world.scales(gid) = g2_scale;
    }

    // update quats (Riemannian)
    if (inside) {
        float4 new_quat = normalize(splat_world.quat);
        v_splat_world.quat -= dot(new_quat, v_splat_world.quat) * new_quat;
        float4 g1_quat = beta1 * g1_splats_world.quats(gid) + (1.f - beta1) * v_splat_world.quat;
        float4 g2_quat = beta2 * g2_splats_world.quats(gid) + (1.f - beta2) * v_splat_world.quat*v_splat_world.quat;
        new_quat -= lr_quats * inv_bias_correction1 * g1_quat / (sqrtf(g2_quat * inv_bias_correction2) + eps);
        splats_world.quats(gid) = normalize(new_quat);
        g1_splats_world.quats(gid) = g1_quat;
        g2_splats_world.quats(gid) = g2_quat;
    }

    // update opacs
    if (inside) {
        float g1_opac = beta1 * g1_splats_world.opacities(gid) + (1.f - beta1) * v_splat_world.opacity;
        float g2_opac = beta2 * g2_splats_world.opacities(gid) + (1.f - beta2) * v_splat_world.opacity*v_splat_world.opacity;
        float new_opac = splat_world.opacity - lr_opacs * inv_bias_correction1 * g1_opac / (sqrtf(g2_opac * inv_bias_correction2) + eps);
        splats_world.opacities(gid) = new_opac;
        g1_splats_world.opacities(gid) = g1_opac;
        g2_splats_world.opacities(gid) = g2_opac;
    }

    // update means (scale agnostic)
    if (inside) {
        float3 g1_mean = g1_splats_world.means(gid);
        float3 g2_mean = g2_splats_world.means(gid);
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
        } else {
            g1_mean = beta1 * g1_mean + (1.f - beta1) * v_splat_world.mean;  // unit: L^-1
            g2_mean = beta2 * g2_mean + (1.f - beta2) * v_splat_world.mean*v_splat_world.mean;  // unit: L^-2
        }
        splats_world.means(gid) = splat_world.mean - lr_means * inv_bias_correction1 * g1_mean /
            (sqrtf(g2_mean * inv_bias_correction2) + eps); // unit: L or dimensionless for dimensionless lr_means
        g1_splats_world.means(gid) = g1_mean;
        g2_splats_world.means(gid) = g2_mean;
    }

    // update features_dc
    if (inside) {
        float3 v_dc = v_splat_world.features_dc;
        // Linear -> sRGB Jacobian inversion on the gradient: divides the
        // working-color-space grad by d(sRGB)/d(linear) so the Adam update
        // is computed in the linear domain (matches fused_adamtr_linear_rgb).
        if constexpr (color_trust_linear) {
            v_dc.x /= SlangPixelWise::linear_rgb_to_srgb_grad(kSh0 * splat_world.features_dc.x + 0.5f);
            v_dc.y /= SlangPixelWise::linear_rgb_to_srgb_grad(kSh0 * splat_world.features_dc.y + 0.5f);
            v_dc.z /= SlangPixelWise::linear_rgb_to_srgb_grad(kSh0 * splat_world.features_dc.z + 0.5f);
        }
        float3 g1_feature_dc = beta1 * g1_splats_world.features_dc(gid) + (1.f - beta1) * v_dc;
        float3 g2_feature_dc = beta2 * g2_splats_world.features_dc(gid) + (1.f - beta2) * v_dc*v_dc;
        float3 denom = sqrtf(g2_feature_dc * inv_bias_correction2) + make_float3(eps);
        float3 delta = make_float3(-lr_features_dc * inv_bias_correction1) * (g1_feature_dc / denom);

        if constexpr (color_trust_linear) {
            // Trust-region clip in the working color space (matches
            // fused_adamtr_rgb_optim_kernel). delta_max = kSh0*sqrt(4*eps_tr*c/opac).
            float opac = sigmoid(splat_world.opacity);
            float3 c = fmaxf(kSh0 * splat_world.features_dc + 0.5f, (1.0f/255.0f)*(1.0f/255.0f));
            float inv_opac = 1.0f / fmaxf(opac, 1e-12f);
            float3 clip = kSh0 * sqrtf(make_float3(
                4.0f * eps_tr * c.x * inv_opac,
                4.0f * eps_tr * c.y * inv_opac,
                4.0f * eps_tr * c.z * inv_opac));
            delta.x = fminf(fmaxf(delta.x, -clip.x), clip.x);
            delta.y = fminf(fmaxf(delta.y, -clip.y), clip.y);
            delta.z = fminf(fmaxf(delta.z, -clip.z), clip.z);
            delta.x = isfinite(delta.x) ? delta.x : 0.0f;
            delta.y = isfinite(delta.y) ? delta.y : 0.0f;
            delta.z = isfinite(delta.z) ? delta.z : 0.0f;
        }
        splats_world.features_dc(gid) = splat_world.features_dc + delta;
        g1_splats_world.features_dc(gid) = g1_feature_dc;
        g2_splats_world.features_dc(gid) = g2_feature_dc;
    }

    // update features_sh
    // TODO: more cache-friendly memory access pattern
    int i = 0;
    float lr_sh = lr_features_sh * inv_bias_correction1;
    float* v_sh_ptr = v_splat_world.features_sh;
    float* sh_ptr = splats_world.features_sh(gid);
    const float reg_weight = sh_reg_weight * (1.0f / (float)num_sh);

    // Per-channel DC color + clip radius for trust-region SH update (mirrors
    // fused_adamtr_rgb_sh_optim_kernel). c_dc is unclamped (used for the
    // linear->sRGB Jacobian inversion); c_dc_clamped is the clamped variant
    // (used for the clip radius). Both are computed once per splat from DC.
    float3 c_dc = make_float3(0.0f), c_dc_clamped = make_float3(0.0f);
    float3 clip_dc = make_float3(0.0f);
    if constexpr (color_trust_linear || color_trust_linear) {
        if (inside) {
            c_dc = kSh0 * splat_world.features_dc + 0.5f;
            c_dc_clamped = fmaxf(c_dc, (1.0f/255.0f)*(1.0f/255.0f));
        }
    }
    if constexpr (color_trust_linear) {
        if (inside) {
            float opac = sigmoid(splat_world.opacity);
            float inv_opac = 1.0f / fmaxf(opac, 1e-12f);
            clip_dc = kSh0 * sqrtf(make_float3(
                4.0f * eps_tr * c_dc_clamped.x * inv_opac,
                4.0f * eps_tr * c_dc_clamped.y * inv_opac,
                4.0f * eps_tr * c_dc_clamped.z * inv_opac));
        }
    }

    if constexpr (BLOCK_SIZE == 0) {
        if (!inside) return;
        float* g1_sh_ptr = g1_splats_world.features_sh(gid);
        float* g2_sh_ptr = g2_splats_world.features_sh(gid);
        // Skip the float4 fast path when color-space trust region is on
        // (per-channel clip / linear correction breaks the lane-coupled fast
        // path). The scalar loop below handles every i.
        if constexpr (!color_trust_linear && !color_trust_linear) {
            #pragma unroll
            for (; i+3 < 3*num_sh; i += 4) {
                float4 sh_coeff = float4{sh_ptr[i], sh_ptr[i+1], sh_ptr[i+2], sh_ptr[i+3]};
                float4 v_sh_coeff = float4{v_sh_ptr[i], v_sh_ptr[i+1], v_sh_ptr[i+2], v_sh_ptr[i+3]} + reg_weight * sh_coeff;
                float4 g1_feature_sh = beta1 * float4{g1_sh_ptr[i], g1_sh_ptr[i+1], g1_sh_ptr[i+2], g1_sh_ptr[i+3]} + (1.f - beta1) * v_sh_coeff;
                float4 g2_feature_sh = beta2 * float4{g2_sh_ptr[i], g2_sh_ptr[i+1], g2_sh_ptr[i+2], g2_sh_ptr[i+3]} + (1.f - beta2) * v_sh_coeff*v_sh_coeff;
                sh_coeff -= lr_sh * g1_feature_sh / (sqrtf(g2_feature_sh * inv_bias_correction2) + eps);
                sh_ptr[i] = sh_coeff.x; sh_ptr[i+1] = sh_coeff.y; sh_ptr[i+2] = sh_coeff.z; sh_ptr[i+3] = sh_coeff.w;
                g1_sh_ptr[i] = g1_feature_sh.x; g1_sh_ptr[i+1] = g1_feature_sh.y; g1_sh_ptr[i+2] = g1_feature_sh.z; g1_sh_ptr[i+3] = g1_feature_sh.w;
                g2_sh_ptr[i] = g2_feature_sh.x; g2_sh_ptr[i+1] = g2_feature_sh.y; g2_sh_ptr[i+2] = g2_feature_sh.z; g2_sh_ptr[i+3] = g2_feature_sh.w;
            }
        }
        #pragma unroll
        for (; i < 3*num_sh; ++i) {
            float sh_coeff = sh_ptr[i];
            float v_sh_coeff = v_sh_ptr[i] + reg_weight * sh_coeff;
            // Channel index 0/1/2 within (band, channel) layout.
            const int ch = i % 3;
            const float c_ch = (ch == 0) ? c_dc.x : (ch == 1) ? c_dc.y : c_dc.z;
            if constexpr (color_trust_linear) {
                v_sh_coeff /= SlangPixelWise::linear_rgb_to_srgb_grad(c_ch);
            }
            float g1_feature_sh = beta1 * g1_sh_ptr[i] + (1.f - beta1) * v_sh_coeff;
            float g2_feature_sh = beta2 * g2_sh_ptr[i] + (1.f - beta2) * v_sh_coeff*v_sh_coeff;
            float delta = -lr_sh * g1_feature_sh / (sqrtf(g2_feature_sh * inv_bias_correction2) + eps);
            if constexpr (color_trust_linear) {
                float clip = (ch == 0) ? clip_dc.x : (ch == 1) ? clip_dc.y : clip_dc.z;
                delta = fminf(fmaxf(delta, -clip), clip);
                delta = isfinite(delta) ? delta : 0.0f;
            }
            sh_ptr[i] = sh_coeff + delta;
            g1_sh_ptr[i] = g1_feature_sh;
            g2_sh_ptr[i] = g2_feature_sh;
        }
    } else if constexpr (num_sh > 0) {
        // Joint (u, sqrt(g2)) AoS quantization for SH momentum, via the shared
        // QuantizedAdamState codec. Each cell = one (splat, coef, channel)
        // triple holds 2 bytes (u, sqrt_g2). The codec dequantizes and rebuilds
        // ordinary (g1, g2) so the Adam math below is unchanged.
        using SHQState = QuantizedAdamState<QUANT_BITS, BLOCK_SIZE>;
        // Index into the packed buffer must use the BUFFER's full SH stride
        // (engine().num_sh = model max). The template `num_sh` here equals
        // SplatPrimitive::num_sh(), which during SH warmup is capped at the
        // runtime sh_degree_to_use and is strictly less than num_sh_buffer --
        // using it as the stride would corrupt neighboring splats' bytes.
        const int64_t sh_base = (int64_t)3 * (int64_t)num_sh_buffer * gid;
        uint8_t* sh_packed_rw = const_cast<uint8_t*>(sh_packed);
        float4 quant_bounds = sh_quant_bounds[blockIdx.x];
        float4 mm = make_float4(1e30f, -1e30f, 1e30f, -1e30f);

        // SH VALUE-quant: per-splat-block layout, one float2 bound per block of
        // 256 splats covering all 3*K cells/splat. Mirrors the optim-state
        // reduction below: read old bound at start, accumulate per-thread
        // min/max during update, block-reduce, encode against new bound.
        using SHValReader8  = QuantizedTensor<8,  BLOCK_SIZE>;
        using SHValReader16 = QuantizedTensor<16, BLOCK_SIZE>;
        uint8_t* sh_value_packed_rw = const_cast<uint8_t*>(sh_value_packed);
        float2 sh_value_old_mm = (VALUE_BITS != 32 && sh_value_bounds != nullptr)
            ? sh_value_bounds[blockIdx.x] : float2{0.f, 0.f};
        float2 sh_value_new_mm = (VALUE_BITS != 32)
            ? float2{1e30f, -1e30f} : float2{0.f, 0.f};
        // Hold per-thread updated SH values until after the value-bounds
        // block reduction, so we encode against the FRESH bounds (mirrors
        // optim-state pattern -- mm computed across all 256 splats, then
        // each thread encodes its 3*K cells).
        float3 sh_updated_vals[num_sh];

        float3 g1_sh_vals[num_sh];
        float3 g2_sh_vals[num_sh];

    if (inside) {
        #pragma unroll
        for (int j = 0; j < num_sh; ++j) {
            // Decode each of the 3 channels via the codec.
            float2 g1g2_x = SHQState::decode_g1g2(sh_packed, sh_base + 3*j + 0, quant_bounds);
            float2 g1g2_y = SHQState::decode_g1g2(sh_packed, sh_base + 3*j + 1, quant_bounds);
            float2 g1g2_z = SHQState::decode_g1g2(sh_packed, sh_base + 3*j + 2, quant_bounds);
            float3 g1_feature_sh = make_float3(g1g2_x.x, g1g2_y.x, g1g2_z.x);
            float3 g2_feature_sh = make_float3(g1g2_x.y, g1g2_y.y, g1g2_z.y);

            // Read source SH value: fp32 directly OR codec-decode from packed
            // bytes using the per-splat-block bound.
            float3 sh_coeff;
            if constexpr (VALUE_BITS == 32) {
                sh_coeff = make_float3(sh_ptr[3*j], sh_ptr[3*j+1], sh_ptr[3*j+2]);
            } else if constexpr (VALUE_BITS == 8) {
                sh_coeff.x = SHValReader8::decode_v(sh_value_packed, sh_base + 3*j + 0, sh_value_old_mm);
                sh_coeff.y = SHValReader8::decode_v(sh_value_packed, sh_base + 3*j + 1, sh_value_old_mm);
                sh_coeff.z = SHValReader8::decode_v(sh_value_packed, sh_base + 3*j + 2, sh_value_old_mm);
            } else /* VALUE_BITS == 16 */ {
                sh_coeff.x = SHValReader16::decode_v(sh_value_packed, sh_base + 3*j + 0, sh_value_old_mm);
                sh_coeff.y = SHValReader16::decode_v(sh_value_packed, sh_base + 3*j + 1, sh_value_old_mm);
                sh_coeff.z = SHValReader16::decode_v(sh_value_packed, sh_base + 3*j + 2, sh_value_old_mm);
            }
            float3 v_sh_coeff = make_float3(
                v_sh_ptr[3*j] + reg_weight * sh_coeff.x,
                v_sh_ptr[3*j+1] + reg_weight * sh_coeff.y,
                v_sh_ptr[3*j+2] + reg_weight * sh_coeff.z
            );
            // Linear -> sRGB Jacobian inversion on the per-channel grad
            // (matches fused_adamtr_linear_rgb_sh_optim_kernel). c_dc is the
            // DC color in [0, 1]; the divisor is d(sRGB)/d(linear) at that c.
            if constexpr (color_trust_linear) {
                v_sh_coeff.x /= SlangPixelWise::linear_rgb_to_srgb_grad(c_dc.x);
                v_sh_coeff.y /= SlangPixelWise::linear_rgb_to_srgb_grad(c_dc.y);
                v_sh_coeff.z /= SlangPixelWise::linear_rgb_to_srgb_grad(c_dc.z);
            }
            float3 g1_updated = beta1 * g1_feature_sh + (1.f - beta1) * v_sh_coeff;
            float3 g2_updated = beta2 * g2_feature_sh + (1.f - beta2) * make_float3(v_sh_coeff.x*v_sh_coeff.x, v_sh_coeff.y*v_sh_coeff.y, v_sh_coeff.z*v_sh_coeff.z);
            float3 denom = sqrtf(g2_updated * inv_bias_correction2) + make_float3(eps);
            float3 sh_delta = make_float3(
                -lr_sh * g1_updated.x / denom.x,
                -lr_sh * g1_updated.y / denom.y,
                -lr_sh * g1_updated.z / denom.z
            );
            if constexpr (color_trust_linear) {
                // Per-channel trust-region clip using the DC's clip radius.
                sh_delta.x = fminf(fmaxf(sh_delta.x, -clip_dc.x), clip_dc.x);
                sh_delta.y = fminf(fmaxf(sh_delta.y, -clip_dc.y), clip_dc.y);
                sh_delta.z = fminf(fmaxf(sh_delta.z, -clip_dc.z), clip_dc.z);
                sh_delta.x = isfinite(sh_delta.x) ? sh_delta.x : 0.0f;
                sh_delta.y = isfinite(sh_delta.y) ? sh_delta.y : 0.0f;
                sh_delta.z = isfinite(sh_delta.z) ? sh_delta.z : 0.0f;
            }
            float3 sh_updated = sh_coeff + sh_delta;
            if constexpr (VALUE_BITS == 32) {
                // fp32 path: write straight to global storage now (no later
                // encode step is needed).
                sh_ptr[3*j] = sh_updated.x;
                sh_ptr[3*j+1] = sh_updated.y;
                sh_ptr[3*j+2] = sh_updated.z;
            } else {
                // Defer encode until after the block-wide value-bounds reduce
                // below -- the bounds need every thread's contribution before
                // we know the encode range. Stash the new values and update
                // this thread's running min/max.
                sh_updated_vals[j] = sh_updated;
                sh_value_new_mm.x = fminf(fminf(fminf(sh_value_new_mm.x,
                    sh_updated.x), sh_updated.y), sh_updated.z);
                sh_value_new_mm.y = fmaxf(fmaxf(fmaxf(sh_value_new_mm.y,
                    sh_updated.x), sh_updated.y), sh_updated.z);
            }

            g1_sh_vals[j] = g1_updated;
            g2_sh_vals[j] = g2_updated;

            // Reduce min/max in the (u, sqrt(g2)) basis (codec helper).
            float2 us_x = SHQState::g1g2_to_us(g1_updated.x, g2_updated.x);
            float2 us_y = SHQState::g1g2_to_us(g1_updated.y, g2_updated.y);
            float2 us_z = SHQState::g1g2_to_us(g1_updated.z, g2_updated.z);
            mm.x = fminf(fminf(fminf(mm.x, us_x.x), us_y.x), us_z.x);
            mm.y = fmaxf(fmaxf(fmaxf(mm.y, us_x.x), us_y.x), us_z.x);
            mm.z = fminf(fminf(fminf(mm.z, us_x.y), us_y.y), us_z.y);
            mm.w = fmaxf(fmaxf(fmaxf(mm.w, us_x.y), us_y.y), us_z.y);
        }
        // SH-degree warmup: when num_sh_buffer > num_sh, the inactive cells
        // (j in [num_sh, num_sh_buffer)) are zero-encoded after the reduce
        // below so they decode to ~0 once they become active in a later
        // step. That zero-encode is endpoint-safe only if BOTH bounds cover
        // 0 -- expand the per-thread bound here so the block reduction
        // produces a bound that includes 0. Done only in inside threads so
        // outside-thread sentinels pass through the reduce unchanged.
        if (num_sh_buffer > (uint32_t)num_sh) {
            if constexpr (VALUE_BITS != 32) {
                sh_value_new_mm.x = fminf(sh_value_new_mm.x, 0.0f);
                sh_value_new_mm.y = fmaxf(sh_value_new_mm.y, 0.0f);
            }
            // (u=0, log_s=0) is the codec's (g1=0, g2=0) fixed point.
            // log_s is non-negative, so forcing mm.z <= 0 is the only side
            // that needs explicit clamping; mm.w already >= 0 from reduce.
            mm.x = fminf(mm.x, 0.0f);
            mm.y = fmaxf(mm.y, 0.0f);
            mm.z = fminf(mm.z, 0.0f);
        }
    }  // inside

        auto warp = cg::tiled_partition<WARP_SIZE>(cg::this_thread_block());
        mm.x = cg::reduce(warp, mm.x, cg::less<float>());
        mm.y = cg::reduce(warp, mm.y, cg::greater<float>());
        mm.z = cg::reduce(warp, mm.z, cg::less<float>());
        mm.w = cg::reduce(warp, mm.w, cg::greater<float>());
        __shared__ float4 shared_reduce[BLOCK_SIZE / WARP_SIZE];
        if (threadIdx.x % WARP_SIZE == 0)
            shared_reduce[threadIdx.x / WARP_SIZE] = mm;
        __syncthreads();
        mm = (threadIdx.x < BLOCK_SIZE / WARP_SIZE) ?
            shared_reduce[threadIdx.x] : make_float4(1e30f, -1e30f, 1e30f, -1e30f);
        mm.x = cg::reduce(warp, mm.x, cg::less<float>());
        mm.y = cg::reduce(warp, mm.y, cg::greater<float>());
        mm.z = cg::reduce(warp, mm.z, cg::less<float>());
        mm.w = cg::reduce(warp, mm.w, cg::greater<float>());
        __syncthreads();
        if (threadIdx.x < BLOCK_SIZE / WARP_SIZE)
            shared_reduce[threadIdx.x] = mm;
        __syncthreads();
        mm = shared_reduce[threadIdx.x / WARP_SIZE];

        if (threadIdx.x == 0)
            sh_quant_bounds[blockIdx.x] = mm;

        // ---- SH VALUE bounds block reduction + encode (VALUE_BITS != 32) ----
        // Mirrors the optim-state reduction above but on a float2 (min, max).
        // Must run BEFORE the early-return for !inside threads -- the bounds
        // reduction needs every thread to participate (out-of-range threads
        // pass sentinel values that don't affect min/max).
        if constexpr (VALUE_BITS != 32) {
            if (!inside) sh_value_new_mm = float2{1e30f, -1e30f};
            auto warp_v = cg::tiled_partition<WARP_SIZE>(cg::this_thread_block());
            sh_value_new_mm.x = cg::reduce(warp_v, sh_value_new_mm.x, cg::less<float>());
            sh_value_new_mm.y = cg::reduce(warp_v, sh_value_new_mm.y, cg::greater<float>());
            __shared__ float2 shared_reduce_v[BLOCK_SIZE / WARP_SIZE];
            if (threadIdx.x % WARP_SIZE == 0)
                shared_reduce_v[threadIdx.x / WARP_SIZE] = sh_value_new_mm;
            __syncthreads();
            sh_value_new_mm = (threadIdx.x < BLOCK_SIZE / WARP_SIZE) ?
                shared_reduce_v[threadIdx.x] : float2{1e30f, -1e30f};
            sh_value_new_mm.x = cg::reduce(warp_v, sh_value_new_mm.x, cg::less<float>());
            sh_value_new_mm.y = cg::reduce(warp_v, sh_value_new_mm.y, cg::greater<float>());
            __syncthreads();
            if (threadIdx.x < BLOCK_SIZE / WARP_SIZE)
                shared_reduce_v[threadIdx.x] = sh_value_new_mm;
            __syncthreads();
            sh_value_new_mm = shared_reduce_v[threadIdx.x / WARP_SIZE];

            if (threadIdx.x == 0)
                sh_value_bounds[blockIdx.x] = sh_value_new_mm;
        }

        if (!inside)
            return;

        // Re-encode the new Adam state via the codec.
        #pragma unroll
        for (int j = 0; j < num_sh; ++j) {
            float3 g1_upd = g1_sh_vals[j];
            float3 g2_upd = g2_sh_vals[j];
            SHQState::encode_g1g2(sh_packed_rw, sh_base + 3*j + 0, g1_upd.x, g2_upd.x, mm);
            SHQState::encode_g1g2(sh_packed_rw, sh_base + 3*j + 1, g1_upd.y, g2_upd.y, mm);
            SHQState::encode_g1g2(sh_packed_rw, sh_base + 3*j + 2, g1_upd.z, g2_upd.z, mm);
        }

        // Re-encode the new SH VALUES via the linear codec against the fresh
        // per-splat-block bound. Only the inside-block (gated above) runs.
        if constexpr (VALUE_BITS == 8) {
            #pragma unroll
            for (int j = 0; j < num_sh; ++j) {
                float3 v = sh_updated_vals[j];
                SHValReader8::encode_v(sh_value_packed_rw, sh_base + 3*j + 0, v.x, sh_value_new_mm);
                SHValReader8::encode_v(sh_value_packed_rw, sh_base + 3*j + 1, v.y, sh_value_new_mm);
                SHValReader8::encode_v(sh_value_packed_rw, sh_base + 3*j + 2, v.z, sh_value_new_mm);
            }
        } else if constexpr (VALUE_BITS == 16) {
            #pragma unroll
            for (int j = 0; j < num_sh; ++j) {
                float3 v = sh_updated_vals[j];
                SHValReader16::encode_v(sh_value_packed_rw, sh_base + 3*j + 0, v.x, sh_value_new_mm);
                SHValReader16::encode_v(sh_value_packed_rw, sh_base + 3*j + 1, v.y, sh_value_new_mm);
                SHValReader16::encode_v(sh_value_packed_rw, sh_base + 3*j + 2, v.z, sh_value_new_mm);
            }
        }

        // SH-degree warmup: zero-encode the inactive cells (j in
        // [num_sh, num_sh_buffer)) so that when the kernel template's
        // num_sh increases at a later warmup step the newly-active bands
        // decode to ~0 instead of the old bound's mm.x (which biases the
        // first few Adam steps after each degree boundary). The bound
        // expansion above guarantees the (0, 0) point is in range, so
        // these encodes land near the codec's true zero fixed point.
        // Runtime-gated -- becomes a no-op once num_sh == num_sh_buffer.
        if (num_sh_buffer > (uint32_t)num_sh) {
            // Adam state: encode (g1=0, g2=0).
            #pragma unroll 1
            for (int j = num_sh; j < (int)num_sh_buffer; ++j) {
                SHQState::encode_g1g2(sh_packed_rw, sh_base + 3*j + 0, 0.0f, 0.0f, mm);
                SHQState::encode_g1g2(sh_packed_rw, sh_base + 3*j + 1, 0.0f, 0.0f, mm);
                SHQState::encode_g1g2(sh_packed_rw, sh_base + 3*j + 2, 0.0f, 0.0f, mm);
            }
            // SH values: encode 0 against the fresh, 0-inclusive bound.
            if constexpr (VALUE_BITS == 8) {
                #pragma unroll 1
                for (int j = num_sh; j < (int)num_sh_buffer; ++j) {
                    SHValReader8::encode_v(sh_value_packed_rw, sh_base + 3*j + 0, 0.0f, sh_value_new_mm);
                    SHValReader8::encode_v(sh_value_packed_rw, sh_base + 3*j + 1, 0.0f, sh_value_new_mm);
                    SHValReader8::encode_v(sh_value_packed_rw, sh_base + 3*j + 2, 0.0f, sh_value_new_mm);
                }
            } else if constexpr (VALUE_BITS == 16) {
                #pragma unroll 1
                for (int j = num_sh; j < (int)num_sh_buffer; ++j) {
                    SHValReader16::encode_v(sh_value_packed_rw, sh_base + 3*j + 0, 0.0f, sh_value_new_mm);
                    SHValReader16::encode_v(sh_value_packed_rw, sh_base + 3*j + 1, 0.0f, sh_value_new_mm);
                    SHValReader16::encode_v(sh_value_packed_rw, sh_base + 3*j + 2, 0.0f, sh_value_new_mm);
                }
            }
        }
    }
}

template<
    typename SplatPrimitive,
    CameraModelType camera_model,
    HessianDiagonalOutputMode hessian_diagonal_output_mode,
    const bool use_scale_agnostic_mean,
    const bool color_trust_linear,
    // SH quantization level. Single int collapses the prior (QUANT_BITS, VALUE_BITS)
    // axes into 2 combos:
    //   LEVEL == 0 : off          -- 32-bit value, 32-bit (fp32) optim state
    //   LEVEL == 1 : light        -- 16-bit value, 8x2-bit optim (2 B / cell)
    // The wrapper derives BLOCK_SIZE / QUANT_BITS / VALUE_BITS internally from
    // LEVEL via constexpr so the kernel template signature is unchanged.
    const int  LEVEL
>
void fused_projection_bwd_optimizer_3dgs_kernel_wrapper(
    cudaStream_t stream,
    // fwd inputs
    const uint32_t C,
    const uint32_t N,
    const uint32_t num_sh_buffer,
    typename SplatPrimitive::WorldBuffer splats_world,
    const float *__restrict__ viewmats, // [C, 4, 4]
    const float4 *__restrict__ intrins,  // [C, 4], fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const uint32_t image_width,
    const uint32_t image_height,
    // fwd outputs
    const int32_t *__restrict__ camera_id_bounds,   // [N+1]
    const int32_t *__restrict__ camera_ids,   // [nnz] -- ORIGINAL (unsorted) order
    const int32_t *__restrict__ perm,         // [nnz] -- sorted_pos -> original_pos
    const float4 *__restrict__ aabb,   // [C, N, 4] or [nnz, 4] (original order)
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
    const uint8_t* __restrict__ sh_packed,      // AoS (u, sqrt_g2) packed SH state
    float4* __restrict__ sh_quant_bounds,
    const uint8_t* __restrict__ sh_value_packed,
    float2* __restrict__ sh_value_bounds,
    // float *__restrict__ v_viewmats // [C, 4, 4] optional
    // optimizer params
    const float* __restrict__ radii,
    const float lr_means,
    const float lr_quats,
    const float lr_scales,
    const float lr_opacs,
    const float lr_features_dc,
    const float lr_features_sh,
    const float max_gauss_ratio,
    const float scale_regularization_weight,
    const float mcmc_opacity_reg_weight,
    const float mcmc_scale_reg_weight,
    const float erank_reg_weight,
    const float erank_reg_weight_s3,
    const float quat_norm_reg_weight,
    const float sh_reg_weight,
    const float eps_tr,
    const int32_t scalar_step,
    const int32_t* __restrict__ steps
) {
    // One kernel instantiation per wrapper instantiation -- the runtime
    // dispatch on the SH quantization level has been lifted into the
    // host-side dispatcher in FusedProjectionBwdOptim.cu so each translation
    // unit compiles a single template specialization.
    constexpr int BLOCK_SIZE_LAUNCH = 256;
    // Derive optim-state and value-quant bit widths from LEVEL:
    //   LEVEL 0: QUANT_BITS=0 (off), VALUE_BITS=32 (off)
    //   LEVEL 1: QUANT_BITS=8       VALUE_BITS=16
    // QUANT_BITS == 0 selects the kernel's no-quant path via its BLOCK_SIZE
    // template (0 = fp32 g1/g2); KERNEL_QUANT_BITS is passed through but
    // unused in that branch (it just needs to satisfy QuantizedAdamState's
    // static_assert).
    constexpr int LEVEL_QUANT_BITS = (LEVEL == 0) ? 0 : 8;
    constexpr int LEVEL_VALUE_BITS = (LEVEL == 0) ? 32 : 16;
    constexpr int KERNEL_BLOCK_SIZE = (LEVEL_QUANT_BITS == 0) ? 0 : BLOCK_SIZE_LAUNCH;
    constexpr int KERNEL_QUANT_BITS = (LEVEL_QUANT_BITS == 0) ? 8 : LEVEL_QUANT_BITS;
    fused_projection_bwd_optimizer_3dgs_kernel<
        SplatPrimitive, camera_model, hessian_diagonal_output_mode,
        use_scale_agnostic_mean, color_trust_linear,
        KERNEL_BLOCK_SIZE, KERNEL_QUANT_BITS, LEVEL_VALUE_BITS
    ><<<_CEIL_DIV(N, BLOCK_SIZE_LAUNCH), BLOCK_SIZE_LAUNCH, 0, stream>>>(
        C, N, num_sh_buffer,
        splats_world, viewmats, intrins, dist_coeffs_buffer, image_width, image_height,
        camera_id_bounds, camera_ids, perm, aabb,
        v_splats_world, vr_splats_world, h_splats_world, v_splats_screen, vr_splats_screen, h_splats_screen,
        g1_splats_world, g2_splats_world, sh_packed, sh_quant_bounds,
        sh_value_packed, sh_value_bounds,
        radii, lr_means, lr_quats, lr_scales, lr_opacs, lr_features_dc, lr_features_sh,
        max_gauss_ratio,
        scale_regularization_weight / (float)N,
        mcmc_opacity_reg_weight / (float)N,
        mcmc_scale_reg_weight / (float)N,
        erank_reg_weight / (float)N,
        erank_reg_weight_s3 / (float)N,
        quat_norm_reg_weight / (float)N,
        2.0f * sh_reg_weight / (float)(3*N),
        eps_tr,
        scalar_step, steps);
}
