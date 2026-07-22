// DensifyScoring.cu -- covariance-based scale initialization and the densification parameter update.
//
// Part of the Densify family -- see DensifyCommon.cuh.

#include "DensifyCommon.cuh"

// ================
// Covariance-Based Scale Initialization
// ================

__global__ void cov_scale_init_kernel(
    int64_t num_points,
    int32_t num_cameras,
    const float3* __restrict__ points,  // [N, 3]
    const bool* __restrict__ is_fisheye,  // [C]; TODO: equisolid
    const int2* __restrict__ sizes,  // [C, 2]
    const float4 *__restrict__ intrins,  // [C, 4], fx, fy, cx, cy
    const float4 *__restrict__ viewmats,  // [C, 4, 4]
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer, // [C]
    float* __restrict__ log_scales  // [N]
) {
    int64_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_points)
        return;

    float3 p_world = points[idx];
    float max_log_scale = __logf(1e-30f);

    for (int32_t i = 0; i < num_cameras; ++i) {
        float4 intrin = intrins[i];
        int width = sizes[i].x, height = sizes[i].y;
        CameraDistortionCoeffs dist_coeffs = dist_coeffs_buffer.load(i);

        float4 p_wh = {p_world.x, p_world.y, p_world.z, 1.0f};
        float3 p_view = {
            dot(viewmats[4*i], p_wh),
            dot(viewmats[4*i+1], p_wh),
            dot(viewmats[4*i+2], p_wh),
        };

        bool valid = false;
        constexpr float eps = 1e-6f;
        float3x3 cov3d = {eps, 0, 0, 0, eps, 0, 0, 0, eps};
        float2x2 cov2d;
        float2 mean2d;
        if (is_fisheye[i]) {
            valid = SlangProjectionUtils::fisheye_proj_3dgs_nav(
                p_view, cov3d, intrin, dist_coeffs, &cov2d, &mean2d
            );
        }
        else {
            valid = SlangProjectionUtils::persp_proj_3dgs_nav(
                p_view, cov3d, intrin, dist_coeffs, width, height, &cov2d, &mean2d
            );
        }

        #pragma unroll
        for (int i = 0; i < 2; ++i) {
            cov2d[i].x = __fmul_rn(cov2d[i].x, 1.0f/eps);
            cov2d[i].y = __fmul_rn(cov2d[i].y, 1.0f/eps);
        }

        float det = cov2d[0].x * cov2d[1].y - cov2d[0].y * cov2d[1].x;
        if (valid && det > 0.0f) {
            float sc = 0.5f * __logf((float)(width * height) / det);
            max_log_scale = fmaxf(max_log_scale, sc);
        }
    }

    log_scales[idx] = max_log_scale;
}

/*[AutoHeaderGeneratorExport]*/
void cov_scale_init_tensor(
    DeviceVector<float3> points,  // [N, 3]
    DeviceVector<bool> is_fisheye,  // [C], bool
    DeviceVector<int2> sizes,  // [C, 2], int32
    DeviceVector<float4> intrins,  // [C, 4]
    DeviceVector<float4> viewmats,  // [C, 4, 4] as 4*C float4 elements
    TorchTensorView dist_coeffs, // [C]
    DeviceVector<float> log_scales  // [N, 1] output
) {
    int64_t N = points.size();
    int64_t C = intrins.size();

    cov_scale_init_kernel<<<_LAUNCH_ARGS_1D(N, 256)>>>(
        N, C,
        points.data_ptr(),
        is_fisheye.data_ptr(),
        sizes.data_ptr(),
        intrins.data_ptr(),
        viewmats.data_ptr(),
        dist_coeffs,
        log_scales.data_ptr()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}



// ================
// Densification Parameter Update
// ================

__global__ void densify_clip_scale_kernel(
    long num_splats,
    float max_scale2d,
    float clip_hardness,
    float max_scale3d,
    const float* __restrict__ radii,
    float3* __restrict__ log_scales,
    float* __restrict__ logit_opacs
) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_splats)
        return;

    if (isfinite(max_scale2d)) {
        float oversize_factor = radii[idx] / max_scale2d;
        if (isfinite(clip_hardness))
            oversize_factor = fminf(oversize_factor, clip_hardness);
        if (oversize_factor > 1.0f) {
            oversize_factor = __logf(oversize_factor);
            log_scales[idx] = log_scales[idx] - make_float3(oversize_factor);
            // this encourages being relocated to, may cause unintended effects in non-MCMC
            if (logit_opacs != nullptr)
                logit_opacs[idx] = fminf(logit_opacs[idx] + 3.0f * oversize_factor, fmaxf(logit_opacs[idx], 5.0f));
        }
    }

    if (isfinite(max_scale3d)) {
        max_scale3d = __logf(max_scale3d);
        log_scales[idx].x = fminf(log_scales[idx].x, max_scale3d);
        log_scales[idx].y = fminf(log_scales[idx].y, max_scale3d);
        log_scales[idx].z = fminf(log_scales[idx].z, max_scale3d);
        // don't touch opacity
    }
}

/*[AutoHeaderGeneratorExport]*/
void densify_clip_scale_tensor(
    int64_t num_splats,
    DeviceVector<float> radii,
    DeviceVector<float3> log_scales,
    float* logit_opacs_ptr,
    float max_scale2d,
    float clip_hardness,
    float max_scale3d
) {
    densify_clip_scale_kernel<<<_LAUNCH_ARGS_1D(num_splats, 256)>>>(
        num_splats, max_scale2d, clip_hardness, max_scale3d,
        radii.data_ptr(),
        log_scales.data_ptr(),
        logit_opacs_ptr
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


// https://www.shadertoy.com/view/4djSRW
__device__ __forceinline__ float hash14(float4 p4) {
    p4 = p4 * float4{.1031, .1030, .0973, .1099};
    p4.x = fmodf(p4.x, 1.0f); p4.y = fmodf(p4.y, 1.0f); p4.z = fmodf(p4.z, 1.0f); p4.w = fmodf(p4.w, 1.0f);
    p4 = p4 + dot(p4, float4{p4.w, p4.z, p4.x, p4.y} + 33.33f);
    return fmodf((p4.x + p4.y) * (p4.z + p4.w), 1.0f);
}

__global__ void densify_update_weight_kernel(
    long num_splats,
    int score_mode,
    const float* __restrict__ radii,  // [N]
    const float3* __restrict__ scales,  // [N, 3], optional
    const float* __restrict__ opacs,  // [N], optional
    const float* __restrict__ accum_weight_scalar,  // [1]
    const float* __restrict__ accum_weight,  // [N]
    const float* __restrict__ accum_weight2,  // [N], optional blend partner
    float blend_w,  // weight of accum_weight2 in the geometric blend
    float2* __restrict__ accum_buffer  // [N, 2]
) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_splats)
        return;

    if (!(radii[idx] > 0))
        return;

    float weight = fabsf(accum_weight[idx]);
    // Geometric blend: weight^(1-w) * |accum_weight2|^w. Invariant (up to a
    // global factor) to each score's overall scale, so relative ranking --
    // all downstream consumers care about -- needs no cross-normalization.
    // The w == 0 / w == 1 cases are short-circuited by the caller passing a
    // single buffer, so no powf is spent there.
    if (accum_weight2 != nullptr) {
        weight = powf(weight, 1.0f - blend_w) * powf(fabsf(accum_weight2[idx]), blend_w);
        if (!isfinite(weight))
            weight = 0.0f;
    }
    if (opacs)
        weight *= sigmoid(opacs[idx]);
    if (accum_weight_scalar != nullptr)
        weight *= accum_weight_scalar[0];
    if (weight == 0.0f)
        return;

    float2 accum = accum_buffer[idx];
    if (score_mode == (int)DensifyScoreMode::Max) {
        accum.x = fmaxf(accum.x, weight);
        accum.y = fmaxf(accum.y, 1.0f);
    } else if (score_mode == (int)DensifyScoreMode::Mean) {
        accum.x = (accum.x * accum.y + weight) / (accum.y + 1.0f);
        accum.y += 1.0f;
    } else if (score_mode == (int)DensifyScoreMode::Median) {
        float rand = hash14(1e2f * float4{weight, accum.x, accum.y, accum.y + 1.5f});
        if (accum.y == 0.0f) {
            accum.x = weight * exp2f(rand - 0.5f);
            accum.y = 1.0f;
        } else if (weight != 0.0f) {
            float s = weight > accum.x ? 1.0f : -1.0f;
            // accum.x *= exp2f(s / sqrtf(accum.y + 1.0f));
            accum.x *= exp2f(s * rand);
            accum.y += 1.0f;
        };
    } else if (score_mode == (int)DensifyScoreMode::Geom) {
        accum.x = __expf((__logf(fmaxf(accum.x, 1e-30f)) * accum.y +  __logf(fmaxf(weight, 1e-30f))) / (accum.y + 1.0f));
        accum.y += 1.0f;
    }
    accum_buffer[idx] = accum;
}


/*[AutoHeaderGeneratorExport]*/
void densify_update_weight(
    int64_t num_splats,
    DeviceVector<float> radii,
    float3* scales_ptr,
    float* opacs_ptr,
    DeviceVector<float> accum_weight,
    DeviceVector<float> accum_weight2,
    float blend_w,
    DeviceVector<float2> accum_buffer,
    int score_mode
) {
    densify_update_weight_kernel<<<_LAUNCH_ARGS_1D(num_splats, 256)>>>(
        num_splats, score_mode,
        radii.data_ptr(),
        scales_ptr,
        opacs_ptr,
        nullptr,
        accum_weight.data_ptr(),
        accum_weight2.data_ptr(),
        blend_w,
        accum_buffer.data_ptr()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


