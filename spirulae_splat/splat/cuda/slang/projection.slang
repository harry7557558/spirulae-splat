#include "projection_utils.slang"

#define _projection_3dgs_ARGS \
    const float3 mean, \
    const float4 quat, \
    const float3 scale, \
    const float in_opacity, \
    const float3x3 R, \
    const float3 t, \
    no_diff const float fx, \
    no_diff const float fy, \
    no_diff const float cx, \
    no_diff const float cy, \
    const uint image_width, \
    const uint image_height, \
    no_diff const float eps2d, \
    no_diff const float near_plane, \
    no_diff const float far_plane, \
    no_diff const float radius_clip, \
    out int2 radii, \
    out float2 mean2d, \
    out float depth, \
    out float3 conic, \
    out float opacity

[ForceInline]
[Differentiable]
void projection_3dgs<let camera_model: CameraModelType, let assume_valid: bool>(
    _projection_3dgs_ARGS
) {

    // transform Gaussian center to camera space
    float3 mean_c;
    posW2C(R, t, mean, mean_c);
    if (!assume_valid && (mean_c.z < near_plane || mean_c.z > far_plane)) {
        radii = {0, 0};
        return;
    }

    // transform Gaussian covariance to camera space
    float3x3 covar;
    {
        // compute from quaternions and scales
        quat_scale_to_covar(quat, scale, covar);
    }
    float3x3 covar_c;
    covarW2C(R, covar, covar_c);

    // perspective projection
    float2x2 covar2d;

    if (camera_model == CameraModelType::PINHOLE)
        persp_proj_3dgs(
            mean_c,
            covar_c,
            fx, fy, cx, cy,
            image_width,
            image_height,
            covar2d,
            mean2d
        );
    else if (camera_model == CameraModelType::FISHEYE)
        fisheye_proj_3dgs(
            mean_c,
            covar_c,
            fx, fy, cx, cy,
            covar2d,
            mean2d
        );
    else if (camera_model == CameraModelType::ORTHO)
        ortho_proj_3dgs(
            mean_c,
            covar_c,
            fx, fy, cx, cy,
            covar2d,
            mean2d
        );

    float compensation;
    float det = add_blur(eps2d, covar2d, compensation);
    if (!assume_valid && det <= 0.f) {
        radii = {0, 0};
        return;
    }

    // compute the inverse of the 2d covariance
    float2x2 covar2d_inv = inverse(covar2d);

    float extend = 3.33f;
    {
        opacity = 1.0 / (1.0 + exp(-in_opacity));
        // opacity *= compensation;  // TODO: mip splatting
        if (!assume_valid && opacity < ALPHA_THRESHOLD) {
            radii = {0, 0};
            return;
        }
        // Compute opacity-aware bounding box.
        // https://arxiv.org/pdf/2402.00525 Section B.2
        extend = min(extend, sqrt(2.0f * log(opacity / ALPHA_THRESHOLD)));
    }

    // compute tight rectangular bounding box (non differentiable)
    // https://arxiv.org/pdf/2402.00525
    float radius_x = ceil(extend * sqrt(covar2d[0].x));
    float radius_y = ceil(extend * sqrt(covar2d[1].y));

    if (!assume_valid && radius_x <= radius_clip && radius_y <= radius_clip) {
        radii = {0, 0};
        return;
    }

    // mask out gaussians outside the image region
    if (!assume_valid && (
        mean2d.x + radius_x <= 0 || mean2d.x - radius_x >= image_width ||
        mean2d.y + radius_y <= 0 || mean2d.y - radius_y >= image_height
    )) {
        radii = {0, 0};
        return;
    }

    radii = {(int)radius_x, (int)radius_y};
    depth = mean_c.z;
    conic = {covar2d_inv[0][0], covar2d_inv[0][1], covar2d_inv[1][1]};
}

[CudaDeviceExport]
void projection_3dgs_persp(_projection_3dgs_ARGS) {
    return projection_3dgs<CameraModelType::PINHOLE, false>(
        mean, quat, scale, in_opacity, R, t, fx, fy, cx, cy,
        image_width, image_height, eps2d, near_plane, far_plane, radius_clip,
        radii, mean2d, depth, conic, opacity
    );
}

[CudaDeviceExport]
void projection_3dgs_fisheye(_projection_3dgs_ARGS) {
    return projection_3dgs<CameraModelType::FISHEYE, false>(
        mean, quat, scale, in_opacity, R, t, fx, fy, cx, cy,
        image_width, image_height, eps2d, near_plane, far_plane, radius_clip,
        radii, mean2d, depth, conic, opacity
    );
}

[CudaDeviceExport]
void projection_3dgs_ortho( _projection_3dgs_ARGS) {
    return projection_3dgs<CameraModelType::ORTHO, false>(
        mean, quat, scale, in_opacity, R, t, fx, fy, cx, cy,
        image_width, image_height, eps2d, near_plane, far_plane, radius_clip,
        radii, mean2d, depth, conic, opacity
    );
}

#define _projection_3dgs_vjp_ARGS \
    const float3 mean, \
    const float4 quat, \
    const float3 scale, \
    const float in_opacity, \
    const float3x3 R, \
    const float3 t, \
    const float fx, \
    const float fy, \
    const float cx, \
    const float cy, \
    const uint image_width, \
    const uint image_height, \
    const float eps2d, \
    const float2 v_mean2d, \
    const float v_depth, \
    const float3 v_conic, \
    const float v_opacity, \
    out float3 v_mean, \
    out float4 v_quat, \
    out float3 v_scale, \
    out float v_in_opacity, \
    out float3x3 v_R, \
    out float3 v_t

[ForceInline]
[Differentiable]
void _projection_3dgs_persp_differentiable(_projection_3dgs_ARGS) {
    return projection_3dgs<CameraModelType::PINHOLE, true>(
        mean, quat, scale, in_opacity, R, t, fx, fy, cx, cy,
        image_width, image_height, eps2d, near_plane, far_plane, radius_clip,
        radii, mean2d, depth, conic, opacity
    );
}

[CudaDeviceExport]
void projection_3dgs_persp_vjp(
    _projection_3dgs_vjp_ARGS
) {
    DifferentialPair<float3> dp_mean = diffPair(mean);
    DifferentialPair<float4> dp_quat = diffPair(quat);
    DifferentialPair<float3> dp_scale = diffPair(scale);
    DifferentialPair<float> dp_in_opacity = diffPair(in_opacity);
    DifferentialPair<float3x3> dp_R = diffPair(R);
    DifferentialPair<float3> dp_t = diffPair(t);
    bwd_diff(_projection_3dgs_persp_differentiable)(
        dp_mean, dp_quat, dp_scale, dp_in_opacity, dp_R, dp_t, fx, fy, cx, cy,
        image_width, image_height, eps2d, 0.0f, 1e30f, -1.0f,
        v_mean2d, v_depth, v_conic, v_opacity
    );
    v_mean = dp_mean.d;
    v_quat = dp_quat.d;
    v_scale = dp_scale.d;
    v_in_opacity = dp_in_opacity.d;
    v_R = dp_R.d;
    v_t = dp_t.d;
}

[ForceInline]
[Differentiable]
void _projection_3dgs_fisheye_differentiable(_projection_3dgs_ARGS) {
    return projection_3dgs<CameraModelType::FISHEYE, true>(
        mean, quat, scale, in_opacity, R, t, fx, fy, cx, cy,
        image_width, image_height, eps2d, near_plane, far_plane, radius_clip,
        radii, mean2d, depth, conic, opacity
    );
}

[CudaDeviceExport]
void projection_3dgs_fisheye_vjp(
    _projection_3dgs_vjp_ARGS
) {
    DifferentialPair<float3> dp_mean = diffPair(mean);
    DifferentialPair<float4> dp_quat = diffPair(quat);
    DifferentialPair<float3> dp_scale = diffPair(scale);
    DifferentialPair<float> dp_in_opacity = diffPair(in_opacity);
    DifferentialPair<float3x3> dp_R = diffPair(R);
    DifferentialPair<float3> dp_t = diffPair(t);
    bwd_diff(_projection_3dgs_fisheye_differentiable)(
        dp_mean, dp_quat, dp_scale, dp_in_opacity, dp_R, dp_t, fx, fy, cx, cy,
        image_width, image_height, eps2d, 0.0f, 1e30f, -1.0f,
        v_mean2d, v_depth, v_conic, v_opacity
    );
    v_mean = dp_mean.d;
    v_quat = dp_quat.d;
    v_scale = dp_scale.d;
    v_in_opacity = dp_in_opacity.d;
    v_R = dp_R.d;
    v_t = dp_t.d;
}

[ForceInline]
[Differentiable]
void _projection_3dgs_ortho_differentiable(_projection_3dgs_ARGS) {
    return projection_3dgs<CameraModelType::ORTHO, true>(
        mean, quat, scale, in_opacity, R, t, fx, fy, cx, cy,
        image_width, image_height, eps2d, near_plane, far_plane, radius_clip,
        radii, mean2d, depth, conic, opacity
    );
}

[CudaDeviceExport]
void projection_3dgs_ortho_vjp(
    _projection_3dgs_vjp_ARGS
) {
    DifferentialPair<float3> dp_mean = diffPair(mean);
    DifferentialPair<float4> dp_quat = diffPair(quat);
    DifferentialPair<float3> dp_scale = diffPair(scale);
    DifferentialPair<float> dp_in_opacity = diffPair(in_opacity);
    DifferentialPair<float3x3> dp_R = diffPair(R);
    DifferentialPair<float3> dp_t = diffPair(t);
    bwd_diff(_projection_3dgs_ortho_differentiable)(
        dp_mean, dp_quat, dp_scale, dp_in_opacity, dp_R, dp_t, fx, fy, cx, cy,
        image_width, image_height, eps2d, 0.0f, 1e30f, -1.0f,
        v_mean2d, v_depth, v_conic, v_opacity
    );
    v_mean = dp_mean.d;
    v_quat = dp_quat.d;
    v_scale = dp_scale.d;
    v_in_opacity = dp_in_opacity.d;
    v_R = dp_R.d;
    v_t = dp_t.d;
}
