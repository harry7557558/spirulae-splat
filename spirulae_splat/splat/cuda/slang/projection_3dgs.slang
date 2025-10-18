#include "projection_utils.slang"


[CudaDeviceExport]
[ForceInline]
[Differentiable]
void persp_proj_3dgs(
    // inputs
    const float3 mean3d,
    const float3x3 cov3d,
    const float fx,
    const float fy,
    const float cx,
    const float cy,
    const uint width,
    const uint height,
    // outputs
    out float2x2 cov2d,
    out float2 mean2d
) {
    float x = mean3d[0], y = mean3d[1], z = mean3d[2];

    float tan_fovx = 0.5f * width / fx;
    float tan_fovy = 0.5f * height / fy;
    float lim_x_pos = (width - cx) / fx + 0.3f * tan_fovx;
    float lim_x_neg = cx / fx + 0.3f * tan_fovx;
    float lim_y_pos = (height - cy) / fy + 0.3f * tan_fovy;
    float lim_y_neg = cy / fy + 0.3f * tan_fovy;

    float rz = 1.f / z;
    float rz2 = rz * rz;
    float tx = z * min(lim_x_pos, max(-lim_x_neg, x * rz));
    float ty = z * min(lim_y_pos, max(-lim_y_neg, y * rz));

    float2x3 J = float2x3(
        fx * rz, 0.f, -fx * tx * rz2,  // 1st row
        0.f, fy * rz, -fy * ty * rz2  // 2nd row
    );
    cov2d = mul(mul(J, cov3d), transpose(J));;
    mean2d = float2(fx * x * rz + cx, fy * y * rz + cy);
}

[CudaDeviceExport]
[ForceInline]
[Differentiable]
void persp_proj_3dgs(
    // inputs
    const float3 mean3d,
    const float3x3 cov3d,
    const float fx,
    const float fy,
    const float cx,
    const float cy,
    // outputs
    out float2x2 cov2d,
    out float2 mean2d
) {
    float x = mean3d[0], y = mean3d[1], z = mean3d[2];

    float rz = 1.f / z;
    float rz2 = rz * rz;
    float tx = x, ty = y;

    float2x3 J = float2x3(
        fx * rz, 0.f, -fx * tx * rz2,  // 1st row
        0.f, fy * rz, -fy * ty * rz2  // 2nd row
    );
    cov2d = mul(mul(J, cov3d), transpose(J));;
    mean2d = float2(fx * x * rz + cx, fy * y * rz + cy);
}

[CudaDeviceExport]
[ForceInline]
[Differentiable]
void fisheye_proj_3dgs(
    // inputs
    const float3 mean3d,
    const float3x3 cov3d,
    const float fx,
    const float fy,
    const float cx,
    const float cy,
    // outputs
    out float2x2 cov2d,
    out float2 mean2d
) {
    float x = mean3d[0], y = mean3d[1], z = mean3d[2];

    float eps = 0.0000001f;
    float xy_len = length(float2(x, y)) + eps;
    float theta = atan2(xy_len, z + eps);
    mean2d = float2(x * fx * theta / xy_len + cx, y * fy * theta / xy_len + cy);

    float x2 = x * x + eps;
    float y2 = y * y;
    float xy = x * y;
    float x2y2 = x2 + y2;
    float x2y2z2_inv = 1.f / (x2y2 + z * z);

    float b = atan2(xy_len, z) / xy_len / x2y2;
    float a = z * x2y2z2_inv / (x2y2);
    float2x3 J = float2x3(
        fx * (x2 * a + y2 * b), fx * xy * (a - b), -fx * x * x2y2z2_inv,
        fy * xy * (a - b), fy * (y2 * a + x2 * b), -fy * y * x2y2z2_inv
    );
    cov2d = mul(mul(J, cov3d), transpose(J));;
}

[CudaDeviceExport]
[ForceInline]
[Differentiable]
void ortho_proj_3dgs(
    // inputs
    const float3 mean3d,
    const float3x3 cov3d,
    const float fx,
    const float fy,
    const float cx,
    const float cy,
    // outputs
    out float2x2 cov2d,
    out float2 mean2d
) {
    float x = mean3d[0], y = mean3d[1], z = mean3d[2];

    float2x3 J = float2x3(
        fx, 0.f, 0.f,
        0.f, fy, 0.f
    );
    cov2d = mul(mul(J, cov3d), transpose(J));;
    mean2d = float2(fx * x + cx, fy * y + cy);
}



[CudaDeviceExport]
[ForceInline]
[Differentiable]
float add_blur(const float eps2d, inout float2x2 covar, out float compensation) {
    float det_orig = covar[0][0] * covar[1][1] - covar[0][1] * covar[1][0];
    covar[0][0] += eps2d;
    covar[1][1] += eps2d;
    float det_blur = covar[0][0] * covar[1][1] - covar[0][1] * covar[1][0];
    compensation = sqrt(max(0.f, det_orig / det_blur));
    return det_blur;
}


#define _projection_3dgs_ARGS \
    const bool antialiased, \
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
    out float depth, \
    out float2 mean2d, \
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
        quat_scale_to_covar(quat, exp(scale), covar);
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
        if (antialiased)
            opacity *= compensation;
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

#define _DEF_projection_3dgs(funname, ...) \
    [CudaDeviceExport] \
    [ForceInline] \
    [Differentiable] \
    void funname(_projection_3dgs_ARGS) { \
        return projection_3dgs<__VA_ARGS__>( \
            antialiased, \
            mean, quat, scale, in_opacity, R, t, fx, fy, cx, cy, \
            image_width, image_height, eps2d, near_plane, far_plane, radius_clip, \
            radii, depth, mean2d, conic, opacity \
        ); \
    }

_DEF_projection_3dgs(projection_3dgs_persp, CameraModelType::PINHOLE, false)
_DEF_projection_3dgs(projection_3dgs_fisheye, CameraModelType::FISHEYE, false)
_DEF_projection_3dgs(projection_3dgs_ortho, CameraModelType::ORTHO, false)


#define _projection_3dgs_vjp_ARGS \
    const bool antialiased, \
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
    const float v_depth, \
    const float2 v_mean2d, \
    const float3 v_conic, \
    const float v_opacity, \
    out float3 v_mean, \
    out float4 v_quat, \
    out float3 v_scale, \
    out float v_in_opacity, \
    out float3x3 v_R, \
    out float3 v_t

_DEF_projection_3dgs(_projection_3dgs_persp_differentiable, CameraModelType::PINHOLE, true)
_DEF_projection_3dgs(_projection_3dgs_fisheye_differentiable, CameraModelType::FISHEYE, true)
_DEF_projection_3dgs(_projection_3dgs_ortho_differentiable, CameraModelType::ORTHO, true)


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
        antialiased,
        dp_mean, dp_quat, dp_scale, dp_in_opacity, dp_R, dp_t, fx, fy, cx, cy,
        image_width, image_height, eps2d, 0.0f, 1e30f, -1.0f,
        v_depth, v_mean2d, v_conic, v_opacity
    );
    v_mean = dp_mean.d;
    v_quat = dp_quat.d;
    v_scale = dp_scale.d;
    v_in_opacity = dp_in_opacity.d;
    v_R = dp_R.d;
    v_t = dp_t.d;
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
        antialiased,
        dp_mean, dp_quat, dp_scale, dp_in_opacity, dp_R, dp_t, fx, fy, cx, cy,
        image_width, image_height, eps2d, 0.0f, 1e30f, -1.0f,
        v_depth, v_mean2d, v_conic, v_opacity
    );
    v_mean = dp_mean.d;
    v_quat = dp_quat.d;
    v_scale = dp_scale.d;
    v_in_opacity = dp_in_opacity.d;
    v_R = dp_R.d;
    v_t = dp_t.d;
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
        antialiased,
        dp_mean, dp_quat, dp_scale, dp_in_opacity, dp_R, dp_t, fx, fy, cx, cy,
        image_width, image_height, eps2d, 0.0f, 1e30f, -1.0f,
        v_depth, v_mean2d, v_conic, v_opacity
    );
    v_mean = dp_mean.d;
    v_quat = dp_quat.d;
    v_scale = dp_scale.d;
    v_in_opacity = dp_in_opacity.d;
    v_R = dp_R.d;
    v_t = dp_t.d;
}
