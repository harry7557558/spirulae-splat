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
// [ForceInline]
[Differentiable]
void persp_proj_3dgs(
    // inputs
    const float3 mean3d,
    const float3x3 cov3d,
    const float fx,
    const float fy,
    const float cx,
    const float cy,
    const float[10] dist_coeffs,
    // outputs
    out float2x2 cov2d,
    out float2 mean2d
) {
    float4 intrins = float4(fx, fy, cx, cy);

    mean2d = persp_proj(mean3d, intrins, dist_coeffs);

    float2x3 J = persp_proj_jac(mean3d, intrins, dist_coeffs);
    cov2d = mul(mul(J, cov3d), transpose(J));
}

[CudaDeviceExport]
// [ForceInline]
[Differentiable]
void fisheye_proj_3dgs(
    // inputs
    const float3 mean3d,
    const float3x3 cov3d,
    const float fx,
    const float fy,
    const float cx,
    const float cy,
    const float[10] dist_coeffs,
    // outputs
    out float2x2 cov2d,
    out float2 mean2d
) {
    float4 intrins = float4(fx, fy, cx, cy);

    mean2d = fisheye_proj(mean3d, intrins, dist_coeffs);

    float2x3 J = fisheye_proj_jac(mean3d, intrins, dist_coeffs);
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
    const float3 sh_coeffs[16], \
    const float3x3 R, \
    const float3 t, \
    no_diff const float fx, \
    no_diff const float fy, \
    no_diff const float cx, \
    no_diff const float cy, \
    no_diff const float[10] dist_coeffs, \
    const uint image_width, \
    const uint image_height, \
    no_diff const float near_plane, \
    no_diff const float far_plane, \
    out int4 aabb_xyxy, \
    out float2 mean2d, \
    out float depth, \
    out float3 conic, \
    out float opacity, \
    out float3 rgb

[ForceInline]
[Differentiable]
void projection_3dgs<let camera_model: CameraModelType, let assume_valid: bool, let eval3d: bool>(
    _projection_3dgs_ARGS
) {

    // transform Gaussian center to camera space
    float3 mean_c;
    posW2C(R, t, mean, mean_c);
    float cull_z = (camera_model == CameraModelType::FISHEYE) ?
        length(mean_c) : mean_c.z;
    if (!assume_valid && (cull_z < near_plane || cull_z > far_plane)) {
        aabb_xyxy = {0, 0, 0, 0};
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
            dist_coeffs,
            covar2d,
            mean2d
        );
    else if (camera_model == CameraModelType::FISHEYE)
        fisheye_proj_3dgs(
            mean_c,
            covar_c,
            fx, fy, cx, cy,
            dist_coeffs,
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
    static const float eps2d = 0.3;
    float det = add_blur(eps2d, covar2d, compensation);
    if (!assume_valid && det <= 0.f) {
        aabb_xyxy = {0, 0, 0, 0};
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
            aabb_xyxy = {0, 0, 0, 0};
            return;
        }
        // Compute opacity-aware bounding box.
        // https://arxiv.org/pdf/2402.00525 Section B.2
        extend = min(extend, sqrt(2.0f * log(opacity / ALPHA_THRESHOLD)));
    }

    // compute tight rectangular bounding box (non differentiable)
    // https://arxiv.org/pdf/2402.00525
    float radius_x = extend * sqrt(covar2d[0].x);
    float radius_y = extend * sqrt(covar2d[1].y);
    float xmin = floor(mean2d.x - radius_x);
    float xmax = ceil(mean2d.x + radius_x);
    float ymin = floor(mean2d.y - radius_y);
    float ymax = ceil(mean2d.y + radius_y);

    // mask out gaussians outside the image region
    if (!assume_valid && (
        xmax <= 0 || xmin >= image_width ||
        ymax <= 0 || ymin >= image_height
    )) {
        aabb_xyxy = {0, 0, 0, 0};
        return;
    }

    aabb_xyxy = {(int)xmin, (int)ymin, (int)xmax, (int)ymax};
    // depth = mean_c.z;
    // depth = length(mean_c);
    depth = 0.5*log(dot(mean_c, mean_c) + 1e-6);
    conic = eval3d ? exp(-scale) :
        float3(covar2d_inv[0][0], covar2d_inv[0][1], covar2d_inv[1][1]);

    float3 campos = -mul(transpose(R), t);
    sh_coeffs_to_color(mean - campos, sh_coeffs, rgb);
    rgb = max(rgb+0.5f, 0.0f);
}

#define _DEF_projection_3dgs(funname, ...) \
    [CudaDeviceExport] \
    [ForceInline] \
    [Differentiable] \
    void funname(_projection_3dgs_ARGS) { \
        return projection_3dgs<__VA_ARGS__>( \
            antialiased, \
            mean, quat, scale, in_opacity, sh_coeffs, R, t, \
            fx, fy, cx, cy, dist_coeffs, \
            image_width, image_height, near_plane, far_plane, \
            aabb_xyxy, mean2d, depth, conic, opacity, rgb \
        ); \
    }

_DEF_projection_3dgs(projection_3dgs_persp, CameraModelType::PINHOLE, false, false)
_DEF_projection_3dgs(projection_3dgs_fisheye, CameraModelType::FISHEYE, false, false)
_DEF_projection_3dgs(projection_3dgs_ortho, CameraModelType::ORTHO, false, false)
_DEF_projection_3dgs(projection_3dgs_eval3d_persp, CameraModelType::PINHOLE, false, true)
_DEF_projection_3dgs(projection_3dgs_eval3d_fisheye, CameraModelType::FISHEYE, false, true)


#define _projection_3dgs_vjp_ARGS \
    const bool antialiased, \
    const float3 mean, \
    const float4 quat, \
    const float3 scale, \
    const float in_opacity, \
    const float3 sh_coeffs[16], \
    const float3x3 R, \
    const float3 t, \
    const float fx, \
    const float fy, \
    const float cx, \
    const float cy, \
    const float[10] dist_coeffs, \
    const uint image_width, \
    const uint image_height, \
    const float2 v_mean2d, \
    const float v_depth, \
    const float3 v_conic, \
    const float v_opacity, \
    const float3 v_rgb, \
    out float3 v_mean, \
    out float4 v_quat, \
    out float3 v_scale, \
    out float v_in_opacity, \
    out float3 v_sh_coeffs[16], \
    out float3x3 v_R, \
    out float3 v_t

_DEF_projection_3dgs(_projection_3dgs_persp_differentiable, CameraModelType::PINHOLE, true, false)
_DEF_projection_3dgs(_projection_3dgs_fisheye_differentiable, CameraModelType::FISHEYE, true, false)
_DEF_projection_3dgs(_projection_3dgs_ortho_differentiable, CameraModelType::ORTHO, true, false)
_DEF_projection_3dgs(_projection_3dgs_eval3d_persp_differentiable, CameraModelType::PINHOLE, true, true)
_DEF_projection_3dgs(_projection_3dgs_eval3d_fisheye_differentiable, CameraModelType::FISHEYE, true, true)


#define _DEF_projection_3dgs_vjp(funname, fwd_funname, ...) \
    [CudaDeviceExport] \
    void funname( \
        _projection_3dgs_vjp_ARGS \
    ) { \
        DifferentialPair<float3> dp_mean = diffPair(mean); \
        DifferentialPair<float4> dp_quat = diffPair(quat); \
        DifferentialPair<float3> dp_scale = diffPair(scale); \
        DifferentialPair<float> dp_in_opacity = diffPair(in_opacity); \
        DifferentialPair<float3[16]> dp_sh_coeffs = diffPair(sh_coeffs); \
        DifferentialPair<float3x3> dp_R = diffPair(R); \
        DifferentialPair<float3> dp_t = diffPair(t); \
        bwd_diff(fwd_funname)( \
            antialiased, \
            dp_mean, dp_quat, dp_scale, dp_in_opacity, dp_sh_coeffs, dp_R, dp_t, \
            fx, fy, cx, cy, dist_coeffs, \
            image_width, image_height, 0.0f, 1e30f, \
            v_mean2d, v_depth, v_conic, v_opacity, v_rgb \
        ); \
        v_mean = dp_mean.d; \
        v_quat = dp_quat.d; \
        v_scale = dp_scale.d; \
        v_in_opacity = dp_in_opacity.d; \
        v_sh_coeffs = dp_sh_coeffs.d; \
        v_R = dp_R.d; \
        v_t = dp_t.d; \
    }

_DEF_projection_3dgs_vjp(projection_3dgs_persp_vjp, _projection_3dgs_persp_differentiable)
_DEF_projection_3dgs_vjp(projection_3dgs_fisheye_vjp, _projection_3dgs_fisheye_differentiable)
_DEF_projection_3dgs_vjp(projection_3dgs_ortho_vjp, _projection_3dgs_ortho_differentiable)
_DEF_projection_3dgs_vjp(projection_3dgs_eval3d_persp_vjp, _projection_3dgs_eval3d_persp_differentiable)
_DEF_projection_3dgs_vjp(projection_3dgs_eval3d_fisheye_vjp, _projection_3dgs_eval3d_fisheye_differentiable)


[CudaDeviceExport]
[Differentiable]
float3x3 compute_3dgut_iscl_rot(float4 quat, float3 scale) {
    float3x3 R = normalized_quat_to_rotmat(quat);
    float3x3 S = float3x3(
        scale.x, 0.f, 0.f,
        0.f, scale.y, 0.f,
        0.f, 0.f, scale.z
    );
    return mul(S, transpose(R));
}

[CudaDeviceExport]
void compute_3dgut_iscl_rot_vjp(
    float4 quat, float3 scale,
    float3x3 v_iscl_rot,
    out float4 v_quat, out float3 v_scale
) {
    DifferentialPair<float4> dp_quat = diffPair(quat);
    DifferentialPair<float3> dp_scale = diffPair(scale);
    bwd_diff(compute_3dgut_iscl_rot)(
        dp_quat, dp_scale, v_iscl_rot
    );
    v_quat = dp_quat.d;
    v_scale = dp_scale.d;
}

[CudaDeviceExport]
[Differentiable]
float evaluate_alpha_3dgs(
    float3 mean, float3x3 iscl_rot, float opacity,
    float3 ray_o, float3 ray_d
) {
    float3 gro = mul(iscl_rot, (ray_o - mean));
    float3 grd = mul(iscl_rot, ray_d);
    float3 gcrod = cross(grd, gro);
    float val = exp(-0.5f * dot(gcrod, gcrod) / dot(grd, grd));
    return opacity * val;
}

[CudaDeviceExport]
void evaluate_alpha_3dgs_vjp(
    float3 mean, float3x3 iscl_rot, float opacity,
    float3 ray_o, float3 ray_d, float v_alpha,
    out float3 v_mean, out float3x3 v_iscl_rot, out float v_opacity,
    out float3 v_ray_o, out float3 v_ray_d
) {
    DifferentialPair<float3> dp_mean = diffPair(mean);
    DifferentialPair<float3x3> dp_iscl_rot = diffPair(iscl_rot);
    DifferentialPair<float> dp_opacity = diffPair(opacity);
    DifferentialPair<float3> dp_ray_o = diffPair(ray_o);
    DifferentialPair<float3> dp_ray_d = diffPair(ray_d);
    bwd_diff(evaluate_alpha_3dgs)(
        dp_mean, dp_iscl_rot, dp_opacity,
        dp_ray_o, dp_ray_d, v_alpha
    );
    v_mean = dp_mean.d;
    v_iscl_rot = dp_iscl_rot.d;
    v_opacity = dp_opacity.d;
    v_ray_o = dp_ray_o.d;
    v_ray_d = dp_ray_d.d;
}

[CudaDeviceExport]
[Differentiable]
void evaluate_color_3dgs(
    float3 mean, float3x3 iscl_rot, float opacity, float3 rgb,
    float3 ray_o, float3 ray_d,
    out float3 out_rgb, out float depth
) {
    out_rgb = rgb;
    // depth = length(mean - ray_o);  // TODO: more accurate depth for large splats
    depth = 0.5*log(dot(mean-ray_o, mean-ray_o) + 1e-6);
}

[CudaDeviceExport]
void evaluate_color_3dgs_vjp(
    float3 mean, float3x3 iscl_rot, float opacity, float3 rgb,
    float3 ray_o, float3 ray_d, float3 v_out_rgb, float v_depth,
    out float3 v_mean, out float3x3 v_iscl_rot, out float v_opacity, out float3 v_rgb,
    out float3 v_ray_o, out float3 v_ray_d
) {
    DifferentialPair<float3> dp_mean = diffPair(mean);
    DifferentialPair<float3x3> dp_iscl_rot = diffPair(iscl_rot);
    DifferentialPair<float> dp_opacity = diffPair(opacity);
    DifferentialPair<float3> dp_rgb = diffPair(rgb);
    DifferentialPair<float3> dp_ray_o = diffPair(ray_o);
    DifferentialPair<float3> dp_ray_d = diffPair(ray_d);
    bwd_diff(evaluate_color_3dgs)(
        dp_mean, dp_iscl_rot, dp_opacity, dp_rgb,
        dp_ray_o, dp_ray_d, v_out_rgb, v_depth
    );
    v_mean = dp_mean.d;
    v_iscl_rot = dp_iscl_rot.d;
    v_opacity = dp_opacity.d;
    v_rgb = dp_rgb.d;
    v_ray_o = dp_ray_o.d;
    v_ray_d = dp_ray_d.d;
}
