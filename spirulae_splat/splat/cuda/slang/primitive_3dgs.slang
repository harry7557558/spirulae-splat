#pragma once
#include "projection_utils.slang"

struct SigmaPoints {
    float3 p[7];
    float w_mean[7];
    float w_cov[7];
};

[ForceInline]
SigmaPoints world_gaussian_sigma_points(
    float3 world_mean,
    float3 world_scale,
    float4 world_quat,
    float3x3 R, float3 t
) {
    static const float alpha = 1.0;
    static const float beta = 2.0;
    static const float kappa = 0.0;
    static const uint dim = 3;
    static const float lambda = alpha * alpha * (dim + kappa) - dim;
    static const float r = sqrt(dim + lambda);

    float3x3 Rq = transpose(normalized_quat_to_rotmat(world_quat));

    SigmaPoints ret;
    ret.p[0] = world_mean;
    ret.w_mean[0] = lambda / (dim + lambda);
    ret.w_cov[0] = lambda / (dim + lambda) + (1.0 - alpha * alpha + beta);

    [ForceUnroll]
    for (uint i = 0; i < dim; i++) {
        float3 delta = sqrt(dim + lambda) * world_scale[i] * Rq[i];
        ret.p[i + 1] = world_mean + delta;
        ret.p[i + 1 + dim] = world_mean - delta;
    }

    [ForceUnroll]
    for (uint i = 1; i <= 2 * dim; i++) {
        ret.w_mean[i] = 0.5 / (dim + lambda);
        ret.w_cov[i] = 0.5 / (dim + lambda);
    }

    [ForceUnroll]
    for (uint i = 0; i < 7; i++) {
        float3 temp = ret.p[i];
        posW2C(R, t, temp, ret.p[i]);
    }

    return ret;
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
[Differentiable]
bool persp_proj_3dgs<let assume_valid : bool>(
    // inputs
    const float3 mean3d,
    const float3x3 cov3d,
    const float4 intrins,
    const float[10] dist_coeffs,
    // outputs
    out float2x2 cov2d,
    out float2 mean2d
) {
    bool valid = persp_proj<assume_valid>(mean3d, intrins, dist_coeffs, mean2d);
    if (!assume_valid && !valid)
        return false;

    float2x3 J = persp_proj_jac(mean3d, intrins, dist_coeffs);
    cov2d = mul(mul(J, cov3d), transpose(J));
    return true;
}

[CudaDeviceExport]
bool persp_proj_3dgs_ut<let assume_valid : bool>(
    // inputs
    const SigmaPoints sigmas,
    const float4 intrins,
    const float[10] dist_coeffs,
    const uint width,
    const uint height,
    // outputs
    out float2x2 cov2d,
    out float2 mean2d
) {
    mean2d = float2(0);
    cov2d = float2x2(0);

    float fx = intrins.x, fy = intrins.y;
    float tan_fovx = 0.5f * width / fx;
    float tan_fovy = 0.5f * height / fy;
    float lim_x_pos = width + 0.3f * tan_fovx * fx;
    float lim_x_neg = 0.3f * tan_fovx * fx;
    float lim_y_pos = height + 0.3f * tan_fovy * fy;
    float lim_y_neg = 0.3f * tan_fovy * fy;

    float2 proj_points[7];
    bool all_valid = true;
    [ForceUnroll] for (int i = 0; i < 7; i++)
        all_valid &= persp_proj<assume_valid>(sigmas.p[i], intrins, dist_coeffs, proj_points[i]);
    if (!assume_valid && !all_valid)
        return false;

    [ForceUnroll] for (int i = 0; i < 7; i++)
        mean2d += sigmas.w_mean[i] * proj_points[i];

    float2 mean2d_clipped = {
        clamp(mean2d.x, -lim_x_neg, lim_x_pos),
        clamp(mean2d.y, -lim_y_neg, lim_y_pos)
    };
    [ForceUnroll] for (int i = 0; i < 7; i++) {
        float2 d = float2(
            clamp(proj_points[i].x, -lim_x_neg, lim_x_pos),
            clamp(proj_points[i].y, -lim_y_neg, lim_y_pos)
        ) - mean2d_clipped;
        cov2d += sigmas.w_cov[i] * float2x2(d.x * d.x, d.x * d.y, d.x * d.y, d.y * d.y);
    }
    return true;
}

[CudaDeviceExport]
[Differentiable]
bool fisheye_proj_3dgs<let assume_valid : bool>(
    // inputs
    const float3 mean3d,
    const float3x3 cov3d,
    const float4 intrins,
    const float[10] dist_coeffs,
    // outputs
    out float2x2 cov2d,
    out float2 mean2d
) {
    bool valid = fisheye_proj<assume_valid>(mean3d, intrins, dist_coeffs, mean2d);
    if (!assume_valid && !valid)
        return false;

    float2x3 J = fisheye_proj_jac(mean3d, intrins, dist_coeffs);
    cov2d = mul(mul(J, cov3d), transpose(J));
    return true;
}

[CudaDeviceExport]
bool fisheye_proj_3dgs_ut<let assume_valid : bool>(
    // inputs
    const SigmaPoints sigmas,
    const float4 intrins,
    const float[10] dist_coeffs,
    // outputs
    out float2x2 cov2d,
    out float2 mean2d
) {
    mean2d = float2(0);
    cov2d = float2x2(0);

    float2 proj_points[7];
    bool all_valid = true;
    [ForceUnroll] for (int i = 0; i < 7; i++)
        all_valid &= fisheye_proj<assume_valid>(sigmas.p[i], intrins, dist_coeffs, proj_points[i]);
    if (!assume_valid && !all_valid)
        return false;

    [ForceUnroll] for (int i = 0; i < 7; i++)
        mean2d += sigmas.w_mean[i] * proj_points[i];

    [ForceUnroll] for (int i = 0; i < 7; i++) {
        float2 d = proj_points[i] - mean2d;
        cov2d += sigmas.w_cov[i] * float2x2(d.x * d.x, d.x * d.y, d.x * d.y, d.y * d.y);
    }
    return true;
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

    float4 intrins = { fx, fy, cx, cy };
    float2x2 covar2d;

    bool is_valid = true;
    if (eval3d) {
        SigmaPoints sigmas = world_gaussian_sigma_points(mean, exp(scale), quat, R, t);

        if (camera_model == CameraModelType::PINHOLE)
            is_valid &= persp_proj_3dgs_ut<assume_valid>(
                sigmas,
                intrins, dist_coeffs, image_width, image_height,
                covar2d, mean2d
            );
        else if (camera_model == CameraModelType::FISHEYE)
            is_valid &= fisheye_proj_3dgs_ut<assume_valid>(
                sigmas,
                intrins, dist_coeffs,
                covar2d, mean2d
            );
    }
    else {
        // transform Gaussian covariance to camera space
        float3x3 covar;
        {
            // compute from quaternions and scales
            quat_scale_to_covar(quat, exp(scale), covar);
        }
        float3x3 covar_c;
        covarW2C(R, covar, covar_c);

        if (camera_model == CameraModelType::PINHOLE)
            // does not support dist_coeffs - asserted on Python side
            persp_proj_3dgs(
                mean_c, covar_c,
                fx, fy, cx, cy, image_width, image_height,
                covar2d, mean2d
            );
        else if (camera_model == CameraModelType::FISHEYE)
            is_valid &= fisheye_proj_3dgs<assume_valid>(
                mean_c, covar_c,
                intrins, dist_coeffs,
                covar2d, mean2d
            );
        else if (camera_model == CameraModelType::ORTHO)
            ortho_proj_3dgs(
                mean_c, covar_c,
                fx, fy, cx, cy,
                covar2d, mean2d
            );
    }
    if (!assume_valid && !is_valid) {
        aabb_xyxy = { 0, 0, 0, 0 };
        return;
    }

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
_DEF_projection_3dgs(projection_3dgut_persp, CameraModelType::PINHOLE, false, true)
_DEF_projection_3dgs(projection_3dgut_fisheye, CameraModelType::FISHEYE, false, true)


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
_DEF_projection_3dgs(_projection_3dgut_persp_differentiable, CameraModelType::PINHOLE, true, true)
_DEF_projection_3dgs(_projection_3dgut_fisheye_differentiable, CameraModelType::FISHEYE, true, true)


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
_DEF_projection_3dgs_vjp(projection_3dgut_persp_vjp, _projection_3dgut_persp_differentiable)
_DEF_projection_3dgs_vjp(projection_3dgut_fisheye_vjp, _projection_3dgut_fisheye_differentiable)


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
    // depth = 0.5 * log(dot(mean - ray_o, mean - ray_o) + 1e-6);
    float3 gro = mul(iscl_rot, (ray_o - mean));
    float3 grd = mul(iscl_rot, ray_d);
    depth = -dot(gro, grd) / dot(grd, grd);
    depth = log(max(depth, 1e-6));
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
