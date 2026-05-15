#pragma once

// suppress CudaDeviceExport to reduce code size
#ifdef CudaDeviceExport
#include "projection_utils.slang"
#else
#define CudaDeviceExport Differentiable
#include "projection_utils.slang"
#undef CudaDeviceExport
#endif

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

[ForceInline]
[Differentiable]
bool persp_proj_3dgs<let assume_valid : bool>(
    // inputs
    const float3 mean3d,
    const float3x3 cov3d,
    const float4 intrins,
    const float[10] dist_coeffs,
    const uint width,
    const uint height,
    // outputs
    out float2x2 cov2d,
    out float2 mean2d
) {
    float fx = intrins.x, fy = intrins.y, cx = intrins.z, cy = intrins.w;

    bool valid = persp_proj<assume_valid>(mean3d, intrins, dist_coeffs, mean2d);
    if (!assume_valid && !valid)
        return false;

    float2x3 J = persp_proj_jac(mean3d, intrins, dist_coeffs);

    float lim_x_pos = width + 0.3f * (0.5f * width);
    float lim_x_neg = 0.3f * (0.5f * width);
    float lim_y_pos = height + 0.3f * (0.5f * height);
    float lim_y_neg = 0.3f * (0.5f * height);

    float rz = 1.f / mean3d.z;
    float max_Jxz = -(-lim_x_neg - cx) * rz;
    float min_Jxz = -(lim_x_pos - cx) * rz;
    float max_Jyz = -(-lim_x_neg - cy) * rz;
    float min_Jyz = -(lim_x_pos - cy) * rz;
    J[0][2] = clamp(J[0][2], min_Jxz, max_Jxz);
    J[1][2] = clamp(J[1][2], min_Jyz, max_Jyz);

    cov2d = mul(mul(J, cov3d), transpose(J));
    return true;
}

[CudaDeviceExport]
bool persp_proj_3dgs_nav(
    const float3 mean3d,
    const float3x3 cov3d,
    const float4 intrins,
    const float[10] dist_coeffs,
    const uint width,
    const uint height,
    out float2x2 cov2d,
    out float2 mean2d
) {
    return persp_proj_3dgs<false>(mean3d, cov3d, intrins, dist_coeffs, width, height, cov2d, mean2d);
}

[ForceInline]
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

[ForceInline]
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

[ForceInline]
[Differentiable]
bool equisolid_proj_3dgs<let assume_valid : bool>(
    // inputs
    const float3 mean3d,
    const float3x3 cov3d,
    const float4 intrins,
    const float[10] dist_coeffs,
    // outputs
    out float2x2 cov2d,
    out float2 mean2d
) {
    bool valid = equisolid_proj<assume_valid>(mean3d, intrins, dist_coeffs, mean2d);
    if (!assume_valid && !valid)
        return false;

    float2x3 J = equisolid_proj_jac(mean3d, intrins, dist_coeffs);
    cov2d = mul(mul(J, cov3d), transpose(J));
    return true;
}

[CudaDeviceExport]
bool fisheye_proj_3dgs_nav(
    const float3 mean3d,
    const float3x3 cov3d,
    const float4 intrins,
    const float[10] dist_coeffs,
    out float2x2 cov2d,
    out float2 mean2d
) {
    return fisheye_proj_3dgs<false>(mean3d, cov3d, intrins, dist_coeffs, cov2d, mean2d);
}

[CudaDeviceExport]
bool equisolid_proj_3dgs_nav(
    const float3 mean3d,
    const float3x3 cov3d,
    const float4 intrins,
    const float[10] dist_coeffs,
    out float2x2 cov2d,
    out float2 mean2d
) {
    return equisolid_proj_3dgs<false>(mean3d, cov3d, intrins, dist_coeffs, cov2d, mean2d);
}

[ForceInline]
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

    // We previously did this in rasterization with dot(raydir, pos) < 0, but it slows down rasterization
    // Here we try a different heuristics
    if (!assume_valid) {
        float max_z = sigmas.p[0].z;
        float2 min_p = proj_points[0], max_p = proj_points[0];
        [ForceUnroll] for (int i = 1; i < 7; i++) {
            max_z = max(max_z, sigmas.p[i].z);
            float2 p = proj_points[0] + (proj_points[i] - proj_points[0]) * 3.329;
            min_p = min(min_p, p);
            max_p = max(max_p, p);
        }
        min_p -= intrins.zw;
        max_p -= intrins.zw;
        if (max_z <= 0.0 && min_p.x * max_p.x < 0.0 && min_p.y * max_p.y < 0.0)
            return false;
    }

    [ForceUnroll] for (int i = 0; i < 7; i++)
        mean2d += sigmas.w_mean[i] * proj_points[i];

    [ForceUnroll] for (int i = 0; i < 7; i++) {
        float2 d = proj_points[i] - mean2d;
        cov2d += sigmas.w_cov[i] * float2x2(d.x * d.x, d.x * d.y, d.x * d.y, d.y * d.y);
    }
    return true;
}

[ForceInline]
bool equisolid_proj_3dgs_ut<let assume_valid : bool>(
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
        all_valid &= equisolid_proj<assume_valid>(sigmas.p[i], intrins, dist_coeffs, proj_points[i]);
    if (!assume_valid && !all_valid)
        return false;

    // We previously did this in rasterization with dot(raydir, pos) < 0, but it slows down rasterization
    // Here we try a different heuristics
    if (!assume_valid) {
        float max_z = sigmas.p[0].z;
        float2 min_p = proj_points[0], max_p = proj_points[0];
        [ForceUnroll] for (int i = 1; i < 7; i++) {
            max_z = max(max_z, sigmas.p[i].z);
            float2 p = proj_points[0] + (proj_points[i] - proj_points[0]) * 3.329;
            min_p = min(min_p, p);
            max_p = max(max_p, p);
        }
        min_p -= intrins.zw;
        max_p -= intrins.zw;
        if (max_z <= 0.0 && min_p.x * max_p.x < 0.0 && min_p.y * max_p.y < 0.0)
            return false;
    }

    [ForceUnroll] for (int i = 0; i < 7; i++)
        mean2d += sigmas.w_mean[i] * proj_points[i];

    [ForceUnroll] for (int i = 0; i < 7; i++) {
        float2 d = proj_points[i] - mean2d;
        cov2d += sigmas.w_cov[i] * float2x2(d.x * d.x, d.x * d.y, d.x * d.y, d.y * d.y);
    }
    return true;
}



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
    const float3 mean, const float4 quat, const float3 scale, const float in_opacity, const float3 sh_coeffs[16], \
    const float3x3 R, const float3 t, \
    no_diff const float fx, no_diff const float fy,  no_diff const float cx, no_diff const float cy, \
    no_diff const float[10] dist_coeffs, const uint image_width, const uint image_height, \
    no_diff out float4 aabb_xyxy, no_diff out float sorting_depth, out float2 mean2d, out float depth, out float3 conic, out float opacity, out float3 rgb

[ForceInline]
[Differentiable]
void projection_3dgs<let camera_model: CameraModelType, let assume_valid: bool, let eval3d: bool>(
    _projection_3dgs_ARGS
) {

    // transform Gaussian center to camera space
    float3 mean_c;
    posW2C(R, t, mean, mean_c);
    float cull_z = (camera_model == CameraModelType::PINHOLE) ?
        mean_c.z : length(mean_c);
    depth = mean_c.z;
    if (!assume_valid && cull_z <= 0.0f) {
        aabb_xyxy = float4(0.0f);
        return;
    }

    float4 intrins = { fx, fy, cx, cy };
    float2x2 covar2d;
    opacity = 1.0 / (1.0 + exp(-in_opacity));

    bool is_valid = true;
    if (eval3d) {
        SigmaPoints sigmas = world_gaussian_sigma_points(mean, exp(scale), normalize(quat), R, t);
        // SigmaPoints sigmas = world_gaussian_sigma_points(mean, exp(scale), quat, R, t);

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
        else if (camera_model == CameraModelType::EQUISOLID)
            is_valid &= equisolid_proj_3dgs_ut<assume_valid>(
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
            quat_scale_to_covar(normalize(quat), exp(scale), covar);
            // quat_scale_to_covar(quat, exp(scale), covar);  // quat is already normalized but can affect gradient
        }
        float3x3 covar_c;
        covarW2C(R, covar, covar_c);

        if (camera_model == CameraModelType::PINHOLE)
            is_valid &= persp_proj_3dgs<assume_valid>(
                mean_c, covar_c,
                intrins, dist_coeffs, image_width, image_height,
                covar2d, mean2d
            );
        else if (camera_model == CameraModelType::FISHEYE)
            is_valid &= fisheye_proj_3dgs<assume_valid>(
                mean_c, covar_c,
                intrins, dist_coeffs,
                covar2d, mean2d
            );
        else if (camera_model == CameraModelType::EQUISOLID)
            is_valid &= equisolid_proj_3dgs<assume_valid>(
                mean_c, covar_c,
                intrins, dist_coeffs,
                covar2d, mean2d
            );

        // do this since sorting for >180deg fisheye is weird
        if (!assume_valid && (camera_model == CameraModelType::FISHEYE || camera_model == CameraModelType::EQUISOLID)) {
            float2 mean2d_c = mean2d - float2(cx, cy);
            float opac = opacity * exp(-0.5f * dot(mul(inverse(covar2d), mean2d_c), mean2d_c));
            if (mean_c.z < 0.0 && opac > 1.0 / 255.0)
                is_valid = false;
        }
    }
    if (!assume_valid && !is_valid) {
        aabb_xyxy = float4(0.0f);
        return;
    }

    float compensation;
    static const float eps2d = 0.3;
    float det = add_blur(eps2d, covar2d, compensation);
    if (!assume_valid && det <= 0.f) {
        aabb_xyxy = float4(0.0f);
        return;
    }

    // compute the inverse of the 2d covariance
    float2x2 covar2d_inv = inverse(covar2d);

    float extend = 3.33f;
    {
        if (antialiased)
            opacity *= compensation;
        if (!assume_valid && opacity < ALPHA_THRESHOLD) {
            aabb_xyxy = float4(0.0f);
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
    float xmin = detach<float>(mean2d.x - radius_x);
    float xmax = detach<float>(mean2d.x + radius_x);
    float ymin = detach<float>(mean2d.y - radius_y);
    float ymax = detach<float>(mean2d.y + radius_y);

    // mask out gaussians outside the image region
    if (!assume_valid && (
        xmax <= 0 || xmin >= image_width ||
        ymax <= 0 || ymin >= image_height
    )) {
        aabb_xyxy = float4(0.0f);
        return;
    }

    aabb_xyxy = float4(xmin, ymin, xmax, ymax);
    sorting_depth = no_diff get_sorting_depth(camera_model, mean_c);
    // conic = eval3d ? exp(-scale) :
    //     float3(covar2d_inv[0][0], covar2d_inv[0][1], covar2d_inv[1][1]);
    conic = float3(covar2d_inv[0][0], covar2d_inv[0][1], covar2d_inv[1][1]);

    float3 campos = -mul(transpose(R), t);
    sh_coeffs_to_color(mean - campos, sh_coeffs, rgb);
    rgb = max(rgb+0.5f, 0.0f);
}

/* Forward Projection */

#define _DEF_projection_3dgs(funname, ...) \
    [CudaDeviceExport] \
    [ForceInline] \
    [Differentiable] \
    void funname(_projection_3dgs_ARGS) { \
        return projection_3dgs<__VA_ARGS__>( \
            antialiased, \
            mean, quat, scale, in_opacity, sh_coeffs, R, t, \
            fx, fy, cx, cy, dist_coeffs, \
            image_width, image_height, \
            aabb_xyxy, sorting_depth, mean2d, depth, conic, opacity, rgb \
        ); \
    }

_DEF_projection_3dgs(projection_3dgs_persp, CameraModelType::PINHOLE, false, false)
_DEF_projection_3dgs(projection_3dgs_fisheye, CameraModelType::FISHEYE, false, false)
_DEF_projection_3dgs(projection_3dgs_equisolid, CameraModelType::EQUISOLID, false, false)
_DEF_projection_3dgs(projection_3dgut_persp, CameraModelType::PINHOLE, false, true)
_DEF_projection_3dgs(projection_3dgut_fisheye, CameraModelType::FISHEYE, false, true)
_DEF_projection_3dgs(projection_3dgut_equisolid, CameraModelType::EQUISOLID, false, true)


/* Backward Projection */

interface _DiffProjection3DGS
{
    [ForceInline]
    [Differentiable]
    static void projection(_projection_3dgs_ARGS);
};

#define _DEF_projection_3dgs_interface(_InterfaceName, ...) \
    struct _InterfaceName : _DiffProjection3DGS { \
        [ForceInline] \
        [Differentiable] \
        static void projection(_projection_3dgs_ARGS) { \
            return projection_3dgs<__VA_ARGS__>( \
                antialiased, \
                mean, quat, scale, in_opacity, sh_coeffs, R, t, \
                fx, fy, cx, cy, dist_coeffs, \
                image_width, image_height, \
                aabb_xyxy, sorting_depth, mean2d, depth, conic, opacity, rgb \
            ); \
        } \
    };

_DEF_projection_3dgs_interface(_DiffProjection3DGS_3dgs_persp, CameraModelType::PINHOLE, true, false)
_DEF_projection_3dgs_interface(_DiffProjection3DGS_3dgs_fisheye, CameraModelType::FISHEYE, true, false)
_DEF_projection_3dgs_interface(_DiffProjection3DGS_3dgs_equisolid, CameraModelType::EQUISOLID, true, false)
_DEF_projection_3dgs_interface(_DiffProjection3DGS_3dgut_persp, CameraModelType::PINHOLE, true, true)
_DEF_projection_3dgs_interface(_DiffProjection3DGS_3dgut_fisheye, CameraModelType::FISHEYE, true, true)
_DEF_projection_3dgs_interface(_DiffProjection3DGS_3dgut_equisolid, CameraModelType::EQUISOLID, true, true)


#define _projection_3dgs_vjp_ARGS \
    const bool antialiased, \
    const float3 mean, const float4 quat, const float3 scale, const float in_opacity, const float3 sh_coeffs[16], \
    const float3x3 R, const float3 t, const float fx, const float fy, const float cx, const float cy, \
    const float[10] dist_coeffs, const uint image_width, const uint image_height, \
    const float2 v_mean2d, const float v_depth, const float3 v_conic, const float v_opacity, const float3 v_rgb, \
    out float3 v_mean, out float4 v_quat, out float3 v_scale, out float v_in_opacity, out float3 v_sh_coeffs[16], \
    out float3x3 v_R, out float3 v_t

[ForceInline]
void _projection_3dgs_vjp(
    _DiffProjection3DGS projectionInterface,
    _projection_3dgs_vjp_ARGS
) {
    DifferentialPair<float3> dp_mean = diffPair(mean);
    DifferentialPair<float4> dp_quat = diffPair(quat);
    DifferentialPair<float3> dp_scale = diffPair(scale);
    DifferentialPair<float> dp_in_opacity = diffPair(in_opacity);
    DifferentialPair<float3[16]> dp_sh_coeffs = diffPair(sh_coeffs);
    DifferentialPair<float3x3> dp_R = diffPair(R);
    DifferentialPair<float3> dp_t = diffPair(t);
    bwd_diff(projectionInterface.projection)(
        antialiased,
        dp_mean, dp_quat, dp_scale, dp_in_opacity, dp_sh_coeffs, dp_R, dp_t,
        fx, fy, cx, cy, dist_coeffs,
        image_width, image_height,
        v_mean2d, v_depth, v_conic, v_opacity, v_rgb
    );
    v_mean = dp_mean.d;
    v_quat = dp_quat.d;
    v_scale = dp_scale.d;
    v_in_opacity = dp_in_opacity.d;
    v_sh_coeffs = dp_sh_coeffs.d;
    v_R = dp_R.d;
    v_t = dp_t.d;
}

#define _DEF_projection_3dgs_vjp(funname, _InterfaceName, ...) \
    [CudaDeviceExport] \
    void funname( \
        _projection_3dgs_vjp_ARGS \
    ) { \
        return _projection_3dgs_vjp( \
            _InterfaceName(), \
            antialiased, \
            mean, quat, scale, in_opacity, sh_coeffs, R, t, \
            fx, fy, cx, cy, dist_coeffs, \
            image_width, image_height, \
            v_mean2d, v_depth, v_conic, v_opacity, v_rgb, \
            v_mean, v_quat, v_scale, v_in_opacity, v_sh_coeffs, v_R, v_t \
        ); \
    }

_DEF_projection_3dgs_vjp(projection_3dgs_persp_vjp, _DiffProjection3DGS_3dgs_persp)
_DEF_projection_3dgs_vjp(projection_3dgs_fisheye_vjp, _DiffProjection3DGS_3dgs_fisheye)
_DEF_projection_3dgs_vjp(projection_3dgs_equisolid_vjp, _DiffProjection3DGS_3dgs_equisolid)
_DEF_projection_3dgs_vjp(projection_3dgut_persp_vjp, _DiffProjection3DGS_3dgut_persp)
_DEF_projection_3dgs_vjp(projection_3dgut_fisheye_vjp, _DiffProjection3DGS_3dgut_fisheye)
_DEF_projection_3dgs_vjp(projection_3dgut_equisolid_vjp, _DiffProjection3DGS_3dgut_equisolid)


/* Backward Projection with Vector-Residual Product and Hessian Diagonal to Mean */

#undef _projection_3dgs_vjp_ARGS
#define _projection_3dgs_vjp_ARGS \
    const bool antialiased, \
    const float3 mean, const float4 quat, const float3 scale, const float in_opacity, const float3 sh_coeffs[16], \
    const float3x3 R, const float3 t, const float fx, const float fy, const float cx, const float cy, \
    const float[10] dist_coeffs, const uint image_width, const uint image_height, \
    const float2 v_mean2d, const float v_depth, const float3 v_conic, const float v_opacity, const float3 v_rgb, \
    const float2 vr_mean2d, const float vr_depth, const float3 vr_conic, const float vr_opacity, const float3 vr_rgb, \
    const float2 h_mean2d, const float h_depth, const float3 h_conic, const float h_opacity, const float3 h_rgb, \
    out float3 v_mean, out float4 v_quat, out float3 v_scale, out float v_in_opacity, out float3 v_sh_coeffs[16], \
    out float3x3 v_R, out float3 v_t, \
    out float3 vr_mean, out float3 h_mean

[ForceInline]
void _projection_3dgs_vjp(
    _DiffProjection3DGS projectionInterface,
    _projection_3dgs_vjp_ARGS
) {
    // gradient
    {
        DifferentialPair<float3> dp_mean = diffPair(mean);
        DifferentialPair<float4> dp_quat = diffPair(quat);
        DifferentialPair<float3> dp_scale = diffPair(scale);
        DifferentialPair<float> dp_in_opacity = diffPair(in_opacity);
        DifferentialPair<float3[16]> dp_sh_coeffs = diffPair(sh_coeffs);
        DifferentialPair<float3x3> dp_R = diffPair(R);
        DifferentialPair<float3> dp_t = diffPair(t);
        bwd_diff(projectionInterface.projection)(
            antialiased,
            dp_mean, dp_quat, dp_scale, dp_in_opacity, dp_sh_coeffs, dp_R, dp_t,
            fx, fy, cx, cy, dist_coeffs,
            image_width, image_height,
            v_mean2d, v_depth, v_conic, v_opacity, v_rgb
        );
        v_mean = dp_mean.d;
        v_quat = dp_quat.d;
        v_scale = dp_scale.d;
        v_in_opacity = dp_in_opacity.d;
        v_sh_coeffs = dp_sh_coeffs.d;
        v_R = dp_R.d;
        v_t = dp_t.d;
    }

    // jacobian residual product
    {
        DifferentialPair<float3> dp_mean = diffPair(mean);
        DifferentialPair<float4> dp_quat = diffPair(quat);
        DifferentialPair<float3> dp_scale = diffPair(scale);
        DifferentialPair<float> dp_in_opacity = diffPair(in_opacity);
        DifferentialPair<float3[16]> dp_sh_coeffs = diffPair(sh_coeffs);
        DifferentialPair<float3x3> dp_R = diffPair(R);
        DifferentialPair<float3> dp_t = diffPair(t);
        bwd_diff(projectionInterface.projection)(
            antialiased,
            dp_mean, dp_quat, dp_scale, dp_in_opacity, dp_sh_coeffs, dp_R, dp_t,
            fx, fy, cx, cy, dist_coeffs,
            image_width, image_height,
            vr_mean2d, vr_depth, vr_conic, vr_opacity, vr_rgb
        );
        vr_mean = dp_mean.d;
    }

    // hessian diagonal (mean)
    [ForceUnroll]
    for (int i = 0; i < 3; ++i) {
        DifferentialPair<float3> dp_mean = diffPair(
            mean, float3(float(i == 0), float(i == 1), float(i == 2)) );
        DifferentialPair<float4> dp_quat = diffPair(quat, float4(0));
        DifferentialPair<float3> dp_scale = diffPair(scale, float3(0));
        DifferentialPair<float> dp_in_opacity = diffPair(in_opacity, 0.0f);
        DifferentialPair<float3[16]> dp_sh_coeffs = diffPair(sh_coeffs, float3[16]( { float3(0) }));
        DifferentialPair<float3x3> dp_R = diffPair(R, float3x3(0));
        DifferentialPair<float3> dp_t = diffPair(t, float3(0));
        DifferentialPair<float2> dp_mean2d;
        DifferentialPair<float> dp_depth;
        DifferentialPair<float3> dp_conic;
        DifferentialPair<float> dp_opacity;
        DifferentialPair<float3> dp_rgb;
        float4 aabb_xyxy;
        float sorting_depth;
        fwd_diff(projectionInterface.projection)(
            antialiased,
            dp_mean, dp_quat, dp_scale, dp_in_opacity, dp_sh_coeffs, dp_R, dp_t,
            fx, fy, cx, cy, dist_coeffs,
            image_width, image_height,
            aabb_xyxy, sorting_depth, dp_mean2d, dp_depth, dp_conic, dp_opacity, dp_rgb
        );
        h_mean[i] =
            dot(h_mean2d, dp_mean2d.d * dp_mean2d.d) +
            h_depth * dp_depth.d * dp_depth.d +
            dot(h_conic, dp_conic.d * dp_conic.d) +
            h_opacity * dp_opacity.d * dp_opacity.d +
            dot(h_rgb, dp_rgb.d * dp_rgb.d);
    }
}

#undef _DEF_projection_3dgs_vjp
#define _DEF_projection_3dgs_vjp(funname, _InterfaceName, ...) \
    [CudaDeviceExport] \
    void funname( \
        _projection_3dgs_vjp_ARGS \
    ) { \
        return _projection_3dgs_vjp( \
            _InterfaceName(), \
            antialiased, \
            mean, quat, scale, in_opacity, sh_coeffs, R, t, \
            fx, fy, cx, cy, dist_coeffs, \
            image_width, image_height, \
            v_mean2d, v_depth, v_conic, v_opacity, v_rgb, \
            vr_mean2d, vr_depth, vr_conic, vr_opacity, vr_rgb, \
            h_mean2d, h_depth, h_conic, h_opacity, h_rgb, \
            v_mean, v_quat, v_scale, v_in_opacity, v_sh_coeffs, v_R, v_t, \
            vr_mean, h_mean \
        ); \
    }

_DEF_projection_3dgs_vjp(projection_3dgs_persp_vjp, _DiffProjection3DGS_3dgs_persp)
_DEF_projection_3dgs_vjp(projection_3dgs_fisheye_vjp, _DiffProjection3DGS_3dgs_fisheye)
_DEF_projection_3dgs_vjp(projection_3dgs_equisolid_vjp, _DiffProjection3DGS_3dgs_equisolid)
_DEF_projection_3dgs_vjp(projection_3dgut_persp_vjp, _DiffProjection3DGS_3dgut_persp)
_DEF_projection_3dgs_vjp(projection_3dgut_fisheye_vjp, _DiffProjection3DGS_3dgut_fisheye)
_DEF_projection_3dgs_vjp(projection_3dgut_equisolid_vjp, _DiffProjection3DGS_3dgut_equisolid)


/* Backward Projection with Vector-Residual Product and Hessian Diagonal to All Parameters except SH */

#undef _projection_3dgs_vjp_ARGS
#define _projection_3dgs_vjp_ARGS \
    const bool antialiased, \
    const float3 mean, const float4 quat, const float3 scale, const float in_opacity, const float3 sh_coeffs[16], \
    const float3x3 R, const float3 t, const float fx, const float fy, const float cx, const float cy, \
    const float[10] dist_coeffs, const uint image_width, const uint image_height, \
    const float2 v_mean2d, const float v_depth, const float3 v_conic, const float v_opacity, const float3 v_rgb, \
    const float2 vr_mean2d, const float vr_depth, const float3 vr_conic, const float vr_opacity, const float3 vr_rgb, \
    const float2 h_mean2d, const float h_depth, const float3 h_conic, const float h_opacity, const float3 h_rgb, \
    out float3 v_mean, out float4 v_quat, out float3 v_scale, out float v_in_opacity, out float3 v_sh_coeffs[16], \
    out float3x3 v_R, out float3 v_t, \
    out float3 vr_mean, out float4 vr_quat, out float3 vr_scale, out float vr_in_opacity, out float3 vr_features_dc, \
    out float3 h_mean, out float4 h_quat, out float3 h_scale, out float h_in_opacity, out float3 h_features_dc

[ForceInline]
void _projection_3dgs_vjp(
    _DiffProjection3DGS projectionInterface,
    _projection_3dgs_vjp_ARGS
) {
    // gradient
    {
        DifferentialPair<float3> dp_mean = diffPair(mean);
        DifferentialPair<float4> dp_quat = diffPair(quat);
        DifferentialPair<float3> dp_scale = diffPair(scale);
        DifferentialPair<float> dp_in_opacity = diffPair(in_opacity);
        DifferentialPair<float3[16]> dp_sh_coeffs = diffPair(sh_coeffs);
        DifferentialPair<float3x3> dp_R = diffPair(R);
        DifferentialPair<float3> dp_t = diffPair(t);
        bwd_diff(projectionInterface.projection)(
            antialiased,
            dp_mean, dp_quat, dp_scale, dp_in_opacity, dp_sh_coeffs, dp_R, dp_t,
            fx, fy, cx, cy, dist_coeffs,
            image_width, image_height,
            v_mean2d, v_depth, v_conic, v_opacity, v_rgb
        );
        v_mean = dp_mean.d;
        v_quat = dp_quat.d;
        v_scale = dp_scale.d;
        v_in_opacity = dp_in_opacity.d;
        v_sh_coeffs = dp_sh_coeffs.d;
        v_R = dp_R.d;
        v_t = dp_t.d;
    }

    // jacobian residual product
    {
        DifferentialPair<float3> dp_mean = diffPair(mean);
        DifferentialPair<float4> dp_quat = diffPair(quat);
        DifferentialPair<float3> dp_scale = diffPair(scale);
        DifferentialPair<float> dp_in_opacity = diffPair(in_opacity);
        DifferentialPair<float3[16]> dp_sh_coeffs = diffPair(sh_coeffs);
        DifferentialPair<float3x3> dp_R = diffPair(R);
        DifferentialPair<float3> dp_t = diffPair(t);
        bwd_diff(projectionInterface.projection)(
            antialiased,
            dp_mean, dp_quat, dp_scale, dp_in_opacity, dp_sh_coeffs, dp_R, dp_t,
            fx, fy, cx, cy, dist_coeffs,
            image_width, image_height,
            vr_mean2d, vr_depth, vr_conic, vr_opacity, vr_rgb
        );
        vr_mean = dp_mean.d;
        vr_quat = dp_quat.d;
        vr_scale = dp_scale.d;
        vr_in_opacity = dp_in_opacity.d;
        vr_features_dc = dp_sh_coeffs.d[0];
    }

    // TODO: some way to refactor this

    // hessian diagonal (mean)
    [ForceUnroll]
    for (int i = 0; i < 3; ++i) {
        DifferentialPair<float3> dp_mean = diffPair(
            mean, float3(float(i == 0), float(i == 1), float(i == 2)) );
        DifferentialPair<float4> dp_quat = diffPair(quat, float4(0));
        DifferentialPair<float3> dp_scale = diffPair(scale, float3(0));
        DifferentialPair<float> dp_in_opacity = diffPair(in_opacity, 0.0f);
        DifferentialPair<float3[16]> dp_sh_coeffs = diffPair(sh_coeffs, float3[16]( { float3(0) }));
        DifferentialPair<float3x3> dp_R = diffPair(R, float3x3(0));
        DifferentialPair<float3> dp_t = diffPair(t, float3(0));
        DifferentialPair<float2> dp_mean2d;
        DifferentialPair<float> dp_depth;
        DifferentialPair<float3> dp_conic;
        DifferentialPair<float> dp_opacity;
        DifferentialPair<float3> dp_rgb;
        float4 aabb_xyxy;
        float sorting_depth;
        fwd_diff(projectionInterface.projection)(
            antialiased,
            dp_mean, dp_quat, dp_scale, dp_in_opacity, dp_sh_coeffs, dp_R, dp_t,
            fx, fy, cx, cy, dist_coeffs,
            image_width, image_height,
            aabb_xyxy, sorting_depth, dp_mean2d, dp_depth, dp_conic, dp_opacity, dp_rgb
        );
        h_mean[i] =
            dot(h_mean2d, dp_mean2d.d * dp_mean2d.d) +
            h_depth * dp_depth.d * dp_depth.d +
            dot(h_conic, dp_conic.d * dp_conic.d) +
            h_opacity * dp_opacity.d * dp_opacity.d +
            dot(h_rgb, dp_rgb.d * dp_rgb.d);
    }

    // hessian diagonal (quat)
    [ForceUnroll]
    for (int i = 0; i < 4; ++i) {
        DifferentialPair<float3> dp_mean = diffPair(mean, float3(0));
        DifferentialPair<float4> dp_quat = diffPair(
            quat, float4(float(i == 0), float(i == 1), float(i == 2), float(i == 3)));
        DifferentialPair<float3> dp_scale = diffPair(scale, float3(0));
        DifferentialPair<float> dp_in_opacity = diffPair(in_opacity, 0.0f);
        DifferentialPair<float3[16]> dp_sh_coeffs = diffPair(sh_coeffs, float3[16]( { float3(0) }));
        DifferentialPair<float3x3> dp_R = diffPair(R, float3x3(0));
        DifferentialPair<float3> dp_t = diffPair(t, float3(0));
        DifferentialPair<float2> dp_mean2d;
        DifferentialPair<float> dp_depth;
        DifferentialPair<float3> dp_conic;
        DifferentialPair<float> dp_opacity;
        DifferentialPair<float3> dp_rgb;
        float4 aabb_xyxy;
        float sorting_depth;
        fwd_diff(projectionInterface.projection)(
            antialiased,
            dp_mean, dp_quat, dp_scale, dp_in_opacity, dp_sh_coeffs, dp_R, dp_t,
            fx, fy, cx, cy, dist_coeffs,
            image_width, image_height,
            aabb_xyxy, sorting_depth, dp_mean2d, dp_depth, dp_conic, dp_opacity, dp_rgb
        );
        h_quat[i] =
            dot(h_mean2d, dp_mean2d.d * dp_mean2d.d) +
            h_depth * dp_depth.d * dp_depth.d +
            dot(h_conic, dp_conic.d * dp_conic.d) +
            h_opacity * dp_opacity.d * dp_opacity.d +
            dot(h_rgb, dp_rgb.d * dp_rgb.d);
    }

    // hessian diagonal (scale)
    [ForceUnroll]
    for (int i = 0; i < 3; ++i) {
        DifferentialPair<float3> dp_mean = diffPair(mean, float3(0));
        DifferentialPair<float4> dp_quat = diffPair(quat, float4(0));
        DifferentialPair<float3> dp_scale = diffPair(
            scale, float3(float(i == 0), float(i == 1), float(i == 2)));
        DifferentialPair<float> dp_in_opacity = diffPair(in_opacity, 0.0f);
        DifferentialPair<float3[16]> dp_sh_coeffs = diffPair(sh_coeffs, float3[16]( { float3(0) }));
        DifferentialPair<float3x3> dp_R = diffPair(R, float3x3(0));
        DifferentialPair<float3> dp_t = diffPair(t, float3(0));
        DifferentialPair<float2> dp_mean2d;
        DifferentialPair<float> dp_depth;
        DifferentialPair<float3> dp_conic;
        DifferentialPair<float> dp_opacity;
        DifferentialPair<float3> dp_rgb;
        float4 aabb_xyxy;
        float sorting_depth;
        fwd_diff(projectionInterface.projection)(
            antialiased,
            dp_mean, dp_quat, dp_scale, dp_in_opacity, dp_sh_coeffs, dp_R, dp_t,
            fx, fy, cx, cy, dist_coeffs,
            image_width, image_height,
            aabb_xyxy, sorting_depth, dp_mean2d, dp_depth, dp_conic, dp_opacity, dp_rgb
        );
        h_scale[i] =
            dot(h_mean2d, dp_mean2d.d * dp_mean2d.d) +
            h_depth * dp_depth.d * dp_depth.d +
            dot(h_conic, dp_conic.d * dp_conic.d) +
            h_opacity * dp_opacity.d * dp_opacity.d +
            dot(h_rgb, dp_rgb.d * dp_rgb.d);
    }

    // hessian diagonal (in_opacity)
    {
        DifferentialPair<float3> dp_mean = diffPair(mean, float3(0));
        DifferentialPair<float4> dp_quat = diffPair(quat, float4(0));
        DifferentialPair<float3> dp_scale = diffPair(scale, float3(0));
        DifferentialPair<float> dp_in_opacity = diffPair(in_opacity, 1.0f);
        DifferentialPair<float3[16]> dp_sh_coeffs = diffPair(sh_coeffs, float3[16]( { float3(0) }));
        DifferentialPair<float3x3> dp_R = diffPair(R, float3x3(0));
        DifferentialPair<float3> dp_t = diffPair(t, float3(0));
        DifferentialPair<float2> dp_mean2d;
        DifferentialPair<float> dp_depth;
        DifferentialPair<float3> dp_conic;
        DifferentialPair<float> dp_opacity;
        DifferentialPair<float3> dp_rgb;
        float4 aabb_xyxy;
        float sorting_depth;
        fwd_diff(projectionInterface.projection)(
            antialiased,
            dp_mean, dp_quat, dp_scale, dp_in_opacity, dp_sh_coeffs, dp_R, dp_t,
            fx, fy, cx, cy, dist_coeffs,
            image_width, image_height,
            aabb_xyxy, sorting_depth, dp_mean2d, dp_depth, dp_conic, dp_opacity, dp_rgb
        );
        h_in_opacity =
            dot(h_mean2d, dp_mean2d.d * dp_mean2d.d) +
            h_depth * dp_depth.d * dp_depth.d +
            dot(h_conic, dp_conic.d * dp_conic.d) +
            h_opacity * dp_opacity.d * dp_opacity.d +
            dot(h_rgb, dp_rgb.d * dp_rgb.d);
    }

    // hessian diagonal (features_dc)
    [ForceUnroll]
    for (int i = 0; i < 3; ++i) {
        DifferentialPair<float3> dp_mean = diffPair(mean, float3(0));
        DifferentialPair<float4> dp_quat = diffPair(quat, float4(0));
        DifferentialPair<float3> dp_scale = diffPair(scale, float3(0));
        DifferentialPair<float> dp_in_opacity = diffPair(in_opacity, 1.0f);
        float3[16] d_sh_coeffs = float3[16]( { float3(0) });
        d_sh_coeffs[0] = float3(float(i == 0), float(i == 1), float(i == 2));
        DifferentialPair<float3[16]> dp_sh_coeffs = diffPair(sh_coeffs, d_sh_coeffs);
        DifferentialPair<float3x3> dp_R = diffPair(R, float3x3(0));
        DifferentialPair<float3> dp_t = diffPair(t, float3(0));
        DifferentialPair<float2> dp_mean2d;
        DifferentialPair<float> dp_depth;
        DifferentialPair<float3> dp_conic;
        DifferentialPair<float> dp_opacity;
        DifferentialPair<float3> dp_rgb;
        float4 aabb_xyxy;
        float sorting_depth;
        fwd_diff(projectionInterface.projection)(
            antialiased,
            dp_mean, dp_quat, dp_scale, dp_in_opacity, dp_sh_coeffs, dp_R, dp_t,
            fx, fy, cx, cy, dist_coeffs,
            image_width, image_height,
            aabb_xyxy, sorting_depth, dp_mean2d, dp_depth, dp_conic, dp_opacity, dp_rgb
        );
        h_features_dc[i] =
            dot(h_mean2d, dp_mean2d.d * dp_mean2d.d) +
            h_depth * dp_depth.d * dp_depth.d +
            dot(h_conic, dp_conic.d * dp_conic.d) +
            h_opacity * dp_opacity.d * dp_opacity.d +
            dot(h_rgb, dp_rgb.d * dp_rgb.d);
    }
}

#undef _DEF_projection_3dgs_vjp
#define _DEF_projection_3dgs_vjp(funname, _InterfaceName, ...) \
    [CudaDeviceExport] \
    void funname( \
        _projection_3dgs_vjp_ARGS \
    ) { \
        return _projection_3dgs_vjp( \
            _InterfaceName(), \
            antialiased, \
            mean, quat, scale, in_opacity, sh_coeffs, R, t, \
            fx, fy, cx, cy, dist_coeffs, \
            image_width, image_height, \
            v_mean2d, v_depth, v_conic, v_opacity, v_rgb, \
            vr_mean2d, vr_depth, vr_conic, vr_opacity, vr_rgb, \
            h_mean2d, h_depth, h_conic, h_opacity, h_rgb, \
            v_mean, v_quat, v_scale, v_in_opacity, v_sh_coeffs, v_R, v_t, \
            vr_mean, vr_quat, vr_scale, vr_in_opacity, vr_features_dc, \
            h_mean, h_quat, h_scale, h_in_opacity, h_features_dc \
        ); \
    }

_DEF_projection_3dgs_vjp(projection_3dgs_persp_vjp, _DiffProjection3DGS_3dgs_persp)
_DEF_projection_3dgs_vjp(projection_3dgs_fisheye_vjp, _DiffProjection3DGS_3dgs_fisheye)
_DEF_projection_3dgs_vjp(projection_3dgs_equisolid_vjp, _DiffProjection3DGS_3dgs_equisolid)
_DEF_projection_3dgs_vjp(projection_3dgut_persp_vjp, _DiffProjection3DGS_3dgut_persp)
_DEF_projection_3dgs_vjp(projection_3dgut_fisheye_vjp, _DiffProjection3DGS_3dgut_fisheye)
_DEF_projection_3dgs_vjp(projection_3dgut_equisolid_vjp, _DiffProjection3DGS_3dgut_equisolid)

