#pragma once
#include "primitive_3dgs.slang"


#define _projection_3dgs_sv_ARGS(num_sv) \
    const bool antialiased, \
    const float3 mean, \
    const float4 quat, \
    const float3 scale, \
    const float in_opacity, \
    const float3 sv_sites[num_sv], \
    const float3 sv_colors[num_sv], \
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
void projection_3dgs_sv<let camera_model : CameraModelType, let assume_valid : bool, let eval3d : bool, let num_sv : int>(
    _projection_3dgs_sv_ARGS(num_sv)
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
    depth = get_sorting_depth(mean_c);
    conic = eval3d ? exp(-scale) :
        float3(covar2d_inv[0][0], covar2d_inv[0][1], covar2d_inv[1][1]);

    float3 campos = -mul(transpose(R), t);
    sv_coeffs_to_color<num_sv>(mean - campos, sv_sites, sv_colors, rgb);
    rgb = max(rgb+0.5f, 0.0f);
}


#define _DEF_projection_3dgs_sv(funname, num_sv, ...) \
    [CudaDeviceExport] \
    [Differentiable] \
    void funname(_projection_3dgs_sv_ARGS(num_sv)) { \
        return projection_3dgs_sv<__VA_ARGS__, num_sv>( \
            antialiased, \
            mean, quat, scale, in_opacity, sv_sites, sv_colors, R, t, \
            fx, fy, cx, cy, dist_coeffs, \
            image_width, image_height, near_plane, far_plane, \
            aabb_xyxy, mean2d, depth, conic, opacity, rgb \
        ); \
    }

// #define _DEF_projection_3dgs_sv_i(num_sv) \
//     _DEF_projection_3dgs_sv(projection_3dgs_sv_persp, num_sv, CameraModelType::PINHOLE, false, false) \
//     _DEF_projection_3dgs_sv(projection_3dgs_sv_fisheye, num_sv, CameraModelType::FISHEYE, false, false) \
//     _DEF_projection_3dgs_sv(projection_3dgs_sv_ortho, num_sv, CameraModelType::ORTHO, false, false) \
//     _DEF_projection_3dgs_sv(projection_3dgut_sv_persp, num_sv, CameraModelType::PINHOLE, false, true) \
//     _DEF_projection_3dgs_sv(projection_3dgut_sv_fisheye, num_sv, CameraModelType::FISHEYE, false, true)

#define _DEF_projection_3dgs_sv_i(num_sv) \
    _DEF_projection_3dgs_sv(projection_3dgut_sv_persp, num_sv, CameraModelType::PINHOLE, false, true) \
    _DEF_projection_3dgs_sv(projection_3dgut_sv_fisheye, num_sv, CameraModelType::FISHEYE, false, true)

// _DEF_projection_3dgs_sv_i(1);
_DEF_projection_3dgs_sv_i(2);
_DEF_projection_3dgs_sv_i(3);
_DEF_projection_3dgs_sv_i(4);
_DEF_projection_3dgs_sv_i(5);
_DEF_projection_3dgs_sv_i(6);
_DEF_projection_3dgs_sv_i(7);
_DEF_projection_3dgs_sv_i(8);


#define _projection_3dgs_sv_vjp_ARGS(num_sv) \
    const bool antialiased, \
    const float3 mean, \
    const float4 quat, \
    const float3 scale, \
    const float in_opacity, \
    const float3 sv_sites[num_sv], \
    const float3 sv_colors[num_sv], \
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
    out float3 v_sv_dirs[num_sv], \
    out float3 v_sv_colors[num_sv], \
    out float3x3 v_R, \
    out float3 v_t

#define _DEF_projection_3dgs_sv_vjp(funname, fwd_funname, num_sv, ...) \
    [CudaDeviceExport] \
    void funname( \
        _projection_3dgs_sv_vjp_ARGS(num_sv) \
    ) { \
        DifferentialPair<float3> dp_mean = diffPair(mean); \
        DifferentialPair<float4> dp_quat = diffPair(quat); \
        DifferentialPair<float3> dp_scale = diffPair(scale); \
        DifferentialPair<float> dp_in_opacity = diffPair(in_opacity); \
        DifferentialPair<float3[num_sv]> dp_sv_dirs = diffPair(sv_sites); \
        DifferentialPair<float3[num_sv]> dp_sv_colors = diffPair(sv_colors); \
        DifferentialPair<float3x3> dp_R = diffPair(R); \
        DifferentialPair<float3> dp_t = diffPair(t); \
        bwd_diff(fwd_funname)( \
            antialiased, \
            dp_mean, dp_quat, dp_scale, dp_in_opacity, dp_sv_dirs, dp_sv_colors, dp_R, dp_t, \
            fx, fy, cx, cy, dist_coeffs, \
            image_width, image_height, 0.0f, 1e30f, \
            v_mean2d, v_depth, v_conic, v_opacity, v_rgb \
        ); \
        v_mean = dp_mean.d; \
        v_quat = dp_quat.d; \
        v_scale = dp_scale.d; \
        v_in_opacity = dp_in_opacity.d; \
        v_sv_dirs = dp_sv_dirs.d; \
        v_sv_colors = dp_sv_colors.d; \
        v_R = dp_R.d; \
        v_t = dp_t.d; \
    }

// #define _DEF_projection_3dgs_sv_vjp_i(num_sv) \
//     _DEF_projection_3dgs_sv(_projection_3dgs_sv_persp_differentiable, num_sv, CameraModelType::PINHOLE, true, false) \
//     _DEF_projection_3dgs_sv(_projection_3dgs_sv_fisheye_differentiable, num_sv, CameraModelType::FISHEYE, true, false) \
//     _DEF_projection_3dgs_sv(_projection_3dgs_sv_ortho_differentiable, num_sv, CameraModelType::ORTHO, true, false) \
//     _DEF_projection_3dgs_sv(_projection_3dgut_sv_persp_differentiable, num_sv, CameraModelType::PINHOLE, true, true) \
//     _DEF_projection_3dgs_sv(_projection_3dgut_sv_fisheye_differentiable, num_sv, CameraModelType::FISHEYE, true, true) \
//     _DEF_projection_3dgs_sv_vjp(projection_3dgs_sv_persp_vjp, _projection_3dgs_sv_persp_differentiable, num_sv) \
//     _DEF_projection_3dgs_sv_vjp(projection_3dgs_sv_fisheye_vjp, _projection_3dgs_sv_fisheye_differentiable, num_sv) \
//     _DEF_projection_3dgs_sv_vjp(projection_3dgs_sv_ortho_vjp, _projection_3dgs_sv_ortho_differentiable, num_sv) \
//     _DEF_projection_3dgs_sv_vjp(projection_3dgut_sv_persp_vjp, _projection_3dgut_sv_persp_differentiable, num_sv) \
//     _DEF_projection_3dgs_sv_vjp(projection_3dgut_sv_fisheye_vjp, _projection_3dgut_sv_fisheye_differentiable, num_sv)

#define _DEF_projection_3dgs_sv_vjp_i(num_sv) \
    _DEF_projection_3dgs_sv(_projection_3dgut_sv_persp_differentiable, num_sv, CameraModelType::PINHOLE, true, true) \
    _DEF_projection_3dgs_sv(_projection_3dgut_sv_fisheye_differentiable, num_sv, CameraModelType::FISHEYE, true, true) \
    _DEF_projection_3dgs_sv_vjp(projection_3dgut_sv_persp_vjp, _projection_3dgut_sv_persp_differentiable, num_sv) \
    _DEF_projection_3dgs_sv_vjp(projection_3dgut_sv_fisheye_vjp, _projection_3dgut_sv_fisheye_differentiable, num_sv)

// _DEF_projection_3dgs_sv_vjp_i(1);
_DEF_projection_3dgs_sv_vjp_i(2);
_DEF_projection_3dgs_sv_vjp_i(3);
_DEF_projection_3dgs_sv_vjp_i(4);
_DEF_projection_3dgs_sv_vjp_i(5);
_DEF_projection_3dgs_sv_vjp_i(6);
_DEF_projection_3dgs_sv_vjp_i(7);
_DEF_projection_3dgs_sv_vjp_i(8);
