#include "projection_utils.slang"


    // const float3 vert0, \
    // const float3 vert1, \
    // const float3 vert2, \

#define _projection_opaque_triangle_ARGS \
    const float3 mean, \
    const float4 quat, \
    const float3 scale, \
    const float2 hardness, \
    const float3x3 R, \
    const float3 t, \
    no_diff const float fx, \
    no_diff const float fy, \
    no_diff const float cx, \
    no_diff const float cy, \
    no_diff const float4 radial_coeffs, \
    no_diff const float2 tangential_coeffs, \
    no_diff const float2 thin_prism_coeffs, \
    const uint image_width, \
    const uint image_height, \
    no_diff const float near_plane, \
    no_diff const float far_plane, \
    out int4 aabb_xyxy, \
    out float depth, \
    out float3 normal, \
    out float2 uv0, \
    out float2 uv1, \
    out float2 uv2, \
    out float2 out_hardness

[ForceInline]
[Differentiable]
void projection_opaque_triangle<let camera_model: CameraModelType, let assume_valid: bool>(
    _projection_opaque_triangle_ARGS
) {
    // float3 mean = (vert0+vert1+vert2)/3.0f;
    float3 mean_c;
    posW2C(R, t, mean, mean_c);
    if (!assume_valid && (mean_c.z < near_plane || mean_c.z > far_plane)) {
        aabb_xyxy = {0, 0, 0, 0};
        return;
    }

    // transform mean+scale+quat to vertices
  #if 0
    float3x3 M;
    quat_scale_to_sqrt_covar(quat, exp(scale), M);
    // float3 vert0 = mul(M, float3(1, 0, 0)) + mean;
    // float3 vert1 = mul(M, float3(-0.5f, sqrt(0.75f), 0)) + mean;
    // float3 vert2 = mul(M, float3(-0.5f, -sqrt(0.75f), 0)) + mean;
    float3 vert0 = mul(float3(1, 0, 0), M) + mean;
    float3 vert1 = mul(float3(-0.5f, sqrt(0.75f), 0), M) + mean;
    float3 vert2 = mul(float3(-0.5f, -sqrt(0.75f), 0), M) + mean;
  #else
    float sx = exp(scale.x),
        sy = exp(scale.y),
        sz = scale.z - 0.5*(scale.x+scale.y);
    float3x3 M = quat_to_rotmat(quat);
    float3 vert0 = mul(M, float3(sx, 0, 0)) + mean;
    float3 vert1 = mul(M, float3(sx*(-0.5+sz), sy, 0)) + mean;
    float3 vert2 = mul(M, float3(sx*(-0.5-sz), -sy, 0)) + mean;
  #endif

    // transform center to camera space
    float3 vert0_c, vert1_c, vert2_c;
    posW2C(R, t, vert0, vert0_c);
    posW2C(R, t, vert1, vert1_c);
    posW2C(R, t, vert2, vert2_c);
    if (!assume_valid && (
        vert0_c.z < near_plane || vert0_c.z > far_plane ||
        vert1_c.z < near_plane || vert1_c.z > far_plane ||
        vert2_c.z < near_plane || vert2_c.z > far_plane
    )) {
        aabb_xyxy = {0, 0, 0, 0};
        return;
    }

    if (camera_model == CameraModelType::PINHOLE) {
        uv0 = vert0_c.xy / vert0_c.z;
        uv1 = vert1_c.xy / vert1_c.z;
        uv2 = vert2_c.xy / vert2_c.z;
        uv0 = float2(fx, fy) * uv0 + float2(cx, cy);
        uv1 = float2(fx, fy) * uv1 + float2(cx, cy);
        uv2 = float2(fx, fy) * uv2 + float2(cx, cy);
    }
    else if (camera_model == CameraModelType::FISHEYE) {
        float4 intrins = float4(fx, fy, cx, cy);
        CameraDistortion dist_coeffs = {radial_coeffs, tangential_coeffs, thin_prism_coeffs};
        // TODO: some approximation better than this
        uv0 = fisheye_proj(vert0_c, intrins, dist_coeffs);
        uv1 = fisheye_proj(vert1_c, intrins, dist_coeffs);
        uv2 = fisheye_proj(vert2_c, intrins, dist_coeffs);
    }
    else if (camera_model == CameraModelType::ORTHO)
        {}  // TODO

    // compute tight rectangular bounding box (non differentiable)

    float2 e0 = uv1-uv0, e1 = uv2-uv1, e2 = uv0-uv2;
    float dmax = abs(det(e0, e1)) / (length(e0) + length(e1) + length(e2));

    float offset = (1.0f / (1.0f - exp2(-1.0f/(1.0f-hardness.y))) - 1.0f) * dmax;

    float xmax = max(max(uv0.x, uv1.x), uv2.x) + offset;
    float xmin = min(min(uv0.x, uv1.x), uv2.x) - offset;
    float ymax = max(max(uv0.y, uv1.y), uv2.y) + offset;
    float ymin = min(min(uv0.y, uv1.y), uv2.y) - offset;

    // mask out splats outside the image region
    if (!assume_valid && (
        xmax <= 0 || xmin >= image_width ||
        ymax <= 0 || ymin >= image_height
    )) {
        aabb_xyxy = {0, 0, 0, 0};
        return;
    }

    aabb_xyxy = {(int)floor(xmin), (int)floor(ymin), (int)ceil(xmax), (int)ceil(ymax)};
    depth = mean_c.z;

    normal = normalize(mul(R, cross(vert1-vert0, vert2-vert0)));
    normal *= -sign(dot(normal, mean_c));

    out_hardness = hardness;
}


#define _DEF_projection_opaque_triangle(funname, ...) \
    [CudaDeviceExport] \
    [ForceInline] \
    [Differentiable] \
    void funname(_projection_opaque_triangle_ARGS) { \
        return projection_opaque_triangle<__VA_ARGS__>( \
            mean, quat, scale, hardness, R, t, \
            /* vert0, vert1, vert2, hardness, R, t, */ \
            fx, fy, cx, cy, radial_coeffs, tangential_coeffs, thin_prism_coeffs, \
            image_width, image_height, near_plane, far_plane, \
            aabb_xyxy, depth, normal, uv0, uv1, uv2, out_hardness \
        ); \
    }

_DEF_projection_opaque_triangle(projection_opaque_triangle_persp, CameraModelType::PINHOLE, false)
_DEF_projection_opaque_triangle(projection_opaque_triangle_fisheye, CameraModelType::FISHEYE, false)

    // const float3 vert0, \
    // const float3 vert1, \
    // const float3 vert2, \

    // out float3 v_vert0, \
    // out float3 v_vert1, \
    // out float3 v_vert2, \

#define _projection_opaque_triangle_vjp_ARGS \
    const float3 mean, \
    const float4 quat, \
    const float3 scale, \
    const float2 hardness, \
    const float3x3 R, \
    const float3 t, \
    const float fx, \
    const float fy, \
    const float cx, \
    const float cy, \
    const float4 radial_coeffs, \
    const float2 tangential_coeffs, \
    const float2 thin_prism_coeffs, \
    const uint image_width, \
    const uint image_height, \
    const float v_depth, \
    const float3 v_normal, \
    const float2 v_uv0, \
    const float2 v_uv1, \
    const float2 v_uv2, \
    const float2 v_out_hardness, \
    out float3 v_mean, \
    out float4 v_quat, \
    out float3 v_scale, \
    out float2 v_hardness, \
    out float3x3 v_R, \
    out float3 v_t

_DEF_projection_opaque_triangle(_projection_opaque_triangle_persp_differentiable, CameraModelType::PINHOLE, true)
_DEF_projection_opaque_triangle(_projection_opaque_triangle_fisheye_differentiable, CameraModelType::FISHEYE, true)


[CudaDeviceExport]
void projection_opaque_triangle_persp_vjp(
    _projection_opaque_triangle_vjp_ARGS
) {
    DifferentialPair<float3> dp_mean = diffPair(mean);
    DifferentialPair<float4> dp_quat = diffPair(quat);
    DifferentialPair<float3> dp_scale = diffPair(scale);
    // DifferentialPair<float3> dp_vert0 = diffPair(vert0);
    // DifferentialPair<float3> dp_vert1 = diffPair(vert1);
    // DifferentialPair<float3> dp_vert2 = diffPair(vert2);
    DifferentialPair<float2> dp_hardness = diffPair(hardness);
    DifferentialPair<float3x3> dp_R = diffPair(R);
    DifferentialPair<float3> dp_t = diffPair(t);
    bwd_diff(_projection_opaque_triangle_persp_differentiable)(
        dp_mean, dp_quat, dp_scale, dp_hardness, dp_R, dp_t,
        // dp_vert0, dp_vert1, dp_vert2, dp_hardness, dp_R, dp_t,
        fx, fy, cx, cy, radial_coeffs, tangential_coeffs, thin_prism_coeffs,
        image_width, image_height, 0.0f, 1e30f,
        v_depth, v_normal, v_uv0, v_uv1, v_uv2, v_out_hardness
    );
    v_mean = dp_mean.d;
    v_quat = dp_quat.d;
    v_scale = dp_scale.d;
    // v_vert0 = dp_vert0.d;
    // v_vert1 = dp_vert1.d;
    // v_vert2 = dp_vert2.d;
    v_hardness = dp_hardness.d;
    v_R = dp_R.d;
    v_t = dp_t.d;
}

[CudaDeviceExport]
void projection_opaque_triangle_fisheye_vjp(
    _projection_opaque_triangle_vjp_ARGS
) {
    DifferentialPair<float3> dp_mean = diffPair(mean);
    DifferentialPair<float4> dp_quat = diffPair(quat);
    DifferentialPair<float3> dp_scale = diffPair(scale);
    // DifferentialPair<float3> dp_vert0 = diffPair(vert0);
    // DifferentialPair<float3> dp_vert1 = diffPair(vert1);
    // DifferentialPair<float3> dp_vert2 = diffPair(vert2);
    DifferentialPair<float2> dp_hardness = diffPair(hardness);
    DifferentialPair<float3x3> dp_R = diffPair(R);
    DifferentialPair<float3> dp_t = diffPair(t);
    bwd_diff(_projection_opaque_triangle_fisheye_differentiable)(
        dp_mean, dp_quat, dp_scale, dp_hardness, dp_R, dp_t,
        // dp_vert0, dp_vert1, dp_vert2, dp_hardness, dp_R, dp_t,
        fx, fy, cx, cy, radial_coeffs, tangential_coeffs, thin_prism_coeffs,
        image_width, image_height, 0.0f, 1e30f,
        v_depth, v_normal, v_uv0, v_uv1, v_uv2, v_out_hardness
    );
    v_mean = dp_mean.d;
    v_quat = dp_quat.d;
    v_scale = dp_scale.d;
    // v_vert0 = dp_vert0.d;
    // v_vert1 = dp_vert1.d;
    // v_vert2 = dp_vert2.d;
    v_hardness = dp_hardness.d;
    v_R = dp_R.d;
    v_t = dp_t.d;
}


[ForceInline]
[Differentiable]
float antialiase_map(float d, float dmax, float hardness) {
    hardness = clamp(hardness, 0.005f, 0.995f);

    // return d/dmax >= 0.0f ? 0.0f : pow(-d/dmax, 1.0f-hardness);

    // float a = d/(2.0f*dmax*(1.0f-hardness)+hardness);
    // // return clamp(-a+0.5f, 0.0f, 0.999f);
    // return 0.999f/(1.0f+exp(4.0f*a));

    // float a = -d / dmax + exp2(-1.0f/(1.0f-hardness));
    // float a = -d / dmax;
    float a = 1.0f - (1.0f + d / dmax) * (1.0f - exp2(-1.0f/(1.0f-hardness)));
    return a <= 0.0f ? 0.0f : min(pow(a, 1.0f-hardness), 0.999f);

    // float a = -d / dmax * (-log(1.0f-hardness));
    // return 0.999f/(1.0f + exp(-4.0f*a));
}

[CudaDeviceExport]
[Differentiable]
float evaluate_alpha_opaque_triangle_fast(
    float2 v0, float2 v1, float2 v2, float2 hardness,
    no_diff float2 p
) {
    float2 e0 = v1-v0, e1 = v2-v1, e2 = v0-v2;
    float se = sign(det(e0, e1));
    float de0 = se * det(p-v0, normalize(e0));
    float de1 = se * det(p-v1, normalize(e1));
    float de2 = se * det(p-v2, normalize(e2));
    float dv = max(max(de0, de1), de2);
    float dmax = abs(det(e0, e1)) / (length(e0) + length(e1) + length(e2));
    return hardness.x * antialiase_map(dv, dmax, hardness.y);
}

[CudaDeviceExport]
void evaluate_alpha_opaque_triangle_fast_vjp(
    float2 v0, float2 v1, float2 v2, float2 hardness,
    float2 p, float v_alpha,
    out float2 v_v0, out float2 v_v1, out float2 v_v2, out float2 v_hardness,
) {
    DifferentialPair<float2> dp_v0 = diffPair(v0);
    DifferentialPair<float2> dp_v1 = diffPair(v1);
    DifferentialPair<float2> dp_v2 = diffPair(v2);
    DifferentialPair<float2> dp_hardness = diffPair(hardness);
    bwd_diff(evaluate_alpha_opaque_triangle_fast)(dp_v0, dp_v1, dp_v2, dp_hardness, p, v_alpha);
    v_v0 = dp_v0.d;
    v_v1 = dp_v2.d;
    v_v2 = dp_v1.d;
    v_hardness = dp_hardness.d;
}

[CudaDeviceExport]
[Differentiable]
float evaluate_alpha_opaque_triangle_precise(
    float2 v0, float2 v1, float2 v2, float2 hardness,
    no_diff float2 p
) {
    float2 e0 = v1-v0, e1 = v2-v1, e2 = v0-v2;
    float dv0 = length((p-v0) - e0 * clamp(dot(p-v0,e0)/dot(e0,e0), 0.f, 1.f));
    float dv1 = length((p-v1) - e1 * clamp(dot(p-v1,e1)/dot(e1,e1), 0.f, 1.f));
    float dv2 = length((p-v2) - e2 * clamp(dot(p-v2,e2)/dot(e2,e2), 0.f, 1.f));
    float se = sign(det(e0, e1));
    float s0 = se * det(p-v0, e0);
    float s1 = se * det(p-v1, e1);
    float s2 = se * det(p-v2, e2);
    float sv = sign(max(max(s0, s1), s2));
    float dv = sv * min(min(dv0, dv1), dv2);
    float dmax = abs(det(e0, e1)) / (length(e0) + length(e1) + length(e2));
    return hardness.x * antialiase_map(dv, dmax, hardness.y);
}

[CudaDeviceExport]
void evaluate_alpha_opaque_triangle_precise_vjp(
    float2 v0, float2 v1, float2 v2, float2 hardness,
    float2 p, float v_alpha,
    out float2 v_v0, out float2 v_v1, out float2 v_v2, out float2 v_hardness,
) {
    DifferentialPair<float2> dp_v0 = diffPair(v0);
    DifferentialPair<float2> dp_v1 = diffPair(v1);
    DifferentialPair<float2> dp_v2 = diffPair(v2);
    DifferentialPair<float2> dp_hardness = diffPair(hardness);
    bwd_diff(evaluate_alpha_opaque_triangle_precise)(dp_v0, dp_v1, dp_v2, dp_hardness, p, v_alpha);
    v_v0 = dp_v0.d;
    v_v1 = dp_v2.d;
    v_v2 = dp_v1.d;
    v_hardness = dp_hardness.d;
}
