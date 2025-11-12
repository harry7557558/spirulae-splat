#include "projection_utils.slang"

[ForceInline]
float hash41(float4 seed) {
    // https://www.shadertoy.com/view/4djSRW
    seed = fract(seed * float4(.1031, .1030, .0973, .1099));
    seed += dot(seed, seed.wzxy + 33.33);
    return fract((seed.x + seed.y) * (seed.z + seed.w));
}

#define _projection_opaque_triangle_eval3d_ARGS \
    const float3 mean, \
    const float4 quat, \
    const float3 scale, \
    const float2 hardness, \
    const float3 sh_coeffs[16], \
    const float3 ch_coeffs[2], \
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
    out float depth, \
    out float3 verts[3], \
    out float3 rgbs[3], \
    out float3 normal

[ForceInline]
[Differentiable]
void projection_opaque_triangle_eval3d<let camera_model: CameraModelType, let assume_valid: bool>(
    _projection_opaque_triangle_eval3d_ARGS
) {
    // float3 mean = (vert0+vert1+vert2)/3.0f;
    float3 mean_c;
    posW2C(R, t, mean, mean_c);
    float cull_zc = (camera_model == CameraModelType::FISHEYE) ?
        length(mean_c) : mean_c.z;
    if (!assume_valid && (cull_zc < near_plane || cull_zc > far_plane)) {
        aabb_xyxy = {0, 0, 0, 0};
        return;
    }

    // transform mean+scale+quat to vertices
    float sx = exp(scale.x),
        sy = exp(scale.y),
        sz = scale.z - 0.5*(scale.x+scale.y);
    float3x3 M = normalized_quat_to_rotmat(quat);
    float3 vert0 = mul(M, float3(sx, 0, 0)) + mean;
    float3 vert1 = mul(M, float3(sx*(-0.5+sz), sy, 0)) + mean;
    float3 vert2 = mul(M, float3(sx*(-0.5-sz), -sy, 0)) + mean;

    // transform center to camera space
    float3 vert0_c, vert1_c, vert2_c;
    posW2C(R, t, vert0, vert0_c);
    posW2C(R, t, vert1, vert1_c);
    posW2C(R, t, vert2, vert2_c);
    float3 cull_z = (camera_model == CameraModelType::FISHEYE) ?
        float3(length(vert0_c), length(vert1_c), length(vert2_c)) :
        float3(vert0_c.z, vert1_c.z, vert2_c.z);
    if (!assume_valid && (
        cull_z.x < near_plane || cull_z.x > far_plane ||
        cull_z.y < near_plane || cull_z.y > far_plane ||
        cull_z.z < near_plane || cull_z.z > far_plane
    )) {
        aabb_xyxy = {0, 0, 0, 0};
        return;
    }

    float2 uv0, uv1, uv2;
    bool all_valid = true;
    if (camera_model == CameraModelType::PINHOLE) {
        float4 intrins = float4(fx, fy, cx, cy);
        // TODO: some approximation better than this
        all_valid &= persp_proj<assume_valid>(vert0_c, intrins, dist_coeffs, uv0);
        all_valid &= persp_proj<assume_valid>(vert1_c, intrins, dist_coeffs, uv1);
        all_valid &= persp_proj<assume_valid>(vert2_c, intrins, dist_coeffs, uv2);
    }
    else if (camera_model == CameraModelType::FISHEYE) {
        float4 intrins = float4(fx, fy, cx, cy);
        // TODO: some approximation better than this
        all_valid &= fisheye_proj<assume_valid>(vert0_c, intrins, dist_coeffs, uv0);
        all_valid &= fisheye_proj<assume_valid>(vert1_c, intrins, dist_coeffs, uv1);
        all_valid &= fisheye_proj<assume_valid>(vert2_c, intrins, dist_coeffs, uv2);
    }
    else if (camera_model == CameraModelType::ORTHO) {} // TODO
    if (!assume_valid && !all_valid) {
        aabb_xyxy = { 0, 0, 0, 0 };
        return;
    }

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
        ymax <= 0 || ymin >= image_height ||
        (mean_c.z <= 0.0 && ((xmin <= 0 && xmax >= image_width) || (ymin <= 0 && ymax >= image_width)))
    )) {
        aabb_xyxy = {0, 0, 0, 0};
        return;
    }

    aabb_xyxy = {(int)floor(xmin), (int)floor(ymin), (int)ceil(xmax), (int)ceil(ymax)};
    // depth = (vert0_c.z + vert1_c.z + vert2_c.z) / 3.0;
    depth = length(vert0_c + vert1_c + vert2_c) / 3.0;
    depth = log(depth + 1e-12);

    float3 campos = -mul(transpose(R), t);
    float3 color;
    sh_coeffs_to_color(mean - campos, sh_coeffs, color);
    #define colormap(c) max((c)+0.5f, 0.0f)
    rgbs[0] = colormap(color + ch_coeffs[0]);
    rgbs[1] = colormap(color - ch_coeffs[0] * 0.5f + ch_coeffs[1] * sqrt(0.75f));
    rgbs[2] = colormap(color - ch_coeffs[0] * 0.5f - ch_coeffs[1] * sqrt(0.75f));
    #undef colormap

    verts[0] = vert0;
    verts[1] = vert1;
    verts[2] = vert2;

    normal = normalize(cross(vert1_c - vert0_c, vert2_c - vert0_c));
    normal *= -sign(dot(normal, mean_c));
}


#define _DEF_projection_opaque_triangle_eval3d(funname, ...) \
    [CudaDeviceExport] \
    [ForceInline] \
    [Differentiable] \
    void funname(_projection_opaque_triangle_eval3d_ARGS) { \
        return projection_opaque_triangle_eval3d<__VA_ARGS__>( \
            mean, quat, scale, hardness, sh_coeffs, ch_coeffs, R, t, \
            fx, fy, cx, cy, dist_coeffs, \
            image_width, image_height, near_plane, far_plane, \
            aabb_xyxy, depth, verts, rgbs, normal \
        ); \
    }

_DEF_projection_opaque_triangle_eval3d(projection_opaque_triangle_eval3d_persp, CameraModelType::PINHOLE, false)
_DEF_projection_opaque_triangle_eval3d(projection_opaque_triangle_eval3d_fisheye, CameraModelType::FISHEYE, false)

#define _projection_opaque_triangle_eval3d_vjp_ARGS \
    const float3 mean, \
    const float4 quat, \
    const float3 scale, \
    const float2 hardness, \
    const float3 sh_coeffs[16], \
    const float3 ch_coeffs[2], \
    const float3x3 R, \
    const float3 t, \
    const float fx, \
    const float fy, \
    const float cx, \
    const float cy, \
    const float[10] dist_coeffs, \
    const uint image_width, \
    const uint image_height, \
    const float v_depth, \
    const float3 v_verts[3], \
    const float3 v_rgbs[3], \
    const float3 v_normal, \
    out float3 v_mean, \
    out float4 v_quat, \
    out float3 v_scale, \
    out float2 v_hardness, \
    out float3 v_sh_coeffs[16], \
    out float3 v_ch_coeffs[2], \
    out float3x3 v_R, \
    out float3 v_t

_DEF_projection_opaque_triangle_eval3d(_projection_opaque_triangle_eval3d_persp_differentiable, CameraModelType::PINHOLE, true)
_DEF_projection_opaque_triangle_eval3d(_projection_opaque_triangle_eval3d_fisheye_differentiable, CameraModelType::FISHEYE, true)


[CudaDeviceExport]
void projection_opaque_triangle_eval3d_persp_vjp(
    _projection_opaque_triangle_eval3d_vjp_ARGS
) {
    DifferentialPair<float3> dp_mean = diffPair(mean);
    DifferentialPair<float4> dp_quat = diffPair(quat);
    DifferentialPair<float3> dp_scale = diffPair(scale);
    DifferentialPair<float2> dp_hardness = diffPair(hardness);
    DifferentialPair<float3[16]> dp_sh_coeffs = diffPair(sh_coeffs);
    DifferentialPair<float3[2]> dp_ch_coeffs = diffPair(ch_coeffs);
    DifferentialPair<float3x3> dp_R = diffPair(R);
    DifferentialPair<float3> dp_t = diffPair(t);
    bwd_diff(_projection_opaque_triangle_eval3d_persp_differentiable)(
        dp_mean, dp_quat, dp_scale, dp_hardness, dp_sh_coeffs, dp_ch_coeffs, dp_R, dp_t,
        fx, fy, cx, cy, dist_coeffs,
        image_width, image_height, 0.0f, 1e30f,
        v_depth, v_verts, v_rgbs, v_normal
    );
    v_mean = dp_mean.d;
    v_quat = dp_quat.d;
    v_scale = dp_scale.d;
    v_hardness = dp_hardness.d;
    v_sh_coeffs = dp_sh_coeffs.d;
    v_ch_coeffs = dp_ch_coeffs.d;
    v_R = dp_R.d;
    v_t = dp_t.d;
}

[CudaDeviceExport]
void projection_opaque_triangle_eval3d_fisheye_vjp(
    _projection_opaque_triangle_eval3d_vjp_ARGS
) {
    DifferentialPair<float3> dp_mean = diffPair(mean);
    DifferentialPair<float4> dp_quat = diffPair(quat);
    DifferentialPair<float3> dp_scale = diffPair(scale);
    DifferentialPair<float2> dp_hardness = diffPair(hardness);
    DifferentialPair<float3[16]> dp_sh_coeffs = diffPair(sh_coeffs);
    DifferentialPair<float3[2]> dp_ch_coeffs = diffPair(ch_coeffs);
    DifferentialPair<float3x3> dp_R = diffPair(R);
    DifferentialPair<float3> dp_t = diffPair(t);
    bwd_diff(_projection_opaque_triangle_eval3d_fisheye_differentiable)(
        dp_mean, dp_quat, dp_scale, dp_hardness, dp_sh_coeffs, dp_ch_coeffs, dp_R, dp_t,
        fx, fy, cx, cy, dist_coeffs,
        image_width, image_height, 0.0f, 1e30f,
        v_depth, v_verts, v_rgbs, v_normal
    );
    v_mean = dp_mean.d;
    v_quat = dp_quat.d;
    v_scale = dp_scale.d;
    v_hardness = dp_hardness.d;
    v_sh_coeffs = dp_sh_coeffs.d;
    v_ch_coeffs = dp_ch_coeffs.d;
    v_R = dp_R.d;
    v_t = dp_t.d;
}

[Differentiable]
[ForceInline]
bool ray_triangle_intersection_uvt(
    float3 ray_o, float3 ray_d,
    float3 verts[3],
    out float u, out float v, out float t
) {
    float3 v1v0 = verts[1] - verts[0];
    float3 v2v0 = verts[2] - verts[0];
    float3 rov0 = ray_o - verts[0];
    float3 n = cross(v1v0, v2v0);
    float3 q = cross(rov0, ray_d);
    float d = 1.0 / dot(ray_d, n);
    u = d * dot(-q, v2v0);
    v = d * dot(q, v1v0);
    t = d * dot(-n, rov0);
    return (u >= 0.0 && v >= 0.0 && (u + v) <= 1.0 && t >= 0.0);
}

[CudaDeviceExport]
[Differentiable]
float evaluate_alpha_opaque_triangle(
    float3 verts[3], float2 hardness,
    float3 ray_o, float3 ray_d
) {
    float u, v, t;
    if (!ray_triangle_intersection_uvt(ray_o, ray_d, verts, u, v, t))
        return 0.0f;
    float opac = min(min(u, v), sqrt(0.5f) * (1.0f - u - v)) * (2.0f + sqrt(2.0f));
    float h = clamp(hardness.y, 0.0f, 0.9999f);
    // float w = pow(opac, 1.0f - h);
    float w = 1.0 - pow(1.0-opac, 1.0/(1.0-h));
    float o = hardness.x;
    return opac < 0.0f ? 0.0f : min(o * w, 0.995f);
}

[CudaDeviceExport]
void evaluate_alpha_opaque_triangle_vjp(
    float3 verts[3], float2 hardness,
    float3 ray_o, float3 ray_d, float v_alpha,
    out float3 v_verts[3], out float2 v_hardness,
    out float3 v_ray_o, out float3 v_ray_d
) {
    DifferentialPair<float3[3]> dp_verts = diffPair(verts);
    DifferentialPair<float2> dp_hardness = diffPair(hardness);
    DifferentialPair<float3> dp_ray_o = diffPair(ray_o);
    DifferentialPair<float3> dp_ray_d = diffPair(ray_d);
    bwd_diff(evaluate_alpha_opaque_triangle)(
        dp_verts, dp_hardness,
        dp_ray_o, dp_ray_d, v_alpha
    );
    v_verts = dp_verts.d;
    v_hardness = dp_hardness.d;
    v_ray_o = dp_ray_o.d;
    v_ray_d = dp_ray_d.d;
}

[CudaDeviceExport]
float evaluate_sorting_depth_opaque_triangle(
    float3 verts[3], float3 rgbs[3],
    float3 ray_o, float3 ray_d
) {
    float u, v, t;
    ray_triangle_intersection_uvt(ray_o, ray_d, verts, u, v, t);
    return t;
}

[CudaDeviceExport]
[Differentiable]
void evaluate_color_opaque_triangle(
    float3 verts[3], float3 rgbs[3],
    float3 ray_o, float3 ray_d,
    out float3 color, out float depth
) {
    float u, v;
    ray_triangle_intersection_uvt(ray_o, ray_d, verts, u, v, depth);
    float w = 1.0f - u - v;

    color = rgbs[0] * w + rgbs[1] * u + rgbs[2] * v;
    // depth = length(verts[0] * w + verts[1] * u + verts[2] * v - ray_o);
    depth = log(max(depth, 1e-12));
}

[CudaDeviceExport]
void evaluate_color_opaque_triangle_vjp(
    float3 verts[3], float3 rgbs[3],
    float3 ray_o, float3 ray_d, float3 v_color, float v_depth,
    out float3 v_verts[3], out float3 v_rgbs[3],
    out float3 v_ray_o, out float3 v_ray_d
) {
    DifferentialPair<float3[3]> dp_verts = diffPair(verts);
    DifferentialPair<float3[3]> dp_rgbs = diffPair(rgbs);
    DifferentialPair<float3> dp_ray_o = diffPair(ray_o);
    DifferentialPair<float3> dp_ray_d = diffPair(ray_d);
    bwd_diff(evaluate_color_opaque_triangle)(
        dp_verts, dp_rgbs, dp_ray_o, dp_ray_d,
        v_color, v_depth
    );
    v_verts = dp_verts.d;
    v_rgbs = dp_rgbs.d;
    v_ray_o = dp_ray_o.d;
    v_ray_d = dp_ray_d.d;
}
