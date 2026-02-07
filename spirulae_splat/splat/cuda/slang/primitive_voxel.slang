#include "projection_utils.slang"

#define _projection_voxel_eval3d_ARGS \
    const no_diff float3 pos, \
    const no_diff float size, \
    const float[8] densities, \
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
    out float depth, \
    out float3 rgbs

[ForceInline]
[Differentiable]
void projection_voxel_eval3d<let camera_model : CameraModelType, let assume_valid : bool>(
    _projection_voxel_eval3d_ARGS
) {
    float3 pos_c[8];
    float min_z = far_plane, max_z = near_plane;
    [ForceUnroll]
    for (int i = 0; i < 8; i++) {
        float3 pos_i = pos + size * float3(i & 1, (i >> 1) & 1, i >> 2);
        posW2C(R, t, pos_i, pos_c[i]);
        float cull_zc = (camera_model == CameraModelType::FISHEYE) ?
            length(pos_c[i]) : pos_c[i].z;
        min_z = min(min_z, cull_zc);
        max_z = max(max_z, cull_zc);
    }
    if (!assume_valid && (min_z < near_plane || max_z > far_plane)) {
        aabb_xyxy = {0, 0, 0, 0};
        return;
    }
    float3 mean_c;
    posW2C(R, t, pos + 0.5 * size, mean_c);

    float2 uv[8];
    bool all_valid = true;
    [ForceUnroll]
    for (int i = 0; i < 8; i++) {
        float4 intrins = float4(fx, fy, cx, cy);
        // TODO: some approximation better than this
        if (camera_model == CameraModelType::PINHOLE)
            all_valid &= persp_proj<assume_valid>(pos_c[i], intrins, dist_coeffs, uv[i]);
        else if (camera_model == CameraModelType::FISHEYE)
            all_valid &= fisheye_proj<assume_valid>(pos_c[i], intrins, dist_coeffs, uv[i]);
    }
    if (!assume_valid && !all_valid) {
        aabb_xyxy = { 0, 0, 0, 0 };
        return;
    }

    // early skip low opacity
    if (!assume_valid) {
        float maxDensity = densities[0];
        [ForceUnroll]
        for (int i = 1; i < 8; i++)
            maxDensity = max(maxDensity, densities[i]);
        float maxPossibleOpacity = 1.0 - exp(-maxDensity * size * sqrt(3.0));
        if (maxPossibleOpacity <= ALPHA_THRESHOLD) {
            aabb_xyxy = { 0, 0, 0, 0 };
            return;
        }
    }

    // compute tight rectangular bounding box (non differentiable)

    float xmax = uv[0].x, xmin = uv[0].x, ymax = uv[0].y, ymin = uv[0].y;
    [ForceUnroll]
    for (int i = 1; i < 8; i++) {
        xmax = max(xmax, uv[i].x);
        xmin = min(xmin, uv[i].x);
        ymax = max(ymax, uv[i].y);
        ymin = min(ymin, uv[i].y);
    }

    // mask out splats outside the image region
    if (!assume_valid && (
        xmax <= 0 || xmin >= image_width ||
        ymax <= 0 || ymin >= image_height ||
        (min_z <= 0.0 && ((xmin <= 0 && xmax >= image_width) || (ymin <= 0 && ymax >= image_width)))
    )) {
        aabb_xyxy = {0, 0, 0, 0};
        return;
    }

    aabb_xyxy = { (int)floor(xmin), (int)floor(ymin), (int)ceil(xmax), (int)ceil(ymax) };
    depth = get_sorting_depth(mean_c);  // TODO: morton code

    float3 campos = -mul(transpose(R), t);
    sh_coeffs_to_color(mean_c - campos, sh_coeffs, rgbs);
    rgbs = max(rgbs + 0.5f, 0.0f);
}


#define _DEF_projection_voxel_eval3d(funname, ...) \
    [CudaDeviceExport] \
    [ForceInline] \
    [Differentiable] \
    void funname(_projection_voxel_eval3d_ARGS) { \
        return projection_voxel_eval3d<__VA_ARGS__>( \
            pos, size, densities, sh_coeffs, R, t, \
            fx, fy, cx, cy, dist_coeffs, \
            image_width, image_height, near_plane, far_plane, \
            aabb_xyxy, depth, rgbs \
        ); \
    }

_DEF_projection_voxel_eval3d(projection_voxel_eval3d_persp, CameraModelType::PINHOLE, false)
_DEF_projection_voxel_eval3d(projection_voxel_eval3d_fisheye, CameraModelType::FISHEYE, false)

#define _projection_voxel_eval3d_vjp_ARGS \
    const no_diff float3 pos, \
    const no_diff float size, \
    const float[8] densities, \
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
    const float3 v_rgb, \
    out float[8] v_densities, \
    out float3 v_sh_coeffs[16], \
    out float3x3 v_R, \
    out float3 v_t

_DEF_projection_voxel_eval3d(_projection_voxel_eval3d_persp_differentiable, CameraModelType::PINHOLE, true)
_DEF_projection_voxel_eval3d(_projection_voxel_eval3d_fisheye_differentiable, CameraModelType::FISHEYE, true)


[CudaDeviceExport]
void projection_voxel_eval3d_persp_vjp(
    _projection_voxel_eval3d_vjp_ARGS
) {
    DifferentialPair<float[8]> dp_densities = diffPair(densities);
    DifferentialPair<float3[16]> dp_sh_coeffs = diffPair(sh_coeffs);
    DifferentialPair<float3x3> dp_R = diffPair(R);
    DifferentialPair<float3> dp_t = diffPair(t);
    bwd_diff(_projection_voxel_eval3d_persp_differentiable)(
        pos, size, dp_densities, dp_sh_coeffs, dp_R, dp_t,
        fx, fy, cx, cy, dist_coeffs,
        image_width, image_height, 0.0f, 1e30f,
        0.0f, v_rgb
    );
    v_densities = dp_densities.d;
    v_sh_coeffs = dp_sh_coeffs.d;
    v_R = dp_R.d;
    v_t = dp_t.d;
}

[CudaDeviceExport]
void projection_voxel_eval3d_fisheye_vjp(
    _projection_voxel_eval3d_vjp_ARGS
) {
    DifferentialPair<float[8]> dp_densities = diffPair(densities);
    DifferentialPair<float3[16]> dp_sh_coeffs = diffPair(sh_coeffs);
    DifferentialPair<float3x3> dp_R = diffPair(R);
    DifferentialPair<float3> dp_t = diffPair(t);
    bwd_diff(_projection_voxel_eval3d_fisheye_differentiable)(
        pos, size, dp_densities, dp_sh_coeffs, dp_R, dp_t,
        fx, fy, cx, cy, dist_coeffs,
        image_width, image_height, 0.0f, 1e30f,
        0.0f, v_rgb
    );
    v_densities = dp_densities.d;
    v_sh_coeffs = dp_sh_coeffs.d;
    v_R = dp_R.d;
    v_t = dp_t.d;
}

[CudaDeviceExport]
[Differentiable]
[ForceInline]
bool ray_aabb_intersection(
    float3 ray_o, float3 ray_d,
    float3 center, float radius,
    out float t0, out float t1
) {
    float3 m = 1.0 / ray_d;
    float3 n = m * (ray_o - center);
    float3 k = abs(m) * radius;
    float3 ta = -n - k, tb = -n + k;
    t0 = max(max(ta.x, ta.y), max(ta.z, 0.0));
    t1 = min(min(tb.x, tb.y), tb.z);
    return t0 < t1;
}

[ForceInline]
[Differentiable]
float explin(float x) {
    return x > 1.1f ? x :
        exp(0.909090909091f*x - 0.904689820196f);
}

[Differentiable]
float interp(const float[8] densities, float3 w) {
    float w0 = (1.0f - w.z) * (1.0f - w.y) * (1.0f - w.x);
    float w1 = (1.0f - w.z) * (1.0f - w.y) * w.x;
    float w2 = (1.0f - w.z) * w.y * (1.0f - w.x);
    float w3 = (1.0f - w.z) * w.y * w.x;
    float w4 = w.z * (1.0f - w.y) * (1.0f - w.x);
    float w5 = w.z * (1.0f - w.y) * w.x;
    float w6 = w.z * w.y * (1.0f - w.x);
    float w7 = w.z * w.y * w.x;
    return w0 * densities[0] + w1 * densities[1] + w2 * densities[2] + w3 * densities[3] +
           w4 * densities[4] + w5 * densities[5] + w6 * densities[6] + w7 * densities[7];
}

[CudaDeviceExport]
[Differentiable]
float evaluate_alpha_voxel(
    const no_diff float3 pos,
    const no_diff float size,
    const float[8] densities,
    float3 ray_o, float3 ray_d
) {
    float t0, t1;
    if (!ray_aabb_intersection(ray_o, ray_d, pos + 0.5 * size, 0.5 * size, t0, t1))
        return 0.0f;
    const int K = 8;
    float accum = 0.0f;
    for (int i = 0; i < K; i++) {
        float t = lerp(t0, t1, ((float)i + 0.5f) / (float)K);
        float a = interp(densities, (ray_o + ray_d * t - pos) / size);
        accum += explin(a);
    }
    return min(1.0 - exp(-(t1 - t0) / K * accum), 0.999f);
}

[CudaDeviceExport]
void evaluate_alpha_voxel_vjp(
    const no_diff float3 pos,
    const no_diff float size,
    const float[8] densities,
    float3 ray_o, float3 ray_d, float v_alpha,
    out float[8] v_densities,
    out float3 v_ray_o, out float3 v_ray_d
) {
    DifferentialPair<float[8]> dp_densities = diffPair(densities);
    DifferentialPair<float3> dp_ray_o = diffPair(ray_o);
    DifferentialPair<float3> dp_ray_d = diffPair(ray_d);
    bwd_diff(evaluate_alpha_voxel)(
        pos, size, dp_densities,
        dp_ray_o, dp_ray_d, v_alpha
    );
    v_densities = dp_densities.d;
    v_ray_o = dp_ray_o.d;
    v_ray_d = dp_ray_d.d;
}

[CudaDeviceExport]
[Differentiable]
void evaluate_color_voxel(
    const no_diff float3 pos,
    const no_diff float size,
    const float[8] densities,
    float3 rgb,
    float3 ray_o, float3 ray_d,
    out float3 out_rgb, out float depth
) {
    out_rgb = rgb;
    float t0, t1;
    ray_aabb_intersection(ray_o, ray_d, pos + 0.5 * size, 0.5 * size, t0, t1);

    const int K = 8;
    float accum = 0.0f, depth_accum = 0.0f;
    for (int i = 0; i < K; i++) {
        float t = lerp(t0, t1, ((float)i + 0.5f) / (float)K);
        float a = interp(densities, (ray_o + ray_d * t - pos) / size);
        a = explin(a);
        accum += a;
        depth_accum += t * a;
    }
    depth = max(depth_accum / accum, 0.0);
}

[CudaDeviceExport]
void evaluate_color_voxel_vjp(
    float3 pos, float size, float[8] densities, float3 rgb,
    float3 ray_o, float3 ray_d, float3 v_out_rgb, float v_depth,
    out float[8] v_densities, out float3 v_rgb,
    out float3 v_ray_o, out float3 v_ray_d
) {
    DifferentialPair<float[8]> dp_densities = diffPair(densities);
    DifferentialPair<float3> dp_rgb = diffPair(rgb);
    DifferentialPair<float3> dp_ray_o = diffPair(ray_o);
    DifferentialPair<float3> dp_ray_d = diffPair(ray_d);
    bwd_diff(evaluate_color_voxel)(
        pos, size, dp_densities, dp_rgb,
        dp_ray_o, dp_ray_d, v_out_rgb, v_depth
    );
    v_densities = dp_densities.d;
    v_rgb = dp_rgb.d;
    v_ray_o = dp_ray_o.d;
    v_ray_d = dp_ray_d.d;
}
