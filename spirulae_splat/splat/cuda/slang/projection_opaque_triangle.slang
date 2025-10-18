#include "projection_utils.slang"


#define _projection_opaque_triangle_ARGS \
    const float3 vert0, \
    const float3 vert1, \
    const float3 vert2, \
    const float hardness, \
    const float3x3 R, \
    const float3 t, \
    no_diff const float fx, \
    no_diff const float fy, \
    no_diff const float cx, \
    no_diff const float cy, \
    const uint image_width, \
    const uint image_height, \
    no_diff const float near_plane, \
    no_diff const float far_plane, \
    no_diff const float radius_clip, \
    out int2 radii, \
    out float depth, \
    out float2 uv0, \
    out float2 uv1, \
    out float2 uv2, \
    out float out_hardness

[ForceInline]
[Differentiable]
void projection_opaque_triangle<let camera_model: CameraModelType, let assume_valid: bool>(
    _projection_opaque_triangle_ARGS
) {

    // transform Gaussian center to camera space
    float3 vert0_c, vert1_c, vert2_c;
    posW2C(R, t, vert0, vert0_c);
    posW2C(R, t, vert1, vert1_c);
    posW2C(R, t, vert2, vert2_c);
    float3 vert_mean = (vert0+vert1+vert2)/3.0;
    float3 vert_mean_c;
    posW2C(R, t, vert_mean, vert_mean_c);
    if (!assume_valid && (
        vert0_c.z < near_plane || vert0_c.z > far_plane ||
        vert1_c.z < near_plane || vert1_c.z > far_plane ||
        vert2_c.z < near_plane || vert2_c.z > far_plane ||
        vert_mean_c.z < near_plane || vert_mean_c.z > far_plane
    )) {
        radii = {0, 0};
        return;
    }

    if (camera_model == CameraModelType::PINHOLE) {
        uv0 = vert0_c.xy / vert0_c.z;
        uv1 = vert1_c.xy / vert1_c.z;
        uv2 = vert2_c.xy / vert2_c.z;
    }
    else if (camera_model == CameraModelType::FISHEYE)
        {}  // TODO
    else if (camera_model == CameraModelType::ORTHO)
        {}  // TODO

    uv0 = float2(fx, fy) * uv0 + float2(cx, cy);
    uv1 = float2(fx, fy) * uv1 + float2(cx, cy);
    uv2 = float2(fx, fy) * uv2 + float2(cx, cy);

    // TODO: extend

    // compute tight rectangular bounding box (non differentiable)
    static const float eps2d = 0.5f;
    float x_max = min(max(max(uv0.x, uv1.x), uv2.x) + eps2d, float(image_width));
    float x_min = max(min(min(uv0.x, uv1.x), uv2.x) - eps2d, 0.0f);
    float y_max = min(max(max(uv0.y, uv1.y), uv2.y) + eps2d, float(image_height));
    float y_min = max(min(min(uv0.y, uv1.y), uv2.y) - eps2d, 0.0f);
    float radius_x = ceil(x_max - x_min);
    float radius_y = ceil(y_max - y_min);

    if (!assume_valid && radius_x <= radius_clip && radius_y <= radius_clip) {
        radii = {0, 0};
        return;
    }

    // mask out splats outside the image region
    if (!assume_valid && (
        x_max <= 0 || x_min >= image_width ||
        y_max <= 0 || y_min >= image_height
    )) {
        radii = {0, 0};
        return;
    }

    radii = {(int)radius_x, (int)radius_y};
    depth = vert_mean_c.z;

    out_hardness = hardness;
}


#define _DEF_projection_opaque_triangle(funname, ...) \
    [CudaDeviceExport] \
    [ForceInline] \
    [Differentiable] \
    void funname(_projection_opaque_triangle_ARGS) { \
        return projection_opaque_triangle<__VA_ARGS__>( \
            vert0, vert1, vert2, hardness, R, t, fx, fy, cx, cy, \
            image_width, image_height, near_plane, far_plane, radius_clip, \
            radii, depth, uv0, uv1, uv2, out_hardness \
        ); \
    }

_DEF_projection_opaque_triangle(projection_opaque_triangle_persp, CameraModelType::PINHOLE, false)


#define _projection_opaque_triangle_vjp_ARGS \
    const float3 vert0, \
    const float3 vert1, \
    const float3 vert2, \
    const float hardness, \
    const float3x3 R, \
    const float3 t, \
    const float fx, \
    const float fy, \
    const float cx, \
    const float cy, \
    const uint image_width, \
    const uint image_height, \
    const float v_depth, \
    const float2 v_uv0, \
    const float2 v_uv1, \
    const float2 v_uv2, \
    const float v_out_hardness, \
    out float3 v_vert0, \
    out float3 v_vert1, \
    out float3 v_vert2, \
    out float v_hardness, \
    out float3x3 v_R, \
    out float3 v_t

_DEF_projection_opaque_triangle(_projection_opaque_triangle_persp_differentiable, CameraModelType::PINHOLE, true)


[CudaDeviceExport]
void projection_opaque_triangle_persp_vjp(
    _projection_opaque_triangle_vjp_ARGS
) {
    DifferentialPair<float3> dp_vert0 = diffPair(vert0);
    DifferentialPair<float3> dp_vert1 = diffPair(vert0);
    DifferentialPair<float3> dp_vert2 = diffPair(vert0);
    DifferentialPair<float> dp_hardness = diffPair(hardness);
    DifferentialPair<float3x3> dp_R = diffPair(R);
    DifferentialPair<float3> dp_t = diffPair(t);
    bwd_diff(_projection_opaque_triangle_persp_differentiable)(
        dp_vert0, dp_vert1, dp_vert2, dp_hardness, dp_R, dp_t, fx, fy, cx, cy,
        image_width, image_height, 0.0f, 1e30f, -1.0f,
        v_depth, v_uv0, v_uv1, v_uv2, v_out_hardness
    );
    v_vert0 = dp_vert0.d;
    v_vert1 = dp_vert1.d;
    v_vert2 = dp_vert2.d;
    v_hardness = dp_hardness.d;
    v_R = dp_R.d;
    v_t = dp_t.d;
}
