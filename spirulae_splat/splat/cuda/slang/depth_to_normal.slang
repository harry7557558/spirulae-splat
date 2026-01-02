#include "projection_utils.slang"

[CudaDeviceExport]
float3 generate_ray_d2n(
    const no_diff float2 pix_pos,
    const no_diff float4 intrins,
    const no_diff float[10] dist_coeffs,
    const bool is_fisheye,
    const bool is_ray_depth
) {
    float2 uv = (pix_pos - intrins.zw) / intrins.xy;
    if (!undistort_point(uv, dist_coeffs, 12, uv)) {
        return float3(0);
    }

    float3 raydir;
    if (is_fisheye) {
        float theta = length(uv);
        raydir = float3(uv / max(theta, 1e-7f) * sin(theta), cos(theta));
        if (!is_ray_depth)
            raydir /= raydir.z;
    }
    else {
        raydir = float3(uv.x, uv.y, 1.0f);
        if (is_ray_depth)
            raydir = normalize(raydir);
    }

    return raydir;
}

[CudaDeviceExport]
[Differentiable]
float3 points_to_normal(
    const float3 points[4] // x-, x+, y-, y+
) {
    float3 dx = points[1] - points[0];
    float3 dy = points[3] - points[2];
    float3 normal = cross(dx, -dy);
    if (dot(normal, normal) != 0.0)
        normal /= length(normal);
    // normal *= -sign(dot(normal, points[0]+points[1]+points[2]+points[3]));
    return normal;
}

[CudaDeviceExport]
void points_to_normal_vjp(
    const float3 points[4], // x-, x+, y-, y+
    const float3 v_normal,
    out float3 v_points[4]
) {
    DifferentialPair<float3[4]> dp_points = diffPair(points);
    bwd_diff(points_to_normal)(
        dp_points, v_normal
    );
    v_points = dp_points.d;
}

[CudaDeviceExport]
[Differentiable]
float3 depth_to_normal(
    const no_diff float2 pix_center,
    const no_diff float4 intrins,  // fx, fy, cx, cy
    const no_diff float[10] dist_coeffs,
    const bool is_fisheye,
    const bool is_ray_depth,
    const float4 depths  // x-, x+, y-, y+
) {
    float3 points[4];

    [ForceUnroll]
    for (int i = 0; i < 4; i++) {
        float2 pix_pos = pix_center + float2(
            (2*i-1)*float(i < 2), (2*(i-2)-1)*float(i>=2)
        );

        float2 uv = (pix_pos - intrins.zw) / intrins.xy;
        if (!no_diff undistort_point(uv, dist_coeffs, 12, uv)) {
            return float3(0.0f);
        }

        float3 raydir;
        if (is_fisheye) {
            float theta = length(uv);
            raydir = float3(uv/max(theta,1e-7f)*sin(theta), cos(theta));
            if (!is_ray_depth)
                raydir /= raydir.z;
        }
        else {
            raydir = float3(uv.x, uv.y, 1.0f);
            if (is_ray_depth)
                raydir = normalize(raydir);
        }

        points[i] = depths[i] * raydir;
    }

    float3 dx = points[1] - points[0];
    float3 dy = points[3] - points[2];
    float3 normal = cross(dx, -dy);
    if (dot(normal, normal) != 0.0)
        normal /= length(normal);
    // normal *= -sign(dot(normal, points[0]+points[1]+points[2]+points[3]));
    return normal;
}

[CudaDeviceExport]
void depth_to_normal_vjp(
    const no_diff float2 pix_center,
    const no_diff float4 intrins,  // fx, fy, cx, cy
    const no_diff float[10] dist_coeffs,
    const bool is_fisheye,
    const bool is_ray_depth,
    const float4 depths,  // x-, x+, y-, y+
    const float3 v_normal,
    out float4 v_depths
) {
    DifferentialPair<float4> dp_depths = diffPair(depths);
    bwd_diff(depth_to_normal)(
        pix_center,
        intrins, dist_coeffs,
        is_fisheye, is_ray_depth,
        dp_depths, v_normal
    );
    v_depths = dp_depths.d;
}


[CudaDeviceExport]
[Differentiable]
float ray_depth_to_linear_depth_factor(
    const no_diff float2 pix_center,
    const no_diff float4 intrins,  // fx, fy, cx, cy
    const no_diff float[10] dist_coeffs,
    const bool is_fisheye
) {

    float2 uv = (pix_center - intrins.zw) / intrins.xy;

    if (!no_diff undistort_point(uv, dist_coeffs, 12, uv))
        return 0.0f;

    float3 raydir;
    if (is_fisheye) {
        float theta = length(uv);
        raydir = float3(uv/max(theta,1e-7f)*sin(theta), cos(theta));
        raydir /= raydir.z;
    }
    else {
        raydir = float3(uv.x, uv.y, 1.0f);
    }

    return sign(raydir.z) / length(raydir);
}
