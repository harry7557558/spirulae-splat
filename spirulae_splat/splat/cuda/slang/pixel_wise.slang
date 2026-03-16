// suppress CudaDeviceExport to reduce code size
#ifdef CudaDeviceExport
#include "projection_utils.slang"
#else
#define CudaDeviceExport Differentiable
#include "projection_utils.slang"
#undef CudaDeviceExport
#endif

[Differentiable]
[CudaDeviceExport]
float rendered_depth_to_expected_depth(float depth, float transmittance) {
    return depth / max(1.0f - transmittance, 1e-10f);
}

[CudaDeviceExport]
void rendered_depth_to_expected_depth_bwd(
    float depth, float transmittance,
    float v_out_depth,
    out float v_depth, out float v_transmittance
) {
    var p_depth = diffPair(depth);
    var p_transmittance = diffPair(transmittance);
    bwd_diff(rendered_depth_to_expected_depth)(
        p_depth, p_transmittance,
        v_out_depth
    );
    v_depth = p_depth.d;
    v_transmittance = p_transmittance.d;
}

[Differentiable]
[CudaDeviceExport]
float3 blend_background(float3 rgb, float transmittance, float3 background) {
    return clamp(rgb + transmittance * background, 0.0, 1.0);
}

[CudaDeviceExport]
void blend_background_bwd(
    float3 rgb, float transmittance, float3 background,
    float3 v_out_rgb,
    out float3 v_rgb, out float v_transmittance, out float3 v_background
) {
    var p_rgb = diffPair(rgb);
    var p_transmittance = diffPair(transmittance);
    var p_background = diffPair(background);
    bwd_diff(blend_background)(
        p_rgb, p_transmittance, p_background,
        v_out_rgb
    );
    v_rgb = p_rgb.d;
    v_transmittance = p_transmittance.d;
    v_background = p_background.d;
}

[ForceInline]
[Differentiable]
[CudaDeviceExport]
float linear_rgb_to_srgb(float x) {
    const float cutoff = 0.0031308f;
    return (x < cutoff ? x * 12.92f : 1.055f * pow(x, 1.0f / 2.4f) - 0.055f);
}

[CudaDeviceExport]
float linear_rgb_to_srgb_grad(float x) {
    var dp_x = diffPair(x, 1.0);
    var dp_y = fwd_diff(linear_rgb_to_srgb)(dp_x);
    return dp_y.d;
}

[Differentiable]
[CudaDeviceExport]
float3 linear_rgb_to_srgb(float3 rgb) {
    return float3(
        linear_rgb_to_srgb(rgb.x),
        linear_rgb_to_srgb(rgb.y),
        linear_rgb_to_srgb(rgb.z)
    );
}

[CudaDeviceExport]
float3 linear_rgb_to_srgb_bwd(
    float3 rgb, float3 v_out_rgb
) {
    var p_rgb = diffPair(rgb);
    bwd_diff(linear_rgb_to_srgb)(
        p_rgb,
        v_out_rgb
    );
    return p_rgb.d;
}



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
