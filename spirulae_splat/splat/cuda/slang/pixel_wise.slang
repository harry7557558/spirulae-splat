[Differentiable]
[CudaDeviceExport]
float3 blend_background(float3 rgb, float alpha, float3 background) {
    return clamp(rgb + (1.0 - alpha) * background, 0.0, 1.0);
}

[CudaDeviceExport]
void blend_background_bwd(
    float3 rgb, float alpha, float3 background,
    float3 v_out_rgb,
    out float3 v_rgb, out float v_alpha, out float3 v_background
) {
    var p_rgb = diffPair(rgb);
    var p_alpha = diffPair(alpha);
    var p_background = diffPair(background);
    bwd_diff(blend_background)(
        p_rgb, p_alpha, p_background,
        v_out_rgb
    );
    v_rgb = p_rgb.d;
    v_alpha = p_alpha.d;
    v_background = p_background.d;
}


[Differentiable]
[CudaDeviceExport]
float3 linear_rgb_to_srgb(float3 rgb) {
    const float cutoff = 0.0031308f;
    rgb.x = (rgb.x < cutoff ? rgb.x * 12.92f : 1.055f * pow(rgb.x, 1.0f / 2.4f) - 0.055f);
    rgb.y = (rgb.y < cutoff ? rgb.y * 12.92f : 1.055f * pow(rgb.y, 1.0f / 2.4f) - 0.055f);
    rgb.z = (rgb.z < cutoff ? rgb.z * 12.92f : 1.055f * pow(rgb.z, 1.0f / 2.4f) - 0.055f);
    return rgb;
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
