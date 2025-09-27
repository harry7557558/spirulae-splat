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
