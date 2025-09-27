// All per-Gaussian losses in one kernel

[ForceInline]
[Differentiable]
float erank_reg_3dgs(float3 scales, float erank_reg, float erank_reg_s3) {
    scales = exp(2.0*scales);
    float x = scales.x, y = scales.y, z = scales.z;
    float s = x + y + z;
    float s1 = max(max(x, y), z) / s;
    float s3 = min(min(x, y), z) / s;
    float s2 = 1.0 - s1 - s3;
    float r = exp(-s1*log(s1) - s2*log(s2) - s3*log(s3));
    float reg = max(-log(r-0.99999), 0.0);
    return erank_reg * reg + erank_reg_s3 * s3;
}

[ForceInline]
[Differentiable]
float erank_reg_2dgs(float2 scales, float erank_reg) {
    scales = exp(2.0*scales);
    float x = scales.x, y = scales.y;
    float s = x + y;
    float s1 = max(x, y) / s;
    float s2 = min(x, y) / s;
    float r = exp(-s1*log(s1) - s2*log(s2));
    float reg = max(-log(r-0.99999), 0.0);
    return erank_reg * reg;
}


[ForceInline]
[Differentiable]
float scale_reg_3dgs(float3 scales, float max_gauss_ratio, float scale_regularization_weight) {
    scales = exp(scales);
    float x = scales.x, y = scales.y, z = scales.z;
    float scale_reg = max(max(max(x, y), z) / min(min(x, y), z), max_gauss_ratio) - max_gauss_ratio;
    return scale_regularization_weight * scale_reg;
}

[ForceInline]
[Differentiable]
float scale_reg_2dgs(float2 scales, float max_gauss_ratio, float scale_regularization_weight) {
    scales = exp(scales);
    float x = scales.x, y = scales.y;
    float scale_reg = max(max(x, y) / min(x, y), max_gauss_ratio) - max_gauss_ratio;
    return scale_regularization_weight * scale_reg;
}


[ForceInline]
[Differentiable]
float mcmc_opacity_reg(float opacity, float mcmc_opacity_reg_weight) {
    float reg = 1.0 / (1.0 + exp(-opacity));
    return mcmc_opacity_reg_weight * reg;
}

[ForceInline]
[Differentiable]
float mcmc_scale_reg_3dgs(float3 scales, float mcmc_scale_reg_weight) {
    scales = exp(scales);
    return mcmc_scale_reg_weight * (scales.x + scales.y + scales.z) / 3.0;
}

[ForceInline]
[Differentiable]
float mcmc_scale_reg_2dgs(float2 scales, float mcmc_scale_reg_weight) {
    scales = exp(scales);
    return mcmc_scale_reg_weight * (scales.x + scales.y) / 2.0;
}

[ForceInline]
[Differentiable]
float quat_norm_reg(float4 quat, float quat_norm_reg_weight) {
    float quat_norm = length(quat);
    float quat_norm_reg = quat_norm - 1.0 - log(quat_norm);
    return quat_norm_reg_weight * quat_norm_reg;
}


[Differentiable]
[CudaDeviceExport]
float[5] per_splat_losses(
    bool is_3dgs,
    float3 scales, float opacity, float4 quat,
    no_diff float mcmc_opacity_reg_weight,
    no_diff float mcmc_scale_reg_weight,
    no_diff float max_gauss_ratio,
    no_diff float scale_regularization_weight,
    no_diff float erank_reg_weight,
    no_diff float erank_reg_weight_s3,
    no_diff float quat_norm_reg_weight
) {
    float losses[5];
    losses[0] = mcmc_opacity_reg(opacity, mcmc_opacity_reg_weight);
    losses[4] = quat_norm_reg(quat, quat_norm_reg_weight);
    if (is_3dgs) {
        losses[1] = mcmc_scale_reg_3dgs(scales, mcmc_scale_reg_weight);
        losses[2] = scale_reg_3dgs(scales, max_gauss_ratio, scale_regularization_weight);
        losses[3] = erank_reg_3dgs(scales, erank_reg_weight, erank_reg_weight_s3);
    } else {
        losses[1] = mcmc_scale_reg_2dgs(scales.xy, mcmc_scale_reg_weight);
        losses[2] = scale_reg_2dgs(scales.xy, max_gauss_ratio, scale_regularization_weight);
        losses[3] = erank_reg_2dgs(scales.xy, erank_reg_weight);
    }
    return losses;
}


[CudaDeviceExport]
void per_splat_losses_bwd(
    bool is_3dgs,
    float3 scales, float opacity, float4 quat,
    float[5] v_loss,
    out float3 v_scales, out float v_opacity, out float4 v_quat,
    float mcmc_opacity_reg_weight,
    float mcmc_scale_reg_weight,
    float max_gauss_ratio,
    float scale_regularization_weight,
    float erank_reg_weight,
    float erank_reg_weight_s3,
    float quat_norm_reg_weight
) {
    var p_scales = diffPair(scales);
    var p_opacity = diffPair(opacity);
    var p_quat = diffPair(quat);
    bwd_diff(per_splat_losses)(
        is_3dgs,
        p_scales, p_opacity, p_quat,
        mcmc_opacity_reg_weight,
        mcmc_scale_reg_weight,
        max_gauss_ratio,
        scale_regularization_weight,
        erank_reg_weight,
        erank_reg_weight_s3,
        quat_norm_reg_weight,
        v_loss
    );
    v_scales = p_scales.d;
    v_opacity = p_opacity.d;
    v_quat = p_quat.d;
}
