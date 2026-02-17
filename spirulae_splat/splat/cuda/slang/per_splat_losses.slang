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


[ForceInline]
[Differentiable]
[CudaDeviceExport]
float[5] per_splat_losses(
    float3 scales, float opacity, float4 quat,
    no_diff float mcmc_opacity_reg_weight,
    no_diff float mcmc_scale_reg_weight,
    no_diff float max_gauss_ratio,
    no_diff float scale_regularization_weight,
    no_diff float erank_reg_weight,
    no_diff float erank_reg_weight_s3,
    no_diff float quat_norm_reg_weight
) {
    constexpr bool is_3dgs = true;
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
[Differentiable]
void per_splat_losses_bwd(
    float3 scales, float opacity, float4 quat,
    no_diff float[5] v_loss,
    out float3 v_scales, out float v_opacity, out float4 v_quat,
    no_diff float mcmc_opacity_reg_weight,
    no_diff float mcmc_scale_reg_weight,
    no_diff float max_gauss_ratio,
    no_diff float scale_regularization_weight,
    no_diff float erank_reg_weight,
    no_diff float erank_reg_weight_s3,
    no_diff float quat_norm_reg_weight
) {
    var p_scales = diffPair(scales);
    var p_opacity = diffPair(opacity);
    var p_quat = diffPair(quat);
    bwd_diff(per_splat_losses)(
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

[CudaDeviceExport]
void per_splat_losses_bwd(
    float3 scales, float opacity, float4 quat,
    float[5] v_loss,
    out float3 v_scales, out float v_opacity, out float4 v_quat,
    out float3 vr_scales, out float vr_opacity, out float4 vr_quat,
    out float3 h_scales, out float h_opacity, out float4 h_quat,
    float mcmc_opacity_reg_weight,
    float mcmc_scale_reg_weight,
    float max_gauss_ratio,
    float scale_regularization_weight,
    float erank_reg_weight,
    float erank_reg_weight_s3,
    float quat_norm_reg_weight
) {
    // gradient
    {
        var p_scales = diffPair(scales);
        var p_opacity = diffPair(opacity);
        var p_quat = diffPair(quat);
        bwd_diff(per_splat_losses)(
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

    // gradient residual product
    {
        float[5] losses = per_splat_losses(
            scales, opacity, quat,
            mcmc_opacity_reg_weight,
            mcmc_scale_reg_weight,
            max_gauss_ratio,
            scale_regularization_weight,
            erank_reg_weight,
            erank_reg_weight_s3,
            quat_norm_reg_weight
        );
        float[5] residuals;
        [ForceUnroll]
        for (int i = 0; i < 5; i++)
            residuals[i] = sqrt(max(losses[i], 0.0));
        var p_scales = diffPair(scales);
        var p_opacity = diffPair(opacity);
        var p_quat = diffPair(quat);
        bwd_diff(per_splat_losses)(
            p_scales, p_opacity, p_quat,
            mcmc_opacity_reg_weight,
            mcmc_scale_reg_weight,
            max_gauss_ratio,
            scale_regularization_weight,
            erank_reg_weight,
            erank_reg_weight_s3,
            quat_norm_reg_weight,
            residuals
        );
        vr_scales = p_scales.d;
        vr_opacity = p_opacity.d;
        vr_quat = p_quat.d;
    }

    // hessian diagonal
    constexpr float[5] loss_weights = { 1, 1, 1, 1, 1 };
    [ForceUnroll]
    for (int i = 0; i < 3; ++i) {
        var dp_scales = diffPair(scales, float3(float(i == 0), float(i == 1), float(i == 2)));
        var dp_opacity = diffPair(opacity, 0.0);
        var dp_quat = diffPair(quat, float4(0));
        DifferentialPair<float3> dp_v_scales;
        DifferentialPair<float> dp_v_opacity;
        DifferentialPair<float4> dp_v_quat;
        fwd_diff(per_splat_losses_bwd)(
            dp_scales, dp_opacity, dp_quat,
            loss_weights,
            dp_v_scales, dp_v_opacity, dp_v_quat,
            mcmc_opacity_reg_weight,
            mcmc_scale_reg_weight,
            max_gauss_ratio,
            scale_regularization_weight,
            erank_reg_weight,
            erank_reg_weight_s3,
            quat_norm_reg_weight
        );
        h_scales[i] = dp_v_scales.d[i];
    }
    {
        var dp_scales = diffPair(scales, float3(0));
        var dp_opacity = diffPair(opacity, 1.0);
        var dp_quat = diffPair(quat, float4(0));
        DifferentialPair<float3> dp_v_scales;
        DifferentialPair<float> dp_v_opacity;
        DifferentialPair<float4> dp_v_quat;
        fwd_diff(per_splat_losses_bwd)(
            dp_scales, dp_opacity, dp_quat,
            loss_weights,
            dp_v_scales, dp_v_opacity, dp_v_quat,
            mcmc_opacity_reg_weight,
            mcmc_scale_reg_weight,
            max_gauss_ratio,
            scale_regularization_weight,
            erank_reg_weight,
            erank_reg_weight_s3,
            quat_norm_reg_weight
        );
        h_opacity = dp_v_opacity.d;
    }
    [ForceUnroll]
    for (int i = 0; i < 4; ++i) {
        var dp_scales = diffPair(scales, float3(0));
        var dp_opacity = diffPair(opacity, 0.0);
        var dp_quat = diffPair(quat, float4(float(i == 0), float(i == 1), float(i == 2), float(i == 3)));
        DifferentialPair<float3> dp_v_scales;
        DifferentialPair<float> dp_v_opacity;
        DifferentialPair<float4> dp_v_quat;
        fwd_diff(per_splat_losses_bwd)(
            dp_scales, dp_opacity, dp_quat,
            loss_weights,
            dp_v_scales, dp_v_opacity, dp_v_quat,
            mcmc_opacity_reg_weight,
            mcmc_scale_reg_weight,
            max_gauss_ratio,
            scale_regularization_weight,
            erank_reg_weight,
            erank_reg_weight_s3,
            quat_norm_reg_weight
        );
        h_quat[i] = dp_v_quat.d[i];
    }
}


// suppress CudaDeviceExport to reduce code size
#define CudaDeviceExport Differentiable
#include "primitive_3dgs.slang"
#undef CudaDeviceExport

[CudaDeviceExport]
void mcmc_add_noise_3dgs(
    float scaler, float min_opacity,
    inout float3 mean, float3 scale, float4 quat, float opac
) {
    float k = 0.5 / min_opacity;
    float x0 = 1.0 - min_opacity;
    float op_sigmoid = 1.0 / (1.0 + exp(-k * ((1.0 - opac) - x0)));

    float3x3 covar;
    quat_scale_to_covar(normalize(quat), exp(scale), covar);

    float3 noise = randn3(float4(dot(mean, mean), dot(mean, scale), dot(scale, scale), dot(quat, float4(opac))));
    noise = noise * scaler * op_sigmoid;
    mean += mul(covar, noise);
}

[CudaDeviceExport]
void mcmc_add_noise_triangle(
    float scaler, float min_opacity,
    inout float3 mean, float3 scale, float4 quat, float opac
) {
    float k = 0.5 / min_opacity;
    float x0 = 1.0 - min_opacity;
    float op_sigmoid = 1.0 / (1.0 + exp(-k * ((1.0 - opac) - x0)));

    float3 vert0, vert1, vert2;
    map_opaque_triangle(mean, normalize(quat), scale, vert0, vert1, vert2);

    float3 vertc = (vert0 + vert1 + vert2) / 3.0f;
    float3 d0 = vert0 - vertc, d1 = vert1 - vertc, d2 = vert2 - vertc;
    float3 dn = 0.5f * min(min(length(d0), length(d1)), length(d2)) * normalize(cross(d0, d1));
    float3x3 covar = 0.5f * (outerProduct(d0) + outerProduct(d1) + outerProduct(d2) + outerProduct(dn)) / 3.5f;

    float3 noise = randn3(float4(dot(mean, mean), dot(mean, scale), dot(scale, scale), dot(quat, float4(opac))));
    noise = noise * scaler * op_sigmoid;
    mean += mul(covar, noise);
}
