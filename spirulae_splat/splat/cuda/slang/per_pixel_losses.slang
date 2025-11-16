// All per-Pixel losses in one kernel

[ForceInline]
[Differentiable]
float mean3(float3 a) {
    return (a.x + a.y + a.z) * (1.0 / 3.0);
}

[ForceInline]
[Differentiable]
float l1_loss(float3 a, float3 b) {
    return mean3(abs(b - a));
}

[ForceInline]
[Differentiable]
float l2_loss(float3 a, float3 b) {
    return dot(b - a, b - a) * (1.0 / 3.0);
}

[ForceInline]
[Differentiable]
float3 normalize_normal(float3 normal, inout bool mask) {
    float norm2 = dot(normal, normal);
    if (norm2 == 0.0) {
        mask = false;
        return float3(0.0);
    }
    return normal * rsqrt(norm2);
}

[ForceInline]
[Differentiable]
float normal_loss(float3 normal, float3 normal_ref) {
    // return 1.0 - dot(normal, normal_ref);
    // return l1_loss(normal, normal_ref);
    float cos_sim_loss = 0.5 - 0.5 * dot(normal, normal_ref);
    return cos_sim_loss + sqrt(max(cos_sim_loss, 1e-12));
}

[ForceInline]
[Differentiable]
float bce_loss(float val, float ref) {
    const float eps = 1e-6;
    return -lerp(log(max(1.0 - val, eps)), log(max(val, eps)), ref);
}

[ForceInline]
[Differentiable]
float alpha_loss(float x, float y) {
    return bce_loss(max(x, y), y);
}

enum RawLossIndex {
    RgbL1,
    RgbL2,
    DepthSupX,
    DepthSupY,
    DepthSupXX,
    DepthSupYY,
    DepthSupXY,
    RenderNormalSup,
    DepthNormalSup,
    AlphaSup,
    AlphaSupUnder,
    NormalReg,
    AlphaReg,
    RgbDistReg,
    DepthDistReg,
    NormalDistReg,
    PixelsTotal,
    MaskTotal,
    DepthMaskTotal,
    RenderNormalMaskTotal,
    DepthNormalMaskTotal,
    NormalRegMaskTotal,
    length
};


[Differentiable]
[CudaDeviceExport]
float[RawLossIndex::length] per_pixel_losses(
    float3 render_rgb,
    float3 ref_rgb,
    float render_depth,
    float ref_depth,
    float3 render_normal,
    float3 depth_normal,
    float3 ref_normal,
    float render_alpha,
    float3 rgb_dist,
    float depth_dist,
    float3 normal_dist,
    bool ref_alpha,
    bool mask,
    bool depth_mask,
    bool normal_mask
) {
    float[RawLossIndex::length] losses;

    // RGB - L1 and L2
    losses[RawLossIndex::RgbL1] = float(mask) * l1_loss(render_rgb, ref_rgb);
    losses[RawLossIndex::RgbL2] = float(mask) *
        clamp(l2_loss(render_rgb, ref_rgb), 0.0f, 1.0f);

    // Depth - Pearson correlation in log space
    depth_mask &= mask;
    render_depth = float(depth_mask) * log(max(render_depth, 1e-4));
    ref_depth = float(depth_mask) * log(max(ref_depth, 1e-4));
    losses[RawLossIndex::DepthSupX] = render_depth;
    losses[RawLossIndex::DepthSupY] = ref_depth;
    losses[RawLossIndex::DepthSupXX] = render_depth * render_depth;
    losses[RawLossIndex::DepthSupYY] = ref_depth * ref_depth;
    losses[RawLossIndex::DepthSupXY] = render_depth * ref_depth;

    // Normal - Supervision and regularization - Pairwise loss between different normals
    normal_mask &= mask;
    bool render_normal_mask = true,
         depth_normal_mask = true;
    render_normal = normalize_normal(render_normal, render_normal_mask);
    depth_normal = normalize_normal(depth_normal, depth_normal_mask);
    ref_normal = normalize_normal(ref_normal, normal_mask);
    losses[RawLossIndex::RenderNormalSup] =
        float(render_normal_mask & normal_mask) * normal_loss(render_normal, ref_normal);
    losses[RawLossIndex::DepthNormalSup] =
        float(depth_normal_mask & normal_mask) * normal_loss(depth_normal, ref_normal);
    losses[RawLossIndex::NormalReg] =
        float(render_normal_mask & depth_normal_mask) * normal_loss(render_normal, depth_normal);

    // Alpha - loss with mask + regularization to go to 0 or 1
    render_alpha = clamp(render_alpha, 0.0f, 1.0f);
    losses[RawLossIndex::AlphaSup] = alpha_loss(render_alpha, float(ref_alpha));
    losses[RawLossIndex::AlphaSupUnder] = alpha_loss(1.0 - render_alpha, 1.0 - float(ref_alpha));
    losses[RawLossIndex::AlphaReg] = 4.0f * render_alpha * (1.0f - render_alpha);  // TODO: support push-to-1 mode

    // Distortion regularization
    losses[RawLossIndex::RgbDistReg] = mean3(rgb_dist) / max(render_alpha, 1e-12);
    losses[RawLossIndex::DepthDistReg] = depth_dist / max(render_alpha, 1e-12);
    losses[RawLossIndex::NormalDistReg] = mean3(normal_dist) / max(render_alpha, 1e-12);

    // Mask total
    losses[RawLossIndex::PixelsTotal] = 1.0f;  // required for patched batching
    losses[RawLossIndex::MaskTotal] = float(mask);
    losses[RawLossIndex::DepthMaskTotal] = float(depth_mask);
    losses[RawLossIndex::RenderNormalMaskTotal] = float(render_normal_mask & normal_mask);
    losses[RawLossIndex::DepthNormalMaskTotal] = float(depth_normal_mask & normal_mask);
    losses[RawLossIndex::NormalRegMaskTotal] = float(render_normal_mask & depth_normal_mask);

    return losses;
}

[CudaDeviceExport]
void per_pixel_losses_bwd(
    float3 render_rgb,
    float3 ref_rgb,
    float render_depth,
    float ref_depth,
    float3 render_normal,
    float3 depth_normal,
    float3 ref_normal,
    float render_alpha,
    float3 rgb_dist,
    float depth_dist,
    float3 normal_dist,
    bool mask,
    bool depth_mask,
    bool normal_mask,
    bool alpha_mask,
    float[RawLossIndex::length] v_losses,
    out float3 v_render_rgb,
    out float3 v_ref_rgb,
    out float v_render_depth,
    out float v_ref_depth,
    out float3 v_render_normal,
    out float3 v_depth_normal,
    out float3 v_ref_normal,
    out float v_render_alpha,
    out float3 v_rgb_dist,
    out float v_depth_dist,
    out float3 v_normal_dist
) {
    DifferentialPair<float3> dp_render_rgb = diffPair(render_rgb);
    DifferentialPair<float3> dp_ref_rgb = diffPair(ref_rgb);
    DifferentialPair<float> dp_render_depth = diffPair(render_depth);
    DifferentialPair<float> dp_ref_depth = diffPair(ref_depth);
    DifferentialPair<float3> dp_render_normal = diffPair(render_normal);
    DifferentialPair<float3> dp_depth_normal = diffPair(depth_normal);
    DifferentialPair<float3> dp_ref_normal = diffPair(ref_normal);
    DifferentialPair<float> dp_render_alpha = diffPair(render_alpha);
    DifferentialPair<float3> dp_rgb_dist = diffPair(rgb_dist);
    DifferentialPair<float> dp_depth_dist = diffPair(depth_dist);
    DifferentialPair<float3> dp_normal_dist = diffPair(normal_dist);
    bwd_diff(per_pixel_losses)(
        dp_render_rgb,
        dp_ref_rgb,
        dp_render_depth,
        dp_ref_depth,
        dp_render_normal,
        dp_depth_normal,
        dp_ref_normal,
        dp_render_alpha,
        dp_rgb_dist,
        dp_depth_dist,
        dp_normal_dist,
        mask,
        depth_mask,
        normal_mask,
        alpha_mask,
        v_losses
    );
    v_render_rgb = dp_render_rgb.d;
    v_ref_rgb = dp_ref_rgb.d;
    v_render_depth = dp_render_depth.d;
    v_ref_depth = dp_ref_depth.d;
    v_render_normal = dp_render_normal.d;
    v_depth_normal = dp_depth_normal.d;
    v_ref_normal = dp_ref_normal.d;
    v_render_alpha = dp_render_alpha.d;
    v_rgb_dist = dp_rgb_dist.d;
    v_depth_dist = dp_depth_dist.d;
    v_normal_dist = dp_normal_dist.d;
}

enum LossWeightIndex {
    RgbSup,
    DepthSup,
    NormalSup,
    AlphaSup,
    AlphaSupUnder,
    NormalReg,
    AlphaReg,
    RgbDistReg,
    DepthDistReg,
    NormalDistReg,
    length
};

enum LossIndex {
    RgbL1,
    RgbPSNR,
    DepthSup,
    NormalSup,
    AlphaSup,
    NormalReg,
    AlphaReg,
    RgbDistReg,
    DepthDistReg,
    NormalDistReg,
    length
};

[Differentiable]
[CudaDeviceExport]
float[LossIndex::length] per_pixel_losses_reduce(
    float[RawLossIndex::length] raw_losses,
    no_diff float[LossWeightIndex::length] weights
) {
    float[LossIndex::length] losses;

    // image loss
    losses[LossIndex::RgbL1] = weights[LossWeightIndex::RgbSup] *
        raw_losses[RawLossIndex::RgbL1] / max(raw_losses[RawLossIndex::MaskTotal], 1);
    losses[LossIndex::RgbPSNR] = -10.0 * log10(
        raw_losses[RawLossIndex::RgbL2] / max(raw_losses[RawLossIndex::MaskTotal], 1));

    // depth supervision loss - Pearson correlation
    losses[LossIndex::DepthSup] = (
        raw_losses[RawLossIndex::DepthMaskTotal] > 0.0 &&
        raw_losses[RawLossIndex::DepthSupY] != 0.0
    ) ?
        weights[LossWeightIndex::DepthSup] * clamp(1.0 -
        (
            raw_losses[RawLossIndex::DepthSupXY] -
                raw_losses[RawLossIndex::DepthSupX] * raw_losses[RawLossIndex::DepthSupY]
                    / raw_losses[RawLossIndex::DepthMaskTotal]
        ) / sqrt(max(1e-12,
            (raw_losses[RawLossIndex::DepthSupXX] -
                raw_losses[RawLossIndex::DepthSupX] * raw_losses[RawLossIndex::DepthSupX]
                    / raw_losses[RawLossIndex::DepthMaskTotal]
            ) * (raw_losses[RawLossIndex::DepthSupYY] -
                raw_losses[RawLossIndex::DepthSupY] * raw_losses[RawLossIndex::DepthSupY]
                    / raw_losses[RawLossIndex::DepthMaskTotal]
            )
         + 1.0))  // numerical issue
    , 0.0, 2.0) : 0.0;

    // normal supervision loss
    losses[LossIndex::NormalSup] = weights[LossWeightIndex::NormalSup] * (
        raw_losses[RawLossIndex::RenderNormalSup] / max(raw_losses[RawLossIndex::RenderNormalMaskTotal], 1) +
        raw_losses[RawLossIndex::DepthNormalSup] / max(raw_losses[RawLossIndex::DepthNormalMaskTotal], 1)
    ) / max(int(raw_losses[RawLossIndex::RenderNormalMaskTotal] > 0.5) +
            int(raw_losses[RawLossIndex::DepthNormalMaskTotal] > 0.5), 1);

    // alpha supervision loss
    losses[LossIndex::AlphaSup] = (
        weights[LossWeightIndex::AlphaSup] * raw_losses[RawLossIndex::AlphaSup] +
        weights[LossWeightIndex::AlphaSupUnder] * raw_losses[RawLossIndex::AlphaSupUnder]
    ) / max(raw_losses[RawLossIndex::PixelsTotal], 1);

    // normal regularization
    losses[LossIndex::NormalReg] = weights[LossWeightIndex::NormalReg] *
        raw_losses[RawLossIndex::NormalReg] / max(raw_losses[RawLossIndex::NormalRegMaskTotal], 1);

    // alpha regularization
    // TODO: push-to-1 mode for random background
    losses[LossIndex::AlphaReg] = weights[LossWeightIndex::AlphaReg] *
        raw_losses[RawLossIndex::AlphaReg] / max(raw_losses[RawLossIndex::PixelsTotal], 1);

    // distortion regularizations
    losses[LossIndex::RgbDistReg] = weights[LossWeightIndex::RgbDistReg] *
        raw_losses[RawLossIndex::RgbDistReg] / max(raw_losses[RawLossIndex::PixelsTotal], 1);
    losses[LossIndex::DepthDistReg] = weights[LossWeightIndex::DepthDistReg] *
        raw_losses[RawLossIndex::DepthDistReg] / max(raw_losses[RawLossIndex::PixelsTotal], 1);
    losses[LossIndex::NormalDistReg] = weights[LossWeightIndex::NormalDistReg] *
        raw_losses[RawLossIndex::NormalDistReg] / max(raw_losses[RawLossIndex::PixelsTotal], 1);

    return losses;
}

[CudaDeviceExport]
float[RawLossIndex::length] per_pixel_losses_reduce_bwd(
    const float[RawLossIndex::length] raw_losses,
    const no_diff float[LossWeightIndex::length] weights,
    float[LossIndex::length] v_losses
) {
    DifferentialPair<float[RawLossIndex::length]> dp_raw_losses = diffPair(raw_losses);
    bwd_diff(per_pixel_losses_reduce)(
        dp_raw_losses, weights,
        v_losses
    );
    return dp_raw_losses.d;
}
