// Densification kernels (Vulkan backend), mirroring csrc/Densify.cu.
// Currently only densify_update_weight (the per-step score accumulator the
// optimizer path calls); splitting/relocation/MCMC land with the
// densify/MCMC port.
//
// Score modes match DensifyScoreMode in csrc/Densify.cuh:
//   0 = Mean, 1 = Max, 2 = Median, 3 = Geom.
// The Median mode's multiplicative random walk uses hash14 (shadertoy
// 4djSRW); its float arithmetic is contraction-sensitive, so CUDA and
// Vulkan can draw different hashes for the same input — parity tests
// compare only the deterministic count channel for Median.

[ForceInline]
float _dn_hash14(float4 p4) {
    p4 = p4 * float4(0.1031f, 0.1030f, 0.0973f, 0.1099f);
    p4 = float4(fmod(p4.x, 1.0f), fmod(p4.y, 1.0f), fmod(p4.z, 1.0f),
                fmod(p4.w, 1.0f));
    p4 = p4 + dot(p4, p4.wzxy + 33.33f);
    return fmod((p4.x + p4.y) * (p4.z + p4.w), 1.0f);
}

[ForceInline]
float _dn_sigmoid(float x) { return 1.0f / (1.0f + exp(-x)); }

[ForceInline]
bool _dn_finitef(float v) {
    return (asuint(v) & 0x7f800000u) != 0x7f800000u;
}

// Optional inputs are flag-gated (use_*) with fallback-backed pointers —
// never null-pointer-gated (see KernelCommon.h null_fallback). The CUDA
// kernel's accum_weight_scalar input is dropped: the launch API always
// passes null for it.
struct DensifyUpdateWeightParams {  // 64 bytes: pushed directly
    float* radii;         // [N]
    float* opacs;         // [N] logit-space (use_opacs gates)
    float* accum_weight;  // [N]
    float* accum_weight2; // [N] blend partner (use_w2 gates)
    float2* accum_buffer; // [N, 2] (score, count)
    float blend_w;
    int score_mode;
    uint use_opacs;
    uint use_w2;
    uint num_splats;
    uint wgs_per_row;
};

[shader("compute")]
[numthreads(256, 1, 1)]
void densify_update_weight(uint3 wg: SV_GroupID, uint tid: SV_GroupThreadID,
                           uniform DensifyUpdateWeightParams p) {
    uint idx = (wg.y * p.wgs_per_row + wg.x) * 256 + tid;
    if (idx >= p.num_splats) return;
    if (!(p.radii[idx] > 0.0f)) return;

    float weight = abs(p.accum_weight[idx]);
    // Geometric blend: weight^(1-w) * |accum_weight2|^w (see Densify.cu).
    if (p.use_w2 != 0) {
        weight = pow(weight, 1.0f - p.blend_w) *
                 pow(abs(p.accum_weight2[idx]), p.blend_w);
        if (!_dn_finitef(weight)) weight = 0.0f;
    }
    if (p.use_opacs != 0) weight *= _dn_sigmoid(p.opacs[idx]);
    if (weight == 0.0f) return;

    float2 accum = p.accum_buffer[idx];
    if (p.score_mode == 1) {  // Max
        accum.x = max(accum.x, weight);
        accum.y = max(accum.y, 1.0f);
    } else if (p.score_mode == 0) {  // Mean
        accum.x = (accum.x * accum.y + weight) / (accum.y + 1.0f);
        accum.y += 1.0f;
    } else if (p.score_mode == 2) {  // Median
        float rand = _dn_hash14(
            1e2f * float4(weight, accum.x, accum.y, accum.y + 1.5f));
        if (accum.y == 0.0f) {
            accum.x = weight * exp2(rand - 0.5f);
            accum.y = 1.0f;
        } else if (weight != 0.0f) {
            float s = weight > accum.x ? 1.0f : -1.0f;
            accum.x *= exp2(s * rand);
            accum.y += 1.0f;
        }
    } else if (p.score_mode == 3) {  // Geom
        accum.x = exp((log(max(accum.x, 1e-30f)) * accum.y +
                       log(max(weight, 1e-30f))) /
                      (accum.y + 1.0f));
        accum.y += 1.0f;
    }
    p.accum_buffer[idx] = accum;
}
