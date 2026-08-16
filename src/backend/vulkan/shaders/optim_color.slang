// Trust-region Adam color kernels (Vulkan backend), mirroring the
// fused_adamtr_rgb_optim / fused_adamtr_rgb_sh_optim kernels in
// kernels/optim/ColorOptim.cu (the non-FPBO DC/SH update path when
// use_color_trust_region is on and no quantization is active). The
// is_linear (sRGB-Jacobian gradient conversion) and zero_grad template
// axes become runtime flags — the kernels are branch-light and have no
// barriers.

#define CudaDeviceExport ForceInline
#include "shaders/pixel_wise.slang"
#undef CudaDeviceExport

static const float kTrBeta1 = 0.9f;
static const float kTrBeta2 = 0.999f;
static const float kTrSh0 = 0.28209479177387814f;  // 1 / (2*sqrt(pi))

[ForceInline]
float _tr_sigmoid(float x) { return 1.0f / (1.0f + exp(-x)); }

struct AdamTrRgbParams {  // 88 bytes: pushed directly
    float* rgbs;        // [N, 3] (features_dc)
    float* grad;        // [N, 3]
    float* exp_avg;     // [N, 3]
    float* exp_avg_sq;  // [N, 3]
    float* opacities;   // [N] logit-space
    float lr;
    float bias_correction1;  // 1 - beta1^step
    float bias_correction2;  // 1 - beta2^step
    float eps;
    float eps_tr;
    float grad_scale;
    uint is_linear;
    uint zero_grad;
    uint num_gs;
    uint wgs_per_row;
};

[shader("compute")]
[numthreads(256, 1, 1)]
void fused_adamtr_rgb(uint3 wg: SV_GroupID, uint tid: SV_GroupThreadID,
                      uniform AdamTrRgbParams p) {
    uint idx = (wg.y * p.wgs_per_row + wg.x) * 256 + tid;
    if (idx >= p.num_gs) return;

    float step_size = p.lr / p.bias_correction1;

    float3 x = float3(p.rgbs[3 * idx], p.rgbs[3 * idx + 1],
                      p.rgbs[3 * idx + 2]);
    float3 v = float3(p.grad[3 * idx], p.grad[3 * idx + 1],
                      p.grad[3 * idx + 2]);
    v = float3(p.grad_scale, p.grad_scale, p.grad_scale) * v;
    float3 m_val = float3(p.exp_avg[3 * idx], p.exp_avg[3 * idx + 1],
                          p.exp_avg[3 * idx + 2]);
    float3 v_val = float3(p.exp_avg_sq[3 * idx], p.exp_avg_sq[3 * idx + 1],
                          p.exp_avg_sq[3 * idx + 2]);

    // Carry the gradient to x = splat_dc_encode(dc), the space Adam and its
    // moments live in when the working colour space is linear.
    if (p.is_linear != 0) {
        v.x /= linear_rgb_to_srgb_grad(kTrSh0 * x.x + 0.5f);
        v.y /= linear_rgb_to_srgb_grad(kTrSh0 * x.y + 0.5f);
        v.z /= linear_rgb_to_srgb_grad(kTrSh0 * x.z + 0.5f);
    }

    m_val = kTrBeta1 * m_val + (1.0f - kTrBeta1) * v;
    v_val = kTrBeta2 * v_val + (1.0f - kTrBeta2) * v * v;

    float3 denom = float3(sqrt(max(v_val.x / p.bias_correction2, 0.0f)),
                          sqrt(max(v_val.y / p.bias_correction2, 0.0f)),
                          sqrt(max(v_val.z / p.bias_correction2, 0.0f))) +
                   float3(p.eps, p.eps, p.eps);
    float3 delta = float3(-step_size, -step_size, -step_size) *
                   (m_val / denom);

    // Trust region. splat_dc_encode already scales the step by brightness, so
    // the linear path clips a bare radius (the 2 makes it identical to the
    // colour-proportional one wherever that binds).
    float opac = _tr_sigmoid(p.opacities[idx]);
    float inv_opac = 1.0f / max(opac, 1e-12f);
    float3 clipv;
    if (p.is_linear != 0) {
        clipv = float3(kTrSh0 * sqrt(max(2.0f * p.eps_tr * inv_opac, 0.0f)));
    } else {
        float3 rgb = max(kTrSh0 * x + float3(0.5f),
                         float3((1.0f / 255.0f) * (1.0f / 255.0f)));
        clipv =
            kTrSh0 * float3(sqrt(max(4.0f * p.eps_tr * rgb.x * inv_opac, 0.0f)),
                            sqrt(max(4.0f * p.eps_tr * rgb.y * inv_opac, 0.0f)),
                            sqrt(max(4.0f * p.eps_tr * rgb.z * inv_opac, 0.0f)));
    }

    delta = clamp(delta, -clipv, clipv);
    delta.x = isfinite(delta.x) ? delta.x : 0.0f;
    delta.y = isfinite(delta.y) ? delta.y : 0.0f;
    delta.z = isfinite(delta.z) ? delta.z : 0.0f;

    float3 new_dc = x + delta;
    if (p.is_linear != 0) {
        new_dc = float3(splat_dc_decode(splat_dc_encode(x.x) + delta.x),
                        splat_dc_decode(splat_dc_encode(x.y) + delta.y),
                        splat_dc_decode(splat_dc_encode(x.z) + delta.z));
    }
    p.rgbs[3 * idx] = new_dc.x;
    p.rgbs[3 * idx + 1] = new_dc.y;
    p.rgbs[3 * idx + 2] = new_dc.z;
    p.exp_avg[3 * idx] = m_val.x;
    p.exp_avg[3 * idx + 1] = m_val.y;
    p.exp_avg[3 * idx + 2] = m_val.z;
    p.exp_avg_sq[3 * idx] = v_val.x;
    p.exp_avg_sq[3 * idx + 1] = v_val.y;
    p.exp_avg_sq[3 * idx + 2] = v_val.z;

    // Last, per optimizer.slang's fused_adam_fwd.
    if (p.zero_grad != 0) {
        p.grad[3 * idx] = 0.0f;
        p.grad[3 * idx + 1] = 0.0f;
        p.grad[3 * idx + 2] = 0.0f;
    }
}

struct AdamTrRgbShParams {  // 96 bytes: pushed directly
    float* param;       // [N, K, 3] (features_sh), scalar-indexed
    float* grad;
    float* exp_avg;
    float* exp_avg_sq;
    float* rgbs;        // [N, 3] (features_dc)
    float* opacities;   // [N] logit-space
    float lr;
    float bias_correction1;
    float bias_correction2;
    float eps;
    float eps_tr;
    float grad_scale;
    uint is_linear;
    uint zero_grad;
    uint num_params;  // num_gs * num_sh * 3, scalar-indexed
    uint num_sh;
    uint wgs_per_row;
};

[shader("compute")]
[numthreads(256, 1, 1)]
void fused_adamtr_rgb_sh(uint3 wg: SV_GroupID, uint tid: SV_GroupThreadID,
                         uniform AdamTrRgbShParams p) {
    uint idx = (wg.y * p.wgs_per_row + wg.x) * 256 + tid;
    if (idx >= p.num_params) return;

    float step_size = p.lr / p.bias_correction1;

    float v = p.grad[idx];
    v *= p.grad_scale;
    float m_val = p.exp_avg[idx];
    float v_val = p.exp_avg_sq[idx];

    // SH stays in coefficient space -- the DC reparameterization does not
    // extend here, because linearizing it about the DC colour needs
    // |sh| << c, which 3DGS does not satisfy.
    float c = p.rgbs[(idx / (3 * p.num_sh)) * 3 + (idx % 3)] * kTrSh0 + 0.5f;
    if (p.is_linear != 0) v /= linear_rgb_to_srgb_grad(c);

    m_val = kTrBeta1 * m_val + (1.0f - kTrBeta1) * v;
    v_val = kTrBeta2 * v_val + (1.0f - kTrBeta2) * v * v;

    float denom = sqrt(max(v_val / p.bias_correction2, 0.0f)) + p.eps;
    float delta = -step_size * (m_val / denom);

    // Trust region
    float opac = _tr_sigmoid(p.opacities[idx / (3 * p.num_sh)]);
    c = max(c, (1.0f / 255.0f) * (1.0f / 255.0f));
    float clip = kTrSh0 * sqrt(max(4.0f * p.eps_tr * c / max(opac, 1e-12f),
                                   0.0f));

    delta = clamp(delta, -clip, clip);
    p.param[idx] += delta;
    p.exp_avg[idx] = m_val;
    p.exp_avg_sq[idx] = v_val;

    // Last, per optimizer.slang's fused_adam_fwd.
    if (p.zero_grad != 0) p.grad[idx] = 0.0f;
}
