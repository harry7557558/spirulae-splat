// Metric3D's model-specific kernels.
//
// Everything else in the forward pass is a general op from nn/ -- convolution,
// transposed convolution, GEMM, LayerNorm, attention, resize, pooling, row
// softmax. What is left here is what could not be anything else: the RAFT
// convex upsampler's weighting rule, the GRU's update blend, and the two
// places where the network's six-channel state is taken apart into a depth and
// a normal.

#include "_common.slang"

// Spec ids 0..1.
[SpecializationConstant] const int kInputF16 = 0;
[SpecializationConstant] const int kUpFactor = 4;

static const int WG = 256;

// ---------------------------------------------------------------------------
// Convex upsampling
//
// RAFT's rule: every output pixel is a convex combination of the 3x3
// neighbourhood of its source pixel, with the nine weights predicted per
// (source pixel, sub-position) and softmaxed.
//
// The mask arrives channel-last as [H, W, 9*f*f] with channel index
// (k*f + i)*f + j -- torch's view(N, 1, 9, f, f, H, W) read through our layout.
// The softmax is over k, whose stride is therefore f*f and NOT contiguous,
// which is why nn::softmax_rows cannot do this half and the whole thing is one
// kernel. Doing it per output pixel rather than per output element also means
// the nine exponentials are computed once for all D channels.
//
// The 0.25 the reference applies to the mask before the softmax is folded into
// the caller's conv, not here -- softmax is not scale-invariant, so it has to
// happen on the logits either way.
// ---------------------------------------------------------------------------

struct ConvexUpParams {
    float* out;    // [H*f, W*f, D]
    uint*  flow;   // [H, W, D]
    uint*  mask;   // [H, W, 9*f*f]
    uint   H;
    uint   W;
    uint   D;
    uint   groups_per_row;
};

[shader("compute")]
[numthreads(WG, 1, 1)]
void convex_upsample(uint3 gid : SV_GroupID, uint3 tid : SV_GroupThreadID,
                     uniform ConvexUpParams p) {
    const uint f = uint(kUpFactor);
    uint i = fold_index(gid, tid, p.groups_per_row, uint(WG));
    if (i >= p.H * p.W * f * f) return;

    // Output pixel (yo, xo) = (y*f + si, x*f + sj).
    uint xo = i % (p.W * f);
    uint yo = i / (p.W * f);
    uint x = xo / f, sj = xo % f;
    uint y = yo / f, si = yo % f;

    const uint mask_stride = f * f;
    const uint mask_base = (y * p.W + x) * 9u * mask_stride + si * f + sj;

    float w[9];
    float hi = -3.402823466e+38f;
    for (uint k = 0; k < 9u; ++k) {
        w[k] = load_elem(p.mask, mask_base + k * mask_stride, kInputF16 != 0);
        hi = max(hi, w[k]);
    }
    float sum = 0.0f;
    for (uint k = 0; k < 9u; ++k) { w[k] = exp(w[k] - hi); sum += w[k]; }
    const float inv = 1.0f / sum;

    // F.unfold(flow, [3,3], padding=1) order: k = ky*3 + kx, offsets -1..1,
    // zero outside.
    for (uint d = 0; d < p.D; ++d) {
        float acc = 0.0f;
        for (uint k = 0; k < 9u; ++k) {
            int sy = int(y) + int(k / 3u) - 1;
            int sx = int(x) + int(k % 3u) - 1;
            if (sy < 0 || sy >= int(p.H) || sx < 0 || sx >= int(p.W)) continue;
            acc += w[k] * inv *
                   load_elem(p.flow, (uint(sy) * p.W + uint(sx)) * p.D + d, kInputF16 != 0);
        }
        p.out[i * p.D + d] = acc;
    }
}

// ---------------------------------------------------------------------------
// The six-channel state
//
// (depth, depth confidence, normal xyz, normal concentration). Both kernels
// below carry norm_normalize, which is Metric3D's own and has two details that
// look like rounding and are not: the epsilon is ADDED to the norm rather than
// clamped against, and the concentration channel is elu(k) + 1 + 0.01 rather
// than a plain activation.
// ---------------------------------------------------------------------------

struct PackParams {
    float* out;      // [H, W, 6]
    uint*  expect;   // [H, W, 1]   depth expectation over the anchors
    uint*  conf;     // [H, W, C]   the FuseBlock tail; last two channels used
    uint*  normal;   // [H, W, 3]   the normal head, before normalization
    uint   n;        // H * W
    uint   conf_ch;  // C
    float  min_depth;
    float  max_depth;
    float  regress_scale;
    uint   groups_per_row;
};

// depth_init = cat(regress_depth, depth_confidence, norm_normalize(normal))
[shader("compute")]
[numthreads(WG, 1, 1)]
void pack_state(uint3 gid : SV_GroupID, uint3 tid : SV_GroupThreadID,
                uniform PackParams p) {
    uint i = fold_index(gid, tid, p.groups_per_row, uint(WG));
    if (i >= p.n) return;

    float d = load_elem(p.expect, i, kInputF16 != 0);
    d = clamp(d, p.min_depth, p.max_depth);

    float nx = load_elem(p.normal, i * 3u + 0u, kInputF16 != 0);
    float ny = load_elem(p.normal, i * 3u + 1u, kInputF16 != 0);
    float nz = load_elem(p.normal, i * 3u + 2u, kInputF16 != 0);
    float inv = 1.0f / (sqrt(nx * nx + ny * ny + nz * nz) + 1e-10f);

    float depth_conf = load_elem(p.conf, i * p.conf_ch + p.conf_ch - 2u, kInputF16 != 0);
    float kappa = load_elem(p.conf, i * p.conf_ch + p.conf_ch - 1u, kInputF16 != 0);

    p.out[i * 6u + 0u] = (d - p.max_depth) / p.regress_scale;
    p.out[i * 6u + 1u] = depth_conf;
    p.out[i * 6u + 2u] = nx * inv;
    p.out[i * 6u + 3u] = ny * inv;
    p.out[i * 6u + 4u] = nz * inv;
    p.out[i * 6u + 5u] = apply_act(kappa, ACT_ELU) + 1.01f;
}

struct UnpackParams {
    float* depth;    // [H, W, 1], or null
    float* normal;   // [H, W, 3], or null
    float* kappa;    // [H, W, 1], or null
    uint*  state;    // [H, W, 6]
    uint   n;
    float  min_depth;
    float  max_depth;
    float  regress_scale;
    uint   groups_per_row;
};

[shader("compute")]
[numthreads(WG, 1, 1)]
void unpack_state(uint3 gid : SV_GroupID, uint3 tid : SV_GroupThreadID,
                  uniform UnpackParams p) {
    uint i = fold_index(gid, tid, p.groups_per_row, uint(WG));
    if (i >= p.n) return;

    float d = load_elem(p.state, i * 6u + 0u, kInputF16 != 0) * p.regress_scale +
              p.max_depth;
    p.depth[i] = clamp(d, p.min_depth, p.max_depth);

    float nx = load_elem(p.state, i * 6u + 2u, kInputF16 != 0);
    float ny = load_elem(p.state, i * 6u + 3u, kInputF16 != 0);
    float nz = load_elem(p.state, i * 6u + 4u, kInputF16 != 0);
    float inv = 1.0f / (sqrt(nx * nx + ny * ny + nz * nz) + 1e-10f);
    p.normal[i * 3u + 0u] = nx * inv;
    p.normal[i * 3u + 1u] = ny * inv;
    p.normal[i * 3u + 2u] = nz * inv;

    p.kappa[i] = apply_act(load_elem(p.state, i * 6u + 5u, kInputF16 != 0), ACT_ELU) +
                 1.01f;
}

// ---------------------------------------------------------------------------
// ConvGRU update
//
// h = (1 - z) * h + z * q, in place. Written as h + z * (q - h) so `h` is read
// once; the two forms agree to fp32 rounding and this one is a third of the
// memory traffic of doing it with three general elementwise passes.
// ---------------------------------------------------------------------------

struct BlendParams {
    float* h;
    uint*  z;
    uint*  q;
    uint   n;
    uint   groups_per_row;
};

[shader("compute")]
[numthreads(WG, 1, 1)]
void gru_blend(uint3 gid : SV_GroupID, uint3 tid : SV_GroupThreadID,
               uniform BlendParams p) {
    uint i = fold_index(gid, tid, p.groups_per_row, uint(WG));
    if (i >= p.n) return;
    float h = p.h[i];
    float z = load_elem(p.z, i, kInputF16 != 0);
    float q = load_elem(p.q, i, kInputF16 != 0);
    p.h[i] = h + z * (q - h);
}
