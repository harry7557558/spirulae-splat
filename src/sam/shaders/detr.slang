// DETR decoder support kernels.
//
// Three small pieces of the reference that are not a matmul, a norm or an
// attention, and would otherwise need a host round-trip per decoder layer:
//
//   sine_embed_boxes  gen_sineembed_for_position on the refined reference boxes
//   rpb_deltas        the log-scaled query-to-column offsets the box-relative
//                     positional bias MLP consumes
//   rpb_combine       the outer sum that turns the two per-axis MLP outputs
//                     into the [head][query][key] attention bias
//   pool_weighted     the masked mean over prompt tokens in DotProductScoring
//
// All four run six times per image (once per decoder layer) on tensors of at
// most 8 x 201 x 5184.

#include "_common.slang"

[SpecializationConstant] const int kUnused = 0;

static const int WG = 256;

// ---------------------------------------------------------------------------
// gen_sineembed_for_position
// ---------------------------------------------------------------------------

struct SineEmbedParams {
    float* out;     // [nq, 512]
    uint*  boxes;   // [nq, 4] f32, (cx, cy, w, h)
    uint*  dim_t;   // [64] f32, the 2*pi / 10000^(2i/128) multipliers
    uint   nq;
    uint   groups_per_row;
};

// The reference emits the four coordinates in the order (cy, cx, w, h) -- NOT
// the storage order -- and interleaves sin/cos within each 128-wide block.
// Both details are load-bearing: the ref_point_head MLP that consumes this is
// trained against exactly that permutation.
[shader("compute")]
[numthreads(WG, 1, 1)]
void sine_embed_boxes(uint3 gid : SV_GroupID, uint3 tid : SV_GroupThreadID,
                      uniform SineEmbedParams p) {
    uint i = fold_index(gid, tid, p.groups_per_row, uint(WG));
    if (i >= p.nq * 512u) return;
    uint q = i / 512u;
    uint c = i % 512u;
    uint slot = c / 128u;          // which coordinate
    uint rest = c % 128u;
    uint freq = rest / 2u;
    bool is_cos = (rest & 1u) != 0u;

    // (cy, cx, w, h) from the stored (cx, cy, w, h).
    uint coord_index = (slot == 0u) ? 1u : ((slot == 1u) ? 0u : slot);
    float v = asfloat(p.boxes[q * 4u + coord_index]);
    float angle = v * asfloat(p.dim_t[freq]);
    p.out[i] = is_cos ? cos(angle) : sin(angle);
}

// ---------------------------------------------------------------------------
// Box-relative positional bias
// ---------------------------------------------------------------------------

struct RpbDeltaParams {
    float* out;     // [nq * S, 2]
    uint*  boxes;   // [nq, 4] f32 (cx, cy, w, h)
    uint*  coords;  // [S] f32, i / S
    uint   nq;
    uint   S;
    uint   axis;    // 0 = x (uses cx, w), 1 = y (uses cy, h)
    uint   groups_per_row;
    uint   _pad0;
};

// delta = coord - box_edge, then sign(8d) * log2(|8d| + 1) / log2(8).
// The log squashing is what lets one small MLP cover the whole image extent;
// the constant 8 and base-8 normalization are the reference's.
[shader("compute")]
[numthreads(WG, 1, 1)]
void rpb_deltas(uint3 gid : SV_GroupID, uint3 tid : SV_GroupThreadID,
                uniform RpbDeltaParams p) {
    uint i = fold_index(gid, tid, p.groups_per_row, uint(WG));
    if (i >= p.nq * p.S * 2u) return;
    uint which = i % 2u;           // 0 -> low edge, 1 -> high edge
    uint rest = i / 2u;
    uint s = rest % p.S;
    uint q = rest / p.S;

    float center = asfloat(p.boxes[q * 4u + p.axis]);
    float extent = asfloat(p.boxes[q * 4u + 2u + p.axis]);
    float edge = (which == 0u) ? (center - 0.5f * extent) : (center + 0.5f * extent);
    float d = (asfloat(p.coords[s]) - edge) * 8.0f;
    float sgn = (d > 0.0f) ? 1.0f : ((d < 0.0f) ? -1.0f : 0.0f);
    float lg = log(abs(d) + 1.0f) / log(2.0f);
    p.out[i] = sgn * lg * (1.0f / 3.0f);   // 1 / log2(8)
}

struct RpbCombineParams {
    float* out;     // [n_heads, nq + 1, S*S]
    uint*  rpb_x;   // [nq * S, n_heads]
    uint*  rpb_y;   // [nq * S, n_heads]
    uint   nq;
    uint   S;
    uint   n_heads;
    uint   groups_per_row;
    uint   _pad0;
};

// bias[h][q][hy * S + hx] = rpb_y[q][hy][h] + rpb_x[q][hx][h], with query row 0
// (the presence token) left at zero -- it has no reference box.
[shader("compute")]
[numthreads(WG, 1, 1)]
void rpb_combine(uint3 gid : SV_GroupID, uint3 tid : SV_GroupThreadID,
                 uniform RpbCombineParams p) {
    uint i = fold_index(gid, tid, p.groups_per_row, uint(WG));
    uint nkv = p.S * p.S;
    uint nq1 = p.nq + 1u;
    if (i >= p.n_heads * nq1 * nkv) return;
    uint kv = i % nkv;
    uint rest = i / nkv;
    uint qi = rest % nq1;
    uint h = rest / nq1;
    if (qi == 0u) {
        p.out[i] = 0.0f;
        return;
    }
    uint q = qi - 1u;
    uint hy = kv / p.S;
    uint hx = kv % p.S;
    float vy = asfloat(p.rpb_y[(q * p.S + hy) * p.n_heads + h]);
    float vx = asfloat(p.rpb_x[(q * p.S + hx) * p.n_heads + h]);
    p.out[i] = vy + vx;
}

// ---------------------------------------------------------------------------
// Masked mean over prompt tokens (DotProductScoring)
// ---------------------------------------------------------------------------

struct PoolParams {
    float* out;      // [C]
    uint*  x;        // [T, C]
    uint*  weights;  // [T] f32; already 1/n_valid for valid tokens, 0 for padding
    uint   T;
    uint   C;
    uint   groups_per_row;
};

[shader("compute")]
[numthreads(WG, 1, 1)]
void pool_weighted(uint3 gid : SV_GroupID, uint3 tid : SV_GroupThreadID,
                   uniform PoolParams p) {
    uint c = fold_index(gid, tid, p.groups_per_row, uint(WG));
    if (c >= p.C) return;
    float acc = 0.0f;
    for (uint t = 0; t < p.T; ++t)
        acc += asfloat(p.x[t * p.C + c]) * asfloat(p.weights[t]);
    p.out[c] = acc;
}

// ---------------------------------------------------------------------------
// Iterative box refinement: boxes = sigmoid(inverse_sigmoid(boxes) + delta)
// ---------------------------------------------------------------------------

struct BoxRefineParams {
    float* boxes;   // [nq, 4], updated in place
    uint*  delta;   // [nq, 4]
    uint   n;
    uint   groups_per_row;
};

[shader("compute")]
[numthreads(WG, 1, 1)]
void box_refine(uint3 gid : SV_GroupID, uint3 tid : SV_GroupThreadID,
                uniform BoxRefineParams p) {
    uint i = fold_index(gid, tid, p.groups_per_row, uint(WG));
    if (i >= p.n) return;
    // eps = 1e-3 on both sides, matching the reference's inverse_sigmoid.
    float x = clamp(p.boxes[i], 1e-3f, 1.0f - 1e-3f);
    float inv = log(x) - log(1.0f - x);
    p.boxes[i] = 1.0f / (1.0f + exp(-(inv + asfloat(p.delta[i]))));
}
