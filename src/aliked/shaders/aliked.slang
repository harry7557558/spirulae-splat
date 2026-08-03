// ALIKED's model-specific kernels: the detector's NMS and sub-pixel
// refinement, and the two coordinate conversions its descriptor head needs.
//
// Everything else in the forward pass is a general op from nn/ -- convolution,
// deformable convolution, pooling, resize, grid sample, GEMM, L2 normalize.
// What is left here is what could not be anything else: a suppression rule, a
// soft-argmax, and the layout conversions between them.

#include "_common.slang"

// Spec ids 0..1.
[SpecializationConstant] const int kInputF16 = 0;
[SpecializationConstant] const int kAct      = 0;

static const int WG = 256;

// ---------------------------------------------------------------------------
// Non-maximum suppression
//
// A faithful port of ALIKED's `simple_nms`, which is NOT "score == maxpool":
//
//   max_mask = scores == max_pool(scores)
//   twice:
//     supp_mask   = max_pool(max_mask) > 0
//     supp_scores = where(supp_mask, 0, scores)
//     max_mask   |= (supp_scores == max_pool(supp_scores)) & ~supp_mask
//
// The two extra rounds recover a peak that a *neighbouring* peak masked, which
// on a smooth score map is a few percent of all keypoints -- not a rounding
// difference, a different feature set. The five max_pool calls this expands to
// are exactly the five MaxPool nodes in the exported graph, which is how the
// port was checked against it.
//
// Split across two entry points per round because supp_mask has to be
// materialized: computing whether a *neighbour* is suppressed needs the whole
// window of max_mask, so it cannot be folded into the comparison that reads it.
// ---------------------------------------------------------------------------

struct NmsParams {
    uint*  mask;      // [H, W] u32, 0/1 -- read and written across rounds
    uint*  supp;      // [H, W] u32, 0/1 scratch
    uint*  scores;    // [H, W] f32
    uint   H;
    uint   W;
    int    radius;
    uint   groups_per_row;
};

// max_mask = (scores == max over the (2r+1)^2 window)
[shader("compute")]
[numthreads(WG, 1, 1)]
void nms_init(uint3 gid : SV_GroupID, uint3 tid : SV_GroupThreadID, uniform NmsParams p) {
    uint i = fold_index(gid, tid, p.groups_per_row, uint(WG));
    if (i >= p.H * p.W) return;
    int x = int(i % p.W), y = int(i / p.W);

    float s = asfloat(p.scores[i]);
    float best = -3.402823466e+38f;
    for (int dy = -p.radius; dy <= p.radius; ++dy) {
        int yy = y + dy;
        if (yy < 0 || yy >= int(p.H)) continue;
        for (int dx = -p.radius; dx <= p.radius; ++dx) {
            int xx = x + dx;
            if (xx < 0 || xx >= int(p.W)) continue;
            best = max(best, asfloat(p.scores[uint(yy) * p.W + uint(xx)]));
        }
    }
    p.mask[i] = (s == best) ? 1u : 0u;
}

// supp = max_pool(max_mask) > 0, i.e. "some accepted peak lies within radius".
[shader("compute")]
[numthreads(WG, 1, 1)]
void nms_suppress(uint3 gid : SV_GroupID, uint3 tid : SV_GroupThreadID,
                  uniform NmsParams p) {
    uint i = fold_index(gid, tid, p.groups_per_row, uint(WG));
    if (i >= p.H * p.W) return;
    int x = int(i % p.W), y = int(i / p.W);

    uint any = 0u;
    for (int dy = -p.radius; dy <= p.radius && any == 0u; ++dy) {
        int yy = y + dy;
        if (yy < 0 || yy >= int(p.H)) continue;
        for (int dx = -p.radius; dx <= p.radius; ++dx) {
            int xx = x + dx;
            if (xx < 0 || xx >= int(p.W)) continue;
            if (p.mask[uint(yy) * p.W + uint(xx)] != 0u) { any = 1u; break; }
        }
    }
    p.supp[i] = any;
}

// max_mask |= (supp_scores == max_pool(supp_scores)) & ~supp
[shader("compute")]
[numthreads(WG, 1, 1)]
void nms_recover(uint3 gid : SV_GroupID, uint3 tid : SV_GroupThreadID,
                 uniform NmsParams p) {
    uint i = fold_index(gid, tid, p.groups_per_row, uint(WG));
    if (i >= p.H * p.W) return;
    if (p.supp[i] != 0u) return;   // ~supp_mask gates the whole term
    int x = int(i % p.W), y = int(i / p.W);

    // supp_scores is `where(supp, 0, scores)`, read rather than materialized.
    float s = asfloat(p.scores[i]);
    float best = -3.402823466e+38f;
    for (int dy = -p.radius; dy <= p.radius; ++dy) {
        int yy = y + dy;
        if (yy < 0 || yy >= int(p.H)) continue;
        for (int dx = -p.radius; dx <= p.radius; ++dx) {
            int xx = x + dx;
            if (xx < 0 || xx >= int(p.W)) continue;
            uint j = uint(yy) * p.W + uint(xx);
            best = max(best, (p.supp[j] != 0u) ? 0.0f : asfloat(p.scores[j]));
        }
    }
    if (s == best) p.mask[i] = 1u;
}

// ---------------------------------------------------------------------------
// Candidate compaction
// ---------------------------------------------------------------------------

struct CollectParams {
    uint*  out;       // [cap, 3]: x, y (u32) and the score's bits
    uint*  counter;   // [1]
    uint*  mask;
    uint*  scores;
    uint   H;
    uint   W;
    int    border;    // DKD zeroes this many rows/cols at the top and left,
                      // and border-1 at the bottom and right (see below)
    float  min_score; // peaks below this never become candidates
    uint   cap;
    uint   groups_per_row;
};

// Appends every surviving peak, with its integer position and score.
//
// The border rule is ALIKED's and is asymmetric -- it zeroes rows [0, r] and
// columns [0, r] but only rows [h-r, h) and columns [w-r, w). So a keypoint
// satisfies r < y < h-r, which is also exactly the condition under which its
// (2r+1)^2 refinement patch is fully inside the map. The refinement kernel
// relies on that and does no bounds testing.
[shader("compute")]
[numthreads(WG, 1, 1)]
void nms_collect(uint3 gid : SV_GroupID, uint3 tid : SV_GroupThreadID,
                 uniform CollectParams p) {
    uint i = fold_index(gid, tid, p.groups_per_row, uint(WG));
    if (i >= p.H * p.W) return;
    if (p.mask[i] == 0u) return;
    // The threshold applies to the PEAK's own score here, and again to the
    // sub-pixel-refined score after refinement. Both, because that is what the
    // reference does: a peak scoring just under the threshold whose refined
    // score lands just over it is dropped, not kept. Filtering here also keeps
    // the candidate list bounded on a textured image.
    if (asfloat(p.scores[i]) < p.min_score) return;
    int x = int(i % p.W), y = int(i / p.W);
    if (x <= p.border || y <= p.border) return;
    if (x >= int(p.W) - p.border || y >= int(p.H) - p.border) return;

    uint slot;
    InterlockedAdd(p.counter[0], 1u, slot);
    if (slot >= p.cap) return;   // the host reads the counter and warns
    p.out[slot * 3u + 0u] = uint(x);
    p.out[slot * 3u + 1u] = uint(y);
    p.out[slot * 3u + 2u] = p.scores[i];
}

// ---------------------------------------------------------------------------
// Sub-pixel refinement
// ---------------------------------------------------------------------------

struct RefineParams {
    float* out_xy;     // [N, 2] refined pixel coordinates, (x, y)
    float* out_score;  // [N]
    uint*  centers;    // [N, 2] i32 integer peak positions, (x, y)
    uint*  scores;     // [H, W] f32
    uint   H;
    uint   W;
    uint   N;
    int    radius;
    float  temperature;
    uint   groups_per_row;
};

// ALIKED's DKD soft-argmax: over the (2r+1)^2 patch of the *score map* around
// the integer peak,
//
//   e   = exp((patch - max(patch)) / T)
//   d   = sum(e * grid) / sum(e)        grid = the (dx, dy) offsets, -r..r
//   xy  = peak + d
//
// and the reported score is the score map resampled *bilinearly at the refined
// position*, not the peak's own value -- COLMAP then thresholds on that one, so
// using the peak value here would filter a different set.
//
// One thread per keypoint: the patch is 25 values and there are a few thousand
// of them, so a workgroup-per-keypoint reduction would leave most lanes idle.
[shader("compute")]
[numthreads(WG, 1, 1)]
void refine_keypoints(uint3 gid : SV_GroupID, uint3 tid : SV_GroupThreadID,
                      uniform RefineParams p) {
    uint n = fold_index(gid, tid, p.groups_per_row, uint(WG));
    if (n >= p.N) return;

    int cx = int(p.centers[2u * n]);
    int cy = int(p.centers[2u * n + 1u]);

    // The border rule guarantees the patch is inside; no clamping, so a
    // violation would be a bug worth crashing on rather than hiding.
    float mx = -3.402823466e+38f;
    for (int dy = -p.radius; dy <= p.radius; ++dy)
        for (int dx = -p.radius; dx <= p.radius; ++dx)
            mx = max(mx, asfloat(p.scores[uint(cy + dy) * p.W + uint(cx + dx)]));

    float sum = 0.0f, sx = 0.0f, sy = 0.0f;
    for (int dy = -p.radius; dy <= p.radius; ++dy)
        for (int dx = -p.radius; dx <= p.radius; ++dx) {
            float e = exp((asfloat(p.scores[uint(cy + dy) * p.W + uint(cx + dx)]) - mx) /
                          p.temperature);
            sum += e;
            sx += e * float(dx);
            sy += e * float(dy);
        }
    float inv = 1.0f / max(sum, 1e-20f);
    float rx = float(cx) + sx * inv;
    float ry = float(cy) + sy * inv;

    p.out_xy[2u * n] = rx;
    p.out_xy[2u * n + 1u] = ry;
    // grid_sample(align_corners=True) at the refined point: in pixel space
    // that is a plain bilinear fetch, and the point is interior by
    // construction, so clamped and zero-padded sampling agree here.
    p.out_score[n] = sample_bilinear(p.scores, false, int(p.H), int(p.W), 1, 0, ry, rx);
}

// ---------------------------------------------------------------------------
// Descriptor-head coordinate conversions
// ---------------------------------------------------------------------------

struct CenterParams {
    uint*  out;    // [N, 2] i32
    uint*  xy;     // [N, 2] f32 pixel coordinates
    uint   N;
    uint   groups_per_row;
};

// The patch centre SDDH gathers around is `kptsi_wh.long()` -- a truncation of
// the *refined* sub-pixel position, which is not always the integer peak the
// refinement started from. Doing this on the host would mean a round trip
// between the refinement and the descriptor head for no reason.
[shader("compute")]
[numthreads(WG, 1, 1)]
void to_int_centers(uint3 gid : SV_GroupID, uint3 tid : SV_GroupThreadID,
                    uniform CenterParams p) {
    uint i = fold_index(gid, tid, p.groups_per_row, uint(WG));
    if (i >= p.N * 2u) return;
    // Coordinates are non-negative here (the border rule keeps them >= 1), so
    // torch's truncation and a floor agree.
    p.out[i] = uint(int(asfloat(p.xy[i])));
}

struct SddhPosParams {
    float* out;        // [N * M, 2] normalized (gx, gy) for grid_sample
    uint*  xy;         // [N, 2] f32 keypoint pixel coordinates
    uint*  offsets;    // [N, 2 * M] f32, raw offset-head output
    uint   N;
    uint   M;
    uint   W;
    uint   H;
    float  max_offset;
    uint   groups_per_row;
};

// SDDH's sampling positions.
//
// The offset head emits 2*M channels which ALIKED reads as `view(N, 2, M)`:
// the FIRST M channels are the x components and the SECOND M the y ones. That
// is not the layout of the deformable convolution's offsets in the same
// network, which are (dy, dx) interleaved per tap. Two conventions, one model;
// swapping either produces descriptors that look normal and match badly.
//
// Positions are then normalized on the align_corners=True convention --
// 2*pos/(w-1, h-1) - 1 -- because that is what ALIKED's grid_sample call uses.
[shader("compute")]
[numthreads(WG, 1, 1)]
void sddh_positions(uint3 gid : SV_GroupID, uint3 tid : SV_GroupThreadID,
                    uniform SddhPosParams p) {
    uint i = fold_index(gid, tid, p.groups_per_row, uint(WG));
    if (i >= p.N * p.M) return;
    uint m = i % p.M;
    uint n = i / p.M;

    float ox = asfloat(p.offsets[n * 2u * p.M + m]);
    float oy = asfloat(p.offsets[n * 2u * p.M + p.M + m]);
    ox = clamp(ox, -p.max_offset, p.max_offset);
    oy = clamp(oy, -p.max_offset, p.max_offset);

    float px = asfloat(p.xy[2u * n]) + ox;
    float py = asfloat(p.xy[2u * n + 1u]) + oy;

    p.out[2u * i] = 2.0f * px / (float(p.W) - 1.0f) - 1.0f;
    p.out[2u * i + 1u] = 2.0f * py / (float(p.H) - 1.0f) - 1.0f;
}

// ---------------------------------------------------------------------------
// LightGlue's assignment tail
//
// The last transformer layer produces a similarity matrix `sim[N0, N1]` and a
// matchability logit per keypoint. LightGlue turns those into an assignment
// with `sigmoid_log_double_softmax`:
//
//   score[i][j] = logsoftmax_j(sim)[i][j] + logsoftmax_i(sim)[i][j]
//                 + logsigmoid(z0[i]) + logsigmoid(z1[j])
//
// and then keeps a pair only if it is the argmax of BOTH its row and its
// column, with exp(score) over a threshold.
//
// Materializing `score` would be a second N0 x N1 matrix -- 16 MB at 2048
// keypoints each, written and read once. It is not needed: every term is
// separable into a row constant and a column constant, so the two log-sum-exps
// are vectors and the arg-max reads `sim` with those applied on the fly. What
// is left is three passes over the matrix instead of five, and no second
// matrix.
// ---------------------------------------------------------------------------

struct AssignParams {
    float* out_val;   // [R] the row's best score, before the row-constant terms
    uint*  out_idx;   // [R] its column
    uint*  sim;       // [N0, N1] f32
    uint*  bias;      // [C] added to every element of a row before the max
    uint   R;         // rows of the logical matrix this pass reduces over
    uint   C;         // columns
    uint   stride;    // element stride between logical rows of `sim`
    uint   step;      // element step between logical columns
    // What `sim` is multiplied by before the bias. TWO for the arg-max, not
    // one: the score is logsoftmax_rows + logsoftmax_cols, and each of those
    // carries its own copy of sim -- so the element being maximized is
    // 2*sim[i][j] - lse_row[i] - lse_col[j] + ls0[i] + ls1[j]. Dropping the
    // factor changes which column wins, not just the value, and leaves every
    // score far enough below zero that exp() is 0 and nothing matches.
    float  scale;
    uint   groups_per_row;
};

// `stride`/`step` are what make one kernel serve both directions: reducing
// over columns of sim is (stride = N1, step = 1), over rows is
// (stride = 1, step = N1). Reading a column is a strided gather either way --
// the matrix is too big to transpose for the sake of one pass over it.

groupshared float g_val[WG];
groupshared uint  g_idx[WG];

// logsumexp over each logical row, in one workgroup per row.
[shader("compute")]
[numthreads(WG, 1, 1)]
void assign_logsumexp(uint3 gid : SV_GroupID, uint3 tid : SV_GroupThreadID,
                      uniform AssignParams p) {
    uint row = gid.y * p.groups_per_row + gid.x;
    uint t = tid.x;

    // Two-pass max-then-sum, both as groupshared trees so wave32 and wave64
    // devices agree bit for bit.
    float m = -3.402823466e+38f;
    if (row < p.R)
        for (uint c = t; c < p.C; c += uint(WG))
            m = max(m, asfloat(p.sim[row * p.stride + c * p.step]));
    g_val[t] = m;
    GroupMemoryBarrierWithGroupSync();
    for (uint s = uint(WG) / 2u; s > 0u; s >>= 1) {
        if (t < s) g_val[t] = max(g_val[t], g_val[t + s]);
        GroupMemoryBarrierWithGroupSync();
    }
    float mx = g_val[0];
    GroupMemoryBarrierWithGroupSync();

    float sum = 0.0f;
    if (row < p.R)
        for (uint c = t; c < p.C; c += uint(WG))
            sum += exp(asfloat(p.sim[row * p.stride + c * p.step]) - mx);
    g_val[t] = sum;
    GroupMemoryBarrierWithGroupSync();
    for (uint s = uint(WG) / 2u; s > 0u; s >>= 1) {
        if (t < s) g_val[t] += g_val[t + s];
        GroupMemoryBarrierWithGroupSync();
    }
    if (row < p.R && t == 0u) p.out_val[row] = mx + log(g_val[0]);
}

// argmax over each logical row of (sim + bias), where `bias` carries every
// term that depends on the column alone.
[shader("compute")]
[numthreads(WG, 1, 1)]
void assign_argmax(uint3 gid : SV_GroupID, uint3 tid : SV_GroupThreadID,
                   uniform AssignParams p) {
    uint row = gid.y * p.groups_per_row + gid.x;
    uint t = tid.x;

    float best = -3.402823466e+38f;
    uint bi = 0u;
    if (row < p.R)
        for (uint c = t; c < p.C; c += uint(WG)) {
            float v = p.scale * asfloat(p.sim[row * p.stride + c * p.step]) +
                      asfloat(p.bias[c]);
            // Strictly greater, so the lowest index wins a tie -- the reduction
            // order below is not the column order, and a mutual-nearest test
            // that disagreed with itself run to run would make the whole
            // reconstruction irreproducible.
            if (v > best) { best = v; bi = c; }
        }
    g_val[t] = best;
    g_idx[t] = bi;
    GroupMemoryBarrierWithGroupSync();
    for (uint s = uint(WG) / 2u; s > 0u; s >>= 1) {
        if (t < s) {
            const bool take = g_val[t + s] > g_val[t] ||
                              (g_val[t + s] == g_val[t] && g_idx[t + s] < g_idx[t]);
            if (take) { g_val[t] = g_val[t + s]; g_idx[t] = g_idx[t + s]; }
        }
        GroupMemoryBarrierWithGroupSync();
    }
    if (row < p.R && t == 0u) {
        p.out_val[row] = g_val[0];
        p.out_idx[row] = g_idx[0];
    }
}
