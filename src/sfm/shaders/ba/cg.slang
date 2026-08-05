// Implicit-Schur preconditioned conjugate gradient: a matrix-free alternative
// to the dense Cholesky for large camera counts. The reduced system S dU = g
// is never formed; S x = B x - sum_p Acp_p W_p (Acp_p^T x) is applied through
// the per-observation A_cp blocks and per-point W, with B the block-diagonal
// damped camera Gram matrix -- "block diagonal" per *image*, which is exact:
// no observation couples two images' poses, and an image's intrinsics columns
// are inside its own block. When a group is shared the per-image blocks overlap
// on those columns, so their contributions are summed atomically instead of
// stored; u1 carries the exclusive/shared choice into the kernels that differ.
// Preconditioner: block Jacobi over u_prec_blocks (sfm/ba/Problem.h), the
// diagonal blocks of S under that partition.
//
// The whole CG loop is recorded into the command buffer: scalar state (rho,
// alpha, beta, tolerance, convergence flag) lives in u_cg_scal on the device
// and every loop kernel early-outs once the flag is set, so the host records
// a fixed iteration count and pays only dispatch overhead past convergence.
// The step accumulates into u_g in place (dp_accum/cam_update read it there,
// same as the dense path). Included at the end of ba.slang; shares its
// bindings, push constants and the g_sh shared pool.

[[vk::binding(29)]] StructuredBuffer<uint> u_cam_ranges;  // (num_images+1) obs grouped by image
[[vk::binding(30)]] StructuredBuffer<uint> u_cam_obs;     // obs indices per image
[[vk::binding(31)]] RWStructuredBuffer<Real> u_cg_r;      // (n_dim) residual
[[vk::binding(32)]] RWStructuredBuffer<Real> u_cg_z;      // (n_dim) M^-1 r
[[vk::binding(33)]] RWStructuredBuffer<Real> u_cg_p;      // (n_dim) search direction
[[vk::binding(34)]] RWStructuredBuffer<Real> u_cg_sp;     // (n_dim) S p
[[vk::binding(35)]] RWStructuredBuffer<Real> u_cg_v;      // (3 per point) W (Acp^T p)
[[vk::binding(36)]] RWStructuredBuffer<Real> u_cg_B;      // per-camera packed damped Jc^T Jc
[[vk::binding(37)]] RWStructuredBuffer<Real> u_cg_M;      // per-camera packed S_cc -> L factor
[[vk::binding(38)]] RWStructuredBuffer<Real> u_cg_scal;   // scalar state (kCg* indices)
[[vk::binding(39)]] RWStructuredBuffer<Real> u_cg_part;   // dot-product partials, 2 slots
[[vk::binding(40)]] RWStructuredBuffer<AWord> u_cg_sp_a;  // atomic views (chunked accumulation)
[[vk::binding(41)]] RWStructuredBuffer<AWord> u_cg_B_a;
[[vk::binding(42)]] RWStructuredBuffer<AWord> u_cg_M_a;
[[vk::binding(43)]] StructuredBuffer<uint> u_cam_chunks;  // 3 per chunk: image, start, count
[[vk::binding(44)]] StructuredBuffer<uint4> u_prec_blocks;  // (col0, len0, col1, len1) per block

static const uint kCgRho = 0;    // r^T z
static const uint kCgPAp = 1;    // p^T S p
static const uint kCgAlpha = 2;
static const uint kCgBeta = 3;
static const uint kCgRR = 4;     // ||r||^2
static const uint kCgTol2 = 5;   // eta^2 ||r0||^2
static const uint kCgConv = 6;   // > 0.5: loop kernels become no-ops
static const uint kCgIters = 7;  // completed CG iterations

// per-camera packed block stride (lower triangle of a kMaxCamDof block)
static const uint kCamBlk = kMaxCamDof * (kMaxCamDof + 1) / 2;

bool cgDone() { return u_cg_scal[kCgConv] > Real(0.5); }

// global columns of camera `img` (pose then intrinsics)
void camCols(uint img, out uint dof, out uint cols[kMaxCamDof]) {
    uint4 gi = u_group_info[u_image_group[img]];
    dof = 6 + gi.z;
    for (int i = 0; i < 6; i++) cols[i] = 6 * img + i;
    for (uint i = 0; i < gi.z; i++) cols[6 + i] = gi.y + i;
}

// global columns of preconditioner block `b` (two contiguous ranges)
void precCols(uint b, out uint dof, out uint cols[kMaxCamDof]) {
    uint4 blk = u_prec_blocks[b];
    dof = blk.y + blk.w;
    for (uint i = 0; i < blk.y; i++) cols[i] = blk.x + i;
    for (uint i = 0; i < blk.w; i++) cols[blk.y + i] = blk.z + i;
}

// ---------------------------------------------------------------------
// per-lambda assembly: B_c, M_c = S_cc, and the reduced RHS
// ---------------------------------------------------------------------

// One warp per chunk of one camera's observation list (u_cam_obs CSR split by
// u_cam_chunks) produces the chunk's share of
//   B_c = Jc^T Jc with (1+lambda) on the diagonal          (packed, u_cg_B)
//   M_c = B_c - sum_obs Acp W Acp^T = S_cc                 (packed, u_cg_M)
//   g_c = Jc^T r - sum_obs Acp W Bp                        (u_g)
// accumulated atomically into the zero-filled outputs (multiplicative damping
// is linear in the sum, so each chunk damps its own partial). This is exactly
// the diagonal-pair part of schur_pair, with the 2x2 middle matrix
// Q = Jp W Jp^T and its contraction Z = Q^T Jc built in shared memory per
// observation (the CG path stores no Y).
//
// M lands in the preconditioner block(s) covering these columns. With exclusive
// groups that is one block per image and the whole dof x dof block goes there;
// with a shared group the pose half goes to the image's block, the intrinsics
// half to the group's, and the cross terms are dropped (they belong to no block
// of the partition).
// u0 = num chunks, u1 = exclusive groups, u2 = group-block base, f0 = lambda
static const uint kShcJ = 0;                    // Jc    2 x dof
static const uint kShcZ = 2 * kMaxCamDof;       // Z = Q^T Jc   2 x dof
static const uint kShcQ = 4 * kMaxCamDof;       // Q = Jp W Jp^T   2 x 2
static const uint kShcY = 4 * kMaxCamDof + 4;   // Y = (Jp W)^T    2 x 3
static const uint kShcR = 4 * kMaxCamDof + 10;  // res   2
static const uint kShcB = 4 * kMaxCamDof + 12;  // Bp    3
static const uint kShcP = 4 * kMaxCamDof + 15;  // Jp    2 x 3

[shader("compute")] [numthreads(kPairWg, 1, 1)]
void cg_cam_diag(uint3 wg : SV_GroupID, uint gi : SV_GroupIndex) {
    if (wg.x >= u0) return;
    uint img = u_cam_chunks[3 * wg.x];
    uint t0 = u_cam_chunks[3 * wg.x + 1];
    uint t1 = t0 + u_cam_chunks[3 * wg.x + 2];
    uint dof;
    uint cols[kMaxCamDof];
    camCols(img, dof, cols);
    uint nel = dof * dof;

    Real accB[kPairEl], accM[kPairEl];
    uint elr[kPairEl], elc[kPairEl];
    [ForceUnroll]
    for (uint k = 0; k < kPairEl; k++) {
        accB[k] = Real(0.0);
        accM[k] = Real(0.0);
        uint el = gi + kPairWg * k;
        elr[k] = el / dof;
        elc[k] = el - elr[k] * dof;
    }
    Real gacc = Real(0.0);

    for (uint t = t0; t < t1; t++) {
        uint o = u_cam_obs[t];
        uint pt = u_obs_point[o];
        uint jb = u_jc_off[o];
        GroupMemoryBarrierWithGroupSync();  // previous obs done being read
        for (uint i = gi; i < 2 * dof; i += kPairWg) g_sh[kShcJ + i] = u_Jc[jb + i];
        if (gi < 6) g_sh[kShcP + gi] = u_Jp[6 * o + gi];
        if (gi < 2) g_sh[kShcR + gi] = u_res[2 * o + gi];
        if (gi < 3) g_sh[kShcB + gi] = u_Bp[3 * pt + gi];
        // Y = (Jp W)^T, then Q = Jp W Jp^T -- one 2x2 middle matrix instead of
        // the dof x 3 block T = Acp W the materialized layout carried per obs
        if (gi < 6) {
            uint j = gi >> 1, row = gi & 1;
            g_sh[kShcY + gi] = u_Jp[6 * o + row] * u_W[9 * pt + j]
                             + u_Jp[6 * o + 2 + row] * u_W[9 * pt + 3 + j]
                             + u_Jp[6 * o + 4 + row] * u_W[9 * pt + 6 + j];
        }
        GroupMemoryBarrierWithGroupSync();
        if (gi < 4) {
            uint qr = gi >> 1, qs = gi & 1;
            g_sh[kShcQ + gi] = g_sh[kShcY + qr] * g_sh[kShcP + qs]
                             + g_sh[kShcY + 2 + qr] * g_sh[kShcP + 2 + qs]
                             + g_sh[kShcY + 4 + qr] * g_sh[kShcP + 4 + qs];
        }
        GroupMemoryBarrierWithGroupSync();
        for (uint i = gi; i < 2 * dof; i += kPairWg) {
            uint a = i >> 1, s = i & 1;
            g_sh[kShcZ + i] = g_sh[kShcJ + 2 * a] * g_sh[kShcQ + s]
                            + g_sh[kShcJ + 2 * a + 1] * g_sh[kShcQ + 2 + s];
        }
        GroupMemoryBarrierWithGroupSync();
        [ForceUnroll]
        for (uint k = 0; k < kPairEl; k++) {
            uint el = gi + kPairWg * k;
            if (el >= nel) continue;
            uint r = elr[k], c = elc[k];
            if (r < c) continue;  // symmetric block: lower triangle only
            accM[k] = accM[k] + g_sh[kShcZ + 2 * r] * g_sh[kShcJ + 2 * c]
                              + g_sh[kShcZ + 2 * r + 1] * g_sh[kShcJ + 2 * c + 1];
            accB[k] = accB[k] + g_sh[kShcJ + 2 * r] * g_sh[kShcJ + 2 * c]
                              + g_sh[kShcJ + 2 * r + 1] * g_sh[kShcJ + 2 * c + 1];
        }
        if (gi < dof) {
            Real yb0 = g_sh[kShcY] * g_sh[kShcB] + g_sh[kShcY + 2] * g_sh[kShcB + 1]
                     + g_sh[kShcY + 4] * g_sh[kShcB + 2];
            Real yb1 = g_sh[kShcY + 1] * g_sh[kShcB] + g_sh[kShcY + 3] * g_sh[kShcB + 1]
                     + g_sh[kShcY + 5] * g_sh[kShcB + 2];
            gacc = gacc + g_sh[kShcJ + 2 * gi] * (yb0 - g_sh[kShcR])
                        + g_sh[kShcJ + 2 * gi + 1] * (yb1 - g_sh[kShcR + 1]);
        }
    }

    uint cb = kCamBlk * img;
    uint gb = kCamBlk * (u2 + u_image_group[img]);
    [ForceUnroll]
    for (uint k = 0; k < kPairEl; k++) {
        uint el = gi + kPairWg * k;
        if (el >= nel) continue;
        uint r = elr[k], c = elc[k];
        if (r < c) continue;
        Real jj = accB[k];
        if (r == c) jj = jj * (Real(1.0) + Real(f0));  // multiplicative damping
        Real m = jj - accM[k];
        if (isfinite(jj)) atomicAccum(u_cg_B_a, cb + pidx(r, c), jj);
        if (!isfinite(m)) continue;
        if (u1 != 0) atomicAccum(u_cg_M_a, cb + pidx(r, c), m);
        else if (r < 6) atomicAccum(u_cg_M_a, cb + pidx(r, c), m);
        else if (c >= 6) atomicAccum(u_cg_M_a, gb + pidx(r - 6, c - 6), m);
    }
    if (gi < dof && isfinite(gacc))
        atomicAccum(u_g_a, cols[gi], -gacc);
}

// Factor each block in place (packed lower Cholesky, guarded like chol_diag).
// One thread per block; blocks are at most kMaxCamDof wide. u0 = num blocks
[shader("compute")] [numthreads(256, 1, 1)]
void cg_prec_fact(uint3 t : SV_DispatchThreadID) {
    uint b = t.x;
    if (b >= u0) return;
    uint4 blk = u_prec_blocks[b];
    uint dof = blk.y + blk.w;
    if (dof == 0) return;
    uint cb = kCamBlk * b;
    Real L[kCamBlk];
    for (uint i = 0; i < dof * (dof + 1) / 2; i++) L[i] = u_cg_M[cb + i];
    for (uint j = 0; j < dof; j++) {
        Real d = sqrt(max(L[pidx(j, j)], Real(1e-30)));
        L[pidx(j, j)] = d;
        for (uint i = j + 1; i < dof; i++) L[pidx(i, j)] = L[pidx(i, j)] / d;
        for (uint c = j + 1; c < dof; c++)
            for (uint i = c; i < dof; i++)
                L[pidx(i, c)] = L[pidx(i, c)] - L[pidx(i, j)] * L[pidx(c, j)];
    }
    for (uint i = 0; i < dof * (dof + 1) / 2; i++) u_cg_M[cb + i] = L[i];
}

// z = M^-1 r, one thread per block (fwd + bwd solve of the small factor).
// u0 = num blocks, u1 = 1 to respect the convergence flag
[shader("compute")] [numthreads(256, 1, 1)]
void cg_prec_apply(uint3 t : SV_DispatchThreadID) {
    if (u1 != 0 && cgDone()) return;
    uint b = t.x;
    if (b >= u0) return;
    uint dof;
    uint cols[kMaxCamDof];
    precCols(b, dof, cols);
    if (dof == 0) return;
    uint cb = kCamBlk * b;
    // Entries below i are always written before they are read; zeroed so
    // slangc can see it.
    Real y[kMaxCamDof];
    for (uint i = 0; i < kMaxCamDof; i++) y[i] = Real(0);
    for (uint i = 0; i < dof; i++) {
        Real v = u_cg_r[cols[i]];
        for (uint j = 0; j < i; j++) v = v - u_cg_M[cb + pidx(i, j)] * y[j];
        y[i] = v / u_cg_M[cb + pidx(i, i)];
    }
    for (int i = int(dof) - 1; i >= 0; i--) {
        uint ui = uint(i);
        Real v = y[ui];
        for (uint j = ui + 1; j < dof; j++) v = v - u_cg_M[cb + pidx(j, ui)] * y[j];
        y[ui] = v / u_cg_M[cb + pidx(ui, ui)];
    }
    for (uint i = 0; i < dof; i++) u_cg_z[cols[i]] = y[i];
}

// ---------------------------------------------------------------------
// CG loop kernels
// ---------------------------------------------------------------------

// r = g (the reduced RHS assembled by cg_cam_diag), x (in u_g) = 0, scalars
// reset. u0 = n_dim
[shader("compute")] [numthreads(256, 1, 1)]
void cg_init(uint3 t : SV_DispatchThreadID) {
    uint i = t.x;
    if (i < 8) u_cg_scal[i] = Real(0.0);
    if (i >= u0) return;
    u_cg_r[i] = u_g[i];
    u_g[i] = Real(0.0);
}

// p = z. u0 = n_dim
[shader("compute")] [numthreads(256, 1, 1)]
void cg_copy(uint3 t : SV_DispatchThreadID) {
    uint i = t.x;
    if (i >= u0) return;
    u_cg_p[i] = u_cg_z[i];
}

// v_p = W_p (Acp_p^T p): one thread per point looping its track (linear in
// track length, like schur_obs). u0 = num_points
[shader("compute")] [numthreads(256, 1, 1)]
void cg_gather(uint3 t : SV_DispatchThreadID) {
    if (cgDone()) return;
    uint p = t.x;
    if (p >= u0) return;
    Real u0v = Real(0.0), u1v = Real(0.0), u2v = Real(0.0);
    uint t0 = u_obs_ranges[p], t1 = u_obs_ranges[p + 1];
    for (uint oi = t0; oi < t1; oi++) {
        uint dof, base;
        uint cols[kMaxCamDof];
        obsCamCols(oi, dof, base, cols);
        // Acp^T p = Jp^T (Jc p): the contraction runs over the 2-vector
        // residual, so only two accumulators cross the dof loop
        Real d0 = Real(0.0), d1 = Real(0.0);
        for (uint a = 0; a < dof; a++) {
            Real xa = u_cg_p[cols[a]];
            d0 = d0 + u_Jc[base + 2 * a] * xa;
            d1 = d1 + u_Jc[base + 2 * a + 1] * xa;
        }
        u0v = u0v + u_Jp[6 * oi] * d0 + u_Jp[6 * oi + 1] * d1;
        u1v = u1v + u_Jp[6 * oi + 2] * d0 + u_Jp[6 * oi + 3] * d1;
        u2v = u2v + u_Jp[6 * oi + 4] * d0 + u_Jp[6 * oi + 5] * d1;
    }
    for (int i = 0; i < 3; i++)
        u_cg_v[3 * p + i] = u_W[9 * p + 3 * i] * u0v + u_W[9 * p + 3 * i + 1] * u1v
                          + u_W[9 * p + 3 * i + 2] * u2v;
}

// Sp_c = B_c p_c: one warp per camera initializes u_cg_sp with the small
// symmetric block multiply; cg_scatter then subtracts the cross term
// atomically. The pose rows belong to this camera alone and are stored; the
// intrinsics rows may be shared with other images, in which case the host has
// pre-zeroed that slice of u_cg_sp and every image adds into it.
// u0 = num_images, u1 = exclusive groups
[shader("compute")] [numthreads(kPairWg, 1, 1)]
void cg_bmul(uint3 wg : SV_GroupID, uint gi : SV_GroupIndex) {
    if (cgDone()) return;
    uint img = wg.x;
    if (img >= u0) return;
    uint dof;
    uint cols[kMaxCamDof];
    camCols(img, dof, cols);
    if (gi >= dof) return;
    uint cb = kCamBlk * img;
    Real acc = Real(0.0);
    for (uint b = 0; b < dof; b++) {
        uint idx = gi >= b ? pidx(gi, b) : pidx(b, gi);
        acc = acc + u_cg_B[cb + idx] * u_cg_p[cols[b]];
    }
    if (u1 != 0 || gi < 6) storeAccum(u_cg_sp_a, cols[gi], acc);
    else if (isfinite(acc)) atomicAccum(u_cg_sp_a, cols[gi], acc);
}

// Sp_c -= sum_{obs of c} Acp v_{pt(obs)} = sum Jc^T (Jp v): one warp per chunk;
// lane i owns Jc element i = 2a+row of every observation in turn, so the block
// loads are fully coalesced, and the per-lane partials reduce as pairs into dof
// atomic adds per chunk. The 2-vector Jp v is small enough that every lane
// recomputes it rather than paying a barrier for it. u0 = num chunks
[shader("compute")] [numthreads(kPairWg, 1, 1)]
void cg_scatter(uint3 wg : SV_GroupID, uint gi : SV_GroupIndex) {
    if (cgDone()) return;
    if (wg.x >= u0) return;
    uint img = u_cam_chunks[3 * wg.x];
    uint t0 = u_cam_chunks[3 * wg.x + 1];
    uint t1 = t0 + u_cam_chunks[3 * wg.x + 2];
    uint dof;
    uint cols[kMaxCamDof];
    camCols(img, dof, cols);

    // 2*dof can reach 36 (> kPairWg), so each lane owns two elements
    uint gi2 = gi + kPairWg;
    uint row = gi & 1, row2 = gi2 & 1;
    Real pacc = Real(0.0), pacc2 = Real(0.0);
    bool live = gi < 2 * dof, live2 = gi2 < 2 * dof;
    for (uint t = t0; t < t1; t++) {
        uint o = u_cam_obs[t];
        uint pt = u_obs_point[o];
        uint base = u_jc_off[o];
        Real d0 = u_Jp[6 * o] * u_cg_v[3 * pt] + u_Jp[6 * o + 2] * u_cg_v[3 * pt + 1]
                + u_Jp[6 * o + 4] * u_cg_v[3 * pt + 2];
        Real d1 = u_Jp[6 * o + 1] * u_cg_v[3 * pt] + u_Jp[6 * o + 3] * u_cg_v[3 * pt + 1]
                + u_Jp[6 * o + 5] * u_cg_v[3 * pt + 2];
        if (live)
            pacc = pacc + u_Jc[base + gi] * (row == 0 ? d0 : d1);
        if (live2)
            pacc2 = pacc2 + u_Jc[base + gi2] * (row2 == 0 ? d0 : d1);
    }
    g_sh[gi] = pacc;
    g_sh[gi2] = pacc2;
    GroupMemoryBarrierWithGroupSync();
    if (gi < dof) {
        Real v = g_sh[2 * gi] + g_sh[2 * gi + 1];
        if (isfinite(v))
            atomicAccum(u_cg_sp_a, cols[gi], -v);
    }
}

// Dot-product partials over the n_dim vectors, one value per workgroup.
// u1 selects the products (0: p.Sp -> slot 0; 1: r.z -> slot 0, r.r -> slot 1)
// u0 = n_dim, u2 = partial stride, u3 = 1 to respect the convergence flag
[shader("compute")] [numthreads(256, 1, 1)]
void cg_red2(uint3 t : SV_DispatchThreadID, uint3 wg : SV_GroupID, uint gi : SV_GroupIndex) {
    if (u3 != 0 && cgDone()) return;
    uint i = t.x;
    Real a = Real(0.0), b = Real(0.0);
    if (i < u0) {
        if (u1 == 0) {
            a = u_cg_p[i] * u_cg_sp[i];
        } else {
            a = u_cg_r[i] * u_cg_z[i];
            b = u_cg_r[i] * u_cg_r[i];
        }
    }
    g_sh[gi] = a;
    g_sh[256 + gi] = b;
    GroupMemoryBarrierWithGroupSync();
    [ForceUnroll]
    for (uint s = 128; s > 0; s >>= 1) {
        if (gi < s) {
            g_sh[gi] = g_sh[gi] + g_sh[gi + s];
            g_sh[256 + gi] = g_sh[256 + gi] + g_sh[256 + gi + s];
        }
        GroupMemoryBarrierWithGroupSync();
    }
    if (gi == 0) {
        u_cg_part[wg.x] = g_sh[0];
        if (u1 != 0) u_cg_part[u2 + wg.x] = g_sh[256];
    }
}

// Final partial reduction + scalar updates; one workgroup.
// u0 = partial count, u1 = mode (0: init rho/tol, 1: alpha, 2: beta/converge),
// u2 = partial stride, f0 = eta (mode 0: relative residual tolerance)
[shader("compute")] [numthreads(256, 1, 1)]
void cg_fin(uint gi : SV_GroupIndex) {
    if (u1 != 0 && cgDone()) return;
    Real a = Real(0.0), b = Real(0.0);
    for (uint i = gi; i < u0; i += 256) {
        a = a + u_cg_part[i];
        if (u1 != 1) b = b + u_cg_part[u2 + i];
    }
    g_sh[gi] = a;
    g_sh[256 + gi] = b;
    GroupMemoryBarrierWithGroupSync();
    [ForceUnroll]
    for (uint s = 128; s > 0; s >>= 1) {
        if (gi < s) {
            g_sh[gi] = g_sh[gi] + g_sh[gi + s];
            g_sh[256 + gi] = g_sh[256 + gi] + g_sh[256 + gi + s];
        }
        GroupMemoryBarrierWithGroupSync();
    }
    if (gi != 0) return;
    Real A = g_sh[0], B = g_sh[256];
    if (u1 == 0) {  // A = r.z, B = ||r||^2 of the initial residual
        u_cg_scal[kCgRho] = A;
        u_cg_scal[kCgRR] = B;
        u_cg_scal[kCgTol2] = Real(f0) * Real(f0) * B;
        u_cg_scal[kCgIters] = Real(0.0);
        // zero RHS (or a broken preconditioner) => keep the zero step
        u_cg_scal[kCgConv] = (isfinite(A) && A > Real(0.0)) ? Real(0.0) : Real(1.0);
    } else if (u1 == 1) {  // A = p.Sp
        Real rho = u_cg_scal[kCgRho];
        bool ok = isfinite(A) && A > Real(0.0);
        u_cg_scal[kCgPAp] = A;
        u_cg_scal[kCgAlpha] = ok ? rho / A : Real(0.0);
        if (!ok) u_cg_scal[kCgConv] = Real(1.0);  // lost positive-definiteness
    } else {  // A = new r.z, B = new ||r||^2
        Real rho = u_cg_scal[kCgRho];
        u_cg_scal[kCgBeta] = A / rho;
        u_cg_scal[kCgRho] = A;
        u_cg_scal[kCgRR] = B;
        u_cg_scal[kCgIters] = u_cg_scal[kCgIters] + Real(1.0);
        if (!(B > u_cg_scal[kCgTol2]) || !isfinite(A) || !(A > Real(0.0)))
            u_cg_scal[kCgConv] = Real(1.0);
    }
}

// x += alpha p (x accumulates in u_g), r -= alpha S p. u0 = n_dim
[shader("compute")] [numthreads(256, 1, 1)]
void cg_axpy(uint3 t : SV_DispatchThreadID) {
    if (cgDone()) return;
    uint i = t.x;
    if (i >= u0) return;
    Real alpha = u_cg_scal[kCgAlpha];
    u_g[i] = u_g[i] + alpha * u_cg_p[i];
    u_cg_r[i] = u_cg_r[i] - alpha * u_cg_sp[i];
}

// p = z + beta p. u0 = n_dim
[shader("compute")] [numthreads(256, 1, 1)]
void cg_updp(uint3 t : SV_DispatchThreadID) {
    if (cgDone()) return;
    uint i = t.x;
    if (i >= u0) return;
    u_cg_p[i] = u_cg_z[i] + u_cg_scal[kCgBeta] * u_cg_p[i];
}
