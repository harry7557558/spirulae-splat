// Implicit-Schur preconditioned conjugate gradient: a matrix-free alternative
// to the dense Cholesky for large camera counts. The reduced system S dU = g
// is never formed; S x = B x - sum_p Acp_p W_p (Acp_p^T x) is applied through
// the per-observation A_cp blocks and per-point W, with B the block-diagonal
// damped camera Gram matrix. Requires per-image intrinsics groups (same
// exclusivity condition as the pair-Schur kernel; checked host-side).
// Preconditioner: block Jacobi with the exact diagonal blocks S_cc.
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

// global columns of camera `img` (pose then intrinsics; exclusive groups)
void camCols(uint img, out uint dof, out uint cols[kMaxCamDof]) {
    uint4 gi = u_group_info[u_image_group[img]];
    dof = 6 + gi.z;
    for (int i = 0; i < 6; i++) cols[i] = 6 * img + i;
    for (uint i = 0; i < gi.z; i++) cols[6 + i] = gi.y + i;
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
// the diagonal-pair part of schur_pair, with T = Acp W recomputed in shared
// memory per observation (the CG path stores no T).
// u0 = num chunks, f0 = lambda
static const uint kShcA = 0;                    // Acp   3 x dof
static const uint kShcT = 3 * kMaxCamDof;       // T     3 x dof
static const uint kShcW = 6 * kMaxCamDof;       // W     9
static const uint kShcJ = 6 * kMaxCamDof + 9;   // Jc    2 x dof
static const uint kShcR = 8 * kMaxCamDof + 9;   // res   2
static const uint kShcB = 8 * kMaxCamDof + 11;  // Bp    3

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
        uint base = u_acp_off[o];
        uint jb = base / 3 * 2;
        GroupMemoryBarrierWithGroupSync();  // previous obs done being read
        for (uint i = gi; i < 3 * dof; i += kPairWg) g_sh[kShcA + i] = u_acp[base + i];
        if (gi < 9) g_sh[kShcW + gi] = u_W[9 * pt + gi];
        for (uint i = gi; i < 2 * dof; i += kPairWg) g_sh[kShcJ + i] = u_Jc[jb + i];
        if (gi < 2) g_sh[kShcR + gi] = u_res[2 * o + gi];
        if (gi < 3) g_sh[kShcB + gi] = u_Bp[3 * pt + gi];
        GroupMemoryBarrierWithGroupSync();
        for (uint i = gi; i < 3 * dof; i += kPairWg) {
            uint r = i / 3, k = i - 3 * r;
            g_sh[kShcT + i] = g_sh[kShcA + 3 * r] * g_sh[kShcW + k]
                            + g_sh[kShcA + 3 * r + 1] * g_sh[kShcW + 3 + k]
                            + g_sh[kShcA + 3 * r + 2] * g_sh[kShcW + 6 + k];
        }
        GroupMemoryBarrierWithGroupSync();
        [ForceUnroll]
        for (uint k = 0; k < kPairEl; k++) {
            uint el = gi + kPairWg * k;
            if (el >= nel) continue;
            uint r = elr[k], c = elc[k];
            if (r < c) continue;  // symmetric block: lower triangle only
            accM[k] = accM[k] + g_sh[kShcT + 3 * r] * g_sh[kShcA + 3 * c]
                              + g_sh[kShcT + 3 * r + 1] * g_sh[kShcA + 3 * c + 1]
                              + g_sh[kShcT + 3 * r + 2] * g_sh[kShcA + 3 * c + 2];
            accB[k] = accB[k] + g_sh[kShcJ + 2 * r] * g_sh[kShcJ + 2 * c]
                              + g_sh[kShcJ + 2 * r + 1] * g_sh[kShcJ + 2 * c + 1];
        }
        if (gi < dof) {
            gacc = gacc + g_sh[kShcT + 3 * gi] * g_sh[kShcB]
                        + g_sh[kShcT + 3 * gi + 1] * g_sh[kShcB + 1]
                        + g_sh[kShcT + 3 * gi + 2] * g_sh[kShcB + 2]
                        - g_sh[kShcJ + 2 * gi] * g_sh[kShcR]
                        - g_sh[kShcJ + 2 * gi + 1] * g_sh[kShcR + 1];
        }
    }

    uint cb = kCamBlk * img;
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
        if (isfinite(m)) atomicAccum(u_cg_M_a, cb + pidx(r, c), m);
    }
    if (gi < dof && isfinite(gacc))
        atomicAccum(u_g_a, cols[gi], -gacc);
}

// Factor each M_c in place (packed lower Cholesky, guarded like chol_diag).
// One thread per camera; blocks are at most kMaxCamDof wide. u0 = num_images
[shader("compute")] [numthreads(256, 1, 1)]
void cg_prec_fact(uint3 t : SV_DispatchThreadID) {
    uint img = t.x;
    if (img >= u0) return;
    uint dof = 6 + u_group_info[u_image_group[img]].z;
    uint cb = kCamBlk * img;
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

// z = M^-1 r, one thread per camera (fwd + bwd solve of the small factor).
// u0 = num_images, u1 = 1 to respect the convergence flag
[shader("compute")] [numthreads(256, 1, 1)]
void cg_prec_apply(uint3 t : SV_DispatchThreadID) {
    if (u1 != 0 && cgDone()) return;
    uint img = t.x;
    if (img >= u0) return;
    uint dof;
    uint cols[kMaxCamDof];
    camCols(img, dof, cols);
    uint cb = kCamBlk * img;
    Real y[kMaxCamDof];
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
        for (uint a = 0; a < dof; a++) {
            Real xa = u_cg_p[cols[a]];
            u0v = u0v + u_acp[base + 3 * a] * xa;
            u1v = u1v + u_acp[base + 3 * a + 1] * xa;
            u2v = u2v + u_acp[base + 3 * a + 2] * xa;
        }
    }
    for (int i = 0; i < 3; i++)
        u_cg_v[3 * p + i] = u_W[9 * p + 3 * i] * u0v + u_W[9 * p + 3 * i + 1] * u1v
                          + u_W[9 * p + 3 * i + 2] * u2v;
}

// Sp_c = B_c p_c: one warp per camera initializes u_cg_sp with the small
// symmetric block multiply; cg_scatter then subtracts the cross term
// atomically. u0 = num_images
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
    storeAccum(u_cg_sp_a, cols[gi], acc);
}

// Sp_c -= sum_{obs of c} Acp v_{pt(obs)}: one warp per chunk; lane i owns Acp
// element i = 3a+k of every observation in turn, so the block loads are fully
// coalesced, and the per-lane partials reduce as triples into dof atomic adds
// per chunk. u0 = num chunks
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

    // 3*dof can reach 36 (> kPairWg), so each lane owns two elements
    uint gi2 = gi + kPairWg;
    uint k3 = gi - 3 * (gi / 3);
    uint k3b = gi2 - 3 * (gi2 / 3);
    Real pacc = Real(0.0), pacc2 = Real(0.0);
    bool live = gi < 3 * dof, live2 = gi2 < 3 * dof;
    for (uint t = t0; t < t1; t++) {
        uint o = u_cam_obs[t];
        uint pt = u_obs_point[o];
        uint base = u_acp_off[o];
        if (live)
            pacc = pacc + u_acp[base + gi] * u_cg_v[3 * pt + k3];
        if (live2)
            pacc2 = pacc2 + u_acp[base + gi2] * u_cg_v[3 * pt + k3b];
    }
    g_sh[gi] = pacc;
    g_sh[gi2] = pacc2;
    GroupMemoryBarrierWithGroupSync();
    if (gi < dof) {
        Real v = g_sh[3 * gi] + g_sh[3 * gi + 1] + g_sh[3 * gi + 2];
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
