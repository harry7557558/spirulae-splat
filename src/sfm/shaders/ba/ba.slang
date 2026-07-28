// Bundle adjustment kernels.
//
// Parameterization: each image has a 6-DOF pose; intrinsics live in "groups"
// so any number of images can share one set of intrinsics, and each group can
// use a different camera model. The camera-side unknown vector u (dimension
// n_dim) is laid out as [pose_0(6), pose_1(6), ..., all group intrinsics].
// Points are eliminated with the Schur complement; the reduced system
//   S dU = g,  S = J_u^T J_u - A_up (App + lambda)^-1 A_up^T
// is assembled directly into a packed lower-triangular dense matrix u_S and
// solved in place by the blocked Cholesky in cholesky.slang.
//
// Push constant slots (meaning per kernel, see host sfm/ba/Solver.h):
//   u0..u3  - counts / offsets / block index
//   f0..f3  - damping and loss parameters

#include "common/camera.slang"
#include "common/loss.slang"

// Bounds on 6 + kNumIntr for the model-agnostic Schur kernels, which size
// register arrays by a compile-time dof. These kernels are *not* per-model
// (a point's observations may span groups with different models), so they take
// the dof as a generic value parameter MD and are instantiated at two tiers:
//   compact -- covers every model up to RADIAL (dof 11); byte-identical to the
//              pre-OpenCV kernels, so problems without OpenCV pay nothing.
//   wide    -- covers OpenCV (dof 14).
// The host (sfm/ba/Solver.h) picks the tier from the problem's actual maximum dof.
// Adding a model with dof > kDofWide means bumping kDofWide (and kMaxCamDof in
// sfm/ba/Problem.h). kDofCompact must stay >= the largest non-OpenCV model's dof.
// Three dof tiers so a problem pays for only the widest camera it actually uses:
//   compact 12 -- up to RADIAL (dof 11); byte-identical to the pre-OpenCV kernels
//   mid     14 -- OpenCV / OpenCV_FISHEYE (dof 14)
//   wide    18 -- FullOpenCV (15) / ThinPrismFisheye (18)
static const int kDofCompact = 12;
static const int kDofMid = 14;
static const int kDofWide = 18;
// The implicit-Schur CG kernels (cg.slang) are not tiered -- they run only for
// very large problems (n_dim > 8192) and keep the full width.
static const int kMaxCamDof = kDofWide;

[[vk::push_constant]] cbuffer PC {
    uint u0; uint u1; uint u2; uint u3;
    float f0; float f1; float f2; float f3;
};

// ---- problem data ----
[[vk::binding(0)]]  StructuredBuffer<Real> u_obs;          // (2 per obs)
[[vk::binding(1)]]  StructuredBuffer<uint> u_obs_image;    // (per obs)
[[vk::binding(2)]]  StructuredBuffer<uint> u_obs_point;    // (per obs)
[[vk::binding(3)]]  StructuredBuffer<uint> u_image_group;  // (per image)
[[vk::binding(4)]]  StructuredBuffer<uint4> u_group_info;  // (per group) x=intr param offset, y=intr column base, z=n_intr, w=model id
[[vk::binding(5)]]  RWStructuredBuffer<Real> u_poses;      // (6 per image)
[[vk::binding(6)]]  RWStructuredBuffer<Real> u_intr;       // (flat intrinsics)
[[vk::binding(7)]]  RWStructuredBuffer<Real> u_points;     // (3 per point)
[[vk::binding(8)]]  StructuredBuffer<uint> u_obs_ranges;   // (num_points+1) obs sorted by point
[[vk::binding(9)]]  StructuredBuffer<uint> u_model_obs;    // per-model obs index lists, concatenated
[[vk::binding(10)]] StructuredBuffer<uint> u_acp_off;      // (per obs) element offset into u_acp

// ---- work buffers (double-bound: Real view + atomic-word view) ----
[[vk::binding(11)]] RWStructuredBuffer<Real> u_acp;        // A_cp blocks, dof x 3 per obs
[[vk::binding(12)]] RWStructuredBuffer<Real> u_S;          // packed lower triangle, n_dim*(n_dim+1)/2
[[vk::binding(13)]] RWStructuredBuffer<AWord> u_S_a;
[[vk::binding(14)]] RWStructuredBuffer<Real> u_g;          // gradient -> in-place solution dU
[[vk::binding(15)]] RWStructuredBuffer<AWord> u_g_a;
[[vk::binding(16)]] RWStructuredBuffer<Real> u_App;        // (9 per point)
[[vk::binding(17)]] RWStructuredBuffer<AWord> u_App_a;
[[vk::binding(18)]] RWStructuredBuffer<Real> u_Bp;         // (3 per point)
[[vk::binding(19)]] RWStructuredBuffer<AWord> u_Bp_a;
[[vk::binding(20)]] RWStructuredBuffer<Real> u_cost;       // (1)
[[vk::binding(21)]] RWStructuredBuffer<AWord> u_cost_a;

// ---- pair-aggregated Schur tables (see sfm/ba/Problem.h buildPairTables) ----
[[vk::binding(22)]] StructuredBuffer<uint2> u_pair_entries; // (obs_i, obs_j) per entry
[[vk::binding(23)]] StructuredBuffer<uint2> u_pair_chunks;  // (offset, count|atomic_flag)
[[vk::binding(24)]] RWStructuredBuffer<Real> u_W;           // (9 per point) (App+l)^-1
[[vk::binding(25)]] RWStructuredBuffer<Real> u_T;           // T = Acp W per obs (Acp layout)
[[vk::binding(26)]] RWStructuredBuffer<Real> u_y;           // (n_dim) triangular-solve scratch
[[vk::binding(27)]] RWStructuredBuffer<Real> u_Jc;          // camera Jacobian rows per obs,
                                                            // [2*a+row], offset acp_off*2/3
[[vk::binding(28)]] RWStructuredBuffer<Real> u_res;         // weighted residual (2 per obs)

// packed lower-triangular index, row >= col
uint pidx(uint r, uint c) { return r * (r + 1) / 2 + c; }

// =====================================================================
// cost
// =====================================================================

static const uint kCostWg = 256;

// One groupshared pool reused by every kernel (the sum of all Workgroup
// variables in a module counts against maxComputeSharedMemorySize on every
// pipeline, so distinct arrays per kernel would blow the limit). The second
// Cholesky tile lives at kSh2 with row pitch 33 to avoid bank conflicts.
static const uint kSh2 = 32 * 32;
static const uint kShPool = 32 * 32 + 32 * 33 + 32;
groupshared Real g_sh[kShPool];

// u0 = obs count for this model, u1 = offset into u_model_obs, f0 = loss param
void cost_impl<M : ICameraModel>(uint tid, uint lid) {
    Real c = Real(0.0);
    if (tid < u0) {
        uint o = u_model_obs[u1 + tid];
        uint img = u_obs_image[o];
        uint pt = u_obs_point[o];
        uint4 gi = u_group_info[u_image_group[img]];

        Real pose[6];
        for (int i = 0; i < 6; i++) pose[i] = u_poses[6 * img + i];
        Real intr[M.kNumIntr];
        for (int i = 0; i < M.kNumIntr; i++) intr[i] = u_intr[gi.x + i];
        RVec3 X = rvec3(u_points[3 * pt + 0], u_points[3 * pt + 1], u_points[3 * pt + 2]);
        RVec2 obs = rvec2(u_obs[2 * o + 0], u_obs[2 * o + 1]);

        RVec2 r = reproject<M>(pose, intr, X, obs);
        Real s = r.x * r.x + r.y * r.y;
        c = Real(0.5) * Loss.cost(s, f0);
    }
    g_sh[lid] = c;
    GroupMemoryBarrierWithGroupSync();
    [ForceUnroll]
    for (uint stride = kCostWg / 2; stride > 0; stride >>= 1) {
        if (lid < stride)
            g_sh[lid] = g_sh[lid] + g_sh[lid + stride];
        GroupMemoryBarrierWithGroupSync();
    }
    if (lid == 0)
        atomicAccum(u_cost_a, 0, g_sh[0]);
}

// =====================================================================
// residual/jacobian evaluation + assembly of S, g, App, Bp, Acp
// =====================================================================

// u0 = obs count for this model, u1 = offset into u_model_obs, f0 = loss param,
// u2 = 1 to defer the S/g accumulation to schur_pair (store Jc and the
// residual per obs instead of issuing heavily contended atomics; for the df
// config the CAS-loop contention dominates the whole kernel)
void jac_impl<M : ICameraModel>(uint tid) {
    if (tid >= u0) return;
    uint o = u_model_obs[u1 + tid];
    uint img = u_obs_image[o];
    uint pt = u_obs_point[o];
    uint4 gi = u_group_info[u_image_group[img]];

    Real pose[6];
    for (int i = 0; i < 6; i++) pose[i] = u_poses[6 * img + i];
    Real intr[M.kNumIntr];
    for (int i = 0; i < M.kNumIntr; i++) intr[i] = u_intr[gi.x + i];
    RVec3 X = rvec3(u_points[3 * pt + 0], u_points[3 * pt + 1], u_points[3 * pt + 2]);
    RVec2 obs = rvec2(u_obs[2 * o + 0], u_obs[2 * o + 1]);

    const int DOF = 6 + M.kNumIntr;    // columns the Jacobian is evaluated over
    // ... and how many of them are assembled. Bundle adjustment refines a
    // leading prefix of a group's intrinsics (gi.z of them); the rest are held
    // constant, own no column of the reduced system, and are simply not written
    // -- COLMAP's refine_principal_point/refine_extra_params, and the whole of
    // EQUIRECTANGULAR (D49/D50). The loops below stay unrolled over the static
    // DOF and predicate on this, so the all-free case costs nothing.
    const uint dofw = 6 + gi.z;

    RVec2 r = reproject<M>(pose, intr, X, obs);
    Real s = r.x * r.x + r.y * r.y;
    Real sw = sqrt(Loss.weight(s, f0));

    // Jacobian rows via two backward passes
    Real Jc[2][DOF];
    Real Jp[2][3];
    {
        var dp_pose = diffPair(pose);
        var dp_intr = diffPair(intr);
        var dp_X = diffPair(X);
        [ForceUnroll]
        for (int row = 0; row < 2; row++) {
            RVec2 seed = rvec2(row == 0 ? Real(1.0) : Real(0.0),
                               row == 0 ? Real(0.0) : Real(1.0));
            bwd_diff(reproject<M>)(dp_pose, dp_intr, dp_X, obs, seed);
            for (int i = 0; i < 6; i++) Jc[row][i] = dp_pose.d[i] * sw;
            for (int i = 0; i < M.kNumIntr; i++) Jc[row][6 + i] = dp_intr.d[i] * sw;
            Jp[row][0] = dp_X.d.x * sw; Jp[row][1] = dp_X.d.y * sw; Jp[row][2] = dp_X.d.z * sw;
        }
    }
    r.x = r.x * sw;
    r.y = r.y * sw;

    // global column indices of this observation's camera-side parameters
    // (strictly increasing: pose columns < intrinsics columns)
    uint cols[DOF];
    for (int i = 0; i < 6; i++) cols[i] = 6 * img + i;
    [ForceUnroll]
    for (int i = 0; i < M.kNumIntr; i++)
        if (uint(i) < gi.z) cols[6 + i] = gi.y + i;

    // A_cp = Jc^T Jp
    uint base = u_acp_off[o];
    [ForceUnroll]
    for (int a = 0; a < DOF; a++)
        if (uint(a) < dofw)
            for (int j = 0; j < 3; j++)
                u_acp[base + 3 * a + j] = Jc[0][a] * Jp[0][j] + Jc[1][a] * Jp[1][j];

    if (u2 != 0) {
        // deferred: schur_pair's diagonal-pair workgroups own the S diagonal
        // blocks and gradient columns exclusively and sum these instead
        uint jb = base / 3 * 2;
        [ForceUnroll]
        for (int a = 0; a < DOF; a++)
            if (uint(a) < dofw) {
                u_Jc[jb + 2 * a] = Jc[0][a];
                u_Jc[jb + 2 * a + 1] = Jc[1][a];
            }
        u_res[2 * o] = r.x;
        u_res[2 * o + 1] = r.y;
    } else {
        // S += Jc^T Jc (lower triangle), g += Jc^T r
        [ForceUnroll]
        for (int a = 0; a < DOF; a++) {
            if (uint(a) >= dofw) continue;
            for (int b = 0; b <= a; b++) {
                Real v = Jc[0][a] * Jc[0][b] + Jc[1][a] * Jc[1][b];
                if (isfinite(v))
                    atomicAccum(u_S_a, pidx(cols[a], cols[b]), v);
            }
            Real gv = Jc[0][a] * r.x + Jc[1][a] * r.y;
            if (isfinite(gv))
                atomicAccum(u_g_a, cols[a], gv);
        }
    }

    // App += Jp^T Jp (full 3x3), Bp += Jp^T r
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            Real v = Jp[0][i] * Jp[0][j] + Jp[1][i] * Jp[1][j];
            if (isfinite(v))
                atomicAccum(u_App_a, 9 * pt + 3 * i + j, v);
        }
        Real bv = Jp[0][i] * r.x + Jp[1][i] * r.y;
        if (isfinite(bv))
            atomicAccum(u_Bp_a, 3 * pt + i, bv);
    }
}

// =====================================================================
// Schur complement (model-agnostic; uses per-obs tables only)
// =====================================================================

// multiply diagonal of packed S by (1 + lambda_cam);  u0 = n_dim, f0 = lambda
[shader("compute")] [numthreads(256, 1, 1)]
void damp_S(uint3 t : SV_DispatchThreadID) {
    uint i = t.x;
    if (i >= u0) return;
    uint idx = pidx(i, i);
    u_S[idx] = u_S[idx] * (Real(1.0) + Real(f0));
}

void obsCamCols<let MD : int>(uint o, out uint dof, out uint base, out uint cols[MD]) {
    uint img = u_obs_image[o];
    uint4 gi = u_group_info[u_image_group[img]];
    dof = 6 + gi.z;
    base = u_acp_off[o];
    for (int i = 0; i < 6; i++) cols[i] = 6 * img + i;
    for (uint i = 0; i < gi.z; i++) cols[6 + i] = gi.y + i;
}

// One thread per observation ("row" i of its point's pair block):
//   S -= A_up (App+l)^-1 A_up^T,  g -= A_up (App+l)^-1 Bp
// Work per thread is linear in the point's track length (one thread per point
// would be quadratic and can exceed the GPU watchdog on long tracks).
// u0 = num_obs, f0 = lambda_point
void schurObsImpl<let MD : int>(uint ti) {
    if (ti >= u0) return;
    uint p = u_obs_point[ti];

    Real A[9];
    for (int i = 0; i < 9; i++) A[i] = u_App[9 * p + i];
    A[0] = A[0] * (Real(1.0) + Real(f0));
    A[4] = A[4] * (Real(1.0) + Real(f0));
    A[8] = A[8] * (Real(1.0) + Real(f0));
    Real Ai[9];
    inverse3(A, Ai);

    uint dof_i, base_i;
    uint cols_i[MD];
    obsCamCols<MD>(ti, dof_i, base_i, cols_i);

    // Mi = Acp_i (App+l)^-1,  dof_i x 3
    Real Mi[3 * MD];
    for (uint a = 0; a < dof_i; a++)
        for (int j = 0; j < 3; j++) {
            Real v = Real(0.0);
            for (int k = 0; k < 3; k++)
                v = v + u_acp[base_i + 3 * a + k] * Ai[3 * k + j];
            Mi[3 * a + j] = v;
        }

    // g -= Mi Bp
    for (uint a = 0; a < dof_i; a++) {
        Real v = Mi[3 * a] * u_Bp[3 * p] + Mi[3 * a + 1] * u_Bp[3 * p + 1]
               + Mi[3 * a + 2] * u_Bp[3 * p + 2];
        if (isfinite(v))
            atomicAccum(u_g_a, cols_i[a], -v);
    }

    // S -= Mi Acp_j^T over all observations j of the same point (ordered
    // pairs across threads); only write row >= col
    uint t0 = u_obs_ranges[p];
    uint t1 = u_obs_ranges[p + 1];
    for (uint tj = t0; tj < t1; tj++) {
        uint dof_j, base_j;
        uint cols_j[MD];
        obsCamCols<MD>(tj, dof_j, base_j, cols_j);
        for (uint c = 0; c < dof_j; c++) {
            uint col = cols_j[c];
            Real b0 = u_acp[base_j + 3 * c];
            Real b1 = u_acp[base_j + 3 * c + 1];
            Real b2 = u_acp[base_j + 3 * c + 2];
            for (uint a = 0; a < dof_i; a++) {
                uint row = cols_i[a];
                if (row < col) continue;
                Real v = Mi[3 * a] * b0 + Mi[3 * a + 1] * b1 + Mi[3 * a + 2] * b2;
                if (isfinite(v))
                    atomicAccum(u_S_a, pidx(row, col), -v);
            }
        }
    }
}

// Precompute per-point W = (App + l)^-1 (used by acp_prep and point_update).
// u0 = num_points, f0 = lambda_point
[shader("compute")] [numthreads(256, 1, 1)]
void point_prep(uint3 t : SV_DispatchThreadID) {
    uint p = t.x;
    if (p >= u0) return;
    Real A[9];
    for (int i = 0; i < 9; i++) A[i] = u_App[9 * p + i];
    A[0] = A[0] * (Real(1.0) + Real(f0));
    A[4] = A[4] * (Real(1.0) + Real(f0));
    A[8] = A[8] * (Real(1.0) + Real(f0));
    Real Ai[9];
    inverse3(A, Ai);
    for (int i = 0; i < 9; i++) u_W[9 * p + i] = Ai[i];
}

// Precompute T = Acp W per observation (T depends only on the obs since its
// point is fixed), so the pair-Schur kernel is a plain dot product per entry.
// u0 = num_obs
[shader("compute")] [numthreads(256, 1, 1)]
void acp_prep(uint3 t : SV_DispatchThreadID) {
    uint o = t.x;
    if (o >= u0) return;
    uint p = u_obs_point[o];
    uint dof = 6 + u_group_info[u_image_group[u_obs_image[o]]].z;
    uint base = u_acp_off[o];
    Real W[9];
    for (int i = 0; i < 9; i++) W[i] = u_W[9 * p + i];
    for (uint r = 0; r < dof; r++) {
        Real a0 = u_acp[base + 3 * r], a1 = u_acp[base + 3 * r + 1], a2 = u_acp[base + 3 * r + 2];
        u_T[base + 3 * r]     = a0 * W[0] + a1 * W[3] + a2 * W[6];
        u_T[base + 3 * r + 1] = a0 * W[1] + a1 * W[4] + a2 * W[7];
        u_T[base + 3 * r + 2] = a0 * W[2] + a1 * W[5] + a2 * W[8];
    }
}

// Pair-aggregated Schur complement: one workgroup (warp) per chunk of one
// unordered image pair {a,b}. Each element of S belongs to exactly one image
// pair (per-image intrinsics groups; checked host-side), so the accumulated
// dof_i x dof_j block is written back without atomics unless the pair was
// split across chunks (flag bit 31). Diagonal pairs (a==b; their entries are
// exactly the observations of that camera) additionally accumulate the
// deferred J^T J block and gradient J^T r - Acp W Bp from u_Jc/u_res, and
// apply the multiplicative damping (1+lambda) to the global diagonal, so
// neither the Jacobian kernel nor damp_S touches S/g on this path.
// u0 = num chunks, f0 = lambda
static const uint kPairWg = 32;
// Wide-tier element count, used by the (untiered) CG pair kernel in cg.slang.
static const uint kPairEl = (kMaxCamDof * kMaxCamDof + kPairWg - 1) / kPairWg;

void schurPairImpl<let MD : int>(uint wgx, uint gi) {
    const uint kEl = (MD * MD + kPairWg - 1) / kPairWg;
    if (wgx >= u0) return;
    uint2 ck = u_pair_chunks[wgx];
    uint off = ck.x;
    uint count = ck.y & 0x7FFFFFFFu;
    bool useAtomic = (ck.y & 0x80000000u) != 0;

    uint2 e0 = u_pair_entries[off];
    uint dof_i, base0_i, dof_j, base0_j;
    uint cols_i[MD], cols_j[MD];
    obsCamCols<MD>(e0.x, dof_i, base0_i, cols_i);
    obsCamCols<MD>(e0.y, dof_j, base0_j, cols_j);
    bool diag = e0.x == e0.y;

    // static element ownership: element el = r*dof_j + c owned by lane el%32
    uint nel = dof_i * dof_j;
    Real acc[kEl];
    Real acc2[kEl];  // deferred J^T J (diagonal pairs only)
    uint elr[kEl], elc[kEl];
    [ForceUnroll]
    for (uint k = 0; k < kEl; k++) {
        acc[k] = Real(0.0);
        acc2[k] = Real(0.0);
        uint el = gi + kPairWg * k;
        elr[k] = el / dof_j;
        elc[k] = el - elr[k] * dof_j;
    }
    Real gacc = Real(0.0);

    for (uint t = 0; t < count; t++) {
        uint2 e = u_pair_entries[off + t];
        uint bi = u_acp_off[e.x];
        uint bj = u_acp_off[e.y];
        uint jb = bi / 3 * 2;
        [ForceUnroll]
        for (uint k = 0; k < kEl; k++) {
            uint el = gi + kPairWg * k;
            if (el >= nel) continue;
            uint r = elr[k], c = elc[k];
            if (diag && r < c) continue;  // symmetric block: lower triangle only
            acc[k] = acc[k] + u_T[bi + 3 * r] * u_acp[bj + 3 * c]
                            + u_T[bi + 3 * r + 1] * u_acp[bj + 3 * c + 1]
                            + u_T[bi + 3 * r + 2] * u_acp[bj + 3 * c + 2];
            if (diag)
                acc2[k] = acc2[k] + u_Jc[jb + 2 * r] * u_Jc[jb + 2 * c]
                                  + u_Jc[jb + 2 * r + 1] * u_Jc[jb + 2 * c + 1];
        }
        if (diag && gi < dof_i) {
            uint p = u_obs_point[e.x];
            gacc = gacc + u_T[bi + 3 * gi] * u_Bp[3 * p]
                        + u_T[bi + 3 * gi + 1] * u_Bp[3 * p + 1]
                        + u_T[bi + 3 * gi + 2] * u_Bp[3 * p + 2]
                        - u_Jc[jb + 2 * gi] * u_res[2 * e.x]
                        - u_Jc[jb + 2 * gi + 1] * u_res[2 * e.x + 1];
        }
    }

    [ForceUnroll]
    for (uint k = 0; k < kEl; k++) {
        uint el = gi + kPairWg * k;
        if (el >= nel) continue;
        uint r = elr[k], c = elc[k];
        if (diag && r < c) continue;
        Real v = -acc[k];
        if (diag) {
            Real jj = acc2[k];
            if (r == c)
                jj = jj * (Real(1.0) + Real(f0));  // damp_S folded in
            v = v + jj;
        }
        if (!isfinite(v)) continue;
        uint R = cols_i[r], C = cols_j[c];
        uint idx = R >= C ? pidx(R, C) : pidx(C, R);
        if (useAtomic)
            atomicAccum(u_S_a, idx, v);
        else
            u_S[idx] = u_S[idx] + v;
    }
    if (diag && gi < dof_i && isfinite(gacc)) {
        uint col = cols_i[gi];
        if (useAtomic)
            atomicAccum(u_g_a, col, -gacc);
        else
            u_g[col] = u_g[col] - gacc;
    }
}

// =====================================================================
// parameter updates (after cholesky solve leaves dU in u_g)
// =====================================================================

// Bp -= Acp^T dU, one thread per observation; u0 = num_obs
void dpAccumImpl<let MD : int>(uint o) {
    if (o >= u0) return;
    uint dof, base;
    uint cols[MD];
    obsCamCols<MD>(o, dof, base, cols);
    uint pt = u_obs_point[o];
    for (int j = 0; j < 3; j++) {
        Real v = Real(0.0);
        for (uint a = 0; a < dof; a++)
            v = v + u_acp[base + 3 * a + j] * u_g[cols[a]];
        if (isfinite(v))
            atomicAccum(u_Bp_a, 3 * pt + j, -v);
    }
}

// dP = (App+l)^-1 Bp; points -= dP.  W was precomputed by point_prep with the
// same lambda; Bp has since been updated in place by dp_accum.  u0 = num_points
[shader("compute")] [numthreads(256, 1, 1)]
void point_update(uint3 t : SV_DispatchThreadID) {
    uint p = t.x;
    if (p >= u0) return;
    for (int i = 0; i < 3; i++) {
        Real d = u_W[9 * p + 3 * i] * u_Bp[3 * p] + u_W[9 * p + 3 * i + 1] * u_Bp[3 * p + 1]
               + u_W[9 * p + 3 * i + 2] * u_Bp[3 * p + 2];
        if (isfinite(d))
            u_points[3 * p + i] = u_points[3 * p + i] - d;
    }
}

// poses/intrinsics -= dU.  u0 = n_dim, u1 = 6*num_images
[shader("compute")] [numthreads(256, 1, 1)]
void cam_update(uint3 t : SV_DispatchThreadID) {
    uint i = t.x;
    if (i >= u0) return;
    Real d = u_g[i];
    if (!isfinite(d)) return;
    u_poses[i] = u_poses[i] - d;
}

// Apply the intrinsics half of the step: one thread per camera group, walking
// its free parameters. Column j of a group maps to its slot j -- but only
// within the group, since a group's stored parameters are contiguous while its
// *columns* are packed against the other groups' (a group whose principal point
// is held owns fewer columns than it stores). Going through the group table is
// what keeps those two layouts independent (D50).  u0 = num groups
[shader("compute")] [numthreads(64, 1, 1)]
void intr_update(uint3 t : SV_DispatchThreadID) {
    uint g = t.x;
    if (g >= u0) return;
    uint4 gi = u_group_info[g];
    for (uint j = 0; j < gi.z; j++) {
        Real d = u_g[gi.y + j];
        if (isfinite(d)) u_intr[gi.x + j] = u_intr[gi.x + j] - d;
    }
}

// =====================================================================
// Schur entry points, two dof tiers (compact / wide); the host dispatches
// the tier that fits the problem's max camera dof (sfm/ba/Solver.h).
// =====================================================================

[shader("compute")] [numthreads(128, 1, 1)]
void schur_obs_c(uint3 t : SV_DispatchThreadID) { schurObsImpl<kDofCompact>(t.x); }
[shader("compute")] [numthreads(128, 1, 1)]
void schur_obs_m(uint3 t : SV_DispatchThreadID) { schurObsImpl<kDofMid>(t.x); }
[shader("compute")] [numthreads(128, 1, 1)]
void schur_obs_w(uint3 t : SV_DispatchThreadID) { schurObsImpl<kDofWide>(t.x); }

[shader("compute")] [numthreads(kPairWg, 1, 1)]
void schur_pair_c(uint3 wg : SV_GroupID, uint gi : SV_GroupIndex) { schurPairImpl<kDofCompact>(wg.x, gi); }
[shader("compute")] [numthreads(kPairWg, 1, 1)]
void schur_pair_m(uint3 wg : SV_GroupID, uint gi : SV_GroupIndex) { schurPairImpl<kDofMid>(wg.x, gi); }
[shader("compute")] [numthreads(kPairWg, 1, 1)]
void schur_pair_w(uint3 wg : SV_GroupID, uint gi : SV_GroupIndex) { schurPairImpl<kDofWide>(wg.x, gi); }

[shader("compute")] [numthreads(256, 1, 1)]
void dp_accum_c(uint3 t : SV_DispatchThreadID) { dpAccumImpl<kDofCompact>(t.x); }
[shader("compute")] [numthreads(256, 1, 1)]
void dp_accum_m(uint3 t : SV_DispatchThreadID) { dpAccumImpl<kDofMid>(t.x); }
[shader("compute")] [numthreads(256, 1, 1)]
void dp_accum_w(uint3 t : SV_DispatchThreadID) { dpAccumImpl<kDofWide>(t.x); }

// =====================================================================
// per-model entry points -- add two lines per new camera model
// =====================================================================

[shader("compute")] [numthreads(256, 1, 1)]
void cost_snavely(uint3 t : SV_DispatchThreadID, uint gi : SV_GroupIndex) { cost_impl<SnavelyModel>(t.x, gi); }
[shader("compute")] [numthreads(128, 1, 1)]
void jac_snavely(uint3 t : SV_DispatchThreadID) { jac_impl<SnavelyModel>(t.x); }

[shader("compute")] [numthreads(256, 1, 1)]
void cost_snavely_f(uint3 t : SV_DispatchThreadID, uint gi : SV_GroupIndex) { cost_impl<SnavelyFModel>(t.x, gi); }
[shader("compute")] [numthreads(128, 1, 1)]
void jac_snavely_f(uint3 t : SV_DispatchThreadID) { jac_impl<SnavelyFModel>(t.x); }

[shader("compute")] [numthreads(256, 1, 1)]
void cost_pinhole_radial(uint3 t : SV_DispatchThreadID, uint gi : SV_GroupIndex) { cost_impl<PinholeRadialModel>(t.x, gi); }
[shader("compute")] [numthreads(128, 1, 1)]
void jac_pinhole_radial(uint3 t : SV_DispatchThreadID) { jac_impl<PinholeRadialModel>(t.x); }

[shader("compute")] [numthreads(256, 1, 1)]
void cost_opencv(uint3 t : SV_DispatchThreadID, uint gi : SV_GroupIndex) { cost_impl<OpenCVModel>(t.x, gi); }
[shader("compute")] [numthreads(128, 1, 1)]
void jac_opencv(uint3 t : SV_DispatchThreadID) { jac_impl<OpenCVModel>(t.x); }

[shader("compute")] [numthreads(256, 1, 1)]
void cost_simple_pinhole(uint3 t : SV_DispatchThreadID, uint gi : SV_GroupIndex) { cost_impl<SimplePinholeModel>(t.x, gi); }
[shader("compute")] [numthreads(128, 1, 1)]
void jac_simple_pinhole(uint3 t : SV_DispatchThreadID) { jac_impl<SimplePinholeModel>(t.x); }

[shader("compute")] [numthreads(256, 1, 1)]
void cost_pinhole(uint3 t : SV_DispatchThreadID, uint gi : SV_GroupIndex) { cost_impl<PinholeModel>(t.x, gi); }
[shader("compute")] [numthreads(128, 1, 1)]
void jac_pinhole(uint3 t : SV_DispatchThreadID) { jac_impl<PinholeModel>(t.x); }

[shader("compute")] [numthreads(256, 1, 1)]
void cost_opencv_fisheye(uint3 t : SV_DispatchThreadID, uint gi : SV_GroupIndex) { cost_impl<FisheyeModel>(t.x, gi); }
[shader("compute")] [numthreads(128, 1, 1)]
void jac_opencv_fisheye(uint3 t : SV_DispatchThreadID) { jac_impl<FisheyeModel>(t.x); }

[shader("compute")] [numthreads(256, 1, 1)]
void cost_full_opencv(uint3 t : SV_DispatchThreadID, uint gi : SV_GroupIndex) { cost_impl<FullOpenCVModel>(t.x, gi); }
[shader("compute")] [numthreads(128, 1, 1)]
void jac_full_opencv(uint3 t : SV_DispatchThreadID) { jac_impl<FullOpenCVModel>(t.x); }

[shader("compute")] [numthreads(256, 1, 1)]
void cost_equirect(uint3 t : SV_DispatchThreadID, uint gi : SV_GroupIndex) { cost_impl<EquirectModel>(t.x, gi); }
[shader("compute")] [numthreads(128, 1, 1)]
void jac_equirect(uint3 t : SV_DispatchThreadID) { jac_impl<EquirectModel>(t.x); }

[shader("compute")] [numthreads(256, 1, 1)]
void cost_thin_prism_fisheye(uint3 t : SV_DispatchThreadID, uint gi : SV_GroupIndex) { cost_impl<ThinPrismFisheyeModel>(t.x, gi); }
[shader("compute")] [numthreads(128, 1, 1)]
void jac_thin_prism_fisheye(uint3 t : SV_DispatchThreadID) { jac_impl<ThinPrismFisheyeModel>(t.x); }

#include "cholesky.slang"
#include "cg.slang"
