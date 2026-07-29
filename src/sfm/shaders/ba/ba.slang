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
// The Jacobian pass stores only what it evaluated: Jc (2 x dof), Jp (2 x 3) and
// the weighted residual, per observation. Every consumer of the camera-point
// block A_cp = Jc^T Jp forms it implicitly instead of reading a materialized
// dof x 3 block -- the products all factor through the 2-vector residual
// dimension rather than the 3-vector point dimension:
//
//   A_cp^T x       = Jp^T (Jc x)                      (cg_gather, dp_accum)
//   A_cp v         = Jc^T (Jp v)                      (cg_scatter)
//   A_cp W A_cp'^T = Jc^T (Jp W Jp'^T) Jc'            (every Schur kernel)
//
// which is both smaller (3*dof scalars per observation replaced by 6, and the
// same again for the prepass product T = A_cp W) and cheaper, since dof > 3.
// `u_Y` holds (Jp W)^T, the one piece the pair kernel wants precomputed.
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
[[vk::binding(10)]] StructuredBuffer<uint> u_jc_off;       // (per obs) element offset into u_Jc

// ---- work buffers (double-bound: Real view + atomic-word view) ----
[[vk::binding(11)]] RWStructuredBuffer<Real> u_Jp;         // point Jacobian, [6*o + 2*j + row]
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
[[vk::binding(25)]] RWStructuredBuffer<Real> u_Y;           // (Jp W)^T per obs, [6*o + 2*j + row]
[[vk::binding(26)]] RWStructuredBuffer<Real> u_y;           // (n_dim) triangular-solve scratch
[[vk::binding(27)]] RWStructuredBuffer<Real> u_Jc;          // camera Jacobian rows per obs,
                                                            // [2*a+row], offset u_jc_off
[[vk::binding(28)]] RWStructuredBuffer<Real> u_res;         // weighted residual (2 per obs)

// packed lower-triangular index, row >= col
uint pidx(uint r, uint c) { return r * (r + 1) / 2 + c; }

// =====================================================================
// cost
// =====================================================================

static const uint kCostWg = 256;

// Three groupshared pools, sized to what actually uses them. The *sum* of a
// module's Workgroup variables counts against maxComputeSharedMemorySize, so
// they cannot be handed out per kernel freely -- but they must not be one pool
// either: a driver allocates each pipeline the shared memory its entry point
// reaches, and a 32-thread workgroup carrying the Cholesky pool's 16.9 kB
// (double) fits about five to an SM. That costs nothing for a kernel that
// never touches shared memory, and 2.2x for schur_pair, which lives on it
// (measured, RTX 4080 Super). The three together are ~26 kB at fp64, inside
// every device's 32 kB floor.
//
// g_chol - the Cholesky tiles: two 32x32 blocks, the second at kSh2 with row
//          pitch 33 to avoid bank conflicts, plus a vector for the solves
// g_pair - schur_pair's per-entry staging (see kPr* below)
// g_sh   - everything small: block reductions, and the CG kernels' per-obs and
//          per-chunk staging
static const uint kSh2 = 32 * 32;
static const uint kShPool = 32 * 32 + 32 * 33 + 32;
groupshared Real g_chol[kShPool];
groupshared Real g_sh[512];

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
// residual/jacobian evaluation: Jc, Jp, weighted residual, App, Bp
// =====================================================================

// The S/g accumulation is *not* done here on any path -- the Schur kernels own
// it, since a workgroup that already holds a camera's whole observation list
// can sum the J^T J block in registers instead of issuing ~dof^2 heavily
// contended atomics per observation (for the df config the CAS-loop contention
// dominated this kernel). It also makes the pass idempotent in the parameters:
// a rejected step restores them and every downstream kernel can simply re-run,
// which is why no path snapshots S.
// u0 = obs count for this model, u1 = offset into u_model_obs, f0 = loss param
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

    uint jb = u_jc_off[o];
    [ForceUnroll]
    for (int a = 0; a < DOF; a++)
        if (uint(a) < dofw) {
            u_Jc[jb + 2 * a] = Jc[0][a];
            u_Jc[jb + 2 * a + 1] = Jc[1][a];
        }
    for (int j = 0; j < 3; j++) {
        u_Jp[6 * o + 2 * j] = Jp[0][j];
        u_Jp[6 * o + 2 * j + 1] = Jp[1][j];
    }
    u_res[2 * o] = r.x;
    u_res[2 * o + 1] = r.y;

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

void obsCamCols<let MD : int>(uint o, out uint dof, out uint base, out uint cols[MD]) {
    uint img = u_obs_image[o];
    uint4 gi = u_group_info[u_image_group[img]];
    dof = 6 + gi.z;
    base = u_jc_off[o];
    for (int i = 0; i < 6; i++) cols[i] = 6 * img + i;
    for (uint i = 0; i < gi.z; i++) cols[6 + i] = gi.y + i;
}

// One thread per observation ("row" i of its point's pair block):
//   S += Jc^T Jc,  g += Jc^T r,  S -= A_up (App+l)^-1 A_up^T,  g -= A_up W Bp
// Work per thread is linear in the point's track length (one thread per point
// would be quadratic and can exceed the GPU watchdog on long tracks).
//
// The atomic fallback for the Schur assembly: the pair kernel below aggregates
// the same sums per image pair and needs no atomics, but it requires exclusive
// column ownership. The (1+lambda) damping of the J^T J diagonal is folded in
// here exactly as it is there.
// u0 = num_obs, f0 = lambda_point
void schurObsImpl<let MD : int>(uint ti) {
    if (ti >= u0) return;
    uint p = u_obs_point[ti];

    // W = (App + lambda)^-1, already inverted per point by point_prep
    Real W[9];
    for (int i = 0; i < 9; i++) W[i] = u_W[9 * p + i];

    uint dof_i, jb_i;
    uint cols_i[MD];
    obsCamCols<MD>(ti, dof_i, jb_i, cols_i);
    Real Jci[2 * MD];
    for (uint a = 0; a < 2 * dof_i; a++) Jci[a] = u_Jc[jb_i + a];

    // Y = (Jp_i W)^T, 2x3 stored column-major like u_Jp
    Real Y[6];
    for (int j = 0; j < 3; j++)
        for (int row = 0; row < 2; row++)
            Y[2 * j + row] = u_Jp[6 * ti + row] * W[j]
                           + u_Jp[6 * ti + 2 + row] * W[3 + j]
                           + u_Jp[6 * ti + 4 + row] * W[6 + j];

    // S += Jc^T Jc (damped diagonal), g += Jc^T r - Jc^T (Y Bp)
    Real yb0 = Y[0] * u_Bp[3 * p] + Y[2] * u_Bp[3 * p + 1] + Y[4] * u_Bp[3 * p + 2];
    Real yb1 = Y[1] * u_Bp[3 * p] + Y[3] * u_Bp[3 * p + 1] + Y[5] * u_Bp[3 * p + 2];
    Real r0 = u_res[2 * ti], r1 = u_res[2 * ti + 1];
    for (uint a = 0; a < dof_i; a++) {
        for (uint b = 0; b <= a; b++) {
            Real v = Jci[2 * a] * Jci[2 * b] + Jci[2 * a + 1] * Jci[2 * b + 1];
            if (a == b) v = v * (Real(1.0) + Real(f0));
            if (isfinite(v))
                atomicAccum(u_S_a, pidx(cols_i[a], cols_i[b]), v);
        }
        Real gv = Jci[2 * a] * (r0 - yb0) + Jci[2 * a + 1] * (r1 - yb1);
        if (isfinite(gv))
            atomicAccum(u_g_a, cols_i[a], gv);
    }

    // S -= Jc_i^T (Jp_i W Jp_j^T) Jc_j over all observations j of the same
    // point (ordered pairs across threads); only write row >= col
    uint t0 = u_obs_ranges[p];
    uint t1 = u_obs_ranges[p + 1];
    for (uint tj = t0; tj < t1; tj++) {
        uint dof_j, jb_j;
        uint cols_j[MD];
        obsCamCols<MD>(tj, dof_j, jb_j, cols_j);
        // Q = Y Jp_j^T, 2x2
        Real q00 = Real(0.0), q01 = Real(0.0), q10 = Real(0.0), q11 = Real(0.0);
        for (int j = 0; j < 3; j++) {
            Real b0 = u_Jp[6 * tj + 2 * j], b1 = u_Jp[6 * tj + 2 * j + 1];
            q00 = q00 + Y[2 * j] * b0;
            q01 = q01 + Y[2 * j] * b1;
            q10 = q10 + Y[2 * j + 1] * b0;
            q11 = q11 + Y[2 * j + 1] * b1;
        }
        for (uint a = 0; a < dof_i; a++) {
            uint row = cols_i[a];
            Real z0 = Jci[2 * a] * q00 + Jci[2 * a + 1] * q10;
            Real z1 = Jci[2 * a] * q01 + Jci[2 * a + 1] * q11;
            for (uint c = 0; c < dof_j; c++) {
                uint col = cols_j[c];
                if (row < col) continue;
                Real v = z0 * u_Jc[jb_j + 2 * c] + z1 * u_Jc[jb_j + 2 * c + 1];
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

// Precompute Y = (Jp W)^T per observation (W depends only on the obs's point),
// so the pair-Schur kernel's per-entry 2x2 middle matrix is one small product.
// u0 = num_obs
[shader("compute")] [numthreads(256, 1, 1)]
void y_prep(uint3 t : SV_DispatchThreadID) {
    uint o = t.x;
    if (o >= u0) return;
    uint p = u_obs_point[o];
    Real W[9];
    for (int i = 0; i < 9; i++) W[i] = u_W[9 * p + i];
    Real Jp[6];
    for (int i = 0; i < 6; i++) Jp[i] = u_Jp[6 * o + i];
    for (int j = 0; j < 3; j++)
        for (int row = 0; row < 2; row++)
            u_Y[6 * o + 2 * j + row] =
                Jp[row] * W[j] + Jp[2 + row] * W[3 + j] + Jp[4 + row] * W[6 + j];
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

// Groupshared staging for the per-entry 2x2 middle matrix Q = Y_i Jp_j^T (plus
// the 2-vector Y_i Bp the diagonal pairs' gradient wants) and the
// half-contracted Z = Q^T Jc_i, so each element of the block costs one 2-term
// dot product (see the A_cp note at the top of the file).
//
// Staged kPrB entries at a time: the cooperative steps need a barrier each, and
// at one entry per barrier the synchronization cost buried the arithmetic
// (measured 2.3x slower than the materialized-A_cp kernel it replaced, which
// had none). Batching amortizes them away while the per-element work stays at
// two multiplies instead of three.
static const uint kPrB = 8;         // entries staged per barrier round
static const uint kPrQ = 0;         // 6 per entry: Q (4) then Y_i Bp (2)
static const uint kPrZ = 6 * kPrB;  // 2*MD per entry, then Jc_j the same again
groupshared Real g_pair[6 * kPrB + 4 * kPrB * kMaxCamDof];

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

    const uint kZs = 2 * MD;              // per-entry stride of the staged blocks
    const uint kPrJ = kPrZ + kPrB * kZs;  // staged Jc_j
    for (uint t0 = 0; t0 < count; t0 += kPrB) {
        uint nb = min(kPrB, count - t0);
        GroupMemoryBarrierWithGroupSync();  // previous batch done being read
        for (uint i = gi; i < nb * 6; i += kPairWg) {
            uint b = i / 6, m = i - 6 * b;
            uint2 e = u_pair_entries[off + t0 + b];
            Real v;
            if (m < 4) {  // Q[qr][qs] = sum_j Y_i[j][qr] Jp_j[j][qs]
                uint qr = m >> 1, qs = m & 1;
                v = u_Y[6 * e.x + qr] * u_Jp[6 * e.y + qs]
                  + u_Y[6 * e.x + 2 + qr] * u_Jp[6 * e.y + 2 + qs]
                  + u_Y[6 * e.x + 4 + qr] * u_Jp[6 * e.y + 4 + qs];
            } else if (diag) {  // (Y_i Bp)[row], only the gradient wants it
                uint row = m - 4;
                uint p = u_obs_point[e.x];
                v = u_Y[6 * e.x + row] * u_Bp[3 * p] + u_Y[6 * e.x + 2 + row] * u_Bp[3 * p + 1]
                  + u_Y[6 * e.x + 4 + row] * u_Bp[3 * p + 2];
            } else {
                continue;
            }
            g_pair[kPrQ + 6 * b + m] = v;
        }
        for (uint i = gi; i < nb * 2 * dof_j; i += kPairWg) {
            uint b = i / (2 * dof_j), m = i - b * 2 * dof_j;
            g_pair[kPrJ + kZs * b + m] = u_Jc[u_jc_off[u_pair_entries[off + t0 + b].y] + m];
        }
        GroupMemoryBarrierWithGroupSync();
        for (uint i = gi; i < nb * 2 * dof_i; i += kPairWg) {
            uint b = i / (2 * dof_i), m = i - b * 2 * dof_i;
            uint a = m >> 1, s = m & 1;
            uint bi = u_jc_off[u_pair_entries[off + t0 + b].x];
            g_pair[kPrZ + kZs * b + m] = u_Jc[bi + 2 * a] * g_pair[kPrQ + 6 * b + s]
                                       + u_Jc[bi + 2 * a + 1] * g_pair[kPrQ + 6 * b + 2 + s];
        }
        GroupMemoryBarrierWithGroupSync();
        for (uint b = 0; b < nb; b++) {
            uint zb = kPrZ + kZs * b, jb = kPrJ + kZs * b;
            // Jc_i stays in global for the deferred J^T J block: a diagonal
            // pair has Jc_i == Jc_j staged right there, but taking it from
            // shared as well costs 45% -- two shared operands per product is
            // one too many, and these reads all hit L1 anyway.
            uint bi = u_jc_off[u_pair_entries[off + t0 + b].x];
            [ForceUnroll]
            for (uint k = 0; k < kEl; k++) {
                uint el = gi + kPairWg * k;
                if (el >= nel) continue;
                uint r = elr[k], c = elc[k];
                if (diag && r < c) continue;  // symmetric block: lower triangle only
                acc[k] = acc[k] + g_pair[zb + 2 * r] * g_pair[jb + 2 * c]
                                + g_pair[zb + 2 * r + 1] * g_pair[jb + 2 * c + 1];
                if (diag)
                    acc2[k] = acc2[k] + u_Jc[bi + 2 * r] * u_Jc[bi + 2 * c]
                                      + u_Jc[bi + 2 * r + 1] * u_Jc[bi + 2 * c + 1];
            }
            if (diag && gi < dof_i) {
                uint ox = u_pair_entries[off + t0 + b].x;
                gacc = gacc + u_Jc[bi + 2 * gi] * (g_pair[kPrQ + 6 * b + 4] - u_res[2 * ox])
                            + u_Jc[bi + 2 * gi + 1] * (g_pair[kPrQ + 6 * b + 5] - u_res[2 * ox + 1]);
            }
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

// Bp -= Acp^T dU = Jp^T (Jc dU), one thread per observation; u0 = num_obs
void dpAccumImpl<let MD : int>(uint o) {
    if (o >= u0) return;
    uint dof, base;
    uint cols[MD];
    obsCamCols<MD>(o, dof, base, cols);
    uint pt = u_obs_point[o];
    Real d0 = Real(0.0), d1 = Real(0.0);
    for (uint a = 0; a < dof; a++) {
        Real x = u_g[cols[a]];
        d0 = d0 + u_Jc[base + 2 * a] * x;
        d1 = d1 + u_Jc[base + 2 * a + 1] * x;
    }
    for (int j = 0; j < 3; j++) {
        Real v = u_Jp[6 * o + 2 * j] * d0 + u_Jp[6 * o + 2 * j + 1] * d1;
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
