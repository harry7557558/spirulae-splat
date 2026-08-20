// Deadzone-quadratic distortion prior, the device half of sfm/ba/CamPrior.h --
// same quadrature, same metric, same Gauss-Newton blocks. The prior sees no
// point, so its contribution is untouched by the Schur complement and simply
// adds onto the group's intrinsics columns after the assembly kernels.
//
// One thread per camera group. Included by ba.slang after cg.slang, since it
// writes the reduced system of whichever path ran.

[[vk::binding(45)]] StructuredBuffer<Real> u_prior;         // kPriorN per group
[[vk::binding(46)]] StructuredBuffer<uint> u_prior_layout;  // kPriorL per group

// Field offsets into the two per-group records; must match packPrior() in
// sfm/ba/Solver.h.
static const uint kPriorN = 4;
static const uint kPrW = 0, kPrDead = 1, kPrRmax = 2;
static const uint kPriorL = 20;
static const uint kPlNrad = 0, kPlNden = 1, kPlNtan = 2, kPlNprism = 3;
static const uint kPlRad = 4, kPlDen = 8, kPlTan = 11, kPlPrism = 13;
static const uint kPlImage = 15, kPlNparam = 16;

// Node and weight as a float pair, the df.slang idiom: a plain float literal
// would leave the double config's quadrature accurate to 1e-8 and every
// CPU/GPU comparison of the prior with it.
static const float kPriorXHi[6] = {0.0337652415f, 0.169395313f, 0.380690396f,
                                   0.619309604f, 0.830604672f, 0.966234744f};
static const float kPriorXLo[6] = {1.39470846e-09f, -5.90002536e-09f, 1.11263398e-08f,
                                   -1.11263398e-08f, 2.08011866e-08f, 1.35064528e-08f};
static const float kPriorWHi[6] = {0.0856622458f, 0.180380791f, 0.233956963f,
                                   0.233956963f, 0.180380791f, 0.0856622458f};
static const float kPriorWLo[6] = {4.39157932e-10f, -4.90163554e-09f, 4.4624775e-09f,
                                   4.4624775e-09f, -4.90163554e-09f, 4.39157932e-10f};

void priorEval(uint gg, out Real cost, out Real go[12], out Real Ho[144]) {
    uint4 gi = u_group_info[gg];
    uint pb = kPriorN * gg, lb = kPriorL * gg;
    cost = Real(0.0);
    for (uint i = 0; i < 12; i++) go[i] = Real(0.0);
    for (uint i = 0; i < 144; i++) Ho[i] = Real(0.0);

    uint nrad = u_prior_layout[lb + kPlNrad], nden = u_prior_layout[lb + kPlNden];
    uint ntan = u_prior_layout[lb + kPlNtan], nprism = u_prior_layout[lb + kPlNprism];
    Real wd = u_prior[pb + kPrW], rmax = u_prior[pb + kPrRmax];
    if (!(wd > Real(0.0) && rmax > Real(0.0) && (nrad + ntan + nprism) > 0)) return;

    uint nfree = min(gi.z, 12u);
    uint np = u_prior_layout[lb + kPlNparam];
    Real p[12];
    for (uint i = 0; i < 12; i++) p[i] = i < np ? u_intr[gi.x + i] : Real(0.0);

    Real u[12], A[144], d[12];
    for (uint i = 0; i < 12; i++) u[i] = Real(0.0);
    for (uint i = 0; i < 144; i++) A[i] = Real(0.0);
    Real s2 = Real(0.0);
    const Real R2 = rmax * rmax;
    for (uint q = 0; q < 6; q++) {
        const Real r2 = R2 * (Real(kPriorXHi[q]) + Real(kPriorXLo[q]));
        const Real w = Real(kPriorWHi[q]) + Real(kPriorWLo[q]);
        Real num = Real(1.0), den = Real(1.0), rp = Real(1.0);
        for (uint i = 0; i < nrad; i++) {
            rp = rp * r2;
            num = num + p[u_prior_layout[lb + kPlRad + i]] * rp;
        }
        rp = Real(1.0);
        for (uint i = 0; i < nden; i++) {
            rp = rp * r2;
            den = den + p[u_prior_layout[lb + kPlDen + i]] * rp;
        }
        for (uint i = 0; i < 12; i++) d[i] = Real(0.0);
        rp = Real(1.0);
        for (uint i = 0; i < nrad; i++) {
            rp = rp * r2;
            d[u_prior_layout[lb + kPlRad + i]] = rp / den;
        }
        rp = Real(1.0);
        for (uint i = 0; i < nden; i++) {
            rp = rp * r2;
            d[u_prior_layout[lb + kPlDen + i]] = -(num * rp) / (den * den);
        }
        const Real dev = num / den - Real(1.0);
        s2 = s2 + w * dev * dev;
        for (uint i = 0; i < 12; i++) {
            u[i] = u[i] + w * d[i] * dev;
            for (uint j = 0; j < 12; j++) A[12 * i + j] = A[12 * i + j] + w * d[i] * d[j];
        }
    }
    const Real ct = Real(2.5) * R2, cp = Real(0.5) * R2;
    for (uint i = 0; i < ntan; i++) {
        uint k = u_prior_layout[lb + kPlTan + i];
        s2 = s2 + ct * p[k] * p[k];
        u[k] = u[k] + ct * p[k];
        A[12 * k + k] = A[12 * k + k] + ct;
    }
    for (uint i = 0; i < nprism; i++) {
        uint k = u_prior_layout[lb + kPlPrism + i];
        s2 = s2 + cp * p[k] * p[k];
        u[k] = u[k] + cp * p[k];
        A[12 * k + k] = A[12 * k + k] + cp;
    }
    const Real s = sqrt(s2), dead = u_prior[pb + kPrDead];
    if (!(s > dead)) return;
    const Real e = s - dead, a = e / s;
    cost = Real(0.5) * wd * e * e;
    const Real rank1 = wd * dead / (s * s * s);
    for (uint i = 0; i < nfree; i++) {
        go[i] = go[i] + wd * a * u[i];
        for (uint j = 0; j < nfree; j++)
            Ho[12 * i + j] = Ho[12 * i + j] + wd * a * A[12 * i + j] + rank1 * u[i] * u[j];
    }
}

// u0 = group count
[shader("compute")] [numthreads(64, 1, 1)]
void cost_prior(uint3 t : SV_DispatchThreadID) {
    uint gg = t.x;
    if (gg >= u0) return;
    Real cost;
    Real go[12], Ho[144];
    priorEval(gg, cost, go, Ho);
    if (isfinite(cost) && cost > Real(0.0)) atomicAccum(u_cost_a, 0, cost);
}

// u0 = group count, u1 = CG path, u2 = group-block base, u3 = exclusive groups,
// f0 = lambda. One thread owns the whole block it writes, so no atomics.
[shader("compute")] [numthreads(64, 1, 1)]
void intr_prior(uint3 t : SV_DispatchThreadID) {
    uint gg = t.x;
    if (gg >= u0) return;
    uint4 gi = u_group_info[gg];
    uint n = min(gi.z, 12u);
    if (n == 0) return;
    Real cost;
    Real go[12], Ho[144];
    priorEval(gg, cost, go, Ho);

    for (uint i = 0; i < n; i++) {
        Ho[12 * i + i] = Ho[12 * i + i] * (Real(1.0) + Real(f0));
        if (isfinite(go[i])) u_g[gi.y + i] = u_g[gi.y + i] + go[i];
    }
    if (u1 == 0) {
        for (uint i = 0; i < n; i++)
            for (uint j = 0; j <= i; j++) {
                Real v = Ho[12 * i + j];
                if (isfinite(v)) u_S[pidx(gi.y + i, gi.y + j)] = u_S[pidx(gi.y + i, gi.y + j)] + v;
            }
        return;
    }
    uint cb = kCamBlk * u_prior_layout[kPriorL * gg + kPlImage];
    uint gb = kCamBlk * (u2 + gg);
    for (uint i = 0; i < n; i++)
        for (uint j = 0; j <= i; j++) {
            Real v = Ho[12 * i + j];
            if (!isfinite(v)) continue;
            u_cg_B[cb + pidx(6 + i, 6 + j)] = u_cg_B[cb + pidx(6 + i, 6 + j)] + v;
            uint m = u3 != 0 ? cb + pidx(6 + i, 6 + j) : gb + pidx(i, j);
            u_cg_M[m] = u_cg_M[m] + v;
        }
}
