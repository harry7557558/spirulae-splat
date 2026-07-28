#pragma once
// Small differentiable vector/matrix helpers over `Real`.
// Custom RVec2/RVec3 are used instead of vector<Real,N> so the emulated DF
// scalar works with the exact same code paths as float/double.

#include "real.slang"

struct RVec2 : IDifferentiable {
    typealias Differential = RVec2;
    Real x; Real y;
    [Differentiable] static RVec2 dzero() { RVec2 v = { Real(0.0), Real(0.0) }; return v; }
    [Differentiable] static RVec2 dadd(RVec2 a, RVec2 b) { RVec2 v = { a.x + b.x, a.y + b.y }; return v; }
};
struct RVec3 : IDifferentiable {
    typealias Differential = RVec3;
    Real x; Real y; Real z;
    [Differentiable] static RVec3 dzero() { RVec3 v = { Real(0.0), Real(0.0), Real(0.0) }; return v; }
    [Differentiable] static RVec3 dadd(RVec3 a, RVec3 b) { RVec3 v = { a.x + b.x, a.y + b.y, a.z + b.z }; return v; }
};

[Differentiable] RVec2 rvec2(Real x, Real y) { RVec2 v = { x, y }; return v; }
[Differentiable] RVec3 rvec3(Real x, Real y, Real z) { RVec3 v = { x, y, z }; return v; }

[Differentiable] RVec3 operator +(RVec3 a, RVec3 b) { return rvec3(a.x + b.x, a.y + b.y, a.z + b.z); }
[Differentiable] RVec3 operator -(RVec3 a, RVec3 b) { return rvec3(a.x - b.x, a.y - b.y, a.z - b.z); }
[Differentiable] RVec3 operator *(RVec3 a, Real s) { return rvec3(a.x * s, a.y * s, a.z * s); }
[Differentiable] RVec2 operator -(RVec2 a, RVec2 b) { return rvec2(a.x - b.x, a.y - b.y); }

// Fused differentiable sincos: one range reduction + polynomial pair per
// call (rsin/rcos each compute both and discard one), with a backward that
// also evaluates it only once. Returns (sin, cos) as an RVec2.
[BackwardDerivative(rsincos_bwd)]
RVec2 rsincos(Real x) {
    Real s, c;
    rsincos_raw(x, s, c);
    return rvec2(s, c);
}
void rsincos_bwd(inout DifferentialPair<Real> x, RVec2 g) {
    Real s, c;
    rsincos_raw(x.p, s, c);
    x = diffPair(x.p, g.x * c - g.y * s);
}

[Differentiable] Real dot3(RVec3 a, RVec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
[Differentiable] RVec3 cross3(RVec3 a, RVec3 b) {
    return rvec3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}
[Differentiable] Real length3(RVec3 a) { return sqrt(dot3(a, a)); }

// Inverse of a symmetric 3x3 given as row-major Real[9] (adjugate formula).
void inverse3(Real m[9], out Real r[9]) {
    Real c00 = m[4] * m[8] - m[5] * m[7];
    Real c01 = m[5] * m[6] - m[3] * m[8];
    Real c02 = m[3] * m[7] - m[4] * m[6];
    Real det = m[0] * c00 + m[1] * c01 + m[2] * c02;
    Real inv_det = Real(1.0) / det;
    r[0] = c00 * inv_det;
    r[1] = (m[2] * m[7] - m[1] * m[8]) * inv_det;
    r[2] = (m[1] * m[5] - m[2] * m[4]) * inv_det;
    r[3] = c01 * inv_det;
    r[4] = (m[0] * m[8] - m[2] * m[6]) * inv_det;
    r[5] = (m[2] * m[3] - m[0] * m[5]) * inv_det;
    r[6] = c02 * inv_det;
    r[7] = (m[1] * m[6] - m[0] * m[7]) * inv_det;
    r[8] = (m[0] * m[4] - m[1] * m[3]) * inv_det;
}
