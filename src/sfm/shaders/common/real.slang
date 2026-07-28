#pragma once
// Scalar type configuration. Exactly one of REAL_FLOAT / REAL_DOUBLE / REAL_DF
// is defined on the slangc command line; everything downstream is written
// against `Real` and the small accumulation API below.
//
//   Real        - the arithmetic scalar type
//   AWord       - element type of buffers that receive atomic accumulation
//                 (bound as a second view of the same VkBuffer)
//   atomicAccum - lock-free `buf[i] += v`
//   loadAccum / storeAccum - non-atomic access through the AWord view
//   rexp/rlog/rsin/rcos - differentiable transcendentals (SPIR-V only has
//   them natively for 32-bit floats; sqrt works for every config)

#include "df.slang"
#include "dmath.slang"

#if defined(REAL_FLOAT)

typealias Real = float;
typealias AWord = float;

void atomicAccum(RWStructuredBuffer<AWord> buf, uint i, Real v) {
    __atomic_add(buf[i], v, MemoryOrder.Relaxed);
}
Real loadAccum(RWStructuredBuffer<AWord> buf, uint i) { return buf[i]; }
void storeAccum(RWStructuredBuffer<AWord> buf, uint i, Real v) { buf[i] = v; }

[Differentiable] Real rexp(Real x) { return exp(x); }
[Differentiable] Real rlog(Real x) { return log(x); }
[Differentiable] Real rsin(Real x) { return sin(x); }
[Differentiable] Real rcos(Real x) { return cos(x); }
void rsincos_raw(Real x, out Real s, out Real c) { s = sin(x); c = cos(x); }
float rToFloat(Real x) { return x; }

#elif defined(REAL_DOUBLE)

typealias Real = double;
typealias AWord = double;

void atomicAccum(RWStructuredBuffer<AWord> buf, uint i, Real v) {
    __atomic_add(buf[i], v, MemoryOrder.Relaxed);
}
Real loadAccum(RWStructuredBuffer<AWord> buf, uint i) { return buf[i]; }
void storeAccum(RWStructuredBuffer<AWord> buf, uint i, Real v) { buf[i] = v; }

[Differentiable] Real rexp(Real x) { return dbl_exp(x); }
[Differentiable] Real rlog(Real x) { return dbl_log(x); }
[Differentiable] Real rsin(Real x) { return dbl_sin(x); }
[Differentiable] Real rcos(Real x) { return dbl_cos(x); }
void rsincos_raw(Real x, out Real s, out Real c) { dbl_sincos(x, s, c); }
float rToFloat(Real x) { return float(x); }

#elif defined(REAL_DF)

typealias Real = DF;
typealias AWord = uint64_t;

// Memory layout of DF is {hi, lo} = 8 bytes; the packed uint64 view puts hi in
// the low word so both views of a buffer see identical bytes (little-endian).
uint64_t packReal(Real v) { return uint64_t(asuint(v.hi)) | (uint64_t(asuint(v.lo)) << 32); }
Real unpackReal(uint64_t w) { return DF(asfloat(uint(w)), asfloat(uint(w >> 32))); }

void atomicAccum(RWStructuredBuffer<AWord> buf, uint i, Real v) {
    uint64_t old_w = buf[i];
    for (;;) {
        uint64_t new_w = packReal(unpackReal(old_w) + v);
        uint64_t got = __atomic_compare_exchange(buf[i], old_w, new_w,
            MemoryOrder.Relaxed, MemoryOrder.Relaxed);
        if (got == old_w) break;
        old_w = got;
    }
}
Real loadAccum(RWStructuredBuffer<AWord> buf, uint i) { return unpackReal(buf[i]); }
void storeAccum(RWStructuredBuffer<AWord> buf, uint i, Real v) { buf[i] = packReal(v); }

[Differentiable] Real rexp(Real x) { return df_exp(x); }
[Differentiable] Real rsin(Real x) { return df_sin(x); }
[Differentiable] Real rcos(Real x) { return df_cos(x); }
void rsincos_raw(Real x, out Real s, out Real c) { df_sincos_impl(x, s, c); }
[BackwardDerivative(df_rlog_bwd)]
Real rlog(Real x) {
    // float log as seed + one Newton step through df_exp
    DF y = DF(log(x.hi));
    DF u = df_mul(x, df_exp(df_neg(y))) - DF(1.0f);
    return y + u - DF(0.5f) * (u * u);
}
void df_rlog_bwd(inout DifferentialPair<DF> x, DF g) {
    x = diffPair(x.p, df_div(g, x.p));
}
float rToFloat(Real x) { return x.hi; }

#else
#error "define one of REAL_FLOAT / REAL_DOUBLE / REAL_DF"
#endif

// atan2(r, z) for r >= 0, returning the angle from the +z axis in [0, pi].
// Generic across all Real configs (SPIR-V has no fp64 inverse-trig): a fp32
// atan2 gives a ~1e-7 seed, then Newton on g(t) = z*sin t - r*cos t = 0 (whose
// root is atan2(r,z), with g' = |(r,z)| != 0 for r>0) reaches full Real
// precision in a few steps, reusing the config's own rsincos. Used by the
// fisheye camera model. A custom backward derivative supplies the exact
// dtheta/dr = z/(r^2+z^2), dtheta/dz = -r/(r^2+z^2), so autodiff never
// differentiates through the Newton loop.
[BackwardDerivative(rAtan2_bwd)]
Real rAtan2(Real r, Real z) {
    Real theta = Real(atan2(rToFloat(r), rToFloat(z)));
    [ForceUnroll]
    for (int i = 0; i < 3; i++) {
        Real s, c;
        rsincos_raw(theta, s, c);
        theta = theta - (z * s - r * c) / (z * c + r * s);
    }
    return theta;
}
void rAtan2_bwd(inout DifferentialPair<Real> r, inout DifferentialPair<Real> z, Real g) {
    Real d = r.p * r.p + z.p * z.p;
    r = diffPair(r.p, g * z.p / d);
    z = diffPair(z.p, -(g * r.p) / d);
}
