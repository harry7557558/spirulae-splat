#pragma once
// Emulated double precision ("double-float" / compensated float pair).
// A value is represented as hi + lo with |lo| <= 0.5 ulp(hi), giving ~49 bits
// of significand using only fp32 hardware. Algorithms follow the QD library
// (Hida, Li, Bailey). All arithmetic used inside differentiable code carries
// custom backward derivatives so Slang autodiff works transparently.
//
// Requires fp32 fma and strict fp evaluation (do not compile with fast-math).

public struct DF : IDifferentiable {
    public typealias Differential = DF;

    public float hi;
    public float lo;

    [Differentiable] public __init(float h) { hi = h; lo = 0.0f; }
    public __init(float h, float l) { hi = h; lo = l; }

    [Differentiable] public static DF dzero() { return DF(0.0f, 0.0f); }
    [Differentiable] public static Differential dadd(Differential a, Differential b) { return df_add(a, b); }
};

// ---------------- error-free transforms ----------------

DF df_two_sum(float a, float b) {
    float s = a + b;
    float bb = s - a;
    float err = (a - (s - bb)) + (b - bb);
    return DF(s, err);
}

// requires |a| >= |b|
DF df_quick_two_sum(float a, float b) {
    float s = a + b;
    return DF(s, b - (s - a));
}

DF df_two_prod(float a, float b) {
    float p = a * b;
    return DF(p, fma(a, b, -p));
}

// ---------------- core arithmetic ----------------

[BackwardDerivative(df_add_bwd)]
public DF df_add(DF a, DF b) {
    DF s = df_two_sum(a.hi, b.hi);
    DF t = df_two_sum(a.lo, b.lo);
    float l = s.lo + t.hi;
    DF r = df_quick_two_sum(s.hi, l);
    l = r.lo + t.lo;
    return df_quick_two_sum(r.hi, l);
}
void df_add_bwd(inout DifferentialPair<DF> a, inout DifferentialPair<DF> b, DF g) {
    a = diffPair(a.p, g);
    b = diffPair(b.p, g);
}

[BackwardDerivative(df_neg_bwd)]
public DF df_neg(DF a) { return DF(-a.hi, -a.lo); }
void df_neg_bwd(inout DifferentialPair<DF> a, DF g) { a = diffPair(a.p, df_neg(g)); }

[BackwardDerivative(df_sub_bwd)]
public DF df_sub(DF a, DF b) { return df_add(a, df_neg(b)); }
void df_sub_bwd(inout DifferentialPair<DF> a, inout DifferentialPair<DF> b, DF g) {
    a = diffPair(a.p, g);
    b = diffPair(b.p, df_neg(g));
}

[BackwardDerivative(df_mul_bwd)]
public DF df_mul(DF a, DF b) {
    DF p = df_two_prod(a.hi, b.hi);
    float l = p.lo + (a.hi * b.lo + a.lo * b.hi);
    return df_quick_two_sum(p.hi, l);
}
void df_mul_bwd(inout DifferentialPair<DF> a, inout DifferentialPair<DF> b, DF g) {
    a = diffPair(a.p, df_mul(g, b.p));
    b = diffPair(b.p, df_mul(g, a.p));
}

DF df_mul_f(DF a, float b) {
    DF p = df_two_prod(a.hi, b);
    float l = p.lo + a.lo * b;
    return df_quick_two_sum(p.hi, l);
}

[BackwardDerivative(df_div_bwd)]
public DF df_div(DF a, DF b) {
    float q1 = a.hi / b.hi;
    DF r = df_sub(a, df_mul_f(b, q1));
    float q2 = r.hi / b.hi;
    r = df_sub(r, df_mul_f(b, q2));
    float q3 = r.hi / b.hi;
    DF q = df_quick_two_sum(q1, q2);
    return df_add(q, DF(q3));
}
void df_div_bwd(inout DifferentialPair<DF> a, inout DifferentialPair<DF> b, DF g) {
    DF inv_b = df_div(DF(1.0f), b.p);
    DF ga = df_mul(g, inv_b);
    a = diffPair(a.p, ga);
    b = diffPair(b.p, df_neg(df_mul(ga, df_mul(a.p, inv_b))));
}

// ---------------- elementary functions ----------------

[BackwardDerivative(df_sqrt_bwd)]
public DF df_sqrt(DF a) {
    if (a.hi <= 0.0f)
        return DF(0.0f, 0.0f);
    float x = 1.0f / sqrt(a.hi);
    float ax = a.hi * x;
    DF ax2 = df_two_prod(ax, ax);
    float err = df_sub(a, ax2).hi * (x * 0.5f);
    return df_quick_two_sum(ax, err);
}
void df_sqrt_bwd(inout DifferentialPair<DF> a, DF g) {
    DF s = df_sqrt(a.p);
    a = diffPair(a.p, df_div(df_mul_f(g, 0.5f), s));
}

// Polynomials are near-minimax (Chebyshev) fits on the reduced range rather
// than Taylor series -- same interval accuracy with several fewer terms.
// The high-order terms contribute < ~6e-9 of the result, far below one df
// ulp after fp32 rounding, so the polynomial tail is evaluated in plain
// fp32 and only the head terms use full DF Horner steps.

static const float kDfInvLn2 = 1.4426950216293335f;
static const DF kDfLn2 = { 0.6931471824645996f, -1.9046542121259336e-09f };
// exp(r), |r| <= ln2/2: degree 10, max rel err 2.9e-16
static const DF kDfExpC[8] = {
    { 1.0f, -0.0f },
    { 1.0f, 6.7423483788263405e-15f },
    { 0.5f, 5.616924818567114e-16f },
    { 0.1666666716337204f, -4.96817653683479e-09f },
    { 0.0416666679084301f, -1.2418569506067456e-09f },
    { 0.008333333767950535f, -3.8228276100227276e-10f },
    { 0.0013888889225199819f, -2.9271123719309955e-11f },
    { 0.00019841169705614448f, 5.648256168333532e-12f },
};
static const float kDfExpTail[3] = { 2.480150396877434e-05f, 2.764018063317053e-06f,
                                     2.762635631370358e-07f };

[BackwardDerivative(df_exp_bwd)]
public DF df_exp(DF a) {
    // exp(a) = 2^k * exp(r),  r = a - k ln2,  |r| <= ln2/2
    float kf = round(a.hi * kDfInvLn2);
    kf = clamp(kf, -126.0f, 126.0f);
    DF r = df_sub(a, df_mul_f(kDfLn2, kf));
    float rf = r.hi;
    float t = (kDfExpTail[2] * rf + kDfExpTail[1]) * rf + kDfExpTail[0];
    DF p = DF(t);
    [ForceUnroll]
    for (int i = 7; i >= 0; i--)
        p = df_add(df_mul(p, r), kDfExpC[i]);
    float scale = exp2(kf);
    return DF(p.hi * scale, p.lo * scale);
}
void df_exp_bwd(inout DifferentialPair<DF> a, DF g) {
    a = diffPair(a.p, df_mul(g, df_exp(a.p)));
}

static const float kDfTwoOverPi = 0.6366197466850281f;
static const DF kDfPiO2 = { 1.5707963705062866f, -4.371138828673793e-08f };
// sin(r)/r and cos(r) as polynomials in s = r^2 <= (pi/4)^2: degree 6,
// max rel err 3.5e-18 / 6.7e-17
static const DF kDfSinC[5] = {
    { 1.0f, -3.147176764917495e-18f },
    { -0.1666666716337204f, 4.9670543234014986e-09f },
    { 0.008333333767950535f, -4.346301651914075e-10f },
    { -0.00019841270113829523f, 2.851792153119037e-12f },
    { 2.7557314297155244e-06f, -9.207551035636183e-14f },
};
static const float kDfSinTail[2] = { -2.5050717056274152e-08f, 1.5894736915544883e-10f };
static const DF kDfCosC[6] = {
    { 1.0f, -4.7157679218211415e-17f },
    { -0.5f, 7.492213094503417e-15f },
    { 0.0416666679084301f, -1.2419577588573816e-09f },
    { -0.0013888889225199819f, 3.552211649826553e-11f },
    { 2.4801578547339886e-05f, -8.776900342307489e-15f },
    { -2.755523382802494e-07f, -5.682247999459365e-16f },
};
static const float kDfCosTail[1] = { 2.063045512201711e-09f };

DF df_sin_poly(DF s, DF r) {
    float sf = s.hi;
    DF p = DF(kDfSinTail[1] * sf + kDfSinTail[0]);
    [ForceUnroll]
    for (int i = 4; i >= 0; i--)
        p = df_add(df_mul(p, s), kDfSinC[i]);
    return df_mul(p, r);
}
DF df_cos_poly(DF s) {
    DF p = DF(kDfCosTail[0]);
    [ForceUnroll]
    for (int i = 5; i >= 0; i--)
        p = df_add(df_mul(p, s), kDfCosC[i]);
    return p;
}

void df_sincos_impl(DF a, out DF s, out DF c) {
    float kf = round(a.hi * kDfTwoOverPi);
    DF r = df_sub(a, df_mul_f(kDfPiO2, kf));
    int j = int(kf) & 3;
    DF s2 = df_mul(r, r);
    DF sr = df_sin_poly(s2, r);
    DF cr = df_cos_poly(s2);
    if (j == 0)      { s = sr;          c = cr; }
    else if (j == 1) { s = cr;          c = df_neg(sr); }
    else if (j == 2) { s = df_neg(sr);  c = df_neg(cr); }
    else             { s = df_neg(cr);  c = sr; }
}

[BackwardDerivative(df_sin_bwd)]
public DF df_sin(DF a) { DF s, c; df_sincos_impl(a, s, c); return s; }
void df_sin_bwd(inout DifferentialPair<DF> a, DF g) {
    a = diffPair(a.p, df_mul(g, df_cos(a.p)));
}

[BackwardDerivative(df_cos_bwd)]
public DF df_cos(DF a) { DF s, c; df_sincos_impl(a, s, c); return c; }
void df_cos_bwd(inout DifferentialPair<DF> a, DF g) {
    a = diffPair(a.p, df_neg(df_mul(g, df_sin(a.p))));
}

// ---------------- operators & misc (overload builtins for DF) ----------------

[Differentiable] public DF operator +(DF a, DF b) { return df_add(a, b); }
[Differentiable] public DF operator -(DF a, DF b) { return df_sub(a, b); }
[Differentiable] public DF operator *(DF a, DF b) { return df_mul(a, b); }
[Differentiable] public DF operator /(DF a, DF b) { return df_div(a, b); }
[Differentiable] public __prefix DF operator -(DF a) { return df_neg(a); }

[Differentiable] public DF sqrt(DF a) { return df_sqrt(a); }
[Differentiable] public DF exp(DF a) { return df_exp(a); }
[Differentiable] public DF sin(DF a) { return df_sin(a); }
[Differentiable] public DF cos(DF a) { return df_cos(a); }

public bool operator <(DF a, DF b) { return a.hi < b.hi || (a.hi == b.hi && a.lo < b.lo); }
public bool operator >(DF a, DF b) { return b < a; }
public bool operator <=(DF a, DF b) { return !(b < a); }
public bool operator >=(DF a, DF b) { return !(a < b); }

public DF abs(DF a) { return a.hi < 0.0f || (a.hi == 0.0f && a.lo < 0.0f) ? df_neg(a) : a; }
public DF max(DF a, DF b) { return a < b ? b : a; }
public DF min(DF a, DF b) { return a < b ? a : b; }
public bool isfinite(DF a) { return isfinite(a.hi) && isfinite(a.lo); }
