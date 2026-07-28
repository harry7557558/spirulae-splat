#pragma once
// Double-precision transcendentals for the SPIR-V target.
//
// GLSL.std.450 only defines Exp/Log/Sin/Cos for 16/32-bit floats (Sqrt and
// Fma are fine on doubles), and slangc would otherwise emit invalid SPIR-V.
// These use standard argument reduction + near-minimax (Chebyshev) fits on
// the reduced range (fewer terms than Taylor at the same interval accuracy;
// a couple of ulps, plenty for BA). All carry custom derivatives so they can
// be used inside [Differentiable] code.

static const double kDInvLn2 = 1.442695040888963387;
static const double kDLn2Hi = 0.69314718048553913832;   // 32 high bits
static const double kDLn2Lo = 7.4406171100123968474e-11;
// exp(r), |r| <= ln2/2: degree 11, max rel err 4.2e-18
static const double kDExpCoef[12] = {
    0.99999999999999999693, 0.99999999999999999976, 0.50000000000000183855,
    0.16666666666666680806, 0.041666666666488095495, 0.0083333333333196006109,
    0.0013888888952314774652, 0.00019841269890047113707, 0.000024801485482328492419,
    2.7557240918578969823e-6, 2.7632639639041029749e-7, 2.5110037605963777712e-8,
};

double dbl_ldexp(double x, int k) {
    // x * 2^k via exponent bits (k clamped well inside range by callers)
    uint hi = uint((1023 + k) << 20);
    return x * asdouble(0u, hi);
}

[BackwardDerivative(dbl_exp_bwd)]
double dbl_exp(double x) {
    if (x > 709.0) return 1.0 / 0.0;
    if (x < -709.0) return 0.0;
    double kf = round(x * kDInvLn2);
    int k = int(kf);
    double r = (x - kf * kDLn2Hi) - kf * kDLn2Lo;
    double p = kDExpCoef[11];
    [ForceUnroll]
    for (int i = 10; i >= 0; i--)
        p = p * r + kDExpCoef[i];
    return dbl_ldexp(p, k);
}
void dbl_exp_bwd(inout DifferentialPair<double> x, double g) {
    x = diffPair(x.p, g * dbl_exp(x.p));
}

[BackwardDerivative(dbl_log_bwd)]
double dbl_log(double x) {
    // float log as seed, one refinement step: log(x) = y + (x e^-y - 1) - (..)^2/2
    double y = double(log(float(x)));
    double u = x * dbl_exp(-y) - 1.0;
    return y + u - 0.5 * u * u;
}
void dbl_log_bwd(inout DifferentialPair<double> x, double g) {
    x = diffPair(x.p, g / x.p);
}

static const double kDTwoOverPi = 0.63661977236758138243;
static const double kDPiO2_1 = 1.5707963267341256142;
static const double kDPiO2_2 = 6.0771005035934605454e-11;
static const double kDPiO2_3 = 2.9127320560933560353e-20;
// sin(r)/r and cos(r) as polynomials in s = r^2 <= (pi/4)^2: degree 6 / 7,
// max rel err 3.5e-18 / 4.3e-20
static const double kDSinCoef[7] = {
    0.99999999999999999685, -0.16666666666666616666, 0.0083333333333203624567,
    -0.00019841269828650300013, 2.7557313376400129128e-6, -2.5050716974102745028e-8,
    1.5894736651849095259e-10,
};
static const double kDCosCoef[8] = {
    0.99999999999999999997, -0.49999999999999999371, 0.041666666666666452389,
    -0.0013888888888861094596, 0.000024801587283881153004, -2.7557313097790086271e-7,
    2.0875582380663953044e-9, -1.1353379638297574126e-11,
};

double dbl_sin_taylor(double r) {
    double s = r * r;
    double p = kDSinCoef[6];
    [ForceUnroll]
    for (int i = 5; i >= 0; i--)
        p = p * s + kDSinCoef[i];
    return p * r;
}
double dbl_cos_taylor(double r) {
    double s = r * r;
    double p = kDCosCoef[7];
    [ForceUnroll]
    for (int i = 6; i >= 0; i--)
        p = p * s + kDCosCoef[i];
    return p;
}

void dbl_sincos(double x, out double sn, out double cs) {
    double kf = round(x * kDTwoOverPi);
    double r = ((x - kf * kDPiO2_1) - kf * kDPiO2_2) - kf * kDPiO2_3;
    int j = int(kf) & 3;
    double sr = dbl_sin_taylor(r);
    double cr = dbl_cos_taylor(r);
    if (j == 0)      { sn = sr;  cs = cr; }
    else if (j == 1) { sn = cr;  cs = -sr; }
    else if (j == 2) { sn = -sr; cs = -cr; }
    else             { sn = -cr; cs = sr; }
}

[BackwardDerivative(dbl_sin_bwd)]
double dbl_sin(double x) { double s, c; dbl_sincos(x, s, c); return s; }
void dbl_sin_bwd(inout DifferentialPair<double> x, double g) {
    x = diffPair(x.p, g * dbl_cos(x.p));
}

[BackwardDerivative(dbl_cos_bwd)]
double dbl_cos(double x) { double s, c; dbl_sincos(x, s, c); return c; }
void dbl_cos_bwd(inout DifferentialPair<double> x, double g) {
    x = diffPair(x.p, -g * dbl_sin(x.p));
}
