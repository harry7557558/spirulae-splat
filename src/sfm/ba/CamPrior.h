// Deadzone-quadratic prior on a camera group's distortion, for bundle
// adjustment. Measured in image space -- the RMS fractional radial displacement
// over the radius the group's picture covers -- so a lens whose coefficients
// are large but cancel is not penalized and one whose coefficients are small
// but fold the image is. See sfm/ba/README.md, "Intrinsics prior".
#pragma once

#include <cmath>
#include <cstdint>

// The weight is per observation of the group, so the prior keeps its share of
// the objective as a model grows. `dead` is the deviation that costs nothing.
struct CamPrior {
    double w = 0, dead = 0;
    double rmax = 0;  // radius of the group's picture: r/f, or radians for a fisheye

    bool active() const { return w > 0 && rmax > 0; }
};

// Which BA parameters are the distortion, per model index of kModels[] in
// sfm/ba/Problem.h. The order is Camera.h packIntrinsics'; `rad` and `den`
// ascend in powers of r^2.
struct CamPriorLayout {
    uint8_t nrad, nden, ntan, nprism;
    uint8_t rad[4], den[3], tan[2], prism[2];
};

static constexpr CamPriorLayout kCamPriorLayout[] = {
    {0, 0, 0, 0, {}, {}, {}, {}},                    // snavely
    {0, 0, 0, 0, {}, {}, {}, {}},                    // snavely_f
    {2, 0, 0, 0, {1, 2}, {}, {}, {}},                // pinhole_radial
    {2, 0, 2, 0, {2, 3}, {}, {4, 5}, {}},            // opencv
    {0, 0, 0, 0, {}, {}, {}, {}},                    // simple_pinhole
    {0, 0, 0, 0, {}, {}, {}, {}},                    // pinhole
    {4, 0, 0, 0, {2, 3, 4, 5}, {}, {}, {}},          // opencv_fisheye
    {3, 3, 2, 0, {2, 3, 6}, {7, 8, 9}, {4, 5}, {}},  // full_opencv
    {4, 0, 2, 2, {2, 3, 6, 7}, {}, {4, 5}, {8, 9}},  // thin_prism_fisheye
    {0, 0, 0, 0, {}, {}, {}, {}},                    // equirect
};

// Six-point Gauss-Legendre in x = (r/rmax)^2, which is the area measure
// 2r dr / rmax^2 and integrates every polynomial model's term exactly.
static constexpr int kCamPriorQuad = 6;
static constexpr double kCamPriorX[kCamPriorQuad] = {
    0.03376524289842399, 0.16939530676686776, 0.38069040695840155,
    0.61930959304159845, 0.83060469323313224, 0.96623475710157601};
static constexpr double kCamPriorW[kCamPriorQuad] = {
    0.08566224618958517, 0.18038078652406930, 0.23395696728634553,
    0.23395696728634553, 0.18038078652406930, 0.08566224618958517};

// Mean square of the distortion's fractional displacement over the disc of
// radius `rmax`, with the gradient pieces `u` = J^T rho and `A` = J^T J when
// they are wanted. Squared because the deadzone works on the root of the sum.
inline double camPriorDist(uint32_t model, const double* p, double rmax, double* u,
                           double (*A)[12]) {
    const CamPriorLayout& L = kCamPriorLayout[model];
    if (!(rmax > 0) || !(L.nrad || L.ntan || L.nprism)) return 0;
    double d[12];
    double s2 = 0;
    const double R2 = rmax * rmax;
    for (int q = 0; q < kCamPriorQuad; q++) {
        const double r2 = R2 * kCamPriorX[q], w = kCamPriorW[q];
        double num = 1, den = 1, rp = 1;
        for (int i = 0; i < L.nrad; i++) { rp *= r2; num += p[L.rad[i]] * rp; }
        rp = 1;
        for (int i = 0; i < L.nden; i++) { rp *= r2; den += p[L.den[i]] * rp; }
        for (int i = 0; i < 12; i++) d[i] = 0;
        rp = 1;
        for (int i = 0; i < L.nrad; i++) { rp *= r2; d[L.rad[i]] = rp / den; }
        rp = 1;
        for (int i = 0; i < L.nden; i++) { rp *= r2; d[L.den[i]] = -num * rp / (den * den); }
        const double dev = num / den - 1.0;
        s2 += w * dev * dev;
        for (int i = 0; u && i < 12; i++) {
            if (d[i] == 0) continue;
            u[i] += w * d[i] * dev;
            for (int j = 0; j < 12; j++) A[i][j] += w * d[i] * d[j];
        }
    }
    // Tangential and thin-prism displacement, averaged over the azimuth in the
    // same fractional metric: 2.5 r^2 per tangential coefficient, 0.5 r^2 per
    // prism one, and both orthogonal to the radial terms.
    const double ct = 2.5 * R2, cp = 0.5 * R2;
    for (int i = 0; i < L.ntan; i++) {
        const uint8_t k = L.tan[i];
        s2 += ct * p[k] * p[k];
        if (!u) continue;
        u[k] += ct * p[k];
        A[k][k] += ct;
    }
    for (int i = 0; i < L.nprism; i++) {
        const uint8_t k = L.prism[i];
        s2 += cp * p[k] * p[k];
        if (!u) continue;
        u[k] += cp * p[k];
        A[k][k] += cp;
    }
    return s2;
}

// How far this lens is from an undistorted one, in the metric above: the RMS
// fractional radial displacement over the picture's radius.
inline double camPriorDeviation(uint32_t model, const double* p, double rmax) {
    return std::sqrt(camPriorDist(model, p, rmax, nullptr, nullptr));
}

// Cost, gradient and Gauss-Newton Hessian of one group's prior, accumulated
// into `g` and `H` over the first `nfree` intrinsics. Parameters past nfree are
// held constant and so contribute to the deviation but own no derivative.
inline double camPriorEval(const CamPrior& pr, uint32_t model, const double* p, uint32_t nfree,
                           double* g, double* H, uint32_t stride) {
    if (!pr.active()) return 0;
    double u[12] = {}, A[12][12] = {};
    const double s = std::sqrt(camPriorDist(model, p, pr.rmax, u, A));
    if (!(s > pr.dead)) return 0;
    const double e = s - pr.dead, a = e / s;
    const double rank1 = pr.w * pr.dead / (s * s * s);
    for (uint32_t i = 0; i < nfree && i < 12; i++) {
        g[i] += pr.w * a * u[i];
        for (uint32_t j = 0; j < nfree && j < 12; j++)
            H[i * stride + j] += pr.w * a * A[i][j] + rank1 * u[i] * u[j];
    }
    return 0.5 * pr.w * e * e;
}
