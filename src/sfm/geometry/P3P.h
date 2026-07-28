// P3P minimal solver (3-point absolute pose), ported from PoseLib's lambdatwist
// implementation (Persson & Nordberg; BSD, the same solver COLMAP uses via
// PoseLib). Unlike DLT PnP this is robust to near-coplanar / elongated point
// configurations, which are common early in SfM, so it is the mapper's PnP
// minimal solver.
//
// Input: unit bearing vectors x[i] (normalized camera rays) and world points
// X[i]; output poses satisfy x_cam = R * X_world + t (world->camera).
#pragma once

#include <array>
#include <cmath>
#include <vector>

#include "sfm/geometry/Essential.h"  // Pose
#include "sfm/geometry/LinAlg.h"

namespace sfm {
namespace p3pdetail {

inline bool root2real(double b, double c, double& r1, double& r2) {
    double v = b * b - 4.0 * c;
    if (v < -1e-12) { r1 = r2 = -0.5 * b; return false; }
    if (v < 0.0) { r1 = -0.5 * b; r2 = -2; return true; }
    double y = std::sqrt(v);
    if (b < 0) { r1 = 0.5 * (-b + y); r2 = 0.5 * (-b - y); }
    else { r1 = 2.0 * c / (-b + y); r2 = 0.5 * (-b - y); }
    return true;
}

// One real root of x^3 + c2 x^2 + c1 x + c0 = 0; returns true if it is the only
// real root (used to short-circuit the two-solution branch).
inline bool solveCubicSingleReal(double c2, double c1, double c0, double& root) {
    double a = c1 - c2 * c2 / 3.0;
    double b = (2.0 * c2 * c2 * c2 - 9.0 * c2 * c1) / 27.0 + c0;
    double c = b * b / 4.0 + a * a * a / 27.0;
    if (c != 0) {
        if (c > 0) {
            c = std::sqrt(c);
            b *= -0.5;
            root = std::cbrt(b + c) + std::cbrt(b - c) - c2 / 3.0;
            return true;
        }
        c = 3.0 * b / (2.0 * a) * std::sqrt(-3.0 / a);
        root = 2.0 * std::sqrt(-a / 3.0) * std::cos(std::acos(c) / 3.0) - c2 / 3.0;
    } else {
        root = -c2 / 3.0 + (a != 0 ? (3.0 * b / a) : 0);
    }
    return false;
}

inline Mat3 colsToMat(const Vec3& c0, const Vec3& c1, const Vec3& c2) {
    return {c0.x, c1.x, c2.x, c0.y, c1.y, c2.y, c0.z, c1.z, c2.z};
}
inline Vec3 col(const Mat3& M, int c) { return {M[c], M[3 + c], M[6 + c]}; }
inline Vec3 row(const Mat3& M, int r) { return {M[3 * r], M[3 * r + 1], M[3 * r + 2]}; }

// pq extraction: two planes whose intersection encodes the depth ratios.
inline std::array<Vec3, 2> computePQ(Mat3 C) {
    auto A = [&](int r, int c) -> double& { return C[3 * r + c]; };
    Mat3 J;
    auto B = [&](int r, int c) -> double& { return J[3 * r + c]; };
    B(0, 0) = A(1, 2) * A(2, 1) - A(1, 1) * A(2, 2);
    B(1, 1) = A(0, 2) * A(2, 0) - A(0, 0) * A(2, 2);
    B(2, 2) = A(0, 1) * A(1, 0) - A(0, 0) * A(1, 1);
    B(0, 1) = A(0, 1) * A(2, 2) - A(0, 2) * A(2, 1);
    B(0, 2) = A(0, 2) * A(1, 1) - A(0, 1) * A(1, 2);
    B(1, 0) = B(0, 1);
    B(1, 2) = A(0, 0) * A(1, 2) - A(0, 2) * A(1, 0);
    B(2, 0) = B(0, 2);
    B(2, 1) = B(1, 2);

    Vec3 v;
    if (B(0, 0) > B(1, 1)) {
        if (B(0, 0) > B(2, 2)) v = col(J, 0) * (1.0 / std::sqrt(B(0, 0)));
        else v = col(J, 2) * (1.0 / std::sqrt(B(2, 2)));
    } else if (B(1, 1) > B(2, 2)) {
        v = col(J, 1) * (1.0 / std::sqrt(B(1, 1)));
    } else {
        v = col(J, 2) * (1.0 / std::sqrt(B(2, 2)));
    }
    A(0, 1) -= v.z; A(0, 2) += v.y; A(1, 2) -= v.x;
    A(1, 0) += v.z; A(2, 0) -= v.y; A(2, 1) += v.x;
    return {col(C, 0), row(C, 0)};
}

inline void refineLambda(double& l1, double& l2, double& l3, double a12, double a13, double a23,
                         double b12, double b13, double b23) {
    for (int it = 0; it < 5; it++) {
        double r1 = l1 * l1 - 2.0 * l1 * l2 * b12 + l2 * l2 - a12;
        double r2 = l1 * l1 - 2.0 * l1 * l3 * b13 + l3 * l3 - a13;
        double r3 = l2 * l2 - 2.0 * l2 * l3 * b23 + l3 * l3 - a23;
        if (std::fabs(r1) + std::fabs(r2) + std::fabs(r3) < 1e-10) return;
        double x11 = l1 - l2 * b12, x12 = l2 - l1 * b12;
        double x21 = l1 - l3 * b13, x23 = l3 - l1 * b13;
        double x32 = l2 - l3 * b23, x33 = l3 - l2 * b23;
        double detJ = 0.5 / (x11 * x23 * x32 + x12 * x21 * x33);
        l1 += (-x23 * x32 * r1 - x12 * x33 * r2 + x12 * x23 * r3) * detJ;
        l2 += (-x21 * x33 * r1 + x11 * x33 * r2 - x11 * x23 * r3) * detJ;
        l3 += (x21 * x32 * r1 - x11 * x32 * r2 - x12 * x21 * r3) * detJ;
    }
}

}  // namespace p3pdetail

inline std::vector<Pose> p3p(std::array<Vec3, 3> x, std::array<Vec3, 3> X) {
    using namespace p3pdetail;
    std::vector<Pose> out;

    Vec3 X01 = X[0] - X[1], X02 = X[0] - X[2], X12 = X[1] - X[2];
    double a01 = X01.dot(X01), a02 = X02.dot(X02), a12 = X12.dot(X12);

    // Order so X12 (=BC) is the largest baseline.
    if (a01 > a02) {
        if (a01 > a12) {
            std::swap(x[0], x[2]); std::swap(X[0], X[2]); std::swap(a01, a12);
            X01 = X12 * -1.0; X02 = X02 * -1.0;
        }
    } else if (a02 > a12) {
        std::swap(x[0], x[1]); std::swap(X[0], X[1]); std::swap(a02, a12);
        X01 = X01 * -1.0; X02 = X12;
    }

    double a12d = 1.0 / a12;
    double a = a01 * a12d, b = a02 * a12d;
    double m01 = x[0].dot(x[1]), m02 = x[0].dot(x[2]), m12 = x[1].dot(x[2]);

    double m12sq = -m12 * m12 + 1.0;
    double m02sq = -1.0 + m02 * m02;
    double m01sq = -1.0 + m01 * m01;
    double ab = a * b, bsq = b * b, asq = a * a;
    double m013 = -2.0 + 2.0 * m01 * m02 * m12;
    double bsqm12sq = bsq * m12sq, asqm12sq = asq * m12sq, abm12sq = 2.0 * ab * m12sq;

    double denom = bsqm12sq + b * m02sq;
    if (std::fabs(denom) < 1e-300) return out;
    double k3_inv = 1.0 / denom;
    double k2 = k3_inv * ((-1.0 + a) * m02sq + abm12sq + bsqm12sq + b * m013);
    double k1 = k3_inv * (asqm12sq + abm12sq + a * m013 + (-1.0 + b) * m01sq);
    double k0 = k3_inv * (asqm12sq + a * m01sq);

    double s;
    bool G = solveCubicSingleReal(k2, k1, k0, s);

    Mat3 C;
    auto CC = [&](int r, int c) -> double& { return C[3 * r + c]; };
    CC(0, 0) = -a + s * (1 - b);
    CC(0, 1) = -m02 * s;
    CC(0, 2) = a * m12 + b * m12 * s;
    CC(1, 0) = CC(0, 1);
    CC(1, 1) = s + 1;
    CC(1, 2) = -m01;
    CC(2, 0) = CC(0, 2);
    CC(2, 1) = CC(1, 2);
    CC(2, 2) = -a - b * s + 1;

    std::array<Vec3, 2> pq = computePQ(C);

    Mat3 XX = colsToMat(X01, X02, X01.cross(X02));
    XX = inverse3(XX);

    for (int i = 0; i < 2; i++) {
        double p0 = pq[i].x, p1 = pq[i].y, p2 = pq[i].z;
        bool switch12 = std::fabs(p0) <= std::fabs(p1);
        double d0, d1, d2;
        if (switch12) {
            double w0 = -p0 / p1, w1 = -p2 / p1;
            double ca = 1.0 / (w1 * w1 - b);
            double cb = 2.0 * (b * m12 - m02 * w1 + w0 * w1) * ca;
            double cc = (w0 * w0 - 2 * m02 * w0 - b + 1.0) * ca;
            double taus[2];
            if (!root2real(cb, cc, taus[0], taus[1])) continue;
            for (double tau : taus) {
                if (tau <= 0) continue;
                d2 = std::sqrt(a12 / (tau * (tau - 2.0 * m12) + 1.0));
                d1 = tau * d2;
                d0 = w0 * d2 + w1 * d1;
                if (d0 < 0) continue;
                refineLambda(d0, d1, d2, a01, a02, a12, m01, m02, m12);
                Vec3 v1 = x[0] * d0 - x[1] * d1, v2 = x[0] * d0 - x[2] * d2;
                Mat3 YY = colsToMat(v1, v2, v1.cross(v2));
                Mat3 R = mul(YY, XX);
                out.push_back({R, x[0] * d0 - mul(R, X[0])});
            }
        } else {
            double w0 = -p1 / p0, w1 = -p2 / p0;
            double ca = 1.0 / (-a * w1 * w1 + 2 * a * m12 * w1 - a + 1);
            double cb = 2 * (a * m12 * w0 - m01 - a * w0 * w1) * ca;
            double cc = (1 - a * w0 * w0) * ca;
            double taus[2];
            if (!root2real(cb, cc, taus[0], taus[1])) continue;
            for (double tau : taus) {
                if (tau <= 0) continue;
                d0 = std::sqrt(a01 / (tau * (tau - 2.0 * m01) + 1.0));
                d1 = tau * d0;
                d2 = w0 * d0 + w1 * d1;
                if (d2 < 0) continue;
                refineLambda(d0, d1, d2, a01, a02, a12, m01, m02, m12);
                Vec3 v1 = x[0] * d0 - x[1] * d1, v2 = x[0] * d0 - x[2] * d2;
                Mat3 YY = colsToMat(v1, v2, v1.cross(v2));
                Mat3 R = mul(YY, XX);
                out.push_back({R, x[0] * d0 - mul(R, X[0])});
            }
        }
        if (!out.empty() && G) break;
    }
    return out;
}

}  // namespace sfm
