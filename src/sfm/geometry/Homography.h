// Homography estimation (src/sfm/README.md).
//
//   estimateHomography  normalized DLT, >= 4 points (minimal sampler = 4,
//                       also the refit solver since DLT takes any N >= 4).
//   symmetricTransferSq symmetric transfer error (px^2), the RANSAC residual.
#pragma once

#include <cmath>
#include <vector>

#include "sfm/geometry/LinAlg.h"

namespace sfm {

// Hartley-normalized DLT homography H mapping image-1 points to image-2 points.
inline std::vector<Mat3> estimateHomography(const std::vector<Vec2>& p1,
                                            const std::vector<Vec2>& p2,
                                            const std::vector<int>& idx) {
    if (idx.size() < 4) return {};
    std::vector<Vec2> a, b;
    Mat3 T1 = hartleyNormalize(p1, idx, a);
    Mat3 T2 = hartleyNormalize(p2, idx, b);

    size_t n = idx.size();
    std::vector<double> A(2 * n * 9, 0.0);
    for (size_t i = 0; i < n; i++) {
        double x1 = a[i].x, y1 = a[i].y, x2 = b[i].x, y2 = b[i].y;
        double* r0 = &A[(2 * i) * 9];
        double* r1 = &A[(2 * i + 1) * 9];
        r0[0] = -x1; r0[1] = -y1; r0[2] = -1;
        r0[6] = x2 * x1; r0[7] = x2 * y1; r0[8] = x2;
        r1[3] = -x1; r1[4] = -y1; r1[5] = -1;
        r1[6] = y2 * x1; r1[7] = y2 * y1; r1[8] = y2;
    }
    auto nv = nullVectors9(A, (int)(2 * n), 1);
    if (nv.empty()) return {};
    Mat3 Hn;
    for (int i = 0; i < 9; i++) Hn[i] = nv[0][i];
    // denormalize: H = T2^-1 Hn T1
    Mat3 H = mul(mul(inverse3(T2), Hn), T1);
    if (std::fabs(H[8]) > 1e-12)
        for (double& v : H) v /= H[8];
    return {H};
}

// Homography on unit bearings (D45). A world plane maps rays to rays --
// b2 ~ H b1 with H = R + t n^T / d -- so the DLT rows are the same, written
// for homogeneous 3-vectors. No Hartley normalization: see estimateEpipolar*.
//
// The overall sign of H is fixed so that most of the fitted sample transfers
// *forwards* (positive dot). Without that, `angularTransferSq` would have to
// treat a ray and its opposite as equal, and a plane behind the camera would
// score as well as the real one.
inline std::vector<Mat3> estimateHomographyBearing(const std::vector<Vec3>& p1,
                                                   const std::vector<Vec3>& p2,
                                                   const std::vector<int>& idx) {
    size_t n = idx.size();
    if (n < 4) return {};
    std::vector<double> A(2 * n * 9, 0.0);
    for (size_t i = 0; i < n; i++) {
        const Vec3& a = p1[idx[i]];
        const Vec3& b = p2[idx[i]];
        double* r0 = &A[(2 * i) * 9];
        double* r1 = &A[(2 * i + 1) * 9];
        r0[0] = -b.z * a.x; r0[1] = -b.z * a.y; r0[2] = -b.z * a.z;
        r0[6] = b.x * a.x;  r0[7] = b.x * a.y;  r0[8] = b.x * a.z;
        r1[3] = -b.z * a.x; r1[4] = -b.z * a.y; r1[5] = -b.z * a.z;
        r1[6] = b.y * a.x;  r1[7] = b.y * a.y;  r1[8] = b.y * a.z;
    }
    auto nv = nullVectors9(A, (int)(2 * n), 1);
    if (nv.empty()) return {};
    Mat3 H;
    for (int i = 0; i < 9; i++) H[i] = nv[0][i];
    int forward = 0;
    for (size_t i = 0; i < n; i++)
        forward += p2[idx[i]].dot(mul(H, p1[idx[i]])) > 0 ? 1 : -1;
    if (forward < 0) for (double& v : H) v = -v;
    return {H};
}

// Symmetric angular transfer error on the sphere, in squared radians. A ray
// transferred to the opposite hemisphere lands past 90 deg and is rejected,
// which is the point of the sign convention above.
inline double angularTransferSq(const Mat3& H, const Mat3& Hinv, const Vec3& a, const Vec3& b) {
    auto ang2 = [](const Vec3& u, const Vec3& v) {
        Vec3 c = u.cross(v);
        double th = std::atan2(c.norm(), u.dot(v));
        return th * th;
    };
    Vec3 hb = mul(H, a), ha = mul(Hinv, b);
    double n1 = hb.norm(), n2 = ha.norm();
    if (n1 < 1e-30 || n2 < 1e-30) return 1e30;
    return ang2(hb * (1.0 / n1), b) + ang2(ha * (1.0 / n2), a);
}

inline Vec2 applyH(const Mat3& H, const Vec2& p) {
    Vec3 q = mul(H, Vec3{p.x, p.y, 1});
    if (std::fabs(q.z) < 1e-30) return {1e15, 1e15};
    return {q.x / q.z, q.y / q.z};
}

// Symmetric transfer error: |b - H a|^2 + |a - H^-1 b|^2.
inline double symmetricTransferSq(const Mat3& H, const Mat3& Hinv, const Vec2& a, const Vec2& b) {
    Vec2 hb = applyH(H, a);
    Vec2 ha = applyH(Hinv, b);
    double d1 = (hb.x - b.x) * (hb.x - b.x) + (hb.y - b.y) * (hb.y - b.y);
    double d2 = (ha.x - a.x) * (ha.x - a.x) + (ha.y - a.y) * (ha.y - a.y);
    return d1 + d2;
}

}  // namespace sfm
