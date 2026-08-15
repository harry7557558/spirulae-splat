#include "data/CameraMath.h"

#include "core/CameraModel.h"

#include <cmath>

namespace camhost {
namespace {

enum { M_PINHOLE = 0, M_FISHEYE = 1, M_EQUISOLID = 2, M_EQUIRECT = 3 };

constexpr double kPi = 3.14159265358979323846;   // MSVC has no M_PI by default

void distort_jac(double qx, double qy, int tier, const float* d, double J[4]) {
    const double e = 1e-5;
    double f[2], fx[2], fy[2];
    distort_lens(qx, qy, tier, d, f);
    distort_lens(qx + e, qy, tier, d, fx);
    distort_lens(qx, qy + e, tier, d, fy);
    J[0] = (fx[0]-f[0])/e; J[1] = (fx[1]-f[1])/e;   // j00 j10
    J[2] = (fy[0]-f[0])/e; J[3] = (fy[1]-f[1])/e;   // j01 j11
}

bool norm3(double x, double y, double z, double out[3]) {
    double l = std::sqrt(x*x + y*y + z*z);
    if (l < 1e-12) l = 1;
    out[0] = x/l; out[1] = y/l; out[2] = z/l;
    return true;
}

}  // namespace

bool has_distortion(int tier, const float* d) {
    if (tier == (int)CameraDistortionType::None) return false;
    for (int i = 0; i < kCameraDistortionParams; i++) if (d[i] != 0.0f) return true;
    return false;
}

// Host mirror of the tier structs in shaders/projection_utils.slang; the slot
// order differs per tier, which is why the tier travels with the coefficients.
void distort_lens(double u, double v, int tier, const float* d, double out[2]) {
    double r2 = u*u + v*v;
    if (tier == (int)CameraDistortionType::None) {
        out[0] = u; out[1] = v;
        return;
    }
    if (tier == (int)CameraDistortionType::OpenCV) {
        double radial = 1 + r2*(d[0] + r2*d[1]);
        out[0] = u*radial + 2*d[2]*u*v + d[3]*(r2 + 2*u*u);
        out[1] = v*radial + 2*d[3]*u*v + d[2]*(r2 + 2*v*v);
        return;
    }
    if (tier == (int)CameraDistortionType::Rational) {
        double radial = (1 + r2*(d[0] + r2*(d[1] + r2*d[2])))
                      / (1 + r2*(d[3] + r2*(d[4] + r2*d[5])));
        out[0] = u*radial + 2*d[6]*u*v + d[7]*(r2 + 2*u*u);
        out[1] = v*radial + 2*d[7]*u*v + d[6]*(r2 + 2*v*v);
        return;
    }
    double radial = 1 + r2*(d[0] + r2*(d[1] + r2*(d[2] + r2*d[3])));
    out[0] = u*radial + 2*d[4]*u*v + d[5]*(r2 + 2*u*u) + d[6]*r2;
    out[1] = v*radial + 2*d[5]*u*v + d[4]*(r2 + 2*v*v) + d[7]*r2;
}

// Newton solve distort(q) = (u,v), with undistort_point_0's failure
// conditions (well-posed forward Jacobian + re-distorts within 0.01).
bool undistort_lens(double u, double v, int tier, const float* d, double out[2]) {
    double qx = u, qy = v;
    for (int it = 0; it < 8; it++) {
        double f[2], J[4];
        distort_lens(qx, qy, tier, d, f);
        distort_jac(qx, qy, tier, d, J);
        double det = J[0]*J[3] - J[2]*J[1];
        if (std::fabs(det) < 1e-12) break;
        double rx = f[0] - u, ry = f[1] - v, inv = 1/det;
        qx -= ( J[3]*rx - J[2]*ry)*inv;
        qy -= (-J[1]*rx + J[0]*ry)*inv;
    }
    if (!std::isfinite(qx) || !std::isfinite(qy)) return false;
    double J[4], f[2];
    distort_jac(qx, qy, tier, d, J);
    double det = J[0]*J[3] - J[2]*J[1];
    if (std::fmin(det, std::fmin(J[0], J[3])) <= 0) return false;   // folded
    distort_lens(qx, qy, tier, d, f);
    if (std::hypot(f[0]-u, f[1]-v) >= 0.01) return false;
    out[0] = qx; out[1] = qy;
    return true;
}

void distort_point(double u, double v, int model, int tier, const float* d,
                   double out[2]) {
    if (model == M_EQUIRECT) {
        out[0] = u; out[1] = v;
        return;
    }
    if (model == M_FISHEYE || model == M_EQUISOLID) {
        const double r = std::hypot(u, v);
        const double theta = std::atan(r);
        double k;
        if (model == M_FISHEYE)
            k = (r < 1e-3) ? (1.0 - theta*theta/6.0) : (theta / r);
        else
            k = (r < 1e-3) ? (1.0 - theta*theta/24.0)
                           : (2.0 * std::sin(0.5*theta) / r);
        u *= k; v *= k;
    }
    distort_lens(u, v, tier, d, out);
}

bool undistort_point(double u, double v, int model, int tier, const float* d,
                     double out[2]) {
    double ray[3];
    if (!generate_ray(u, v, model, tier, d, ray)) return false;
    if (ray[2] <= 1e-12) return false;
    out[0] = ray[0] / ray[2];
    out[1] = ray[1] / ray[2];
    return true;
}

bool generate_ray(double u, double v, int model, int tier, const float* d,
                  double out[3]) {
    if (model == M_EQUIRECT) {
        if (std::fabs(u) > kPi || std::fabs(v) > kPi/2) return false;
        double cl = std::cos(v);
        out[0] = cl*std::sin(u); out[1] = std::sin(v); out[2] = cl*std::cos(u);
        return true;
    }
    double uu = u, vv = v;
    if (has_distortion(tier, d)) {
        double q[2];
        if (!undistort_lens(u, v, tier, d, q)) return false;
        uu = q[0]; vv = q[1];
    }
    double r = std::hypot(uu, vv);
    if (model == M_FISHEYE) {                        // equidistant: theta = r
        if (r >= kPi) return false;
        double s = r < 1e-3 ? (1 - r*r/6) : (std::sin(r)/r);
        return norm3(uu*s, vv*s, std::cos(r), out);
    }
    if (model == M_EQUISOLID) {                      // r = 2 sin(theta/2)
        if (r >= 2.0) return false;
        double k = std::sqrt(std::fmax(0.0, 1 - 0.25*r*r));
        return norm3(uu*k, vv*k, 1 - 0.5*r*r, out);
    }
    return norm3(uu, vv, 1.0, out);                  // pinhole
}

bool project_ray(const double ray[3], int model, int tier, const float* d,
                 double out[2]) {
    if (model == M_EQUIRECT) {
        out[0] = std::atan2(ray[0], ray[2]);
        out[1] = std::atan2(ray[1], std::hypot(ray[0], ray[2]));
        return true;
    }
    const double rz = ray[2];
    const double rxy = std::hypot(ray[0], ray[1]);
    if (model == M_PINHOLE) {
        if (rz <= 1e-12) return false;
        distort_point(ray[0]/rz, ray[1]/rz, model, tier, d, out);
        return true;
    }
    // Both wide models are angular, so theta comes straight off the ray. This
    // does NOT go through distort_point: that one recovers theta as atan(r),
    // which is the wrong branch past 90 degrees -- exactly the half of a wide
    // fisheye the split exists to keep.
    const double theta = std::atan2(rxy, rz);
    if (theta >= kPi) return false;
    const double r_ang = (model == M_FISHEYE) ? theta : 2.0 * std::sin(0.5 * theta);
    const double s = (rxy < 1e-12) ? 0.0 : r_ang / rxy;
    distort_lens(ray[0] * s, ray[1] * s, tier, d, out);
    return true;
}

double pinhole_coverage(int model, int width, int height, double fx, double fy) {
    if (model == M_EQUIRECT) return 0.0;
    if (model != M_FISHEYE && model != M_EQUISOLID) return 1.0;
    if (!(fx > 0) || !(fy > 0) || width <= 0 || height <= 0) return 1.0;

    // Both images share an axis and intrinsics, so at azimuth phi the frame
    // reaches radius rect(phi) in the wide one and tan(theta) = rect(phi) in
    // the pinhole -- the same expression, hence no root finding.
    const double U = 0.5 * width / fx, V = 0.5 * height / fy;
    const auto kept = [&](double phi) {
        const double c = std::cos(phi), s = std::sin(phi);
        const double rect = std::fmin(U / std::fmax(c, 1e-300),
                                      V / std::fmax(s, 1e-300));
        const double theta = std::atan(rect);
        const double r = (model == M_FISHEYE) ? theta : 2.0 * std::sin(0.5 * theta);
        return r * r;
    };
    // Simpson on each side of the corner, where the frame's radius has a kink.
    const double corner = std::atan2(V, U);
    double area = 0.0;
    for (int half = 0; half < 2; ++half) {
        const double a = half ? corner : 0.0, b = half ? 0.5 * kPi : corner;
        const int n = 128;                       // even, for Simpson
        const double h = (b - a) / n;
        double s = kept(a) + kept(b);
        for (int i = 1; i < n; ++i) s += kept(a + i * h) * (i & 1 ? 4 : 2);
        area += s * h / 3.0;
    }
    // Against 2 * (the quadrant's area): both are integrals of r^2/2 over a
    // quarter turn, so the leading half cancels.
    return std::fmin(1.0, area / (2.0 * U * V));
}

}  // namespace camhost
