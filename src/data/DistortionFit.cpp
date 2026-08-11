#include "data/DistortionFit.h"

#include <algorithm>
#include <cmath>
#include <vector>

// Method, in one place:
//
//  1. Scan the source's field of view by bisecting theta along 128 azimuths of
//     the target's own axis, giving the valid domain as a star region.
//  2. Sample a tensor-product Chebyshev-Gauss grid over that region's bounding
//     box in the TARGET's undistorted normalized coordinates -- theta-space for
//     FISHEYE/EQUISOLID, tan-space for PINHOLE (projection_utils.slang). Equal
//     weights at those nodes integrate against 1/sqrt(1-x^2), which is what
//     makes the L2 solution track the L-infinity one.
//  3. The residual is bilinear -- linear in k1..k4/p1/p2/sx1/sy1 at fixed
//     intrinsics, linear in fx/fy/cx/cy at fixed coefficients -- so a linear
//     solve with c = 0 gives the intrinsics, then Levenberg-Marquardt over all
//     12 unknowns converges to machine precision. Alternating the two exact
//     blocks instead looks tempting and is not: fx and k1..k4 are strongly
//     coupled and it stalls near 2 px on a source the tier represents exactly.
//  4. A few Lawson rounds (w <- w * |r|) push that toward minimax; the best
//     max-error iterate is kept.
//  5. If the fitted distortion turns out not to be invertible over the domain,
//     the high-order coefficients are masked off and 1-4 repeated. A slightly
//     worse camera the kernels accept beats an exact one they reject.

namespace dsfit {
namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr int kNC = 8;   // ThinPrism coefficients
constexpr int kNU = 12;  // ... plus fx, cx, fy, cy

struct Model {
    double fx = 1.0, fy = 1.0, cx = 0.0, cy = 0.0;
    double c[kNC] = {};
};

struct Sample {
    double u = 0.0, v = 0.0;    // target undistorted normalized
    double us = 0.0, vs = 0.0;  // source pixel
    double w = 1.0;
    double gx[kNC] = {}, gy[kNC] = {};
};

// Past 175 degrees a lens is not something this fits, and tan(85 deg) = 11.4
// is already as far as a polynomial in tan-space can be pushed.
constexpr double kWidestTheta = 175.0 * kPi / 180.0;

double theta_cap(CameraModelType m) {
    return m == CameraModelType::PINHOLE ? 85.0 * kPi / 180.0 : kWidestTheta;
}

double rho_of_theta(CameraModelType m, double th) {
    switch (m) {
        case CameraModelType::PINHOLE:   return std::tan(th);
        case CameraModelType::FISHEYE:   return th;
        case CameraModelType::EQUISOLID: return 2.0 * std::sin(0.5 * th);
        default:                         return 0.0;
    }
}

bool unproject(CameraModelType m, double u, double v, double d[3]) {
    if (m == CameraModelType::PINHOLE) {
        d[0] = u; d[1] = v; d[2] = 1.0;
        return true;
    }
    double r = std::sqrt(u * u + v * v), th;
    if (m == CameraModelType::FISHEYE) {
        th = r;
    } else if (m == CameraModelType::EQUISOLID) {
        if (r > 2.0) return false;
        th = 2.0 * std::asin(0.5 * r);
    } else {
        return false;
    }
    if (th > kPi) return false;
    double s = r > 1e-12 ? std::sin(th) / r : 1.0;
    d[0] = u * s; d[1] = v * s; d[2] = std::cos(th);
    return true;
}

// Column basis of DistThinPrism::distort, matching the slot order
// k1 k2 k3 k4 p1 p2 sx1 sy1.
void distort_basis(double u, double v, double gx[kNC], double gy[kNC]) {
    double r2 = u * u + v * v, r4 = r2 * r2, r6 = r4 * r2, r8 = r4 * r4;
    gx[0] = u * r2; gx[1] = u * r4; gx[2] = u * r6; gx[3] = u * r8;
    gx[4] = 2.0 * u * v; gx[5] = r2 + 2.0 * u * u; gx[6] = r2; gx[7] = 0.0;
    gy[0] = v * r2; gy[1] = v * r4; gy[2] = v * r6; gy[3] = v * r8;
    gy[4] = r2 + 2.0 * v * v; gy[5] = 2.0 * u * v; gy[6] = 0.0; gy[7] = r2;
}

// Mirrors is_valid_distortion<DistThinPrism>: the distortion must stay locally
// orientation-preserving or the kernels drop the primitive.
bool distortion_invertible(double u, double v, const double c[kNC]) {
    double r2 = u * u + v * v, r4 = r2 * r2, r6 = r4 * r2;
    double R = 1.0 + r2 * (c[0] + r2 * (c[1] + r2 * (c[2] + r2 * c[3])));
    double Rp = c[0] + 2.0 * c[1] * r2 + 3.0 * c[2] * r4 + 4.0 * c[3] * r6;
    double p1 = c[4], p2 = c[5], sx = c[6], sy = c[7];
    double j00 = R + 2.0 * u * u * Rp + 2.0 * p1 * v + 6.0 * p2 * u + 2.0 * sx * u;
    double j01 = 2.0 * u * v * Rp + 2.0 * p1 * u + 2.0 * p2 * v + 2.0 * sx * v;
    double j10 = 2.0 * u * v * Rp + 2.0 * p2 * v + 2.0 * p1 * u + 2.0 * sy * u;
    double j11 = R + 2.0 * v * v * Rp + 2.0 * p2 * u + 6.0 * p1 * v + 2.0 * sy * v;
    return std::min(j00 * j11 - j01 * j10, std::min(j00, j11)) > 0.0;
}

// Cholesky of the lower triangle of A (row major, n <= kNU), then solve into b.
bool solve_spd(int n, double* A, double* b) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double s = A[i * n + j];
            for (int k = 0; k < j; k++) s -= A[i * n + k] * A[j * n + k];
            if (i == j) {
                if (!(s > 0.0)) return false;
                A[i * n + i] = std::sqrt(s);
            } else {
                A[i * n + j] = s / A[j * n + j];
            }
        }
    }
    for (int i = 0; i < n; i++) {
        double s = b[i];
        for (int k = 0; k < i; k++) s -= A[i * n + k] * b[k];
        b[i] = s / A[i * n + i];
    }
    for (int i = n - 1; i >= 0; i--) {
        double s = b[i];
        for (int k = i + 1; k < n; k++) s -= A[k * n + i] * b[k];
        b[i] = s / A[i * n + i];
    }
    return true;
}

// Undistorted-then-distorted normalized point of `m` at sample `p`.
void distorted(const Model& m, const Sample& p, double* dx, double* dy) {
    double x = p.u, y = p.v;
    for (int j = 0; j < kNC; j++) {
        x += m.c[j] * p.gx[j];
        y += m.c[j] * p.gy[j];
    }
    *dx = x;
    *dy = y;
}

double cost(const std::vector<Sample>& s, const Model& m) {
    double f = 0.0;
    for (const Sample& p : s) {
        double dx, dy;
        distorted(m, p, &dx, &dy);
        double ru = m.fx * dx + m.cx - p.us;
        double rv = m.fy * dy + m.cy - p.vs;
        f += p.w * (ru * ru + rv * rv);
    }
    return f;
}

// Levenberg-Marquardt on the weighted least squares. Unknown order:
// c[0..7], then fx, cx, fy, cy; `act` lists the ones that may move, so a
// masked-off coefficient keeps whatever it was initialized to (zero).
void optimize(const std::vector<Sample>& s, Model& m,
              const std::vector<int>& act) {
    const int nu = (int)act.size();
    if (nu == 0) return;
    double lambda = 1e-6;
    double f0 = cost(s, m);
    for (int it = 0; it < 200 && lambda < 1e14; it++) {
        double jtj[kNU * kNU] = {}, jtr[kNU] = {};
        for (const Sample& p : s) {
            double dx, dy;
            distorted(m, p, &dx, &dy);
            double ru = m.fx * dx + m.cx - p.us;
            double rv = m.fy * dy + m.cy - p.vs;
            double ju[kNU] = {}, jv[kNU] = {};
            for (int j = 0; j < kNC; j++) {
                ju[j] = m.fx * p.gx[j];
                jv[j] = m.fy * p.gy[j];
            }
            ju[8] = dx; ju[9] = 1.0;
            jv[10] = dy; jv[11] = 1.0;
            for (int i = 0; i < nu; i++) {
                const int ai = act[i];
                jtr[i] += p.w * (ju[ai] * ru + jv[ai] * rv);
                for (int j = 0; j <= i; j++) {
                    const int aj = act[j];
                    jtj[i * nu + j] += p.w * (ju[ai] * ju[aj] + jv[ai] * jv[aj]);
                }
            }
        }
        bool stepped = false;
        for (int attempt = 0; attempt < 24 && !stepped && lambda < 1e14; attempt++) {
            double A[kNU * kNU], b[kNU];
            std::copy(jtj, jtj + nu * nu, A);
            for (int i = 0; i < nu; i++) {
                // Marquardt scaling: fx ~ 1e3 and k1 ~ 1e-1 in one system.
                A[i * nu + i] += lambda * std::max(jtj[i * nu + i], 1e-300);
                b[i] = -jtr[i];
            }
            if (!solve_spd(nu, A, b)) { lambda *= 10.0; continue; }
            Model t = m;
            for (int i = 0; i < nu; i++) {
                switch (act[i]) {
                    case 8:  t.fx += b[i]; break;
                    case 9:  t.cx += b[i]; break;
                    case 10: t.fy += b[i]; break;
                    case 11: t.cy += b[i]; break;
                    default: t.c[act[i]] += b[i]; break;
                }
            }
            double f1 = cost(s, t);
            if (f1 < f0) {
                bool done = f0 - f1 <= 1e-14 * f0;
                m = t;
                f0 = f1;
                lambda = std::max(lambda * 0.1, 1e-12);
                stepped = true;
                if (done) return;
            } else {
                lambda *= 10.0;
            }
        }
        if (!stepped) return;
    }
}

void solve_intrinsics(const std::vector<Sample>& s, Model& m) {
    double sw = 0.0;
    double sxx = 0.0, sx = 0.0, sxu = 0.0, su = 0.0;
    double syy = 0.0, sy = 0.0, syv = 0.0, sv = 0.0;
    for (const Sample& p : s) {
        double dx, dy;
        distorted(m, p, &dx, &dy);
        sw += p.w;
        sxx += p.w * dx * dx; sx += p.w * dx; sxu += p.w * dx * p.us; su += p.w * p.us;
        syy += p.w * dy * dy; sy += p.w * dy; syv += p.w * dy * p.vs; sv += p.w * p.vs;
    }
    double dx_det = sxx * sw - sx * sx;
    if (std::fabs(dx_det) > 1e-30 * (sxx * sw + 1.0)) {
        m.fx = (sxu * sw - sx * su) / dx_det;
        m.cx = (sxx * su - sx * sxu) / dx_det;
    }
    double dy_det = syy * sw - sy * sy;
    if (std::fabs(dy_det) > 1e-30 * (syy * sw + 1.0)) {
        m.fy = (syv * sw - sy * sv) / dy_det;
        m.cy = (syy * sv - sy * syv) / dy_det;
    }
}

void errors(const std::vector<Sample>& s, const Model& m,
            std::vector<double>* e, double* rms, double* mx) {
    if (e) e->resize(s.size());
    double sum2 = 0.0, worst = 0.0;
    for (size_t i = 0; i < s.size(); i++) {
        const Sample& p = s[i];
        double dx, dy;
        distorted(m, p, &dx, &dy);
        double ru = m.fx * dx + m.cx - p.us;
        double rv = m.fy * dy + m.cy - p.vs;
        double d = std::sqrt(ru * ru + rv * rv);
        if (e) (*e)[i] = d;
        sum2 += ru * ru + rv * rv;
        worst = std::max(worst, d);
    }
    *rms = s.empty() ? 0.0 : std::sqrt(sum2 / (double)s.size());
    *mx = worst;
}

// The source's valid domain as theta_max per azimuth: star-shaped about the
// target axis, which every lens model satisfies on its FIRST branch. Only the
// first branch counts -- EUCM and the unified models fold back and re-enter the
// image mirrored, and fitting both branches at once is meaningless.
struct Domain {
    static constexpr int kAz = 128;
    double theta[kAz] = {};

    double limit(double phi) const {
        double t = phi * (kAz / (2.0 * kPi));
        t -= std::floor(t / kAz) * kAz;
        int i = (int)t;
        double f = t - i;
        return theta[i % kAz] * (1.0 - f) + theta[(i + 1) % kAz] * f;
    }
};

void scan_domain(const SourceProject& src, double w, double h, double cap,
                 Domain* dom) {
    // Where the optical axis lands, so "outward" can be measured.
    double cu = 0.5 * w, cv = 0.5 * h;
    src(0.0, 0.0, 1.0, &cu, &cv);
    for (int a = 0; a < Domain::kAz; a++) {
        double phi = 2.0 * kPi * a / Domain::kAz;
        const double cp = std::cos(phi), sp = std::sin(phi);
        auto radius = [&](double th, double* r) {
            double st = std::sin(th), ct = std::cos(th), u, v;
            if (!src(st * cp, st * sp, ct, &u, &v)) return false;
            if (!(u >= 0.0 && u <= w && v >= 0.0 && v <= h)) return false;
            *r = std::hypot(u - cu, v - cv);
            return true;
        };
        // The domain ends where the source leaves the image AND where its
        // image radius stops growing fast enough. Without the second test a
        // model that folds back INSIDE the frame -- r_d peaks and comes back,
        // which the division and unified models do -- is fitted across the
        // fold, and fitting a two-to-one map is meaningless.
        //
        // "Fast enough" is is_valid_distortion's Jdet > 0.25 read in
        // theta-space: the derivative is zero AT the fold, so any coefficient
        // error flips its sign and the ladder discards coefficients the source
        // needs. The kernels drop those pixels either way.
        const double back = cap * 1e-3;
        double r_ax = 0.0;
        const bool have_ax = radius(back, &r_ax);  // radius(0) is 0 by
        const double min_growth = 0.25 * r_ax;     //   construction
        auto inside = [&](double th) {
            double r1, r0;
            if (!radius(th, &r1)) return false;
            if (th <= back) return true;
            if (!radius(th - back, &r0)) return false;
            return have_ax ? (r1 - r0 > min_growth) : (r1 > r0);
        };
        const int kSteps = 48;
        double lo = -1.0, hi = -1.0;
        for (int i = 0; i <= kSteps; i++) {
            double th = cap * i / kSteps;
            if (inside(th)) lo = th;
            else if (lo >= 0.0) { hi = th; break; }
        }
        if (lo >= 0.0 && hi >= 0.0) {
            for (int i = 0; i < 40; i++) {
                double mid = 0.5 * (lo + hi);
                if (inside(mid)) lo = mid; else hi = mid;
            }
        }
        dom->theta[a] = std::max(lo, 0.0);
    }
}

double domain_max(const Domain& d) {
    double m = 0.0;
    for (double t : d.theta) m = std::max(m, t);
    return m;
}

std::vector<Sample> build_samples(const SourceProject& src, CameraModelType m,
                                  double w, double h, const double box[4],
                                  const Domain& dom, int n, bool chebyshev) {
    std::vector<double> xs(n), ys(n);
    for (int i = 0; i < n; i++) {
        double t = chebyshev ? std::cos(kPi * (2 * i + 1) / (2.0 * n))
                             : (n == 1 ? 0.0 : -1.0 + 2.0 * i / (n - 1));
        xs[i] = 0.5 * (box[0] + box[1]) + 0.5 * (box[1] - box[0]) * t;
        ys[i] = 0.5 * (box[2] + box[3]) + 0.5 * (box[3] - box[2]) * t;
    }
    std::vector<Sample> out;
    out.reserve((size_t)n * n);
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            Sample p;
            p.u = xs[i];
            p.v = ys[j];
            double d[3], us, vs;
            if (!unproject(m, p.u, p.v, d)) continue;
            double th = std::atan2(std::sqrt(d[0] * d[0] + d[1] * d[1]), d[2]);
            if (th > dom.limit(std::atan2(p.v, p.u))) continue;
            if (!src(d[0], d[1], d[2], &us, &vs)) continue;
            if (!(us >= 0.0 && us <= w && vs >= 0.0 && vs <= h)) continue;
            p.us = us;
            p.vs = vs;
            distort_basis(p.u, p.v, p.gx, p.gy);
            out.push_back(p);
        }
    }
    return out;
}

// The rim of the valid domain falls between grid nodes, and a near-minimax
// residual peaks there, so the reported max is only honest with it sampled.
void append_rim(const SourceProject& src, CameraModelType m, double w, double h,
                const Domain& dom, std::vector<Sample>* out) {
    const int kN = 4 * Domain::kAz;
    for (int a = 0; a < kN; a++) {
        double phi = 2.0 * kPi * a / kN;
        double rho = rho_of_theta(m, dom.limit(phi) * 0.999);
        Sample p;
        p.u = rho * std::cos(phi);
        p.v = rho * std::sin(phi);
        double d[3], us, vs;
        if (!unproject(m, p.u, p.v, d)) continue;
        if (!src(d[0], d[1], d[2], &us, &vs)) continue;
        if (!(us >= 0.0 && us <= w && vs >= 0.0 && vs <= h)) continue;
        p.us = us;
        p.vs = vs;
        distort_basis(p.u, p.v, p.gx, p.gy);
        out->push_back(p);
    }
}

// Coefficient sets tried in turn until the fitted distortion is invertible:
// everything, then no thin prism, then no k3/k4, then radial only, then k1
// alone, then none. Each rung is a subset of the one above it.
constexpr unsigned kCoeffLadder[] = {0xFFu, 0x3Fu, 0x33u, 0x03u, 0x01u, 0x00u};

std::vector<int> active_unknowns(unsigned mask, bool refine) {
    std::vector<int> act;
    for (int j = 0; j < kNC; j++)
        if (mask & (1u << j)) act.push_back(j);
    if (refine) { act.push_back(8); act.push_back(9); act.push_back(10); act.push_back(11); }
    return act;
}

// One (model, coefficient mask) attempt over an already-scanned domain.
FitResult fit_one(const SourceProject& src, double w, double h,
                  CameraModelType model, const FitOptions& opt, unsigned mask,
                  const Domain& dom) {
    FitResult res;
    res.target.model = model;
    res.target.distortion = CameraDistortionType::ThinPrism;

    double box[4] = {0.0, 0.0, 0.0, 0.0};  // xmin xmax ymin ymax
    for (int a = 0; a < Domain::kAz; a++) {
        double phi = 2.0 * kPi * a / Domain::kAz;
        double rho = rho_of_theta(model, dom.theta[a]);
        double x = rho * std::cos(phi), y = rho * std::sin(phi);
        box[0] = std::min(box[0], x); box[1] = std::max(box[1], x);
        box[2] = std::min(box[2], y); box[3] = std::max(box[3], y);
    }
    if (!(box[1] > box[0]) || !(box[3] > box[2]))
        return res;

    const int n = std::min(512, std::max(8, opt.grid));
    std::vector<Sample> fit = build_samples(src, model, w, h, box, dom, n, true);
    if ((int)fit.size() < 4 * kNC)
        return res;

    Model m;
    solve_intrinsics(fit, m);
    if (!(m.fx > 0.0) || !(m.fy > 0.0))
        return res;

    const std::vector<int> act = active_unknowns(mask, opt.refine_intrinsics);
    Model best = m;
    double best_max = 1e300;
    std::vector<double> e;
    const int lawson = std::max(0, opt.lawson_iters);
    for (int round = 0; round <= lawson; round++) {
        optimize(fit, m, act);
        double rms, mx;
        errors(fit, m, &e, &rms, &mx);
        if (mx < best_max) { best_max = mx; best = m; }
        if (round == lawson) break;
        // Lawson: reweight by residual, floored so a node cannot be retired for
        // good, then renormalized to keep the ridge scale meaningful.
        double floor_e = 1e-3 * mx, sum = 0.0;
        for (size_t i = 0; i < fit.size(); i++) {
            fit[i].w *= e[i] + floor_e;
            sum += fit[i].w;
        }
        if (!(sum > 0.0)) break;
        double scale = (double)fit.size() / sum;
        for (Sample& p : fit) p.w *= scale;
    }

    std::vector<Sample> check = build_samples(src, model, w, h, box, dom,
                                              std::min(1024, 2 * n + 1), false);
    append_rim(src, model, w, h, dom, &check);
    if (check.empty()) check = fit;
    double rms = 0.0, mx = 0.0;
    errors(check, best, nullptr, &rms, &mx);

    bool invertible = true;
    for (const Sample& p : check)
        invertible = invertible && distortion_invertible(p.u, p.v, best.c);

    res.target.fx = (float)best.fx;
    res.target.fy = (float)best.fy;
    res.target.cx = (float)best.cx;
    res.target.cy = (float)best.cy;
    for (int j = 0; j < kNC; j++) res.target.coeffs[j] = (float)best.c[j];
    res.rms_px = rms;
    res.max_px = mx;
    res.samples = (int)check.size();
    res.invertible = invertible;
    res.ok = invertible && rms <= opt.max_rms_px;
    return res;
}

// Walk the ladder until the distortion is invertible, keeping the best rung
// reached. `wide` is the domain scanned at the widest cap; it is narrowed to
// what `model` can address.
FitResult fit_model(const SourceProject& src, double w, double h,
                    CameraModelType model, const FitOptions& opt,
                    const Domain& wide) {
    FitResult res;
    res.target.model = model;
    res.target.distortion = CameraDistortionType::ThinPrism;
    if (!camera_distortion_is_compiled(model, CameraDistortionType::ThinPrism))
        return res;

    Domain dom = wide;
    const double cap = theta_cap(model);
    for (int a = 0; a < Domain::kAz; a++) dom.theta[a] = std::min(dom.theta[a], cap);

    unsigned prev = ~0u;
    for (unsigned rung : kCoeffLadder) {
        unsigned mask = rung & opt.coeff_mask;
        if (mask == prev) continue;
        prev = mask;
        FitResult r = fit_one(src, w, h, model, opt, mask, dom);
        // A rung that failed to solve at all says nothing about the next one.
        if (r.samples > 0 && (r.invertible || res.samples == 0)) res = r;
        if (r.invertible) break;
    }
    return res;
}

}  // namespace

FitResult fit_camera(const SourceProject& src, int width, int height,
                     CameraModelType target_model, const FitOptions& opt) {
    FitResult res;
    res.target.model = target_model;
    res.target.distortion = CameraDistortionType::ThinPrism;
    if (!src || width <= 0 || height <= 0)
        return res;
    const double w = (double)width, h = (double)height;
    Domain dom;
    scan_domain(src, w, h, theta_cap(target_model), &dom);
    res = fit_model(src, w, h, target_model, opt, dom);
    res.fov_deg = 2.0 * domain_max(dom) * 180.0 / kPi;
    return res;
}

FitResult fit_camera_auto(const SourceProject& src, int width, int height,
                          const FitOptions& opt) {
    FitResult res;
    res.target.distortion = CameraDistortionType::ThinPrism;
    if (!src || width <= 0 || height <= 0)
        return res;

    const double w = (double)width, h = (double)height;
    Domain dom;
    scan_domain(src, w, h, kWidestTheta, &dom);
    const double half_fov = domain_max(dom);
    res.fov_deg = 2.0 * half_fov * 180.0 / kPi;

    // Preference order, widest usable model last. A perspective target keeps
    // the straight lines straight but its tan() blows the polynomial up past
    // ~60 degrees, and equisolid only earns its keep on a true ultra-wide.
    CameraModelType cands[3];
    int nc = 0;
    if (half_fov < 60.0 * kPi / 180.0) cands[nc++] = CameraModelType::PINHOLE;
    cands[nc++] = CameraModelType::FISHEYE;
    if (half_fov > 70.0 * kPi / 180.0) cands[nc++] = CameraModelType::EQUISOLID;

    // Sub-tenth-of-a-pixel is below what bilinear resampling can resolve, so
    // the first candidate that reaches it wins outright. Past that, a later
    // candidate has to be clearly better to displace the preferred one --
    // swapping a perspective camera for a fisheye to shave 10% off the
    // residual is not worth the change of character.
    constexpr double kGoodPx = 0.1;
    bool have = false;
    for (int i = 0; i < nc; i++) {
        FitResult r = fit_model(src, w, h, cands[i], opt, dom);
        if (r.samples == 0) continue;
        bool better = !have || (r.invertible && !res.invertible) ||
                      (r.invertible == res.invertible &&
                       r.max_px < 0.5 * res.max_px);
        if (better) {
            double fov = res.fov_deg;
            res = r;
            res.fov_deg = fov;
            have = true;
        }
        if (r.invertible && r.max_px <= kGoodPx) break;
    }
    if (have)
        return res;

    // Nothing projected. Rather than fail the load, hand back a plain fisheye
    // spanning the image; the caller warns and the dataset still opens.
    res.target.model = CameraModelType::FISHEYE;
    res.target.distortion = CameraDistortionType::None;
    res.target.fx = res.target.fy = (float)(0.5 * std::max(w, h));
    res.target.cx = (float)(0.5 * w);
    res.target.cy = (float)(0.5 * h);
    for (int j = 0; j < kNC; j++) res.target.coeffs[j] = 0.0f;
    res.invertible = true;
    res.ok = false;
    return res;
}

}  // namespace dsfit

// Self-check against two lens models no tier represents exactly, at 1600x1200:
//   g++ -std=c++17 -I src -DDSFIT_SELFTEST src/data/DistortionFit.cpp -o /tmp/dfit
#ifdef DSFIT_SELFTEST
#include <cstdio>

int main() {
    auto run = [](const char* name, const dsfit::SourceProject& s, CameraModelType m) {
        dsfit::FitResult r = dsfit::fit_camera(s, 1600, 1200, m);
        std::printf("%-22s -> %-9s rms %.3e px  max %.3e px  n %d  inv %d  ok %d\n",
                    name, camera_model_to_string(m), r.rms_px, r.max_px,
                    r.samples, (int)r.invertible, (int)r.ok);
        std::printf("    f %.3f %.3f  c %.3f %.3f  k %.5g %.5g %.5g %.5g"
                    "  p %.5g %.5g  s %.5g %.5g\n",
                    r.target.fx, r.target.fy, r.target.cx, r.target.cy,
                    r.target.coeffs[0], r.target.coeffs[1], r.target.coeffs[2],
                    r.target.coeffs[3], r.target.coeffs[4], r.target.coeffs[5],
                    r.target.coeffs[6], r.target.coeffs[7]);
    };

    // COLMAP FOV: r_d = atan(2 r tan(omega/2)) / omega on the pinhole plane.
    const double omega = 0.9;
    dsfit::SourceProject fov = [&](double x, double y, double z, double* u, double* v) {
        if (z <= 1e-9) return false;
        double a = x / z, b = y / z, r = std::sqrt(a * a + b * b);
        double f = r < 1e-9 ? 1.0
                            : std::atan(2.0 * r * std::tan(0.5 * omega)) / (r * omega);
        *u = 1000.0 * a * f + 800.0;
        *v = 1000.0 * b * f + 600.0;
        return true;
    };
    run("COLMAP FOV w=0.9", fov, CameraModelType::PINHOLE);
    run("COLMAP FOV w=0.9", fov, CameraModelType::FISHEYE);

    // EUCM: d = sqrt(beta*(x^2+y^2) + z^2), denom = alpha*d + (1-alpha)*z.
    // Deliberately unguarded past its fold-back angle; Domain must find it.
    const double alpha = 0.6, beta = 1.1;
    dsfit::SourceProject eucm = [&](double x, double y, double z, double* u, double* v) {
        double d = std::sqrt(beta * (x * x + y * y) + z * z);
        double den = alpha * d + (1.0 - alpha) * z;
        if (den <= 1e-9) return false;
        *u = 600.0 * x / den + 800.0;
        *v = 600.0 * y / den + 600.0;
        return true;
    };
    run("EUCM a=0.6 b=1.1", eucm, CameraModelType::FISHEYE);
    run("EUCM a=0.6 b=1.1", eucm, CameraModelType::EQUISOLID);
    run("EUCM a=0.6 b=1.1", eucm, CameraModelType::PINHOLE);
    return 0;
}
#endif
