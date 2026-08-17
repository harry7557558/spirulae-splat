#include "moge/model/Recover.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace moge {
namespace {

// The residual at one shift, and the focal that minimizes it there. With
// `fixed` positive the focal is not free and this is just the sum of squares.
double residual_at(const float* p, const float* uv, int n, double shift, double fixed,
                   double* focal_out) {
    double num = 0.0, den = 0.0;
    for (int i = 0; i < n; ++i) {
        const double z = (double)p[i * 3 + 2] + shift;
        if (!(z > 1e-12)) return std::numeric_limits<double>::infinity();
        const double px = (double)p[i * 3 + 0] / z, py = (double)p[i * 3 + 1] / z;
        num += px * (double)uv[i * 2 + 0] + py * (double)uv[i * 2 + 1];
        den += px * px + py * py;
    }
    const double focal = fixed > 0.0 ? fixed : (den > 1e-30 ? num / den : 1.0);

    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        const double z = (double)p[i * 3 + 2] + shift;
        const double du = focal * (double)p[i * 3 + 0] / z - (double)uv[i * 2 + 0];
        const double dv = focal * (double)p[i * 3 + 1] / z - (double)uv[i * 2 + 1];
        sum += du * du + dv * dv;
    }
    if (focal_out) *focal_out = focal;
    return sum;
}

}  // namespace

Recovered recover_shift(const float* points, const float* uv, int n, float focal_in) {
    Recovered r;
    r.focal = focal_in > 0.0f ? focal_in : 1.0f;
    if (n < 2) return r;

    double zmin = std::numeric_limits<double>::infinity();
    std::vector<float> zs((size_t)n);
    for (int i = 0; i < n; ++i) {
        zs[(size_t)i] = points[i * 3 + 2];
        zmin = std::min(zmin, (double)points[i * 3 + 2]);
    }
    if (!(zmin > 0.0)) return r;   // the exp remap makes this impossible
    std::nth_element(zs.begin(), zs.begin() + n / 2, zs.end());
    const double zmed = std::max((double)zs[(size_t)(n / 2)], 1e-9);

    // Quadratic spacing over (-zmin, 8*zmed]: the interesting region is the one
    // just above the barrier, where a small change in shift is a large change in
    // the implied field of view.
    constexpr int kSteps = 192;
    const double lo = -zmin * (1.0 - 1e-4), span = zmin + 8.0 * zmed;
    auto at = [&](int k) { return lo + span * ((double)k / kSteps) * ((double)k / kSteps); };

    int best = 0;
    double best_r = std::numeric_limits<double>::infinity();
    for (int k = 0; k <= kSteps; ++k) {
        const double v = residual_at(points, uv, n, at(k), focal_in, nullptr);
        if (v < best_r) { best_r = v; best = k; }
    }
    if (!std::isfinite(best_r)) return r;

    // Golden section on the bracketing interval. 80 steps takes a bracket of
    // a few percent of the scene scale to fp32's last digit.
    double a = at(std::max(best - 1, 0)), b = at(std::min(best + 1, kSteps));
    const double phi = 0.6180339887498949;
    double c = b - phi * (b - a), d = a + phi * (b - a);
    double fc = residual_at(points, uv, n, c, focal_in, nullptr);
    double fd = residual_at(points, uv, n, d, focal_in, nullptr);
    for (int it = 0; it < 80 && b - a > 1e-9 * (1.0 + std::fabs(a)); ++it) {
        if (fc < fd) {
            b = d; d = c; fd = fc;
            c = b - phi * (b - a);
            fc = residual_at(points, uv, n, c, focal_in, nullptr);
        } else {
            a = c; c = d; fc = fd;
            d = a + phi * (b - a);
            fd = residual_at(points, uv, n, d, focal_in, nullptr);
        }
    }

    double focal = r.focal;
    r.shift = (float)(0.5 * (a + b));
    residual_at(points, uv, n, r.shift, focal_in, &focal);
    r.focal = (float)focal;
    r.solved = true;
    return r;
}

}  // namespace moge
