// EvalMetrics.cpp -- see EvalMetrics.h.

#include "app/EvalMetrics.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace ssplat {

float image_l1(const float* a, const float* b, int64_t n) {
    double s = 0.0;
    for (int64_t i = 0; i < n; i++) s += std::abs((double)a[i] - (double)b[i]);
    return (float)(s / (double)std::max<int64_t>(n, 1));
}

float image_psnr(const float* a, const float* b, int64_t n) {
    double s = 0.0;
    for (int64_t i = 0; i < n; i++) {
        double d = (double)a[i] - (double)b[i];
        s += d * d;
    }
    double mse = s / (double)std::max<int64_t>(n, 1);
    if (mse <= 0.0) return std::numeric_limits<float>::infinity();
    return (float)(10.0 * std::log10(1.0 / mse));   // data_range = 1
}

namespace {

// torchmetrics' SSIM constants for data_range = 1.
constexpr double kK1 = 0.01, kK2 = 0.03;
constexpr int    kWin = 11;
constexpr double kSigma = 1.5;

// Separable gaussian, normalized -- torchmetrics builds the 2D window as the
// outer product of this with itself.
std::vector<double> gaussian_window() {
    std::vector<double> g(kWin);
    const int half = kWin / 2;
    double sum = 0.0;
    for (int i = 0; i < kWin; i++) {
        double x = (double)(i - half);
        g[i] = std::exp(-(x * x) / (2.0 * kSigma * kSigma));
        sum += g[i];
    }
    for (double& v : g) v /= sum;
    return g;
}

// torch's mode="reflect": index -1 maps to 1, not to 0 (the border pixel is
// not repeated).
int reflect(int i, int n) {
    if (i < 0) return -i;
    if (i >= n) return 2 * (n - 1) - i;
    return i;
}

// Separable gaussian blur of one channel with reflect padding: [H, W] -> [H, W].
void conv_reflect(const std::vector<double>& src, int H, int W,
                  const std::vector<double>& g, std::vector<double>& dst,
                  std::vector<double>& tmp) {
    const int p = kWin / 2;
    tmp.assign((size_t)H * W, 0.0);
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++) {
            double s = 0.0;
            for (int k = 0; k < kWin; k++)
                s += src[(size_t)y * W + reflect(x + k - p, W)] * g[k];
            tmp[(size_t)y * W + x] = s;
        }
    dst.assign((size_t)H * W, 0.0);
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++) {
            double s = 0.0;
            for (int k = 0; k < kWin; k++)
                s += tmp[(size_t)reflect(y + k - p, H) * W + x] * g[k];
            dst[(size_t)y * W + x] = s;
        }
}

}  // namespace

float image_ssim(const float* a, const float* b, int H, int W, int C) {
    if (H < kWin || W < kWin)
        throw std::runtime_error("image_ssim: image smaller than the 11x11 window");
    const std::vector<double> g = gaussian_window();
    const double C1 = (kK1 * 1.0) * (kK1 * 1.0);
    const double C2 = (kK2 * 1.0) * (kK2 * 1.0);

    std::vector<double> x((size_t)H * W), y((size_t)H * W),
                        xx((size_t)H * W), yy((size_t)H * W), xy((size_t)H * W);
    std::vector<double> mx, my, mxx, myy, mxy, tmp;

    double total = 0.0;
    for (int c = 0; c < C; c++) {
        for (int64_t p = 0; p < (int64_t)H * W; p++) {
            double va = a[p * C + c], vb = b[p * C + c];
            x[(size_t)p] = va;  y[(size_t)p] = vb;
            xx[(size_t)p] = va * va;
            yy[(size_t)p] = vb * vb;
            xy[(size_t)p] = va * vb;
        }
        conv_reflect(x,  H, W, g, mx,  tmp);
        conv_reflect(y,  H, W, g, my,  tmp);
        conv_reflect(xx, H, W, g, mxx, tmp);
        conv_reflect(yy, H, W, g, myy, tmp);
        conv_reflect(xy, H, W, g, mxy, tmp);

        double acc = 0.0;
        for (size_t i = 0; i < mx.size(); i++) {
            double ux = mx[i], uy = my[i];
            double sxx = mxx[i] - ux * ux;
            double syy = myy[i] - uy * uy;
            double sxy = mxy[i] - ux * uy;
            double num = (2.0 * ux * uy + C1) * (2.0 * sxy + C2);
            double den = (ux * ux + uy * uy + C1) * (sxx + syy + C2);
            acc += num / den;
        }
        total += acc / (double)((size_t)H * W);
    }
    return (float)(total / (double)C);
}


// --- colour correction -------------------------------------------------------
//
// The design matrix per pixel is [quadratic terms | linear | bias]:
//   c=0 -> r*r, r*g, r*b ; c=1 -> g*g, g*b ; c=2 -> b*b   (upper triangle)
//   then r, g, b, then 1.
// For C=3 that is 10 columns, so each channel solve is a 10x10 normal-equation
// system -- small enough for a plain Cholesky in double.

namespace {

int design_cols(int C) { return C * (C + 1) / 2 + C + 1; }

void build_row(const float* px, int C, double* row) {
    int k = 0;
    for (int c = 0; c < C; c++)
        for (int d = c; d < C; d++) row[k++] = (double)px[c] * (double)px[d];
    for (int c = 0; c < C; c++) row[k++] = (double)px[c];
    row[k++] = 1.0;
}

// Solve (A^T A) w = A^T b for a symmetric positive-(semi)definite normal
// matrix. Falls back to a tiny ridge if the system is singular, which happens
// when a channel's mask leaves too few unsaturated pixels.
bool solve_spd(std::vector<double>& M, std::vector<double>& v, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double s = M[(size_t)i * n + j];
            for (int k = 0; k < j; k++)
                s -= M[(size_t)i * n + k] * M[(size_t)j * n + k];
            if (i == j) {
                if (s <= 0.0) return false;
                M[(size_t)i * n + j] = std::sqrt(s);
            } else {
                M[(size_t)i * n + j] = s / M[(size_t)j * n + j];
            }
        }
    }
    for (int i = 0; i < n; i++) {                     // forward
        double s = v[i];
        for (int k = 0; k < i; k++) s -= M[(size_t)i * n + k] * v[k];
        v[i] = s / M[(size_t)i * n + i];
    }
    for (int i = n - 1; i >= 0; i--) {                // back
        double s = v[i];
        for (int k = i + 1; k < n; k++) s -= M[(size_t)k * n + i] * v[k];
        v[i] = s / M[(size_t)i * n + i];
    }
    return true;
}

}  // namespace

std::vector<float> color_correct(const float* img, const float* ref,
                                 int64_t num_pixels, int C,
                                 int num_iters, float eps) {
    const int n = design_cols(C);
    const double lo = (double)eps, hi = 1.0 - (double)eps;
    auto unclipped = [&](double z) { return z >= lo && z <= hi; };

    std::vector<float> cur(img, img + (size_t)num_pixels * C);

    // Saturation of the ORIGINAL input, fixed across iterations (mask0).
    std::vector<uint8_t> mask0((size_t)num_pixels * C);
    for (int64_t i = 0; i < num_pixels * C; i++)
        mask0[(size_t)i] = unclipped((double)img[i]) ? 1 : 0;

    std::vector<double> row((size_t)n);
    std::vector<double> warp((size_t)n * C);
    std::vector<double> M((size_t)n * n), v((size_t)n);

    for (int it = 0; it < num_iters; it++) {
        for (int c = 0; c < C; c++) {
            std::fill(M.begin(), M.end(), 0.0);
            std::fill(v.begin(), v.end(), 0.0);
            for (int64_t p = 0; p < num_pixels; p++) {
                const float* px = &cur[(size_t)p * C];
                double b = (double)ref[(size_t)p * C + c];
                // Rows saturated in the input, in the current estimate, or in
                // the reference contribute nothing.
                if (!mask0[(size_t)p * C + c] || !unclipped((double)px[c]) ||
                    !unclipped(b))
                    continue;
                build_row(px, C, row.data());
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j <= i; j++)
                        M[(size_t)i * n + j] += row[i] * row[j];
                    v[i] += row[i] * b;
                }
            }
            for (int i = 0; i < n; i++)                 // mirror to full
                for (int j = i + 1; j < n; j++)
                    M[(size_t)i * n + j] = M[(size_t)j * n + i];
            if (!solve_spd(M, v, n)) {
                // Degenerate channel: leave it untouched (identity warp).
                std::fill(v.begin(), v.end(), 0.0);
                v[(size_t)(C * (C + 1) / 2 + c)] = 1.0;
            }
            for (int i = 0; i < n; i++) warp[(size_t)i * C + c] = v[i];
        }
        // Apply the warp to update the running estimate.
        std::vector<float> next((size_t)num_pixels * C);
        for (int64_t p = 0; p < num_pixels; p++) {
            build_row(&cur[(size_t)p * C], C, row.data());
            for (int c = 0; c < C; c++) {
                double s = 0.0;
                for (int i = 0; i < n; i++) s += row[i] * warp[(size_t)i * C + c];
                next[(size_t)p * C + c] = (float)std::min(std::max(s, 0.0), 1.0);
            }
        }
        cur.swap(next);
    }
    return cur;
}

}  // namespace ssplat
