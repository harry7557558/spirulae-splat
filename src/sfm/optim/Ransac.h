// Generic LO-RANSAC with MSAC scoring (src/sfm/README.md).
//
// Estimator-agnostic: the caller supplies a minimal solver (`fit`, exactly
// `min_samples` points, may return several candidates), a non-minimal solver
// (`refit`, all current inliers, for the local-optimization step), and a
// squared-residual function. Scoring is MSAC (truncated quadratic), inliers are
// residual^2 < max_error^2, and the trial count adapts to the best inlier ratio
// (COLMAP's defaults: confidence 0.999, 100..10000 trials).
//
// SPRT early termination (roadmap) is not implemented; the adaptive stopping
// rule below is the MVP.
#pragma once

#include <algorithm>
#include <cmath>
#include <functional>
#include <random>
#include <vector>

namespace sfm {

struct RansacOptions {
    double max_error = 4.0;          // inlier threshold in the residual's units (px)
    double confidence = 0.999;
    int min_num_trials = 100;
    int max_num_trials = 10000;
    double min_inlier_ratio = 0.0;   // give up early if clearly below this
    unsigned seed = 0;
    int lo_iters = 10;               // local-optimization refit rounds
};

template <class Model>
struct RansacReport {
    Model model{};
    std::vector<char> inlier_mask;   // per point
    int num_inliers = 0;
    double score = 0;                // MSAC (lower is better)
    bool success = false;
    int trials = 0;
};

template <class Model>
using FitFn = std::function<std::vector<Model>(const std::vector<int>&)>;
template <class Model>
using ResidualFn = std::function<double(const Model&, int)>;

// Score a model over all points; fills inliers, returns MSAC score (lower
// better) and inlier count.
template <class Model>
static double scoreModel(const Model& m, int n, double thr2, const ResidualFn<Model>& res,
                         std::vector<char>& inliers, int& count) {
    inliers.assign(n, 0);
    count = 0;
    double score = 0;
    for (int i = 0; i < n; i++) {
        double r2 = res(m, i);
        if (r2 < thr2) {
            inliers[i] = 1;
            count++;
            score += r2;
        } else {
            score += thr2;
        }
    }
    return score;
}

template <class Model>
RansacReport<Model> loransac(int n, int min_samples, const FitFn<Model>& fit,
                             const FitFn<Model>& refit, const ResidualFn<Model>& res,
                             const RansacOptions& opt) {
    RansacReport<Model> best;
    best.score = 1e300;
    if (n < min_samples) return best;

    const double thr2 = opt.max_error * opt.max_error;
    std::mt19937 rng(opt.seed);
    std::uniform_int_distribution<int> uni(0, n - 1);

    auto consider = [&](const Model& m) {
        std::vector<char> inl;
        int cnt = 0;
        double sc = scoreModel(m, n, thr2, res, inl, cnt);
        if (cnt > best.num_inliers || (cnt == best.num_inliers && sc < best.score)) {
            best.model = m;
            best.inlier_mask = inl;
            best.num_inliers = cnt;
            best.score = sc;
            best.success = cnt >= min_samples;
        }
    };

    int max_trials = opt.max_num_trials;
    int trial = 0;
    for (; trial < max_trials; trial++) {
        // random distinct minimal sample
        std::vector<int> sample;
        sample.reserve(min_samples);
        while ((int)sample.size() < min_samples) {
            int idx = uni(rng);
            if (std::find(sample.begin(), sample.end(), idx) == sample.end()) sample.push_back(idx);
        }
        for (const Model& m : fit(sample)) consider(m);

        // adaptive stopping from the best inlier ratio so far
        if (best.num_inliers > 0) {
            double w = (double)best.num_inliers / n;
            double denom = std::log(1.0 - std::pow(w, min_samples));
            if (denom < 0) {
                double need = std::log(1.0 - opt.confidence) / denom;
                int dyn = (int)std::ceil(need);
                int cap = std::max(opt.min_num_trials, std::min(opt.max_num_trials, dyn));
                if (trial + 1 >= cap) { trial++; break; }
            }
        }
    }

    // Local optimization: iteratively refit on the current inliers.
    if (best.success && refit) {
        for (int it = 0; it < opt.lo_iters; it++) {
            std::vector<int> inl;
            for (int i = 0; i < n; i++)
                if (best.inlier_mask[i]) inl.push_back(i);
            if ((int)inl.size() < min_samples) break;
            int before = best.num_inliers;
            for (const Model& m : refit(inl)) consider(m);
            if (best.num_inliers <= before) break;  // converged
        }
    }
    best.trials = trial;
    if (opt.min_inlier_ratio > 0 && (double)best.num_inliers / n < opt.min_inlier_ratio)
        best.success = false;
    return best;
}

}  // namespace sfm
