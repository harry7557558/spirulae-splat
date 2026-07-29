// Generic LO-RANSAC with MSAC scoring (src/sfm/README.md).
//
// Estimator-agnostic: the caller supplies a minimal solver (`fit`, exactly
// `min_samples` points, may return several candidates), a non-minimal solver
// (`refit`, all current inliers, for the local-optimization step), and a
// squared-residual function. Scoring is MSAC (truncated quadratic), inliers are
// residual^2 < max_error^2, and the trial count adapts to the best inlier ratio
// (COLMAP's defaults: confidence 0.999, 100..10000 trials).
//
// SPRT was measured and rejected (D26): residual evaluation is a few percent of
// RANSAC's cost, so a statistical test that stops scoring early buys almost
// nothing. What does buy something, and is exact rather than statistical, is
// the bail in `scoreModel` -- a model that cannot reach the incumbent's inlier
// count stops being scored -- and running the local optimization *inside* the
// trial loop, which raises the inlier ratio the stopping rule reads.
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
// better) and inlier count. `res` is a template parameter, not a
// std::function: this runs `trials * n` times per estimate and an indirect
// call per residual was pure overhead against a few flops of work.
//
// `need` is the inlier count a model must be able to reach to be worth
// finishing: once the points left cannot lift `count` to it, the model can
// neither beat nor tie the incumbent, so scoring stops and `count` is returned
// as -1. Most trials of a RANSAC are such models, and each of them was
// previously scored against every correspondence. The saving is exact -- a
// model that would have won is never cut, because the bail needs
// count + remaining < need, not <=.
template <class Model, class Res>
static double scoreModel(const Model& m, int n, double thr2, const Res& res,
                         std::vector<char>& inliers, int& count, int need = 0) {
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
            if (count + (n - 1 - i) < need) { count = -1; return score; }
        }
    }
    return score;
}

// `fit` / `refit` / `res` are deduced, so callers pass lambdas and everything
// inlines; the FitFn / ResidualFn aliases above still work for a caller that
// wants a type-erased one. Model stays explicit at the call sites.
template <class Model, class Fit, class Refit, class Res>
RansacReport<Model> loransac(int n, int min_samples, const Fit& fit, const Refit& refit,
                             const Res& res, const RansacOptions& opt) {
    RansacReport<Model> best;
    best.score = 1e300;
    if (n < min_samples) return best;

    const double thr2 = opt.max_error * opt.max_error;
    std::mt19937 rng(opt.seed);
    std::uniform_int_distribution<int> uni(0, n - 1);

    // Scratch reused across trials. `consider` ran thousands of times per
    // estimate and allocated a fresh mask each time; now the winning mask is
    // swapped into `best` and the loser's storage comes back as scratch.
    std::vector<char> inl;
    auto consider = [&](const Model& m) {
        int cnt = 0;
        double sc = scoreModel(m, n, thr2, res, inl, cnt, best.num_inliers);
        if (cnt < 0) return;  // bailed: cannot reach the incumbent's inlier count
        if (cnt > best.num_inliers || (cnt == best.num_inliers && sc < best.score)) {
            best.model = m;
            best.inlier_mask.swap(inl);
            best.num_inliers = cnt;
            best.score = sc;
            best.success = cnt >= min_samples;
        }
    };

    // Local optimization on the current inliers, `rounds` times or until it
    // stops helping. Used both inside the trial loop (where a better model
    // raises the inlier ratio and so lowers the trial count the stopping rule
    // demands) and once more at the end.
    std::vector<int> lo_idx;
    auto localOptimize = [&](int rounds) {
        for (int it = 0; it < rounds; it++) {
            lo_idx.clear();
            for (int i = 0; i < n; i++)
                if (best.inlier_mask[i]) lo_idx.push_back(i);
            if ((int)lo_idx.size() < min_samples) break;
            const int before = best.num_inliers;
            for (const Model& m : refit(lo_idx)) consider(m);
            if (best.num_inliers <= before) break;  // converged
        }
    };

    int max_trials = opt.max_num_trials;
    int trial = 0;
    std::vector<int> sample;
    sample.reserve(min_samples);
    for (; trial < max_trials; trial++) {
        // random distinct minimal sample
        sample.clear();
        while ((int)sample.size() < min_samples) {
            int idx = uni(rng);
            if (std::find(sample.begin(), sample.end(), idx) == sample.end()) sample.push_back(idx);
        }
        const int before = best.num_inliers;
        for (const Model& m : fit(sample)) consider(m);
        // Textbook LO-RANSAC: optimize as soon as the incumbent improves, not
        // only at the end. The refit model explains more of the data, which
        // both raises the inlier ratio the stopping rule reads (fewer trials
        // for the same confidence) and tightens the early-bail threshold
        // (cheaper trials). One round is enough here -- the full ladder runs
        // after the loop.
        if (opt.lo_iters > 0 && best.num_inliers > before && best.success) localOptimize(1);

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

    // Final local optimization. (`refit` is a deduced callable rather than a
    // std::function, so there is no "empty" state to test for -- every caller
    // passes a real one, and an estimator with no non-minimal solver passes its
    // minimal one twice.)
    if (best.success) localOptimize(opt.lo_iters);
    best.trials = trial;
    if (opt.min_inlier_ratio > 0 && (double)best.num_inliers / n < opt.min_inlier_ratio)
        best.success = false;
    return best;
}

}  // namespace sfm
