#pragma once
// Online contextual bandit that chooses the bilateral-grid grid-gradient
// backward implementation at runtime, per (family, resolution, grid shape),
// with no offline calibration phase. See BilagridBackwardSelection.md for the
// full design + rationale.
//
// Shape ported from vksplat's Gaussian Thompson-sampling scheduler (per-arm
// running stats + exponential forgetting + warmup-noise annealing + forced
// init + staleness re-probe), but made CONTEXTUAL: every context key owns its
// own arm table, so the selector does not average a fast/slow crossover across
// a run's resolutions. Because the reward (GPU time) is nearly a deterministic
// function of the key plus small noise, each key converges in a handful of
// samples, and that exploration is amortized into normal training steps -- no
// separate warmup wall-clock.
//
// Pure host bookkeeping: no CUDA / Vulkan / Slang dependency, header-only.
// Reward unit is GPU milliseconds (matches backend::event_elapsed_ms).

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <random>
#include <string>
#include <tuple>
#include <vector>

namespace ssplat {
namespace bilagrid {

// ---------------------------------------------------------------------------
// Arms
// ---------------------------------------------------------------------------

enum class BwdImpl : int { V1Gather = 0, V2Scatter = 1 };

// One arm = one concrete backward launch config.
struct BwdArm {
    BwdImpl impl;
    // Footprint-split tuning for V1Gather (feeds mult_x/mult_y in the launcher;
    // small -> aggressive split / more parallelism + atomics, large -> collapses
    // to the no-split path / minimal atomics + max serial work). Ignored (-1)
    // for V2Scatter.
    int target_tile_size;
};

// Arm sets are PER BILAGRID TYPE and TUNABLE: edit these tables after profiling
// on real hardware/datasets.
//
// v1 tile sizes span the useful range 2..8. The optimum is BOTH resolution- and
// hardware-dependent and is often NON-power-of-2: benchmarks show tile=3 winning
// at low resolution (everywhere), tile=5 at mid, tile=5/8 at high (NVIDIA/AMD
// differ) -- so 4 and 6 are included to fill the 3-5 and 5-8 gaps. Tiles >=13
// collapse to mult=1 (one thread serially scans a whole cell footprint) at all
// but extreme resolutions and never won, so they are dropped.

// PPISP (9-channel, nonlinear): tile-only. v2's atomic scatter was measured
// slower than v1 at every tested resolution/GPU (catastrophic on AMD: up to
// 384 ms), so it is intentionally excluded -- exploring it only wastes warmup.
inline const std::vector<BwdArm>& ppisp_arms() {
    static const std::vector<BwdArm> arms = {
        { BwdImpl::V1Gather, 2 }, { BwdImpl::V1Gather, 3 },
        { BwdImpl::V1Gather, 4 }, { BwdImpl::V1Gather, 5 },
        { BwdImpl::V1Gather, 6 }, { BwdImpl::V1Gather, 8 },
        { BwdImpl::V2Scatter, -1 },
    };
    return arms;
}

// Affine (12-channel, LINEAR): tiles + v2. Unlike PPISP, the affine transform
// is cheap per pixel, so v2's fused per-pixel scatter can beat v1 on some
// datasets/hardware -- worth keeping as an arm. Index 6 = v2.
inline const std::vector<BwdArm>& affine_arms() {
    static const std::vector<BwdArm> arms = {
        { BwdImpl::V1Gather, 2 }, { BwdImpl::V1Gather, 3 },
        { BwdImpl::V1Gather, 4 }, { BwdImpl::V1Gather, 5 },
        { BwdImpl::V1Gather, 6 }, { BwdImpl::V1Gather, 8 },
        { BwdImpl::V2Scatter, -1 },
    };
    return arms;
}

// Log-linear (9-channel, nonlinear exp-diagonal): tile-only. Its v2 scatter is
// IMPLEMENTED and reachable, but dropped from the default arm set -- like PPISP
// (RGB nonlinear) v2 loses (8 ms NV / 173 ms AMD vs tile ~0.6-0.9 ms). To
// re-enable, add `{ BwdImpl::V2Scatter, -1 }` here.
inline const std::vector<BwdArm>& loglinear_arms() { return ppisp_arms(); }

// Depth (2-channel) / normal (3-channel): GT-side grids, needs_image_grad=false.
// Their v2 (pure scatter) is IMPLEMENTED and reachable, but is dropped from the
// default arm set: benchmarks show it loses badly on the tiny geometry grids
// (4x8x8 -> thousands of pixels/cell -> atomic contention, ~50x the v1 gather).
// To re-enable, add `{ BwdImpl::V2Scatter, -1 }` to these tables.
inline const std::vector<BwdArm>& depth_arms() { return ppisp_arms(); }
inline const std::vector<BwdArm>& normal_arms() { return ppisp_arms(); }

// Generic default (used by unit tests): tiles + v2.
inline const std::vector<BwdArm>& default_arms() { return affine_arms(); }

// ---------------------------------------------------------------------------
// Context key: what the fast/slow choice depends on
// ---------------------------------------------------------------------------

// family: 0 affine, 1 ppisp, 2 loglinear, 3 depth, 4 normal (caller's mapping).
struct ContextKey {
    int family = 0;
    int H = 0, W = 0;        // sample resolution
    int L = 0, gH = 0, gW = 0;  // grid shape
    bool operator<(const ContextKey& o) const {
        return std::tie(family, H, W, L, gH, gW) <
               std::tie(o.family, o.H, o.W, o.L, o.gH, o.gW);
    }
};

// ---------------------------------------------------------------------------
// Selector
// ---------------------------------------------------------------------------

class BilagridBwdSelector {
public:
    struct Params {
        double initial_noise = 1.0;  // ms; injected exploration std, decays
        // Re-probe a (competitive) arm at least this often. The reward is
        // near-stationary per (device, resolution), so re-probing is mostly
        // insurance against drift -- keep it INFREQUENT so a losing-but-
        // competitive arm (e.g. v2 on Intel iGPU, ~2.25x the best -> within
        // reprobe_factor) does not tax the steady state. At 1500 an Intel v2
        // re-probe (~3 s) amortizes to ~2 ms/iter.
        int max_no_update = 1500;
        double adapt_tau = 1000.0;   // forgetting time constant (updates)
        double warmup_tau = 8.0;     // per-key noise decay constant (updates)
        int min_count = 2;           // forced samples per arm before TS
        // Only re-probe (staleness) arms whose trusted mean is within this
        // factor of the current best. Clearly-inferior arms (e.g. v2's atomic
        // scatter on AMD, ~10-60x the best v1) are measured once during forced
        // init and then never re-probed -- otherwise a periodic re-probe of a
        // catastrophic arm dominates the runtime. Generous enough (8x) to keep
        // re-probing merely-slower arms so a genuine ranking change is caught.
        double reprobe_factor = 8.0;
        uint64_t seed = 42;
    };

    // Delegating ctors rather than `= Params()` default args: a default
    // argument is not a complete-class context, so it cannot see Params's
    // in-class member initializers here (GCC rejects it).
    BilagridBwdSelector() : BilagridBwdSelector(default_arms(), Params{}) {}

    explicit BilagridBwdSelector(std::vector<BwdArm> arms)
        : BilagridBwdSelector(std::move(arms), Params{}) {}

    BilagridBwdSelector(std::vector<BwdArm> arms, Params p)
        : arms_(std::move(arms)),
          params_(p),
          adapt_beta_(std::exp(-1.0 / p.adapt_tau)),
          warmup_beta_(std::exp(-1.0 / p.warmup_tau)),
          rng_(p.seed) {
        parse_override();
    }

    const std::vector<BwdArm>& arms() const { return arms_; }
    int num_arms() const { return (int)arms_.size(); }

    // Pick an arm index for `key`. Honors the SSPLAT_BILAGRID_BWD override.
    int sample(const ContextKey& key) {
        Bucket& bk = bucket(key);
        if (forced_arm_ >= 0) return forced_arm_;

        const int n = num_arms();

        // 1) forced init: every arm gets >= min_count samples first.
        {
            int n_new = 0, pick = -1;
            for (int i = 0; i < n; i++) {
                if (bk.arms[i].count < params_.min_count) {
                    // reservoir pick among under-sampled arms (uniform)
                    std::uniform_int_distribution<int> d(0, n_new);
                    if (d(rng_) == 0) pick = i;
                    n_new++;
                }
            }
            if (pick >= 0) return pick;
        }

        // 2) staleness: re-probe the most-stale arm, but only among arms still
        //    COMPETITIVE with the best (within reprobe_factor). A clearly-worse
        //    arm with a trusted mean (>= min_count samples) is left alone, so a
        //    catastrophic arm never eats runtime on periodic re-probes.
        {
            bool have_best = false;
            double best_mean = 0.0;
            for (int i = 0; i < n; i++) {
                if (bk.arms[i].count <= 0.0) continue;
                double m = bk.arms[i].mean();
                if (!have_best || m < best_mean) { best_mean = m; have_best = true; }
            }
            int stale = -1, worst = params_.max_no_update - 1;
            for (int i = 0; i < n; i++) {
                if (bk.arms[i].num_no_update <= worst) continue;
                if (have_best && bk.arms[i].count >= params_.min_count &&
                    bk.arms[i].mean() > best_mean * params_.reprobe_factor)
                    continue;  // trusted-inferior -> don't re-probe
                worst = bk.arms[i].num_no_update;
                stale = i;
            }
            if (stale >= 0) return stale;
        }

        // 3) Thompson draw: one Gaussian per arm, pick the min predicted cost.
        int best = 0;
        double best_val = 0.0;
        for (int i = 0; i < n; i++) {
            const Stat& s = bk.arms[i];
            std::normal_distribution<double> d(s.mean(),
                                               s.stdev_mean(bk.warmup_noise));
            double v = d(rng_);
            if (i == 0 || v < best_val) {
                best_val = v;
                best = i;
            }
        }
        return best;
    }

    // Feed back the measured GPU time (ms) for the arm returned by sample().
    void update(const ContextKey& key, int arm, double gpu_ms) {
        if (arm < 0 || arm >= num_arms()) return;
        Bucket& bk = bucket(key);
        const int n = num_arms();

        // Staleness bookkeeping runs regardless of whether we keep the sample:
        // the chosen arm just ran; the others got staler.
        for (int i = 0; i < n; i++)
            if (i != arm) bk.arms[i].num_no_update++;
        bk.arms[arm].num_no_update = 0;

        // DISCARD the first measurement of each arm. The first dispatch of a
        // kernel variant pays a one-time pipeline-compile / kernel-load cost
        // (observed on Intel iGPU: first v2 ~6900 ms vs ~2000 ms steady) that
        // would otherwise inflate the arm's mean and can mislead the choice.
        // The forced-init phase keeps sampling this arm (count stays < min_count)
        // until it has min_count *warm* samples.
        if (!bk.arms[arm].warmed) {
            bk.arms[arm].warmed = true;
            return;
        }

        // Exponential forgetting of all arms + warmup-noise decay.
        for (int i = 0; i < n; i++) {
            Stat& s = bk.arms[i];
            if (s.count > params_.min_count) {
                double factor = adapt_beta_;
                double floor = params_.min_count / s.count;
                if (floor > factor) factor = floor;  // never forget below min_count
                s.count *= factor;
                s.sum *= factor;
                s.sum2 *= factor;
            }
        }
        bk.warmup_noise *= warmup_beta_;
        Stat& c = bk.arms[arm];
        c.count += 1.0;
        c.sum += gpu_ms;
        c.sum2 += gpu_ms * gpu_ms;
    }

    // --- introspection (for tests / logging) ---
    struct Stat {
        double count = 0.0, sum = 0.0, sum2 = 0.0;
        int num_no_update = 0;
        bool warmed = false;  // has this arm been dispatched at least once?
        double mean() const { return count > 0.0 ? sum / count : 0.0; }
        double var_mean() const {
            if (count < 2.0) return 0.0;
            double v = (sum2 / count - mean() * mean()) / (count - 1.0);
            return v > 0.0 ? v : 0.0;
        }
        double stdev_mean(double noise) const {
            return std::sqrt(var_mean() + noise * noise);
        }
    };

    // Current best-known arm for a key by empirical mean (no exploration draw).
    // Returns -1 if the key was never seen. Under a pin, returns the pinned arm.
    int best_known(const ContextKey& key) const {
        if (forced_arm_ >= 0) return forced_arm_;
        auto it = table_.find(key);
        if (it == table_.end()) return -1;
        const Bucket& bk = it->second;
        int best = -1;
        double best_mean = 0.0;
        for (int i = 0; i < (int)bk.arms.size(); i++) {
            if (bk.arms[i].count <= 0.0) continue;
            double m = bk.arms[i].mean();
            if (best < 0 || m < best_mean) {
                best_mean = m;
                best = i;
            }
        }
        return best;
    }

    double arm_mean(const ContextKey& key, int arm) const {
        auto it = table_.find(key);
        if (it == table_.end() || arm < 0 || arm >= (int)it->second.arms.size())
            return 0.0;
        return it->second.arms[arm].mean();
    }

    bool pinned() const { return forced_arm_ >= 0; }
    int forced_arm() const { return forced_arm_; }
    size_t num_keys() const { return table_.size(); }

    // Human-readable dump of every context key's per-arm mean cost (ms), best
    // arm marked. For SSPLAT_BILAGRID_PROFILE benchmarking.
    void dump(std::FILE* f) const {
        for (const auto& kv : table_) {
            const ContextKey& key = kv.first;
            const Bucket& bk = kv.second;
            int best = -1;
            double bestm = 0.0;
            for (int i = 0; i < (int)bk.arms.size(); i++) {
                if (bk.arms[i].count <= 0.0) continue;
                double m = bk.arms[i].mean();
                if (best < 0 || m < bestm) { bestm = m; best = i; }
            }
            std::fprintf(f, "[bilagrid-sel] fam=%d res=%dx%d grid=%dx%dx%d\n",
                         key.family, key.W, key.H, key.L, key.gH, key.gW);
            for (int i = 0; i < (int)bk.arms.size(); i++) {
                const BwdArm& a = arms_[i];
                std::fprintf(f,
                    "    arm%d %-2s tile=%-3d  mean=%.4f ms  n=%.0f%s\n", i,
                    a.impl == BwdImpl::V2Scatter ? "v2" : "v1",
                    a.target_tile_size, bk.arms[i].mean(), bk.arms[i].count,
                    i == best ? "   <== best" : "");
            }
        }
    }

private:
    struct Bucket {
        std::vector<Stat> arms;
        double warmup_noise;
    };

    Bucket& bucket(const ContextKey& key) {
        auto it = table_.find(key);
        if (it != table_.end()) return it->second;
        Bucket bk;
        bk.arms.resize(arms_.size());
        bk.warmup_noise = params_.initial_noise;
        return table_.emplace(key, std::move(bk)).first->second;
    }

    // SSPLAT_BILAGRID_BWD = auto (default) | v1 | v2 | v1:<tile> | #<arm-index>
    void parse_override() {
        forced_arm_ = -1;
        const char* e = std::getenv("SSPLAT_BILAGRID_BWD");
        if (!e || !*e) return;
        std::string s(e);
        if (s == "auto") return;
        if (!s.empty() && s[0] == '#') {  // raw arm index
            int idx = std::atoi(s.c_str() + 1);
            if (idx >= 0 && idx < num_arms()) forced_arm_ = idx;
            return;
        }
        if (s == "v2") {
            forced_arm_ = find_impl(BwdImpl::V2Scatter, -1);
            return;
        }
        if (s.rfind("v1", 0) == 0) {  // "v1" or "v1:<tile>"
            int tile = -1;  // -1 = first v1 arm
            auto colon = s.find(':');
            if (colon != std::string::npos) tile = std::atoi(s.c_str() + colon + 1);
            forced_arm_ = find_impl(BwdImpl::V1Gather, tile);
            return;
        }
        // unknown token -> leave auto
    }

    int find_impl(BwdImpl impl, int tile) const {
        int first = -1;
        for (int i = 0; i < num_arms(); i++) {
            if (arms_[i].impl != impl) continue;
            if (first < 0) first = i;
            if (tile < 0 || arms_[i].target_tile_size == tile) return i;
        }
        return first;  // requested tile absent -> first arm of that impl
    }

    std::vector<BwdArm> arms_;
    Params params_;
    double adapt_beta_;
    double warmup_beta_;
    std::map<ContextKey, Bucket> table_;
    std::mt19937_64 rng_;
    int forced_arm_ = -1;
};

}  // namespace bilagrid
}  // namespace ssplat
