// Host-only unit test for BilagridBwdSelector (Phase 1 of the bilagrid
// backward-selection plan; see BilagridBackwardSelection.md). No GPU / backend
// calls -- it simulates per-(key,arm) GPU times and checks the bandit's
// behavior. Auto-globbed into the backend test group; also builds standalone:
//   g++ -std=c++17 -O2 bilagrid_selector_test.cpp -o /tmp/sel && /tmp/sel
//
// Exit code 0 = all checks pass, 1 = a check failed.

#include "kernels/bilagrid/BilagridBwdSelector.h"

#include <cstdio>
#include <random>
#include <string>
#include <vector>

using namespace ssplat::bilagrid;

static int g_fail = 0;
#define CHECK(cond, msg)                                                       \
    do {                                                                       \
        if (!(cond)) {                                                         \
            std::fprintf(stderr, "  FAIL: %s\n", msg);                         \
            g_fail++;                                                          \
        } else {                                                              \
            std::fprintf(stderr, "  ok:   %s\n", msg);                         \
        }                                                                      \
    } while (0)

// Deterministic simulated GPU time (ms) for an arm at a key, given per-arm
// base costs + small multiplicative noise. Separate RNG from the selector's.
struct Sim {
    std::mt19937_64 rng{12345};
    double time(const std::vector<double>& base, int arm) {
        std::uniform_real_distribution<double> jit(0.97, 1.03);
        return base[arm] * jit(rng);
    }
    int argmin(const std::vector<double>& base) {
        int b = 0;
        for (int i = 1; i < (int)base.size(); i++)
            if (base[i] < base[b]) b = i;
        return b;
    }
};

// ---------------------------------------------------------------------------
// Test A: forced init covers every arm at least min_count times on a new key.
// ---------------------------------------------------------------------------
static void test_forced_init() {
    std::fprintf(stderr, "[forced-init]\n");
    BilagridBwdSelector sel;  // default 5 arms
    ContextKey k{1, 840, 1296, 8, 16, 16};
    Sim sim;
    std::vector<double> base = {2.0, 1.95, 1.9, 1.8, 1.75, 1.7, 1.0};
    std::vector<int> picks(sel.num_arms(), 0);
    // Forced-init samples every arm to min_count WARM samples; with the first
    // (cold) sample per arm discarded, that is min_count+1 dispatches per arm.
    // Run generously so every arm is covered, then check the selector has
    // enough warm data to identify the true optimum (arm 6 = the 1.0 entry).
    for (int i = 0; i < 4 * sel.num_arms(); i++) {
        int a = sel.sample(k);
        picks[a]++;
        sel.update(k, a, sim.time(base, a));
    }
    bool all_covered = true;
    for (int i = 0; i < sel.num_arms(); i++)
        if (picks[i] < 2) all_covered = false;
    CHECK(all_covered, "every arm sampled during warmup (incl. discarded first)");
    CHECK(sel.best_known(k) == 6, "warmup identifies the true optimum (arm 6)");
}

// ---------------------------------------------------------------------------
// Test B: contextual separation -- two keys with DIFFERENT fastest arms each
// converge to their own optimum (the whole point of keying on resolution).
// ---------------------------------------------------------------------------
static void test_contextual_separation() {
    std::fprintf(stderr, "[contextual-separation]\n");
    BilagridBwdSelector sel;
    Sim sim;
    ContextKey lo{1, 200, 300, 8, 16, 16};    // low-res: scatter (arm 6) wins
    ContextKey hi{1, 2000, 3000, 8, 16, 16};  // high-res: v1 tile=5 (arm 3) wins
    std::vector<double> base_lo = {2.0, 1.95, 1.9, 1.8, 1.75, 1.7, 1.0};
    std::vector<double> base_hi = {1.5, 1.3, 1.1, 1.0, 1.15, 1.2, 2.5};

    int hit_lo = 0, hit_hi = 0, tail = 0;
    const int N = 3000, tail_start = 2500;
    for (int i = 0; i < N; i++) {
        int a_lo = sel.sample(lo);
        sel.update(lo, a_lo, sim.time(base_lo, a_lo));
        int a_hi = sel.sample(hi);
        sel.update(hi, a_hi, sim.time(base_hi, a_hi));
        if (i >= tail_start) {
            tail++;
            if (a_lo == sim.argmin(base_lo)) hit_lo++;
            if (a_hi == sim.argmin(base_hi)) hit_hi++;
        }
    }
    CHECK(sel.num_keys() == 2, "two distinct context keys tracked");
    CHECK(sel.best_known(lo) == 6, "low-res best_known == scatter (arm 6)");
    CHECK(sel.best_known(hi) == 3, "high-res best_known == v1 tile=5 (arm 3)");
    std::fprintf(stderr, "  (tail exploitation: lo=%.0f%% hi=%.0f%%)\n",
                 100.0 * hit_lo / tail, 100.0 * hit_hi / tail);
    CHECK(hit_lo > 0.7 * tail, "low-res exploits its optimum in the tail");
    CHECK(hit_hi > 0.7 * tail, "high-res exploits its optimum in the tail");
}

// ---------------------------------------------------------------------------
// Test C: elasticity -- when the fastest arm changes mid-run, the selector
// re-probes (staleness) and migrates to the new optimum.
// ---------------------------------------------------------------------------
static void test_elasticity() {
    std::fprintf(stderr, "[elasticity]\n");
    // Faster forgetting + staleness so the flip is caught within the test.
    BilagridBwdSelector::Params p;
    p.adapt_tau = 200.0;
    p.max_no_update = 64;
    BilagridBwdSelector sel(default_arms(), p);
    Sim sim;
    ContextKey k{0, 500, 500, 8, 16, 16};
    std::vector<double> before = {0.5, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0};  // arm 0 fast
    std::vector<double> after  = {2.0, 2.0, 2.0, 0.5, 2.0, 2.0, 2.0};  // arm 3 fast

    for (int i = 0; i < 1500; i++) {
        int a = sel.sample(k);
        sel.update(k, a, sim.time(before, a));
    }
    CHECK(sel.best_known(k) == 0, "before flip: converged to arm 0");
    for (int i = 0; i < 4000; i++) {
        int a = sel.sample(k);
        sel.update(k, a, sim.time(after, a));
    }
    CHECK(sel.best_known(k) == 3, "after flip: migrated to arm 3");
}

// ---------------------------------------------------------------------------
// Test D: env override pins the arm regardless of measurements.
// ---------------------------------------------------------------------------

#if defined(_WIN32) || defined(_WIN64)
inline int setenv(const char *name, const char *value, int overwrite) {
    if (!overwrite) {
        if (getenv(name) != nullptr) {
            return 0;
        }
    }
    return _putenv_s(name, value);
}
inline int unsetenv(const char *name) {
    return _putenv_s(name, "");
}
#endif

static void test_override(const char* val, int expect, const char* label) {
    setenv("SSPLAT_BILAGRID_BWD", val, 1);
    BilagridBwdSelector sel;
    unsetenv("SSPLAT_BILAGRID_BWD");  // selector already parsed it in ctor
    ContextKey k{2, 100, 100, 8, 16, 16};
    Sim sim;
    std::vector<double> base = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7};  // arm 0 fastest
    bool pinned = true;
    for (int i = 0; i < 50; i++) {
        int a = sel.sample(k);
        if (a != expect) pinned = false;
        sel.update(k, a, sim.time(base, a));
    }
    CHECK(sel.pinned() && pinned, label);
}

int main() {
    std::fprintf(stderr, "== BilagridBwdSelector unit test ==\n");
    test_forced_init();
    test_contextual_separation();
    test_elasticity();
    std::fprintf(stderr, "[override]\n");
    test_override("v2", 6, "SSPLAT_BILAGRID_BWD=v2 pins scatter (arm 6)");
    test_override("v1", 0, "SSPLAT_BILAGRID_BWD=v1 pins first v1 (arm 0)");
    test_override("v1:8", 5, "SSPLAT_BILAGRID_BWD=v1:8 pins tile=8 (arm 5)");
    test_override("#3", 3, "SSPLAT_BILAGRID_BWD=#3 pins raw index 3");

    if (g_fail) {
        std::fprintf(stderr, "\n== %d CHECK(s) FAILED ==\n", g_fail);
        return 1;
    }
    std::fprintf(stderr, "\n== all checks passed ==\n");
    return 0;
}
