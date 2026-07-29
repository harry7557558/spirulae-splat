// Lightweight wall-clock accumulators for mapper profiling (host-time work,
// src/sfm/README.md "The mapper is the stage that is left"). Zero-dependency,
// always compiled; reporting is opt-in via SSPLAT_SFM_MAP_PROF=1 so normal runs
// produce identical output.
#pragma once

#include <chrono>
#include <cstdio>
#include <cstdlib>

namespace sfm {

struct MapProf {
    // mapper host phases
    double init_seed = 0;   // initialize(): seed search incl. two-view RANSAC
    double choose = 0;      // chooseNextImages(): ranking candidates
    double reg = 0;         // registerImage(): PnP RANSAC + refine + recount
    double tri = 0;         // triangulateForImage() during growth
    double retri = 0;       // completeAndRetriangulate()
    double merge = 0;       // mergeTracks(), inside the above
    double filter = 0;      // sanitizeCameras + filterPoints + filterImages
    double snapshot = 0;    // checkedRefine() reconstruction copy
    // global BA, split (sfm/map/Bundle.h)
    double ba_build = 0;    // BAProblem assembly from the Reconstruction
    double ba_init = 0;     // BundleSolver ctor + init (device + pipelines + upload)
    double ba_solve = 0;    // solver.solve()
    double ba_write = 0;    // download + writeback
    long n_ba = 0, n_ba_iters = 0, n_choose = 0, n_reg_try = 0, n_reg_ok = 0;
    long n_merged = 0;      // observations absorbed by track merging

    static bool enabled() {
        static bool e = std::getenv("SSPLAT_SFM_MAP_PROF") != nullptr;
        return e;
    }

    void report(double total_s) const {
        if (!enabled()) return;
        double ba = ba_build + ba_init + ba_solve + ba_write;
        double accounted = init_seed + choose + reg + tri + retri + filter + snapshot + ba;
        fprintf(stderr,
                "[prof] mapper total %.2f s, accounted %.2f s (%.0f%%)\n"
                "[prof]   seed search   %8.2f s\n"
                "[prof]   choose-next   %8.2f s  (%ld calls)\n"
                "[prof]   register      %8.2f s  (%ld tries, %ld ok)\n"
                "[prof]   triangulate   %8.2f s\n"
                "[prof]   retriangulate %8.2f s  (of which merge %.2f s, %ld obs absorbed)\n"
                "[prof]   filter        %8.2f s\n"
                "[prof]   snapshot      %8.2f s\n"
                "[prof]   BA            %8.2f s  (%ld calls, %ld LM iters): "
                "build %.2f + init %.2f + solve %.2f + write %.2f\n",
                total_s, accounted, total_s > 0 ? 100.0 * accounted / total_s : 0.0,
                init_seed, choose, n_choose, reg, n_reg_try, n_reg_ok, tri, retri, merge,
                n_merged, filter, snapshot, ba, n_ba, n_ba_iters, ba_build, ba_init, ba_solve,
                ba_write);
    }
};

inline MapProf g_map_prof;

// RAII accumulator: adds elapsed seconds to `acc` at scope exit.
class ProfTimer {
public:
    explicit ProfTimer(double& acc)
        : acc_(acc), t0_(std::chrono::steady_clock::now()) {}
    ~ProfTimer() {
        acc_ += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0_).count();
    }

private:
    double& acc_;
    std::chrono::steady_clock::time_point t0_;
};

}  // namespace sfm
