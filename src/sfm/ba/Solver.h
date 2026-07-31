// Levenberg-Marquardt driver: owns GPU buffers, records the per-iteration
// command stream (assembly -> linear solve -> updates), and runs the
// accept/reject loop on the host (one small readback per iteration).
//
// Two linear solvers for the reduced camera system S dU = g:
//   dense - packed in-place blocked Cholesky (cholesky.slang), with S built
//           by the pair-aggregated (or atomic fallback) Schur kernels
//   cg    - implicit-Schur block-Jacobi PCG (cg.slang); S is never formed,
//           so the packed S, Y and pair-entry buffers are not allocated and
//           peak VRAM stays linear in the observation count
// Selection is automatic by problem size and a VRAM budget (see decidePaths),
// or forced with SolverOptions::solver. Optionally the dense machinery is
// kept allocated as a fallback: if CG fails to converge, the iteration is
// re-solved densely (assembly reuse, like a rejected step).
#pragma once

#include <chrono>
#include <cmath>
#include <cstdint>
#include <string>
#include <memory>
#include <vector>

#include "sfm/ba/Problem.h"
#include "sfm/vk/EmbeddedSpirv.h"
#include "sfm/vk/VkContext.h"

enum class RealCfg { F32, F64, DF64 };

inline const char* realCfgName(RealCfg c) {
    switch (c) {
        case RealCfg::F32: return "float";
        case RealCfg::F64: return "double";
        default: return "df";
    }
}
inline size_t realSize(RealCfg c) { return c == RealCfg::F32 ? 4 : 8; }

inline void packReals(std::vector<uint8_t>& out, const double* v, size_t n, RealCfg cfg) {
    out.resize(n * realSize(cfg));
    if (cfg == RealCfg::F32) {
        float* p = (float*)out.data();
        for (size_t i = 0; i < n; i++) p[i] = (float)v[i];
    } else if (cfg == RealCfg::F64) {
        memcpy(out.data(), v, n * 8);
    } else {
        float* p = (float*)out.data();
        for (size_t i = 0; i < n; i++) {
            float hi = (float)v[i];
            p[2 * i] = hi;
            p[2 * i + 1] = (float)(v[i] - hi);
        }
    }
}

inline void unpackReals(std::vector<double>& out, const uint8_t* v, size_t n, RealCfg cfg) {
    out.resize(n);
    if (cfg == RealCfg::F32) {
        const float* p = (const float*)v;
        for (size_t i = 0; i < n; i++) out[i] = p[i];
    } else if (cfg == RealCfg::F64) {
        memcpy(out.data(), v, n * 8);
    } else {
        const float* p = (const float*)v;
        for (size_t i = 0; i < n; i++) out[i] = (double)p[2 * i] + (double)p[2 * i + 1];
    }
}

enum class SolverSel { Auto, Dense, CG };
enum class CgFallback { Auto, On, Off };

struct SolverOptions {
    RealCfg real = RealCfg::F64;
    float loss_param = 1.0f;      // Huber delta / Cauchy c (unused by trivial loss)
    int max_iters = 50;
    double init_damping = 1e-2;
    double rtol = 1e-6;
    int patience = 10;
    SolverSel solver = SolverSel::Auto;
    double vram_budget_mb = 0;    // 0 = 90% of the device-local heap
    int cg_max_iters = 100;       // CG iteration cap per LM step
    double cg_tol = 0.1;          // relative residual tolerance eta
    CgFallback cg_fallback = CgFallback::Auto;
    // The kernels are compiled per (real, loss); `loss` selects the embedded
    // blob "ba_<real>_<loss>". spv_path overrides it with a module from disk
    // (a hand-compiled shader, for iteration without relinking).
    std::string loss = "trivial";
    std::string spv_path;
    int device = -1;
    bool validate = false;
    bool verbose = true;
    bool profile = false;
};

struct SolverStats {
    double initial_cost = 0, final_cost = 0;
    int iterations = 0, accepted = 0;
    double solve_seconds = 0;
    double vram_mb = 0;
    const char* solver = "dense";
    double cg_iters_total = 0;    // CG iterations summed over LM solves
    int cg_solves = 0;
    int cg_fallbacks = 0;         // LM iterations re-solved densely
};

class BundleSolver {
    enum class LinSolve { DenseObs, DensePair, CG };
    static constexpr uint32_t kCamBlk = kMaxCamDof * (kMaxCamDof + 1) / 2;  // matches cg.slang

public:
    // With `shared` the solver runs on a caller-owned persistent context:
    // device, pipelines and descriptor machinery are created once and reused
    // across solver instances (the mapper's periodic global BAs), and only the
    // problem-sized buffers are per-instance -- created in init(), freed in
    // the destructor. Without it the solver owns a scoped context, the old
    // behavior (ba_main, selftests).
    BundleSolver(BAProblem& P, const SolverOptions& opt, VkContext* shared = nullptr)
        : P_(P), opt_(opt), owned_(shared ? nullptr : new VkContext),
          ctx_(shared ? *shared : *owned_) {}

    ~BundleSolver() {
        // Persistent context: return this problem's VRAM now, not at ctx
        // teardown. (Owned context frees everything in ~VkContext anyway.)
        if (!owned_)
            for (GpuBuffer* b : ownBufs_) ctx_.destroyBuffer(*b);
    }

    void init() {
        auto prof_t0 = std::chrono::steady_clock::now();
        auto prof_lap = [&prof_t0] {
            auto t1 = std::chrono::steady_clock::now();
            double dt = std::chrono::duration<double>(t1 - prof_t0).count();
            prof_t0 = t1;
            return dt;
        };
        VkContextOptions vopt;
        vopt.needFloat64 = opt_.real == RealCfg::F64;
        vopt.needFloatAtomics = opt_.real != RealCfg::DF64;
        vopt.needInt64Atomics = opt_.real == RealCfg::DF64;
        vopt.deviceIndex = opt_.device;
        vopt.validate = opt_.validate;
        vopt.profile = opt_.profile;
        if (!ctx_.initialized()) ctx_.init(vopt);
        double t_ctx = prof_lap();

        decidePaths();

        const size_t rs = realSize(opt_.real);
        const uint64_t packed = (uint64_t)P_.n_dim * (P_.n_dim + 1) / 2;
        const bool needS = !useCG_ || haveFallback_;
        const uint32_t npart = (P_.n_dim + 255) / 256;

        auto mkReal = [&](uint64_t n) { return ctx_.createBuffer(std::max<uint64_t>(n, 1) * rs); };
        auto mkUint = [&](uint64_t n) { return ctx_.createBuffer(std::max<uint64_t>(n, 1) * 4); };

        bObs_ = mkReal(2 * (uint64_t)P_.num_obs);
        bObsImage_ = mkUint(P_.num_obs);
        bObsPoint_ = mkUint(P_.num_obs);
        bImageGroup_ = mkUint(P_.num_images);
        bGroupInfo_ = ctx_.createBuffer(std::max<size_t>(P_.groups.size(), 1) * 16);
        bPoses_ = mkReal(P_.pose_dim);
        bIntr_ = mkReal(P_.total_intr);
        bPoints_ = mkReal(3 * (uint64_t)P_.num_points);
        bObsRanges_ = mkUint(P_.num_points + 1);
        bModelObs_ = mkUint(P_.model_obs.size());
        bJcOff_ = mkUint(P_.num_obs);
        bJp_ = mkReal(6 * (uint64_t)P_.num_obs);
        bS_ = mkReal(needS ? packed : 1);
        bG_ = mkReal(P_.n_dim);
        bApp_ = mkReal(9 * (uint64_t)P_.num_points);
        bBp_ = mkReal(3 * (uint64_t)P_.num_points);
        bCost_ = mkReal(4);
        bPosesBak_ = mkReal(P_.pose_dim);
        bIntrBak_ = mkReal(P_.total_intr);
        bPointsBak_ = mkReal(3 * (uint64_t)P_.num_points);
        bPairEntries_ = mkUint(std::max<size_t>(P_.pair_entries.size(), 2));
        bPairChunks_ = mkUint(std::max<size_t>(P_.pair_chunks.size(), 2));
        // S/g are rebuilt from the per-observation Jacobians on every path, so
        // a rejected step needs no snapshot of them -- only Bp, which the
        // point back-substitution overwrites in place.
        bBp0_ = mkReal(3 * (uint64_t)P_.num_points);
        bJc_ = mkReal(P_.jc_total);
        bRes_ = mkReal(2 * (uint64_t)P_.num_obs);
        bW_ = mkReal(9 * (uint64_t)P_.num_points);
        bYp_ = mkReal(needS && P_.use_pair_schur ? 6 * (uint64_t)P_.num_obs : 1);
        bY_ = mkReal(P_.n_dim);
        // CG path buffers (1-element dummies when unused)
        bCamRanges_ = mkUint(cgAllocated_ ? P_.num_images + 1 : 1);
        bCamObs_ = mkUint(cgAllocated_ ? P_.num_obs : 1);
        bCamChunks_ = mkUint(cgAllocated_ ? P_.cam_chunks.size() : 3);
        bCgR_ = mkReal(cgAllocated_ ? P_.n_dim : 1);
        bCgZ_ = mkReal(cgAllocated_ ? P_.n_dim : 1);
        bCgP_ = mkReal(cgAllocated_ ? P_.n_dim : 1);
        bCgSp_ = mkReal(cgAllocated_ ? P_.n_dim : 1);
        bCgV_ = mkReal(cgAllocated_ ? 3 * (uint64_t)P_.num_points : 1);
        bCgB_ = mkReal(cgAllocated_ ? (uint64_t)kCamBlk * P_.num_images : 1);
        bCgM_ = mkReal(cgAllocated_ ? (uint64_t)kCamBlk * P_.num_prec_blocks : 1);
        bCgScal_ = mkReal(8);
        bCgPart_ = mkReal(cgAllocated_ ? 2 * (uint64_t)npart : 1);
        bPrecBlocks_ = mkUint(cgAllocated_ ? P_.prec_blocks.size() : 4);

        // binding order must match sfm/shaders/ba/ba.slang + cg.slang; atomic views
        // alias the same VkBuffer at the odd bindings
        std::vector<VkBuffer> binds = {
            bObs_.buf, bObsImage_.buf, bObsPoint_.buf, bImageGroup_.buf, bGroupInfo_.buf,
            bPoses_.buf, bIntr_.buf, bPoints_.buf, bObsRanges_.buf, bModelObs_.buf, bJcOff_.buf,
            bJp_.buf,
            bS_.buf, bS_.buf, bG_.buf, bG_.buf, bApp_.buf, bApp_.buf, bBp_.buf, bBp_.buf,
            bCost_.buf, bCost_.buf,
            bPairEntries_.buf, bPairChunks_.buf, bW_.buf, bYp_.buf, bY_.buf,
            bJc_.buf, bRes_.buf,
            bCamRanges_.buf, bCamObs_.buf, bCgR_.buf, bCgZ_.buf, bCgP_.buf, bCgSp_.buf,
            bCgV_.buf, bCgB_.buf, bCgM_.buf, bCgScal_.buf, bCgPart_.buf,
            bCgSp_.buf, bCgB_.buf, bCgM_.buf, bCamChunks_.buf, bPrecBlocks_.buf,
        };
        ownBufs_ = {&bObs_, &bObsImage_, &bObsPoint_, &bImageGroup_, &bGroupInfo_,
                    &bPoses_, &bIntr_, &bPoints_, &bObsRanges_, &bModelObs_, &bJcOff_,
                    &bJp_, &bS_, &bG_, &bApp_, &bBp_, &bCost_,
                    &bPosesBak_, &bIntrBak_, &bPointsBak_,
                    &bBp0_, &bJc_, &bRes_,
                    &bPairEntries_, &bPairChunks_, &bW_, &bYp_, &bY_,
                    &bCamRanges_, &bCamObs_, &bCamChunks_, &bCgR_, &bCgZ_, &bCgP_, &bCgSp_,
                    &bCgV_, &bCgB_, &bCgM_, &bCgScal_, &bCgPart_, &bPrecBlocks_};
        double t_buf = prof_lap();
        ctx_.createDescriptors(binds);
        // Fresh-from-the-driver allocations happen to arrive zeroed; memory
        // recycled within a persistent context does not -- it still holds the
        // previous solve's data, and any kernel that accumulates into a
        // buffer it assumes zero would inherit garbage. Make the guarantee
        // explicit instead of allocation-dependent. (After createDescriptors:
        // begin() binds the descriptor set, which must not still reference
        // the previous solve's freed buffers.)
        {
            VkCommandBuffer cb = ctx_.begin();
            for (GpuBuffer* b : ownBufs_) ctx_.fillZero(cb, *b);
            ctx_.submit(cb);
        }

        // Pick the Schur-kernel dof tier from the problem's largest camera dof:
        // compact (<=12) is byte-identical to the pre-OpenCV kernels, mid (<=14)
        // fits OpenCV/fisheye, wide (<=18) fits FullOpenCV/ThinPrismFisheye. A
        // problem pays only for the widest camera it uses (must match
        // kDofCompact/kDofMid/kDofWide in ba.slang).
        uint32_t maxDof = 0;
        for (const BAProblem::Group& g : P_.groups) maxDof = std::max(maxDof, 6 + g.n_intr);
        schurSuffix_ = maxDof <= 12 ? "_c" : maxDof <= 14 ? "_m" : "_w";

        // Only the entry points this problem dispatches: a cold driver cache
        // spends ~90 ms compiling each, and the module has thirty-odd, most of
        // them camera models the problem does not use and dof tiers it does not
        // need. loadPipelines skips what a shared context already built, so the
        // mapper's later solves add at most the kernels a new path wants.
        std::vector<std::string> entries = {
            "point_prep", "point_update", "cam_update", "intr_update",
            std::string("dp_accum") + schurSuffix_,
        };
        if (useCG_) {
            const char* cg[] = {"cg_cam_diag", "cg_prec_fact", "cg_prec_apply", "cg_init",
                                "cg_copy", "cg_gather", "cg_bmul", "cg_scatter", "cg_red2",
                                "cg_fin", "cg_axpy", "cg_updp"};
            entries.insert(entries.end(), std::begin(cg), std::end(cg));
        }
        if (!useCG_ || haveFallback_) {
            const char* ch[] = {"chol_diag", "chol_panel", "chol_update", "tri_fwd", "tri_bwd"};
            entries.insert(entries.end(), std::begin(ch), std::end(ch));
            if (P_.use_pair_schur) {
                entries.push_back("y_prep");
                entries.push_back(std::string("schur_pair") + schurSuffix_);
            } else {
                entries.push_back(std::string("schur_obs") + schurSuffix_);
            }
        }
        for (const BAProblem::ModelRange& mr : P_.model_ranges) {
            entries.push_back(kModels[mr.model].cost_entry);
            entries.push_back(kModels[mr.model].jac_entry);
        }
        // A shared context keeps its pipelines across solver instances; the
        // caller owning it must keep (real, loss) fixed, since the module is
        // compiled per pair (runGlobalBA's cache does -- the mapper never
        // varies them mid-run).
        if (!opt_.spv_path.empty()) {
            ctx_.loadPipelines(opt_.spv_path, entries);
        } else {
            std::string blob =
                std::string("ba_") + realCfgName(opt_.real) + "_" + opt_.loss;
            size_t words = 0;
            const uint32_t* code = sfm::findSpirv(blob.c_str(), &words);
            if (!code) {
                std::string have;
                for (size_t i = 0; const char* n = sfm::spirvBlobName(i); i++)
                    have += (i ? ", " : "") + std::string(n);
                throw std::runtime_error("shader variant '" + blob +
                                         "' is not built into this binary "
                                         "(configured with SSPLAT_SFM_REALS / "
                                         "SSPLAT_SFM_LOSSES); available: " + have);
            }
            ctx_.loadPipelines(code, words * 4, entries);
        }
        double t_pipe = prof_lap();

        // Upload static data + initial parameters, in one submit. Seventeen
        // fenced copies were most of a small solve's cost once several solvers
        // share a device (see VkContext::uploadMany).
        std::vector<uint8_t> obs, poses, intr, points;
        packReals(obs, P_.obs_xy.data(), P_.obs_xy.size(), opt_.real);
        packReals(poses, P_.poses.data(), P_.poses.size(), opt_.real);
        packReals(intr, P_.intr.data(), P_.intr.size(), opt_.real);
        packReals(points, P_.points.data(), P_.points.size(), opt_.real);
        std::vector<uint32_t> gi;
        for (auto& g : P_.groups) {
            gi.push_back(g.intr_offset);
            gi.push_back(g.intr_col);
            gi.push_back(g.n_intr);
            gi.push_back(g.model);
        }
        std::vector<VkContext::UploadItem> up = {
            {&bObs_, obs.data(), obs.size()},
            {&bObsImage_, P_.obs_image.data(), P_.obs_image.size() * 4},
            {&bObsPoint_, P_.obs_point.data(), P_.obs_point.size() * 4},
            {&bImageGroup_, P_.image_group.data(), P_.image_group.size() * 4},
            {&bGroupInfo_, gi.data(), gi.size() * 4},
            {&bObsRanges_, P_.obs_ranges.data(), P_.obs_ranges.size() * 4},
            {&bModelObs_, P_.model_obs.data(), P_.model_obs.size() * 4},
            {&bJcOff_, P_.jc_off.data(), P_.jc_off.size() * 4},
            {&bPoses_, poses.data(), poses.size()},
            {&bIntr_, intr.data(), intr.size()},
            {&bPoints_, points.data(), points.size()},
        };
        if (P_.use_pair_schur) {
            up.push_back({&bPairEntries_, P_.pair_entries.data(), P_.pair_entries.size() * 4});
            up.push_back({&bPairChunks_, P_.pair_chunks.data(), P_.pair_chunks.size() * 4});
        }
        if (cgAllocated_) {
            up.push_back({&bCamRanges_, P_.cam_obs_ranges.data(), P_.cam_obs_ranges.size() * 4});
            up.push_back({&bCamObs_, P_.cam_obs.data(), P_.cam_obs.size() * 4});
            up.push_back({&bCamChunks_, P_.cam_chunks.data(), P_.cam_chunks.size() * 4});
            up.push_back({&bPrecBlocks_, P_.prec_blocks.data(), P_.prec_blocks.size() * 4});
        }
        ctx_.uploadMany(up.data(), up.size());

        double t_upload = prof_lap();
        if (std::getenv("SSPLAT_SFM_MAP_PROF"))
            fprintf(stderr, "[prof]   solver init: ctx %.3f buf %.3f pipe %.3f upload %.3f s\n",
                    t_ctx, t_buf, t_pipe, t_upload);
        stats_.vram_mb = ctx_.totalAllocatedMB();
        stats_.solver = useCG_ ? (haveFallback_ ? "cg+fallback" : "cg") : "dense";
        if (opt_.verbose)
            fprintf(stderr, "[vk] n_dim = %u, solver = %s, VRAM allocated = %.1f MB\n",
                    P_.n_dim, stats_.solver, stats_.vram_mb);
    }

    void uploadParams() {
        std::vector<uint8_t> poses, intr, points;
        packReals(poses, P_.poses.data(), P_.poses.size(), opt_.real);
        packReals(intr, P_.intr.data(), P_.intr.size(), opt_.real);
        packReals(points, P_.points.data(), P_.points.size(), opt_.real);
        const VkContext::UploadItem up[] = {
            {&bPoses_, poses.data(), poses.size()},
            {&bIntr_, intr.data(), intr.size()},
            {&bPoints_, points.data(), points.size()},
        };
        ctx_.uploadMany(up, 3);
    }

    void downloadParams() {
        std::vector<uint8_t> tmp(std::max({bPoses_.size, bIntr_.size, bPoints_.size}));
        ctx_.download(bPoses_, tmp.data(), P_.poses.size() * realSize(opt_.real));
        unpackReals(P_.poses, tmp.data(), P_.poses.size(), opt_.real);
        ctx_.download(bIntr_, tmp.data(), P_.intr.size() * realSize(opt_.real));
        unpackReals(P_.intr, tmp.data(), P_.intr.size(), opt_.real);
        ctx_.download(bPoints_, tmp.data(), P_.points.size() * realSize(opt_.real));
        unpackReals(P_.points, tmp.data(), P_.points.size(), opt_.real);
    }

    double computeCost() {
        VkCommandBuffer cb = ctx_.begin();
        recordCost(cb);
        ctx_.barrier(cb);
        ctx_.recordDownload(cb, bCost_, realSize(opt_.real), 0, kDlCost);
        ctx_.submit(cb);
        return readCost();
    }

    void solve() {
        auto t0 = std::chrono::high_resolution_clock::now();
        double damping = opt_.init_damping;
        double cost = computeCost();
        stats_.initial_cost = cost;
        int noimprov = 0;

        bool reuse = false;  // after a reject, the assembly still matches the params
        double reject_mult = 2.0;
        int consec_fallbacks = 0;
        for (int it = 0; it < opt_.max_iters; it++) {
            if (opt_.verbose)
                fprintf(stderr, "iter %3d: cost = %.9e, damping = %.3g%s\n", it, cost, damping,
                        reuse ? " (reuse)" : "");

            LinSolve path = useCG_ ? LinSolve::CG : densePath_;
            VkCommandBuffer cb = ctx_.begin();
            recordIteration(cb, (float)damping, reuse, path);
            ctx_.submit(cb);
            double newCost = readCost();
            stats_.iterations = it + 1;

            if (path == LinSolve::CG) {
                bool conv;
                double cg_iters;
                readCgStatus(conv, cg_iters);
                stats_.cg_solves++;
                stats_.cg_iters_total += cg_iters;
                if (conv) {
                    consec_fallbacks = 0;
                    // adapt the recorded iteration cap to the observed count
                    cgMaxit_ = std::min<uint32_t>(
                        std::max<uint32_t>((uint32_t)(1.5 * cg_iters) + 8, 16),
                        (uint32_t)opt_.cg_max_iters);
                } else {
                    uint32_t usedCap = cgMaxit_;
                    cgMaxit_ = (uint32_t)opt_.cg_max_iters;
                    // a truncated-CG step is still a damped descent step; keep
                    // it if it improved the cost and only pay for the dense
                    // re-solve when the step would be rejected anyway
                    bool stepOk = std::isfinite(newCost) && newCost <= cost * (1.0 + opt_.rtol);
                    if (stepOk) consec_fallbacks = 0;
                    if (haveFallback_ && !stepOk) {
                        // discard the step and redo this iteration with the
                        // dense solver, reusing the assembly (reject flow)
                        if (opt_.verbose)
                            fprintf(stderr, "iter %3d: CG hit %u-iteration cap, dense fallback\n",
                                    it, usedCap);
                        VkCommandBuffer rb = ctx_.begin();
                        ctx_.copy(rb, bPosesBak_, bPoses_, bPoses_.size);
                        ctx_.copy(rb, bIntrBak_, bIntr_, bIntr_.size);
                        ctx_.copy(rb, bPointsBak_, bPoints_, bPoints_.size);
                        ctx_.submit(rb);
                        cb = ctx_.begin();
                        recordIteration(cb, (float)damping, true, densePath_);
                        ctx_.submit(cb);
                        newCost = readCost();
                        stats_.cg_fallbacks++;
                        if (++consec_fallbacks >= 3) {
                            useCG_ = false;  // CG is not paying off; stay dense
                            stats_.solver = "cg->dense";
                            if (opt_.verbose)
                                fprintf(stderr, "[vk] repeated CG stalls, switching to dense\n");
                        }
                    }
                }
            }

            if (std::isfinite(newCost) && newCost <= cost * (1.0 + opt_.rtol)) {
                if (newCost / cost >= 1.0 - opt_.rtol) {
                    // tie-zone accept: count toward patience and nudge lambda
                    // UP -- shrinking it on micro-improvements lets it
                    // collapse and the solver stall in an ill-conditioned
                    // plateau (observed with the df config on 871)
                    if (++noimprov >= opt_.patience) { cost = newCost; break; }
                } else {
                    noimprov = 0;
                    damping *= 1.0 / 3.0;
                }
                cost = newCost;
                stats_.accepted++;
                reuse = false;
                reject_mult = 2.0;
            } else {
                // reject: restore parameters; the assembly snapshot stays valid
                VkCommandBuffer rb = ctx_.begin();
                ctx_.copy(rb, bPosesBak_, bPoses_, bPoses_.size);
                ctx_.copy(rb, bIntrBak_, bIntr_, bIntr_.size);
                ctx_.copy(rb, bPointsBak_, bPoints_, bPoints_.size);
                ctx_.submit(rb);
                if (!std::isfinite(newCost)) {
                    if (++noimprov >= opt_.patience) break;
                } else
                    noimprov = 0;
                // rejects are cheap (assembly is reused), so search lambda
                // finely at first, but escalate on consecutive rejects
                damping *= reject_mult;
                reject_mult = std::min(reject_mult * 2.0, 32.0);
                reuse = true;
            }
        }
        stats_.final_cost = cost;
        auto t1 = std::chrono::high_resolution_clock::now();
        stats_.solve_seconds = std::chrono::duration<double>(t1 - t0).count();
        ctx_.printProfile();
    }

    const SolverStats& stats() const { return stats_; }
    VkContext& ctx() { return ctx_; }
    GpuBuffer& bufS() { return bS_; }
    GpuBuffer& bufG() { return bG_; }

    // debug: run one full assembly (no factor/solve) so S and g can be dumped
    void debugAssemble(float damping) {
        VkCommandBuffer cb = ctx_.begin();
        recordAssembly(cb, damping, false, densePath_);
        ctx_.submit(cb);
    }

    // debug: solve the same assembly with both paths and report the step
    // difference (requires the cg+fallback configuration)
    double debugCompareStep(float damping) {
        if (!useCG_ || !haveFallback_)
            throw std::runtime_error("step comparison needs --solver cg + fallback on");
        VkCommandBuffer cb = ctx_.begin();
        recordAssembly(cb, damping, false, LinSolve::CG);
        recordPCG(cb, (uint32_t)opt_.cg_max_iters);
        ctx_.barrier(cb);
        ctx_.recordDownload(cb, bCgScal_, 8 * realSize(opt_.real), 0, kDlCgScal);
        ctx_.submit(cb);
        std::vector<double> xcg = downloadG();
        bool conv;
        double cg_iters;
        readCgStatus(conv, cg_iters);
        cb = ctx_.begin();
        recordAssembly(cb, damping, true, densePath_);  // reuse the same assembly
        recordCholesky(cb);
        ctx_.submit(cb);
        std::vector<double> xd = downloadG();
        double dmax = 0, xmax = 0;
        for (uint32_t i = 0; i < P_.n_dim; i++) {
            dmax = std::max(dmax, std::fabs(xcg[i] - xd[i]));
            xmax = std::max(xmax, std::fabs(xd[i]));
        }
        double rel = dmax / std::max(xmax, 1e-300);
        printf("cmp-step lambda=%g: cg %s in %.0f iters, |dx_cg - dx_dense|_inf/|dx|_inf = %.3e\n",
               damping, conv ? "converged" : "hit cap", cg_iters, rel);
        return rel;
    }

    // record the in-place dense factor + triangular solves (public for
    // selftest). chol_update also factors the next diagonal tile, so
    // chol_diag proper only runs for the first block; the triangular solves
    // are one fused dispatch per block (see cholesky.slang).
    void recordCholesky(VkCommandBuffer cb) {
        const uint32_t n = P_.n_dim, bs = 32;
        const uint32_t nb = (n + bs - 1) / bs;
        Push p;
        p.u0 = n;
        p.u1 = 0;
        ctx_.dispatch(cb, "chol_diag", 1, p);
        ctx_.barrier(cb);
        for (uint32_t k = 0; k + 1 < nb; k++) {
            p.u1 = k;
            uint32_t below = nb - 1 - k;
            ctx_.dispatch(cb, "chol_panel", below, p);
            ctx_.barrier(cb);
            p.u2 = below * (below + 1) / 2;
            ctx_.dispatch(cb, "chol_update", p.u2, p);
            ctx_.barrier(cb);
        }
        for (uint32_t k = 0; k < nb; k++) {
            p.u1 = k;
            p.u2 = (k + 1) * bs < n ? n - (k + 1) * bs : 0;
            ctx_.dispatch(cb, "tri_fwd", std::max(1u, (p.u2 + 255) / 256), p);
            ctx_.barrier(cb);
        }
        for (int k = (int)nb - 1; k >= 0; k--) {
            p.u1 = (uint32_t)k;
            p.u2 = (uint32_t)k * bs;
            ctx_.dispatch(cb, "tri_bwd", std::max(1u, (p.u2 + 255) / 256), p);
            ctx_.barrier(cb);
        }
    }

private:
    // Choose the linear solver path from the problem shape, the options and
    // the VRAM budget, and build the host-side tables the choice needs.
    void decidePaths() {
        const bool exclusive = exclusiveGroups(P_);
        const uint64_t packed = (uint64_t)P_.n_dim * (P_.n_dim + 1) / 2;
        const bool denseOk = packed <= 0x7FFFFFFFull;  // 32-bit packed indexing
        const uint64_t pairEntries = pairEntryCount(P_);
        const bool cgOk = P_.num_obs > 0;
        const double budget = opt_.vram_budget_mb > 0 ? opt_.vram_budget_mb
                                                      : 0.9 * ctx_.deviceLocalHeapMB();
        // The pair-aggregated Schur assembly is the faster of the two dense
        // assemblies but the only one that costs memory (the entry lists and
        // the per-observation Y). Treat it as what it is -- an optional
        // accelerator -- so a dense problem that no longer fits with it
        // degrades to the atomic kernel instead of jumping to CG.
        const double denseObsMB = estimateMB(true, false, false, 0);
        const double densePairMB = estimateMB(true, false, true, pairEntries);
        const bool pairOk = exclusive && P_.num_obs > 0 && pairEntries <= kMaxPairEntries &&
                            densePairMB <= budget;
        const double denseMB = pairOk ? densePairMB : denseObsMB;
        const double cgMB = estimateMB(false, true, false, 0);
        const double bothMB = estimateMB(true, true, pairOk, pairEntries);

        // n_dim above which the cubic dense factor loses to CG. Measured on an
        // RTX 4080 Super at fp64: 871 cams (n=7839, ill-conditioned BAL set,
        // ~56 CG iters/solve) still favors dense; 1936 cams (n=17424, ~9 CG
        // iters/solve) favors CG by ~50x. The dense factor costs ~(n/8192)^3
        // seconds per iteration here, so n=8192 keeps dense under ~0.3 s/iter.
        const uint32_t kDenseMaxDim = 8192;

        switch (opt_.solver) {
            case SolverSel::Dense:
                useCG_ = false;
                if (!denseOk)
                    throw std::runtime_error(
                        "reduced system too large for the dense solver (n_dim*(n_dim+1)/2 "
                        "exceeds 32-bit packed indexing); use --solver cg");
                break;
            case SolverSel::CG:
                useCG_ = cgOk;
                if (!cgOk) {
                    fprintf(stderr, "[vk] warning: no observations, falling back to dense\n");
                    if (!denseOk) throw std::runtime_error("no usable solver path");
                }
                break;
            case SolverSel::Auto:
                useCG_ = cgOk && (P_.n_dim > kDenseMaxDim || !denseOk || denseMB > budget);
                if (!useCG_ && !denseOk)
                    throw std::runtime_error("reduced system too large for the dense solver");
                break;
        }

        haveFallback_ = false;
        if (useCG_) {
            bool want = opt_.cg_fallback == CgFallback::On ||
                        (opt_.cg_fallback == CgFallback::Auto && bothMB <= 0.5 * budget);
            haveFallback_ = want && pairOk && denseOk;
            if (opt_.cg_fallback == CgFallback::On && !haveFallback_)
                fprintf(stderr, "[vk] warning: dense fallback unavailable "
                                "(pair-Schur or packed-index limits)\n");
        }

        if (opt_.verbose)
            fprintf(stderr,
                    "[vk] VRAM estimates: dense %.0f MB (%s Schur), cg %.0f MB (budget %.0f MB)\n",
                    denseMB, pairOk ? "pair" : "per-obs", cgMB, budget);
        // Say so before the driver does. There is nothing below CG to fall back
        // to -- its footprint is the problem data plus a few vectors -- so this
        // is the point at which the answer is a smaller problem or more memory.
        const double needMB = (useCG_ ? cgMB : denseMB) + (haveFallback_ ? denseMB : 0);
        if (needMB > budget)
            fprintf(stderr,
                    "[vk] warning: the %s solver needs ~%.0f MB and the budget is %.0f MB; "
                    "this may run out of device memory\n",
                    useCG_ ? "cg" : "dense", needMB, budget);

        // host tables for the chosen paths
        P_.use_pair_schur = false;
        if ((!useCG_ || haveFallback_) && pairOk) buildPairTables(P_);
        densePath_ = P_.use_pair_schur ? LinSolve::DensePair : LinSolve::DenseObs;
        cgAllocated_ = useCG_;
        if (cgAllocated_) {
            buildCamTables(P_);
            buildPrecBlocks(P_, exclusive);
        }
        cgMaxit_ = (uint32_t)opt_.cg_max_iters;
    }

    // device-buffer footprint of a path combination, in MB (mirrors init())
    double estimateMB(bool withDense, bool withCG, bool pairTables, uint64_t pairEntries) const {
        const double rs = (double)realSize(opt_.real);
        const double n = P_.n_dim, no = P_.num_obs, np = P_.num_points, ni = P_.num_images;
        const double packed = n * (n + 1) / 2;
        double b = 0;
        b += no * (2 * rs + 12) + 4 * (double)P_.model_obs.size();  // obs + index tables
        b += 4 * ni + 16 * (double)P_.groups.size();
        b += 2 * (6 * ni + P_.total_intr) * rs;                    // poses/intr + backups
        b += 3 * np * rs * 3;                                      // points, backup, Bp0
        b += 4 * (np + 1);
        b += ((double)P_.jc_total + 8 * no) * rs;                  // Jc, Jp, res
        b += (9 + 3 + 9) * np * rs;                                // App, Bp, W
        b += 2 * n * rs;                                           // g, y
        if (withDense) {
            b += packed * rs;
            if (pairTables)
                b += 8.0 * pairEntries * 1.01 + 6 * no * rs;       // pair entries + Y
        }
        if (withCG)
            b += (4 * n + 3 * np + (2.0 * ni + (double)P_.groups.size()) * kCamBlk) * rs +
                 4 * (ni + 1) + 4 * no + 12 * (no / 1024 + ni) +
                 16 * (ni + (double)P_.groups.size());  // chunk + prec-block tables
        return b / (1024.0 * 1024.0);
    }

    void recordCost(VkCommandBuffer cb) {
        ctx_.fillZero(cb, bCost_);
        ctx_.barrier(cb);
        for (auto& mr : P_.model_ranges) {
            Push p;
            p.u0 = mr.count;
            p.u1 = mr.offset;
            p.f0 = opt_.loss_param;
            ctx_.dispatch(cb, kModels[mr.model].cost_entry, (mr.count + 255) / 256, p);
        }
        ctx_.barrier(cb);
    }

    void recordAssembly(VkCommandBuffer cb, float damping, bool reuse, LinSolve path) {
        const bool dense = path != LinSolve::CG;
        if (reuse) {
            // Params were restored after a reject; the per-observation
            // Jacobians (Jc/Jp/res/App) still match them, so skip the Jacobian
            // pass. S and g are rebuilt from those by the Schur kernels -- on
            // every path, which is why none of them is snapshotted. Bp is the
            // one thing the back-substitution overwrote, so it is.
            ctx_.copy(cb, bBp0_, bBp_, bBp_.size);
            if (dense) {
                ctx_.fillZero(cb, bS_);
                ctx_.fillZero(cb, bG_);
            }
            ctx_.barrier(cb);
        } else {
            // backup params for possible reject
            ctx_.copy(cb, bPoses_, bPosesBak_, bPoses_.size);
            ctx_.copy(cb, bIntr_, bIntrBak_, bIntr_.size);
            ctx_.copy(cb, bPoints_, bPointsBak_, bPoints_.size);

            if (dense) {
                ctx_.fillZero(cb, bS_);
                ctx_.fillZero(cb, bG_);
            }
            ctx_.fillZero(cb, bApp_);
            ctx_.fillZero(cb, bBp_);
            ctx_.barrier(cb);

            for (auto& mr : P_.model_ranges) {
                Push p;
                p.u0 = mr.count;
                p.u1 = mr.offset;
                p.f0 = opt_.loss_param;
                ctx_.dispatch(cb, kModels[mr.model].jac_entry, (mr.count + 127) / 128, p);
            }
            ctx_.barrier(cb);

            ctx_.copy(cb, bBp_, bBp0_, bBp_.size);
            ctx_.barrier(cb);
        }

        {
            Push q;
            q.u0 = P_.num_points;
            q.f0 = damping;
            ctx_.dispatch(cb, "point_prep", (P_.num_points + 255) / 256, q);
            if (path == LinSolve::CG) {  // chunked cg_cam_diag accumulates atomically
                ctx_.fillZero(cb, bCgB_);
                ctx_.fillZero(cb, bCgM_);
                ctx_.fillZero(cb, bG_);
            }
        }
        ctx_.barrier(cb);

        {
            Push p;
            p.f0 = damping;
            if (path == LinSolve::DensePair) {
                p.u0 = P_.num_obs;
                ctx_.dispatch(cb, "y_prep", (P_.num_obs + 255) / 256, p);
                ctx_.barrier(cb);
                p.u0 = P_.num_pair_chunks;
                ctx_.dispatch(cb, std::string("schur_pair") + schurSuffix_, P_.num_pair_chunks, p);
            } else if (path == LinSolve::DenseObs) {
                p.u0 = P_.num_obs;
                ctx_.dispatch(cb, std::string("schur_obs") + schurSuffix_, (P_.num_obs + 127) / 128, p);
            } else {
                p.u0 = P_.num_cam_chunks;
                p.u1 = P_.prec_exclusive ? 1 : 0;
                p.u2 = P_.num_images;  // group blocks follow the per-image ones
                ctx_.dispatch(cb, "cg_cam_diag", P_.num_cam_chunks, p);
                ctx_.barrier(cb);
                p.u0 = P_.num_prec_blocks;
                ctx_.dispatch(cb, "cg_prec_fact", (P_.num_prec_blocks + 255) / 256, p);
            }
        }
        ctx_.barrier(cb);
    }

    // Record the device-side PCG loop (see cg.slang). A fixed iteration count
    // is recorded; every kernel no-ops once the convergence flag is set.
    void recordPCG(VkCommandBuffer cb, uint32_t maxit) {
        const uint32_t n = P_.n_dim;
        const uint32_t ng = (n + 255) / 256;
        const uint32_t npart = ng;
        const uint32_t nib = (P_.num_prec_blocks + 255) / 256;
        // shared intrinsics groups: cg_bmul accumulates into the intrinsics
        // slice of Sp instead of storing it, so that slice must start at zero
        const VkDeviceSize intrOff = (VkDeviceSize)P_.pose_dim * realSize(opt_.real);
        const VkDeviceSize intrSize = (VkDeviceSize)(n - P_.pose_dim) * realSize(opt_.real);
        const bool zeroIntr = !P_.prec_exclusive && intrSize > 0;
        Push pn;
        pn.u0 = n;
        ctx_.dispatch(cb, "cg_init", ng, pn);
        ctx_.barrier(cb);
        Push pc;
        pc.u0 = P_.num_prec_blocks;
        pc.u1 = 0;  // flag was just cleared
        ctx_.dispatch(cb, "cg_prec_apply", nib, pc);
        ctx_.barrier(cb);
        Push pr;
        pr.u0 = n;
        pr.u1 = 1;
        pr.u2 = npart;
        pr.u3 = 0;
        ctx_.dispatch(cb, "cg_red2", ng, pr);
        ctx_.barrier(cb);
        Push pf;
        pf.u0 = npart;
        pf.u1 = 0;
        pf.u2 = npart;
        pf.f0 = (float)opt_.cg_tol;
        ctx_.dispatch(cb, "cg_fin", 1, pf);
        ctx_.barrier(cb);
        ctx_.dispatch(cb, "cg_copy", ng, pn);
        ctx_.barrier(cb);
        pc.u1 = 1;
        pr.u3 = 1;
        for (uint32_t it = 0; it < maxit; it++) {
            Push pg;
            pg.u0 = P_.num_points;
            ctx_.dispatch(cb, "cg_gather", (P_.num_points + 255) / 256, pg);
            ctx_.barrier(cb);
            Push ps;
            ps.u0 = P_.num_images;
            ps.u1 = P_.prec_exclusive ? 1 : 0;
            if (zeroIntr) {
                ctx_.fillZero(cb, bCgSp_, intrOff, intrSize);
                ctx_.barrier(cb);
            }
            ctx_.dispatch(cb, "cg_bmul", P_.num_images, ps);
            ctx_.barrier(cb);
            ps.u0 = P_.num_cam_chunks;
            ctx_.dispatch(cb, "cg_scatter", P_.num_cam_chunks, ps);
            ctx_.barrier(cb);
            pr.u1 = 0;
            ctx_.dispatch(cb, "cg_red2", ng, pr);
            ctx_.barrier(cb);
            pf.u1 = 1;
            ctx_.dispatch(cb, "cg_fin", 1, pf);
            ctx_.barrier(cb);
            ctx_.dispatch(cb, "cg_axpy", ng, pn);
            ctx_.barrier(cb);
            ctx_.dispatch(cb, "cg_prec_apply", nib, pc);
            ctx_.barrier(cb);
            pr.u1 = 1;
            ctx_.dispatch(cb, "cg_red2", ng, pr);
            ctx_.barrier(cb);
            pf.u1 = 2;
            ctx_.dispatch(cb, "cg_fin", 1, pf);
            ctx_.barrier(cb);
            ctx_.dispatch(cb, "cg_updp", ng, pn);
            ctx_.barrier(cb);
        }
    }

    void recordIteration(VkCommandBuffer cb, float damping, bool reuse, LinSolve path) {
        recordAssembly(cb, damping, reuse, path);

        if (path == LinSolve::CG)
            recordPCG(cb, cgMaxit_);
        else
            recordCholesky(cb);

        {
            Push p;
            p.u0 = P_.num_obs;
            ctx_.dispatch(cb, std::string("dp_accum") + schurSuffix_, (P_.num_obs + 255) / 256, p);
        }
        ctx_.barrier(cb);

        {
            Push p;
            p.u0 = P_.num_points;
            ctx_.dispatch(cb, "point_update", (P_.num_points + 255) / 256, p);
            Push q;
            q.u0 = P_.pose_dim;
            ctx_.dispatch(cb, "cam_update", (P_.pose_dim + 255) / 256, q);
            if (!P_.groups.empty()) {
                Push r;
                r.u0 = (uint32_t)P_.groups.size();
                ctx_.dispatch(cb, "intr_update", ((uint32_t)P_.groups.size() + 63) / 64, r);
            }
        }
        ctx_.barrier(cb);

        recordCost(cb);
        // Fold the two readbacks the LM loop needs into this command buffer.
        // A separate download() is its own fenced submit, so taking them here
        // halves the submits per iteration -- and a submit's latency, not its
        // arithmetic, is what a forty-image solve costs. The atom phase of a
        // bottom-up run spends five thousand iterations on problems that size,
        // with several solvers sharing the device.
        ctx_.barrier(cb);
        ctx_.recordDownload(cb, bCost_, realSize(opt_.real), 0, kDlCost);
        if (path == LinSolve::CG)
            ctx_.recordDownload(cb, bCgScal_, 8 * realSize(opt_.real), 0, kDlCgScal);
    }

    // Offsets into the download staging buffer for the folded readbacks above.
    static constexpr VkDeviceSize kDlCost = 0, kDlCgScal = 64;

    double readCost() {
        std::vector<double> v;
        unpackReals(v, (const uint8_t*)ctx_.stagingDownloadPtr() + kDlCost, 1, opt_.real);
        return v[0];
    }

    void readCgStatus(bool& converged, double& iters) {
        std::vector<double> v;
        unpackReals(v, (const uint8_t*)ctx_.stagingDownloadPtr() + kDlCgScal, 8, opt_.real);
        // flag set AND tolerance met (flag alone can also mean breakdown;
        // flag unset means the recorded iteration cap was hit)
        converged = v[6] > 0.5 && v[4] <= v[5];
        iters = v[7];
    }

    std::vector<double> downloadG() {
        std::vector<uint8_t> raw(P_.n_dim * realSize(opt_.real));
        ctx_.download(bG_, raw.data(), raw.size());
        std::vector<double> v;
        unpackReals(v, raw.data(), P_.n_dim, opt_.real);
        return v;
    }

    BAProblem& P_;
    SolverOptions opt_;
    const char* schurSuffix_ = "_c";  // "_c"/"_w" dof tier for the Schur kernels
    SolverStats stats_;
    std::unique_ptr<VkContext> owned_;      // null when running on a shared context
    VkContext& ctx_;
    std::vector<GpuBuffer*> ownBufs_;   // this instance's buffers (freed in dtor if shared)

    LinSolve densePath_ = LinSolve::DenseObs;
    bool useCG_ = false;       // CG is the active path (may demote to dense)
    bool cgAllocated_ = false; // CG buffers/tables exist
    bool haveFallback_ = false;
    uint32_t cgMaxit_ = 100;

    GpuBuffer bObs_, bObsImage_, bObsPoint_, bImageGroup_, bGroupInfo_;
    GpuBuffer bPoses_, bIntr_, bPoints_, bObsRanges_, bModelObs_, bJcOff_;
    GpuBuffer bJp_, bS_, bG_, bApp_, bBp_, bCost_;
    GpuBuffer bPosesBak_, bIntrBak_, bPointsBak_;
    GpuBuffer bBp0_, bJc_, bRes_;
    GpuBuffer bPairEntries_, bPairChunks_, bW_, bYp_, bY_;
    GpuBuffer bCamRanges_, bCamObs_, bCamChunks_, bCgR_, bCgZ_, bCgP_, bCgSp_;
    GpuBuffer bCgV_, bCgB_, bCgM_, bCgScal_, bCgPart_, bPrecBlocks_;
};
