// Multi-model orchestration (docs/notes/sfm-design.md D44).
//
// The mapper produces a set of models (D41) and the merger fuses the ones that
// share images (D43). Neither alone gets a fragmented capture back to one
// reconstruction, because the two operations feed each other:
//
//   * merging two models gives the result more structure, so images that could
//     not be registered into either half may register into the whole;
//   * registering more images into a model creates overlap with models it
//     shared nothing with, which is exactly what a merge needs.
//
// So this drives them in rounds -- merge, grow, prune, repeat -- until a round
// changes nothing. Everything it calls is a public operation on Mapper
// (continueFrom / refine / seedFurtherModels) or on MergeSession, so the same
// loop is available to a future bottom-up hierarchical mapper, which needs the
// identical vocabulary over models built from image clusters.
//
// Rounds are cheap when there is nothing to do: `continueFrom` returns its
// input untouched when it registers nothing, and a single-model capture skips
// the whole thing.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "sfm/core/Model.h"
#include "sfm/map/Mapper.h"
#include "sfm/map/Merge.h"
#include "sfm/map/ModelOps.h"

namespace sfm {

class ModelManager {
public:
    ModelManager(Mapper& mapper, ManagerOptions opt = {}) : mapper_(mapper), opt_(opt) {}

    const ManagerStats& stats() const { return st_; }

    std::vector<Reconstruction> run(std::vector<Reconstruction> models) {
        st_ = ManagerStats();
        st_.models_before = models.size();
        st_.covered_before = coveredImages(models).size();
        // The incoming models are the mapper's own work, checked by the rules
        // that built them; auditing one before anything has touched it can only
        // spend a bundle adjustment to confirm what it already knows. They
        // become auditable the moment a merge or a growth pass changes them,
        // which changes their shape and so their signature.
        // Before anything expensive touches them. The mapper's seed retries and
        // its sub-model search can both hand over a model that a bigger one
        // already contains (D19, D41), and until this ran only after a growth
        // pass, the first round bundle-adjusted, split, merged, audited and
        // grew every one of those copies first.
        models = dropRedundant(std::move(models));
        for (const Reconstruction& m : models) memo_.audited.insert(ModelMemo::of(m));
        if (models.size() < 2 && !opt_.do_grow && !opt_.do_reseed) {
            st_.models_after = models.size();
            st_.covered_after = st_.covered_before;
            return models;
        }
        for (int round = 0; round < std::max(1, opt_.max_rounds); round++) {
            st_.rounds = (size_t)round + 1;
            const size_t before_models = models.size();
            const size_t before_cov = coveredImages(models).size();

            if (opt_.focal_consensus_tol > 0) fixOutlierCameras(models);
            if (opt_.do_joint_ba && models.size() > 1) jointRefinePass(models);
            if (opt_.do_split) models = splitPass(models);
            if (opt_.do_merge && models.size() > 1) models = mergePass(models);
            models = auditPass(models);   // the seam, before anything is built on it
            if (opt_.do_grow) models = growPass(models);
            models = dropRedundant(models);
            if (opt_.do_reseed) models = reseedPass(models, round);

            const size_t cov = coveredImages(models).size();
            if (opt_.verbose)
                fprintf(stderr, "[mgr] round %d: %zu -> %zu models, %zu -> %zu images covered\n",
                        round + 1, before_models, models.size(), before_cov, cov);
            if (models.size() == before_models && cov == before_cov) break;
        }
        // Growth registers against the model as it stood, so the last round's
        // registrations have never been checked against the finished one.
        models = auditPass(models);
        // ... and the reseed at the end of the last round, plus that audit,
        // both change what the models cover after the loop's own prune ran.
        models = dropRedundant(models);
        // Folds are cut *last*, once, and never inside the loop.
        if (opt_.do_duplicate_split) {
            models = duplicatePass(std::move(models));
            for (Reconstruction& m : models) {
                if (m.numRegistered() < 2) continue;
                if (memo_.audited.insert(ModelMemo::of(m)).second) m = mapper_.refine(m);
            }
        }
        sortModels(models);
        st_.models_after = models.size();
        st_.covered_after = coveredImages(models).size();
        return models;
    }

private:
    // ---- rounds ----------------------------------------------------------

    std::vector<Reconstruction> mergePass(std::vector<Reconstruction> models) {
        MergeOptions mo = opt_.merge;
        mo.duplicate = opt_.duplicate;
        mo.validate = seamValidator(mapper_, opt_, &st_);
        MergeSession s(std::move(models), mo);
        const size_t merges = s.runAuto();
        st_.merges += merges;
        for (const MergeAttempt& a : s.log())
            if (!a.merged) st_.merges_refused++;
        return s.take();
    }

    std::vector<Reconstruction> splitPass(std::vector<Reconstruction> models) {
        return splitInconsistentModels(mapper_, std::move(models), opt_, memo_, st_);
    }

    std::vector<Reconstruction> duplicatePass(std::vector<Reconstruction> models) {
        return splitFoldedModels(mapper_, std::move(models), opt_, memo_, st_);
    }

    // One bundle adjustment over every component, intrinsics shared per camera
    // group. Runs before the merge pass: alignment is measured in pixels of
    // reprojection error, so a component whose focal is its own invention
    // fails merge tests it should pass.
    void jointRefinePass(std::vector<Reconstruction>& models) {
        std::map<uint32_t, std::pair<double, double>> before;  // id -> (min f, max f)
        std::map<uint32_t, double> radius;                     // ... and its frame corner
        for (const Reconstruction& m : models)
            for (const auto& kv : m.cameras) {
                auto it = before.find(kv.first);
                if (it == before.end()) {
                    before[kv.first] = {kv.second.focal(), kv.second.focal()};
                    radius[kv.first] = std::hypot(kv.second.cx, kv.second.cy);
                } else {
                    it->second.first = std::min(it->second.first, kv.second.focal());
                    it->second.second = std::max(it->second.second, kv.second.focal());
                }
            }
        // What this pass is for is a component that fitted its own focal, and
        // the disagreement between the components is how much of that there is
        // to fix. When they already agree to within a pixel at the frame
        // corner, it is a bundle adjustment over the whole capture that ends
        // where it started -- and the manage loop runs it every round.
        const double spread_px = focalSpreadPx(models);
        if (spread_px <= kJointBaMovedPx) {
            if (opt_.verbose)
                fprintf(stderr, "[mgr] joint BA skipped: the components' focals already agree to "
                        "%.2f px at the frame corner\n", spread_px);
            return;
        }
        mapper_.jointRefine(models);
        st_.joint_ba_rounds++;
        // How far the pass moved a feature at the corner of the frame, in
        // pixels: that is the unit every decision the memos hold was made in.
        double moved_px = 0;
        for (const auto& kv : before) {
            double f = 0;
            for (const Reconstruction& m : models) {
                auto it = m.cameras.find(kv.first);
                if (it == m.cameras.end()) continue;
                f = it->second.focal();
                break;
            }
            const double span = std::max(kv.second.second - f, f - kv.second.first);
            if (f > 0) moved_px = std::max(moved_px, span / f * radius[kv.first]);
            if (opt_.verbose)
                fprintf(stderr, "[mgr] joint BA: camera %u focal %.0f..%.0f -> %.0f (shared)\n",
                        kv.first, kv.second.first, kv.second.second, f);
        }
        // Every model's shape is unchanged but its geometry may not be, and a
        // model that could not grow before may grow now -- so the memos have to
        // forget them. Only when the pass actually moved something, though: the
        // memos are keyed by shape, so a model a merge or a split changed is
        // already outside them, and the pass that has to be repeated is the one
        // whose *inputs* changed while its shape did not. Round one pulls the
        // components' focals together and genuinely changes them; round two
        // onward starts from the solution round one converged to, and
        // re-growing and re-auditing every model because a focal moved by less
        // than a pixel at the frame corner is most of what the manage loop
        // spends on a large fragmented capture.
        if (moved_px <= kJointBaMovedPx) return;
        memo_.barren.clear();
        memo_.audited.clear();
        for (const Reconstruction& m : models) memo_.audited.insert(ModelMemo::of(m));
    }

    std::vector<Reconstruction> auditPass(std::vector<Reconstruction> models) {
        return auditModels(mapper_, std::move(models), opt_, memo_, st_);
    }

    // Grow every model into whatever it can still register. Largest first, so
    // the strongest model gets first claim on the contested images and the
    // smaller ones grow into what is left (and, by re-registering images the
    // big model now holds, become mergeable with it).
    std::vector<Reconstruction> growPass(std::vector<Reconstruction> models) {
        sortModels(models);
        for (size_t i = 0; i < models.size(); i++) {
            // A model that grew nothing last round will grow nothing this
            // round either, unless something changed it in between. Models are
            // moved and re-indexed by every pass, so they are remembered by
            // shape rather than by index; a collision only costs one skipped
            // attempt that was overwhelmingly likely to fail anyway.
            const ModelMemo::Signature sig = ModelMemo::of(models[i]);
            if (memo_.barren.count(sig)) continue;
            Mapper::GrowStats gs;
            std::vector<const Reconstruction*> others;
            for (size_t j = 0; j < models.size(); j++)
                if (j != i) others.push_back(&models[j]);
            Reconstruction grown = mapper_.continueFrom(models[i], &gs, others);
            if (!gs.registered) {
                memo_.barren.insert(sig);
                continue;
            }
            st_.grown_images += gs.registered;
            if (opt_.verbose)
                fprintf(stderr, "[mgr] model %zu grew %u -> %u images (%u registered)\n", i,
                        gs.before, gs.after, gs.registered);
            models[i] = std::move(grown);
        }
        return models;
    }

    std::vector<Reconstruction> dropRedundant(std::vector<Reconstruction> models) {
        return dropRedundantModels(std::move(models), opt_, st_);
    }

    // Look for models among the images nothing covers. After a growth round
    // the claimed set has changed, so seeds that were ineligible before (a
    // pair with one claimed image) may be available now -- hence the
    // relaxation ladder is re-armed.
    std::vector<Reconstruction> reseedPass(std::vector<Reconstruction> models, int round) {
        const size_t before = models.size();
        mapper_.claimAll(models);
        if (!mapper_.unclaimed()) return models;
        mapper_.seedFurtherModels(models, round > 0);
        st_.reseeded_models += models.size() - before;
        return models;
    }

    Mapper& mapper_;
    ManagerOptions opt_;
    ManagerStats st_;
    ModelMemo memo_;

    // ---- camera consensus ------------------------------------------------
    //
    // Every model of one capture is looking at the same physical cameras, so a
    // camera id that one model puts at a wildly different focal from the rest
    // is that model's mistake -- the focal search landing in a bad basin, kept
    // by a BA that had too few images to contradict it.
    //
    // The fix is to give them the consensus intrinsics and re-run the mapper's
    // refinement, which re-fits the poses and structure to them. If that was
    // the only thing wrong, the model comes back consistent with everything
    // else; if it was not, refinement filters it down and the merge tests
    // still refuse it. Weighted by image count, so the consensus comes from
    // the models with the evidence.
    void fixOutlierCameras(std::vector<Reconstruction>& models) {
        if (models.size() < 3) return;  // no population to take a consensus from
        std::map<uint32_t, std::vector<std::pair<double, uint32_t>>> focals;  // id -> (f, weight)
        std::map<uint32_t, uint32_t> widest;  // id -> most images any model gives it
        for (const Reconstruction& m : models) {
            std::map<uint32_t, uint32_t> used;
            for (const auto& kv : m.images)
                if (kv.second.registered) used[kv.second.camera_id]++;
            for (const auto& kv : m.cameras) {
                auto u = used.find(kv.first);
                if (u == used.end()) continue;
                focals[kv.first].push_back({kv.second.focal(), u->second});
                widest[kv.first] = std::max(widest[kv.first], u->second);
            }
        }
        std::map<uint32_t, double> consensus;
        for (const auto& kv : focals) {
            // A camera group of one image has no population to be an outlier
            // of: its focal is fitted from that image alone in every model that
            // holds it, so another model's value is not better evidence, and
            // "refitting" costs a full bundle adjustment of the whole model to
            // swap one guess for another. On an internet collection, where
            // every image is its own camera (D20), that fired on hundreds of
            // cameras and made the manage loop the longest stage of the run.
            if (widest[kv.first] < 2 || kv.second.size() < 3) continue;
            // Weighted median: the model with the most images decides ties,
            // and a couple of runaway small models cannot move it.
            std::vector<std::pair<double, uint32_t>> v = kv.second;
            std::sort(v.begin(), v.end());
            uint64_t total = 0;
            for (const auto& p : v) total += p.second;
            uint64_t acc = 0;
            double med = v.empty() ? 0 : v.front().first;
            for (const auto& p : v) {
                acc += p.second;
                med = p.first;
                if (acc * 2 >= total) break;
            }
            consensus[kv.first] = med;
        }
        for (Reconstruction& m : models) {
            std::vector<uint32_t> bad;
            for (const auto& kv : m.cameras) {
                auto c = consensus.find(kv.first);
                if (c == consensus.end()) continue;
                const double f = kv.second.focal(), ref = c->second;
                if (ref > 0 && std::fabs(f - ref) > opt_.focal_consensus_tol * ref)
                    bad.push_back(kv.first);
            }
            if (bad.empty()) continue;
            // Take the whole camera from the model that owns the consensus
            // focal, not just its focal: the distortion terms of a camera that
            // ran away are fitted to the same mistake.
            for (uint32_t id : bad) {
                const Camera* donor = nullptr;
                for (const Reconstruction& other : models) {
                    auto it = other.cameras.find(id);
                    if (it == other.cameras.end()) continue;
                    if (std::fabs(it->second.focal() - consensus[id]) < 1e-6 * consensus[id]) {
                        donor = &it->second;
                        break;
                    }
                }
                if (!donor) continue;
                if (opt_.verbose)
                    fprintf(stderr,
                            "[mgr] camera %u of a %u-image model: focal %.0f vs the %.0f the "
                            "other models agree on; refitting\n", id, m.numRegistered(),
                            m.cameras[id].focal(), donor->focal());
                m.cameras[id] = *donor;
                st_.cameras_refit++;
            }
            const uint32_t before = m.numRegistered();
            m = mapper_.refine(m);
            if (opt_.verbose)
                fprintf(stderr, "[mgr] refit model: %u -> %u images, %zu points\n", before,
                        m.numRegistered(), m.points3D.size());
        }
    }
};

}  // namespace sfm
