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

namespace sfm {

struct ManagerOptions {
    MergeOptions merge;
    int max_rounds = 4;
    bool do_merge = true;
    bool do_grow = true;      // register unregistered images into existing models
    bool do_reseed = true;    // look for models among images still uncovered
    // Audit a merged model's poses against the correspondence graph and
    // re-register what it cannot support (Mapper::audit). Off means a merged
    // model is only bundle-adjusted, which cannot undo a gross misplacement.
    bool do_audit = true;
    // A model this fraction of whose images another (larger) model also holds
    // carries nothing of its own: growth has made it a duplicate. It is
    // dropped rather than written as a reconstruction in its own right.
    double redundant_ratio = 0.9;
    // Cameras whose focal deviates from the population consensus by more than
    // this factor are refit before merging (see fixOutlierCameras). A model
    // built on a runaway focal cannot align with anything, and on a rig every
    // model is looking at the same physical cameras.
    double focal_consensus_tol = 0.15;  // 0 disables
    // Bundle-adjust every component in one problem with intrinsics shared per
    // camera group (Mapper::jointRefine, D45). This is what stops a small
    // component from fitting its own focal to its own noise; the outlier refit
    // above still runs first, because a camera that has already run away needs
    // replacing, not averaging.
    bool do_joint_ba = true;
    // Cross-seam agreement required of a merge (Mapper::checkSeam). A verified
    // pair whose two images end up on opposite sides of the seam is evidence
    // the alignment never saw; the merged model has to reproduce most of them.
    // 0 disables the check.
    double seam_min_agreement = 0.6;   // fraction of tested cross-seam pairs
    double seam_min_pair_fraction = 0.5;  // ... of a pair's matches, to call it holding
    int seam_min_pairs = 10;           // below this there is nothing to judge on
    // Break a model whose own verified pairs do not agree with it
    // (Mapper::splitInconsistent). Only pairs with at least
    // `split_min_matches` are trusted to vote, and a resulting group smaller
    // than `split_min_group` is de-registered rather than written out.
    bool do_split = true;
    int split_min_matches = 30;
    size_t split_min_group = 10;
    // Folded models: two similar-looking parts of the capture written on top of
    // each other (findDuplicateStructure). Detected from missing structure, not
    // from disagreement, because a fold agrees with itself.
    //
    // On by default since D46, when the verdict stopped being the conflict
    // *rate* -- which does not separate a fold from a dense capture -- and
    // became what the implied split would cost. Measured with the same options
    // on real-world datasets, the rate is *higher* on the sound model; the cut
    // differs by a hundredfold, because two places written on top of each other
    // never shared structure to  begin with, while cutting a sound model always
    // runs through structure it genuinely shares.
    // See DuplicateOptions::max_cut_fraction.
    bool do_duplicate_split = true;
    DuplicateOptions duplicate;
    bool verbose = true;
};

struct ManagerStats {
    size_t rounds = 0;
    size_t merges = 0, merges_refused = 0;
    size_t grown_images = 0;       // images registered by growth passes
    size_t reseeded_models = 0;
    size_t dropped_redundant = 0;
    size_t audited_out = 0, audited_repaired = 0;
    size_t cameras_refit = 0;
    size_t joint_ba_rounds = 0;
    size_t seam_checked = 0, seam_refused = 0;
    size_t splits = 0, duplicate_splits = 0, split_dropped = 0;
    size_t models_before = 0, models_after = 0;
    size_t covered_before = 0, covered_after = 0;
};

// The cross-seam agreement test, as a merge validator (D45).
//
// Every test MergeSession can run by itself is computed from the evidence the
// alignment already used, and repeated structure satisfies all of it -- a fold
// agrees with itself. This one asks the correspondence graph instead: do the
// verified two-view geometries that cross the seam still hold in the merged
// model? It is the only test that separates "these two places look the same"
// from "these two places are the same".
//
// A free function because the bottom-up merge tree needs it as much as the
// manage loop does. Its levels perform hundreds of merges, and without this
// they run with nothing watching for repeated structure at all.
inline std::function<std::string(const Reconstruction&, const Reconstruction&, const Sim3&)>
seamValidator(Mapper& mapper, const ManagerOptions& opt, ManagerStats* st = nullptr) {
    if (opt.seam_min_agreement <= 0) return {};
    Mapper* mp = &mapper;
    const double max_err = opt.merge.max_reproj_error;
    const double pair_frac = opt.seam_min_pair_fraction;
    const double min_agree = opt.seam_min_agreement;
    const size_t min_pairs = (size_t)std::max(0, opt.seam_min_pairs);
    return [mp, max_err, pair_frac, min_agree, min_pairs, st](
               const Reconstruction& merged, const Reconstruction& src, const Sim3&) -> std::string {
        std::set<uint32_t> src_side;
        for (const auto& kv : src.images)
            if (kv.second.registered) src_side.insert(kv.first);
        Mapper::SeamCheck sc = mp->checkSeam(merged, src_side, max_err, pair_frac);
        if (sc.tested < min_pairs) return "";  // nothing to judge on
        const double agree = (double)sc.agree / (double)sc.tested;
        if (st) st->seam_checked++;
        if (agree >= min_agree) return "";
        if (st) st->seam_refused++;
        char buf[224];
        snprintf(buf, sizeof buf,
                 "only %zu of the %zu verified pairs that cross the seam still hold "
                 "(median %.0f%% of their matches explained)",
                 sc.agree, sc.tested, 100.0 * sc.median_frac);
        return std::string(buf);
    };
}

// Distinct images any model registers.
inline std::set<uint32_t> coveredImages(const std::vector<Reconstruction>& models) {
    std::set<uint32_t> ids;
    for (const Reconstruction& m : models)
        for (const auto& kv : m.images)
            if (kv.second.registered) ids.insert(kv.first);
    return ids;
}

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
        for (const Reconstruction& m : models) audited_.insert(signature(m));
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
                if (audited_.insert(signature(m)).second) m = mapper_.refine(m);
            }
        }
        sortModels(models);
        st_.models_after = models.size();
        st_.covered_after = coveredImages(models).size();
        return models;
    }

private:
    static void sortModels(std::vector<Reconstruction>& models) {
        std::stable_sort(models.begin(), models.end(),
                         [](const Reconstruction& a, const Reconstruction& b) {
                             return a.points3D.size() > b.points3D.size();
                         });
    }

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

    // Break models the correspondence graph says are not one thing.
    //
    // The merge pass refuses a bad merge before it commits, but models arrive
    // here that were never merged and are wrong anyway -- a chain of
    // registrations through repeated structure does the same damage, and D41's
    // own output has models like that. Splitting them is what makes the rest of
    // the loop work: each piece is a sound reconstruction, and the merge pass
    // gets to try again with the seam test watching.
    //
    // Only models that changed since they were last split are re-examined, so
    // a split and a merge cannot chase each other round after round.
    std::vector<Reconstruction> splitPass(std::vector<Reconstruction> models) {
        std::vector<Reconstruction> out;
        for (Reconstruction& m : models) {
            const Signature sig = signature(m);
            if (m.numRegistered() < 2 * opt_.split_min_group || !split_.insert(sig).second) {
                out.push_back(std::move(m));
                continue;
            }
            Mapper::SplitStats ss;
            std::vector<Reconstruction> parts = mapper_.splitInconsistent(
                m, opt_.merge.max_reproj_error, opt_.seam_min_pair_fraction,
                opt_.split_min_matches, opt_.split_min_group, &ss);
            if (parts.size() <= 1) {
                out.push_back(std::move(m));
                continue;
            }
            st_.splits++;
            st_.split_dropped += ss.dropped_images;
            if (opt_.verbose)
                fprintf(stderr,
                        "[mgr] split a %u-image model into %zu (%zu of %zu inner pairs hold; "
                        "largest group %zu, %zu images dropped)\n",
                        m.numRegistered(), parts.size(), ss.pairs_agree, ss.pairs_tested,
                        ss.largest, ss.dropped_images);
            for (Reconstruction& p : parts) {
                split_.insert(signature(p));
                out.push_back(std::move(p));
            }
        }
        return out;
    }

    // Break a model that has written two places on top of each other. This is
    // the failure the agreement split above cannot see -- a fold agrees with
    // itself -- and it is detected from what the model is *missing*: images it
    // places in the same spot looking the same way, with no structure in
    // common (findDuplicateStructure).
    std::vector<Reconstruction> duplicatePass(std::vector<Reconstruction> models) {
        std::vector<Reconstruction> out;
        MatchedFn matched = mapper_.matchedPredicate();
        for (Reconstruction& m : models) {
            if (m.numRegistered() < 2 * opt_.split_min_group) {
                out.push_back(std::move(m));
                continue;
            }
            DuplicateReport dr = findDuplicateStructure(m, opt_.duplicate, matched);
            if (!dr.duplicated(opt_.duplicate)) {
                out.push_back(std::move(m));
                continue;
            }
            size_t dropped = 0;
            DuplicateCut cut;
            std::vector<Reconstruction> parts =
                splitDuplicateStructure(m, dr, opt_.split_min_group, &dropped, &cut);
            // The conflicts say a fold is possible; the cut says whether the
            // split is cheap. Only a genuine fold is both (D46) -- cutting a
            // sound model always runs through structure it really does share,
            // and that is what saves drjohnson here.
            if (parts.size() <= 1 || !foldSplitAccepted(dr, cut, opt_.duplicate)) {
                if (opt_.verbose && parts.size() > 1)
                    fprintf(stderr,
                            "[mgr] a %u-image model has %zu of %zu co-located pairs with nothing "
                            "in common, but splitting it would sever %.1f%% of its co-visibility "
                            "(>%.1f%%): keeping it whole\n",
                            m.numRegistered(), dr.conflicts, dr.colocated, 100.0 * cut.fraction(),
                            100.0 * opt_.duplicate.max_cut_fraction);
                out.push_back(std::move(m));
                continue;
            }
            st_.duplicate_splits++;
            st_.split_dropped += dropped;
            if (opt_.verbose) {
                fprintf(stderr,
                        "[mgr] a %u-image model has %zu of %zu co-located image pairs with no "
                        "structure in common and no match either (%zu more share nothing but "
                        "were matched) and a cut that severs %.2f%% of its co-visibility: two "
                        "places written on top of each other. Splitting into",
                        m.numRegistered(), dr.conflicts, dr.colocated, dr.unmatched_but_seen,
                        100.0 * cut.fraction());
                for (const Reconstruction& p : parts) fprintf(stderr, " %u", p.numRegistered());
                if (dropped) fprintf(stderr, " (%zu images dropped)", dropped);
                fprintf(stderr, "\n");
            }
            for (Reconstruction& p : parts) {
                split_.insert(signature(p));
                out.push_back(std::move(p));
            }
        }
        return out;
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
        double spread_px = 0;
        for (const auto& kv : before) {
            const double f = 0.5 * (kv.second.first + kv.second.second);
            if (f > 0) spread_px = std::max(spread_px,
                                            (kv.second.second - kv.second.first) / f *
                                                radius[kv.first]);
        }
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
        barren_.clear();
        audited_.clear();
        for (const Reconstruction& m : models) audited_.insert(signature(m));
    }

    // Check every model that has changed, and refine it.
    //
    // A merged model is two halves glued along a seam that has never been
    // optimized as one, and a single similarity cannot place both halves of a
    // drifted model correctly -- it fits the overlap and misplaces what is far
    // from it. Bundle adjustment cannot walk a misplacement back, so the
    // mapper first asks the correspondence graph whether every image belongs
    // where the model puts it (Mapper::audit) and re-registers the ones that
    // do not. Models nothing has touched keep their exact shape and are
    // skipped, so this costs nothing on a dataset with nothing to fix.
    std::vector<Reconstruction> auditPass(std::vector<Reconstruction> models) {
        for (Reconstruction& m : models) {
            const Signature sig = signature(m);
            if (!audited_.insert(sig).second) continue;
            const size_t pts = m.points3D.size();
            const uint32_t imgs = m.numRegistered();
            Mapper::AuditStats as;
            m = opt_.do_audit ? mapper_.audit(m, &as) : mapper_.refine(m);
            st_.audited_out += as.deregistered;
            st_.audited_repaired += as.reregistered;
            audited_.insert(signature(m));
            if (opt_.verbose && (as.unsupported || m.numRegistered() != imgs))
                fprintf(stderr,
                        "[mgr] audited a model: %u -> %u images, %zu -> %zu points "
                        "(%u unsupported, %u re-registered, %u dropped)\n",
                        imgs, m.numRegistered(), pts, m.points3D.size(), as.unsupported,
                        as.reregistered, as.deregistered);
        }
        return models;
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
            const Signature sig = signature(models[i]);
            if (barren_.count(sig)) continue;
            Mapper::GrowStats gs;
            std::vector<const Reconstruction*> others;
            for (size_t j = 0; j < models.size(); j++)
                if (j != i) others.push_back(&models[j]);
            Reconstruction grown = mapper_.continueFrom(models[i], &gs, others);
            if (!gs.registered) {
                barren_.insert(sig);
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

    // Models that growth has turned into copies of what is already kept.
    //
    // Against the *union* of the models kept before it, not against one of them
    // at a time. Testing one at a time misses the case that actually occurs: on
    // a 5356-image capture the run ended with models of 5249, 5168, 3770 and
    // 1389 images where each of the last three was 98-100% inside the union of
    // the others and none was 90% inside any single one, because their images
    // were split between two bigger models. Largest first, so what a model is
    // measured against is only ever models that outrank it.
    std::vector<Reconstruction> dropRedundant(std::vector<Reconstruction> models) {
        sortModels(models);
        std::vector<char> keep(models.size(), 1);
        std::set<uint32_t> covered;
        for (size_t i = 0; i < models.size(); i++) {
            const uint32_t n = models[i].numRegistered();
            size_t shared = 0;
            for (const auto& kv : models[i].images)
                if (kv.second.registered && covered.count(kv.first)) shared++;
            if (i && n && (double)shared >= opt_.redundant_ratio * (double)n) {
                keep[i] = 0;
                st_.dropped_redundant++;
                if (opt_.verbose)
                    fprintf(stderr,
                            "[mgr] dropped a %u-image model: %zu of its images are already in "
                            "larger ones\n", n, shared);
                continue;
            }
            for (const auto& kv : models[i].images)
                if (kv.second.registered) covered.insert(kv.first);
        }
        std::vector<Reconstruction> out;
        for (size_t i = 0; i < models.size(); i++)
            if (keep[i]) out.push_back(std::move(models[i]));
        return out;
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

    // Movement at the frame corner, in pixels, below which a joint refinement
    // is taken to have changed nothing anyone downstream can act on -- well
    // under the registration and filtering tolerances, which are a few pixels.
    static constexpr double kJointBaMovedPx = 1.0;

    using Signature = std::pair<uint32_t, size_t>;  // registered images, 3D points
    static Signature signature(const Reconstruction& m) {
        return {m.numRegistered(), m.points3D.size()};
    }

    Mapper& mapper_;
    ManagerOptions opt_;
    ManagerStats st_;
    std::set<Signature> barren_;   // models a growth pass could not extend
    std::set<Signature> audited_;  // models already checked in this shape
    std::set<Signature> split_;    // models already examined for a bad seam

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
