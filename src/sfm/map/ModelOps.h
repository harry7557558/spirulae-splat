// The vocabulary of operations on a *set* of reconstructions.
//
// A capture that does not come back as one model needs the same five things
// whatever produced the pieces: merge what belongs together, break what does
// not, register the images nobody claimed, drop the copies, and cut a model
// that has written two places on top of each other. The manage loop (D44) drives
// them in rounds over the flat mapper's output; the bottom-up tree (D57) drives
// them level by level over its own. Both need the identical passes with the
// identical thresholds, so they live here and neither owns them.
//
// Everything here is a free function over `Mapper`'s public operations. The
// memo is the one piece of state: passes are skipped for models they have
// already seen in this shape, which is what keeps a loop that runs them
// repeatedly from paying for them repeatedly.
#pragma once

#include <algorithm>
#include <cmath>
#include <numeric>
#include <cstdio>
#include <functional>
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
    // What an audit's repair pass may grow a model to, as a fraction of what
    // it already held (see Mapper::audit). 0 leaves it uncapped, which is what
    // the manage loop wants -- it has no growth schedule of its own to defer
    // to. The bottom-up tree does, and caps it there.
    double audit_growth_frac = 0;
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
    size_t seam_checked = 0, seam_refused = 0, seam_skipped = 0;
    size_t splits = 0, duplicate_splits = 0, split_dropped = 0;
    size_t models_before = 0, models_after = 0;
    size_t covered_before = 0, covered_after = 0;
};

// Models are moved and re-indexed by every pass, so a pass that wants to skip
// what it has already done cannot remember an index. It remembers the shape --
// registered images and 3D points -- which changes whenever anything touched
// the model, and whose collisions cost one skipped attempt that was
// overwhelmingly likely to be a no-op anyway.
struct ModelMemo {
    using Signature = std::pair<uint32_t, size_t>;
    static Signature of(const Reconstruction& m) {
        return {m.numRegistered(), m.points3D.size()};
    }
    std::set<Signature> barren;   // a growth pass could not extend these
    std::set<Signature> audited;  // already checked in this shape
    std::set<Signature> split;    // already examined for a bad seam
    void clear() {
        barren.clear();
        audited.clear();
        split.clear();
    }
};

// Movement at the frame corner, in pixels, below which a joint refinement is
// taken to have changed nothing anyone downstream can act on -- well under the
// registration and filtering tolerances, which are a few pixels.
inline constexpr double kJointBaMovedPx = 1.0;

// How far apart the components' intrinsics are, in pixels at the frame corner.
//
// This is the quantity a joint refinement exists to close (D45): a component
// that fitted its own focal to its own noise. When the components already agree
// to under a pixel there, the solve is a bundle adjustment over the whole
// capture that ends where it started -- and both the manage loop and the
// bottom-up tree would otherwise run one every round.
inline double focalSpreadPx(const std::vector<Reconstruction>& models) {
    std::map<uint32_t, std::pair<double, double>> range;  // id -> (min f, max f)
    std::map<uint32_t, double> radius;                    // ... and its frame corner
    for (const Reconstruction& m : models)
        for (const auto& kv : m.cameras) {
            auto it = range.find(kv.first);
            if (it == range.end()) {
                range[kv.first] = {kv.second.focal(), kv.second.focal()};
                radius[kv.first] = std::hypot(kv.second.cx, kv.second.cy);
            } else {
                it->second.first = std::min(it->second.first, kv.second.focal());
                it->second.second = std::max(it->second.second, kv.second.focal());
            }
        }
    double px = 0;
    for (const auto& kv : range) {
        const double f = 0.5 * (kv.second.first + kv.second.second);
        if (f > 0)
            px = std::max(px, (kv.second.second - kv.second.first) / f * radius[kv.first]);
    }
    return px;
}

// The cross-seam agreement test, as a merge validator (D45).
//
// Every test MergeSession can run by itself is computed from the evidence the
// alignment already used, and repeated structure satisfies all of it -- a fold
// agrees with itself. This one asks the correspondence graph instead: do the
// verified two-view geometries that cross the seam still hold in the merged
// model? It is the only test that separates "these two places look the same"
// from "these two places are the same".
//
// The bottom-up tree needs it as much as the manage loop does. Its levels
// perform hundreds of merges, and without this they run with nothing watching
// for repeated structure at all.
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
        // Nothing to judge on. Counted, because "the test passed" and "the test
        // could not run" are different facts and a merge tree that accumulates
        // the second one is merging unwatched.
        if (sc.tested < min_pairs) {
            if (st) st->seam_skipped++;
            return "";
        }
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

// COLMAP orders written models by 3D point count, descending
// (ReconstructionManager::Write), so sparse/0 is the model with the most
// structure. Every pass that ranks models uses the same order.
inline void sortModels(std::vector<Reconstruction>& models) {
    std::stable_sort(models.begin(), models.end(),
                     [](const Reconstruction& a, const Reconstruction& b) {
                         return a.points3D.size() > b.points3D.size();
                     });
}

// ---- the passes ---------------------------------------------------------

// Models that growth has turned into copies of what is already kept.
//
// Against the *union* of the models kept before it, not against one of them at
// a time. Testing one at a time misses the case that actually occurs: on a
// 5356-image capture the run ended with models of 5249, 5168, 3770 and 1389
// images where each of the last three was 98-100% inside the union of the
// others and none was 90% inside any single one, because their images were
// split between two bigger models. Largest first, so what a model is measured
// against is only ever models that outrank it.
// `follow` arrays are reordered and filtered with the models. The bottom-up
// tree carries a per-model "this one merged" flag across this pass, and losing
// it would make its final audit re-check every model instead of the ones with a
// new seam.
inline std::vector<Reconstruction> dropRedundantModels(std::vector<Reconstruction> models,
                                                       const ManagerOptions& opt,
                                                       ManagerStats& st,
                                                       std::vector<std::vector<char>*> follow = {}) {
    std::vector<size_t> order(models.size());
    std::iota(order.begin(), order.end(), (size_t)0);
    std::stable_sort(order.begin(), order.end(), [&models](size_t a, size_t b) {
        return models[a].points3D.size() > models[b].points3D.size();
    });
    std::vector<size_t> keep;
    std::set<uint32_t> covered;
    for (size_t rank = 0; rank < order.size(); rank++) {
        const Reconstruction& m = models[order[rank]];
        const uint32_t n = m.numRegistered();
        size_t shared = 0;
        for (const auto& kv : m.images)
            if (kv.second.registered && covered.count(kv.first)) shared++;
        if (rank && n && (double)shared >= opt.redundant_ratio * (double)n) {
            st.dropped_redundant++;
            if (opt.verbose)
                fprintf(stderr,
                        "[mgr] dropped a %u-image model: %zu of its images are already in "
                        "larger ones\n", n, shared);
            continue;
        }
        for (const auto& kv : m.images)
            if (kv.second.registered) covered.insert(kv.first);
        keep.push_back(order[rank]);
    }
    for (std::vector<char>* f : follow) {
        std::vector<char> next;
        next.reserve(keep.size());
        for (size_t i : keep) next.push_back(i < f->size() ? (*f)[i] : 0);
        *f = std::move(next);
    }
    std::vector<Reconstruction> out;
    out.reserve(keep.size());
    for (size_t i : keep) out.push_back(std::move(models[i]));
    return out;
}

// Break models the correspondence graph says are not one thing.
//
// The merge pass refuses a bad merge before it commits, but models arrive here
// that were never merged and are wrong anyway -- a chain of registrations
// through repeated structure does the same damage, and D41's own output has
// models like that. Splitting them is what makes the rest of the loop work:
// each piece is a sound reconstruction, and the merge pass gets to try again
// with the seam test watching.
//
// Only models that changed since they were last split are re-examined, so a
// split and a merge cannot chase each other round after round.
inline std::vector<Reconstruction> splitInconsistentModels(Mapper& mapper,
                                                           std::vector<Reconstruction> models,
                                                           const ManagerOptions& opt,
                                                           ModelMemo& memo, ManagerStats& st) {
    std::vector<Reconstruction> out;
    for (Reconstruction& m : models) {
        const ModelMemo::Signature sig = ModelMemo::of(m);
        if (m.numRegistered() < 2 * opt.split_min_group || !memo.split.insert(sig).second) {
            out.push_back(std::move(m));
            continue;
        }
        Mapper::SplitStats ss;
        std::vector<Reconstruction> parts = mapper.splitInconsistent(
            m, opt.merge.max_reproj_error, opt.seam_min_pair_fraction, opt.split_min_matches,
            opt.split_min_group, &ss);
        if (parts.size() <= 1) {
            out.push_back(std::move(m));
            continue;
        }
        st.splits++;
        st.split_dropped += ss.dropped_images;
        if (opt.verbose)
            fprintf(stderr,
                    "[mgr] split a %u-image model into %zu (%zu of %zu inner pairs hold; "
                    "largest group %zu, %zu images dropped)\n",
                    m.numRegistered(), parts.size(), ss.pairs_agree, ss.pairs_tested, ss.largest,
                    ss.dropped_images);
        for (Reconstruction& p : parts) {
            memo.split.insert(ModelMemo::of(p));
            out.push_back(std::move(p));
        }
    }
    return out;
}

// Break a model that has written two places on top of each other. This is the
// failure the agreement split above cannot see -- a fold agrees with itself --
// and it is detected from what the model is *missing*: images it places in the
// same spot looking the same way, with no structure in common
// (findDuplicateStructure).
inline std::vector<Reconstruction> splitFoldedModels(Mapper& mapper,
                                                     std::vector<Reconstruction> models,
                                                     const ManagerOptions& opt, ModelMemo& memo,
                                                     ManagerStats& st) {
    std::vector<Reconstruction> out;
    MatchedFn matched = mapper.matchedPredicate();
    for (Reconstruction& m : models) {
        if (m.numRegistered() < 2 * opt.split_min_group) {
            out.push_back(std::move(m));
            continue;
        }
        DuplicateReport dr = findDuplicateStructure(m, opt.duplicate, matched);
        if (!dr.duplicated(opt.duplicate)) {
            out.push_back(std::move(m));
            continue;
        }
        size_t dropped = 0;
        DuplicateCut cut;
        std::vector<Reconstruction> parts =
            splitDuplicateStructure(m, dr, opt.split_min_group, &dropped, &cut);
        // The conflicts say a fold is possible; the cut says whether the split
        // is cheap. Only a genuine fold is both (D46) -- cutting a sound model
        // always runs through structure it really does share.
        if (parts.size() <= 1 || !foldSplitAccepted(dr, cut, opt.duplicate)) {
            if (opt.verbose && parts.size() > 1)
                fprintf(stderr,
                        "[mgr] a %u-image model has %zu of %zu co-located pairs with nothing "
                        "in common, but splitting it would sever %.1f%% of its co-visibility "
                        "(>%.1f%%): keeping it whole\n",
                        m.numRegistered(), dr.conflicts, dr.colocated, 100.0 * cut.fraction(),
                        100.0 * opt.duplicate.max_cut_fraction);
            out.push_back(std::move(m));
            continue;
        }
        st.duplicate_splits++;
        st.split_dropped += dropped;
        if (opt.verbose) {
            fprintf(stderr,
                    "[mgr] a %u-image model has %zu of %zu co-located image pairs with no "
                    "structure in common and no match either (%zu more share nothing but were "
                    "matched) and a cut that severs %.2f%% of its co-visibility: two places "
                    "written on top of each other. Splitting into",
                    m.numRegistered(), dr.conflicts, dr.colocated, dr.unmatched_but_seen,
                    100.0 * cut.fraction());
            for (const Reconstruction& p : parts) fprintf(stderr, " %u", p.numRegistered());
            if (dropped) fprintf(stderr, " (%zu images dropped)", dropped);
            fprintf(stderr, "\n");
        }
        for (Reconstruction& p : parts) {
            memo.split.insert(ModelMemo::of(p));
            out.push_back(std::move(p));
        }
    }
    return out;
}

// Check every model that has changed, and refine it.
//
// A merged model is two halves glued along a seam that has never been optimized
// as one, and a single similarity cannot place both halves of a drifted model
// correctly -- it fits the overlap and misplaces what is far from it. Bundle
// adjustment cannot walk a misplacement back, so the mapper first asks the
// correspondence graph whether every image belongs where the model puts it
// (Mapper::audit) and re-registers the ones that do not. Models nothing has
// touched keep their exact shape and are skipped, so this costs nothing on a
// dataset with nothing to fix.
inline std::vector<Reconstruction> auditModels(Mapper& mapper, std::vector<Reconstruction> models,
                                               const ManagerOptions& opt, ModelMemo& memo,
                                               ManagerStats& st) {
    for (Reconstruction& m : models) {
        if (!memo.audited.insert(ModelMemo::of(m)).second) continue;
        const size_t pts = m.points3D.size();
        const uint32_t imgs = m.numRegistered();
        Mapper::AuditStats as;
        const uint32_t cap =
            opt.audit_growth_frac > 0
                ? (uint32_t)std::ceil((1.0 + opt.audit_growth_frac) * (double)imgs)
                : 0;
        m = opt.do_audit ? mapper.audit(m, &as, cap) : mapper.refine(m);
        st.audited_out += as.deregistered;
        st.audited_repaired += as.reregistered;
        memo.audited.insert(ModelMemo::of(m));
        if (opt.verbose && (as.unsupported || m.numRegistered() != imgs))
            fprintf(stderr,
                    "[mgr] audited a model: %u -> %u images, %zu -> %zu points "
                    "(%u unsupported, %u re-registered, %u dropped)\n",
                    imgs, m.numRegistered(), pts, m.points3D.size(), as.unsupported,
                    as.reregistered, as.deregistered);
    }
    return models;
}

}  // namespace sfm
