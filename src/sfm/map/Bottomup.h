// Bottom-up reconstruction: many small models, merged upwards (D57).
//
// The incremental mapper plus the manage loop already alternate between
// registering and merging -- the manager's rounds exist precisely because a
// merge makes more registrations possible and a registration makes more merges
// possible. What that alternation never does is start small. It builds one
// large model first and only then discovers what it could not reach, so the
// expensive whole-model passes (global BA, retriangulation, filtering) all run
// at full size, and every repair runs at full size too.
//
// This turns the schedule around. The view graph is cut into *atoms* of a few
// dozen images (sfm/map/Partition.h), each reconstructed by the ordinary
// incremental mapper -- at that size its passes are trivial and its failure
// modes are local. The atoms are then merged into each other, largest overlap
// first, until one model is left or nothing else will merge.
//
// Two things make the merging work where a single pass over the mapper's own
// output does not:
//
//   * **Shared intrinsics throughout.** Every model in flight is bundle-
//     adjusted in one problem with the intrinsics shared per camera group
//     (Mapper::jointRefine). A forty-image atom cannot determine its own focal
//     and must not try; when each model keeps its own answer, the merger is
//     asked to align two reconstructions of the same place in two different
//     gauges, and the pixel-space tests it uses to accept a merge are exactly
//     what that breaks. There is nothing to average at merge time because
//     nothing ever diverged.
//   * **A bundle adjustment between rounds.** Merging chains: A absorbs B,
//     then AB absorbs C across a seam that has never been optimized. Errors
//     accumulate along the chain until an acceptance test refuses, which is
//     what "no two models can be merged" usually means. So a round that
//     merged anything is followed by one joint refinement, and the next round
//     gets to try the refused pairs again on a model that has been made
//     consistent.
//
// The third piece is the hybrid, and it is what stops the tree ending in a pile
// of fragments: **every level grows the models that did not merge, before the
// joint refinement** (detail::growModels). Merging is good at joining models
// that overlap and cannot invent overlap that is not there, so a model that
// failed to merge is short of images, not short of a better merge test --
// registration is the one operation that creates the shared set the next Sim(3)
// aligns on, and images another model already holds are the target rather than
// a side effect.
//
// Growing every level rather than only once the tree has stalled completely is
// what makes it cheap: a stall is reached one refused pair at a time, so
// waiting for it spends a level, and a joint bundle adjustment, per pair. It
// also ramps by itself -- at the first levels almost everything merges and
// there is little to grow, and by the time a few large models are left it is
// nearly all of them.
//
// Growth is audited (Mapper::audit), and that is not optional: unaudited it
// does not merely add bad poses, it makes the cross-seam test agree with the
// merges built on them. Whatever the tree still cannot reach is left to the
// manage loop, which does the same two things without the level structure.
#pragma once

#include <algorithm>
#include <chrono>
#include <numeric>
#include <cstdio>
#include <map>
#include <string>
#include <set>
#include <vector>

#include "sfm/core/Model.h"
#include "sfm/map/Mapper.h"
#include "sfm/map/Merge.h"
#include "sfm/map/Partition.h"

namespace sfm {

struct BottomUpOptions {
    // Atom size, and how many images neighbouring atoms share. The overlap is
    // what a Sim(3) merge aligns on, so it is not optional: two atoms with
    // nothing in common can only be joined by growth, which is the slow path.
    // Overlap below min_part on purpose: it is what keeps a cut part strictly
    // smaller than what it was cut from (see bisect).
    PartitionOptions partition{48, 12, 16};
    // Below this many images the flat mapper is both faster and better: the
    // whole-model passes this exists to shrink are already small, and cutting a
    // capture that reconstructs in one piece only creates seams to repair.
    size_t min_images = 300;
    // Levels of the merge tree. Each halves the model count, so the ceiling is
    // log2(atoms) plus the extra levels a refused pair earns after a
    // refinement; hitting it leaves the rest to the manage loop.
    int max_rounds = 16;
    // One bundle adjustment across every live model with intrinsics shared per
    // camera group, after the atoms and after every level that merged.
    bool joint_intrinsics = true;
    // Bundle-adjustment cadence while an atom is being built. Far looser than
    // the mapper's global default, because a forty-image solve does not fill
    // the hardware: at 1.1 the atom phase ran ~5500 of them and spent half its
    // time there, and what they were protecting against is a bad registration
    // in a model small enough for the joint refinement after the atom phase --
    // one solve over every atom at once, with the intrinsics shared -- to
    // absorb. The per-call overhead is most of a small solve's cost, so the
    // saving is close to the ratio.
    double atom_ba_growth = 2.0;
    // One joint bundle adjustment over every atom, with intrinsics shared per
    // camera group, before any merging. This is the solve the per-atom ones
    // above are traded for: the same work as hundreds of small problems, in one
    // that saturates the device, and it is where the atoms stop each having
    // their own opinion about the focal.
    bool joint_after_atoms = true;
    size_t joint_min_models = 32;
    // Seed attempts an atom may spend looking for further components inside
    // itself. An atom that fragments is usually dust, and the images it leaves
    // behind are picked up by growth later.
    int atom_model_trials = 4;
    // Levels between growth passes over the models that did not merge; 0
    // disables in-tree growth and hands every stall to the manage loop, which
    // does the same thing with the same audit. `max_grow_passes` bounds the
    // total so a large tree cannot spend all its time registering.
    //
    // Growth *must* be followed by the audit growModels runs: unaudited, it
    // does not merely add bad poses, it makes the cross-seam test agree with
    // the merges built on them (measured: 21 points of AUC@5, and the same
    // whole-session displacement the flat mapper produces).
    int grow_every = 1;
    int max_grow_passes = 8;
    // What one growth pass may add to a model, as a fraction of what it already
    // holds and never fewer than `grow_budget_min`. The pass is asked for a
    // bridge to its neighbour, not for a reconstruction: the overlap bound in
    // continueFrom limits only images another model holds, so without this a
    // model beside a large uncovered region grows into all of it.
    double grow_budget_frac = 0.25;
    size_t grow_budget_min = 25;
    bool verbose = true;
};

struct BottomUpStats {
    size_t atoms = 0;
    size_t atom_images = 0;        // summed over atoms, so overlap counts twice
    size_t models_from_atoms = 0;
    size_t rounds = 0;
    size_t merges = 0, merges_refused = 0;
    size_t joint_ba = 0;
    size_t grown_images = 0;       // registered by a growth pass
    size_t grown_rejected = 0;     // ... and dropped again by the pose check
    double t_atoms = 0, t_merge = 0, t_ba = 0, t_grow = 0;
};

namespace detail {

// One level of the tree: every model absorbs at most one other, best overlap
// first and smallest pair first, so the models grow *together*.
//
// The alternative -- MergeSession::runAuto, which chains until nothing else
// merges -- is both slower and worse here. Slower because every attempt copies
// the anchor (the merge has no cheaper inverse), and an anchor that absorbs a
// hundred atoms in turn is copied at every size it passes through: quadratic,
// where a matching copies the smaller half and doubles the model size per
// level. Worse because each link in the chain is a seam that has never been
// optimized, and the ones at the end of the chain are asked to align across
// all of them at once.
//
// Carries a per-model "has absorbed something" flag through the reindexing, so
// the caller can refine only what changed.
inline size_t mergeLevel(std::vector<Reconstruction>& models, const MergeOptions& opt,
                         std::vector<char>& dirty, std::vector<char>& stalled, size_t& refused,
                         std::map<std::string, size_t>* why = nullptr) {
    MergeSession s(std::move(models), opt);
    std::vector<MergeCandidate> cands = s.candidates();
    std::stable_sort(cands.begin(), cands.end(),
                     [&s](const MergeCandidate& a, const MergeCandidate& b) {
                         if (a.common_images != b.common_images)
                             return a.common_images > b.common_images;
                         const uint32_t sa = s.model(a.dst).numRegistered() +
                                             s.model(a.src).numRegistered();
                         const uint32_t sb = s.model(b.dst).numRegistered() +
                                             s.model(b.src).numRegistered();
                         return sa < sb;
                     });
    std::vector<char> busy(s.numModels(), 0);
    size_t merges = 0;
    for (const MergeCandidate& c : cands) {
        if (busy[c.dst] || busy[c.src]) continue;
        // A refusal does not retire either model for the level: it may still
        // have a partner it can align with.
        const MergeAttempt a = s.tryMerge(c.dst, c.src);
        if (a.merged) {
            merges++;
            busy[c.dst] = busy[c.src] = 1;
        } else {
            refused++;
            // Reasons carry counts, which would make every one unique; keep the
            // leading words so a level's refusals group into kinds.
            if (why) {
                std::string k = a.reason.substr(0, a.reason.find_first_of("0123456789"));
                while (!k.empty() && (k.back() == ' ' || k.back() == '(')) k.pop_back();
                (*why)[k.empty() ? a.reason : k]++;
            }
        }
    }
    const std::set<size_t> touched = s.changed();
    const std::vector<char> was = dirty;
    models.clear();
    dirty.clear();
    stalled.clear();
    for (size_t i = 0; i < s.numModels(); i++) {
        if (!s.alive(i)) continue;
        stalled.push_back(busy[i] ? 0 : 1);  // survived the level without merging
        models.push_back(std::move(s.modelMut(i)));
        dirty.push_back((i < was.size() && was[i]) || touched.count(i) ? 1 : 0);
    }
    return merges;
}

// Register whatever each model can still take, largest first.
//
// This is the incremental half of the hybrid, and it is what a stalled level
// needs. A level stops for one of two reasons and growth answers both: either
// too few pairs still share `min_common_images`, in which case the models are
// not wrong but simply do not touch, and only registering more images can make
// them touch; or a merge was refused on thin evidence, in which case more
// shared images give the alignment more to fit and the cross-seam test more to
// judge on -- so a correct merge that was refused gets another chance, and a
// wrong one is refused with more confidence.
//
// An image another model already holds is a legitimate target, and in fact the
// point: that overlap *is* what the next Sim(3) aligns on. `continueFrom`
// bounds it to max_model_overlap of them, which is enough to determine a
// transform without letting one model quietly absorb its neighbour.
inline size_t growModels(Mapper& mapper, std::vector<Reconstruction>& models,
                         std::vector<char>& dirty, const std::vector<char>& which,
                         double budget_frac, size_t budget_min, size_t& rejected) {
    std::vector<size_t> order(models.size());
    std::iota(order.begin(), order.end(), (size_t)0);
    std::stable_sort(order.begin(), order.end(), [&models](size_t a, size_t b) {
        return models[a].numRegistered() > models[b].numRegistered();
    });
    size_t registered = 0;
    for (size_t i : order) {
        if (i < which.size() && !which[i]) continue;
        const uint32_t have = models[i].numRegistered();
        if (have < 2) continue;
        std::vector<const Reconstruction*> others;
        for (size_t j = 0; j < models.size(); j++)
            if (j != i) others.push_back(&models[j]);
        // Bridge, not reconstruct. The pass is asked for enough images to make
        // a merge possible, scaled to the model so a big one may reach further
        // than a small one, and it is stopped well before it could grow into a
        // whole uncovered region on its own -- which would be the flat mapper
        // again, at the flat mapper's cost.
        const uint32_t budget =
            (uint32_t)std::max((double)budget_min, budget_frac * (double)have);
        Mapper::GrowStats gs;
        uint32_t bad = 0;
        Reconstruction grown =
            mapper.growByPnP(models[i], &gs, others, have + budget, &bad);
        rejected += bad;
        if (!gs.registered) continue;
        registered += gs.registered;
        models[i] = std::move(grown);
        dirty[i] = 1;
    }
    return registered;
}

}  // namespace detail

// Reconstruct bottom-up. `mapper` must already be set up for the whole
// database; it is restricted per atom and released before returning.
inline std::vector<Reconstruction> bottomUpReconstruct(Mapper& mapper, const MatchesDatabase& db,
                                                       const BottomUpOptions& opt,
                                                       const MergeOptions& merge_opt,
                                                       BottomUpStats& st) {
    auto clk = [] { return std::chrono::steady_clock::now(); };
    auto secs = [](auto a, auto b) { return std::chrono::duration<double>(b - a).count(); };

    // The starting intrinsics are chosen once, over the whole database, and
    // every atom inherits them. Left to the atoms, the first one built would
    // pick for all the others from a few dozen images (D48).
    mapper.bootstrapCameras();

    ViewGraph g = buildViewGraph(db);
    std::vector<std::vector<uint32_t>> atoms = partitionViewGraph(g, opt.partition);
    st.atoms = atoms.size();
    for (const std::vector<uint32_t>& a : atoms) st.atom_images += a.size();
    if (opt.verbose) {
        size_t small = SIZE_MAX, big = 0;
        for (const std::vector<uint32_t>& a : atoms) {
            small = std::min(small, a.size());
            big = std::max(big, a.size());
        }
        fprintf(stderr, "[bup] %zu image(s) -> %zu atom(s) of %zu..%zu images (%.2fx cover)\n",
                db.images.size(), atoms.size(), atoms.empty() ? 0 : small, big,
                db.images.empty() ? 0.0 : (double)st.atom_images / (double)db.images.size());
    }

    // A capture that does not split into at least two atoms has nothing to
    // merge, and the flat mapper is what a single atom would have run anyway.
    if (atoms.size() < 2) {
        if (opt.verbose) fprintf(stderr, "[bup] one atom: reconstructing it flat\n");
        return mapper.run();
    }

    // ---- the atoms -------------------------------------------------------
    MapperOptions saved = mapper.options();
    mapper.options().ba_growth_ratio = std::max(1.0 + 1e-9, opt.atom_ba_growth);
    mapper.options().max_model_trials = opt.atom_model_trials;
    mapper.options().verbose = false;
    std::vector<Reconstruction> models;
    auto t0 = clk();
    for (size_t i = 0; i < atoms.size(); i++) {
        mapper.restrictTo(atoms[i]);
        mapper.claimAll({});
        uint32_t reg = 0;
        size_t kept = 0;
        for (Reconstruction& m : mapper.run()) {
            if (m.numRegistered() < 2) continue;
            reg += m.numRegistered();
            kept++;
            models.push_back(std::move(m));
        }
        if (opt.verbose)
            fprintf(stderr, "[bup] atom %zu/%zu: %zu images -> %zu model(s), %u registered\n",
                    i + 1, atoms.size(), atoms[i].size(), kept, reg);
    }
    mapper.restrictTo({});
    mapper.claimAll({});
    mapper.options() = saved;
    st.t_atoms = secs(t0, clk());
    st.models_from_atoms = models.size();
    // Every atom failed. The flat mapper will fail the same way and fail fast,
    // and it is what gives the caller a model to report the failure on.
    if (models.empty()) return mapper.run();

    // The solve the cheap per-atom cadence is traded for: every atom in one
    // problem, intrinsics shared per camera group. Hundreds of forty-image
    // solves do not fill the device; this is the same work in one that does,
    // and it is where the atoms stop each having their own focal.
    // Only when there are enough atoms for it to be the trade it is meant to
    // be. With a dozen atoms it is one extra full solve buying back a handful
    // of tiny ones -- measured, that cost 19 % on a 480-image capture.
    if (opt.joint_after_atoms && models.size() >= opt.joint_min_models) {
        t0 = clk();
        mapper.jointRefine(models);
        st.joint_ba++;
        st.t_ba += secs(t0, clk());
        if (opt.verbose)
            fprintf(stderr, "[bup] joint refinement over %zu atom model(s): %.1f s\n",
                    models.size(), st.t_ba);
    }

    // ---- upwards ---------------------------------------------------------
    std::vector<char> dirty(models.size(), 0);    // models that absorbed another
    std::vector<char> stalled(models.size(), 1);  // ... and those that did not
    int grow_passes = 0;
    for (int round = 0; round < std::max(1, opt.max_rounds); round++) {
        st.rounds = (size_t)round + 1;
        size_t refused = 0;
        std::map<std::string, size_t> why;
        t0 = clk();
        const size_t merges =
            detail::mergeLevel(models, merge_opt, dirty, stalled, refused, &why);
        st.t_merge += secs(t0, clk());
        st.merges += merges;
        st.merges_refused += refused;
        if (opt.verbose) {
            fprintf(stderr, "[bup] level %zu: %zu merge(s), %zu refused, %zu model(s) left\n",
                    st.rounds, merges, refused, models.size());
            // Hundreds of merges a level, so the reasons are summarized by kind
            // rather than printed one by one -- but they have to be visible, or
            // a level that refuses everything looks the same as one with
            // nothing left to merge.
            for (const auto& kv : why)
                fprintf(stderr, "[bup]   %4zu x %s\n", kv.second, kv.first.c_str());
        }
        if (models.size() <= 1) break;

        // Grow the models that did *not* merge, before the joint refinement.
        //
        // A model that failed to merge is short of overlap, and registration is
        // the only thing that creates overlap. Waiting for a level to stall
        // completely before growing spends a level per refused pair discovering
        // that one at a time; growing every level lets the next one merge
        // instead. Models that just merged are left alone -- they have a fresh
        // seam to settle and already have work at the next level.
        //
        // This ramps by itself: at the first levels almost everything merges so
        // there is little to grow, and by the time the tree is down to a few
        // large models it is nearly all of them.
        size_t reg = 0;
        if (opt.grow_every > 0 && grow_passes < opt.max_grow_passes &&
            (round % opt.grow_every) == 0) {
            size_t want = 0;
            for (char c : stalled) want += c ? 1 : 0;
            if (want) {
                t0 = clk();
                size_t rejected = 0;
                reg = detail::growModels(mapper, models, dirty, stalled, opt.grow_budget_frac,
                                         opt.grow_budget_min, rejected);
                st.grown_rejected += rejected;
                st.t_grow += secs(t0, clk());
                st.grown_images += reg;
                if (reg) grow_passes++;
                if (opt.verbose)
                    fprintf(stderr, "[bup]   growth: %zu image(s) into %zu of %zu model(s) that "
                            "did not merge (%zu rejected by the pose check)\n", reg, want,
                            models.size(), st.grown_rejected);
            }
        }
        // Nothing merged and nothing grew: the tree is finished, and what is
        // left is for the manage loop.
        if (merges == 0 && reg == 0) break;

        t0 = clk();
        if (opt.joint_intrinsics && models.size() > 1) {
            mapper.jointRefine(models);
            st.joint_ba++;
        } else {
            for (size_t i = 0; i < models.size(); i++)
                if (dirty[i] && models[i].numRegistered() >= 2)
                    models[i] = mapper.refine(models[i]);
        }
        st.t_ba += secs(t0, clk());
    }

    // The joint pass optimizes but does not filter: it is one bundle
    // adjustment, with none of the mapper's retriangulation, track completion
    // or de-registration. Every model that absorbed another gets that treatment
    // once, here, before the manage loop builds anything else on it. An atom
    // that never merged already had it, from the run that built it.
    t0 = clk();
    for (size_t i = 0; i < models.size(); i++)
        if (dirty[i] && models[i].numRegistered() >= 2) models[i] = mapper.refine(models[i]);
    st.t_ba += secs(t0, clk());

    std::stable_sort(models.begin(), models.end(),
                     [](const Reconstruction& a, const Reconstruction& b) {
                         return a.points3D.size() > b.points3D.size();
                     });
    if (opt.verbose)
        fprintf(stderr,
                "[bup] %zu atom model(s) -> %zu after %zu level(s): %zu merged, %zu refused, "
                "%zu grown (atoms %.1f s, merge %.1f s, grow %.1f s, BA %.1f s)\n",
                st.models_from_atoms, models.size(), st.rounds, st.merges, st.merges_refused,
                st.grown_images, st.t_atoms, st.t_merge, st.t_grow, st.t_ba);
    mapper.claimAll(models);
    return models;
}

}  // namespace sfm
