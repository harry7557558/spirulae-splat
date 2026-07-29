// Bottom-up hierarchical mapping (src/sfm/README.md).
//
// The incremental mapper is the right algorithm and the wrong schedule for a
// large capture: every global bundle adjustment, retriangulation and filter
// pass walks the whole model, and it runs dozens of them, so the last hundred
// images cost far more than the first hundred. Reconstructing overlapping
// *clusters* of the view graph and gluing them together pays the same total
// number of registrations at a fraction of the model size, and the glue is
// already here -- the clusters are cut to overlap, so the Sim(3) merge (D43)
// always has shared poses to align on.
//
// Shape:
//   partition (sfm/map/Partition.h)  ->  one model per leaf, built by the same
//   Mapper restricted to that leaf   ->  pairwise merge, largest first
//   ->  the manager's usual merge / grow / audit rounds on the result.
//
// What it is *not*: a different reconstruction algorithm. Every model here is
// built by Mapper::run's own rules, and everything after the leaves is
// MergeSession and ModelManager. That is deliberate -- the hierarchy is a
// schedule, so a regression in it cannot be a regression in the geometry (D55).
#pragma once

#include <algorithm>
#include <cstdio>
#include <set>
#include <vector>

#include "sfm/core/Model.h"
#include "sfm/map/Manager.h"
#include "sfm/map/Mapper.h"
#include "sfm/map/Merge.h"
#include "sfm/map/Partition.h"

namespace sfm {

// Measured on three ~1100-1200 image captures, flat schedule first:
//
//   vicon_room    1003 s / 1100 images / 89.9 AUC@5  ->  761 s / 1120 / 93.6
//   windmill       981 s / 1194        / 92.7        -> 1144 s / 1194 / 90.5
//   zipnerf_nyc    246 s / 1001        / 96.7        ->  281 s /  998 / 96.5
//
// Structurally it does what it is for: 12-14 clusters of 77-189 images, every
// cluster reconstructed, every leaf merged and none refused. What it does *not*
// do is make an already-easy capture cheaper -- a flat run's bundle adjustments
// are mostly small (the model grows from two images) and its cost sits in the
// last few full-size ones, which the hierarchical run pays as well, on top of
// reconstructing the clusters' overlap twice. Where it pays is the capture the
// flat schedule finds hard: vicon_room is 24% faster, 20 images fuller and 3.7
// AUC better, because each cluster is a problem small enough to get right.
//
// So it is an option and not the default, and the two things that would make it
// the default are reconstructing clusters in parallel (a `rec_` per worker) and
// a capture size where the flat model's tail dominates outright.
struct HierarchicalOptions {
    PartitionOptions partition;
    // Below this many images the flat incremental mapper is both faster and
    // better: clustering buys nothing when the model never gets big.
    size_t min_images = 400;
    bool verbose = true;
};

struct HierarchicalStats {
    size_t clusters = 0;
    size_t cluster_images = 0;    // summed over clusters, overlap counted twice
    size_t models_from_leaves = 0;
    size_t merges = 0, merges_refused = 0;
    double seconds_leaves = 0, seconds_merge = 0;
};

// Reconstruct by clusters and merge. `mapper` must already be set up for the
// whole database; it is restricted per cluster and released at the end.
inline std::vector<Reconstruction> hierarchicalReconstruct(Mapper& mapper,
                                                           const MatchesDatabase& db,
                                                           const HierarchicalOptions& opt,
                                                           const MergeOptions& merge_opt,
                                                           HierarchicalStats& st) {
    ViewGraph g = buildViewGraph(db);
    std::vector<std::vector<uint32_t>> clusters = partitionViewGraph(g, opt.partition);
    st.clusters = clusters.size();
    for (const std::vector<uint32_t>& c : clusters) st.cluster_images += c.size();
    if (opt.verbose) {
        fprintf(stderr, "[hier] %zu image(s) -> %zu cluster(s) of", db.images.size(),
                clusters.size());
        for (size_t i = 0; i < clusters.size() && i < 12; i++)
            fprintf(stderr, " %zu", clusters[i].size());
        if (clusters.size() > 12) fprintf(stderr, " ...");
        fprintf(stderr, " images\n");
    }

    // ---- leaves ----------------------------------------------------------
    std::vector<Reconstruction> models;
    for (size_t i = 0; i < clusters.size(); i++) {
        mapper.restrictTo(clusters[i]);
        mapper.claimAll({});
        std::vector<Reconstruction> part = mapper.run();
        size_t kept = 0;
        for (Reconstruction& m : part) {
            if (m.numRegistered() < 2) continue;
            kept++;
            models.push_back(std::move(m));
        }
        if (opt.verbose)
            fprintf(stderr, "[hier] cluster %zu/%zu (%zu images) -> %zu model(s)\n", i + 1,
                    clusters.size(), clusters[i].size(), kept);
    }
    mapper.restrictTo({});
    mapper.claimAll({});
    st.models_from_leaves = models.size();

    // ---- glue ------------------------------------------------------------
    // Straight into the merger: sibling clusters share `overlap` images by
    // construction, which is exactly what alignReconstructions needs. The
    // manager's rounds run afterwards, on whatever is left over.
    if (models.size() > 1) {
        MergeSession s(std::move(models), merge_opt);
        st.merges = s.runAuto();
        for (const MergeAttempt& a : s.log())
            if (!a.merged) st.merges_refused++;
        models = s.take();
    }
    if (opt.verbose)
        fprintf(stderr, "[hier] %zu leaf model(s) -> %zu after merging (%zu merge(s), %zu "
                "refused)\n", st.models_from_leaves, models.size(), st.merges,
                st.merges_refused);
    return models;
}

}  // namespace sfm
