// GPU pair selection: decide which pairs are worth full matching by matching a
// small top-scale descriptor subset of one image against the other image's
// descriptors, then keeping only each image's best-scoring partners. This is
// the "pair selection" step commercial SfM pipelines run instead of exhaustive
// matching; it fills the roadmap's retrieval slot without a vocabulary tree
// (docs/notes/sfm-design.md D35).
//
// Why this form: exhaustive full matching is quadratic in *pair count* -- at
// the measured 3.5 ms/pair, 1000 images is 500k pairs, ~30 min -- and on
// sparse collections most of those pairs die in verification anyway.
// The score used here is a subsample of the true objective:
// the number of ratio-test survivors when the K largest-scale features of one
// image are matched against the other image's features, computed by the same
// GPU brute-force matcher. Scoring stays O(N^2) but with a constant small
// enough that full matching, now O(N*k), dominates again. A vocabulary tree
// only becomes interesting when even the reduced O(N^2) sweep hurts (~10k
// images); its score is also a strictly weaker pair-relevance signal than
// actually matching descriptors.
//
// The scoring is deliberately asymmetric (mini query vs full train). The
// symmetric mini-vs-mini variant is ~30x cheaper again but was measured and
// rejected: requiring *both* endpoints of a correspondence to fall in their
// image's top-K squares the selection loss, and a 256-descriptor train side
// makes the ratio test toothless (the second-best distance is far), so real
// pairs score low while junk scores nonzero.
#pragma once

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <functional>
#include <utility>
#include <vector>

#include "sfm/core/Features.h"
#include "sfm/feature/Matcher.h"
#include "sfm/feature/Pairing.h"

namespace sfm {

struct PairSelectionOptions {
    uint32_t num_features = 512;  // K: per-image top-scale query subset
    // Train-side cap. Scoring cost is num_features * train_features; the
    // default full-resolution train side is what makes the score reliable.
    // 0 = no cap (use every descriptor).
    uint32_t train_features = 0;
    uint32_t num_neighbors = 32;  // k: keep each image's k best-scoring partners
    uint32_t min_score = 4;       // mini-matches below this never qualify a pair
    // Lowe ratio for the scoring pass. Independent of the full matcher's
    // ratio: scoring only ranks pairs, so it can afford to be looser.
    float ratio = 0.8f;
    // Scoring problems are ~1/32 the size of full matching, so batch more
    // pairs per submit than the full matcher does.
    int batch_pairs = 256;
    int device = -1;
};

// Gather f's K largest-scale features (K = 0 or >= count keeps everything,
// but still gathered in scale order). The canonical feature order is by
// position, deliberately not scale (D16), so this is an explicit host-side
// gather; scale ties break by index.
inline FeatureSet topScaleSubset(const FeatureSet& f, uint32_t K) {
    if (K == 0 || K > f.count()) K = f.count();
    std::vector<uint32_t> idx(f.count());
    for (uint32_t j = 0; j < f.count(); j++) idx[j] = j;
    std::partial_sort(idx.begin(), idx.begin() + K, idx.end(), [&](uint32_t a, uint32_t b) {
        if (f.keypoints[a].scale != f.keypoints[b].scale)
            return f.keypoints[a].scale > f.keypoints[b].scale;
        return a < b;
    });
    FeatureSet m;
    m.width = f.width;
    m.height = f.height;
    m.dim = f.dim;
    m.dtype = f.dtype;
    const uint32_t dsz = f.dim * dtypeSize(f.dtype);
    m.keypoints.reserve(K);
    m.descriptors.resize((size_t)K * dsz);
    for (uint32_t k = 0; k < K; k++) {
        m.keypoints.push_back(f.keypoints[idx[k]]);
        std::memcpy(&m.descriptors[(size_t)k * dsz], &f.descriptors[(size_t)idx[k] * dsz], dsz);
    }
    return m;
}

// Score every i<j pair and return the union of each image's top-k partners,
// in lexicographic pair order. Deterministic: the matcher is exact, subset
// selection and every tie-break are index-ordered. `progress(done, total)` is
// optional.
inline std::vector<std::pair<uint32_t, uint32_t>> prefilterPairs(
    const std::vector<FeatureSet>& feats, const PairSelectionOptions& opt,
    const std::function<void(size_t, size_t)>& progress = nullptr) {
    const uint32_t n = (uint32_t)feats.size();

    // One index space drives the matcher: [0, n) are the query subsets, [n, 2n)
    // the train sides, and every scoring pair is (query i, train j). The train
    // side is the memory cost (the queries are ~32 KB each), so it gets the
    // optional cap -- and when it is uncapped the view points straight at the
    // caller's feature sets instead of copying a gigabyte of descriptors.
    std::vector<FeatureSet> owned(n + (opt.train_features == 0 ? 0 : (size_t)n));
    std::vector<const FeatureSet*> sets(2 * (size_t)n);
    for (uint32_t i = 0; i < n; i++) {
        owned[i] = topScaleSubset(feats[i], opt.num_features);
        sets[i] = &owned[i];
        if (opt.train_features == 0) {
            sets[(size_t)n + i] = &feats[i];
        } else {
            owned[(size_t)n + i] = topScaleSubset(feats[i], opt.train_features);
            sets[(size_t)n + i] = &owned[(size_t)n + i];
        }
    }

    // Every ordered pair is scored -- (i queries, j trains) and the reverse --
    // and a pair's score is the max of its two directions. One direction alone
    // misses pairs whose shared content is top-scale-prominent in only one of
    // the two images. Train-major order: all (i, j) for one train j
    // before moving to the next, so the matcher's resident-descriptor
    // allocator streams each train set through VRAM once instead of thrashing
    // when the dataset outgrows the budget.
    std::vector<std::pair<uint32_t, uint32_t>> all;
    all.reserve((size_t)n * (n - 1));
    for (uint32_t j = 0; j < n; j++)
        for (uint32_t i = 0; i < n; i++)
            if (i != j) all.emplace_back(i, (uint32_t)(n + j));

    MatchOptions mo;
    mo.device = opt.device;
    mo.batch_pairs = opt.batch_pairs;
    mo.max_num_matches = 0;
    mo.max_ratio = opt.ratio;
    // No cross-check: the ratio test against a full-size train side is already
    // selective, and skipping it means the matcher never dispatches the column
    // reduction or downloads the train-side results -- at K=256 queries that
    // cuts the readback per pair from (K + n_train) to K entries, ~33x here.
    mo.cross_check = false;
    BruteForceMatcher matcher(mo);

    // Score on the GPU, chunked. Only the counts are ever wanted, so the
    // matcher counts on the spot rather than building (and immediately
    // discarding) half a million match lists.
    std::vector<uint32_t> score(all.size(), 0);
    // Only a bound on how often progress is reported and how much of `all` the
    // matcher sees at once -- it splits the range into device-sized chunks
    // itself, and a scoring pair's results are small enough that thousands fit
    // in one submit.
    const size_t chunk = (size_t)std::max(1, opt.batch_pairs) * 64;
    std::vector<uint32_t> out;
    for (size_t b = 0; b < all.size(); b += chunk) {
        const size_t e = std::min(b + chunk, all.size());
        matcher.countBatch(sets, all, b, e, out);
        for (size_t p = b; p < e; p++) score[p] = out[p - b];
        if (progress) progress(e, all.size());
    }

    // Fold the two directions into one unordered score. In the train-major
    // order above, (query i, train j) sits at index j*(n-1) + (i < j ? i : i-1).
    auto dirIndex = [&](uint32_t q, uint32_t t) {
        return (size_t)t * (n - 1) + (q < t ? q : q - 1);
    };
    // Union of per-image top-k: a pair survives if either endpoint ranks it.
    std::vector<std::vector<std::pair<uint32_t, uint32_t>>> adj(n);  // (score, partner)
    for (uint32_t i = 0; i < n; i++)
        for (uint32_t j = i + 1; j < n; j++) {
            uint32_t s = std::max(score[dirIndex(i, j)], score[dirIndex(j, i)]);
            if (s < opt.min_score) continue;
            adj[i].push_back({s, j});
            adj[j].push_back({s, i});
        }
    std::vector<std::pair<uint32_t, uint32_t>> sel;
    for (uint32_t i = 0; i < n; i++) {
        auto& a = adj[i];
        const size_t k = std::min<size_t>(opt.num_neighbors, a.size());
        std::partial_sort(a.begin(), a.begin() + k, a.end(),
                          [](const std::pair<uint32_t, uint32_t>& x,
                             const std::pair<uint32_t, uint32_t>& y) {
                              if (x.first != y.first) return x.first > y.first;
                              return x.second < y.second;
                          });
        for (size_t t = 0; t < k; t++) {
            uint32_t j = a[t].second;
            sel.emplace_back(std::min(i, j), std::max(i, j));
        }
    }
    std::sort(sel.begin(), sel.end());
    sel.erase(std::unique(sel.begin(), sel.end()), sel.end());
    return sel;
}

}  // namespace sfm
