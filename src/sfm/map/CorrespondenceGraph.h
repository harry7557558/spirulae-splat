// Correspondence graph: for every (image, feature), the set of features in
// other images it is matched to (COLMAP's scene/correspondence_graph).
// Built from the verified two-view matches; the incremental mapper
// queries it to find 2D-3D correspondences (register-next) and to grow tracks
// (triangulation).
//
// Stored per image as CSR -- one offsets array over the image's features plus
// one flat correspondence array -- rather than a vector per feature. A large
// capture has millions of features and tens of millions of correspondences, so
// the vector-per-feature layout spent a heap allocation and 24 bytes of header
// on each of them, and made the mapper's innermost loops chase a pointer per
// feature. `at()` returns a view over the flat array; every caller iterates.
#pragma once

#include <cstdint>
#include <vector>

#include "sfm/core/Matches.h"

namespace sfm {

struct Correspondence {
    uint32_t image_id;
    uint32_t feature_idx;
};

// Contiguous run of correspondences for one (image, feature).
struct CorrespondenceView {
    const Correspondence* first = nullptr;
    const Correspondence* last = nullptr;
    const Correspondence* begin() const { return first; }
    const Correspondence* end() const { return last; }
    size_t size() const { return (size_t)(last - first); }
    bool empty() const { return first == last; }
};

class CorrespondenceGraph {
public:
    // `num_features[i]` = feature count of image i (matching db.images order).
    void build(const MatchesDatabase& db, const std::vector<uint32_t>& num_features) {
        const size_t n = num_features.size();
        starts_.assign(n, {});
        data_.assign(n, {});
        for (size_t i = 0; i < n; i++) starts_[i].assign((size_t)num_features[i] + 1, 0);

        // Count, prefix-sum, scatter. The +1 offset in the counting pass makes
        // the prefix sum land directly in the final offsets array.
        for (const TwoViewMatches& p : db.pairs)
            for (const FeatureMatch& m : p.matches) {
                if (m.idx1 + 1 >= starts_[p.image1].size() ||
                    m.idx2 + 1 >= starts_[p.image2].size())
                    continue;
                starts_[p.image1][m.idx1 + 1]++;
                starts_[p.image2][m.idx2 + 1]++;
            }
        std::vector<std::vector<uint32_t>> fill(n);
        for (size_t i = 0; i < n; i++) {
            std::vector<uint32_t>& s = starts_[i];
            for (size_t f = 1; f < s.size(); f++) s[f] += s[f - 1];
            data_[i].resize(s.empty() ? 0 : s.back());
            fill[i].assign(s.begin(), s.end() - (s.empty() ? 0 : 1));
        }
        for (const TwoViewMatches& p : db.pairs)
            for (const FeatureMatch& m : p.matches) {
                if (m.idx1 + 1 >= starts_[p.image1].size() ||
                    m.idx2 + 1 >= starts_[p.image2].size())
                    continue;
                data_[p.image1][fill[p.image1][m.idx1]++] = {p.image2, m.idx2};
                data_[p.image2][fill[p.image2][m.idx2]++] = {p.image1, m.idx1};
            }
    }

    CorrespondenceView at(uint32_t image, uint32_t feature) const {
        const std::vector<uint32_t>& s = starts_[image];
        const Correspondence* base = data_[image].data();
        return {base + s[feature], base + s[feature + 1]};
    }

    size_t numImages() const { return starts_.size(); }
    uint32_t numFeatures(uint32_t image) const {
        return starts_[image].empty() ? 0u : (uint32_t)(starts_[image].size() - 1);
    }

private:
    // Per image: offsets over the image's features (num_features+1 entries) and
    // the correspondences they index into.
    std::vector<std::vector<uint32_t>> starts_;
    std::vector<std::vector<Correspondence>> data_;
};

}  // namespace sfm
