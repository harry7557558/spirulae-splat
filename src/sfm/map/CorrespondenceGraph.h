// Correspondence graph: for every (image, feature), the set of features in
// other images it is matched to (COLMAP's scene/correspondence_graph).
// Built from the verified two-view matches; the incremental mapper
// queries it to find 2D-3D correspondences (register-next) and to grow tracks
// (triangulation).
#pragma once

#include <cstdint>
#include <vector>

#include "sfm/core/Matches.h"

namespace sfm {

struct Correspondence {
    uint32_t image_id;
    uint32_t feature_idx;
};

class CorrespondenceGraph {
public:
    // `num_features[i]` = feature count of image i (matching db.images order).
    void build(const MatchesDatabase& db, const std::vector<uint32_t>& num_features) {
        adj_.assign(num_features.size(), {});
        for (size_t i = 0; i < num_features.size(); i++) adj_[i].resize(num_features[i]);
        for (const TwoViewMatches& p : db.pairs) {
            for (const FeatureMatch& m : p.matches) {
                if (m.idx1 >= adj_[p.image1].size() || m.idx2 >= adj_[p.image2].size()) continue;
                adj_[p.image1][m.idx1].push_back({p.image2, m.idx2});
                adj_[p.image2][m.idx2].push_back({p.image1, m.idx1});
            }
        }
    }

    const std::vector<Correspondence>& at(uint32_t image, uint32_t feature) const {
        return adj_[image][feature];
    }

    size_t numImages() const { return adj_.size(); }
    uint32_t numFeatures(uint32_t image) const { return (uint32_t)adj_[image].size(); }

private:
    // adj_[image][feature] -> matched (image, feature) across all pairs
    std::vector<std::vector<std::vector<Correspondence>>> adj_;
};

}  // namespace sfm
