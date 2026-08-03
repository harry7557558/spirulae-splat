// Which image pairs to match. Kept separate from the matcher so the pairing
// strategy and the descriptor-matching algorithm compose independently: any
// IFeatureMatcher can be driven over any pair list.
//
// Exhaustive (all i<j) and sequential (a sliding window over a video) are
// generated here; Prefilter (GPU pair selection by mini-matching top-scale
// subsets) needs the features and lives in sfm/feature/PairSelection.h.
//
// Sequential is a chain and nothing more: it has no link between two passes
// over the same place, which is what `--loop-closure` (on by default) adds by
// unioning the Prefilter shortlist into the window -- see matchFeatureDir.
#pragma once

#include <cstdint>
#include <utility>
#include <vector>

namespace sfm {

enum class PairMode { Exhaustive, Sequential, Prefilter };

// Generate ordered image pairs (i < j) for `n` images.
//   Exhaustive: every pair.
//   Sequential: each image with the next `overlap` images (video-style);
//     `overlap` <= 0 falls back to 10 (COLMAP's default).
inline std::vector<std::pair<uint32_t, uint32_t>> generatePairs(uint32_t n, PairMode mode,
                                                                int overlap = 10) {
    std::vector<std::pair<uint32_t, uint32_t>> pairs;
    if (mode == PairMode::Exhaustive) {
        pairs.reserve((size_t)n * (n - 1) / 2);
        for (uint32_t i = 0; i < n; i++)
            for (uint32_t j = i + 1; j < n; j++) pairs.emplace_back(i, j);
    } else {  // Sequential
        int ov = overlap > 0 ? overlap : 10;
        for (uint32_t i = 0; i < n; i++)
            for (uint32_t k = 1; k <= (uint32_t)ov && i + k < n; k++)
                pairs.emplace_back(i, i + k);
    }
    return pairs;
}

}  // namespace sfm
