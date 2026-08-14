#pragma once

// The reader for what a running reconstruction writes about itself.
//
// The other half of sfm/core/Progress.h: the GUI passes `--progress-dir` to its
// own child and polls the two files that appear there. Both are written whole
// and renamed, so a read either gets the previous snapshot or this one.
//
// Polling rather than watching on purpose -- two files at 2 Hz is a pair of
// stats, and an inotify/ReadDirectoryChangesW pair is a portability problem for
// nothing.

#include "data/DatasetParser.h"

#include <cstdint>
#include <string>
#include <vector>

namespace gui {

// A reconstruction in progress, in the shape ViewportPanel::attach_preview_data
// takes. `ds` carries only what the preview reads -- poses, intrinsics and
// points; there is no dataset behind this and no images to load.
struct LiveModel {
    ParsedDataset ds;
    PostSplitCameras post;
    uint32_t n_images = 0;      // in the capture
    uint32_t n_registered = 0;  // ... and in the model
    uint64_t n_points = 0;      // before subsampling
    bool empty() const { return n_registered == 0 && ds.points.num() == 0; }
};

// Inlier counts between image ranges. `bins` is `n_images` for a capture small
// enough and 512 above that, so a cell is a block of images rather than one.
struct PairMatrix {
    uint32_t n_images = 0, bins = 0;
    std::vector<uint32_t> counts;   // bins*bins, symmetric
    uint32_t peak = 0;
    bool empty() const { return bins == 0 || counts.empty(); }
    uint32_t at(uint32_t r, uint32_t c) const {
        return r < bins && c < bins ? counts[(size_t)r * bins + c] : 0;
    }
};

// `dir` is the --progress-dir the child was given. Both return false when the
// file is absent, unfinished or not newer than `mtime` -- which the caller
// keeps, so a poll that finds nothing new costs one stat.
bool read_live_model(const std::string& dir, int64_t& mtime, LiveModel& out);
bool read_pair_matrix(const std::string& dir, int64_t& mtime, PairMatrix& out);

// The same matrix from a finished `matches.bin`, which is what a run leaves
// behind: the live file is deleted with the rest of the intermediates.
bool read_pair_matrix_from_matches(const std::string& matches_path,
                                   int64_t& mtime, PairMatrix& out);

// Keypoints of one extracted image, in its own pixels, or false when the
// feature file is not there yet. Descriptors are skipped -- they are the file,
// and nothing here looks at them.
//
// `rel_stem` is the image's path under the image directory without its
// extension ("cam0/01973"): the feature files mirror the image tree, so that
// stems from different camera folders cannot collide.
bool read_keypoints(const std::string& features_dir, const std::string& rel_stem,
                    int width, int height, std::vector<float>& xy_out);

}  // namespace gui
