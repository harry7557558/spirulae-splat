#pragma once

// FrameSelect -- motion-blur-aware frame selection for video datasets.
// C++ multithreaded port of scripts/extract_frames.py's FrameSelector: the
// sharpness metric is the variance of the 3x3 Laplacian of the mean-
// subtracted 512x512 grayscale image; within each group of `group`
// consecutive candidate frames the sharpest one is kept.
//
// The GUI extracts candidates with ffmpeg at (target fps x group), then
// calls this to keep the best frame per group -- same output rate as a
// plain fps filter, but each kept frame is the least blurry of its
// neighborhood.

#include <atomic>
#include <functional>
#include <string>

namespace gui {

// Scores every image in cand_dir (sorted by filename) across all hardware
// threads, keeps the sharpest of each consecutive `group`, and MOVES the
// keepers to out_dir as <prefix>NNNNN.<ext> (contiguous numbering). Losers
// are deleted. If keeping one per group still exceeds max_frames, the group
// size is enlarged. Returns the number of frames kept, or -1 on error /
// cancellation.
int select_sharpest_frames(const std::string& cand_dir,
                           const std::string& out_dir,
                           const std::string& prefix,
                           int group, int max_frames,
                           const std::function<void(const std::string&)>& log,
                           const std::atomic<bool>& cancel);

}  // namespace gui
