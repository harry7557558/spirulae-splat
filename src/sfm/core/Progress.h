#pragma once

// Snapshots of a run in progress, for a front end that is watching it.
//
// Off unless `--progress-dir` was given, so an ordinary CLI run writes nothing
// and pays for nothing. The GUI passes one and polls the files; that is the
// whole of the channel, because the reconstruction is a child process and its
// stdout already belongs to the user (SfmRunner.h says why it still is one).
//
// Every file is written whole to a `.tmp` and renamed, so a reader never sees
// half of one, and every write is rate-limited by wall clock -- the mapper
// registers an image every few milliseconds on a small capture and every few
// seconds on a large one, and a screen wants the same cadence from both.

#include <cstdint>
#include <string>

namespace sfm {

struct Reconstruction;

namespace progress {

// Longest side of the pair matrix, and the point budget of a model snapshot.
// Both are what a screen can show rather than what the run holds: 512^2 u32 is
// 1 MB, and 50k points is already more than a preview resolves.
inline constexpr uint32_t kMatrixBins = 512;
inline constexpr uint32_t kMaxPoints = 50000;

// Where snapshots go. Empty (the default) disables everything below.
void set_dir(const std::string& dir);
bool enabled();

// model.bin: "VKPM", u32 version=2, u32 images, u32 registered, u64 points,
// then per registered image { f32 c2w[12] (OpenGL camera-to-world), u32 w, u32
// h, u32 colmap_model_id, u32 num_params, f64 params[num_params] }, then u32
// count and that many { f32 xyz[3], u8 rgb[3] }.
//
// The model as it stands, subsampled to kMaxPoints. Call it as often as is
// convenient; it returns immediately until the interval has passed, unless
// `force` says this is the last word on a stage.
void model(const Reconstruction& rec, bool force = false);

// Matching is about to verify pairs among `n_images` images.
void begin_matching(uint32_t n_images);
// One verified pair. Safe to call from the verification workers.
void pair(uint32_t image1, uint32_t image2, uint32_t inliers);

// Write whatever is buffered, whatever the clock says. Call at the end of a
// stage so the last state on screen is the final one.
void flush();

}  // namespace progress
}  // namespace sfm
