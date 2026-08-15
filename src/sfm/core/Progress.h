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
#include <functional>
#include <string>
#include <utility>
#include <vector>

namespace sfm {

struct Point3D;
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
//
// `color` fills one point's rgb, and is called only for the points a snapshot
// actually writes: the mapper colours a model when it finishes one, so without
// it every preview of a model under construction is neutral grey.
using PointColor = std::function<void(const Point3D&, uint8_t rgb[3])>;
void model(const Reconstruction& rec, bool force = false,
           const PointColor& color = {});

// pairs.bin: "VKPP", u32 version=2, u32 images, u32 bins, then three
// bins*bins u32 planes -- summed inliers, candidate pairs, verified pairs.
// The last two are what tell a cell nothing has reached it yet from a cell
// pairing was never going to try.
//
// Matching is about to verify `pairs` among `n_images` images.
void begin_matching(uint32_t n_images,
                    const std::vector<std::pair<uint32_t, uint32_t>>& pairs);
// One verified pair, `inliers` of 0 meaning it did not survive verification.
// Safe to call from the verification workers.
void pair(uint32_t image1, uint32_t image2, uint32_t inliers);

// Write whatever is buffered, whatever the clock says. Call at the end of a
// stage so the last state on screen is the final one.
void flush();

}  // namespace progress
}  // namespace sfm
