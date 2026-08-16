#pragma once

// The depth-and-normal step of a dataset run: `spirula geometry` over the
// dataset the reconstruction just wrote, or over one that was already there.
//
// A child process for the reasons SfmRunner gives, and a free function rather
// than a class because both runners end with it -- one copy of the argument
// list, so a flag cannot reach only the built-in path.

#include "app/gui/FilmReel.h"
#include "app/gui/PrepProgress.h"
#include "i18n/Message.h"

#include <atomic>
#include <cstdint>
#include <string>
#include <vector>

namespace gui {

struct GeometryJob {
    bool enable = false;
    // A model id (geometry_models()) or a path to an .onnx file.
    std::string model = "metric3d-vit-large";
    int  max_size = 1064;         // longest side of one face the network runs
    bool want_normal = true;
    bool want_depth = false;
    bool normal_jpg = false;
    int  jpeg_quality = 95;
    bool depth_mm = false;
    // 0 auto, 1 yes, 2 no -- `spirula geometry --ray-depth` and `--split`.
    int  ray_depth = 0;
    int  split = 0;
    bool overwrite = false;       // recompute maps that are already on disk
};

// One checkpoint the screen offers, largest last. Index 1 is the default: it
// is what the reference pipeline runs and what the timings in
// src/metric3d/README.md were measured on.
struct GeometryModel {
    const char* id;               // what --model spells
    const spirula::i18n::Msg* label;
    const spirula::i18n::Msg* blurb;
    uint64_t bytes;               // the whole download, both files
};
const std::vector<GeometryModel>& geometry_models();
bool geometry_model_cached(const std::string& id);

// What is left to fetch for `id`, in order; empty when it is all there or
// `id` is a path the user pointed at. vit-giant2 is two files -- the reader
// resolves the second one by name, so both land in the same directory.
struct GeometryDownload {
    std::string url, dest;
    uint64_t bytes = 0;
};
std::vector<GeometryDownload> geometry_model_downloads(const std::string& id);

// "" when this build can estimate geometry, otherwise why it cannot.
std::string geometry_availability();

// Run the step over `dataset`, reporting into `prog` and showing what it
// writes on `reel` (null for no screen). `images` is the folder the maps
// describe, so the reel can show each beside its photograph.
bool run_geometry_step(const GeometryJob& job, const std::string& dataset,
                       const std::string& images, RunProgress& prog,
                       FilmReel* reel, const std::atomic<bool>& cancel,
                       std::string& error);

}  // namespace gui
