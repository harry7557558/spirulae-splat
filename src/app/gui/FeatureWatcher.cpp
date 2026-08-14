// FeatureWatcher.cpp -- see FeatureWatcher.h.

#include "app/gui/FeatureWatcher.h"

#include "app/FrameMask.h"          // app::load_rgb
#include "app/gui/SfmProgress.h"

#include <algorithm>
#include <chrono>
#include <filesystem>

namespace fs = std::filesystem;

namespace gui {

namespace {

bool is_image(const fs::path& p) {
    std::string e = p.extension().string();
    for (char& c : e) c = (char)std::tolower((unsigned char)c);
    return e == ".jpg" || e == ".jpeg" || e == ".png" || e == ".webp" ||
           e == ".tif" || e == ".tiff" || e == ".bmp";
}

}  // namespace

FeatureWatcher::~FeatureWatcher() { stop(); }

void FeatureWatcher::start(const std::string& image_dir,
                           const std::string& features_dir, FilmStrip* film) {
    if (_worker.joinable() && _image_dir == image_dir &&
        _features_dir == features_dir)
        return;
    stop();
    if (image_dir.empty() || features_dir.empty() || !film) return;
    _image_dir = image_dir;
    _features_dir = features_dir;
    _stop = false;
    _worker = std::thread([this, image_dir, features_dir, film] {
        run(image_dir, features_dir, film);
    });
}

void FeatureWatcher::stop() {
    _stop = true;
    if (_worker.joinable()) _worker.join();
    _image_dir.clear();
    _features_dir.clear();
}

void FeatureWatcher::run(std::string image_dir, std::string features_dir,
                         FilmStrip* film) {
    std::error_code ec;
    std::vector<fs::path> files;
    for (fs::recursive_directory_iterator
             it(image_dir, fs::directory_options::skip_permission_denied |
                               fs::directory_options::follow_directory_symlink, ec),
         end;
         !ec && it != end; it.increment(ec))
        if (it->is_regular_file(ec) && is_image(it->path()))
            files.push_back(it->path());
    std::sort(files.begin(), files.end());

    // A fixed number of images spread across the capture, rather than whatever
    // a wall-clock rate limit happened to catch: forty images and four thousand
    // should both end up showing the same span of the dataset.
    const size_t stride = files.size() / (size_t)kFilmSlots + 1;
    const fs::path root(image_dir);
    size_t index = 0;
    for (const fs::path& f : files) {
        const bool show = index++ % stride == 0;
        // The feature files mirror the image tree, so a camera folder is part
        // of the name. Lexical: fs::relative resolves symlinks, and an images/
        // of links into a raw capture would relativize right out of the tree.
        fs::path rel = f.lexically_relative(root);
        if (rel.empty() || *rel.begin() == "..") rel = f.filename();
        const std::string rel_stem =
            (rel.parent_path() / rel.stem()).generic_string();

        std::vector<float> xy;
        // Wait for this image's turn rather than skipping ahead: out of order,
        // the strip would show whatever the disk happened to have and none of
        // the sense of the stage working through the capture.
        while (!_stop.load() && !read_keypoints(features_dir, rel_stem, 0, 0, xy))
            std::this_thread::sleep_for(std::chrono::milliseconds(120));
        if (_stop.load()) return;
        if (!show) continue;

        int w = 0, h = 0;
        std::vector<uint8_t> rgb;
        if (!app::load_rgb(f.string(), w, h, rgb)) continue;
        // read_keypoints was asked for the file's own coordinates above (0, 0);
        // now that the size is known, ask again so they land on the picture.
        if (!read_keypoints(features_dir, rel_stem, w, h, xy)) continue;
        film->offer(rgb.data(), w, h, rel.generic_string(), nullptr,
                    FramePoints{xy.data(), xy.size() / 2});
    }
}

}  // namespace gui
