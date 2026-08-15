// FeatureWatcher.cpp -- see FeatureWatcher.h.

#include "app/gui/FeatureWatcher.h"

#include "app/FrameMask.h"          // app::load_rgb, app::load_stencil
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
                           const std::string& mask_dir,
                           const std::string& features_dir, FilmReel* film) {
    if (_worker.joinable() && _image_dir == image_dir &&
        _mask_dir == mask_dir && _features_dir == features_dir)
        return;
    stop();
    if (image_dir.empty() || features_dir.empty() || !film) return;
    _image_dir = image_dir;
    _mask_dir = mask_dir;
    _features_dir = features_dir;
    _stop = false;
    _worker = std::thread([this, image_dir, mask_dir, features_dir, film] {
        run(image_dir, mask_dir, features_dir, film);
    });
}

void FeatureWatcher::stop() {
    _stop = true;
    if (_worker.joinable()) _worker.join();
    _image_dir.clear();
    _mask_dir.clear();
    _features_dir.clear();
}

void FeatureWatcher::run(std::string image_dir, std::string mask_dir,
                         std::string features_dir, FilmReel* film) {
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

    const fs::path root(image_dir);
    for (const fs::path& f : files) {
        // The feature files mirror the image tree, so a camera folder is part
        // of the name. Lexical: fs::relative resolves symlinks, and an images/
        // of links into a raw capture would relativize right out of the tree.
        fs::path rel = f.lexically_relative(root);
        if (rel.empty() || *rel.begin() == "..") rel = f.filename();
        const std::string rel_stem =
            (rel.parent_path() / rel.stem()).generic_string();
        const std::string feat =
            (fs::path(features_dir) / (rel_stem + ".bin")).string();

        std::vector<KeyPoint2D> kp;
        // Wait for this image's turn rather than skipping ahead: out of order,
        // the reel would show whatever the disk happened to have and none of
        // the sense of the stage working through the capture.
        while (!_stop.load() && !read_keypoints_file(feat, 0, 0, kp))
            std::this_thread::sleep_for(std::chrono::milliseconds(120));
        if (_stop.load()) return;

        FilmFrame frame;
        frame.name = rel.generic_string();
        frame.image_path = f.string();
        frame.points_path = feat;
        if (!mask_dir.empty()) {
            const fs::path m = fs::path(mask_dir) / (rel_stem + ".png");
            if (fs::exists(m, ec)) frame.mask_path = m.string();
        }
        // Every image goes on the reel, so the slider covers the whole
        // capture; only as many as the screen keeps up with are decoded here,
        // and the rest are read back from these paths if it is dragged to them.
        int w = 0, h = 0;
        std::vector<uint8_t> rgb, mask;
        // read_keypoints_file was asked for the file's own coordinates above;
        // now that the size is known, ask again so they land on the picture.
        if (film->wants() && app::load_rgb(f.string(), w, h, rgb) &&
            read_keypoints_file(feat, w, h, kp)) {
            int mw = 0, mh = 0;
            if (!frame.mask_path.empty())
                app::load_stencil(frame.mask_path, mw, mh, mask);
            film->add(frame, rgb.data(), w, h,
                      mask.size() == (size_t)w * h ? mask.data() : nullptr,
                      FramePoints{kp.data(), kp.size()});
        } else {
            film->add(frame);
        }
    }
}

}  // namespace gui
