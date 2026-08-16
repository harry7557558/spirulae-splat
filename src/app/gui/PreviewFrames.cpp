// PreviewFrames.cpp -- see PreviewFrames.h.

#include "app/gui/PreviewFrames.h"

#include "app/FrameMask.h"
#include "app/gui/DatasetPrep.h"
#include "app/gui/Subprocess.h"
#include "i18n/catalog/Dataset.h"

#ifdef SS_HAVE_VIDEO
#include "video/Video.h"
#endif

#include <algorithm>
#include <cctype>
#include <chrono>
#include <filesystem>

namespace fs = std::filesystem;
namespace dmsg = spirula::i18n::msg::dataset;

namespace gui {

namespace {

bool is_image_file(const fs::path& p) {
    std::string e = p.extension().string();
    for (auto& c : e) c = (char)std::tolower((unsigned char)c);
    return e == ".jpg" || e == ".jpeg" || e == ".png" || e == ".webp" ||
           e == ".tif" || e == ".tiff" || e == ".bmp";
}

fs::path temp_still(const char* tag) {
    return fs::temp_directory_path() /
           (std::string("spirula-") + tag + "-" +
            std::to_string(
                std::chrono::steady_clock::now().time_since_epoch().count()) +
            ".jpg");
}

}  // namespace

void collect_preview_frames(PreviewSource& src, int offers,
                            std::vector<PreviewFrame>& frames,
                            std::vector<std::string>& all_files,
                            const std::atomic<bool>& cancel) {
    frames.clear();
    all_files.clear();
    src.video_seconds = 0.0;
    offers = std::max(1, offers);
    std::error_code ec;

    if (src.is_video) {
        // Frames are decoded on demand; the list is just how many offers the
        // slider makes. Seeking is not supported by the decoder, so "frame N"
        // means "the Nth frame we sample while reading forward".
        long long total = 0;
#ifdef SS_HAVE_VIDEO
        if (src.builtin_decode) {
            video::VideoReader r;
            if (r.open(src.input)) {
                total = r.info().frame_count;
                // Free here, and the only thing a later fallback would be
                // missing: ffmpeg seeks by time, not by frame.
                if (r.info().fps > 0.0)
                    src.video_seconds = (double)total / r.info().fps;
            }
        }
#endif
        if (total <= 0) {
            // No in-process decoding on this machine, or it could not open the
            // file: ffmpeg is what will read this video, in the preview and in
            // the run, so ffmpeg is who to ask how long it is.
            VideoFacts facts;
            if (ffmpeg_probe_video(src.ffmpeg_exe, src.input, facts, cancel)) {
                src.video_seconds = facts.duration;
                total = facts.frames;
            }
        }
        total = std::max(total, (long long)offers);
        for (int i = 0; i < offers; i++) {
            PreviewFrame f;
            f.index = std::min(total - 1, (long long)i * (total / offers));
            f.position = (float)((double)f.index / (double)std::max(1LL, total - 1));
            frames.push_back(f);
        }
        return;
    }

    // follow_directory_symlink for the same reason DatasetPrep walks that way:
    // a prepared capture's images/ is often a link into the raw one, and the
    // default iterator quietly returns nothing for it.
    for (fs::recursive_directory_iterator it(
             src.input, fs::directory_options::skip_permission_denied |
                            fs::directory_options::follow_directory_symlink, ec), end;
         !ec && it != end; it.increment(ec))
        if (it->is_regular_file(ec) && is_image_file(it->path()))
            all_files.push_back(it->path().string());
    std::sort(all_files.begin(), all_files.end());

    // A spread through the capture is enough to judge a setting, and keeps the
    // slider meaningful on a 3000-photo folder. The index kept is the one in
    // the FULL list, because that is what a run counts.
    const size_t n = all_files.size();
    const size_t want = std::min<size_t>(n, (size_t)offers);
    for (size_t i = 0; i < want; i++) {
        PreviewFrame f;
        f.index = want > 1 ? (long long)(i * (n - 1) / (want - 1)) : 0;
        f.path = all_files[(size_t)f.index];
        f.position = n > 1 ? (float)((double)f.index / (double)(n - 1)) : 0.0f;
        frames.push_back(f);
    }
}

bool load_preview_frame(const PreviewSource& src, const PreviewFrame& frame,
                        int& w, int& h, std::vector<uint8_t>& rgb,
                        std::string& error, const std::atomic<bool>& cancel) {
    w = h = 0;
    rgb.clear();
    if (!src.is_video) {
        if (frame.path.empty() || !app::load_rgb(frame.path, w, h, rgb)) {
            error = dmsg::preview_frame_unreadable.get();
            return false;
        }
        return true;
    }

#ifdef SS_HAVE_VIDEO
    if (src.builtin_decode) {
        // No seek: read forward to the frame this offer names.
        video::VideoReader r;
        if (r.open(src.input)) {
            for (long long i = 0; i <= frame.index; i++) {
                if (cancel.load()) return false;
                nn::Image f = r.readFrame();
                if (f.empty()) break;
                w = f.width;
                h = f.height;
                rgb = std::move(f.data);
            }
        }
    }
#endif
    // The same fallback the dataset run takes when the driver cannot decode
    // (DatasetPrep::extract_video). A whole subprocess for one still, which is
    // why it is not the first choice.
    if (rgb.empty()) {
        if (!command_exists(src.ffmpeg_exe)) {
            error = spirula::i18n::format(dmsg::preview_needs_ffmpeg,
                                          {src.ffmpeg_exe});
            return false;
        }
        // ffmpeg seeks by time, and the very end of a file is past the last
        // frame often enough to be worth backing away from.
        const double at = std::min((double)frame.position * src.video_seconds,
                                   std::max(0.0, src.video_seconds - 0.1));
        const fs::path tmp = temp_still("preview");
        const bool ok =
            ffmpeg_extract_frame(src.ffmpeg_exe, src.input, at, tmp.string(), cancel);
        if (cancel.load()) return false;
        if (ok) app::load_rgb(tmp.string(), w, h, rgb);
        std::error_code rm;
        fs::remove(tmp, rm);
    }
    if (rgb.empty()) {
        error = dmsg::preview_frame_unreadable.get();
        return false;
    }
    return true;
}

}  // namespace gui
