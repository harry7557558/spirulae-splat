// DatasetPrep.cpp -- see DatasetPrep.h.

#include "app/gui/DatasetPrep.h"

#include "app/gui/FrameSelect.h"
#include "app/gui/Subprocess.h"

#include "app_generated/mask_py.h"   // kMaskPy[] (CMake-embedded scripts/mask.py)

#include "external/stb_image.h"      // stbi_info (image size probe)

#ifdef SSPLAT_BUILD_SAM
#include "nn/Device.h"
#include "nn/io/Image.h"
#include "sam/Masking.h"
#endif
#ifdef SSPLAT_HAVE_VIDEO
#include "app/FrameExtract.h"
#include "video/Video.h"
#endif

#ifndef _WIN32
#include <ftw.h>
#endif

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>

namespace fs = std::filesystem;

namespace gui {

const char* const kVideoExtensions[kNumVideoExtensions] = {
    ".mp4", ".mov", ".mkv", ".webm", ".m4v", ".insv", ".avi",
    ".mts", ".m2ts", ".360", ".ts", ".wmv",
};

namespace {

// NOT std::filesystem::remove_all -- on the torch build libtorch.so interposes
// an ABI-incompatible copy (see AGENTS.md gotchas).
void remove_tree(const fs::path& p) {
#ifndef _WIN32
    nftw(p.string().c_str(),
         [](const char* f, const struct stat*, int, struct FTW*) {
             return ::remove(f);
         }, 16, FTW_DEPTH | FTW_PHYS);
#else
    std::error_code ec;
    std::filesystem::remove_all(p, ec);
#endif
}

bool is_image_file(const fs::path& p) {
    std::string e = p.extension().string();
    for (auto& c : e) c = (char)std::tolower((unsigned char)c);
    return e == ".jpg" || e == ".jpeg" || e == ".png" || e == ".webp" ||
           e == ".tif" || e == ".tiff" || e == ".bmp";
}

std::string lower_ext(const std::string& path) {
    std::string e = fs::path(path).extension().string();
    for (auto& c : e) c = (char)std::tolower((unsigned char)c);
    return e;
}

std::string human_duration(double seconds) {
    if (seconds < 1.0) return "a moment";
    char b[64];
    if (seconds < 90.0) std::snprintf(b, sizeof b, "%.0f s", seconds);
    else if (seconds < 5400.0) std::snprintf(b, sizeof b, "%.0f min", seconds / 60.0);
    else std::snprintf(b, sizeof b, "%.1f h", seconds / 3600.0);
    return b;
}

// Progress that answers "how long is this going to take", which is the only
// question a user has during a twenty-minute masking pass, and the one a
// counter that ticks every tenth frame does not answer.
//
// Rate-limited by wall clock rather than by a frame count: the same call site
// serves a decode running at a thousand frames a second and a segmentation
// running at one every two seconds.
class RateLimitedProgress {
public:
    using Clock = std::chrono::steady_clock;

    RateLimitedProgress(std::function<void(const std::string&)> log,
                        std::string noun, int64_t total)
        : _log(std::move(log)), _noun(std::move(noun)), _total(total),
          _start(Clock::now()), _last(_start) {}

    void update(int64_t done, bool force = false) {
        const auto now = Clock::now();
        const double since =
            std::chrono::duration<double>(now - _last).count();
        if (!force && (done == _reported || since < 2.0)) return;
        _reported = done;
        _last = now;
        // The rate is measured from the first item, not from the start: the
        // several seconds a checkpoint takes to reach the GPU would otherwise
        // be spread over every frame and put the first estimate out by 3x.
        if (_anchor_done < 0) {
            _anchor_done = done;
            _start = now;
        }
        const double elapsed = std::chrono::duration<double>(now - _start).count();
        const int64_t measured = done - _anchor_done;
        std::string line = "  " + std::to_string(done);
        if (_total > 0) line += " / " + std::to_string(_total);
        line += " " + _noun;
        if (measured > 0 && elapsed > 0.5) {
            const double per = elapsed / (double)measured;
            char rate[64];
            if (per >= 0.5) std::snprintf(rate, sizeof rate, "%.1f s each", per);
            else            std::snprintf(rate, sizeof rate, "%.0f/s", 1.0 / per);
            line += std::string("  (") + rate;
            if (_total > done)
                line += ", about " + human_duration(per * (double)(_total - done)) +
                        " left";
            line += ")";
        }
        _log(line);
    }

private:
    std::function<void(const std::string&)> _log;
    std::string _noun;
    int64_t     _total = 0;
    int64_t     _reported = -1;
    int64_t     _anchor_done = -1;   // count at the first report; see update()
    Clock::time_point _start, _last;
};

}  // namespace

const Backends& backends() {
    static const Backends probed = [] {
        Backends b;
#ifdef SSPLAT_HAVE_VIDEO
        b.video_reason = app::video_decode_availability();
        b.builtin_video = b.video_reason.empty();
        if (!b.builtin_video)
            b.video_note = "this graphics driver cannot decode video, so "
                           "frames are extracted with ffmpeg";
#else
        b.video_reason = "built without the video decoder "
                         "(-DSSPLAT_ENABLE_PATENTED=OFF)";
        b.video_note = "frames are extracted with ffmpeg";
#endif
#ifdef SSPLAT_BUILD_SAM
        b.builtin_masking = true;
#else
        b.masking_reason = "built without the segmentation module "
                           "(-DSSPLAT_BUILD_SAM=OFF)";
        b.masking_note =
            "Masks are made by an external Python script (scripts/mask.py with "
            "lang-segment-anything, which needs a CUDA PyTorch). Set the Python "
            "path under Tool locations if it is not on PATH.";
#endif
        return b;
    }();
    return probed;
}

int DatasetPrep::count_images(const std::string& dir) {
    int n = 0;
    std::error_code ec;
    for (fs::recursive_directory_iterator it(
             dir, fs::directory_options::skip_permission_denied, ec), end;
         !ec && it != end; it.increment(ec))
        if (it->is_regular_file(ec) && is_image_file(it->path())) n++;
    return n;
}

bool DatasetPrep::first_image_dims(const std::string& dir, int& W, int& H) {
    std::error_code ec;
    for (fs::recursive_directory_iterator it(
             dir, fs::directory_options::skip_permission_denied, ec), end;
         !ec && it != end; it.increment(ec))
        if (it->is_regular_file(ec) && is_image_file(it->path())) {
            int c = 0;
            return stbi_info(it->path().string().c_str(), &W, &H, &c) != 0;
        }
    return false;
}

int DatasetPrep::exec(const std::vector<std::string>& argv) {
    std::string cmd;
    for (const auto& a : argv) cmd += (cmd.empty() ? "$ " : " ") + a;
    _log(cmd);
    return run_process(argv, "", _log, _cancel);
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

bool DatasetPrep::run(const PrepJob& job, PrepResult& out, std::string& error) {
#ifdef SSPLAT_BUILD_SAM
    // Hand the GPU back on the way out, by whichever of the dozen exits is
    // taken. A SAM 3 checkpoint is about 2 GB of VRAM and the inference layer's
    // pool is process-wide and grow-only, so without this it stays resident
    // for the life of the GUI -- through the reconstruction and the training
    // run that follow, which are exactly what wants the memory back.
    //
    // Safe because the mask preview owns the only other Session, and the
    // dataset screen closes it before starting a job.
    struct ReleaseDevice {
        ~ReleaseDevice() { nn::shutdown(); }
    } release_device;
#endif
    const fs::path ws = job.workspace;
    std::error_code ec;
    fs::create_directories(ws, ec);

    if (job.is_video) {
        if (!extract_video(job, out, error)) return false;
    } else {
        out.image_dir = fs::absolute(job.input_path).string();
        // Photos are referenced where they are, not copied: a 40 GB folder of
        // raw captures does not want a second copy, and the parsers accept an
        // absolute image_dir.
        out.image_dir_cfg = out.image_dir;
    }

    out.n_images = count_images(out.image_dir);
    _log("Found " + std::to_string(out.n_images) + " images in " + out.image_dir);
    if (out.n_images < 3) {
        error = "need at least 3 images (found " + std::to_string(out.n_images) +
                ")";
        return false;
    }

    if (job.mask_enable) {
        if (job.mask_prompt.empty()) {
            error = "masking is on but the prompt is empty (e.g. \"people; cars\")";
            return false;
        }
        // A built-in video run already produced masks in the same pass.
        if (out.mask_dir.empty()) {
            const std::string rel = job.is_video ? "images" : out.image_dir_cfg;
            if (!generate_masks(job, out.image_dir, rel, error)) return false;
            out.mask_dir = (ws / "masks").string();
        }
    }
    return true;
}

// ---------------------------------------------------------------------------
// Video -> frames
// ---------------------------------------------------------------------------

bool DatasetPrep::extract_video(const PrepJob& job, PrepResult& out,
                                std::string& error) {
    const fs::path ws = job.workspace;
    out.image_dir = (ws / "images").string();
    out.image_dir_cfg = "images";

    // Resume: frames are moved into place in one batch after selection, so a
    // non-empty images/ means a previous extraction finished.
    if (job.resume) {
        const int have = count_images(out.image_dir);
        if (have > 0) {
            _log("Resume: keeping " + std::to_string(have) +
                 " extracted frames in " + out.image_dir +
                 " (delete the folder to re-extract)");
            std::error_code ec;
            out.multi_track = fs::is_directory(fs::path(out.image_dir) / "cam1", ec);
            if (job.mask_enable) {
                const fs::path masks = ws / "masks";
                std::error_code mec;
                if (fs::is_directory(masks, mec) && !fs::is_empty(masks, mec)) {
                    _log("Resume: keeping the masks in " + masks.string());
                    out.mask_dir = masks.string();
                }
            }
            return true;
        }
    }

    const bool want_builtin = !job.force_external_decode && backends().builtin_video;
    if (want_builtin) {
        if (extract_video_builtin(job, out, error)) return true;
        if (_cancel.load()) return false;
        // A container or profile the driver cannot decode is exactly what the
        // fallback is for, and the user should not have to know which is which.
        _log("Built-in decoding could not handle this file (" + error +
             "); falling back to ffmpeg");
    }
    return extract_video_ffmpeg(job, out, error);
}

bool DatasetPrep::extract_video_builtin(const PrepJob& job, PrepResult& out,
                                        std::string& error) {
#ifndef SSPLAT_HAVE_VIDEO
    (void)job; (void)out;
    error = backends().video_reason;
    return false;
#else
    _stage("Extracting frames (GPU decode)");

    // fps -> "one frame every N source frames". The source rate is what the
    // container states; a variable-rate file is close enough for this.
    std::string probe_err;
    const int tracks = app::video_track_count(job.input_path, probe_err);
    if (tracks <= 0) {
        error = probe_err.empty() ? "no video track" : probe_err;
        return false;
    }
    out.multi_track = tracks > 1;

    double src_fps = 30.0;
    int64_t src_frames = 0;
    {
        video::VideoReader r;
        if (r.open(job.input_path)) {
            if (r.info().fps > 1.0) src_fps = r.info().fps;
            src_frames = r.info().frame_count;
        }
    }
    const int window = std::max(job.sharp_window, 1);
    int skip = (int)std::lround(src_fps / std::max(job.video_fps, 0.01f));
    skip = std::max(skip, 1);

    // What the container claims, for the estimate; the loop stops on the real
    // end of stream either way.
    int64_t expect = src_frames > 0 ? (src_frames / skip) * (int64_t)tracks : 0;
    if (job.max_frames > 0 && (expect == 0 || expect > (int64_t)job.max_frames * tracks))
        expect = (int64_t)job.max_frames * tracks;

    app::FrameExtractJob fx;
    fx.input = job.input_path;
    fx.image_dir = out.image_dir;
    fx.mask_dir = (fs::path(job.workspace) / "masks").string();
    fx.skip = skip;
    fx.keep = window > 1 ? window : 0;
    fx.max_frames = job.max_frames;
    fx.quality = 95;
    if (job.mask_enable && !job.force_external_masking &&
        !job.mask_model_path.empty()) {
        // Masking rides along on the decode: the frame is already on the
        // device, so this is far cheaper than a second pass over the JPEGs.
        _stage("Extracting frames and masking (GPU)");
        fx.mask.model = job.mask_model_path;
        fx.mask.text = job.mask_prompt;
        fx.mask.neg_text = job.mask_negative_prompt;
        fx.mask.keep_prompted = job.mask_keep_subject;
        fx.mask.max_size = job.mask_max_image_size;
        fx.mask.video = true;
    }

    app::FrameExtractSinks sinks;
    sinks.log = _log;
    sinks.cancel = &_cancel;
    RateLimitedProgress progress(
        _log, fx.mask.model.empty() ? "frames written" : "frames written and masked",
        expect);
    sinks.progress = [&](int64_t written, int64_t decoded) {
        (void)decoded;
        progress.update(written);
    };

    app::FrameExtractStats stats;
    if (!app::extract_frames(fx, sinks, stats, error)) return false;
    _log(app::format_extract_stats(stats, out.image_dir, !fx.mask.model.empty()));
    if (stats.written == 0) {
        error = "no frames were extracted";
        return false;
    }
    if (!fx.mask.model.empty()) out.mask_dir = fx.mask_dir;
    return true;
#endif
}

bool DatasetPrep::extract_video_ffmpeg(const PrepJob& job, PrepResult& out,
                                       std::string& error) {
    const fs::path ws = job.workspace;
    if (!command_exists(job.ffmpeg_exe)) {
        error = "ffmpeg not found ('" + job.ffmpeg_exe +
                "'). Install it, set its path under Tool locations, or build "
                "with -DSSPLAT_ENABLE_PATENTED=ON for in-process decoding.";
        return false;
    }

    // Multi-track videos (Insta360 .insv): one folder per track, one camera
    // per folder, as in scripts/extract_frames.py.
    std::vector<int> streams = {0};
    if (lower_ext(job.input_path) == ".insv") {
        std::vector<int> found;
        run_process({"ffprobe", "-v", "error", "-select_streams", "v",
                     "-show_entries", "stream=index", "-of", "csv=p=0",
                     job.input_path}, "",
                    [&](const std::string& l) {
                        try { found.push_back(std::stoi(l)); } catch (...) {}
                    }, _cancel);
        if (found.size() > 1) streams = found;
    }
    out.multi_track = streams.size() > 1;

    const int window = std::max(job.sharp_window, 1);
    for (size_t tr = 0; tr < streams.size(); tr++) {
        std::string track_path = job.input_path;
        const fs::path out_dir = streams.size() > 1
            ? ws / "images" / ("cam" + std::to_string(tr))
            : ws / "images";
        if (job.resume && count_images(out_dir.string()) > 0) {
            _log("Resume: keeping the frames already in " + out_dir.string());
            continue;
        }
        if (streams.size() > 1) {
            _stage("Splitting video track " + std::to_string(tr));
            const fs::path tmp_track =
                ws / ("track_cam" + std::to_string(tr) + ".mp4");
            int rc = exec({job.ffmpeg_exe, "-nostdin", "-y", "-i", job.input_path,
                           "-map", "0:v:" + std::to_string(tr), "-c", "copy",
                           tmp_track.string()});
            if (rc == kCancelled) { error = "cancelled"; return false; }
            if (rc != 0) { error = "ffmpeg track split failed (see log)"; return false; }
            track_path = tmp_track.string();
        }

        _stage(window > 1 ? "Extracting candidate frames (ffmpeg)"
                          : "Extracting frames (ffmpeg)");
        const fs::path cand = ws / "frames_tmp";
        remove_tree(cand);
        std::error_code ec;
        fs::create_directories(cand, ec);
        char vf[64];
        std::snprintf(vf, sizeof vf, "fps=%g", (double)job.video_fps * window);
        int rc = exec({job.ffmpeg_exe, "-nostdin", "-y", "-i", track_path,
                       "-vf", vf, "-qscale:v", "2",
                       (cand / "c_%06d.jpg").string()});
        if (rc == kCancelled) { error = "cancelled"; return false; }
        if (rc != 0) { error = "ffmpeg frame extraction failed (see log)"; return false; }

        if (window > 1) _stage("Selecting sharpest frames (multithreaded)");
        fs::create_directories(out_dir, ec);
        const int kept = select_sharpest_frames(cand.string(), out_dir.string(), "",
                                                window, job.max_frames, _log, _cancel);
        remove_tree(cand);
        if (streams.size() > 1) fs::remove(track_path, ec);
        if (kept < 0) {
            error = _cancel.load() ? "cancelled" : "frame selection failed";
            return false;
        }
        _log("Kept " + std::to_string(kept) + " frames -> " + out_dir.string());
    }
    return true;
}

// ---------------------------------------------------------------------------
// Masks
// ---------------------------------------------------------------------------

bool DatasetPrep::generate_masks(const PrepJob& job, const std::string& images,
                                 const std::string& images_rel,
                                 std::string& error) {
    const bool want_builtin = !job.force_external_masking &&
                              backends().builtin_masking &&
                              !job.mask_model_path.empty();
    if (want_builtin) return generate_masks_builtin(job, images, error);
    return generate_masks_python(job, images_rel, error);
}

bool DatasetPrep::generate_masks_builtin(const PrepJob& job,
                                         const std::string& images,
                                         std::string& error) {
#ifndef SSPLAT_BUILD_SAM
    (void)job; (void)images;
    error = backends().masking_reason;
    return false;
#else
    _stage("Generating masks (segmentation)");

    std::vector<fs::path> files;
    std::error_code ec;
    for (fs::recursive_directory_iterator it(
             images, fs::directory_options::skip_permission_denied, ec), end;
         !ec && it != end; it.increment(ec))
        if (it->is_regular_file(ec) && is_image_file(it->path()))
            files.push_back(it->path());
    std::sort(files.begin(), files.end());
    if (files.empty()) { error = "no images to mask"; return false; }

    sam::MaskOptions mo;
    mo.model = job.mask_model_path;
    mo.text = job.mask_prompt;
    mo.neg_text = job.mask_negative_prompt;
    mo.keep_prompted = job.mask_keep_subject;
    mo.max_size = job.mask_max_image_size;
    // A folder of photos is not a video: consecutive files may be anywhere in
    // the scene, so tracking across them would carry a memory bank that does
    // not apply. Frames extracted from a video are the ordered case, and they
    // are masked during extraction instead.
    mo.video = false;

    sam::Masker masker;
    if (!masker.init(mo, error)) return false;

    const fs::path mask_root = fs::path(job.workspace) / "masks";
    fs::create_directories(mask_root, ec);
    const fs::path image_root(images);

    int done = 0;
    RateLimitedProgress progress(_log, "images masked", (int64_t)files.size());
    for (const fs::path& f : files) {
        if (_cancel.load()) { error = "cancelled"; return false; }
        // Masks mirror the image tree, so cam0/ and cam1/ keep their names.
        fs::path rel = fs::relative(f, image_root, ec);
        if (ec || rel.empty()) rel = f.filename();
        const fs::path dst = mask_root / rel.parent_path() /
                             (rel.stem().string() + ".png");
        if (job.resume && fs::exists(dst, ec)) { ++done; continue; }
        fs::create_directories(dst.parent_path(), ec);

        nn::Image img = nn::load_image(f.string());
        if (img.empty()) {
            _log("warning: could not read " + f.string() + "; skipped");
            continue;
        }
        sam::Mask mask;
        if (!masker.run(img, mask, nullptr)) {
            error = "masking failed on " + f.filename().string() + ": " +
                    masker.lastError();
            return false;
        }
        sam::save_mask_png(mask, dst.string());
        progress.update(++done, done == (int)files.size());
    }
    return true;
#endif
}

// The Python fallback: the embedded scripts/mask.py run through an external
// interpreter with lang-segment-anything. It prints an install hint and exits
// 0 when the packages are missing, so that is detected from its output rather
// than its exit code.
bool DatasetPrep::generate_masks_python(const PrepJob& job,
                                        const std::string& images_rel,
                                        std::string& error) {
    _stage("Generating masks (external Python)");
    if (!command_exists(job.python_exe)) {
        error = "Python not found ('" + job.python_exe +
                "'); external masking needs Python with the "
                "lang-segment-anything package. Set the Python path under "
                "Tool locations, or use the built-in segmentation.";
        return false;
    }
    const fs::path ws = job.workspace;
    const fs::path script = ws / ".ssplat_mask.py";
    {
        FILE* f = std::fopen(script.string().c_str(), "wb");
        if (!f) { error = "cannot write " + script.string(); return false; }
        std::fwrite(kMaskPy, 1, kMaskPySize, f);
        std::fclose(f);
    }
    std::vector<std::string> argv = {
        job.python_exe, script.string(), ws.string(),
        "--prompt", job.mask_prompt,
        "--images", images_rel,
        "--masks", "masks",
        "--max_image_size", std::to_string(job.mask_max_image_size),
        "--model", job.mask_model_name,
    };
    if (!job.mask_negative_prompt.empty()) {
        argv.push_back("--negative_prompt");
        argv.push_back(job.mask_negative_prompt);
    }
    // Without this the external path silently ignores the polarity and always
    // removes what the prompt named -- the exact opposite of what an object
    // capture asked for.
    if (job.mask_keep_subject) argv.push_back("--keep_prompted");
    std::string install_hint;
    std::string cmd;
    for (const auto& a : argv) cmd += (cmd.empty() ? "$ " : " ") + a;
    _log(cmd);
    const int rc = run_process(argv, "", [&](const std::string& l) {
        _log(l);
        if (l.find("not found or not installed properly") != std::string::npos ||
            l.find("ModuleNotFoundError") != std::string::npos)
            install_hint = l;
    }, _cancel);
    if (rc == kCancelled) { error = "cancelled"; return false; }

    std::error_code ec;
    const bool have_masks = fs::is_directory(ws / "masks", ec) &&
                            !fs::is_empty(ws / "masks", ec);
    if (!install_hint.empty() || rc != 0 || !have_masks) {
        error = "Mask generation failed";
        if (!install_hint.empty())
            error += " -- missing Python packages. Install "
                     "lang-segment-anything (pip install git+https://github.com/"
                     "luca-medeiros/lang-segment-anything, needs CUDA PyTorch)"
                     + std::string(job.mask_model_name == "sam3"
                         ? ", or for SAM 3: https://github.com/facebookresearch/sam3"
                         : "");
        else
            error += " (see log)";
        return false;
    }
    return true;
}

}  // namespace gui
