// DatasetPrep.cpp -- see DatasetPrep.h.

#include "app/gui/DatasetPrep.h"

#include "app/gui/FrameSelect.h"
#include "app/gui/Subprocess.h"

#include "app_generated/mask_py.h"   // kMaskPy[] (CMake-embedded scripts/mask.py)

#include "external/stb_image.h"      // stbi_info (image size probe)

#ifdef SS_BUILD_SAM
#include "app/WriterPool.h"
#include "nn/Device.h"
#include "nn/io/Image.h"
#include "sam/Masking.h"
#endif
#ifdef SS_HAVE_VIDEO
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

// Every walk over an image tree follows directory symlinks. A prepared capture
// whose images/ and masks/ are links into the raw one is an ordinary layout,
// and the default iterator returns nothing at all for it -- which reads as "no
// images here" rather than as "not looked".
constexpr fs::directory_options kWalk =
    fs::directory_options::skip_permission_denied |
    fs::directory_options::follow_directory_symlink;

// Images under `dir`, sorted, never descending into `skip` (see count_images).
std::vector<fs::path> walk_images(const fs::path& dir, const fs::path& skip = {}) {
    std::vector<fs::path> out;
    std::error_code ec;
    const bool guard = !skip.empty() && fs::is_directory(skip, ec);
    for (fs::recursive_directory_iterator it(dir, kWalk, ec), end; !ec && it != end;
         it.increment(ec)) {
        if (guard && it->is_directory(ec) &&
            fs::equivalent(it->path(), skip, ec)) {
            it.disable_recursion_pending();
            continue;
        }
        if (it->is_regular_file(ec) && is_image_file(it->path()))
            out.push_back(it->path());
    }
    std::sort(out.begin(), out.end());
    return out;
}

#ifdef SS_BUILD_SAM
// Preview clicks -> the seeds the masker takes. Clicks of one object made on
// one frame become ONE prompt (several positive points describe one thing);
// clicks of the same object on another frame become a second prompt, which the
// masker applies as a correction when it gets there.
//
// `total_frames` is how many frames the run will see, and is only needed when
// the click's own frame numbering does not survive to it -- see MaskClick.
// Pass 0 to keep the recorded index.
std::vector<sam::SeedPrompt> seeds_from_clicks(const std::vector<MaskClick>& clicks,
                                               int64_t total_frames) {
    std::vector<sam::SeedPrompt> seeds;
    for (const MaskClick& c : clicks) {
        const int64_t frame =
            total_frames > 1
                ? (int64_t)std::llround((double)c.position * (double)(total_frames - 1))
                : (total_frames == 1 ? 0 : c.frame);
        sam::SeedPrompt* seed = nullptr;
        for (sam::SeedPrompt& s : seeds)
            if (s.object == c.object && s.frame == frame) seed = &s;
        if (!seed) {
            seeds.push_back({});
            seed = &seeds.back();
            seed->object = c.object;
            seed->frame = frame;
        }
        if (c.positive) seed->prompt.pos_points.push_back({c.x, c.y});
        else            seed->prompt.neg_points.push_back({c.x, c.y});
    }
    return seeds;
}
#endif

std::string lower_ext(const std::string& path) {
    std::string e = fs::path(path).extension().string();
    for (auto& c : e) c = (char)std::tolower((unsigned char)c);
    return e;
}

// base / sub, where an empty `sub` means base itself (and not "base/").
fs::path under(const std::string& base, const std::string& sub) {
    return sub.empty() ? fs::path(base) : fs::path(base) / sub;
}

// Are the job's clicked objects prompts for THIS input? See
// PrepJob::mask_clicks_source: an unnamed source means the clicks were drawn
// before there was anything to distinguish, which is the single-input case.
bool clicks_apply_to(const PrepJob& job, const PrepInput& in) {
    return job.mask_clicks.empty() || job.mask_clicks_source.empty() ||
           job.mask_clicks_source == in.path;
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

bool is_video_path(const std::string& path) {
    const std::string e = lower_ext(path);
    for (const char* v : kVideoExtensions)
        if (e == v) return true;
    return false;
}

bool is_dual_fisheye_path(const std::string& path) {
    return lower_ext(path) == ".insv";
}

namespace {

// The last component, without a trailing separator ("a/b/" -> "b").
std::string leaf_name(const fs::path& p) {
    return p.filename().empty() ? p.parent_path().filename().string()
                                : p.filename().string();
}

bool named(const fs::path& p, const char* what) {
    std::string n = leaf_name(p);
    for (auto& c : n) c = (char)std::tolower((unsigned char)c);
    return n == what;
}

// A folder that exists and has at least one image in it, at any depth. Stops at
// the first one: this runs on the UI thread, and the tree can hold thousands.
bool any_image(const fs::path& p) {
    std::error_code ec;
    if (!fs::is_directory(p, ec)) return false;
    for (fs::recursive_directory_iterator it(p, kWalk, ec), end; !ec && it != end;
         it.increment(ec))
        if (it->is_regular_file(ec) && is_image_file(it->path())) return true;
    return false;
}

}  // namespace

bool folder_has_images(const std::string& dir) { return any_image(dir); }

WorkspaceState probe_workspace(const std::string& workspace,
                               const std::vector<PrepInput>& inputs) {
    WorkspaceState st;
    std::error_code ec;
    const fs::path ws(workspace);
    if (workspace.empty() || !fs::is_directory(ws, ec)) return st;

    // A folder is the run's leftover only if the run would write it. When
    // images/ or masks/ under the output IS an input, it is the capture.
    auto is_input = [&](const fs::path& dir, bool masks) {
        for (const PrepInput& in : inputs) {
            const std::string& p = masks ? in.mask_dir : in.path;
            if (!p.empty() && fs::equivalent(dir, p, ec)) return true;
        }
        return false;
    };
    auto has_content = [&](const fs::path& p) {
        return fs::is_directory(p, ec) && !fs::is_empty(p, ec);
    };

    st.frames = has_content(ws / "images") && !is_input(ws / "images", false);
    st.masks = has_content(ws / "masks") && !is_input(ws / "masks", true);
    st.features = has_content(ws / "features") || fs::exists(ws / "matches.bin", ec) ||
                  fs::exists(ws / "database.db", ec);
    st.model = has_content(ws / "sparse");
    return st;
}

bool is_mask_folder(const std::string& path) {
    return named(fs::path(path), "masks");
}

void resolve_photo_folder(const std::string& picked, std::string& images,
                          std::string& masks) {
    std::error_code ec;
    const fs::path p = fs::absolute(picked, ec);
    images = p.string();
    masks.clear();
    // A dataset folder: index images/, not the folder holding it (which also
    // holds the masks, the point cloud, and whatever else was left there).
    // This is the same probe `spirula sfm auto` prints as "<dir> contains
    // images/, using <dir>/images as the image directory".
    if (any_image(p / "images")) images = (p / "images").string();
    // The masks belong beside the images: under the folder that holds them, or
    // -- when the folder IS the images/ of a dataset -- next to it.
    std::vector<fs::path> candidates{fs::path(images) / "masks"};
    if (named(images, "images"))
        candidates.push_back(fs::path(images).parent_path() / "masks");
    for (const fs::path& cand : candidates) {
        if (fs::is_directory(cand, ec) && any_image(cand)) {
            masks = cand.string();
            break;
        }
    }
}

const Backends& backends() {
    static const Backends probed = [] {
        Backends b;
#ifdef SS_HAVE_VIDEO
        b.video_reason = app::video_decode_availability();
        b.builtin_video = b.video_reason.empty();
        if (!b.builtin_video)
            b.video_note = "this graphics driver cannot decode video, so "
                           "frames are extracted with ffmpeg";
#else
        b.video_reason = "built without the video decoder "
                         "(-DSS_ENABLE_PATENTED=OFF)";
        b.video_note = "frames are extracted with ffmpeg";
#endif
#ifdef SS_BUILD_SAM
        b.builtin_masking = true;
#else
        b.masking_reason = "built without the segmentation module "
                           "(-DSS_BUILD_SAM=OFF)";
        b.masking_note =
            "Masks are made by an external Python script (scripts/mask.py with "
            "lang-segment-anything, which needs a CUDA PyTorch). Set the Python "
            "path under Tool locations if it is not on PATH.";
#endif
        return b;
    }();
    return probed;
}

int DatasetPrep::count_images(const std::string& dir, const std::string& skip) {
    return (int)walk_images(dir, skip).size();
}

bool DatasetPrep::first_image_dims(const std::string& dir, int& W, int& H) {
    for (const fs::path& f : walk_images(dir)) {
        int c = 0;
        if (stbi_info(f.string().c_str(), &W, &H, &c) != 0) return true;
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
#ifdef SS_BUILD_SAM
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
    if (job.inputs.empty()) {
        error = "nothing to prepare: no video or photo folder was picked";
        return false;
    }

    // Where each input's images and masks ended up, so the masking pass below
    // can run per input (see generate_masks) instead of over one flat tree.
    struct Prepared {
        std::string images, masks;      // absolute
        std::string images_rel, masks_rel;  // relative to the workspace
        // This input's masks already exist: written on the way out of the
        // decoder, brought along by the input, or kept by a resumed run.
        bool have_masks = false;
    };
    std::vector<Prepared> per(job.inputs.size());
    // A masks/ nested under the images is not a folder of views (see
    // count_images); only the in-place case can have one, since everything else
    // writes into the dataset's own images/.
    std::string skip_dir;
    // Is the dataset's own masks/ this run's, rather than an earlier one's? A
    // stale masks/ from a previous run with masking on must not be picked up by
    // a run that turned it off -- that is what --no-masks exists to say.
    bool want_masks = job.mask_enable;

    if (reads_photos_in_place(job.inputs)) {
        // Photos are referenced where they are, not copied: a 40 GB folder of
        // raw captures does not want a second copy, and the parsers accept an
        // absolute image_dir.
        const PrepInput& in = job.inputs[0];
        out.image_dir = fs::absolute(in.path).string();
        out.image_dir_cfg = out.image_dir;
        per[0].images = per[0].images_rel = out.image_dir;
        if (in.mask_dir.empty()) {
            // Nothing came with them, so anything generated goes in the
            // dataset, next to the reconstruction rather than next to the
            // photos -- a folder we were only asked to read.
            per[0].masks = (ws / "masks").string();
            per[0].masks_rel = "masks";
        } else {
            per[0].masks = per[0].masks_rel = fs::absolute(in.mask_dir, ec).string();
            per[0].have_masks = true;
            out.mask_dir = out.mask_dir_cfg = per[0].masks;
            skip_dir = per[0].masks;
            _log("Using the masks that came with the photos: " + out.mask_dir);
        }
    } else {
        out.image_dir = (ws / "images").string();
        out.image_dir_cfg = "images";
        for (size_t i = 0; i < job.inputs.size(); i++) {
            const PrepInput& in = job.inputs[i];
            Prepared& p = per[i];
            p.images = under(out.image_dir, in.subdir).string();
            p.masks = under((ws / "masks").string(), in.subdir).string();
            p.images_rel = under("images", in.subdir).generic_string();
            p.masks_rel = under("masks", in.subdir).generic_string();
            // Several inputs share one image tree only by living in their own
            // folders, and a folder is what makes them separate cameras.
            if (!in.subdir.empty()) out.per_folder_cameras = true;
            if (in.is_video) {
                if (!extract_video(job, in, p.images, p.masks, out, p.have_masks,
                                   error))
                    return false;
            } else if (!gather_photos(job, in, p.images, p.masks, p.have_masks,
                                      error)) {
                return false;
            }
            // Masks an input brought with it are in the dataset now, so they
            // count even when nothing asked for masking.
            if (p.have_masks && !in.mask_dir.empty()) want_masks = true;
        }
    }

    out.n_images = count_images(out.image_dir, skip_dir);
    _log("Found " + std::to_string(out.n_images) + " images in " + out.image_dir);
    if (out.n_images < 3) {
        error = "need at least 3 images (found " + std::to_string(out.n_images) +
                ")";
        return false;
    }

    if (job.mask_enable) {
        if (job.mask_prompt.empty() && job.mask_clicks.empty()) {
            error = "masking is on but nothing says what to mask: type a prompt "
                    "(e.g. \"people; cars\"), or open \"Try the mask\" and click "
                    "the object";
            return false;
        }
        for (size_t i = 0; i < job.inputs.size(); i++) {
            // Already masked: by the decoder in the same pass, or by whoever
            // made the masks this input arrived with. Segmenting over those
            // would replace an answer the user already has.
            if (per[i].have_masks) continue;
            if (!generate_masks(job, job.inputs[i], per[i].images,
                                per[i].images_rel, per[i].masks, per[i].masks_rel,
                                error))
                return false;
        }
    }
    // Masks are whatever ended up in the dataset's own masks/ -- generated
    // here, written by the decoder, or linked in beside gathered photos. The
    // in-place case has already named the folder it reads them from.
    std::error_code mec;
    if (out.mask_dir.empty() && want_masks && fs::is_directory(ws / "masks", mec) &&
        !fs::is_empty(ws / "masks", mec)) {
        out.mask_dir = (ws / "masks").string();
        out.mask_dir_cfg = "masks";
    }
    return true;
}

// ---------------------------------------------------------------------------
// Video -> frames
// ---------------------------------------------------------------------------

bool DatasetPrep::extract_video(const PrepJob& job, const PrepInput& in,
                                const std::string& images,
                                const std::string& masks, PrepResult& out,
                                bool& masked, std::string& error) {
    // Resume: frames are moved into place in one batch after selection, so a
    // non-empty folder means a previous extraction of THIS input finished.
    if (job.resume) {
        const int have = count_images(images);
        if (have > 0) {
            _log("Resume: keeping " + std::to_string(have) +
                 " extracted frames in " + images +
                 " (delete the folder to re-extract)");
            std::error_code ec;
            if (fs::is_directory(fs::path(images) / "cam1", ec))
                out.per_folder_cameras = true;
            if (job.mask_enable) {
                std::error_code mec;
                if (fs::is_directory(masks, mec) && !fs::is_empty(masks, mec)) {
                    _log("Resume: keeping the masks in " + masks);
                    masked = true;
                }
            }
            return true;
        }
    }

    const bool want_builtin = !job.force_external_decode && backends().builtin_video;
    if (want_builtin) {
        if (extract_video_builtin(job, in, images, masks, out, masked, error))
            return true;
        if (_cancel.load()) return false;
        // A container or profile the driver cannot decode is exactly what the
        // fallback is for, and the user should not have to know which is which.
        _log("Built-in decoding could not handle this file (" + error +
             "); falling back to ffmpeg");
    }
    return extract_video_ffmpeg(job, in, images, out, error);
}

bool DatasetPrep::extract_video_builtin(const PrepJob& job, const PrepInput& in,
                                        const std::string& images,
                                        const std::string& masks, PrepResult& out,
                                        bool& masked, std::string& error) {
#ifndef SS_HAVE_VIDEO
    (void)job; (void)in; (void)images; (void)masks; (void)out; (void)masked;
    error = backends().video_reason;
    return false;
#else
    _stage("Extracting frames (GPU decode)");
    _log("Video: " + in.path);

    // fps -> "one frame every N source frames". The source rate is what the
    // container states; a variable-rate file is close enough for this.
    std::string probe_err;
    const int tracks = app::video_track_count(in.path, probe_err);
    if (tracks <= 0) {
        error = probe_err.empty() ? "no video track" : probe_err;
        return false;
    }
    if (tracks > 1) out.per_folder_cameras = true;

    double src_fps = 30.0;
    int64_t src_frames = 0;
    {
        video::VideoReader r;
        if (r.open(in.path)) {
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
    fx.input = in.path;
    fx.image_dir = images;
    fx.mask_dir = masks;
    fx.skip = skip;
    fx.keep = window > 1 ? window : 0;
    fx.max_frames = job.max_frames;
    fx.quality = 95;
    const bool clicks_here = clicks_apply_to(job, in);
    if (job.mask_enable && !job.force_external_masking &&
        !job.mask_model_path.empty() &&
        (!job.mask_prompt.empty() || clicks_here)) {
        // Masking rides along on the decode: the frame is already on the
        // device, so this is far cheaper than a second pass over the JPEGs.
        _stage("Extracting frames and masking (GPU)");
        fx.mask.model = job.mask_model_path;
        fx.mask.text = job.mask_prompt;
        fx.mask.neg_text = job.mask_negative_prompt;
        fx.mask.keep_prompted = job.mask_keep_subject;
        fx.mask.max_size = job.mask_max_image_size;
        fx.mask.video = true;
        // This path reads the video itself, so a click's decoded index means
        // the same thing here as it did in the preview.
        if (clicks_here) fx.mask.seeds = seeds_from_clicks(job.mask_clicks, 0);
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
    _log(app::format_extract_stats(stats, images, !fx.mask.model.empty()));
    if (stats.written == 0) {
        error = "no frames were extracted";
        return false;
    }
    if (!fx.mask.model.empty()) masked = true;
    return true;
#endif
}

bool DatasetPrep::extract_video_ffmpeg(const PrepJob& job, const PrepInput& in,
                                       const std::string& images,
                                       PrepResult& out, std::string& error) {
    const fs::path ws = job.workspace;
    if (!command_exists(job.ffmpeg_exe)) {
        error = "ffmpeg not found ('" + job.ffmpeg_exe +
                "'). Install it, set its path under Tool locations, or build "
                "with -DSS_ENABLE_PATENTED=ON for in-process decoding.";
        return false;
    }
    _log("Video: " + in.path);

    // Multi-track videos (Insta360 .insv): one folder per track, one camera
    // per folder, as in scripts/extract_frames.py.
    std::vector<int> streams = {0};
    if (is_dual_fisheye_path(in.path)) {
        std::vector<int> found;
        run_process({"ffprobe", "-v", "error", "-select_streams", "v",
                     "-show_entries", "stream=index", "-of", "csv=p=0",
                     in.path}, "",
                    [&](const std::string& l) {
                        try { found.push_back(std::stoi(l)); } catch (...) {}
                    }, _cancel);
        if (found.size() > 1) streams = found;
    }
    if (streams.size() > 1) out.per_folder_cameras = true;

    const int window = std::max(job.sharp_window, 1);
    for (size_t tr = 0; tr < streams.size(); tr++) {
        std::string track_path = in.path;
        const fs::path out_dir = streams.size() > 1
            ? fs::path(images) / ("cam" + std::to_string(tr))
            : fs::path(images);
        if (job.resume && count_images(out_dir.string()) > 0) {
            _log("Resume: keeping the frames already in " + out_dir.string());
            continue;
        }
        if (streams.size() > 1) {
            _stage("Splitting video track " + std::to_string(tr));
            const fs::path tmp_track =
                ws / ("track_cam" + std::to_string(tr) + ".mp4");
            int rc = exec({job.ffmpeg_exe, "-nostdin", "-y", "-i", in.path,
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
// Photos -> the dataset's own images/
// ---------------------------------------------------------------------------

// Only when the photos cannot be read where they are: a job with more than one
// input reconstructs from ONE image tree, and the folders under it are what
// make its inputs separate cameras.
//
// Hard links where the filesystem gives one, a copy otherwise. A link costs a
// directory entry, which matters when the alternative is a second copy of a
// folder of raw captures; falling back is what makes it work across devices
// (and on a filesystem that has no links at all).
bool DatasetPrep::gather_photos(const PrepJob& job, const PrepInput& in,
                                const std::string& images,
                                const std::string& masks, bool& have_masks,
                                std::string& error) {
    std::error_code ec;
    const fs::path src = fs::absolute(in.path, ec);
    if (!fs::is_directory(src, ec)) {
        error = "not a folder: " + in.path;
        return false;
    }
    _stage("Collecting photos");

    // The masks come across too, into the folder that mirrors the images -- a
    // mask is found by its image's relative name, so the two trees have to move
    // together or every mask stops matching.
    const bool with_masks = !in.mask_dir.empty();
    struct Tree { fs::path from, to; const char* noun; };
    std::vector<Tree> trees{{src, fs::path(images), "photos"}};
    if (with_masks)
        trees.push_back({fs::absolute(in.mask_dir, ec), fs::path(masks), "masks"});

    for (const Tree& t : trees) {
        _log(std::string(t.noun) + ": " + t.from.string() + " -> " + t.to.string());
        const std::vector<fs::path> files = walk_images(t.from);
        if (files.empty()) {
            error = std::string("no ") + t.noun + " in " + t.from.string();
            return false;
        }
        int linked = 0, copied = 0, kept = 0;
        RateLimitedProgress progress(_log, std::string(t.noun) + " collected",
                                     (int64_t)files.size());
        for (const fs::path& f : files) {
            if (_cancel.load()) { error = "cancelled"; return false; }
            fs::path rel = f.lexically_relative(t.from);
            if (rel.empty() || *rel.begin() == "..") rel = f.filename();
            const fs::path out_path = t.to / rel;
            if (fs::exists(out_path, ec)) { kept++; continue; }
            fs::create_directories(out_path.parent_path(), ec);
            std::error_code link_ec;
            fs::create_hard_link(f, out_path, link_ec);
            if (!link_ec) {
                linked++;
            } else {
                std::error_code copy_ec;
                fs::copy_file(f, out_path, copy_ec);
                if (copy_ec) {
                    error = "could not put " + f.string() + " into " +
                            t.to.string() + " (" + copy_ec.message() + ")";
                    return false;
                }
                copied++;
            }
            progress.update(linked + copied + kept);
        }
        _log("  " + std::to_string(linked) + " linked, " + std::to_string(copied) +
             " copied, " + std::to_string(kept) + " already there");
    }
    if (with_masks) have_masks = true;
    return true;
}

// ---------------------------------------------------------------------------
// Masks
// ---------------------------------------------------------------------------

bool DatasetPrep::generate_masks(const PrepJob& job, const PrepInput& in,
                                 const std::string& images,
                                 const std::string& images_rel,
                                 const std::string& masks,
                                 const std::string& masks_rel,
                                 std::string& error) {
    // A click describes one frame of one capture (see PrepJob::mask_clicks_
    // source), so an input it does not belong to has only the text prompt --
    // and if there is none, nothing at all to mask by. Saying so beats writing
    // masks from a prompt the user never gave for these images.
    if (!clicks_apply_to(job, in) && job.mask_prompt.empty()) {
        _log("Note: " + in.path + " is not the input the clicked objects were "
             "drawn on and there is no text prompt, so its frames are left "
             "unmasked. Add a prompt, or click the object on this input too.");
        return true;
    }
    const bool want_builtin = !job.force_external_masking &&
                              backends().builtin_masking &&
                              !job.mask_model_path.empty();
    if (want_builtin) return generate_masks_builtin(job, in, images, masks, error);
    if (!job.mask_clicks.empty()) {
        // The Python fallback is lang-segment-anything: text in, masks out. It
        // has no way to take a click, so saying so beats writing masks that
        // quietly ignore half of what the user asked for.
        error = "clicked objects need the built-in segmentation; the external "
                "Python masker only understands text prompts. Turn off "
                "\"external masking\", or describe the object in words.";
        return false;
    }
    return generate_masks_python(job, images_rel, masks_rel, error);
}

bool DatasetPrep::generate_masks_builtin(const PrepJob& job, const PrepInput& in,
                                         const std::string& images,
                                         const std::string& masks,
                                         std::string& error) {
#ifndef SS_BUILD_SAM
    (void)job; (void)in; (void)images; (void)masks;
    error = backends().masking_reason;
    return false;
#else
    _stage("Generating masks (segmentation)");

    std::error_code ec;
    // Never the masks themselves: a nested masks/ is only possible for photos
    // read where they are, and those already have masks -- but the guard costs
    // nothing and the alternative is masking a folder of masks.
    const std::vector<fs::path> files = walk_images(images, masks);
    if (files.empty()) { error = "no images to mask"; return false; }

    sam::MaskOptions mo;
    mo.model = job.mask_model_path;
    mo.text = job.mask_prompt;
    mo.neg_text = job.mask_negative_prompt;
    mo.keep_prompted = job.mask_keep_subject;
    mo.max_size = job.mask_max_image_size;
    // A folder of photos is not a video: consecutive files may be anywhere in
    // the scene, so tracking across them would carry a memory bank that does
    // not apply. Clicks are the exception and force it on -- a click says
    // nothing about any frame but its own, so without a memory bank to carry
    // the object forward there is nothing to propagate, and the input in that
    // case is an ordered capture (this path also masks frames that ffmpeg
    // extracted, which is where a clicked object usually arrives).
    const bool clicks_here = clicks_apply_to(job, in);
    mo.video = clicks_here && !job.mask_clicks.empty();
    // For a folder of photos the preview counted the same files in the same
    // order, so a click's frame index is exact. For frames ffmpeg extracted it
    // is not -- the preview measured against the video and ffmpeg resampled it
    // to a different rate -- and the fraction through the capture is the part
    // that survives, so the index is recomputed from it here.
    if (clicks_here)
        mo.seeds = seeds_from_clicks(job.mask_clicks,
                                     in.is_video ? (int64_t)files.size() : 0);

    sam::Masker masker;
    if (!masker.init(mo, error)) return false;

    const fs::path mask_root(masks);
    fs::create_directories(mask_root, ec);
    const fs::path image_root(images);

    // Encoding a 1080p mask costs about a third of what the model costs to
    // produce it, and none of it needs the GPU.
    app::WriterPool writers;
    int done = 0;
    int64_t index = -1;
    RateLimitedProgress progress(_log, "images masked", (int64_t)files.size());
    for (const fs::path& f : files) {
        ++index;
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
        if (!masker.run(img, mask, nullptr, index)) {
            error = "masking failed on " + f.filename().string() + ": " +
                    masker.lastError();
            return false;
        }
        app::WriteJob wj;
        wj.mask = std::move(mask);
        wj.path = dst.string();
        writers.submit(std::move(wj));
        progress.update(++done, done == (int)files.size());
    }
    // The last few masks are still in the queue; the caller goes straight on to
    // structure from motion, which will look for them.
    writers.finish();
    if (writers.failures() > 0) {
        error = std::to_string(writers.failures()) + " mask(s) could not be "
                "written to " + mask_root.string();
        return false;
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
                                        const std::string& masks_rel,
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
    const fs::path script = ws / ".spirula_mask.py";
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
        "--masks", masks_rel,
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
    const bool have_masks = fs::is_directory(ws / masks_rel, ec) &&
                            !fs::is_empty(ws / masks_rel, ec);
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
