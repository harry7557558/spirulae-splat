// SfmRunner.cpp -- see SfmRunner.h.

#include "app/gui/SfmRunner.h"

#include "i18n/Locale.h"
#include "i18n/catalog/Log.h"

#include "app/gui/AppPaths.h"
#include "app/gui/Subprocess.h"
#ifdef SS_TOOL_SFM
// For the stage tags the child prints; a build without the module has no child
// to read (see availability()).
#include "sfm/core/Log.h"
#endif

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <sstream>

namespace fs = std::filesystem;
namespace lmsg = spirula::i18n::msg::log;
using spirula::i18n::format;

// Shorthand: every log line that carries a path or a count is a format()
// call, and there are enough of them here to be worth a short name.
inline std::string fmt(const spirula::i18n::Msg& m,
                       std::initializer_list<spirula::i18n::Arg> a) {
    return format(m, a);
}

namespace gui {

namespace {

const char* kQuality[] = {"low", "medium", "high", "extreme"};
const char* kDataType[] = {"individual", "video", "internet"};
const char* kCameraMode[] = {"single", "folder", "image"};
const char* kPairs[] = {"auto", "exhaustive", "sequential", "prefilter"};
const char* kMapper[] = {"flat", "bottom-up"};
const char* kFeatures[] = {"sift", "aliked-n16rot", "aliked-n32"};
const char* kMatcher[] = {"bruteforce", "lightglue"};

template <int N>
const char* pick(const char* const (&table)[N], int i, int fallback = 0) {
    return table[(i >= 0 && i < N) ? i : fallback];
}

// Splits a free-form flag string the way a shell would for the simple cases:
// whitespace separated, with "quoted runs" kept together. Enough for
// `--max-error 2 --masks "/path/with space"`.
std::vector<std::string> split_args(const std::string& s) {
    std::vector<std::string> out;
    std::string cur;
    bool in_quote = false, have = false;
    char quote = 0;
    for (char c : s) {
        if (in_quote) {
            if (c == quote) in_quote = false;
            else cur += c;
        } else if (c == '"' || c == '\'') {
            in_quote = true;
            quote = c;
            have = true;
        } else if (std::isspace((unsigned char)c)) {
            if (have) { out.push_back(cur); cur.clear(); have = false; }
        } else {
            cur += c;
            have = true;
        }
    }
    if (have) out.push_back(cur);
    return out;
}

// The mapper only writes a model when it finishes, so any model on disk is
// from a completed run.
bool has_model(const fs::path& sparse) {
    std::error_code ec;
    for (fs::directory_iterator it(sparse, ec), end; !ec && it != end;
         it.increment(ec))
        if (fs::exists(it->path() / "cameras.bin", ec)) return true;
    return false;
}

}  // namespace

std::string SfmRunner::availability() {
#ifndef SS_TOOL_SFM
    return "this build has no structure-from-motion module "
           "(-DSS_BUILD_SFM=OFF); use COLMAP instead.";
#else
    if (exe_path().empty())
        return "this program could not work out its own path, so it cannot "
               "run the reconstruction step.";
    return "";
#endif
}

SfmRunner::~SfmRunner() {
    cancel();
    if (_worker.joinable()) _worker.join();
}

void SfmRunner::start(const SfmJob& job) {
    if (_state.load() == State::Running) return;
    if (_worker.joinable()) _worker.join();
    _cancel = false;
    _progress = -1.0f;
    _partial = false;
    {
        std::lock_guard<std::mutex> lk(_mu);
        _error.clear();
        _dataset_dir.clear();
        _image_dir.clear();
        _mask_dir.clear();
    }
    _state = State::Running;
    _worker = std::thread([this, job] { run(job); });
}

void SfmRunner::cancel() { _cancel = true; }

std::string SfmRunner::stage() {
    std::lock_guard<std::mutex> lk(_mu);
    return _stage;
}
std::string SfmRunner::error() {
    std::lock_guard<std::mutex> lk(_mu);
    return _error;
}
std::string SfmRunner::dataset_dir() {
    std::lock_guard<std::mutex> lk(_mu);
    return _dataset_dir;
}
std::string SfmRunner::image_dir() {
    std::lock_guard<std::mutex> lk(_mu);
    return _image_dir;
}
std::string SfmRunner::mask_dir() {
    std::lock_guard<std::mutex> lk(_mu);
    return _mask_dir;
}
std::vector<std::string> SfmRunner::drain_log() {
    std::lock_guard<std::mutex> lk(_mu);
    std::vector<std::string> out;
    out.swap(_log);
    return out;
}

void SfmRunner::log(const std::string& line) {
    std::lock_guard<std::mutex> lk(_mu);
    _log.push_back(line);
    if (_log.size() > 5000) _log.erase(_log.begin(), _log.begin() + 1000);
}

void SfmRunner::set_stage_if_new(const char* s) {
    {
        std::lock_guard<std::mutex> lk(_mu);
        if (_stage == s) return;
    }
    set_stage(s);
}

void SfmRunner::set_stage(const std::string& s) {
    {
        std::lock_guard<std::mutex> lk(_mu);
        _stage = s;
    }
    _progress = -1.0f;
    log("==== " + s + " ====");
}

// spirula-sfm's own progress lines. Reading them is a little grubby, but it is
// the one thing a child process cannot hand over structurally.
//
// What it keys on is the stage TAG, and the tag is translated -- so this asks
// the same catalog the child printed from (sfm/core/Log.h), rather than
// matching English text that a `--lang ja` run will never emit. The numbers
// are the first "a/b" on the line, which every progress message puts there and
// no other message on those stages has.
void SfmRunner::note_progress(const std::string& l) {
#ifdef SS_TOOL_SFM
    using sfm::slog::Tag;
    auto tagged = [&l](Tag t) {
        const std::string p = sfm::slog::prefix(t);
        return l.compare(0, p.size(), p) == 0;
    };
    // First "<digits>/<digits>" anywhere in the line.
    auto fraction = [&l](unsigned long& a, unsigned long& b) {
        for (size_t i = 0; i + 2 < l.size(); i++) {
            if (!std::isdigit((unsigned char)l[i])) continue;
            if (i && std::isdigit((unsigned char)l[i - 1])) continue;
            if (std::sscanf(l.c_str() + i, "%lu/%lu", &a, &b) == 2 && b > 0) return true;
            while (i + 1 < l.size() && std::isdigit((unsigned char)l[i + 1])) i++;
        }
        return false;
    };
    unsigned long a = 0, b = 0;
    if (tagged(Tag::Extract)) {
        if (!fraction(a, b)) return;
        set_stage_if_new(lmsg::stage_finding_features.get());
        _progress = (float)((double)a / (double)b);
        return;
    }
    if (tagged(Tag::Match)) {
        if (!fraction(a, b)) return;
        set_stage_if_new(lmsg::stage_matching_images.get());
        _progress = (float)((double)a / (double)b);
        return;
    }
    if (tagged(Tag::Map)) {
        // The mapper counts registrations rather than working towards a
        // total, so this is a stage change and a spinner, not a fraction.
        set_stage_if_new(lmsg::stage_reconstructing.get());
        _progress = -1.0f;
    }
#else
    (void)l;
#endif
}

// A capture shot as several inputs is several cameras, and the lens is a
// property of the input rather than of the dataset. `spirula sfm` says that with
// PREFIX=VALUE, where the prefix matches an image name by path prefix -- which
// is exactly the sub-folder each input's frames went into. A dual-lens video's
// cam0/ and cam1/ both sit under that folder, so one prefix covers both, and
// they still group as two cameras because grouping is per folder.
//
// The focal length is carried as a fraction of the image width rather than in
// pixels, because the width is not known until the frames exist (an Insta360
// X5 is fx = fy ~ 0.269 * width whatever it was shot at). A lone input has no
// prefix, so its value IS the dataset-wide one -- and an explicit focal typed
// into the panel wins over a preset-derived one.
void SfmRunner::append_camera_overrides(const SfmJob& job, const PrepResult& prep,
                                        std::vector<std::string>& argv) {
    for (const PrepInput& in : job.prep.inputs) {
        const std::string prefix = in.subdir.empty() ? "" : in.subdir + "=";
        // For a lone input the panel's own "Camera / lens" is the single source
        // of truth and has already been passed; only a named group adds one.
        if (!in.subdir.empty() && !in.camera_model.empty()) {
            argv.push_back("--camera-model");
            argv.push_back(prefix + in.camera_model);
        }
        if (!(in.focal_factor > 0)) continue;
        if (in.subdir.empty() && job.init_focal_px > 0) continue;
        const std::string dir =
            (in.subdir.empty() ? fs::path(prep.image_dir)
                               : fs::path(prep.image_dir) / in.subdir).string();
        int W = 0, H = 0;
        if (!DatasetPrep::first_image_dims(dir, W, H)) {
            log(fmt(lmsg::sfm_focal_unreadable, {dir}));
            continue;
        }
        char buf[32];
        std::snprintf(buf, sizeof buf, "%g", (double)in.focal_factor * W);
        argv.push_back("--focal");
        argv.push_back(prefix + buf);
        log(fmt(lmsg::sfm_initial_focal,
                {in.subdir.empty() ? lmsg::sfm_the_capture.get() : in.subdir.c_str(),
                 buf, in.focal_factor, (long long)W}));
    }
}

void SfmRunner::run(SfmJob job) {
    auto fail = [&](const std::string& why) {
        std::lock_guard<std::mutex> lk(_mu);
        _error = why;
        _state = _cancel.load() ? State::Cancelled : State::Failed;
    };

    try {
        const fs::path ws = job.prep.workspace;
        std::error_code ec;
        fs::create_directories(ws, ec);

        if (std::string why = availability(); !why.empty()) return fail(why);

        // What was there before this run touched anything. Not a reason to
        // refuse: an existing sparse/ (see WorkspaceState::model), nor the
        // input's own images -- writing the dataset next to the images it was
        // built from is the layout every parser expects.
        const WorkspaceState prior = probe_workspace(ws.string(), job.prep.inputs);
        if (prior.resumable() && !job.prep.resume)
            return fail(lmsg::err_unfinished_run.get());
        if (prior.resumable())
            log(fmt(lmsg::sfm_resuming, {ws.string()}));
        if (prior.model)
            log(fmt(lmsg::sfm_will_overwrite, {(ws / "sparse").string()}));

        // ---- 1. frames and masks ------------------------------------------
        PrepResult prep;
        {
            DatasetPrep dp([this](const std::string& l) { log(l); },
                           [this](const std::string& s) { set_stage(s); },
                           _cancel);
            std::string err;
            if (!dp.run(job.prep, prep, err)) return fail(err);
        }
        if (prep.per_folder_cameras && job.camera_mode == 0) {
            log(lmsg::one_camera_per_folder.get());
            job.camera_mode = 1;
        }

        // ---- 2. reconstruction --------------------------------------------
        // Resume: a completed model means the mapper already ran -- but only
        // when this run's own leftovers are there to say the model is ours. A
        // sparse/ on its own is somebody else's dataset (or a finished one),
        // and skipping the mapper for it would "reconstruct" by doing nothing.
        if (job.prep.resume && prior.features && has_model(ws / "sparse")) {
            log(lmsg::sfm_resume_skip_recon.get());
        } else {
            set_stage(lmsg::stage_reconstructing_features.get());
            std::vector<std::string> argv = {
                // The child is this same executable, so it has the same
                // thirteen languages -- tell it which one, or its output
                // lands in the log in whatever the machine's locale is.
                exe_path(), "--lang", spirula::i18n::code(spirula::i18n::current()),
                "sfm", "auto", prep.image_dir,
                "-o", ws.string(),
                "--quality", pick(kQuality, job.quality, 2),
                "--data-type", pick(kDataType, job.data_type),
                "--camera-model", job.camera_model,
                "--camera-mode", pick(kCameraMode, job.camera_mode, 1),
                "--mapper", pick(kMapper, job.mapper),
                "--features", pick(kFeatures, job.features),
                // LightGlue only exists for the learned descriptors; asking
                // for it with SIFT selected is a usage error, so do not.
                "--matcher",
                pick(kMatcher, job.features == 0 ? 0 : job.matcher),
            };
            if (job.pairs > 0) {
                argv.push_back("--pairs");
                argv.push_back(pick(kPairs, job.pairs));
                if (job.pairs == 2) {
                    argv.push_back("--overlap");
                    argv.push_back(std::to_string(job.overlap));
                }
            }
            // Sequential is what `auto` resolves to for video, so this has to
            // be passed whenever sequential is reachable, not only when it was
            // named. It is a no-op under the other pair modes.
            if (!job.loop_closure) argv.push_back("--no-loop-closure");
            if (job.init_focal_px > 0) {
                char buf[32];
                std::snprintf(buf, sizeof buf, "%g", job.init_focal_px);
                argv.push_back("--focal");
                argv.push_back(buf);
            }
            append_camera_overrides(job, prep, argv);
            if (job.max_features > 0) {
                // Each frontend has its own count flag, because their budgets
                // are not comparable -- a learned detector emits a few
                // thousand better-localized points where SIFT wants tens of
                // thousands. One spinner, routed to whichever is running.
                argv.push_back(job.features == 0 ? "--max-features"
                                                 : "--aliked-max-features");
                argv.push_back(std::to_string(job.max_features));
            }
            if (job.max_image_size > 0) {
                argv.push_back("--max-image-size");
                argv.push_back(std::to_string(job.max_image_size));
            }
            if (!prep.mask_dir.empty()) {
                argv.push_back("--masks");
                argv.push_back(prep.mask_dir);
            } else {
                // Otherwise `auto` picks up a stale masks/ sitting beside the
                // images from an earlier run with masking on.
                argv.push_back("--no-masks");
            }
            for (const std::string& a : split_args(job.extra_args))
                argv.push_back(a);

            std::string cmd;
            for (const auto& a : argv) cmd += (cmd.empty() ? "$ " : " ") + a;
            log(cmd);
            const int rc = run_process(argv, "", [this](const std::string& l) {
                log(l);
                note_progress(l);
            }, _cancel);
            if (rc == kCancelled) return fail(lmsg::err_cancelled.get());
            if (rc == kSpawnFailed)
                return fail(fmt(lmsg::err_spawn_recon, {argv[0]}));
            // `auto` spends exit code 2 on "nothing reconstructed" and 3 on
            // "reconstructed, but under half the images registered or the
            // reprojection error is high" (src/sfm/README.md). 3 is a warning
            // here, not a failure: a partial model still trains, and throwing
            // it away over a threshold would be worse than saying so.
            if (rc == 3) {
                _partial = true;
                log(lmsg::sfm_partial.get());
            } else if (rc != 0) {
                return fail(lmsg::err_recon_failed.get());
            }
        }

        if (!has_model(ws / "sparse"))
            return fail(lmsg::err_no_reconstruction.get());

        // ---- 3. tidy up ----------------------------------------------------
        if (!job.keep_intermediate) {
            set_stage(lmsg::stage_cleaning_up.get());
            for (const char* name : {"features", "matches.bin"}) {
                const fs::path p = ws / name;
                if (!fs::exists(p, ec)) continue;
                if (fs::is_directory(p, ec)) {
                    for (fs::directory_iterator it(p, ec), end; !ec && it != end;
                         it.increment(ec))
                        fs::remove(it->path(), ec);
                }
                fs::remove(p, ec);
            }
        }

        if (reads_photos_in_place(job.prep.inputs))
            log(fmt(lmsg::photos_referenced_in_place, {prep.image_dir_cfg}));

        set_stage(lmsg::stage_done.get());
        _progress = 1.0f;
        {
            std::lock_guard<std::mutex> lk(_mu);
            _dataset_dir = ws.string();
            _image_dir = prep.image_dir_cfg;
            _mask_dir = prep.mask_dir_cfg;
        }
        _state = State::Done;
    } catch (const std::exception& e) {
        fail(e.what());
    }
}

}  // namespace gui
