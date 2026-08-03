// SfmRunner.cpp -- see SfmRunner.h.

#include "app/gui/SfmRunner.h"

#include "app/gui/AppPaths.h"
#include "app/gui/Subprocess.h"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <sstream>

namespace fs = std::filesystem;

namespace gui {

namespace {

const char* kQuality[] = {"low", "medium", "high", "extreme"};
const char* kDataType[] = {"individual", "video", "internet"};
const char* kCameraMode[] = {"single", "folder", "image"};
const char* kPairs[] = {"auto", "exhaustive", "sequential", "prefilter"};
const char* kMapper[] = {"auto", "flat", "bottom-up"};
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
#ifndef SSPLAT_TOOL_SFM
    return "this build has no structure-from-motion module "
           "(-DSSPLAT_BUILD_SFM=OFF); use COLMAP instead.";
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

// ssplat-sfm's own progress lines. Reading them is a little grubby, but it is
// the one thing a child process cannot hand over structurally, and the formats
// below are the ones its stages print on every run.
void SfmRunner::note_progress(const std::string& l) {
    auto frac = [&](size_t a, size_t b) {
        if (b > 0) _progress = (float)((double)a / (double)b);
    };
    unsigned long a = 0, b = 0;
    if (l.size() > 1 && l[0] == '[' && std::isdigit((unsigned char)l[1]) &&
        std::sscanf(l.c_str(), "[%lu/%lu]", &a, &b) == 2) {
        set_stage_if_new("Finding features");
        frac(a, b);
        return;
    }
    if (l.rfind("[match]", 0) == 0 &&
        std::sscanf(l.c_str(), "[match] %lu/%lu pairs", &a, &b) == 2) {
        set_stage_if_new("Matching images");
        frac(a, b);
        return;
    }
    if (l.rfind("[map]", 0) == 0) {
        set_stage_if_new("Reconstructing cameras (the slow part)");
        // "[map] registered image 42 (118/193 PnP inliers), 57 total"
        const char* p = std::strstr(l.c_str(), " total");
        if (p) _progress = -1.0f;   // count, not a fraction: show a spinner
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

        {
            const bool prior = fs::is_directory(ws / "sparse", ec) ||
                               fs::is_directory(ws / "features", ec) ||
                               (fs::is_directory(ws / "images", ec) &&
                                !fs::is_empty(ws / "images", ec));
            if (prior && !job.prep.resume)
                return fail("the output folder already holds a previous run "
                            "(images / features / sparse); tick \"Resume "
                            "previous run\" to reuse it, or pick an empty "
                            "folder");
            if (prior)
                log("Resuming the previous run in " + ws.string());
        }

        // ---- 1. frames and masks ------------------------------------------
        PrepResult prep;
        {
            DatasetPrep dp([this](const std::string& l) { log(l); },
                           [this](const std::string& s) { set_stage(s); },
                           _cancel);
            std::string err;
            if (!dp.run(job.prep, prep, err)) return fail(err);
        }
        if (prep.multi_track && job.camera_mode == 0) {
            log("Multi-track video: switching to one camera per folder");
            job.camera_mode = 1;
        }

        // ---- 2. reconstruction --------------------------------------------
        // Resume: a completed model means the mapper already ran.
        if (job.prep.resume && has_model(ws / "sparse")) {
            log("Resume: a reconstruction already exists under sparse/; "
                "skipping (delete it to reconstruct again)");
        } else {
            set_stage("Reconstructing (finding features)");
            std::vector<std::string> argv = {
                exe_path(), "sfm", "auto", prep.image_dir,
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
            if (job.init_focal_px > 0) {
                char buf[32];
                std::snprintf(buf, sizeof buf, "%g", job.init_focal_px);
                argv.push_back("--focal");
                argv.push_back(buf);
            }
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
            if (rc == kCancelled) return fail("cancelled");
            if (rc == kSpawnFailed)
                return fail("could not start the reconstruction (" + argv[0] + ")");
            // `auto` spends exit code 2 on "nothing reconstructed" and 3 on
            // "reconstructed, but under half the images registered or the
            // reprojection error is high" (src/sfm/README.md). 3 is a warning
            // here, not a failure: a partial model still trains, and throwing
            // it away over a threshold would be worse than saying so.
            if (rc == 3) {
                _partial = true;
                log("Note: only part of the capture reconstructed. It will "
                    "still train, but expect gaps.");
            } else if (rc != 0) {
                return fail("reconstruction failed (see the log). Common "
                            "causes: too few overlapping images, not enough "
                            "overlap between them, or the wrong camera model "
                            "for the lens.");
            }
        }

        if (!has_model(ws / "sparse"))
            return fail("no reconstruction was produced -- the images may not "
                        "overlap enough. Try more photos around the subject, "
                        "or a higher quality setting.");

        // ---- 3. tidy up ----------------------------------------------------
        if (!job.keep_intermediate) {
            set_stage("Cleaning up");
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

        if (!job.prep.is_video)
            log("Note: the photos are referenced where they are. If you reopen "
                "this dataset later, set image_dir to " + prep.image_dir_cfg +
                " under the dataset options.");

        set_stage("Done");
        _progress = 1.0f;
        {
            std::lock_guard<std::mutex> lk(_mu);
            _dataset_dir = ws.string();
            _image_dir = prep.image_dir_cfg;
        }
        _state = State::Done;
    } catch (const std::exception& e) {
        fail(e.what());
    }
}

}  // namespace gui
