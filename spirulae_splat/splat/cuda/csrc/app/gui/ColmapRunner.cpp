// ColmapRunner.cpp -- see ColmapRunner.h. CLI flags mirror
// scripts/run_colmap.bash (COLMAP >= 4.x; use_gpu-style flags are gone).

#include "ColmapRunner.h"
#include "FrameSelect.h"
#include "Subprocess.h"

#include "app_generated/mask_py.h"   // kMaskPy[] (CMake-embedded scripts/mask.py)

#ifndef _WIN32
#include <ftw.h>
#endif

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>

namespace fs = std::filesystem;

namespace gui {

namespace {

// NOT std::filesystem::remove_all -- on the torch build libtorch.so
// interposes an ABI-incompatible copy (see app/README.md gotchas).
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

const char* kVocabTreeName = "vocab_tree_faiss_flickr100K_words256K.bin";
const char* kVocabTreeUrl =
    "https://github.com/colmap/colmap/releases/download/3.11.1/"
    "vocab_tree_faiss_flickr100K_words256K.bin";

bool is_image_file(const fs::path& p) {
    std::string e = p.extension().string();
    for (auto& c : e) c = (char)std::tolower((unsigned char)c);
    return e == ".jpg" || e == ".jpeg" || e == ".png" || e == ".webp" ||
           e == ".tif" || e == ".tiff" || e == ".bmp";
}

// Recursive: COLMAP's feature_extractor indexes image_path recursively, so
// the sanity count must match.
int count_images(const std::string& dir) {
    int n = 0;
    std::error_code ec;
    for (fs::recursive_directory_iterator it(
             dir, fs::directory_options::skip_permission_denied, ec), end;
         !ec && it != end; it.increment(ec))
        if (it->is_regular_file(ec) && is_image_file(it->path())) n++;
    return n;
}

fs::path cache_dir() {
#ifdef _WIN32
    const char* base = std::getenv("LOCALAPPDATA");
    fs::path dir = base ? fs::path(base) : fs::path(".");
    dir /= "spirulae-splat";
#else
    fs::path dir;
    if (const char* x = std::getenv("XDG_CACHE_HOME")) dir = x;
    else if (const char* h = std::getenv("HOME")) dir = fs::path(h) / ".cache";
    else dir = ".";
    dir /= "spirulae-splat";
#endif
    std::error_code ec;
    fs::create_directories(dir, ec);
    return dir;
}

}  // namespace

ColmapRunner::~ColmapRunner() {
    cancel();
    if (_worker.joinable()) _worker.join();
}

void ColmapRunner::start(const ColmapJob& job) {
    if (_state.load() == State::Running) return;
    if (_worker.joinable()) _worker.join();
    _cancel = false;
    {
        std::lock_guard<std::mutex> lk(_mu);
        _error.clear();
        _dataset_dir.clear();
        _image_dir.clear();
    }
    _state = State::Running;
    _worker = std::thread([this, job] { run(job); });
}

void ColmapRunner::cancel() { _cancel = true; }

std::string ColmapRunner::stage() {
    std::lock_guard<std::mutex> lk(_mu);
    return _stage;
}
std::string ColmapRunner::error() {
    std::lock_guard<std::mutex> lk(_mu);
    return _error;
}
std::string ColmapRunner::dataset_dir() {
    std::lock_guard<std::mutex> lk(_mu);
    return _dataset_dir;
}
std::string ColmapRunner::image_dir() {
    std::lock_guard<std::mutex> lk(_mu);
    return _image_dir;
}
std::vector<std::string> ColmapRunner::drain_log() {
    std::lock_guard<std::mutex> lk(_mu);
    std::vector<std::string> out;
    out.swap(_log);
    return out;
}

void ColmapRunner::log(const std::string& line) {
    std::lock_guard<std::mutex> lk(_mu);
    _log.push_back(line);
    if (_log.size() > 5000) _log.erase(_log.begin(), _log.begin() + 1000);
}

void ColmapRunner::set_stage(const std::string& s) {
    {
        std::lock_guard<std::mutex> lk(_mu);
        _stage = s;
    }
    log("==== " + s + " ====");
}

int ColmapRunner::exec(const std::vector<std::string>& argv) {
    std::string cmd;
    for (const auto& a : argv) cmd += (cmd.empty() ? "$ " : " ") + a;
    log(cmd);
    return run_process(argv, "", [this](const std::string& l) { log(l); }, _cancel);
}

// COLMAP >= 4 required: 3.x still had the SiftExtraction.use_gpu-era CLI and
// misses flags this pipeline passes (run_colmap.bash targets 4.x).
bool ColmapRunner::check_colmap_version(const ColmapJob& job, std::string& err) {
    if (!command_exists(job.colmap_exe)) {
        err = "colmap not found ('" + job.colmap_exe +
              "'); install COLMAP 4.x or set its path under Tool Locations";
        return false;
    }
    int major = -1, minor = -1;
    run_process({job.colmap_exe, "help"}, "",
                [&](const std::string& l) {
                    if (major >= 0) return;
                    const char* p = std::strstr(l.c_str(), "COLMAP ");
                    if (p) std::sscanf(p, "COLMAP %d.%d", &major, &minor);
                },
                _cancel);
    if (major < 0) {
        log("warning: could not detect the COLMAP version; continuing anyway");
        return true;
    }
    log("Detected COLMAP " + std::to_string(major) + "." + std::to_string(minor));
    if (major < 4) {
        err = "COLMAP " + std::to_string(major) + "." + std::to_string(minor) +
              " found, but version 4.x or newer is required (the command-line "
              "options changed). Please upgrade COLMAP.";
        return false;
    }
    return true;
}

// Resolve the vocab-tree file: explicit path > found near the workspace or
// in the app cache > downloaded into the cache (curl).
std::string ColmapRunner::resolve_vocab_tree(const ColmapJob& job) {
    std::error_code ec;
    if (!job.vocab_tree_path.empty()) {
        if (fs::exists(job.vocab_tree_path, ec)) return job.vocab_tree_path;
        log("vocab tree not found at '" + job.vocab_tree_path + "'");
        return "";
    }
    fs::path ws = job.workspace;
    for (const fs::path& dir : {ws, ws.parent_path(), cache_dir()}) {
        if (!fs::is_directory(dir, ec)) continue;
        for (fs::directory_iterator it(dir, ec), end; !ec && it != end;
             it.increment(ec)) {
            std::string name = it->path().filename().string();
            if (name.rfind("vocab_tree", 0) == 0 &&
                it->path().extension() == ".bin") {
                log("Found vocabulary tree: " + it->path().string());
                return it->path().string();
            }
        }
    }
    // Download into the cache.
    fs::path dst = cache_dir() / kVocabTreeName;
    set_stage("Downloading vocabulary tree (one-time, ~150 MB)");
    if (!command_exists("curl")) {
        log("curl not found -- download it manually:");
        log(std::string("  ") + kVocabTreeUrl + " -> " + dst.string());
        return "";
    }
    fs::path tmp = dst;
    tmp += ".part";
    int rc = exec({"curl", "-L", "-f", "--progress-bar",
                   "-o", tmp.string(), kVocabTreeUrl});
    if (rc != 0) {
        fs::remove(tmp, ec);
        log("vocabulary tree download failed");
        return "";
    }
    fs::rename(tmp, dst, ec);
    return ec ? "" : dst.string();
}

// AI masking via the embedded scripts/mask.py. The script prints an install
// hint and exits 0 when lang-sam / SAM-3 is missing, so detect that from
// its output rather than the exit code.
bool ColmapRunner::run_masking(const ColmapJob& job, const std::string& images,
                               std::string& err) {
    set_stage("Generating masks (AI segmentation)");
    if (!command_exists(job.python_exe)) {
        err = "Python not found ('" + job.python_exe + "'); masking needs "
              "Python with the lang-segment-anything package. Set the Python "
              "path under Tool Locations, or disable masking.";
        return false;
    }
    fs::path ws = job.workspace;
    fs::path script = ws / ".ssplat_mask.py";
    {
        FILE* f = std::fopen(script.string().c_str(), "wb");
        if (!f) { err = "cannot write " + script.string(); return false; }
        std::fwrite(kMaskPy, 1, kMaskPySize, f);
        std::fclose(f);
    }
    std::vector<std::string> argv = {
        job.python_exe, script.string(), ws.string(),
        "--prompt", job.mask_prompt,
        "--images", images,
        "--masks", "masks",
        "--max_image_size", std::to_string(job.mask_max_image_size),
        "--model", job.mask_model,
    };
    if (!job.mask_negative_prompt.empty()) {
        argv.push_back("--negative_prompt");
        argv.push_back(job.mask_negative_prompt);
    }
    std::string install_hint;
    std::string cmd;
    for (const auto& a : argv) cmd += (cmd.empty() ? "$ " : " ") + a;
    log(cmd);
    int rc = run_process(argv, "", [&](const std::string& l) {
        log(l);
        if (l.find("not found or not installed properly") != std::string::npos ||
            l.find("ModuleNotFoundError") != std::string::npos)
            install_hint = l;
    }, _cancel);
    if (rc == kCancelled) { err = "cancelled"; return false; }
    std::error_code ec;
    bool have_masks = fs::is_directory(ws / "masks", ec) &&
                      !fs::is_empty(ws / "masks", ec);
    if (!install_hint.empty() || rc != 0 || !have_masks) {
        err = "Mask generation failed";
        if (!install_hint.empty())
            err += " -- missing Python packages. Install lang-segment-anything "
                   "(pip install git+https://github.com/luca-medeiros/"
                   "lang-segment-anything, needs CUDA PyTorch)"
                   + std::string(job.mask_model == "sam3"
                       ? ", or for SAM-3: https://github.com/facebookresearch/sam3"
                       : "");
        else
            err += " (see log)";
        return false;
    }
    return true;
}

void ColmapRunner::run(ColmapJob job) {
    auto fail = [&](const std::string& why) {
        std::lock_guard<std::mutex> lk(_mu);
        _error = why;
        _state = _cancel.load() ? State::Cancelled : State::Failed;
    };
    try {
        const fs::path ws = job.workspace;
        fs::create_directories(ws);

        std::string err;
        if (!check_colmap_version(job, err)) return fail(err);

        // ---- 1. resolve the image directory --------------------------------
        std::string images;
        std::string image_dir_cfg;   // what the trainer's image_dir should be
        if (job.is_video) {
            if (!command_exists(job.ffmpeg_exe))
                return fail("ffmpeg not found ('" + job.ffmpeg_exe +
                            "'); install it or set its path under Tool Locations");

            // Multi-track videos (Insta360 .insv): one folder per track,
            // one COLMAP camera per folder (extract_frames.py behavior).
            std::vector<int> streams = {0};
            {
                std::string ext = fs::path(job.input_path).extension().string();
                for (auto& c : ext) c = (char)std::tolower((unsigned char)c);
                if (ext == ".insv") {
                    std::vector<int> found;
                    run_process({"ffprobe", "-v", "error", "-select_streams", "v",
                                 "-show_entries", "stream=index", "-of", "csv=p=0",
                                 job.input_path}, "",
                                [&](const std::string& l) {
                                    try { found.push_back(std::stoi(l)); } catch (...) {}
                                }, _cancel);
                    if (found.size() > 1) streams = found;
                }
            }
            if (streams.size() > 1 && job.camera_mode == 0) {
                log("Multi-track video: switching to one camera per folder");
                job.camera_mode = 1;
            }

            int window = std::max(job.sharp_window, 1);
            for (size_t tr = 0; tr < streams.size(); tr++) {
                std::string track_path = job.input_path;
                fs::path out_dir = ws / "images";
                if (streams.size() > 1) {
                    set_stage("Splitting video track " + std::to_string(tr));
                    track_path = (ws / ("track_cam" + std::to_string(tr) + ".mp4")).string();
                    int rc = exec({job.ffmpeg_exe, "-nostdin", "-y",
                                   "-i", job.input_path,
                                   "-map", "0:v:" + std::to_string(tr),
                                   "-c", "copy", track_path});
                    if (rc == kCancelled) return fail("cancelled");
                    if (rc != 0) return fail("ffmpeg track split failed (see log)");
                    out_dir = ws / "images" / ("cam" + std::to_string(tr));
                }

                set_stage(window > 1
                    ? "Extracting candidate frames (ffmpeg)"
                    : "Extracting frames (ffmpeg)");
                fs::path cand = ws / "frames_tmp";
                remove_tree(cand);
                fs::create_directories(cand);
                char vf[64];
                std::snprintf(vf, sizeof vf, "fps=%g",
                              (double)job.video_fps * window);
                int rc = exec({job.ffmpeg_exe, "-nostdin", "-y", "-i", track_path,
                               "-vf", vf, "-qscale:v", "2",
                               (cand / "c_%06d.jpg").string()});
                if (rc == kCancelled) return fail("cancelled");
                if (rc != 0) return fail("ffmpeg frame extraction failed (see log)");

                if (window > 1)
                    set_stage("Selecting sharpest frames (multithreaded)");
                fs::create_directories(out_dir);
                int kept = select_sharpest_frames(
                    cand.string(), out_dir.string(), "",
                    window, job.max_frames,
                    [this](const std::string& l) { log(l); }, _cancel);
                remove_tree(cand);
                if (streams.size() > 1) {
                    std::error_code ec;
                    fs::remove(track_path, ec);
                }
                if (kept < 0)
                    return fail(_cancel.load() ? "cancelled"
                                               : "frame selection failed");
                log("Kept " + std::to_string(kept) + " frames -> " +
                    out_dir.string());
            }
            images = (ws / "images").string();
            image_dir_cfg = "images";
        } else {
            images = job.input_path;
            // Reference the source images in place (no copy); the parsers
            // accept an absolute image_dir.
            image_dir_cfg = fs::absolute(job.input_path).string();
        }
        int n_images = count_images(images);
        log("Found " + std::to_string(n_images) + " images in " + images);
        if (n_images < 3)
            return fail("need at least 3 images (found " +
                        std::to_string(n_images) + ")");

        // ---- 2. optional AI masking ----------------------------------------
        bool have_masks = false;
        if (job.mask_enable) {
            if (job.mask_prompt.empty())
                return fail("masking is enabled but the prompt is empty "
                            "(e.g. \"people; cars\")");
            if (!run_masking(job, job.is_video ? "images" : image_dir_cfg, err))
                return fail(err);
            have_masks = true;
        }

        // Quality knobs (per run_colmap.bash: fewer features = much faster
        // matching; O(n^2) in feature count).
        int features = job.max_num_features > 0 ? job.max_num_features
                     : job.quality == 0 ? 4096
                     : job.quality == 1 ? 8192 : 16384;
        const std::string db = (ws / "database.db").string();

        // ---- 3. feature extraction -----------------------------------------
        set_stage("Extracting features (colmap)");
        std::vector<std::string> fe = {job.colmap_exe, "feature_extractor",
            "--database_path", db,
            "--image_path", images,
            "--ImageReader.camera_model", job.camera_model,
            "--SiftExtraction.max_num_features", std::to_string(features)};
        if (job.camera_mode == 0) {
            fe.push_back("--ImageReader.single_camera");
            fe.push_back("1");
        } else if (job.camera_mode == 1) {
            fe.push_back("--ImageReader.single_camera_per_folder");
            fe.push_back("1");
        }
        int size_cap = job.max_image_size > 0 ? job.max_image_size
                     : job.quality == 0 ? 2000 : 0;
        if (size_cap > 0) {
            fe.push_back("--FeatureExtraction.max_image_size");
            fe.push_back(std::to_string(size_cap));
        }
        if (job.estimate_affine_shape) {
            fe.push_back("--SiftExtraction.estimate_affine_shape");
            fe.push_back("1");
        }
        if (have_masks) {
            fe.push_back("--ImageReader.mask_path");
            fe.push_back((ws / "masks").string());
        }
        int rc = exec(fe);
        if (rc == kCancelled) return fail("cancelled");
        if (rc != 0) return fail("colmap feature_extractor failed (see log)");

        // ---- 4. matching -----------------------------------------------------
        // auto: video -> sequential (frames are ordered); photo sets match
        // exhaustively up to ~400 images, then via vocabulary tree.
        int matcher = job.matcher;
        if (matcher == 0)
            matcher = job.is_video ? 2 : (n_images <= 400 ? 1 : 3);
        std::vector<std::string> ma;
        if (matcher == 3) {
            std::string vt = resolve_vocab_tree(job);
            if (vt.empty()) {
                if (job.matcher == 3)
                    return fail("no vocabulary tree available (see log)");
                log("warning: no vocabulary tree; falling back to exhaustive "
                    "matching (slow for large sets)");
                matcher = 1;
            } else {
                set_stage("Matching features (vocabulary tree)");
                ma = {job.colmap_exe, "vocab_tree_matcher",
                      "--database_path", db,
                      "--VocabTreeMatching.vocab_tree_path", vt};
            }
        }
        if (matcher == 2) {
            set_stage("Matching features (sequential)");
            ma = {job.colmap_exe, "sequential_matcher",
                  "--database_path", db,
                  "--SequentialMatching.overlap", std::to_string(job.seq_overlap),
                  "--SequentialMatching.quadratic_overlap",
                  job.seq_quadratic_overlap ? "1" : "0"};
        } else if (matcher == 1) {
            set_stage("Matching features (exhaustive)");
            ma = {job.colmap_exe, "exhaustive_matcher", "--database_path", db};
        }
        if (job.estimate_affine_shape) {
            ma.push_back("--FeatureMatching.guided_matching");
            ma.push_back("1");
        }
        rc = exec(ma);
        if (rc == kCancelled) return fail("cancelled");
        if (rc != 0) return fail("colmap matcher failed (see log)");

        // ---- 5. sparse reconstruction ----------------------------------------
        set_stage("Reconstructing cameras (mapper; this is the slow part)");
        fs::create_directories(ws / "sparse");
        rc = exec({job.colmap_exe, "mapper",
                   "--database_path", db,
                   "--image_path", images,
                   "--output_path", (ws / "sparse").string(),
                   "--Mapper.ba_use_gpu", job.ba_use_gpu ? "1" : "0",
                   "--Mapper.structure_less_registration_fallback", "0"});
        if (rc == kCancelled) return fail("cancelled");
        if (rc != 0) return fail("colmap mapper failed (see log)");
        // COLMAP may emit several models (sparse/0, sparse/1, ...); the
        // trainer auto-picks the one with the most registered images.
        int n_models = 0;
        {
            std::error_code ec;
            for (fs::directory_iterator it(ws / "sparse", ec), end;
                 !ec && it != end; it.increment(ec))
                if (fs::exists(it->path() / "cameras.bin")) n_models++;
        }
        if (n_models == 0)
            return fail("mapper produced no reconstruction (too few matches? "
                        "try higher quality or more overlapping images)");
        if (n_models > 1)
            log("Note: mapper produced " + std::to_string(n_models) +
                " partial models; the trainer uses the largest one");

        // ---- 6. bundle-adjustment refinement ---------------------------------
        if (job.final_bundle_adjust) {
            set_stage("Refining cameras (bundle adjustment)");
            rc = exec({job.colmap_exe, "bundle_adjuster",
                       "--input_path", (ws / "sparse" / "0").string(),
                       "--output_path", (ws / "sparse" / "0").string(),
                       "--BundleAdjustment.refine_focal_length", "1",
                       "--BundleAdjustment.refine_principal_point", "1",
                       "--BundleAdjustment.refine_extra_params", "1"});
            if (rc == kCancelled) return fail("cancelled");
            if (rc != 0) log("warning: bundle_adjuster failed; keeping the "
                             "mapper result");
        }

        // ---- 7. dataset marker ------------------------------------------------
        // Records the image dir so re-opening the dataset later needs no
        // manual configuration.
        {
            FILE* f = std::fopen((ws / "ssplat_dataset.json").string().c_str(), "w");
            if (f) {
                std::string esc;
                for (char c : image_dir_cfg) {
                    if (c == '"' || c == '\\') esc += '\\';
                    esc += c;
                }
                std::fprintf(f, "{\"image_dir\": \"%s\"}\n", esc.c_str());
                std::fclose(f);
            }
        }

        set_stage("Done");
        {
            std::lock_guard<std::mutex> lk(_mu);
            _dataset_dir = ws.string();
            _image_dir = image_dir_cfg;
        }
        _state = State::Done;
    } catch (const std::exception& e) {
        fail(e.what());
    }
}

}  // namespace gui
