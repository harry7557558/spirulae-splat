// ColmapRunner.cpp -- see ColmapRunner.h. CLI flags mirror
// scripts/run_colmap.bash (COLMAP >= 4.x; use_gpu-style flags are gone).

#include "app/gui/ColmapRunner.h"
#include "app/gui/FrameSelect.h"
#include "app/gui/Subprocess.h"

#include "app_generated/mask_py.h"   // kMaskPy[] (CMake-embedded scripts/mask.py)

#include "external/stb_image.h"   // stbi_info (image size probe)

#ifndef _WIN32
#include <ftw.h>
#endif

#include <algorithm>
#include <cctype>
#include <cmath>
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

bool is_fisheye_model(const std::string& m) {
    return m.find("FISHEYE") != std::string::npos;
}

// First image found under dir (recursive), for probing dimensions.
bool first_image_dims(const std::string& dir, int& W, int& H) {
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

// Initial ImageReader.camera_params for a model: given focal length and a
// centered principal point, all distortion coefficients zero. Order follows
// COLMAP's camera model definitions.
std::string compose_camera_params(const std::string& model, double f,
                                  double cx, double cy) {
    struct M { const char* name; int n_focal; int n_extra; };
    static const M kModels[] = {
        {"SIMPLE_PINHOLE", 1, 0},        {"PINHOLE", 2, 0},
        {"SIMPLE_RADIAL", 1, 1},         {"RADIAL", 1, 2},
        {"OPENCV", 2, 4},                {"FULL_OPENCV", 2, 8},
        {"OPENCV_FISHEYE", 2, 4},        {"THIN_PRISM_FISHEYE", 2, 8},
        {"SIMPLE_RADIAL_FISHEYE", 1, 1}, {"RADIAL_FISHEYE", 1, 2},
    };
    for (const M& m : kModels) {
        if (model != m.name) continue;
        char buf[64];
        std::string out;
        for (int i = 0; i < m.n_focal; i++) {
            std::snprintf(buf, sizeof buf, "%.4f,", f);
            out += buf;
        }
        std::snprintf(buf, sizeof buf, "%.4f,%.4f", cx, cy);
        out += buf;
        for (int i = 0; i < m.n_extra; i++) out += ",0";
        return out;
    }
    return "";
}

// Registered-image count of a COLMAP model dir (uint64 head of images.bin;
// same trick as ColmapParser's largest-model pick).
int64_t model_num_images(const fs::path& dir) {
    FILE* f = std::fopen((dir / "images.bin").string().c_str(), "rb");
    if (!f) return 0;
    uint64_t n = 0;
    size_t got = std::fread(&n, sizeof n, 1, f);
    std::fclose(f);
    return got == 1 ? (int64_t)n : 0;
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

double ColmapRunner::model_reproj_error(const ColmapJob& job,
                                        const std::string& model) {
    double err = 1e30;
    run_process({job.colmap_exe, "model_analyzer", "--path", model}, "",
                [&](const std::string& l) {
                    const char* key = "Mean reprojection error:";
                    size_t p = l.find(key);
                    if (p == std::string::npos) return;
                    try { err = std::stod(l.substr(p + std::strlen(key))); }
                    catch (...) {}
                },
                _cancel);
    return err;
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

        // Interrupted-run handling: either resume (reuse everything the
        // previous run completed) or insist on a clean folder -- never
        // silently mix a fresh run into stale artifacts.
        {
            std::error_code ec;
            bool prior = fs::exists(ws / "database.db", ec) ||
                         (fs::is_directory(ws / "images", ec) &&
                          !fs::is_empty(ws / "images", ec)) ||
                         fs::is_directory(ws / "sparse", ec);
            if (prior && !job.resume)
                return fail("the workspace already contains a previous run "
                            "(database.db / images / sparse); enable "
                            "\"Resume previous run\" to reuse it, or choose "
                            "an empty folder");
            if (prior)
                log("Resuming previous run in " + ws.string() +
                    " (completed stages are reused)");
        }

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
                fs::path out_dir = streams.size() > 1
                    ? ws / "images" / ("cam" + std::to_string(tr))
                    : ws / "images";
                // Resume: a non-empty output folder means this track's
                // extraction completed (frames are moved in one batch after
                // selection, and a cancelled selection leaves it empty).
                if (job.resume) {
                    int have = count_images(out_dir.string());
                    if (have > 0) {
                        log("Resume: keeping " + std::to_string(have) +
                            " extracted frames in " + out_dir.string() +
                            " (delete the folder to re-extract)");
                        continue;
                    }
                }
                if (streams.size() > 1) {
                    set_stage("Splitting video track " + std::to_string(tr));
                    int rc = exec({job.ffmpeg_exe, "-nostdin", "-y",
                                   "-i", job.input_path,
                                   "-map", "0:v:" + std::to_string(tr),
                                   "-c", "copy",
                                   (ws / ("track_cam" + std::to_string(tr) + ".mp4")).string()});
                    if (rc == kCancelled) return fail("cancelled");
                    if (rc != 0) return fail("ffmpeg track split failed (see log)");
                    track_path = (ws / ("track_cam" + std::to_string(tr) + ".mp4")).string();
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
        // matching; O(n^2) in feature count). ALIKED extracts far fewer,
        // higher-quality keypoints than SIFT.
        bool aliked = job.feature_type == 1;
        int features = job.max_num_features > 0 ? job.max_num_features
                     : aliked ? (job.quality == 0 ? 1024
                               : job.quality == 1 ? 2048 : 4096)
                     : (job.quality == 0 ? 4096
                      : job.quality == 1 ? 8192 : 16384);
        const std::string db = (ws / "database.db").string();

        // Initial camera parameters: an explicit ImageReader.camera_params
        // string wins; otherwise compose one from the focal-length factor
        // (fx = fy = factor * width, centered principal point, zero
        // distortion). A good initial focal length stabilizes mapper
        // initialization a lot, especially for fisheye lenses.
        std::string cam_params = job.camera_params;
        if (cam_params.empty() && job.init_focal_factor > 0) {
            int W = 0, H = 0;
            if (first_image_dims(images, W, H)) {
                cam_params = compose_camera_params(
                    job.camera_model, (double)job.init_focal_factor * W,
                    0.5 * W, 0.5 * H);
                if (cam_params.empty())
                    log("warning: no camera_params template for " +
                        job.camera_model + "; skipping the initial focal length");
                else
                    log("Initial camera (" + job.camera_model + ", " +
                        std::to_string(W) + "x" + std::to_string(H) +
                        "): " + cam_params);
            } else {
                log("warning: could not read an image size; skipping the "
                    "initial focal length");
            }
        }

        // ---- 3. feature extraction -----------------------------------------
        set_stage(aliked ? "Extracting features (colmap, ALIKED)"
                         : "Extracting features (colmap)");
        std::vector<std::string> fe = {job.colmap_exe, "feature_extractor",
            "--database_path", db,
            "--image_path", images,
            "--ImageReader.camera_model", job.camera_model};
        if (aliked) {
            fe.push_back("--FeatureExtraction.type");
            fe.push_back("ALIKED");
            fe.push_back("--AlikedExtraction.max_num_features");
            fe.push_back(std::to_string(features));
        } else {
            fe.push_back("--SiftExtraction.max_num_features");
            fe.push_back(std::to_string(features));
        }
        if (!cam_params.empty()) {
            fe.push_back("--ImageReader.camera_params");
            fe.push_back(cam_params);
        }
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
        if (job.estimate_affine_shape && !aliked) {
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
        // The matcher is an explicit choice (the GUI presets sequential for
        // video / exhaustive for photos; stale configs fall back the same
        // way). The vocabulary tree (matcher + sequential loop detection)
        // indexes SIFT descriptors only.
        int matcher = job.matcher;
        if (matcher < 1 || matcher > 3)
            matcher = job.is_video ? 2 : (n_images <= 400 ? 1 : 3);
        std::string match_type;   // FeatureMatching.type ("" = COLMAP default)
        if (aliked)
            match_type = job.lightglue ? "ALIKED_LIGHTGLUE" : "ALIKED_BRUTEFORCE";
        else if (job.lightglue)
            match_type = "SIFT_LIGHTGLUE";
        std::vector<std::string> ma;
        if (matcher == 3) {
            if (aliked)
                return fail("vocabulary-tree matching requires SIFT features "
                            "(the tree indexes SIFT descriptors); pick the "
                            "sequential or exhaustive matcher for ALIKED");
            std::string vt = resolve_vocab_tree(job);
            if (vt.empty())
                return fail("no vocabulary tree available (see log)");
            set_stage("Matching features (vocabulary tree)");
            ma = {job.colmap_exe, "vocab_tree_matcher",
                  "--database_path", db,
                  "--VocabTreeMatching.vocab_tree_path", vt};
        }
        if (matcher == 2) {
            std::string vt;
            if (job.seq_loop_closure) {
                if (aliked)
                    log("note: loop detection needs SIFT descriptors "
                        "(vocabulary tree); relying on quadratic overlap "
                        "for loop closure");
                else if ((vt = resolve_vocab_tree(job)).empty())
                    log("warning: no vocabulary tree; loop detection disabled");
            }
            set_stage("Matching features (sequential)");
            ma = {job.colmap_exe, "sequential_matcher",
                  "--database_path", db,
                  "--SequentialMatching.overlap", std::to_string(job.seq_overlap),
                  "--SequentialMatching.quadratic_overlap",
                  job.seq_quadratic_overlap ? "1" : "0"};
            if (!vt.empty()) {
                ma.push_back("--SequentialMatching.loop_detection");
                ma.push_back("1");
                ma.push_back("--SequentialMatching.vocab_tree_path");
                ma.push_back(vt);
            }
        } else if (matcher == 1) {
            set_stage("Matching features (exhaustive)");
            ma = {job.colmap_exe, "exhaustive_matcher", "--database_path", db};
        }
        if (!match_type.empty()) {
            ma.push_back("--FeatureMatching.type");
            ma.push_back(match_type);
        }
        if (job.estimate_affine_shape && !aliked) {
            ma.push_back("--FeatureMatching.guided_matching");
            ma.push_back("1");
        }
        // Wrong-match suppression (repetitive scenes).
        if (job.match_max_ratio > 0 && !aliked) {
            char b[16];
            std::snprintf(b, sizeof b, "%g", job.match_max_ratio);
            ma.push_back("--SiftMatching.max_ratio");
            ma.push_back(b);
        }
        if (job.min_inliers_per_pair > 0) {
            ma.push_back("--TwoViewGeometry.min_num_inliers");
            ma.push_back(std::to_string(job.min_inliers_per_pair));
        }
        rc = exec(ma);
        if (rc == kCancelled) return fail("cancelled");
        if (rc != 0) return fail("colmap matcher failed (see log)");

        // ---- 5. sparse reconstruction ----------------------------------------
        set_stage("Reconstructing cameras (mapper; this is the slow part)");
        fs::create_directories(ws / "sparse");
        bool fisheye = is_fisheye_model(job.camera_model);
        bool ba_gpu = job.ba_use_gpu;
        if (ba_gpu && fisheye) {
            ba_gpu = false;
            log("note: COLMAP's GPU bundle adjustment does not support "
                "fisheye camera models; using the CPU backend");
        }
        // Low-distortion perspective models map more stably with the
        // distortion coefficients FIXED at their initial values; the final
        // refinement pass then recovers them (run_colmap.bash's advice).
        // Fisheye models need their distortion refined during mapping.
        bool no_extra = job.camera_model == "PINHOLE" ||
                        job.camera_model == "SIMPLE_PINHOLE";
        bool fix_extra = !no_extra &&
                         (job.mapper_extra_params == 2 ||
                          (job.mapper_extra_params == 0 && !fisheye));
        std::vector<std::string> mp = {job.colmap_exe, "mapper",
                   "--database_path", db,
                   "--image_path", images,
                   "--output_path", (ws / "sparse").string(),
                   "--Mapper.ba_use_gpu", ba_gpu ? "1" : "0",
                   "--Mapper.structure_less_registration_fallback", "0"};
        if (fix_extra) {
            mp.push_back("--Mapper.ba_refine_extra_params");
            mp.push_back("0");
            log(job.final_bundle_adjust
                ? "Distortion held fixed during mapping; the final "
                  "refinement pass recovers it"
                : "warning: distortion held fixed during mapping but the "
                  "final refinement pass is disabled -- coefficients stay "
                  "at their initial values");
        }
        if (job.min_num_matches > 0) {
            mp.push_back("--Mapper.min_num_matches");
            mp.push_back(std::to_string(job.min_num_matches));
        }
        // Stricter image registration (repetitive scenes).
        char fbuf[16];
        if (job.abs_pose_min_num_inliers > 0) {
            mp.push_back("--Mapper.abs_pose_min_num_inliers");
            mp.push_back(std::to_string(job.abs_pose_min_num_inliers));
        }
        if (job.abs_pose_min_inlier_ratio > 0) {
            std::snprintf(fbuf, sizeof fbuf, "%g", job.abs_pose_min_inlier_ratio);
            mp.push_back("--Mapper.abs_pose_min_inlier_ratio");
            mp.push_back(fbuf);
        }
        if (job.abs_pose_max_error > 0) {
            std::snprintf(fbuf, sizeof fbuf, "%g", job.abs_pose_max_error);
            mp.push_back("--Mapper.abs_pose_max_error");
            mp.push_back(fbuf);
        }
        // COLMAP may emit several models (sparse/0, sparse/1, ...); the
        // trainer auto-picks the one with the most registered images.
        auto enumerate_models = [&]() {
            std::vector<std::pair<int64_t, fs::path>> ms;   // (-count, dir)
            std::error_code ec;
            for (fs::directory_iterator it(ws / "sparse", ec), end;
                 !ec && it != end; it.increment(ec))
                if (it->path().filename().string().rfind(".", 0) != 0 &&
                    fs::exists(it->path() / "cameras.bin"))
                    ms.push_back({-model_num_images(it->path()), it->path()});
            std::sort(ms.begin(), ms.end());
            return ms;
        };
        // Resume: the mapper only writes models on completion, so any
        // existing model is from a FINISHED mapper run -- reuse it. (An
        // interrupted mapper leaves nothing and simply reruns; features and
        // matches were reused from the database either way.)
        std::vector<std::pair<int64_t, fs::path>> models;
        if (job.resume && !(models = enumerate_models()).empty()) {
            log("Resume: " + std::to_string(models.size()) +
                " existing model(s) under sparse/; skipping the mapper "
                "(delete sparse/ to re-reconstruct)");
        } else {
            rc = exec(mp);
            if (rc == kCancelled) return fail("cancelled");
            if (rc != 0) return fail("colmap mapper failed (see log)");
            models = enumerate_models();
        }
        if (models.empty())
            return fail("mapper produced no reconstruction (too few matches? "
                        "try higher quality or more overlapping images)");
        fs::path best = models[0].second;
        if (models.size() > 1) {
            std::string counts;
            for (auto& m : models)
                counts += (counts.empty() ? "" : ", ") + std::to_string(-m.first);
            log("Note: mapper produced " + std::to_string(models.size()) +
                " partial models (" + counts + " images)");
        }

        // ---- 6. merge partial models (best effort) ---------------------------
        // model_merger only succeeds when two models share registered
        // frames, which mapper partials usually don't -- but when it works
        // it recovers a broken-apart reconstruction, so it's always worth
        // the few seconds. A merge is kept only when the merged model
        // registers more images than the base.
        if (models.size() > 1 && job.merge_models) {
            set_stage("Merging partial models");
            fs::path cur = best;
            int64_t cur_n = -models[0].first;
            fs::path acc = ws / "sparse" / ".merge_acc";
            fs::path tmp = ws / "sparse" / ".merge_tmp";
            int merged = 0;
            for (size_t i = 1; i < models.size(); i++) {
                remove_tree(tmp);
                fs::create_directories(tmp);
                rc = exec({job.colmap_exe, "model_merger",
                           "--input_path1", cur.string(),
                           "--input_path2", models[i].second.string(),
                           "--output_path", tmp.string()});
                if (rc == kCancelled) return fail("cancelled");
                int64_t n = rc == 0 ? model_num_images(tmp) : 0;
                if (n > cur_n) {
                    remove_tree(acc);
                    std::error_code ec;
                    fs::rename(tmp, acc, ec);
                    if (ec) break;
                    cur = acc;
                    cur_n = n;
                    merged++;
                    log("Merged " + models[i].second.filename().string() +
                        " -> " + std::to_string(n) + " images");
                } else {
                    log("Could not merge " +
                        models[i].second.filename().string() +
                        " (models share no registered frames); skipped");
                }
            }
            remove_tree(tmp);
            if (merged > 0) {
                // Persist as the next sparse/<N> so the trainer's
                // largest-model auto-pick finds it.
                int next = 0;
                while (fs::exists(ws / "sparse" / std::to_string(next))) next++;
                fs::path dst = ws / "sparse" / std::to_string(next);
                std::error_code ec;
                fs::rename(acc, dst, ec);
                if (!ec) {
                    best = dst;
                    log("Merged model written to sparse/" +
                        std::to_string(next) + " (" + std::to_string(cur_n) +
                        " images); the trainer auto-picks it");
                }
            }
        }

        // ---- 7. bundle-adjustment refinement ---------------------------------
        // Runs on the LARGEST (or merged) model, releasing what the mapper
        // held fixed: focal length, distortion, and (perspective models
        // only) the principal point -- releasing pp together with the 8
        // thin-prism coefficients on a ~200-degree fisheye can diverge to
        // astronomic reprojection errors. The refined model is checked
        // against the input (mean reprojection error via model_analyzer)
        // and REVERTED when it got worse or non-finite.
        if (job.final_bundle_adjust) {
            set_stage("Refining cameras (bundle adjustment)");
            fs::path tmp = ws / "sparse" / ".ba_tmp";
            remove_tree(tmp);
            fs::create_directories(tmp);
            rc = exec({job.colmap_exe, "bundle_adjuster",
                       "--input_path", best.string(),
                       "--output_path", tmp.string(),
                       "--BundleAdjustment.refine_focal_length", "1",
                       "--BundleAdjustment.refine_principal_point",
                       fisheye ? "0" : "1",
                       "--BundleAdjustment.refine_extra_params", "1"});
            if (rc == kCancelled) return fail("cancelled");
            double before = model_reproj_error(job, best.string());
            double after = rc == 0 ? model_reproj_error(job, tmp.string()) : 1e30;
            char msg[160];
            std::snprintf(msg, sizeof msg,
                          "Mean reprojection error: %.4g px -> %.4g px",
                          before, after);
            log(msg);
            if (rc == 0 && std::isfinite(after) &&
                (after <= before || !std::isfinite(before))) {
                std::error_code ec;
                for (const char* f : {"cameras.bin", "images.bin",
                                      "points3D.bin", "frames.bin",
                                      "rigs.bin"}) {
                    if (!fs::exists(tmp / f, ec)) continue;
                    fs::rename(tmp / f, best / f, ec);
                }
                log("Refinement kept");
            } else {
                log("warning: bundle_adjuster " +
                    std::string(rc != 0 ? "failed" : "made the model worse") +
                    "; keeping the mapper result");
            }
            remove_tree(tmp);
        }

        // No marker file is written: the image dir is handed to the GUI
        // in-memory for the immediate open; on later re-opens the parser
        // default applies (video datasets use images/ anyway) and photo-in-
        // place datasets need data.image_dir set in the dataparser options.
        if (!job.is_video)
            log("Note: images are referenced in place; when re-opening this "
                "dataset later, set image_dir to " + image_dir_cfg +
                " under the dataset-parsing options");

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
