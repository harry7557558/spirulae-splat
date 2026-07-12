// ColmapRunner.cpp -- see ColmapRunner.h.

#include "ColmapRunner.h"
#include "Subprocess.h"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <filesystem>

namespace fs = std::filesystem;

namespace gui {

namespace {

bool is_image_file(const fs::path& p) {
    std::string e = p.extension().string();
    for (auto& c : e) c = (char)std::tolower((unsigned char)c);
    return e == ".jpg" || e == ".jpeg" || e == ".png" || e == ".webp" ||
           e == ".tif" || e == ".tiff" || e == ".bmp";
}

int count_images(const std::string& dir) {
    int n = 0;
    std::error_code ec;
    for (fs::directory_iterator it(dir, ec), end; !ec && it != end; it.increment(ec))
        if (it->is_regular_file(ec) && is_image_file(it->path())) n++;
    return n;
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

void ColmapRunner::run(ColmapJob job) {
    auto fail = [&](const std::string& why) {
        std::lock_guard<std::mutex> lk(_mu);
        _error = why;
        _state = _cancel.load() ? State::Cancelled : State::Failed;
    };
    try {
        const fs::path ws = job.workspace;
        fs::create_directories(ws);

        // ---- 1. resolve the image directory --------------------------------
        std::string images;
        std::string image_dir_cfg;   // what the trainer's image_dir should be
        if (job.is_video) {
            set_stage("Extracting frames (ffmpeg)");
            if (!command_exists(job.ffmpeg_exe))
                return fail("ffmpeg not found ('" + job.ffmpeg_exe +
                            "'); install it or set its path in the form above");
            fs::create_directories(ws / "images");
            char vf[64];
            std::snprintf(vf, sizeof vf, "fps=%g", (double)job.video_fps);
            int rc = exec({job.ffmpeg_exe, "-y", "-i", job.input_path,
                           "-vf", vf, "-qscale:v", "2",
                           (ws / "images" / "frame_%05d.jpg").string()});
            if (rc == kCancelled) return fail("cancelled");
            if (rc != 0) return fail("ffmpeg frame extraction failed (see log)");
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

        if (!command_exists(job.colmap_exe))
            return fail("colmap not found ('" + job.colmap_exe +
                        "'); install COLMAP or set its path in the form above");

        // Quality knobs (COLMAP's own low/medium/high presets, trimmed).
        const char* max_image_size  = job.quality == 0 ? "1600"
                                    : job.quality == 1 ? "2400" : "3200";
        const char* max_num_features = job.quality == 0 ? "4096"
                                     : job.quality == 1 ? "8192" : "16384";
        const std::string db = (ws / "database.db").string();
        const char* gpu = job.use_gpu ? "1" : "0";

        // ---- 2. feature extraction -----------------------------------------
        set_stage("Extracting features (colmap)");
        int rc = exec({job.colmap_exe, "feature_extractor",
                       "--database_path", db,
                       "--image_path", images,
                       "--ImageReader.camera_model", job.camera_model,
                       "--ImageReader.single_camera", job.single_camera ? "1" : "0",
                       "--SiftExtraction.use_gpu", gpu,
                       "--SiftExtraction.max_image_size", max_image_size,
                       "--SiftExtraction.max_num_features", max_num_features});
        if (rc == kCancelled) return fail("cancelled");
        if (rc != 0) return fail("colmap feature_extractor failed (see log)");

        // ---- 3. matching -----------------------------------------------------
        // auto: video or large image sets match sequentially (frames are
        // ordered), photo sets match exhaustively.
        bool sequential = job.matcher == 2 ||
                          (job.matcher == 0 && (job.is_video || n_images > 400));
        if (sequential) {
            set_stage("Matching features (sequential)");
            rc = exec({job.colmap_exe, "sequential_matcher",
                       "--database_path", db,
                       "--SequentialMatching.overlap", "10",
                       "--SiftMatching.use_gpu", gpu});
        } else {
            set_stage("Matching features (exhaustive)");
            rc = exec({job.colmap_exe, "exhaustive_matcher",
                       "--database_path", db,
                       "--SiftMatching.use_gpu", gpu});
        }
        if (rc == kCancelled) return fail("cancelled");
        if (rc != 0) return fail("colmap matcher failed (see log)");

        // ---- 4. sparse reconstruction ----------------------------------------
        set_stage("Reconstructing cameras (mapper; this is the slow part)");
        fs::create_directories(ws / "sparse");
        rc = exec({job.colmap_exe, "mapper",
                   "--database_path", db,
                   "--image_path", images,
                   "--output_path", (ws / "sparse").string()});
        if (rc == kCancelled) return fail("cancelled");
        if (rc != 0) return fail("colmap mapper failed (see log)");
        if (!fs::exists(ws / "sparse" / "0" / "cameras.bin"))
            return fail("mapper produced no reconstruction (too few matches? "
                        "try higher quality or more overlapping images)");

        // ---- 5. dataset marker ------------------------------------------------
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
