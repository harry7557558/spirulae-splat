// MeshRunner.cpp -- see MeshRunner.h.

#include "app/gui/MeshRunner.h"

#include "app/gui/AppPaths.h"
#include "app/gui/Subprocess.h"
#include "i18n/Locale.h"
#include "i18n/catalog/Log.h"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>

namespace fs = std::filesystem;
namespace lmsg = spirula::i18n::msg::log;
using spirula::i18n::format;

namespace gui {

namespace {

// The pipeline's stages, in the order the child prints them, and how far
// through the run each one is. The child's own `[meshing] <name> (...)` lines
// drive this -- no separate protocol to keep in sync.
// Each phase announces itself when it STARTS, so `at` is where the bar sits
// on entry and the next row's `at` is where it may creep to.
struct StageMark {
    const char* key;    // what the child prints after "[meshing] "
    float at;           // fraction complete when this line appears
};
const StageMark kStages[] = {
    {"loading", 0.00f},        {"point cloud", 0.04f},
    {"Delaunay", 0.10f},       {"occupancy field", 0.34f},
    {"cut edges", 0.60f},      {"bisection", 0.64f},
    {"marching tets", 0.74f},  {"merge", 0.78f},
    {"cleanup", 0.82f},        {"cull unseen", 0.85f},
    {"quality", 0.88f},        {"orient", 0.90f},
    {"color", 0.91f},          {"uv", 0.94f},
    {"bake", 0.96f},           {"texture", 0.98f},
    {"stats", 0.99f},          {"wrote", 0.995f},
};
constexpr int kNumStages = (int)(sizeof(kStages) / sizeof(kStages[0]));

// How long a phase is assumed to take when it reports nothing (seconds). The
// creep is 1 - exp(-t/tau), so it approaches the next mark without ever
// reaching or passing it -- a bar that moves but never lies about finishing.
constexpr double kCreepTau = 25.0;

// Whichever of the requested formats carries the most: a textured GLB is
// self-contained, a glTF needs its sidecars, OBJ has no vertex colors, PLY has
// no texture. Preview quality follows that order.
const char* kPreviewOrder[] = {"glb", "gltf", "obj", "ply"};

double now_seconds() {
    using namespace std::chrono;
    return duration<double>(steady_clock::now().time_since_epoch()).count();
}

}  // namespace

std::string MeshJob::preview_path() const {
    if (output.empty()) return {};
    for (const char* want : kPreviewOrder)
        for (int i = 0; i < kNumMeshFormats; ++i)
            if (formats[i] && std::strcmp(kMeshFormats[i], want) == 0)
                return output + "." + want;
    return {};
}

MeshRunner::~MeshRunner() {
    _cancel = true;
    if (_worker.joinable()) _worker.join();
}

void MeshRunner::start(const MeshJob& job) {
    if (_worker.joinable()) _worker.join();
    {
        std::lock_guard<std::mutex> lk(_mu);
        _error.clear();
        _stage.clear();
        _output.clear();
        _log.clear();
    }
    _verts = 0;
    _faces = 0;
    _stage_lo = -1.0f;
    _stage_hi = -1.0f;
    _stage_frac = -1.0f;
    _stage_at = now_seconds();
    _cancel = false;
    _run_id.fetch_add(1);
    _state = State::Running;
    _worker = std::thread([this, job] { run(job); });
}

void MeshRunner::cancel() { _cancel = true; }

std::string MeshRunner::stage() {
    std::lock_guard<std::mutex> lk(_mu);
    return _stage;
}
std::string MeshRunner::error() {
    std::lock_guard<std::mutex> lk(_mu);
    return _error;
}
std::string MeshRunner::output_path() {
    std::lock_guard<std::mutex> lk(_mu);
    return _output;
}

void MeshRunner::log(const std::string& line) {
    std::lock_guard<std::mutex> lk(_mu);
    _log.push_back(line);
    if (_log.size() > 4000) _log.erase(_log.begin(), _log.begin() + 1000);
}

std::vector<std::string> MeshRunner::drain_log() {
    std::lock_guard<std::mutex> lk(_mu);
    std::vector<std::string> out;
    out.swap(_log);
    return out;
}

float MeshRunner::progress() const {
    const float lo = _stage_lo.load();
    if (lo < 0.0f) return -1.0f;
#if 1
    return std::powf(lo, 2.0f);
#else
    const float hi = _stage_hi.load();
    const float sub = _stage_frac.load();
    if (sub >= 0.0f) return lo + (hi - lo) * std::min(std::max(sub, 0.0f), 1.0f);
    const double dt = std::max(0.0, now_seconds() - _stage_at.load());
    return lo + (hi - lo) * (float)(1.0 - std::exp(-dt / kCreepTau));
#endif
}

void MeshRunner::note_line(const std::string& line) {
    const std::string tag = "[meshing] ";
    const size_t at = line.find(tag);
    if (at == std::string::npos) return;
    const std::string rest = line.substr(at + tag.size());
    for (int i = 0; i < kNumStages; ++i) {
        if (rest.rfind(kStages[i].key, 0) != 0) continue;
        _stage_lo = kStages[i].at;
        _stage_hi = i + 1 < kNumStages ? kStages[i + 1].at : 1.0f;
        _stage_frac = -1.0f;
        _stage_at = now_seconds();
        std::lock_guard<std::mutex> lk(_mu);
        _stage = kStages[i].key;
        break;
    }
    // "<stage>: rendered 12/120 cameras" -- a phase reporting its own
    // fraction. Anchored on the " cameras" suffix rather than on "any a/b in
    // the line", so a path like /data/12/34/model.ply cannot be read as
    // progress. Monotone within a bracket: the bisection loop sweeps the
    // camera set once per iteration, and a bar that runs backwards reads as
    // a bug.
    const std::string kUnit = " cameras";
    const size_t unit = rest.rfind(kUnit);
    if (unit != std::string::npos && unit + kUnit.size() >= rest.size() - 1 &&
        _stage_lo.load() >= 0.0f) {
        size_t b = unit;
        while (b > 0 && std::isdigit((unsigned char)rest[b - 1])) --b;
        if (b > 0 && b < unit && rest[b - 1] == '/') {
            const size_t sl = b - 1;
            size_t a = sl;
            while (a > 0 && std::isdigit((unsigned char)rest[a - 1])) --a;
            if (a < sl) {
                const double done = std::atof(rest.c_str() + a);
                const double total = std::atof(rest.c_str() + b);
                if (total > 0.0) {
                    const float f = (float)(done / total);
                    if (f > _stage_frac.load()) _stage_frac = f;
                }
            }
        }
    }
    // "wrote <path>: N vertices, M faces" -- the file and its size, from the
    // one line that knows both.
    if (rest.rfind("wrote ", 0) == 0) {
        const size_t colon = rest.find(':', 6);
        if (colon != std::string::npos) {
            std::lock_guard<std::mutex> lk(_mu);
            if (_output.empty()) _output = rest.substr(6, colon - 6);
        }
    }
    long long v = 0, f = 0;
    if (std::sscanf(rest.c_str(), "done: %lld vertices, %lld faces", &v, &f) == 2) {
        _verts = v;
        _faces = f;
    }
}

void MeshRunner::run(MeshJob job) {
    auto fail = [&](const std::string& why) {
        std::lock_guard<std::mutex> lk(_mu);
        _error = why;
        _state = _cancel.load() ? State::Cancelled : State::Failed;
    };

    std::string formats;
    for (int i = 0; i < kNumMeshFormats; ++i)
        if (job.formats[i])
            formats += (formats.empty() ? "" : ",") + std::string(kMeshFormats[i]);
    if (formats.empty()) formats = "ply";

    if (job.output.empty()) {
        // Beside the checkpoint, as the CLI's own default does -- and for a
        // FILE source, `<name>_mesh` rather than `<name>`: the latter writes
        // `splat.ply` when meshing `splat.ply`. The child refuses that
        // outright (check_mesh_outputs_safe); this is the reason it never
        // has to.
        std::error_code ec;
        fs::path base(job.checkpoint);
        if (fs::is_regular_file(base, ec))
            base.replace_filename(base.stem().string() + "_mesh");
        else
            base = base / "mesh";
        job.output = base.string();
    }

    std::vector<std::string> argv = {
        // The child is this same executable, so it has the same thirteen
        // languages -- tell it which one, or its output lands in the log in
        // whatever the machine's locale is.
        exe_path(), "--lang",
        spirula::i18n::code(spirula::i18n::current()),
        "mesh", job.checkpoint,
        "--format", formats,
        "--color", kMeshColorModes[job.color],
        "--output", job.output,
    };
    if (!job.use_data) {
        argv.push_back("--no-data");
    } else if (!job.data_dir.empty()) {
        argv.push_back("--data");
        argv.push_back(job.data_dir);
    }
    if (job.max_cameras > 0) {
        argv.push_back("--max-cameras");
        argv.push_back(std::to_string(job.max_cameras));
    }
    if (job.texture_size > 0 && job.color == 2) {
        argv.push_back("--texture-size");
        argv.push_back(std::to_string(job.texture_size));
    }
    if (job.iso > 0.0f) {
        char buf[32];
        std::snprintf(buf, sizeof buf, "%g", job.iso);
        argv.push_back("--iso");
        argv.push_back(buf);
    }
    if (job.bisection_iters != 3) {
        argv.push_back("--bisection-iters");
        argv.push_back(std::to_string(job.bisection_iters));
    }
    if (job.merge_factor != 1.0f) {
        char buf[32];
        std::snprintf(buf, sizeof buf, "%g", job.merge_factor);
        argv.push_back("--merge-factor");
        argv.push_back(buf);
    }
    if (job.quality_iters != 3) {
        argv.push_back("--quality-iters");
        argv.push_back(std::to_string(job.quality_iters));
    }
    if (job.floater_min_faces != 10) {
        argv.push_back("--floater-min-faces");
        argv.push_back(std::to_string(job.floater_min_faces));
    }
    if (job.carve_k != 1) {
        argv.push_back("--carve-k");
        argv.push_back(std::to_string(job.carve_k));
    }
    if (!job.cull_unseen) argv.push_back("--no-cull-unseen");
    for (const std::string& a : split_args(job.extra_args)) argv.push_back(a);

    std::string cmd;
    for (const auto& a : argv) cmd += (cmd.empty() ? "$ " : " ") + a;
    log(cmd);

    const int rc = run_process(argv, "", [this](const std::string& l) {
        log(l);
        note_line(l);
    }, _cancel);

    if (rc == kCancelled) {
        _state = State::Cancelled;
        std::lock_guard<std::mutex> lk(_mu);
        _error = lmsg::err_cancelled.get();
        return;
    }
    if (rc == kSpawnFailed)
        return fail(format(lmsg::err_spawn_recon, {argv[0]}));
    if (rc != 0) return fail(lmsg::err_mesh_failed.get());

    {
        // The child writes the formats in checkbox order and its first
        // "wrote" line is whichever came first -- but the preview wants the
        // one that carries the most (a textured GLB over a bare PLY), so the
        // requested-format preference wins when that file is really there.
        std::error_code ec;
        const std::string best = job.preview_path();
        std::lock_guard<std::mutex> lk(_mu);
        if (!best.empty() && fs::is_regular_file(best, ec)) _output = best;
        else if (_output.empty()) _output = best;
    }
    _stage_lo = 1.0f;
    _stage_hi = 1.0f;
    _stage_frac = -1.0f;
    _state = State::Done;
}

}  // namespace gui
