#pragma once

// TrainerCore -- the engine-side training session shared by the standalone
// CLI (main.cpp) and the native GUI (gui/). Extracted from main.cpp; the
// per-step plumbing remains a direct port of the Python managed-path trainer
// (see the mapping table in app/README.md):
//   model.py  _build_loss_weights / engine_train_step_managed
//   core.py   _build_optim_config / _build_densify_config
//   optimizer.py get_scheduled_lr
//   model.py  populate_modules (seeding) / _maybe_init_bilagrid /
//             background + color-space init
//   trainer.py train() save cadence / _setup_cpp_data_manager (non-warp)
//
// TrainerSession splits the run into phases so front-ends can interleave
// their own UI between them:
//   check_config()  -> throws for unported features, logs warnings
//   load_dataset()  -> parse + post-split bake (no GPU work)
//   setup_engine()  -> output dir, seeding, engine + DataManager init
//   train()         -> the step loop; pause / stop / render-fairness via
//                      the public atomics, per-step callback for progress
//
// The engine is a process-global singleton: one live session at a time.
// setup_engine() calls engine_reset(), so a fresh session can follow a
// finished one in the same process (the GUI's "train again" path).

#include "../Engine.h"
#include "DatasetParser.h"
#include "RenderWorker.h"
#include "generated/cli_config.h"

#include <array>
#include <atomic>
#include <chrono>
#include <deque>
#include <filesystem>
#include <functional>
#include <map>
#include <mutex>
#include <optional>
#include <string>
#include <vector>

namespace ssplat {

// ===========================================================================
// Color-space handling (port of _wrapper_per_pixel.get_color_transform_matrix
// + the resolution logic in model.py populate_modules:612-631 / __init__:532)
// ===========================================================================

using Mat3f = std::array<float, 9>;

Mat3f gamut_to_rec709(const std::string& name);
Mat3f invert3x3(const Mat3f& m);

struct ColorResolution {
    std::string splat_gamut;   // "" = Rec.709 / none
    std::string image_gamut;
    bool splat_linear  = false;
    bool image_linear  = false;
    bool convert_seed  = false;  // convert_initial_point_cloud_color resolved
};

ColorResolution resolve_color(const SsplatConfig& c);


// ===========================================================================
// LR schedule -- port of OptimizerConfig.get_scheduled_lr (optimizer.py:54).
// ===========================================================================

float scheduled_lr(int step, int max_steps, float lr,
                   std::optional<float> lr_final = std::nullopt,
                   std::optional<int> warmup = std::nullopt);


// ===========================================================================
// Splat seeding -- port of model.py populate_modules:546-651 (3dgs branch)
// ===========================================================================

struct SeedSplats {
    int64_t num = 0;                 // live splats (<= cap rows)
    std::vector<float> means;        // [cap, 3]
    std::vector<float> quats;        // [cap, 4]
    std::vector<float> scales;       // [cap, 3]  (log)
    std::vector<float> opacities;    // [cap, 1]  (logit)
    std::vector<float> features_dc;  // [cap, 3]
    std::vector<float> features_sh;  // [cap, dim_sh-1, 3]  (zeros)
};

SeedSplats seed_splats(const ColmapPoints3D& pts, const SsplatConfig& cfg,
                       const ColorResolution& color);


// ===========================================================================
// Per-step EngineStepConfig -- ports of model.py _build_loss_weights /
// engine_train_step_managed and core.py _build_optim_config /
// _build_densify_config.
// ===========================================================================

struct RunState {
    float train_frame_scale = 1.0f;
    bool  splat_linear = false;
    bool  bilagrid_rgb_init    = false;
    bool  bilagrid_depth_init  = false;
    bool  bilagrid_normal_init = false;
    bool  ppisp_init           = false;
};

std::array<float, (int)LossWeightIndex::length>
build_loss_weights(const SsplatConfig& c, int step);

EngineStepConfig build_step_config(const SsplatConfig& c, const RunState& st,
                                   int step);

// Nested-by-group config.json dump using the original Python field names.
void save_config_json(const SsplatConfig& c, const std::filesystem::path& out_dir,
                      const std::string& preset);


// ===========================================================================
// TrainerSession
// ===========================================================================

struct TrainerProgress {
    int step = 0;              // 0-based step that just finished
    int total_steps = 0;
    double step_latency = 0.0; // seconds, this step
    int64_t num_splats = 0;
    std::map<std::string, float> losses;
};

struct TrainerCallbacks {
    // Called after every completed step, engine mutex released.
    std::function<void(const TrainerProgress&)> on_step;
};

class TrainerSession {
public:
    TrainerSession() : _start_time(std::chrono::steady_clock::now()) {}
    TrainerSession(const TrainerSession&) = delete;
    TrainerSession& operator=(const TrainerSession&) = delete;

    // Inputs. Set before check_config().
    SsplatConfig cfg;
    std::string preset = "3dgs";
    // Human-readable progress/warning messages. Default (unset) = stdout.
    std::function<void(const std::string&)> log_fn;

    // Filled by load_dataset().
    ParsedDataset ds;
    PostSplitCameras post;
    bool has_depth = false;
    bool has_normal = false;

    // Filled by setup_engine().
    std::filesystem::path out_dir;
    RunState st;

    // Coordination between the train loop, viewer render workers, and
    // front-end controls (same dance as trainer.py:146-179).
    std::mutex engine_mutex;
    std::atomic<bool> paused{false};
    std::atomic<bool> stop_requested{false};
    std::atomic<bool> render_pending{false};
    std::atomic<int>  cur_step{0};

    // Throws std::runtime_error for features the managed C++ path does not
    // support; logs warnings for approximated ones.
    void check_config();

    // Parse the dataset + bake POST-split cameras. No GPU work.
    void load_dataset();

    // Create the output dir, dump config.json, reset + seed the engine,
    // set up the DataManager and bilagrid/PPISP. Requires load_dataset().
    void setup_engine();

    // The training loop. Returns when all steps ran or stop_requested was
    // set (a final checkpoint is saved either way unless steps_per_save==0).
    void train(const TrainerCallbacks& cb = {});

    void save_checkpoint(int step);

    // Trainer.get_progress port (trainer.py:181-205); the /progress body.
    std::string progress_json();

    // Viewer wiring shared by the web viewer and the GUI viewport.
    ViewerRenderConfig make_viewer_config() const;
    ViewerHooks make_viewer_hooks();

    void log(const std::string& msg);

private:
    std::chrono::steady_clock::time_point _start_time;
    std::mutex _progress_mutex;            // guards the latency window
    std::deque<double> _step_latencies;    // last 100, seconds
};

}  // namespace ssplat
