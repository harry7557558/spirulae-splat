#pragma once

// Viewer -- C++ port of the Python web viewer server (viewer/server.py +
// http_server.py + render_worker.py + annotation.py) for the standalone CLI
// trainer. Serves the SAME single-file viewer.html (embedded at build time)
// over a dependency-free HTTP server, so the browser client is unchanged --
// ssh port-forward the viewer port from a headless cloud box and open
// http://localhost:<port>/.
//
// Endpoints (unchanged from the Python server):
//   GET /              viewer.html
//   GET /render?...    JPEG of the requested buffer from the requested camera
//   GET /buffers       JSON list of viewable buffer keys
//   GET /progress      JSON training progress (step / eta / latency / paused)
//   GET /pause-toggle  toggle training pause, returns {"paused": bool}
//
// Threading model (parity with the Python implementation):
//   - HTTP thread: parses requests; /render submits to the render worker
//     (latest-wins pending slot) and waits up to 10 s for the matching
//     result.
//   - Render worker thread: owns the actual engine render. Before touching
//     the engine it takes hooks.engine_mutex -- the same mutex the training
//     loop holds during each train step. hooks.set_render_pending(true) is
//     called synchronously at submit time so the training loop can yield
//     the mutex at the very next iteration boundary (the Python
//     _render_pending fairness dance, trainer.py:146-153).

#include "DatasetParser.h"

#include <array>
#include <atomic>
#include <functional>
#include <memory>
#include <mutex>
#include <optional>
#include <string>

// Static per-run info the render path needs (subset of model config).
struct ViewerRenderConfig {
    std::string primitive = "3dgs";
    bool  packed = false;
    int   sh_degree_warmup_every = 1000;
    std::optional<float> relative_scale;
    // Median-depth buffers are rendered only when a median loss is
    // configured (model.py get_outputs:1103-1108); requesting them
    // otherwise 400s, matching Python.
    bool  output_median = false;
    // Distortion render enabled when a distortion regularizer is configured;
    // also enabled on demand when a *_distortion buffer is requested.
    bool  distortion_reg_on = false;
    // Viewer-client c2w remap into the training frame (trainer.py:557-576).
    float train_frame_scale = 1.0f;
    std::array<float, 16> train_to_normalized{1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
};

struct ViewerHooks {
    // Serializes engine access between the training loop and viewer renders.
    std::mutex* engine_mutex = nullptr;
    // Current training step (for the SH-degree warmup schedule).
    std::function<int()> current_step;
    // JSON body for /progress (built by the training loop, which owns the
    // timing state).
    std::function<std::string()> progress_json;
    // Toggle pause; returns the new paused state.
    std::function<bool()> pause_toggle;
    // "Render desired" flag: called with true synchronously from the HTTP
    // thread at submit, with false from the worker once the queue drains.
    std::function<void(bool)> set_render_pending;
};

class ViewerServer {
public:
    ViewerServer();
    ~ViewerServer();

    // Uploads the POST-split camera table to the engine's viewer cache
    // (engine_viewer_init: frustum annotation + thumbnail slots; thumbnails
    // fill as training batches flow through), then starts the render worker
    // and HTTP threads. Call after engine_setup_data_manager.
    void start(const std::string& host, int port,
               ViewerRenderConfig cfg, ViewerHooks hooks,
               const PostSplitCameras& post);
    void stop();

private:
    struct Impl;
    std::unique_ptr<Impl> _impl;
};
