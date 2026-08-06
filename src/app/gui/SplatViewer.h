#pragma once

// SplatViewer -- a splat file open for looking at, with no dataset and no
// training behind it.
//
// The engine renders from `world` splats and a camera; a training run is only
// one way to fill those in. This fills them from a file instead, so the same
// RenderWorker the training viewport uses (and the same one the web viewer
// uses) can show a finished model -- ours, or anyone else's 3D Gaussian
// Splatting PLY.
//
// Lifetime, and why it matters here more than usual: THE ENGINE IS A PROCESS-
// GLOBAL SINGLETON. Opening a file calls engine_reset() and takes it over, so
// a training run and an open file cannot coexist; GuiApp routes an open
// through the same "stop training?" confirmation as any other session-
// destroying action, and closes the viewer before setting a run up. Detach
// anything rendering from this (the viewport's RenderWorker) BEFORE close().
//
// Loading runs on a worker thread: a large model is a few hundred MB of PLY
// and the window has to stay responsive while it lands.

#include "app/webviewer/RenderWorker.h"
#include "data/DatasetParser.h"

#include <atomic>
#include <cstdint>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

namespace gui {

class SplatViewer {
public:
    enum class State { Idle, Loading, Ready, Failed };
    // What the file turned out to hold. A `.ply` is a container, not a
    // format: the same extension carries Gaussians, an SfM point cloud and a
    // mesh, and only the property list says which. Points still open -- in the
    // GL preview the trainer uses before a run -- because that is what a user
    // who dropped `points3D.ply` on the window meant.
    enum class Kind { Splats, Points };

    ~SplatViewer();

    // Accepts a splat .ply, a step-*.ckpt directory, or a run directory (whose
    // newest checkpoint is used). Asynchronous; watch state().
    void open(const std::string& path);
    // Give the engine back. Safe to call when nothing is open.
    void close();

    State state() const { return _state.load(); }
    bool ready() const { return _state.load() == State::Ready; }
    Kind kind() const { return _kind.load(); }

    std::string error();
    std::string path();          // what was asked for
    std::string file();          // the .ply actually opened
    int64_t num_splats() const { return _num_splats.load(); }
    int sh_degree() const { return _sh_degree.load(); }

    // Valid once ready(). The client navigates in a normalized frame centred
    // on the model and sized so it is about one unit across, exactly as it
    // does for a dataset.
    ViewerRenderConfig render_config();
    ViewerHooks make_hooks();
    // Identity of what is loaded, for the viewport's framing memory.
    std::string scene_key();

    // Valid once ready() with kind() == Points: a dataset with no cameras and
    // the cloud in it, which is all PreviewRenderer needs.
    const ParsedDataset& points() const { return _points; }
    const PostSplitCameras& post() const { return _post; }

    std::vector<std::string> drain_log();

private:
    void run(std::string path);
    void log(const std::string& s);

    std::thread _worker;
    std::atomic<State> _state{State::Idle};
    std::atomic<Kind> _kind{Kind::Splats};
    std::atomic<int64_t> _num_splats{0};
    std::atomic<int> _sh_degree{0};
    // True once this object has called engine_reset(): what close() has to
    // undo, and what says a half-finished load still left the engine dirty.
    std::atomic<bool> _owns_engine{false};

    std::mutex _mu;                    // guards everything below
    std::string _error, _path, _file;
    ViewerRenderConfig _cfg;
    std::vector<std::string> _log;
    // Points only. Written by the worker before it publishes Ready, read by
    // the GUI thread after; the state flag is the handoff.
    ParsedDataset _points;
    PostSplitCameras _post;

    // Handed to the RenderWorker as ViewerHooks::engine_mutex. Nothing else
    // touches the engine while a file is open, but the worker still takes it:
    // the hook is not optional, and it is what makes a future second reader
    // safe by construction.
    std::mutex _engine_mutex;
};

}  // namespace gui
