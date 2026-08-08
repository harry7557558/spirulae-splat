#pragma once

// SplatViewer -- a FILE open for looking at, with no dataset and no training
// behind it: a 3D Gaussian Splatting model, an SfM point cloud, or an
// extracted triangle mesh. What it turns out to be decides which renderer
// the viewport uses; opening it is the same action either way, which is what
// makes "drop a file on the window" work.
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
#include "mesh/MeshExport.h"   // meshing::MeshData

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
    // mesh, and only the property list says which. Points and meshes open in
    // the GL preview renderer the trainer uses before a run; only Splats need
    // the engine (and so only Splats take it over).
    enum class Kind { Splats, Points, Mesh };

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

    // Valid once ready() with kind() == Mesh: the triangles, and the
    // similarity that maps them into the normalized frame the viewport
    // navigates (row-major 3x4, the layout PreviewRenderer takes).
    const meshing::MeshData& mesh() const { return _mesh; }
    const float* mesh_to_normalized() const { return _mesh_t2n; }
    int64_t num_faces() const { return _num_faces.load(); }

    // Re-run the engine's linear / wide-gamut setup for the OPEN model.
    // The viewer offers it because a PLY records no color space: a model
    // trained in ACEScg looks wrong until it is told so, and the only way to
    // find out is to try. Takes the engine lock; call from the GUI thread.
    // `gamut` is a gamut_to_rec709 name, "" for Rec.709.
    void set_color_space(const char* gamut, bool linear);
    // What the run's config.json said, so the controls open on the right
    // setting rather than on "none".
    std::string gamut() const;
    bool linear_color() const;

    // Hand the previous primitive's screen buffers back after a primitive
    // switch (see engine_release_screen_buffers). Takes the engine lock;
    // call from the GUI thread, between renders.
    void release_screen_buffers();

    std::vector<std::string> drain_log();

private:
    void run(std::string path);
    void log(const std::string& s);

    std::thread _worker;
    std::atomic<State> _state{State::Idle};
    std::atomic<Kind> _kind{Kind::Splats};
    std::atomic<int64_t> _num_splats{0};
    std::atomic<int> _sh_degree{0};
    std::atomic<int64_t> _num_faces{0};
    // True once this object has called engine_reset(): what close() has to
    // undo, and what says a half-finished load still left the engine dirty.
    std::atomic<bool> _owns_engine{false};

    std::mutex _mu;                    // guards everything below
    std::string _error, _path, _file;
    std::string _gamut;                // "" = Rec.709
    bool _linear = false;
    ViewerRenderConfig _cfg;
    std::vector<std::string> _log;
    // Points only. Written by the worker before it publishes Ready, read by
    // the GUI thread after; the state flag is the handoff.
    ParsedDataset _points;
    PostSplitCameras _post;
    // Mesh only, same handoff.
    meshing::MeshData _mesh;
    float _mesh_t2n[12] = {1,0,0,0, 0,1,0,0, 0,0,1,0};

    // Handed to the RenderWorker as ViewerHooks::engine_mutex. Nothing else
    // touches the engine while a file is open, but the worker still takes it:
    // the hook is not optional, and it is what makes a future second reader
    // safe by construction.
    std::mutex _engine_mutex;
};

}  // namespace gui
