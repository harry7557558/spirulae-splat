#pragma once

// Viewer -- C++ port of the Python web viewer server (viewer/server.py +
// http_server.py) for the standalone CLI trainer. Serves the SAME
// single-file viewer.html (embedded at build time) over a dependency-free
// HTTP server, so the browser client is unchanged -- ssh port-forward the
// viewer port from a headless cloud box and open http://localhost:<port>/.
//
// Endpoints (unchanged from the Python server):
//   GET /              viewer.html
//   GET /render?...    JPEG of the requested buffer from the requested camera
//   GET /buffers       JSON list of viewable buffer keys
//   GET /progress      JSON training progress (step / eta / latency / paused)
//   GET /pause-toggle  toggle training pause, returns {"paused": bool}
//
// The actual render path (latest-wins worker thread + engine render +
// annotation blit) lives in RenderWorker.h -- shared with the native GUI
// viewport. This file adds only the HTTP layer + JPEG encode.

#include "data/DatasetParser.h"
#include "app/webviewer/RenderWorker.h"

#include <memory>
#include <string>

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
