#pragma once

// ViewportPanel -- the native 3D viewport. Drives the shared RenderWorker
// (the same latest-wins engine render path the web viewer uses) and shows
// the result in an OpenGL texture inside ImGui, with turntable navigation
// matching the web client (Z-up normalized frame):
//   left-drag orbit, right/middle-drag pan, wheel zoom.
//
// GUI thread only. Attach after the engine is set up (TrainRunner
// engine_ready()); detach BEFORE the session it renders from is destroyed.

#include "../RenderWorker.h"

#include <cstdint>
#include <string>
#include <vector>

namespace ssplat { class TrainerSession; }

namespace gui {

class ViewportPanel {
public:
    void attach(ssplat::TrainerSession& session);
    void detach();
    bool attached() const { return _attached; }

    // Draw the viewport (controls row + image) into the current ImGui
    // window/child. `training` enables continuous refresh.
    void draw(bool training);

    // Free the GL texture. Call while the GL context is still current
    // (before ImGui/GLFW shutdown).
    void destroy_gl();

private:
    void reset_view();
    void build_request(ViewRequest& q, int W, int H) const;
    void upload(const ViewResult& res);

    RenderWorker _worker;
    bool _attached = false;

    // Scene framing (normalized frame, Z-up).
    float _target[3] = {0, 0, 0};
    float _home_target[3] = {0, 0, 0};
    float _yaw = 0.0f, _pitch = 0.6f;
    float _dist = 3.0f, _home_dist = 3.0f;
    float _fov_deg = 60.0f;

    // Render options.
    int _buffer_idx = 0;
    std::vector<std::string> _buffer_keys;
    bool _show_cams = false;
    int _scale_idx = 1;              // 0 = 50%, 1 = 75%, 2 = 100%
    bool _auto_refresh = true;

    // In-flight request / result texture.
    uint64_t _pending = 0;
    double _last_submit = 0.0;
    bool _dirty = true;
    unsigned int _tex = 0;
    int _tex_w = 0, _tex_h = 0;
    std::string _last_error;
};

}  // namespace gui
