#pragma once

// ViewportPanel -- the native 3D viewport, with two backends behind the
// browser-identical NavCamera navigation (see NavCamera.h: four modes,
// mouse / keyboard / touch-via-OS-pointer / gamepad input):
//
//   Preview: pure-GL sparse point cloud + camera frusta (PreviewRenderer),
//            shown as soon as a dataset is parsed, before training.
//   Engine:  the shared RenderWorker (same latest-wins engine render path
//            as the web viewer) once training has set up the engine; all
//            web-viewer camera models (pinhole / fisheye-equidistant /
//            fisheye-equisolid / equirectangular) with per-model FOV.
//
// GUI thread only. attach*() after the corresponding phase; detach() BEFORE
// the session it renders from is destroyed.

#include "../RenderWorker.h"
#include "NavCamera.h"
#include "PreviewRenderer.h"

#include <cstdint>
#include <string>
#include <vector>

struct ImVec2;
namespace ssplat { class TrainerSession; }

namespace gui {

class ViewportPanel {
public:
    // Dataset preview (needs load_dataset() done; no GPU engine).
    void attach_preview(ssplat::TrainerSession& session);
    // Engine renderer (needs engine_ready).
    void attach(ssplat::TrainerSession& session);
    void detach();
    bool attached() const { return _mode == Mode::Engine; }
    bool preview_active() const { return _mode == Mode::Preview; }

    // Draw the viewport (controls rows + image) into the current ImGui
    // window/child. `training` enables continuous refresh.
    void draw(bool training);

    // Free GL resources. Call while the GL context is still current
    // (before ImGui/GLFW shutdown).
    void destroy_gl();

private:
    enum class Mode { None, Preview, Engine };

    void compute_framing(const ssplat::TrainerSession& session);
    // Frame the scene only when a different dataset arrives; a preview ->
    // engine transition on the same dataset keeps the navigated pose and
    // intrinsics (no jump when training starts).
    void maybe_frame(const ssplat::TrainerSession& session);
    void reset_view();
    const char* camera_model_name() const;
    float fov_min() const;
    float fov_max() const;
    void compute_intrinsics(int W, int H, float& fx, float& fy) const;
    // Camera-to-orbit-target distance (normalized frame); drives the
    // zoom-adaptive grid cell size.
    float nav_dist() const;
    void build_request(ViewRequest& q, int W, int H) const;
    void view_matrix(float out[16]) const;   // row-major world-to-view
    void upload(const ViewResult& res);
    void draw_controls(bool engine);
    void handle_input(float item_h);         // mouse/keys/gamepad on last item
    void draw_engine(bool training, const ImVec2& avail);
    void draw_preview(const ImVec2& avail);

    Mode _mode = Mode::None;
    RenderWorker _worker;
    PreviewRenderer _preview;

    // Browser-identical camera + navigation.
    NavCamera _cam;
    NavCamera _home;                 // pose restored by Reset view
    float _home_dist = 3.0f;         // scene radius (near/far for preview)
    std::string _framed_key;         // dataset identity the pose belongs to
    bool _dragging = false;
    int _drag_button = 0;            // ImGuiMouseButton at drag start

    // Camera model + per-model FOV memory (browser _fovMemory equivalent).
    int _cam_model = 0;              // index into kViewerCameraModels
    float _fov_deg[4] = {90, 90, 90, 0};

    // Render options.
    int _buffer_idx = 0;
    std::vector<std::string> _buffer_keys;
    bool _show_cams = false;
    bool _show_grid = false;         // axes + ground-plane grid overlay
    float _frustum_scale = 1.0f;     // camera-frustum size multiplier
    int _scale_idx = 1;              // 0 = 50%, 1 = 75%, 2 = 100%
    bool _auto_refresh = true;

    // In-flight request / result texture (engine mode).
    uint64_t _pending = 0;
    double _last_submit = 0.0;
    bool _dirty = true;
    unsigned int _tex = 0;
    int _tex_w = 0, _tex_h = 0;
    std::string _last_error;
};

}  // namespace gui
