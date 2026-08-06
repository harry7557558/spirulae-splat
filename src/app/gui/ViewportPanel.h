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

#include "app/webviewer/RenderWorker.h"
#include "app/gui/NavCamera.h"
#include "app/gui/PreviewRenderer.h"

#include <cstdint>
#include <string>
#include <vector>

struct ImVec2;
namespace spirula { class TrainerSession; }

namespace gui {

class ViewportPanel {
public:
    // Dataset preview (needs load_dataset() done; no GPU engine).
    void attach_preview(spirula::TrainerSession& session);
    // The same GL preview over a bare point cloud -- a .ply that turned out to
    // hold points rather than Gaussians. No cameras, so no frusta.
    void attach_preview_data(const ParsedDataset& ds, const PostSplitCameras& post,
                             const std::string& key, float radius = 1.0f);
    // Engine renderer (needs engine_ready).
    void attach(spirula::TrainerSession& session);
    // Engine renderer over something that is not a training session -- a splat
    // file opened in the viewer (SplatViewer). `key` identifies the scene, so
    // reopening the same file keeps the pose; `radius` is the scene radius in
    // the client's normalized frame, which is 1 for anything normalized.
    // There are no training cameras behind this, so the frustum controls are
    // not offered.
    void attach_scene(const ViewerRenderConfig& cfg, const ViewerHooks& hooks,
                      const std::string& key, float radius = 1.0f);
    void detach();
    bool attached() const { return _mode == Mode::Engine; }
    bool preview_active() const { return _mode == Mode::Preview; }

    // Draw the viewport (controls rows + image) into the current ImGui
    // window/child. `training` enables continuous refresh; `step` is the
    // training iteration that has just finished (< 0 when unknown), which is
    // what paces the refresh while nobody is steering the camera.
    void draw(bool training, int step = -1);

    // Free GL resources. Call while the GL context is still current
    // (before ImGui/GLFW shutdown).
    void destroy_gl();

private:
    enum class Mode { None, Preview, Engine };

    void compute_framing(const spirula::TrainerSession& session);
    // The client-frame default pose (web viewer cam.reset() + orbit(0,-250)).
    void reset_pose(float radius);
    // Frame the scene only when a different dataset arrives; a preview ->
    // engine transition on the same dataset keeps the navigated pose and
    // intrinsics (no jump when training starts).
    void maybe_frame(const spirula::TrainerSession& session);
    void reset_view();
    // Update _moving / _last_move from the camera pose. Runs every frame in
    // every scale mode: what it feeds is no longer only the adaptive scale.
    void note_motion(double now);
    float render_scale();
    // Training iterations between refreshes while the camera is still, sized
    // from the measured render and step costs so the viewport keeps costing
    // about the same slice of the run whatever the dataset.
    int idle_step_interval() const;
    // Render resolution for `avail` at the current scale, aspect PRESERVED.
    void render_size(const ImVec2& avail, int& W, int& H) const;
    // Double-click centering (webgl viewer's recenterAt): pan laterally so
    // p (normalized frame) sits on the optical axis and make it the orbit
    // pivot; if p is behind the camera plane (>180-degree models), rotate
    // toward it instead.
    void recenter_at(const float p[3]);
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
    void draw_engine(bool training, const ImVec2& avail, int step);
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
    // Pending double-click pick, as fractional viewport-image coordinates
    // (set in handle_input; consumed by draw_preview / draw_engine).
    bool _dbl_pending = false;
    float _dbl_u = 0.0f, _dbl_v = 0.0f;

    // Camera model + per-model FOV memory (browser _fovMemory equivalent).
    int _cam_model = 0;              // index into kViewerCameraModels
    float _fov_deg[4] = {90, 180, 180, 0};

    // Whether there are training cameras behind what is being rendered. A
    // splat file has none, so the frustum controls are hidden rather than
    // shown doing nothing.
    bool _has_cameras = true;

    // Render options.
    int _buffer_idx = 0;
    std::vector<std::string> _buffer_keys;
    bool _show_cams = false;
    bool _show_grid = false;         // axes + ground-plane grid overlay
    float _frustum_scale = 1.0f;     // camera-frustum size multiplier
    // 0 = auto (see render_scale), 1 = 50%, 2 = 75%, 3 = 100%
    int _scale_idx = 0;
    float _last_pose[10] = {};       // pos + rot + target, to spot motion
    double _last_move = -1e9;
    bool _moving = false;
    bool _auto_refresh = true;

    // Refresh pacing while training (see idle_step_interval): the step the
    // last render was submitted at, and exponential averages of what a render
    // and a training step each cost.
    int _last_render_step = -1;
    int _last_seen_step = -1;
    double _last_step_time = -1.0;
    double _render_secs = 0.0;
    double _step_secs = 0.0;

    // In-flight request / result texture (engine mode).
    uint64_t _pending = 0;
    double _last_submit = 0.0;
    bool _dirty = true;
    unsigned int _tex = 0;
    int _tex_w = 0, _tex_h = 0;
    std::string _last_error;
};

}  // namespace gui
