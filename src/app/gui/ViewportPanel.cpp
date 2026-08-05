// ViewportPanel.cpp -- see ViewportPanel.h.

#include "app/gui/ViewportPanel.h"

#include "app/TrainerCore.h"
#include "app/gui/Ui.h"

#include "i18n/catalog/Gui.h"
#include "app/gui/GlLoader.h"   // GL types + 1.1 entry points

#include "imgui.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "core/Env.h"

namespace msg = spirula::i18n::msg::gui;

namespace gui {

namespace {

// Camera models offered by the web client (viewer.html #camera-model), with
// per-model FOV ranges (viewer.html FOV_RANGE; EQUIRECTANGULAR is a fixed
// full-sphere view with no FOV).
struct ViewerCamModel {
    const char* name;      // engine camera_model string
    const char* label;
    float fov_min, fov_max;
};
const ViewerCamModel kViewerCameraModels[] = {
    {"PINHOLE",         "Pinhole",                 10, 150},
    {"FISHEYE",         "Fisheye (Equidistant)",   10, 360},
    {"EQUISOLID",       "Fisheye (Equisolid)",     10, 360},
    {"EQUIRECTANGULAR", "Equirectangular (360\xc2\xb0)", 0, 0},
};

// fovToIntrinsics port: inverting the actual projection (not the pinhole tan
// map) makes the slider a true field of view for the fisheye models.
void fov_to_intrinsics(float fov_deg, int w, int h, const char* model,
                       float& fx, float& fy) {
    float fov_rad = fov_deg * 3.14159265358979f / 180.0f;
    if (!std::strcmp(model, "FISHEYE"))
        fx = (w / 2.0f) / (fov_rad / 2.0f);                     // r = f*theta
    else if (!std::strcmp(model, "EQUISOLID"))
        fx = (w / 2.0f) / (2.0f * std::sin(fov_rad / 4.0f));    // r = 2f sin(theta/2)
    else
        fx = (w / 2.0f) / std::tan(fov_rad / 2.0f);             // pinhole
    fy = fx;
}

}  // namespace

// ---------------------------------------------------------------------------
// Framing + camera math
// ---------------------------------------------------------------------------

void ViewportPanel::compute_framing(const spirula::TrainerSession& session) {
    // Match the web viewer's cam.reset() exactly: target = client-frame
    // origin (the normalized frame is centered on the CAMERA POSES via
    // center_method="poses", i.e. the captured object for object-centric
    // datasets -- NOT the point-cloud centroid, which distant background
    // points drag far away), pos = [0,0,1], then orbit(0, -250).
    _cam.pos[0] = 0; _cam.pos[1] = 0; _cam.pos[2] = 1;
    _cam.rot[0] = _cam.rot[1] = _cam.rot[2] = 0; _cam.rot[3] = 1;
    _cam.target[0] = _cam.target[1] = _cam.target[2] = 0;
    _cam.orbit(0, -250);
    _home = _cam;

    // Scene radius (drives only the preview depth range): spread of the
    // camera positions in the client frame.
    double A[16] = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
    const auto& ds = session.ds;
    if (ds.train_frame_scale != 1.0f) {
        double T[16];
        for (int i = 0; i < 16; i++) T[i] = ds.train_to_normalized[i];
        dsparse::invert_affine4x4(T, A);
    }
    double radius = 1.0;
    for (int64_t i = 0; i < ds.num_cameras; i++) {
        float p[3] = {ds.c2w[i*12 + 3], ds.c2w[i*12 + 7], ds.c2w[i*12 + 11]};
        double q[3];
        for (int r = 0; r < 3; r++)
            q[r] = A[r*4+0]*p[0] + A[r*4+1]*p[1] + A[r*4+2]*p[2] + A[r*4+3];
        radius = std::max(radius, std::sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]));
    }
    _home_dist = (float)radius;
    _dirty = true;
}

void ViewportPanel::maybe_frame(const spirula::TrainerSession& session) {
    std::string key = session.cfg.data + ":" +
        std::to_string(session.ds.num_cameras) + ":" +
        std::to_string(session.ds.points.num());
    if (key == _framed_key) return;   // same dataset: keep the current pose
    _framed_key = key;
    compute_framing(session);
    _show_cams = true;   // default on for a fresh dataset preview
}

void ViewportPanel::reset_view() {
    _cam = _home;
    _dirty = true;
}

void ViewportPanel::recenter_at(const float p[3]) {
    // viewer.html recenterAt: lateral pan keeps the view orientation; the
    // point becomes the orbit pivot.
    float fwd[3];
    _cam.axis_forward(fwd);
    float rel[3] = {p[0] - _cam.pos[0], p[1] - _cam.pos[1],
                    p[2] - _cam.pos[2]};
    float d = rel[0]*fwd[0] + rel[1]*fwd[1] + rel[2]*fwd[2];
    if (d > 0.0f) {
        for (int i = 0; i < 3; i++)
            _cam.pos[i] += rel[i] - d * fwd[i];
    } else {
        float up[3], eye[3] = {_cam.pos[0], _cam.pos[1], _cam.pos[2]};
        _cam.axis_up(up);
        _cam.look_at(eye, p, up);
    }
    for (int i = 0; i < 3; i++) _cam.target[i] = p[i];
    _dirty = true;
}

const char* ViewportPanel::camera_model_name() const {
    return kViewerCameraModels[_cam_model].name;
}
float ViewportPanel::fov_min() const { return kViewerCameraModels[_cam_model].fov_min; }
float ViewportPanel::fov_max() const { return kViewerCameraModels[_cam_model].fov_max; }

void ViewportPanel::compute_intrinsics(int W, int H, float& fx, float& fy) const {
    const char* model = camera_model_name();
    if (!std::strcmp(model, "EQUIRECTANGULAR")) {
        // Full 360x180 panorama across the image (viewer.html
        // sendRenderRequest): fx = w/2pi, fy = h/pi.
        fx = (float)W / (2.0f * 3.14159265358979f);
        fy = (float)H / 3.14159265358979f;
    } else {
        fov_to_intrinsics(_fov_deg[_cam_model], W, H, model, fx, fy);
    }
}

float ViewportPanel::nav_dist() const {
    float dx = _cam.pos[0] - _cam.target[0];
    float dy = _cam.pos[1] - _cam.target[1];
    float dz = _cam.pos[2] - _cam.target[2];
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

void ViewportPanel::build_request(ViewRequest& q, int W, int H) const {
    _cam.c2w(q.c2w);
    compute_intrinsics(W, H, q.fx, q.fy);
    q.cx = 0.5f * (float)W;
    q.cy = 0.5f * (float)H;
    q.W = W;
    q.H = H;
    q.model = camera_model_name();
    q.key = _buffer_keys.empty() ? "rgb" : _buffer_keys[_buffer_idx];
    q.show_cams = _show_cams;
    q.show_grid = _show_grid;
    q.grid_dist = nav_dist();
    for (int i = 0; i < 3; i++) q.grid_target[i] = _cam.target[i];
    q.cam_size_scale = _frustum_scale;
}

// Row-major world-to-view for the GL preview (inverse of the camera c2w).
void ViewportPanel::view_matrix(float out[16]) const {
    float m[12];
    _cam.c2w(m);
    // c2w columns are the view axes; view = [R^T | -R^T pos].
    for (int r = 0; r < 3; r++) {
        float ax = m[0*4 + r], ay = m[1*4 + r], az = m[2*4 + r];
        out[r*4 + 0] = ax;
        out[r*4 + 1] = ay;
        out[r*4 + 2] = az;
        out[r*4 + 3] = -(ax*m[3] + ay*m[7] + az*m[11]);
    }
    out[12] = out[13] = out[14] = 0;
    out[15] = 1;
}

// ---------------------------------------------------------------------------
// Attach / detach
// ---------------------------------------------------------------------------

void ViewportPanel::attach_preview(spirula::TrainerSession& session) {
    detach();
    if (!_preview.build(session)) {
        _last_error = "preview renderer unavailable (OpenGL 3.2 required)";
        return;
    }
    maybe_frame(session);
    _last_error.clear();
    _mode = Mode::Preview;
}

void ViewportPanel::attach(spirula::TrainerSession& session) {
    detach();
    _worker.start(session.make_viewer_config(), session.make_viewer_hooks());
    _buffer_keys = _worker.buffer_keys();
    _buffer_idx = std::min<int>(_buffer_idx, (int)_buffer_keys.size() - 1);
    maybe_frame(session);
    _pending = 0;
    _last_error.clear();
    _mode = Mode::Engine;
}

void ViewportPanel::detach() {
    if (_mode == Mode::Engine) _worker.stop();
    _pending = 0;
    _mode = Mode::None;
}

void ViewportPanel::destroy_gl() {
    _preview.destroy_gl();
    if (_tex) {
        GLuint t = (GLuint)_tex;
        glDeleteTextures(1, &t);
        _tex = 0;
    }
}

void ViewportPanel::upload(const ViewResult& res) {
    if (!_tex) {
        GLuint t = 0;
        glGenTextures(1, &t);
        _tex = t;
        glBindTexture(GL_TEXTURE_2D, (GLuint)_tex);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    }
    glBindTexture(GL_TEXTURE_2D, (GLuint)_tex);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, res.W, res.H, 0,
                 GL_RGB, GL_UNSIGNED_BYTE, res.rgb8.data());
    _tex_w = res.W;
    _tex_h = res.H;
}

// ---------------------------------------------------------------------------
// Input (browser-identical: viewer.html event wiring)
// ---------------------------------------------------------------------------

void ViewportPanel::handle_input(float /*item_h*/) {
    ImGuiIO& io = ImGui::GetIO();
    bool hovered = ImGui::IsItemHovered();

    // Pointer (mouse; single-touch and OS touch/trackpad gestures arrive as
    // emulated mouse + wheel events -- one finger orbits, two-finger
    // pan/pinch maps to the wheel). Drag mapping per viewer.html:
    //   MMB / RMB / Shift+drag -> pan
    //   LMB: orbit (turntable/trackball) or look (fps/fly)
    if (hovered && !_dragging) {
        for (int b : {ImGuiMouseButton_Left, ImGuiMouseButton_Right,
                      ImGuiMouseButton_Middle}) {
            if (ImGui::IsMouseClicked(b)) {
                _dragging = true;
                _drag_button = b;
                break;
            }
        }
    }
    if (_dragging) {
        if (!ImGui::IsMouseDown(_drag_button)) {
            _dragging = false;
        } else {
            float dx = io.MouseDelta.x, dy = io.MouseDelta.y;
            if (dx != 0 || dy != 0) {
                bool is_pan = _drag_button == ImGuiMouseButton_Right ||
                              _drag_button == ImGuiMouseButton_Middle ||
                              io.KeyShift;
                if (spirula::env("NAV_DEBUG"))
                    std::fprintf(stderr,
                        "[nav] btn=%d pan=%d shift=%d d=(%.0f,%.0f) tgt=(%.3f,%.3f,%.3f) pos=(%.3f,%.3f,%.3f)\n",
                        _drag_button, (int)is_pan, (int)io.KeyShift, dx, dy,
                        _cam.target[0], _cam.target[1], _cam.target[2],
                        _cam.pos[0], _cam.pos[1], _cam.pos[2]);
                if (is_pan)
                    _cam.pan(dx, dy);
                else if (_cam.mode == NavCamera::Turntable ||
                         _cam.mode == NavCamera::Trackball)
                    _cam.orbit(dx, dy);
                else
                    _cam.look(dx, dy);
                _dirty = true;
            }
        }
    }

    // Double-click = center the view on the 3D content under the cursor
    // (viewer.html pickRecenter). Recorded here as fractional image
    // coordinates; resolved by the mode-specific draw (preview: CPU point
    // pick, engine: depth readback on the next render).
    if (hovered && ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) {
        ImVec2 mn = ImGui::GetItemRectMin();
        ImVec2 sz = ImGui::GetItemRectSize();
        if (sz.x > 0 && sz.y > 0) {
            _dbl_u = (io.MousePos.x - mn.x) / sz.x;
            _dbl_v = (io.MousePos.y - mn.y) / sz.y;
            _dbl_pending = true;
            _dirty = true;
        }
    }

    // Scroll = dolly (browser wheel deltaY is ~+-100 per notch, ImGui is
    // +-1 with the opposite sign convention).
    if (hovered && io.MouseWheel != 0.0f) {
        _cam.dolly(-io.MouseWheel * 100.0f);
        _dirty = true;
    }

    // Keyboard (viewer.html listens on the window; here: while the pointer
    // is over the viewport or dragging, and no text field wants input).
    if ((hovered || _dragging) && !io.WantTextInput) {
        NavCamera::Keys k;
        k.w = ImGui::IsKeyDown(ImGuiKey_W);
        k.a = ImGui::IsKeyDown(ImGuiKey_A);
        k.s = ImGui::IsKeyDown(ImGuiKey_S);
        k.d = ImGui::IsKeyDown(ImGuiKey_D);
        k.e = ImGui::IsKeyDown(ImGuiKey_E);
        k.q = ImGui::IsKeyDown(ImGuiKey_Q);
        k.up = ImGui::IsKeyDown(ImGuiKey_UpArrow);
        k.down = ImGui::IsKeyDown(ImGuiKey_DownArrow);
        k.left = ImGui::IsKeyDown(ImGuiKey_LeftArrow);
        k.right = ImGui::IsKeyDown(ImGuiKey_RightArrow);
        float dt = std::min(io.DeltaTime, 0.1f);
        if (_cam.keyboard_tick(dt, k)) _dirty = true;
    }

    // Gamepad: always active, like the browser's gamepadTick loop.
    {
        float dt = std::min(io.DeltaTime, 0.1f);
        if (_cam.gamepad_tick(dt)) _dirty = true;
    }
}

// ---------------------------------------------------------------------------
// Draw
// ---------------------------------------------------------------------------

void ViewportPanel::draw_controls(bool engine) {
    // All controls fit on one row on a wide enough viewport; otherwise the
    // navigation/camera group wraps to a second row.
    bool one_row = ImGui::GetContentRegionAvail().x >= 1190.0f;

    // ---- display group ----
    if (engine && !_buffer_keys.empty()) {
        ImGui::SetNextItemWidth(130);
        // Buffer names ("color", "depth") come from the engine.
        if (ui::BeginComboRaw("##buffer", _buffer_keys[_buffer_idx].c_str())) {
            for (int i = 0; i < (int)_buffer_keys.size(); i++)
                if (ui::SelectableRaw(_buffer_keys[i], i == _buffer_idx)) {
                    _buffer_idx = i;
                    _dirty = true;
                }
            ImGui::EndCombo();
        }
        ui::help_on_hover(msg::viewport_buffer_help);
        ImGui::SameLine();
    }
    if (!engine) {
        ui::TextDisabled(msg::viewport_dataset_preview);
        ui::help_on_hover(msg::viewport_dataset_preview_help);
        ImGui::SameLine();
    }
    if (ui::Checkbox(msg::viewport_cameras, &_show_cams)) _dirty = true;
    if (_show_cams) {
        ImGui::SameLine();
        ImGui::SetNextItemWidth(110);
        if (ui::SliderFloatRaw("##fsize", &_frustum_scale, 0.1f, 10.0f,
                               "size x%.2f", ImGuiSliderFlags_Logarithmic))
            _dirty = true;
        ui::help_on_hover(msg::viewport_frustum_size_help);
    }
    ImGui::SameLine();
    if (ui::Checkbox(msg::viewport_grid, &_show_grid)) _dirty = true;
    ui::help_on_hover(msg::viewport_cameras_help);
    if (engine) {
        ImGui::SameLine();
        ImGui::SetNextItemWidth(66);
        const char* scales[] = {"50%", "75%", "100%"};
        if (ui::ComboRaw("##scale", &_scale_idx, scales, 3)) _dirty = true;
        ui::help_on_hover(msg::viewport_scale_help);
        ImGui::SameLine();
        ui::Checkbox(msg::viewport_live, &_auto_refresh);
        ui::help_on_hover(msg::viewport_live_help);
    }
    ImGui::SameLine();
    if (ui::Button(msg::viewport_reset_view)) reset_view();

    // ---- navigation + camera group ----
    // The navigation modes and camera models are named exactly as the web
    // viewer names them, and the tooltip below says so; renaming them per
    // language would break that correspondence.
    if (one_row) ImGui::SameLine();
    const char* nav_modes[] = {"Turntable", "Trackball", "First Person", "Free Fly"};
    int nav = (int)_cam.mode;
    ImGui::SetNextItemWidth(110);
    if (ui::ComboRaw("##navmode", &nav, nav_modes, 4))
        _cam.mode = (NavCamera::Mode)nav;
    ui::help_on_hover(msg::viewport_nav_help);
    ImGui::SameLine();
    ImGui::SetNextItemWidth(100);
    ui::SliderFloatRaw("##speed", &_cam.speed_exp, -2.0f, 2.0f, "spd 10^%.1f");
    ui::help_on_hover(msg::viewport_speed_help);
    ImGui::SameLine();
    ImGui::SetNextItemWidth(160);
    if (ui::BeginComboRaw("##cammodel", kViewerCameraModels[_cam_model].label)) {
        for (int i = 0; i < 4; i++)
            if (ui::SelectableRaw(kViewerCameraModels[i].label, i == _cam_model)) {
                _cam_model = i;
                // Clamp the remembered FOV into the new model's range
                // (viewer.html _syncCameraModelUI).
                if (fov_max() > 0)
                    _fov_deg[i] = std::clamp(_fov_deg[i], fov_min(), fov_max());
                _dirty = true;
            }
        ImGui::EndCombo();
    }
    ui::help_on_hover(msg::viewport_projection_help);
    ImGui::SameLine();
    if (fov_max() > 0) {
        ImGui::SetNextItemWidth(120);
        if (ui::SliderFloatRaw("##fov", &_fov_deg[_cam_model],
                               fov_min(), fov_max(), "fov %.0f\xc2\xb0"))
            _dirty = true;
    } else {
        ui::TextDisabledRaw("360\xc2\xb0 x 180\xc2\xb0");
    }
}

void ViewportPanel::draw(bool training) {
    draw_controls(_mode == Mode::Engine);

    ImVec2 avail = ImGui::GetContentRegionAvail();
    avail.x = std::max(avail.x, 64.0f);
    avail.y = std::max(avail.y, 64.0f);

    if (_mode == Mode::None) {
        ImGui::Dummy(ImVec2(avail.x, avail.y * 0.4f));
        const char* line = _last_error.empty()
            ? msg::viewport_open_a_dataset.get() : _last_error.c_str();
        float tw = ImGui::CalcTextSize(line).x;
        ImGui::SetCursorPosX(std::max(0.0f, (ImGui::GetWindowWidth() - tw) * 0.5f));
        ui::TextDisabledRaw(line);
        return;
    }

    if (_mode == Mode::Preview) draw_preview(avail);
    else                        draw_engine(training, avail);
}

void ViewportPanel::draw_preview(const ImVec2& avail) {
    int W = (int)std::clamp(avail.x, 64.0f, 4096.0f);
    int H = (int)std::clamp(avail.y, 64.0f, 4096.0f);
    float view[16];
    view_matrix(view);
    // Same intrinsics as the engine render request -> seamless transition.
    float fx, fy;
    compute_intrinsics(W, H, fx, fy);
    unsigned tex = _preview.render(W, H, view,
                                   (PreviewProjection)_cam_model,
                                   fx / (0.5f * W), fy / (0.5f * H),
                                   _home_dist, nav_dist(), _cam.target,
                                   _show_cams, _frustum_scale, _show_grid);
    if (!tex) {
        ui::TextDisabled(msg::viewport_render_failed);
        return;
    }
    // FBO textures are bottom-up; flip V.
    ImGui::Image((ImTextureID)(intptr_t)tex, ImVec2((float)W, (float)H),
                 ImVec2(0, 1), ImVec2(1, 0));
    handle_input((float)H);

    // Double-click centering: CPU pick against the displayed point cloud
    // (nearest point along the cursor ray, 3% angular cone).
    if (_dbl_pending) {
        _dbl_pending = false;
        float u = (_dbl_u - 0.5f) * (float)W / fx;
        float v = (_dbl_v - 0.5f) * (float)H / fy;
        float dcv[3];
        if (viewer_pixel_ray(_cam_model, u, v, dcv)) {
            float m[12];
            _cam.c2w(m);
            float ro[3] = {m[3], m[7], m[11]}, rd[3], p[3];
            // CV ray -> world through the OpenGL c2w (y/z column flip).
            for (int r = 0; r < 3; r++)
                rd[r] = m[r*4+0]*dcv[0] - m[r*4+1]*dcv[1] - m[r*4+2]*dcv[2];
            if (_preview.pick_point(ro, rd, p)) recenter_at(p);
        }
    }

    char info[96];
    std::snprintf(info, sizeof info, "%lld points",
                  (long long)_preview.num_points());
    ImVec2 p = ImGui::GetItemRectMin();
    ImGui::GetWindowDrawList()->AddText(ImVec2(p.x + 8, p.y + 6),
                                        IM_COL32(200, 200, 200, 180), info);
}

void ViewportPanel::draw_engine(bool training, const ImVec2& avail) {
    // Poll the in-flight render.
    if (_pending) {
        ViewResult res;
        if (_worker.try_get_result(_pending, res)) {
            _pending = 0;
            if (res.error.empty()) {
                upload(res);
                _last_error.clear();
                // Double-click centering result (background clicks miss).
                if (res.pick_hit) recenter_at(res.pick_point);
            } else {
                _last_error = res.error;
            }
        }
    }

    // Submit a fresh render when needed.
    float scale = _scale_idx == 0 ? 0.5f : _scale_idx == 1 ? 0.75f : 1.0f;
    int W = (int)std::clamp(avail.x * scale, 64.0f, 1920.0f);
    int H = (int)std::clamp(avail.y * scale, 64.0f, 1920.0f);
    double now = ImGui::GetTime();
    bool want = _dirty || (training && _auto_refresh && now - _last_submit > 0.15);
    if (!_pending && want) {
        ViewRequest q;
        build_request(q, W, H);
        // Attach a pending double-click pick to this render; the picked
        // point comes back with the result (depth readback, no extra pass).
        if (_dbl_pending) {
            q.pick_px = (int)(_dbl_u * (float)W);
            q.pick_py = (int)(_dbl_v * (float)H);
            _dbl_pending = false;
        }
        _pending = _worker.submit(q);
        _last_submit = now;
        _dirty = false;
    }

    // Draw the last result, aspect-fit.
    if (_tex && _tex_w > 0) {
        float ar_img = (float)_tex_w / (float)_tex_h;
        float ar_avail = avail.x / avail.y;
        ImVec2 size = ar_img > ar_avail
            ? ImVec2(avail.x, avail.x / ar_img)
            : ImVec2(avail.y * ar_img, avail.y);
        ImVec2 pad((avail.x - size.x) * 0.5f, (avail.y - size.y) * 0.5f);
        ImGui::SetCursorPos(ImVec2(ImGui::GetCursorPosX() + pad.x,
                                   ImGui::GetCursorPosY() + pad.y));
        ImGui::Image((ImTextureID)(intptr_t)_tex, size);
        handle_input(size.y);
    } else {
        ImGui::Dummy(ImVec2(avail.x, avail.y * 0.4f));
        const char* line = _last_error.empty() ? msg::viewport_rendering.get()
                                               : _last_error.c_str();
        float tw = ImGui::CalcTextSize(line).x;
        ImGui::SetCursorPosX(std::max(0.0f, (ImGui::GetWindowWidth() - tw) * 0.5f));
        ui::TextDisabledRaw(line);
    }
    if (!_last_error.empty() && _tex)
        ui::TextColored(ImVec4(1, 0.4f, 0.4f, 1), msg::viewport_render_error,
                        {_last_error});
}

}  // namespace gui
