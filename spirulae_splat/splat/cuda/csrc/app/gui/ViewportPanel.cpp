// ViewportPanel.cpp -- see ViewportPanel.h.

#include "ViewportPanel.h"

#include "../TrainerCore.h"
#include "ConfigUI.h"   // help_tooltip_on_hover

#include "imgui.h"

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <GL/gl.h>
#elif defined(__APPLE__)
#define GL_SILENCE_DEPRECATION
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

#include <algorithm>
#include <cmath>
#include <cstring>

namespace gui {

namespace {

void v3_normalize(float v[3]) {
    float n = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    if (n < 1e-12f) { v[0] = 1; v[1] = 0; v[2] = 0; return; }
    v[0] /= n; v[1] /= n; v[2] /= n;
}
void v3_cross(const float a[3], const float b[3], float out[3]) {
    out[0] = a[1]*b[2] - a[2]*b[1];
    out[1] = a[2]*b[0] - a[0]*b[2];
    out[2] = a[0]*b[1] - a[1]*b[0];
}

}  // namespace

void ViewportPanel::attach(ssplat::TrainerSession& session) {
    detach();
    _worker.start(session.make_viewer_config(), session.make_viewer_hooks());
    _buffer_keys = _worker.buffer_keys();
    _buffer_idx = std::min<int>(_buffer_idx, (int)_buffer_keys.size() - 1);

    // ---- frame the scene in the client (normalized) frame ------------------
    // The render path remaps client c2w through train_to_normalized when
    // train_frame_scale != 1, so positions must be framed in that space:
    // p_normalized = inv(train_to_normalized) * p_train.
    double A[16] = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
    const auto& ds = session.ds;
    if (ds.train_frame_scale != 1.0f) {
        double T[16];
        for (int i = 0; i < 16; i++) T[i] = ds.train_to_normalized[i];
        dsparse::invert_affine4x4(T, A);
    }
    auto map_pt = [&](const float p[3], double out[3]) {
        for (int r = 0; r < 3; r++)
            out[r] = A[r*4+0]*p[0] + A[r*4+1]*p[1] + A[r*4+2]*p[2] + A[r*4+3];
    };

    // Target = centroid of the seed points (the scene), sampled.
    double centroid[3] = {0, 0, 0};
    int64_t n_pts = ds.points.num();
    if (n_pts > 0) {
        int64_t stride = std::max<int64_t>(1, n_pts / 10000);
        int64_t cnt = 0;
        for (int64_t i = 0; i < n_pts; i += stride, cnt++) {
            double q[3];
            map_pt(&ds.points.xyz[i*3], q);
            for (int r = 0; r < 3; r++) centroid[r] += q[r];
        }
        for (int r = 0; r < 3; r++) centroid[r] /= std::max<int64_t>(cnt, 1);
    }

    // Distance / start direction from the camera ring.
    double cam_c[3] = {0, 0, 0};
    std::vector<double> dists;
    dists.reserve(ds.num_cameras);
    for (int64_t i = 0; i < ds.num_cameras; i++) {
        float p[3] = {ds.c2w[i*12 + 3], ds.c2w[i*12 + 7], ds.c2w[i*12 + 11]};
        double q[3];
        map_pt(p, q);
        for (int r = 0; r < 3; r++) cam_c[r] += q[r];
        double dx = q[0]-centroid[0], dy = q[1]-centroid[1], dz = q[2]-centroid[2];
        dists.push_back(std::sqrt(dx*dx + dy*dy + dz*dz));
    }
    double dist = 3.0;
    if (!dists.empty()) {
        for (int r = 0; r < 3; r++) cam_c[r] /= (double)ds.num_cameras;
        std::nth_element(dists.begin(), dists.begin() + dists.size()/2, dists.end());
        dist = std::max(1e-3, 1.4 * dists[dists.size()/2]);
    }

    for (int r = 0; r < 3; r++)
        _home_target[r] = _target[r] = (float)centroid[r];
    _home_dist = _dist = (float)dist;
    double dir[3] = {cam_c[0]-centroid[0], cam_c[1]-centroid[1], cam_c[2]-centroid[2]};
    double dn = std::sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
    if (dn > 1e-9) {
        _yaw = (float)std::atan2(dir[1], dir[0]);
        _pitch = (float)std::asin(std::clamp(dir[2] / dn, -0.95, 0.95));
    } else {
        _yaw = 0.0f;
        _pitch = 0.6f;
    }

    _pending = 0;
    _dirty = true;
    _last_error.clear();
    _attached = true;
}

void ViewportPanel::detach() {
    if (!_attached) return;
    _worker.stop();
    _pending = 0;
    _attached = false;
}

void ViewportPanel::reset_view() {
    std::memcpy(_target, _home_target, sizeof _target);
    _dist = _home_dist;
    _dirty = true;
}

void ViewportPanel::destroy_gl() {
    if (_tex) {
        GLuint t = (GLuint)_tex;
        glDeleteTextures(1, &t);
        _tex = 0;
    }
}

void ViewportPanel::build_request(ViewRequest& q, int W, int H) const {
    // Orbit camera -> OpenGL-convention c2w in the Z-up client frame.
    float cp = std::cos(_pitch), sp = std::sin(_pitch);
    float cy = std::cos(_yaw), sy = std::sin(_yaw);
    float eye[3] = {_target[0] + _dist * cp * cy,
                    _target[1] + _dist * cp * sy,
                    _target[2] + _dist * sp};
    float f[3] = {_target[0] - eye[0], _target[1] - eye[1], _target[2] - eye[2]};
    v3_normalize(f);
    const float up_w[3] = {0, 0, 1};
    float r[3], u[3];
    v3_cross(f, up_w, r);
    v3_normalize(r);
    v3_cross(r, f, u);
    for (int i = 0; i < 3; i++) {
        q.c2w[i*4 + 0] = r[i];
        q.c2w[i*4 + 1] = u[i];
        q.c2w[i*4 + 2] = -f[i];
        q.c2w[i*4 + 3] = eye[i];
    }
    float fy = 0.5f * (float)H / std::tan(0.5f * _fov_deg * 3.14159265f / 180.0f);
    q.fx = fy;
    q.fy = fy;
    q.cx = 0.5f * (float)W;
    q.cy = 0.5f * (float)H;
    q.W = W;
    q.H = H;
    q.model = "PINHOLE";
    q.key = _buffer_keys.empty() ? "rgb" : _buffer_keys[_buffer_idx];
    q.show_cams = _show_cams;
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

void ViewportPanel::draw(bool training) {
    // ---- controls row -------------------------------------------------------
    if (!_buffer_keys.empty()) {
        ImGui::SetNextItemWidth(140);
        if (ImGui::BeginCombo("##buffer", _buffer_keys[_buffer_idx].c_str())) {
            for (int i = 0; i < (int)_buffer_keys.size(); i++)
                if (ImGui::Selectable(_buffer_keys[i].c_str(), i == _buffer_idx)) {
                    _buffer_idx = i;
                    _dirty = true;
                }
            ImGui::EndCombo();
        }
        help_tooltip_on_hover("Which render buffer to display (color, depth, "
                              "alpha, normals from depth, ...).");
        ImGui::SameLine();
    }
    if (ImGui::Checkbox("cameras", &_show_cams)) _dirty = true;
    help_tooltip_on_hover("Overlay the training camera frusta (with image "
                          "thumbnails once training has visited them).");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(70);
    const char* scales[] = {"50%", "75%", "100%"};
    if (ImGui::Combo("##scale", &_scale_idx, scales, 3)) _dirty = true;
    help_tooltip_on_hover("Render resolution relative to the viewport size. "
                          "Lower is faster and steals less time from training.");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(110);
    if (ImGui::SliderFloat("##fov", &_fov_deg, 20.0f, 120.0f, "fov %.0f\xc2\xb0"))
        _dirty = true;
    ImGui::SameLine();
    if (ImGui::Button("Reset view")) reset_view();
    ImGui::SameLine();
    ImGui::Checkbox("live", &_auto_refresh);
    help_tooltip_on_hover("Continuously re-render while training so the "
                          "viewport follows the optimization.");

    // ---- image region --------------------------------------------------------
    ImVec2 avail = ImGui::GetContentRegionAvail();
    avail.x = std::max(avail.x, 64.0f);
    avail.y = std::max(avail.y, 64.0f);

    if (!_attached) {
        ImGui::Dummy(ImVec2(avail.x, avail.y * 0.4f));
        const char* msg = "Viewport activates when training starts";
        float tw = ImGui::CalcTextSize(msg).x;
        ImGui::SetCursorPosX((ImGui::GetWindowWidth() - tw) * 0.5f);
        ImGui::TextDisabled("%s", msg);
        return;
    }

    // Poll the in-flight render.
    if (_pending) {
        ViewResult res;
        if (_worker.try_get_result(_pending, res)) {
            _pending = 0;
            if (res.error.empty()) {
                upload(res);
                _last_error.clear();
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

        // ---- navigation ------------------------------------------------------
        if (ImGui::IsItemHovered()) {
            ImGuiIO& io = ImGui::GetIO();
            if (io.MouseWheel != 0.0f) {
                _dist *= std::exp(-io.MouseWheel * 0.15f);
                _dist = std::clamp(_dist, 1e-4f, 1e5f);
                _dirty = true;
            }
            if (ImGui::IsMouseDragging(ImGuiMouseButton_Left, 0.0f)) {
                _yaw   -= io.MouseDelta.x * 0.008f;
                _pitch += io.MouseDelta.y * 0.008f;
                _pitch = std::clamp(_pitch, -1.55f, 1.55f);
                if (io.MouseDelta.x || io.MouseDelta.y) _dirty = true;
            }
            bool pan = ImGui::IsMouseDragging(ImGuiMouseButton_Right, 0.0f) ||
                       ImGui::IsMouseDragging(ImGuiMouseButton_Middle, 0.0f);
            if (pan && (io.MouseDelta.x || io.MouseDelta.y)) {
                // Pixel-accurate pan at the target depth.
                float ppu = 2.0f * _dist *
                    std::tan(0.5f * _fov_deg * 3.14159265f / 180.0f) / size.y;
                float cp = std::cos(_pitch), sp = std::sin(_pitch);
                float cy = std::cos(_yaw), sy = std::sin(_yaw);
                float f[3] = {-cp * cy, -cp * sy, -sp};
                const float up_w[3] = {0, 0, 1};
                float r[3], u[3];
                v3_cross(f, up_w, r);
                v3_normalize(r);
                v3_cross(r, f, u);
                for (int i = 0; i < 3; i++)
                    _target[i] += (-io.MouseDelta.x * r[i] + io.MouseDelta.y * u[i]) * ppu;
                _dirty = true;
            }
        }
    } else {
        ImGui::Dummy(ImVec2(avail.x, avail.y * 0.4f));
        const char* msg = _last_error.empty() ? "Rendering..." : _last_error.c_str();
        float tw = ImGui::CalcTextSize(msg).x;
        ImGui::SetCursorPosX(std::max(0.0f, (ImGui::GetWindowWidth() - tw) * 0.5f));
        ImGui::TextDisabled("%s", msg);
    }
    if (!_last_error.empty() && _tex)
        ImGui::TextColored(ImVec4(1, 0.4f, 0.4f, 1), "render error: %s",
                           _last_error.c_str());
}

}  // namespace gui
