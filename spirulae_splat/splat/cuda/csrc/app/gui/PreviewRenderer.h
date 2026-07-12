#pragma once

// PreviewRenderer -- pure-OpenGL dataset preview: the SfM sparse point cloud
// (vertex-colored) + training-camera frusta, rendered into an offscreen FBO
// texture. Used by the viewport between "dataset loaded" and "training
// started", when the CUDA engine has nothing to render yet. Geometry is
// built in the same Z-up normalized frame the viewport navigates, and the
// vertex shader implements the same camera models as the engine viewer
// (pinhole / fisheye-equidistant / fisheye-equisolid / equirectangular), so
// the hand-off to the engine render at training start is seamless.

#include <cstdint>

namespace ssplat { class TrainerSession; }

namespace gui {

// Projection selector for render() (matches ViewportPanel's camera-model
// dropdown order / the web client's #camera-model).
enum class PreviewProjection {
    Pinhole = 0, Fisheye, Equisolid, Equirectangular
};

class PreviewRenderer {
public:
    // Build GL buffers from the parsed dataset. Requires a current GL
    // context (GUI thread) and the session's load_dataset() to have
    // completed. Returns false when GL init fails (missing functions).
    bool build(const ssplat::TrainerSession& session);
    bool built() const { return _built; }

    // Render into the internal FBO; returns the color texture (0 on error).
    // view is a row-major 4x4 world-to-view matrix. sx/sy are the engine
    // intrinsics normalized to NDC: fx/(W/2), fy/(H/2). frustum_scale
    // multiplies the base camera size; scene_radius drives the depth range.
    unsigned render(int W, int H, const float view[16],
                    PreviewProjection proj, float sx, float sy,
                    float scene_radius, bool show_cams, float frustum_scale);

    // Base frustum size (kNN heuristic, normalized frame).
    float base_camera_size() const { return _base_cam_size; }
    int64_t num_points() const { return _num_points; }

    void destroy_gl();

private:
    bool ensure_program();
    bool ensure_fbo(int W, int H);

    bool _built = false;
    bool _gl_ok = false;
    unsigned _prog = 0;
    int _u_view = -1, _u_scale = -1, _u_color = -1;
    int _u_model = -1, _u_s = -1, _u_zrange = -1;
    unsigned _vao_pts = 0, _vbo_pts = 0;
    unsigned _vao_cam = 0, _vbo_cam = 0;
    int64_t _num_points = 0;
    int64_t _num_cam_verts = 0;
    float _base_cam_size = 0.1f;

    unsigned _fbo = 0, _color_tex = 0, _depth_rb = 0;
    int _fbo_w = 0, _fbo_h = 0;
};

}  // namespace gui
