// PreviewRenderer.cpp -- see PreviewRenderer.h.

#include "PreviewRenderer.h"

#include "GlLoader.h"
#include "../TrainerCore.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace gui {

namespace {

// One shader for both draws. Vertices are (pos, aux); for the point cloud
// aux is the RGB color (u_scale = 0), for camera frusta aux is the offset
// from the camera center, scaled by the live frustum-size uniform
// (u_color.a > 0 selects the flat line color).
// Projection happens in the shader so the preview supports the same camera
// models as the engine viewer: u_model selects pinhole(0) / fisheye-
// equidistant(1) / equisolid(2) / equirectangular(3); u_s is the engine
// intrinsics normalized to NDC (fx/(W/2), fy/(H/2)). Depth is the linear
// view-space distance mapped over u_zrange (fine for points/lines and valid
// for the >180-degree projections where a perspective z is meaningless).
const char* kVert = R"(#version 150
in vec3 a_pos;
in vec3 a_aux;
uniform mat4 u_view;
uniform float u_scale;
uniform vec4 u_color;
uniform int u_model;
uniform vec2 u_s;
uniform vec2 u_zrange;
out vec4 v_col;
void main() {
    vec3 p = a_pos + u_scale * a_aux;
    vec3 v = (u_view * vec4(p, 1.0)).xyz;
    float dist = max(length(v), 1e-9);
    vec2 ndc;
    bool clipped = false;
    if (u_model == 0) {                     // pinhole
        if (v.z > -1e-6) clipped = true;
        ndc = u_s * (v.xy / -v.z);
    } else if (u_model == 3) {              // equirectangular
        float lon = atan(v.x, -v.z);
        float lat = asin(clamp(v.y / dist, -1.0, 1.0));
        ndc = u_s * vec2(lon, lat);
    } else {                                // fisheye / equisolid
        float theta = acos(clamp(-v.z / dist, -1.0, 1.0));
        float rlen = length(v.xy);
        vec2 dir2 = rlen > 1e-9 ? v.xy / rlen : vec2(0.0);
        float r = (u_model == 1) ? theta : 2.0 * sin(0.5 * theta);
        ndc = u_s * dir2 * r;
    }
    float z = (dist - u_zrange.x) / (u_zrange.y - u_zrange.x) * 2.0 - 1.0;
    if (clipped) z = 3.0;
    v_col = (u_color.a > 0.0) ? u_color : vec4(a_aux, 1.0);
    gl_Position = vec4(ndc, z, 1.0);
}
)";

const char* kFrag = R"(#version 150
in vec4 v_col;
out vec4 frag;
void main() { frag = vec4(v_col.rgb, 1.0); }
)";

GLuint compile(GLenum type, const char* src) {
    GLuint sh = glx::CreateShader(type);
    glx::ShaderSource(sh, 1, &src, nullptr);
    glx::CompileShader(sh);
    GLint ok = 0;
    glx::GetShaderiv(sh, GL_COMPILE_STATUS, &ok);
    if (!ok) {
        char log[512];
        glx::GetShaderInfoLog(sh, sizeof log, nullptr, log);
        std::fprintf(stderr, "[preview] shader error: %s\n", log);
        glx::DeleteShader(sh);
        return 0;
    }
    return sh;
}

struct V { float px, py, pz, ax, ay, az; };

}  // namespace

bool PreviewRenderer::ensure_program() {
    if (_prog) return true;
    if (!glx::init()) return false;
    GLuint vs = compile(GL_VERTEX_SHADER, kVert);
    GLuint fs = compile(GL_FRAGMENT_SHADER, kFrag);
    if (!vs || !fs) return false;
    _prog = glx::CreateProgram();
    glx::AttachShader(_prog, vs);
    glx::AttachShader(_prog, fs);
    glx::BindAttribLocation(_prog, 0, "a_pos");
    glx::BindAttribLocation(_prog, 1, "a_aux");
    glx::LinkProgram(_prog);
    glx::DeleteShader(vs);
    glx::DeleteShader(fs);
    GLint ok = 0;
    glx::GetProgramiv(_prog, GL_LINK_STATUS, &ok);
    if (!ok) {
        char log[512];
        glx::GetProgramInfoLog(_prog, sizeof log, nullptr, log);
        std::fprintf(stderr, "[preview] link error: %s\n", log);
        glx::DeleteProgram(_prog);
        _prog = 0;
        return false;
    }
    _u_view = glx::GetUniformLocation(_prog, "u_view");
    _u_scale = glx::GetUniformLocation(_prog, "u_scale");
    _u_color = glx::GetUniformLocation(_prog, "u_color");
    _u_model = glx::GetUniformLocation(_prog, "u_model");
    _u_s = glx::GetUniformLocation(_prog, "u_s");
    _u_zrange = glx::GetUniformLocation(_prog, "u_zrange");
    return true;
}

bool PreviewRenderer::build(const ssplat::TrainerSession& session) {
    destroy_gl();
    if (!ensure_program()) return false;
    const auto& ds = session.ds;

    // train -> normalized frame similarity (identity when scale == 1),
    // matching how the viewport frames the scene.
    double A[16] = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
    if (ds.train_frame_scale != 1.0f) {
        double T[16];
        for (int i = 0; i < 16; i++) T[i] = ds.train_to_normalized[i];
        dsparse::invert_affine4x4(T, A);
    }
    auto map_pt = [&](const float p[3], float out[3]) {
        for (int r = 0; r < 3; r++)
            out[r] = (float)(A[r*4+0]*p[0] + A[r*4+1]*p[1] + A[r*4+2]*p[2] + A[r*4+3]);
    };
    // Rotation+scale part for direction vectors, and the similarity scale.
    double sA = std::sqrt(A[0]*A[0] + A[4]*A[4] + A[8]*A[8]);
    if (sA < 1e-12) sA = 1.0;
    auto map_dir = [&](const double d[3], float out[3]) {
        for (int r = 0; r < 3; r++)
            out[r] = (float)((A[r*4+0]*d[0] + A[r*4+1]*d[1] + A[r*4+2]*d[2]) / sA);
    };

    // ---- point cloud (capped; stride-sampled) -------------------------------
    int64_t n_src = ds.points.num();
    int64_t stride = std::max<int64_t>(1, n_src / 4000000);
    std::vector<V> pts;
    pts.reserve(n_src / stride + 1);
    for (int64_t i = 0; i < n_src; i += stride) {
        V v;
        float p[3];
        map_pt(&ds.points.xyz[i*3], p);
        v.px = p[0]; v.py = p[1]; v.pz = p[2];
        v.ax = ds.points.rgb[i*3 + 0] / 255.0f;
        v.ay = ds.points.rgb[i*3 + 1] / 255.0f;
        v.az = ds.points.rgb[i*3 + 2] / 255.0f;
        pts.push_back(v);
    }
    _num_points = (int64_t)pts.size();

    // ---- camera frusta (per-INPUT cameras, 8 line segments each) ------------
    std::vector<V> cams;
    cams.reserve(ds.num_cameras * 16 * 8);
    for (int64_t i = 0; i < ds.num_cameras; i++) {
        const float* M = &ds.c2w[i*12];
        float c[3];
        float t[3] = {M[3], M[7], M[11]};
        map_pt(t, c);
        // Column basis (OpenGL convention; -Z forward), mapped + unit-scaled.
        float X[3], Y[3], Z[3];
        double dx[3] = {M[0], M[4], M[8]};
        double dy[3] = {M[1], M[5], M[9]};
        double dz[3] = {M[2], M[6], M[10]};
        map_dir(dx, X); map_dir(dy, Y); map_dir(dz, Z);
        float fx = ds.intrins[i*4 + 0], fy = ds.intrins[i*4 + 1];
        float tanx = std::min(0.5f * ds.widths[i]  / std::max(fx, 1.0f), 1.5f);
        float tany = std::min(0.5f * ds.heights[i] / std::max(fy, 1.0f), 1.5f);
        float corner[4][3];
        for (int k = 0; k < 4; k++) {
            float sx = (k & 1) ? tanx : -tanx;
            float sy = (k & 2) ? tany : -tany;
            for (int r = 0; r < 3; r++)
                corner[k][r] = sx * X[r] + sy * Y[r] - Z[r];
        }
        // Subdivided so the edges curve correctly under the fisheye /
        // equirectangular preview projections.
        constexpr int kSub = 8;
        auto emit = [&](const float* a, const float* b) {
            float za[3] = {0, 0, 0};
            if (!a) a = za;
            if (!b) b = za;
            for (int t = 0; t < kSub; t++) {
                float f0 = (float)t / kSub, f1 = (float)(t + 1) / kSub;
                V v;
                v.px = c[0]; v.py = c[1]; v.pz = c[2];
                v.ax = a[0] + (b[0]-a[0])*f0;
                v.ay = a[1] + (b[1]-a[1])*f0;
                v.az = a[2] + (b[2]-a[2])*f0;
                cams.push_back(v);
                v.ax = a[0] + (b[0]-a[0])*f1;
                v.ay = a[1] + (b[1]-a[1])*f1;
                v.az = a[2] + (b[2]-a[2])*f1;
                cams.push_back(v);
            }
        };
        emit(nullptr, corner[0]);
        emit(nullptr, corner[1]);
        emit(nullptr, corner[2]);
        emit(nullptr, corner[3]);
        emit(corner[0], corner[1]);
        emit(corner[1], corner[3]);
        emit(corner[3], corner[2]);
        emit(corner[2], corner[0]);
    }
    _num_cam_verts = (int64_t)cams.size();

    // Base frustum size: the engine's kNN heuristic (train frame) mapped to
    // normalized units.
    _base_cam_size = viewer_camera_size_heuristic(session.post) * (float)sA;

    auto make_vao = [&](unsigned& vao, unsigned& vbo, const std::vector<V>& data) {
        GLuint va = 0, vb = 0;
        glx::GenVertexArrays(1, &va);
        glx::GenBuffers(1, &vb);
        glx::BindVertexArray(va);
        glx::BindBuffer(GL_ARRAY_BUFFER, vb);
        glx::BufferData(GL_ARRAY_BUFFER, (glx::glSizeiptr)(data.size() * sizeof(V)),
                        data.data(), GL_STATIC_DRAW);
        glx::EnableVertexAttribArray(0);
        glx::VertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(V), (void*)0);
        glx::EnableVertexAttribArray(1);
        glx::VertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(V),
                                 (void*)(3 * sizeof(float)));
        glx::BindVertexArray(0);
        vao = va;
        vbo = vb;
    };
    make_vao(_vao_pts, _vbo_pts, pts);
    make_vao(_vao_cam, _vbo_cam, cams);

    _gl_ok = true;
    _built = true;
    return true;
}

bool PreviewRenderer::ensure_fbo(int W, int H) {
    if (_fbo && W == _fbo_w && H == _fbo_h) return true;
    if (!_fbo) {
        GLuint fbo = 0, tex = 0, rb = 0;
        glx::GenFramebuffers(1, &fbo);
        glGenTextures(1, &tex);
        glx::GenRenderbuffers(1, &rb);
        _fbo = fbo; _color_tex = tex; _depth_rb = rb;
    }
    glBindTexture(GL_TEXTURE_2D, (GLuint)_color_tex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, W, H, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
    glx::BindRenderbuffer(GL_RENDERBUFFER, (GLuint)_depth_rb);
    glx::RenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, W, H);
    glx::BindFramebuffer(GL_FRAMEBUFFER, (GLuint)_fbo);
    glx::FramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                              GL_TEXTURE_2D, (GLuint)_color_tex, 0);
    glx::FramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,
                                 GL_RENDERBUFFER, (GLuint)_depth_rb);
    bool ok = glx::CheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE;
    glx::BindFramebuffer(GL_FRAMEBUFFER, 0);
    _fbo_w = W;
    _fbo_h = H;
    return ok;
}

unsigned PreviewRenderer::render(int W, int H, const float view[16],
                                 PreviewProjection proj, float sx, float sy,
                                 float scene_radius, bool show_cams,
                                 float frustum_scale) {
    if (!_built || !_gl_ok || W < 1 || H < 1) return 0;
    if (!ensure_fbo(W, H)) return 0;

    float zn = std::max(1e-5f, 0.002f * scene_radius);
    float zf = std::max(10.0f * zn, 500.0f * scene_radius);

    glx::BindFramebuffer(GL_FRAMEBUFFER, (GLuint)_fbo);
    glViewport(0, 0, W, H);
    glClearColor(0.05f, 0.055f, 0.065f, 1.0f);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glx::UseProgram(_prog);
    glx::UniformMatrix4fv(_u_view, 1, GL_TRUE, view);
    glx::Uniform1i(_u_model, (int)proj);
    glx::Uniform2f(_u_s, sx, sy);
    glx::Uniform2f(_u_zrange, zn, zf);

    // Point cloud (aux = vertex color).
    glPointSize(2.0f);
    glx::Uniform1f(_u_scale, 0.0f);
    glx::Uniform4f(_u_color, 0, 0, 0, 0);
    glx::BindVertexArray(_vao_pts);
    glDrawArrays(GL_POINTS, 0, (GLsizei)_num_points);

    // Camera frusta (aux = offset).
    if (show_cams && _num_cam_verts > 0) {
        glx::Uniform1f(_u_scale, _base_cam_size * std::max(frustum_scale, 1e-3f));
        glx::Uniform4f(_u_color, 1.0f, 0.62f, 0.25f, 1.0f);
        glx::BindVertexArray(_vao_cam);
        glDrawArrays(GL_LINES, 0, (GLsizei)_num_cam_verts);
    }

    glx::BindVertexArray(0);
    glx::UseProgram(0);
    glDisable(GL_DEPTH_TEST);
    glx::BindFramebuffer(GL_FRAMEBUFFER, 0);
    return _color_tex;
}

void PreviewRenderer::destroy_gl() {
    if (_vbo_pts) { GLuint b[2] = {(GLuint)_vbo_pts, (GLuint)_vbo_cam}; glx::DeleteBuffers(2, b); }
    if (_vao_pts) { GLuint a[2] = {(GLuint)_vao_pts, (GLuint)_vao_cam}; glx::DeleteVertexArrays(2, a); }
    _vbo_pts = _vbo_cam = _vao_pts = _vao_cam = 0;
    if (_fbo) {
        GLuint fbo = (GLuint)_fbo;
        glx::DeleteFramebuffers(1, &fbo);
        GLuint tex = (GLuint)_color_tex;
        glDeleteTextures(1, &tex);
        GLuint rb = (GLuint)_depth_rb;
        glx::DeleteRenderbuffers(1, &rb);
        _fbo = _color_tex = _depth_rb = 0;
        _fbo_w = _fbo_h = 0;
    }
    _built = false;
    _num_points = _num_cam_verts = 0;
}

}  // namespace gui
