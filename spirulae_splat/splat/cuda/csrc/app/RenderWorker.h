#pragma once

// RenderWorker -- the interactive render path shared by the web-viewer HTTP
// server (Viewer.cpp, which JPEG-encodes the result) and the native GUI
// viewport (gui/, which uploads it to an OpenGL texture). Extracted from
// Viewer.cpp; ports:
//   viewer/render_worker.py latest-wins pending slot + result handoff
//   trainer.py _render/render  c2w remap + engine render + annotation call
//   model.py get_outputs    (viewer subset: rgb/depth/alpha/normals/median/
//                            distortion display transforms)
//   viewer/annotation.py    engine_viewer_init upload + engine_blit_view call
//                           + the camera_size kNN heuristic
//
// Threading model: submit() stores the request in a single latest-wins slot
// and synchronously flips hooks.set_render_pending(true) so the training
// loop yields the engine mutex at its next iteration boundary (the Python
// _render_pending fairness dance, trainer.py:146-153). The worker thread
// takes hooks.engine_mutex around all engine work and publishes a ViewResult;
// callers either block (wait_result, HTTP thread) or poll (try_get_result,
// GUI frame loop).

#include "DatasetParser.h"

#include <array>
#include <atomic>
#include <cstdint>
#include <functional>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <vector>

// Static per-run info the render path needs (subset of model config).
struct ViewerRenderConfig {
    std::string primitive = "3dgs";
    bool  packed = false;
    int   sh_degree_warmup_every = 1000;
    std::optional<float> relative_scale;
    // Median-depth buffers are rendered only when a median loss is
    // configured (model.py get_outputs:1103-1108); requesting them
    // otherwise 400s, matching Python.
    bool  output_median = false;
    // Distortion render enabled when a distortion regularizer is configured;
    // also enabled on demand when a *_distortion buffer is requested.
    bool  distortion_reg_on = false;
    // Viewer-client c2w remap into the training frame (trainer.py:557-576).
    float train_frame_scale = 1.0f;
    std::array<float, 16> train_to_normalized{1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
    // Unscaled frustum size from the camera-kNN heuristic (the return value
    // of viewer_upload_cameras). Multiplied by ViewRequest::cam_size_scale
    // per render; 0 disables live frustum-size updates.
    float base_camera_size = 0.0f;
};

struct ViewerHooks {
    // Serializes engine access between the training loop and viewer renders.
    std::mutex* engine_mutex = nullptr;
    // Current training step (for the SH-degree warmup schedule).
    std::function<int()> current_step;
    // JSON body for /progress (built by the training loop, which owns the
    // timing state).
    std::function<std::string()> progress_json;
    // Toggle pause; returns the new paused state.
    std::function<bool()> pause_toggle;
    // "Render desired" flag: called with true synchronously from the
    // submitting thread, with false from the worker once the queue drains.
    std::function<void(bool)> set_render_pending;
};

// One interactive render request. c2w is a row-major 3x4 camera-to-world in
// the OpenGL/nerfstudio convention, expressed in the client's normalized
// frame (the worker remaps through cfg.train_to_normalized).
struct ViewRequest {
    float c2w[12] = {1,0,0,0, 0,1,0,0, 0,0,1,0};
    float fx = 500, fy = 500, cx = 256, cy = 256;
    int W = 512, H = 512;
    std::string model = "PINHOLE";
    std::string key = "rgb";
    bool show_cams = false;
    bool show_grid = false;   // axes/grid overlay (viewer_upload_grid)
    // Nav view distance (camera to orbit target) in the client's normalized
    // frame; drives the zoom-adaptive grid cell size. 0 = unknown (the grid
    // keeps its current spacing and center).
    float grid_dist = 0.0f;
    // Nav orbit target (client's normalized frame); the grid's line patch
    // recenters on it so it stays under the viewer at any zoom.
    float grid_target[3] = {0, 0, 0};
    // Frustum-size multiplier applied to cfg.base_camera_size when
    // show_cams is on (user-facing "camera size" slider).
    float cam_size_scale = 1.0f;
    // Double-click pick: when >= 0, also return the 3D point under this
    // pixel (render resolution), read from the ray-depth channel the render
    // already downloads -- no extra VRAM or render pass.
    int pick_px = -1, pick_py = -1;
};

struct ViewResult {
    uint64_t id = 0;
    int W = 0, H = 0;
    std::vector<uint8_t> rgb8;   // [H, W, 3]
    std::string error;           // non-empty on failure
    // Pick result (ViewRequest::pick_px/py): point under the pixel in the
    // client's normalized frame; pick_hit false = background / invalid ray.
    bool pick_hit = false;
    float pick_point[3] = {0, 0, 0};
};

class RenderWorker {
public:
    RenderWorker();
    ~RenderWorker();

    void start(ViewerRenderConfig cfg, ViewerHooks hooks);
    void stop();

    // Latest-wins submit; returns the request id.
    uint64_t submit(const ViewRequest& q);
    // Block up to timeout_s for the result of request `id`.
    bool wait_result(uint64_t id, ViewResult& out, double timeout_s);
    // Non-blocking: true when the result of request `id` is available.
    bool try_get_result(uint64_t id, ViewResult& out);

    const ViewerRenderConfig& config() const;

    // Viewable buffer keys for this config (the Python /buffers set; the
    // distortion buffers appear only when a distortion regularizer is
    // configured -- they cost full-resolution VRAM to produce).
    std::vector<std::string> buffer_keys() const;

private:
    struct Impl;
    std::unique_ptr<Impl> _impl;
};

// One-shot upload of the post-split camera table for frustum annotation +
// thumbnails (annotation.ensure_viewer_initialized port: engine_viewer_init
// + the camera_size kNN heuristic). Call after engine_setup_data_manager,
// before rendering with show_cams. Returns the kNN-based camera size, for
// ViewerRenderConfig::base_camera_size.
float viewer_upload_cameras(const PostSplitCameras& post);

// The kNN camera-size heuristic on its own (0.2 * median 4-NN camera
// distance), for renderers that draw frusta without the engine (the GUI's
// dataset preview).
float viewer_camera_size_heuristic(const PostSplitCameras& post);

// Unit ray direction (CV convention: x right, y down, z forward) for a
// viewer pixel, u = (px - cx)/fx, v = (py - cy)/fy. Host port of
// projection_utils generate_ray for the distortion-free display models
// 0 pinhole / 1 fisheye-equidistant / 2 equisolid / 3 equirectangular;
// false when the pixel is outside the model's domain.
bool viewer_pixel_ray(int camera_model, float u, float v, float dir[3]);

// One-shot upload of the axes/grid overlay (engine_viewer_set_grid). The
// grid is axis-aligned in the engine's training frame -- the frame splats
// are saved in -- so grid lines mark round coordinates of the exported
// model; its cell size then adapts to ViewRequest::grid_dist per render.
// Call alongside viewer_upload_cameras; drawn when ViewRequest::show_grid
// is set.
void viewer_upload_grid(const PostSplitCameras& post);
