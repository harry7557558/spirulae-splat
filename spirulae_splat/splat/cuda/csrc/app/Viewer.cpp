// Viewer.cpp -- see Viewer.h. Ports:
//   viewer/http_server.py   endpoint handlers (route bodies below)
//   viewer/render_worker.py latest-wins pending slot + result handoff
//   trainer.py _render/render  c2w remap + engine render + annotation call
//   model.py get_outputs    (viewer subset: rgb/depth/alpha/normals/median/
//                            distortion display transforms)
//   viewer/annotation.py    engine_viewer_init upload + engine_blit_view call
//                           + the camera_size kNN heuristic

#include "Viewer.h"
#include "HttpServer.h"

#include "../Engine.h"    // engine render + viewer entry points
#include "../Camera.h"    // camera_model_from_name

#include <cuda_runtime.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#define STBI_WRITE_NO_STDIO
#include "../external/stb_image_write.h"

#include "app_generated/viewer_html.h"   // kViewerHtml[] (CMake-embedded)

#include <algorithm>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <thread>
#include <unordered_set>
#include <vector>


namespace {

// ---------------------------------------------------------------------------
// Small helpers
// ---------------------------------------------------------------------------

TorchTensorView tvp(const void* p, uint32_t elem_size, std::vector<int64_t> shape) {
    return {(uint64_t)(uintptr_t)p, elem_size, std::move(shape)};
}
TorchTensorView tv_null() { return {0, 0, {}}; }

// Grow-only device buffer (viewer scratch; a handful per worker, so no pool).
struct DevBuf {
    void* ptr = nullptr;
    size_t cap = 0;
    void* ensure(size_t bytes) {
        if (bytes > cap) {
            if (ptr) cudaFree(ptr);
            if (cudaMalloc(&ptr, bytes) != cudaSuccess) {
                ptr = nullptr; cap = 0;
                throw std::runtime_error("viewer: cudaMalloc failed");
            }
            cap = bytes;
        }
        return ptr;
    }
    void* upload(const void* host, size_t bytes) {
        void* p = ensure(bytes);
        if (cudaMemcpy(p, host, bytes, cudaMemcpyHostToDevice) != cudaSuccess)
            throw std::runtime_error("viewer: H2D copy failed");
        return p;
    }
    ~DevBuf() { if (ptr) cudaFree(ptr); }
};

// Median distance to the 4th-nearest unique camera position, for the frustum
// scale (port of annotation.py _knn_dist; O(N^2) is fine at camera counts).
double camera_knn_dist(const std::vector<float>& c2w_flip, int64_t n, int k = 4) {
    if (n <= 1) return 1.0;
    std::vector<std::array<double, 3>> pos;
    pos.reserve(n);
    double scale = 1e-9;
    for (int64_t i = 0; i < n; i++)
        for (int r = 0; r < 3; r++)
            scale = std::max(scale, (double)std::abs(c2w_flip[i*12 + r*4 + 3]));
    std::unordered_set<uint64_t> seen;   // hash of quantized position
    for (int64_t i = 0; i < n; i++) {
        std::array<double, 3> p;
        uint64_t h = 1469598103934665603ull;
        for (int r = 0; r < 3; r++) {
            p[r] = c2w_flip[i*12 + r*4 + 3];
            int64_t q = (int64_t)std::llround(p[r] / scale * 1e6);
            h = (h ^ (uint64_t)q) * 1099511628211ull;
        }
        if (seen.insert(h).second) pos.push_back(p);
    }
    int64_t m = (int64_t)pos.size();
    if (m <= 1) return 1.0;
    std::vector<double> kth;
    std::vector<double> d2(m);
    for (int64_t i = 0; i < m; i++) {
        for (int64_t j = 0; j < m; j++) {
            double dx = pos[j][0]-pos[i][0], dy = pos[j][1]-pos[i][1], dz = pos[j][2]-pos[i][2];
            d2[j] = dx*dx + dy*dy + dz*dz;
        }
        int kk = (int)std::min<int64_t>(k, m - 1);
        std::nth_element(d2.begin(), d2.begin() + kk, d2.end());
        kth.push_back(std::sqrt(d2[kk]));   // max over the k nearest (self = 0)
    }
    std::nth_element(kth.begin(), kth.begin() + kth.size()/2, kth.end());
    return kth[kth.size()/2];
}

void jpeg_write_cb(void* ctx, void* data, int size) {
    auto* out = (std::vector<uint8_t>*)ctx;
    out->insert(out->end(), (uint8_t*)data, (uint8_t*)data + size);
}

std::string json_escape(const std::string& s) {
    std::string o;
    for (char c : s) { if (c == '"' || c == '\\') o += '\\'; o += c; }
    return o;
}

struct RenderReq {
    uint64_t id = 0;
    float c2w[12];
    float fx, fy, cx, cy;
    int W = 512, H = 512;
    std::string model = "PINHOLE";
    std::string key = "rgb";
    bool show_cams = false;
    int quality = 75;
};

struct RenderRes {
    uint64_t id = 0;
    std::vector<uint8_t> jpeg;
    std::string error;
};

}  // namespace


// ===========================================================================
// Impl
// ===========================================================================

struct ViewerServer::Impl {
    ViewerRenderConfig cfg;
    ViewerHooks hooks;
    HttpServer http;

    std::thread worker;
    std::atomic<bool> running{false};
    std::mutex mu;
    std::condition_variable cv_work, cv_result;
    bool has_pending = false;
    RenderReq pending;
    RenderRes result;
    uint64_t next_id = 1;

    std::vector<std::string> buffer_keys;

    // Device scratch (worker thread only).
    DevBuf d_buffer, d_depth, d_alpha, d_out;
    DevBuf d_intr, d_vm, d_dist, d_normals, d_depth_in;

    // ---- render worker -----------------------------------------------------

    uint64_t submit(RenderReq q) {
        std::unique_lock<std::mutex> lk(mu);
        q.id = next_id++;
        pending = q;
        has_pending = true;
        uint64_t id = q.id;
        lk.unlock();
        // Synchronously flip the "render desired" flag so the training loop
        // yields at its next iteration boundary (render_worker.py:44-53).
        if (hooks.set_render_pending) hooks.set_render_pending(true);
        cv_work.notify_one();
        return id;
    }

    bool wait_result(uint64_t id, RenderRes& out, double timeout_s) {
        std::unique_lock<std::mutex> lk(mu);
        bool ok = cv_result.wait_for(lk, std::chrono::duration<double>(timeout_s),
                                     [&] { return result.id == id; });
        if (ok) out = result;
        return ok;
    }

    void worker_loop() {
        while (running) {
            RenderReq q;
            {
                std::unique_lock<std::mutex> lk(mu);
                cv_work.wait(lk, [&] { return has_pending || !running; });
                if (!running) break;
                q = pending;
                has_pending = false;
            }
            RenderRes res;
            res.id = q.id;
            try {
                res.jpeg = render_once(q);
            } catch (const std::exception& e) {
                std::fprintf(stderr, "[viewer] render error: %s\n", e.what());
                res.error = e.what();
            }
            {
                std::lock_guard<std::mutex> lk(mu);
                result = std::move(res);
            }
            cv_result.notify_all();
            {
                std::lock_guard<std::mutex> lk(mu);
                if (!has_pending && hooks.set_render_pending)
                    hooks.set_render_pending(false);
            }
        }
    }

    // ---- the actual render (trainer._render + get_outputs viewer subset) ---

    std::vector<uint8_t> render_once(const RenderReq& q) {
        const int W = q.W, H = q.H;
        const int64_t npx = (int64_t)W * H;
        static const float D[3] = {1.f, -1.f, -1.f};

        // Client c2w is in the legacy normalized frame; remap to the training
        // frame (trainer.py:557-576). T is a similarity (scale*R | t): full
        // similarity on the position, unit rotation on the basis.
        float c2w[12];
        std::memcpy(c2w, q.c2w, sizeof c2w);
        if (cfg.train_frame_scale != 1.0f) {
            const auto& T = cfg.train_to_normalized;
            double s = std::sqrt((double)T[0]*T[0] + (double)T[4]*T[4] + (double)T[8]*T[8]);
            float out[12];
            for (int r = 0; r < 3; r++) {
                for (int c = 0; c < 3; c++) {
                    double v = 0.0;
                    for (int m = 0; m < 3; m++)
                        v += (double)T[r*4 + m] / s * q.c2w[m*4 + c];
                    out[r*4 + c] = (float)v;
                }
                double t = T[r*4 + 3];
                for (int m = 0; m < 3; m++) t += (double)T[r*4 + m] * q.c2w[m*4 + 3];
                out[r*4 + 3] = (float)t;
            }
            std::memcpy(c2w, out, sizeof c2w);
        }

        // c2w -> engine viewmat (model.py get_outputs:1016-1025).
        double t[3] = {c2w[3], c2w[7], c2w[11]};
        if (cfg.relative_scale.has_value())
            for (auto& v : t) v *= *cfg.relative_scale;
        double Rf[3][3];
        for (int r = 0; r < 3; r++)
            for (int c = 0; c < 3; c++)
                Rf[r][c] = c2w[r*4 + c] * D[c];
        float vm[16] = {0};
        for (int r = 0; r < 3; r++) {
            double ti = 0.0;
            for (int c = 0; c < 3; c++) {
                vm[r*4 + c] = (float)Rf[c][r];
                ti -= Rf[c][r] * t[c];
            }
            vm[r*4 + 3] = (float)ti;
        }
        vm[15] = 1.f;

        float intr[4] = {q.fx, q.fy, q.cx, q.cy};
        float dist0[10] = {0};
        bool want_dist = cfg.distortion_reg_on ||
                         q.key.find("distortion") != std::string::npos;
        bool want_median = cfg.output_median;
        int sh_deg = hooks.current_step
            ? hooks.current_step() / std::max(cfg.sh_degree_warmup_every, 1) : 100;

        std::vector<float> rgb(npx * 3), depth(npx), Ts(npx);
        std::vector<float> median(want_median ? npx : 0);
        std::vector<float> rgbd, depthd;

        std::vector<uint8_t> out_host(npx * 3);
        {
            // Engine access: exclusive with the training step. Everything up
            // to the annotated-u8 download stays under the lock, mirroring
            // Trainer.render's "post-processor + D->H inside the lock".
            std::lock_guard<std::mutex> lk(*hooks.engine_mutex);

            set_camera_params(W, H, q.model,
                              tvp(vm, 4, {1, 4, 4}),
                              tvp(intr, 4, {1, 4}),
                              tvp(dist0, 4, {1, 10}));
            forward_3dgs(cfg.primitive, sh_deg, cfg.packed,
                         want_median, want_dist ? 2 : 0);
            engine_copy_render_to_host(
                tvp(rgb.data(), 4, {1, H, W, 3}),
                tvp(depth.data(), 4, {1, H, W, 1}),
                tvp(Ts.data(), 4, {1, H, W, 1}),
                tv_null(),
                want_median ? tvp(median.data(), 4, {1, H, W, 1}) : tv_null());
            if (want_dist) {
                rgbd.resize(npx * 3);
                depthd.resize(npx);
                engine_copy_distortion_to_host(
                    tvp(rgbd.data(), 4, {1, H, W, 3}),
                    tvp(depthd.data(), 4, {1, H, W, 1}));
            }

            // Display transforms (model.py get_outputs:1223-1254).
            float inv_rs = cfg.relative_scale.has_value() ? 1.0f / *cfg.relative_scale : 1.0f;
            for (auto& v : depth) v *= inv_rs;
            if (want_median) for (auto& v : median) v *= inv_rs;
            auto dist_display = [](std::vector<float>& v) {
                constexpr float e = 1.0f / 255.0f;
                for (auto& x : v) x = std::sqrt(x + e * e) - e;
            };
            if (want_dist) { dist_display(rgbd); dist_display(depthd); }

            std::vector<float> alpha3(npx * 3);
            for (int64_t i = 0; i < npx; i++) {
                float a = 1.0f - Ts[i];
                alpha3[i*3] = alpha3[i*3+1] = alpha3[i*3+2] = a;
            }

            // depth -> normal display (0.5 + 0.5n), via the engine kernel on
            // device copies (model.py:1141-1161).
            auto depth_normal_display = [&](const std::vector<float>& d_host) {
                float* dd = (float*)d_depth_in.upload(d_host.data(), npx * 4);
                float* dn = (float*)d_normals.ensure(npx * 12);
                float* di = (float*)d_intr.upload(intr, sizeof intr);
                float* dc = (float*)d_dist.upload(dist0, sizeof dist0);
                depth_to_normal_forward_tv(
                    q.model, tvp(di, 4, {1, 4}), tvp(dc, 4, {1, 10}),
                    /*is_ray_depth=*/true,
                    tvp(dd, 4, {1, H, W, 1}), tvp(dn, 4, {1, H, W, 3}));
                std::vector<float> n_host(npx * 3);
                if (cudaMemcpy(n_host.data(), dn, npx * 12,
                               cudaMemcpyDeviceToHost) != cudaSuccess)
                    throw std::runtime_error("viewer: normal D2H failed");
                for (auto& v : n_host) v = 0.5f + 0.5f * v;
                return n_host;
            };

            // Select the display buffer for the requested key.
            std::vector<float> local;
            const std::vector<float>* buf = nullptr;
            int C = 3;
            if (q.key == "rgb")               { buf = &rgb; }
            else if (q.key == "depth")        { buf = &depth; C = 1; }
            else if (q.key == "alpha")        { buf = &alpha3; }
            else if (q.key == "depth_normal") { local = depth_normal_display(depth); buf = &local; }
            else if (q.key == "depth_median" && want_median)
                { buf = &median; C = 1; }
            else if (q.key == "normal_median" && want_median)
                { local = depth_normal_display(median); buf = &local; }
            else if (q.key == "rgb_distortion" && want_dist)   { buf = &rgbd; }
            else if (q.key == "depth_distortion" && want_dist) { buf = &depthd; C = 1; }
            else throw std::runtime_error("unsupported buffer_key: " + q.key);

            // Annotate + colormap on GPU (annotation.annotate_train_cameras).
            float* db = (float*)d_buffer.upload(buf->data(), npx * C * 4);
            float* dz = (float*)d_depth.upload(depth.data(), npx * 4);
            float* da = (float*)d_alpha.upload(alpha3.data(), npx * 12);
            float* di = (float*)d_intr.upload(intr, sizeof intr);
            float* dv = (float*)d_vm.upload(vm, sizeof vm);
            float* dc = (float*)d_dist.upload(dist0, sizeof dist0);
            uint8_t* dout = (uint8_t*)d_out.ensure(npx * 3);
            engine_blit_view(
                q.key,
                tvp(db, 4, {H, W, C}),
                tvp(dz, 4, {H, W, 1}),
                tvp(da, 4, {H, W, 3}),
                (int)camera_model_from_name(q.model),
                tvp(di, 4, {1, 4}),
                tvp(dv, 4, {4, 4}),
                tvp(dc, 4, {1, 10}),
                q.show_cams,
                tvp(dout, 1, {H, W, 3}));
            if (cudaMemcpy(out_host.data(), dout, npx * 3,
                           cudaMemcpyDeviceToHost) != cudaSuccess)
                throw std::runtime_error("viewer: blit D2H failed");
        }

        // JPEG encode outside the engine lock.
        std::vector<uint8_t> jpeg;
        jpeg.reserve(npx / 4);
        if (!stbi_write_jpg_to_func(jpeg_write_cb, &jpeg, W, H, 3,
                                    out_host.data(), q.quality))
            throw std::runtime_error("viewer: JPEG encode failed");
        return jpeg;
    }

    // ---- endpoint handlers (http_server.py port) ----------------------------

    HttpResponse handle_index() {
        // Dev override so viewer.html edits don't need a rebuild.
        if (const char* p = std::getenv("SSPLAT_VIEWER_HTML")) {
            std::ifstream f(p, std::ios::binary);
            if (f) {
                HttpResponse r;
                r.content_type = "text/html; charset=utf-8";
                r.body.assign(std::istreambuf_iterator<char>(f), {});
                return r;
            }
            std::fprintf(stderr, "[viewer] SSPLAT_VIEWER_HTML=%s not readable; "
                         "serving embedded copy\n", p);
        }
        HttpResponse r;
        r.content_type = "text/html; charset=utf-8";
        r.body.assign(kViewerHtml, kViewerHtml + kViewerHtmlSize);
        return r;
    }

    HttpResponse handle_render(const HttpRequest& req) {
        RenderReq q;
        q.key = req.get("buffer_key", "rgb");
        std::string c2w_str = req.get("c2w");
        {
            size_t pos = 0;
            int i = 0;
            while (i < 12 && pos <= c2w_str.size()) {
                size_t comma = c2w_str.find(',', pos);
                try {
                    q.c2w[i++] = std::stof(c2w_str.substr(
                        pos, comma == std::string::npos ? std::string::npos : comma - pos));
                } catch (...) { return HttpResponse::text(400, "bad c2w"); }
                if (comma == std::string::npos) break;
                pos = comma + 1;
            }
            if (i != 12) return HttpResponse::text(400, "c2w must have 12 values");
        }
        q.fx = (float)req.get_double("fx", 500);
        q.fy = (float)req.get_double("fy", 500);
        q.cx = (float)req.get_double("cx", 256);
        q.cy = (float)req.get_double("cy", 256);
        q.W  = req.get_int("width", 512);
        q.H  = req.get_int("height", 512);
        q.model = req.get("camera_model", "PINHOLE");
        q.quality = req.get_int("jpeg_quality", 75);
        q.show_cams = req.get_bool("show_training_cameras", false);

        constexpr int MAX_DIM = 2160;   // prevent OOM (http_server.py:74)
        if (q.W <= 0 || q.H <= 0 || std::max(q.W, q.H) > MAX_DIM)
            return HttpResponse::text(400, "image too large");
        if (std::find(buffer_keys.begin(), buffer_keys.end(), q.key) == buffer_keys.end())
            return HttpResponse::text(400, "unsupported buffer_key");
        if ((int)camera_model_from_name(q.model) < 0)
            return HttpResponse::text(400, "unsupported camera_model");

        uint64_t id = submit(q);
        RenderRes res;
        if (!wait_result(id, res, 10.0))
            return HttpResponse::text(500, "render timeout");
        if (!res.error.empty())
            return HttpResponse::text(500, res.error);
        HttpResponse r;
        r.content_type = "image/jpeg";
        r.body = std::move(res.jpeg);
        return r;
    }

    HttpResponse handle_buffers() {
        std::string body = "[";
        for (size_t i = 0; i < buffer_keys.size(); i++)
            body += (i ? ", \"" : "\"") + buffer_keys[i] + "\"";
        body += "]";
        return HttpResponse::json(body);
    }
};


// ===========================================================================
// ViewerServer
// ===========================================================================

ViewerServer::ViewerServer() : _impl(new Impl) {}
ViewerServer::~ViewerServer() { stop(); }

void ViewerServer::start(const std::string& host, int port,
                         ViewerRenderConfig cfg, ViewerHooks hooks,
                         const PostSplitCameras& post) {
    Impl& im = *_impl;
    im.cfg = std::move(cfg);
    im.hooks = std::move(hooks);

    // Viewable buffers (the Python /buffers set minus the debug-only
    // sh / refinement_score renders, which the CLI doesn't produce).
    im.buffer_keys = {"rgb", "depth", "alpha", "depth_normal"};
    if (im.cfg.output_median) {
        im.buffer_keys.push_back("depth_median");
        im.buffer_keys.push_back("normal_median");
    }
    im.buffer_keys.push_back("rgb_distortion");
    im.buffer_keys.push_back("depth_distortion");

    // One-shot upload of the post-split camera table for frustum annotation
    // + thumbnails (annotation.ensure_viewer_initialized). Host views are
    // fine -- the engine copies.
    float camera_size = 0.2f * (float)camera_knn_dist(post.c2w_flip, post.n_post);
    engine_viewer_init(
        tvp(post.post_models.data(), 4, {post.n_post}),
        tvp(post.intrins.data(), 4, {post.n_post, 4}),
        tvp(post.dist_coeffs.data(), 4, {post.n_post, 10}),
        tvp(post.c2w_flip.data(), 4, {post.n_post, 3, 4}),
        tvp(post.post_widths.data(), 4, {post.n_post}),
        tvp(post.post_heights.data(), 4, {post.n_post}),
        camera_size);

    im.running = true;
    im.worker = std::thread([&im] { im.worker_loop(); });

    im.http.route("/",           [&im](const HttpRequest&) { return im.handle_index(); });
    im.http.route("/index.html", [&im](const HttpRequest&) { return im.handle_index(); });
    im.http.route("/render",     [&im](const HttpRequest& r) { return im.handle_render(r); });
    im.http.route("/buffers",    [&im](const HttpRequest&) { return im.handle_buffers(); });
    im.http.route("/progress",   [&im](const HttpRequest&) {
        return HttpResponse::json(im.hooks.progress_json ? im.hooks.progress_json()
                                                         : "{}");
    });
    im.http.route("/pause-toggle", [&im](const HttpRequest&) {
        if (!im.hooks.pause_toggle)
            return HttpResponse::json(
                "{\"paused\": false, \"error\": \"pause_toggle not available\"}");
        bool paused = im.hooks.pause_toggle();
        return HttpResponse::json(std::string("{\"paused\": ") +
                                  (paused ? "true" : "false") + "}");
    });
    im.http.start(host, port);
}

void ViewerServer::stop() {
    Impl& im = *_impl;
    if (!im.running.exchange(false)) return;
    im.cv_work.notify_all();
    if (im.worker.joinable()) im.worker.join();
    im.http.stop();
}
