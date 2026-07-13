// Viewer.cpp -- see Viewer.h. HTTP endpoint handlers are ports of
// viewer/http_server.py; the render path itself lives in RenderWorker.cpp.

#include "Viewer.h"
#include "HttpServer.h"

// Implementation lives in csrc/stb_image_write_impl.cpp (shared with the
// mesh texture export); declarations only here.
#include "../external/stb_image_write.h"

#include "app_generated/viewer_html.h"   // kViewerHtml[] (CMake-embedded)

#include "../Camera.h"    // camera_model_from_name

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <vector>


namespace {

void jpeg_write_cb(void* ctx, void* data, int size) {
    auto* out = (std::vector<uint8_t>*)ctx;
    out->insert(out->end(), (uint8_t*)data, (uint8_t*)data + size);
}

}  // namespace


// ===========================================================================
// Impl
// ===========================================================================

struct ViewerServer::Impl {
    ViewerHooks hooks;
    HttpServer http;
    RenderWorker worker;
    std::vector<std::string> buffer_keys;
    bool started = false;

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
        ViewRequest q;
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
        int quality = req.get_int("jpeg_quality", 75);
        q.show_cams = req.get_bool("show_training_cameras", false);
        q.show_grid = req.get_bool("show_grid", false);
        q.grid_dist = (float)req.get_double("grid_dist", 0.0);
        q.grid_target[0] = (float)req.get_double("grid_tx", 0.0);
        q.grid_target[1] = (float)req.get_double("grid_ty", 0.0);
        q.grid_target[2] = (float)req.get_double("grid_tz", 0.0);
        q.cam_size_scale = (float)req.get_double("camera_size_scale", 1.0);

        constexpr int MAX_DIM = 2160;   // prevent OOM (http_server.py:74)
        if (q.W <= 0 || q.H <= 0 || std::max(q.W, q.H) > MAX_DIM)
            return HttpResponse::text(400, "image too large");
        if (std::find(buffer_keys.begin(), buffer_keys.end(), q.key) == buffer_keys.end())
            return HttpResponse::text(400, "unsupported buffer_key");
        if ((int)camera_model_from_name(q.model) < 0)
            return HttpResponse::text(400, "unsupported camera_model");

        uint64_t id = worker.submit(q);
        ViewResult res;
        if (!worker.wait_result(id, res, 10.0))
            return HttpResponse::text(500, "render timeout");
        if (!res.error.empty())
            return HttpResponse::text(500, res.error);

        // JPEG encode outside the engine lock.
        std::vector<uint8_t> jpeg;
        jpeg.reserve((size_t)res.W * res.H / 4);
        if (!stbi_write_jpg_to_func(jpeg_write_cb, &jpeg, res.W, res.H, 3,
                                    res.rgb8.data(), quality))
            return HttpResponse::text(500, "JPEG encode failed");
        HttpResponse r;
        r.content_type = "image/jpeg";
        r.body = std::move(jpeg);
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
    im.hooks = hooks;

    // One-shot upload of the post-split camera table for frustum annotation
    // + thumbnails (annotation.ensure_viewer_initialized). Host views are
    // fine -- the engine copies.
    cfg.base_camera_size = viewer_upload_cameras(post);
    viewer_upload_grid(post);

    im.worker.start(std::move(cfg), std::move(hooks));
    im.buffer_keys = im.worker.buffer_keys();
    im.started = true;

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
    if (!im.started) return;
    im.started = false;
    im.worker.stop();
    im.http.stop();
}
