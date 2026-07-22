// Python bindings for the native web-viewer server (app/webviewer/) and the
// post-split camera bake (data/DatasetParser.h).
//
// This is the C++ half of docs/restructure-proposal.md §4.2: Python stops
// carrying a second HTTP server + render worker and drives the same one the
// CLI trainer and the GUI use. viewer.html is already single-source (embedded
// at build time), so the browser client is bit-identical either way.
//
// GIL discipline
// --------------
// The render worker thread must never call into Python -- doing so would
// deadlock against a training loop that holds the engine mutex while holding
// the GIL. So the hooks handed to ViewerServer read plain state (atomics and a
// mutex-guarded string) that Python *pushes* into this object, rather than
// calling Python back. The only shared lock is `engine_mutex`, and Python
// acquires it through `engine_lock()`, which releases the GIL while blocking.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "app/webviewer/Viewer.h"
#include "app/webviewer/RenderWorker.h"
#include "data/DatasetParser.h"

#include <atomic>
#include <memory>
#include <mutex>
#include <string>

namespace py = pybind11;

namespace {

template <typename T>
py::array_t<T> as_array(const std::vector<T>& v, std::vector<py::ssize_t> shape) {
    py::ssize_t total = 1;
    for (auto s : shape) total *= s;
    if (total != (py::ssize_t)v.size())
        throw std::runtime_error("bind_viewer: shape does not match vector size");
    py::array_t<T> out(shape);
    if (!v.empty())
        std::memcpy(out.mutable_data(), v.data(), v.size() * sizeof(T));
    return out;
}

template <typename T>
py::array_t<T> as_array(const std::vector<T>& v) {
    return as_array(v, {(py::ssize_t)v.size()});
}

// ---------------------------------------------------------------------------
// WebViewer -- ViewerServer plus the state its hooks read.
// ---------------------------------------------------------------------------
class WebViewer {
public:
    WebViewer() = default;
    ~WebViewer() { stop(); }

    WebViewer(const WebViewer&) = delete;
    WebViewer& operator=(const WebViewer&) = delete;

    void start(const std::string& host, int port,
               ViewerRenderConfig cfg, const PostSplitCameras& post) {
        if (_started) throw std::runtime_error("WebViewer already started");

        ViewerHooks hooks;
        hooks.engine_mutex = &_engine_mutex;
        hooks.current_step = [this] { return _step.load(); };
        hooks.progress_json = [this] {
            std::lock_guard<std::mutex> lk(_progress_mu);
            return _progress_json;
        };
        hooks.pause_toggle = [this] {
            bool now = !_paused.load();
            _paused.store(now);
            return now;
        };
        hooks.set_render_pending = [this](bool v) { _render_pending.store(v); };

        _server.start(host, port, std::move(cfg), std::move(hooks), post);
        _started = true;
    }

    void stop() {
        if (!_started) return;
        _server.stop();
        _started = false;
    }

    bool started() const { return _started; }

    void set_step(int s) { _step.store(s); }
    int  step() const    { return _step.load(); }

    void set_progress_json(const std::string& j) {
        std::lock_guard<std::mutex> lk(_progress_mu);
        _progress_json = j;
    }

    bool paused() const        { return _paused.load(); }
    void set_paused(bool v)    { _paused.store(v); }
    bool render_pending() const { return _render_pending.load(); }

    std::mutex& engine_mutex() { return _engine_mutex; }

private:
    ViewerServer _server;
    bool _started = false;

    std::mutex _engine_mutex;
    std::atomic<int>  _step{0};
    std::atomic<bool> _paused{false};
    std::atomic<bool> _render_pending{false};

    std::mutex  _progress_mu;
    std::string _progress_json = "{}";
};

// Context manager over WebViewer::engine_mutex. Releasing the GIL while
// blocking is what keeps a render (worker thread, holds the mutex, no GIL)
// from deadlocking against the training loop (holds the GIL, wants the mutex).
class EngineLock {
public:
    explicit EngineLock(WebViewer& v) : _v(v) {}

    void enter() {
        py::gil_scoped_release release;
        _v.engine_mutex().lock();
        _held = true;
    }

    void exit(const py::object&, const py::object&, const py::object&) {
        if (!_held) return;
        _held = false;
        _v.engine_mutex().unlock();
    }

private:
    WebViewer& _v;
    bool _held = false;
};

}  // namespace

void bind_viewer(py::module_& m) {
    // -----------------------------------------------------------------
    // PostSplitCameras + bake_post_split
    // -----------------------------------------------------------------
    // The cubemap split for fisheye / equisolid / equirectangular inputs.
    // This is the C++ original that trainer.py::_setup_cpp_data_manager was
    // ported from; the arrays feed engine_setup_data_manager directly.
    py::class_<PostSplitCameras>(m, "PostSplitCameras")
        .def_readonly("n_post", &PostSplitCameras::n_post)
        .def_readonly("any_warp", &PostSplitCameras::any_warp)
        .def_readonly("any_fisheye_warp", &PostSplitCameras::any_fisheye_warp)
        .def_readonly("direct_equirect", &PostSplitCameras::direct_equirect)
        // Kept as plain lists: engine_setup_data_manager takes std::vector.
        .def_readonly("K_per_camera", &PostSplitCameras::K_per_camera)
        .def_readonly("post_offsets", &PostSplitCameras::post_offsets)
        .def_readonly("viewmats", &PostSplitCameras::viewmats)
        .def_readonly("intrins", &PostSplitCameras::intrins)
        .def_readonly("dist_coeffs", &PostSplitCameras::dist_coeffs)
        .def_readonly("input_intrins", &PostSplitCameras::input_intrins)
        .def_readonly("input_dist_coeffs", &PostSplitCameras::input_dist_coeffs)
        .def_readonly("post_widths", &PostSplitCameras::post_widths)
        .def_readonly("post_heights", &PostSplitCameras::post_heights)
        .def_readonly("post_models", &PostSplitCameras::post_models)
        // Shaped views for the arrays callers actually inspect.
        .def_property_readonly("viewmats_array", [](const PostSplitCameras& p) {
            return as_array(p.viewmats, {p.n_post, 4, 4});
        })
        .def_property_readonly("intrins_array", [](const PostSplitCameras& p) {
            return as_array(p.intrins, {p.n_post, 4});
        })
        .def_property_readonly("c2w_flip", [](const PostSplitCameras& p) {
            return as_array(p.c2w_flip, {p.n_post, 3, 4});
        })
        .def("__repr__", [](const PostSplitCameras& p) {
            return "<PostSplitCameras n_post=" + std::to_string(p.n_post) +
                   (p.any_warp ? " warped>" : ">");
        });

    m.def("bake_post_split", &bake_post_split,
          py::arg("dataset"), py::arg("warp_to_pinhole") = false,
          py::arg("warp_spherical_to_pinhole") = false,
          "Expand a ParsedDataset into the POST-split camera arrays "
          "engine_setup_data_manager consumes (fisheye/equisolid -> 5 cubemap "
          "faces, equirectangular -> 6, everything else K=1).");

    // -----------------------------------------------------------------
    // ViewerRenderConfig
    // -----------------------------------------------------------------
    py::class_<ViewerRenderConfig>(m, "ViewerRenderConfig")
        .def(py::init<>())
        .def_readwrite("primitive", &ViewerRenderConfig::primitive)
        .def_readwrite("packed", &ViewerRenderConfig::packed)
        .def_readwrite("sh_degree_warmup_every",
                       &ViewerRenderConfig::sh_degree_warmup_every)
        .def_readwrite("relative_scale", &ViewerRenderConfig::relative_scale)
        .def_readwrite("output_median", &ViewerRenderConfig::output_median)
        .def_readwrite("distortion_reg_on", &ViewerRenderConfig::distortion_reg_on)
        .def_readwrite("train_frame_scale", &ViewerRenderConfig::train_frame_scale)
        .def_property("train_to_normalized",
            [](const ViewerRenderConfig& c) {
                py::array_t<float> out({(py::ssize_t)4, (py::ssize_t)4});
                std::memcpy(out.mutable_data(), c.train_to_normalized.data(),
                            16 * sizeof(float));
                return out;
            },
            [](ViewerRenderConfig& c, py::array_t<float> v) {
                if (v.size() != 16)
                    throw std::runtime_error("train_to_normalized needs 16 floats");
                auto buf = v.template unchecked<>();
                (void)buf;
                std::memcpy(c.train_to_normalized.data(), v.data(),
                            16 * sizeof(float));
            });

    // -----------------------------------------------------------------
    // WebViewer
    // -----------------------------------------------------------------
    py::class_<EngineLock>(m, "_EngineLock")
        .def("__enter__", &EngineLock::enter)
        .def("__exit__", &EngineLock::exit);

    py::class_<WebViewer>(m, "WebViewer", R"doc(
The native web-viewer server (same HTTP endpoints and same viewer.html as
ssplat-train and ssplat-gui).

Wrap every engine call in the training loop with ``engine_lock()`` so the
render worker gets a turn:

    viewer = WebViewer()
    viewer.start("0.0.0.0", 7007, cfg, post)
    for step in range(n):
        viewer.set_step(step)
        viewer.set_progress_json(json.dumps({...}))
        with viewer.engine_lock():
            engine_train_step_managed(...)
    viewer.stop()

``render_pending`` is set by the worker when a browser request is waiting; a
loop that does long uninterrupted work should yield on it.
)doc")
        .def(py::init<>())
        .def("start",
             [](WebViewer& v, const std::string& host, int port,
                ViewerRenderConfig cfg, const PostSplitCameras& post) {
                 py::gil_scoped_release release;
                 v.start(host, port, std::move(cfg), post);
             },
             py::arg("host") = "0.0.0.0", py::arg("port") = 7007,
             py::arg("config") = ViewerRenderConfig(),
             py::arg("post") = PostSplitCameras())
        .def("stop", [](WebViewer& v) {
                 py::gil_scoped_release release;
                 v.stop();
             })
        .def_property_readonly("started", &WebViewer::started)
        .def("set_step", &WebViewer::set_step, py::arg("step"))
        .def_property_readonly("step", &WebViewer::step)
        .def("set_progress_json", &WebViewer::set_progress_json, py::arg("json"))
        .def_property("paused", &WebViewer::paused, &WebViewer::set_paused)
        .def_property_readonly("render_pending", &WebViewer::render_pending)
        .def("engine_lock", [](WebViewer& v) { return EngineLock(v); },
             py::keep_alive<0, 1>(),
             "Context manager over the engine mutex shared with the render "
             "worker. Releases the GIL while blocking.");
}
