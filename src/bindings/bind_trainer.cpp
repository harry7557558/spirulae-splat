// Python bindings for the native training session (app/TrainerCore.{h,cpp}).
//
// This is the C++ half of docs/restructure-proposal.md §4.3, the last of the
// three duplicated subsystems. TrainerCore.h's own header comment describes
// itself as a line-by-line port of the Python managed path -- a comment doing
// the job a shared implementation should do. Binding it makes the C++ side
// the single implementation and lets `model.py::engine_train_step_managed`
// plus `core.py::_build_{optim,densify}_config` go away.
//
// What is bound, and why each piece:
//
//   SsplatConfig          the flattened training config. Generated from the
//                         Python dataclasses by tools/codegen/generate_cli_config.py,
//                         so this binding is written against the generated
//                         X-macro rather than a hand-listed field set -- add a
//                         field to a Python dataclass, re-run the generator,
//                         and it appears here with no edit to this file.
//   ssplat_config_fields() the same X-macro as data, so Python can walk its
//                         own dataclasses and fill SsplatConfig without a
//                         second hand-maintained name mapping.
//   build_step_config()   the per-step EngineStepConfig builder. The actual
//                         drift bomb, and what the parity test compares.
//   TrainerSession        the phased session (check_config / load_dataset /
//                         setup_engine / train), plus train_step() for a
//                         front-end that keeps its own loop.
//
// GIL discipline matches bind_viewer.cpp: every call that blocks or takes the
// engine mutex releases the GIL, and the C++ side never calls into Python
// except through the explicit on_step callback (which re-acquires).

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>

#include "app/TrainerCore.h"
#include "bindings/EngineLock.h"

#include <string>
#include <vector>

namespace py = pybind11;
using namespace ssplat;

void bind_trainer(py::module_& m) {
    // -----------------------------------------------------------------
    // SsplatConfig -- generated field set, bound mechanically
    // -----------------------------------------------------------------
    py::class_<SsplatConfig> cfg_cls(m, "SsplatConfig", R"doc(
The flattened native training config (src/app/generated/cli_config.h).

Field names are the CLI keys with '_' separators; defaults are the Python
dataclass defaults, baked in by the generator. Build one from a Python
TrainerConfig with spirulae_splat.modules.native_trainer.to_native_config().
)doc");
    cfg_cls.def(py::init<>());
#define X(member, key, pyname, group, choices, help) \
    cfg_cls.def_readwrite(key, &SsplatConfig::member, help);
    SSPLAT_CONFIG_FIELDS(X)
#undef X

    // The preset overrides ss_trainer.py applies via its tyro subcommands.
    m.def("ssplat_apply_preset",
          [](SsplatConfig& c, const std::string& name) {
              if (!ssplat_apply_preset(c, name))
                  throw std::runtime_error("unknown preset: " + name);
          },
          py::arg("config"), py::arg("name"),
          "Apply a preset's field overrides in place (3dgs, 360-camera, "
          "in-the-wild, linear-color, synthetic, meshing, academic-baseline).");

    // (cli_key, pyname, group, choices, help) per field, straight off the
    // generated X-macro. Python uses (group, pyname) to find the value on its
    // own dataclass tree and cli_key to set it here -- so the name mapping
    // exists once, in the generator, instead of once more in Python.
    m.def("ssplat_config_fields", []() {
        py::list out;
#define X(member, key, pyname, group, choices, help) \
        out.append(py::make_tuple(key, pyname, group, choices, help));
        SSPLAT_CONFIG_FIELDS(X)
#undef X
        return out;
    }, "Field table for SsplatConfig: (cli_key, pyname, group, choices, help).");

    // -----------------------------------------------------------------
    // The ported per-step builders
    // -----------------------------------------------------------------
    py::class_<RunState>(m, "RunState", R"doc(
Setup-time state build_step_config() needs beyond the config itself: the
training-frame scale and which of the optional modules actually got
initialized. TrainerSession.run_state carries the real one.
)doc")
        .def(py::init<>())
        .def_readwrite("train_frame_scale",    &RunState::train_frame_scale)
        .def_readwrite("splat_linear",         &RunState::splat_linear)
        .def_readwrite("bilagrid_rgb_init",    &RunState::bilagrid_rgb_init)
        .def_readwrite("bilagrid_depth_init",  &RunState::bilagrid_depth_init)
        .def_readwrite("bilagrid_normal_init", &RunState::bilagrid_normal_init)
        .def_readwrite("ppisp_init",           &RunState::ppisp_init);

    m.def("build_step_config", &build_step_config,
          py::arg("config"), py::arg("run_state"), py::arg("step"),
          "The per-step EngineStepConfig (loss weights, LR schedules, densify "
          "/ bilagrid / PPISP / background args). Replaces the Python "
          "model.py::engine_train_step_managed config half.");

    m.def("build_loss_weights",
          [](const SsplatConfig& c, int step) {
              auto w = build_loss_weights(c, step);
              return std::vector<float>(w.begin(), w.end());
          },
          py::arg("config"), py::arg("step"),
          "Per-step loss weights, indexed by LossWeightIndex.");

    m.def("scheduled_lr",
          [](int step, int max_steps, float lr,
             std::optional<float> lr_final, std::optional<int> warmup) {
              return scheduled_lr(step, max_steps, lr, lr_final, warmup);
          },
          py::arg("step"), py::arg("max_steps"), py::arg("lr"),
          py::arg("lr_final") = py::none(), py::arg("warmup") = py::none(),
          "OptimizerConfig.get_scheduled_lr, natively.");

    // -----------------------------------------------------------------
    // TrainerProgress
    // -----------------------------------------------------------------
    py::class_<TrainerProgress>(m, "TrainerProgress")
        .def_readonly("step",         &TrainerProgress::step)
        .def_readonly("total_steps",  &TrainerProgress::total_steps)
        .def_readonly("step_latency", &TrainerProgress::step_latency)
        .def_readonly("num_splats",   &TrainerProgress::num_splats)
        .def_readonly("losses",       &TrainerProgress::losses)
        .def("__repr__", [](const TrainerProgress& p) {
            return "<TrainerProgress step=" + std::to_string(p.step) + "/" +
                   std::to_string(p.total_steps) + ">";
        });

    // -----------------------------------------------------------------
    // TrainerSession
    // -----------------------------------------------------------------
    py::class_<TrainerSession>(m, "TrainerSession", R"doc(
The native training session -- the same one ssplat-train and ssplat-gui run.

Phases, in order:

    sess = TrainerSession()
    sess.config = to_native_config(trainer_config)
    sess.front_end_handles_resume = True   # Python does its own resume
    sess.check_config()
    sess.load_dataset()      # parse + POST-split bake, no GPU work
    sess.setup_engine()      # output dir, seeding, engine + DataManager

Then either hand the loop over:

    sess.train(on_step=lambda p: print(p.step, p.losses))

or keep your own, which is what a front-end with resume / eval / profiling
wants:

    for step in range(start, sess.config.num_iterations):
        with sess.engine_lock():
            losses = sess.train_step(step)

The engine is a process-global singleton, so only one session may be live at
a time; setup_engine() calls engine_reset() first.
)doc")
        .def(py::init<>())
        // The config is a value member, so expose it by reference: callers
        // mutate `sess.config.num_iterations` in place.
        .def_readwrite("config", &TrainerSession::cfg,
                       py::return_value_policy::reference_internal)
        .def_readwrite("preset", &TrainerSession::preset)
        .def_readwrite("front_end_handles_resume",
                       &TrainerSession::front_end_handles_resume)
        .def_readwrite("front_end_handles_eval",
                       &TrainerSession::front_end_handles_eval)
        .def_readwrite("out_dir_override", &TrainerSession::out_dir_override)
        .def_readwrite("write_config_json", &TrainerSession::write_config_json)

        // Filled by load_dataset() / setup_engine().
        .def_readonly("dataset",    &TrainerSession::ds)
        .def_readonly("post",       &TrainerSession::post)
        .def_readonly("has_depth",  &TrainerSession::has_depth)
        .def_readonly("has_normal", &TrainerSession::has_normal)
        .def_readonly("run_state",  &TrainerSession::st)
        .def_property_readonly("output_dir", [](const TrainerSession& s) {
            return s.out_dir.string();
        })
        .def_readwrite("viewer_base_camera_size",
                       &TrainerSession::viewer_base_camera_size)

        // Run control. These are atomics on the C++ side, shared with the
        // viewer render workers.
        .def_property("paused",
            [](const TrainerSession& s) { return s.paused.load(); },
            [](TrainerSession& s, bool v) { s.paused = v; })
        .def_property("stop_requested",
            [](const TrainerSession& s) { return s.stop_requested.load(); },
            [](TrainerSession& s, bool v) { s.stop_requested = v; })
        .def_property_readonly("render_pending",
            [](const TrainerSession& s) { return s.render_pending.load(); })
        .def_property("cur_step",
            [](const TrainerSession& s) { return s.cur_step.load(); },
            [](TrainerSession& s, int v) { s.cur_step = v; })

        // Progress / warning sink. Called from phases that have released the
        // GIL, so the wrapper takes it back.
        .def("set_log_fn",
             [](TrainerSession& s, std::function<void(const std::string&)> fn) {
                 if (!fn) { s.log_fn = nullptr; return; }
                 s.log_fn = [fn](const std::string& msg) {
                     py::gil_scoped_acquire acquire;
                     fn(msg);
                 };
             },
             py::arg("fn"), py::keep_alive<1, 2>(),
             "Route log output to a Python callable (default: stdout).")

        .def("check_config", [](TrainerSession& s) { s.check_config(); },
             "Throw for unported features; warn for approximated ones.")
        .def("load_dataset", [](TrainerSession& s) {
                 py::gil_scoped_release release;
                 s.load_dataset();
             }, "Parse the dataset and bake the POST-split cameras. No GPU work.")
        .def("setup_engine", [](TrainerSession& s) {
                 py::gil_scoped_release release;
                 s.setup_engine();
             }, "Create the output dir, dump config.json, reset + seed the "
                "engine, install the DataManager and bilagrid/PPISP.")

        .def("train",
             [](TrainerSession& s,
                std::function<void(const TrainerProgress&)> on_step) {
                 TrainerCallbacks cb;
                 if (on_step) {
                     // Called from the (GIL-free) training thread, so take
                     // the GIL back for the duration of the Python call.
                     cb.on_step = [&on_step](const TrainerProgress& p) {
                         py::gil_scoped_acquire acquire;
                         on_step(p);
                     };
                 }
                 py::gil_scoped_release release;
                 s.train(cb);
             },
             py::arg("on_step") = py::none(),
             "Run the whole loop, honouring paused / stop_requested and "
             "saving checkpoints on the configured cadence.")

        .def("train_step",
             [](TrainerSession& s, int step) {
                 py::gil_scoped_release release;
                 return s.train_step(step);
             },
             py::arg("step"),
             "One step: build this step's EngineStepConfig and run it. "
             "Returns the loss dict. Caller must hold engine_lock().")

        .def("save_checkpoint",
             [](TrainerSession& s, int step) {
                 py::gil_scoped_release release;
                 s.save_checkpoint(step);
             }, py::arg("step"))

        .def("progress_json", &TrainerSession::progress_json,
             "The /progress body: step, total_steps, elapsed_time, eta, "
             "latency_ms, paused.")
        .def("make_viewer_config", &TrainerSession::make_viewer_config,
             "ViewerRenderConfig for this run (also used by the GUI viewport).")

        .def("engine_lock",
             [](TrainerSession& s) {
                 return ssplat_bindings::EngineLock(&s.engine_mutex);
             },
             py::keep_alive<0, 1>(),
             "Context manager over the engine mutex shared with the viewer "
             "render workers. Releases the GIL while blocking.");
}
