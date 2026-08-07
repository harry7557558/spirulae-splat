# Testing

Three suites, in descending order of how much you should trust them.

## 1. Native cross-backend parity tests (the important ones)

`src/backend/tests/*.cpp` — currently 17 tools covering projection (fwd, bwd,
quant-grad), rasterization bwd, tile intersect, warp, FPBO, optimizer (general
+ geometry), densify, per-pixel train, PPISP, bilagrid, multi-scale loss, plus
`backend/tests/engine/` which drives the *real* engine end to end
(render parity, train parity). `backend/vulkan/tests/` adds 3 Vulkan-only
smoke tests (runtime, pipeline, sort/scan).

**The same source builds under both backends.** Each file only touches
`backend::`, the generated launch declarations, and `Tensor.h`. The workflow
is dump-then-compare:

```bash
# on the CUDA machine
./build_cuda/projection_parity dump ref.bin
# on the target machine / device
./build/projection_parity compare ref.bin
```

Inputs are deterministic, and comparison is tolerance-based — fast-math
exp/sqrt chains legitimately differ across compilers, and borderline-cull
flips change whole rows, so a small allowance for those is built in.

Several tools also take a `*_DUMP_GOT` environment variable
(`FPBO_DUMP_GOT`, `PPISP_DUMP_GOT`, `MSLOSS_DUMP_GOT`, `PWTRAIN_DUMP_GOT`,
`BILAGRID_DUMP_GOT`, `DENSIFY_DUMP_GOT`) to write the *actual* values
alongside the reference, which is how you diff a mismatch numerically instead
of guessing.

### Building them

```bash
# CUDA branch: opt-in
bash build_develop.bash -B build_cuda -DSS_BACKEND=cuda -DSS_BUILD_BACKEND_TESTS=ON
# Vulkan branch: built unconditionally
bash build_develop.bash -DSS_BACKEND=vulkan
```

Each `.cpp` becomes an executable of the same base name in the build dir.

### Cross-machine / cross-vendor runs

The comparison target is often a different machine (e.g. an AMD GPU box), and
often offline. The pattern that works:

1. Transfer a matching `slangc` to the target and point `-DSS_SLANGC=` at
   it — SPIR-V is compiled at build time and never committed, so the target
   needs a compiler, and the version is pinned.
2. Dump references on the CUDA host.
3. Copy the `.bin` files over and run `compare` on the target.

Keep reference dumps out of git (`parity_refs/` is gitignored).

## 2. GUI / viewer checks

The web viewer can be driven headlessly over the Chrome DevTools Protocol.
Headless defaults to SwiftShader; to exercise a real GPU, run against a real
display (`DISPLAY=:0`).

If you script a Python training run that serves the viewer, pass
**`--no-keep-viewer-alive`** or the process hangs at exit waiting on it.

## 3. Python tests (`tests/python/`)

```bash
pytest tests/python/
```

`test_per_pixel.py`, `test_ppisp.py`, `test_rasterization.py`,
`test_render_background_sh.py`, `test_ssim.py`, plus `utils.py`.

**These are unmaintained and may fail.** They date from the PyTorch-first era,
compare CUDA kernels against `_torch_impl.py` reference implementations, and
have not been run recently. Treat a failure here as "needs investigation",
not as a release blocker — but do investigate, since they are the only
remaining check on the `_torch_impl.py` reference math.

## 4. Dataparser gate

`tests/python/test_dataparser_parity.py` runs COLMAP (text + binary),
Nerfstudio and Metashape fixtures through the native parsers and checks the
frame set, poses, intrinsics, distortion, seed cloud, train-frame scalars and
both split sides against `dataparser_golden.json`.

Those golden values are the frozen result of the §4.1 parity comparison: they
were captured while `modules/dataparser.py` still had its implementation, and
certified field-by-field against it (4 formats x 4 config variants x 2 splits)
before that implementation was deleted. Regenerate only when a parser change
is intentional, with `make_dataparser_golden.py`, and say in the commit
message what moved.

Two checks do not depend on the golden and are worth keeping in mind: the
train and eval splits must partition the frames (a bug that dropped frames
from *both* sides would leave each side self-consistent), and all four
fixture formats must describe the same scene (otherwise the golden agrees on
garbage).

Fixtures come from a fixed seed (`tests/python/dataset_fixtures.py`), so the
test depends on no dataset present on the machine.

```bash
pytest tests/python/test_dataparser_parity.py -q
```

Unlike the rest of `tests/python/`, **this one is maintained**.

## 5. Web-viewer binding smoke test

`tests/python/test_webviewer.py` brings up the **native** viewer server from
Python against a real engine (fixture -> `parse_dataset` -> `bake_post_split`
-> `engine_setup_data_manager`), checks `/`, `/buffers`, `/progress` and
`/pause-toggle`, exercises `engine_lock()`, and asserts `stop()` returns and
releases the port. Needs a CUDA device.

```bash
pytest tests/python/test_webviewer.py -q
```

It also covers `WebViewer.start_for_session()`, the wiring that lets the
viewer read its step counter / pause flag / progress JSON straight off a
`TrainerSession` instead of having Python push them.

## 6. Trainer gate

`tests/python/test_trainer_parity.py` is the §4.3 gate. Three parts:

1. **Config conversion** — `to_native_config(PresetClass())` must equal
   `TrainConfig()` + `train_apply_preset(name)` for all seven presets. These
   are *two live representations*: `src/config/TrainConfig.h` is the source of
   truth and the Python dataclasses are the downstream copy, so this test is
   what catches the copy drifting until the dataclasses are deleted.
2. **Per-step `EngineStepConfig`** — 8 config variants × 4 run states × 20
   steps chosen to straddle every warmup/decay boundary, checked against
   `step_config_golden.json`. Those 640 configs are the frozen result of the
   parity comparison against `model.py::engine_train_step_managed` +
   `core.py::_build_{optim,densify}_config`: all 42,880 fields matched, and
   the certification was re-run against the Python implementation immediately
   before deleting it. Drift in a ported LR schedule now fails here instead of
   showing up as a quality regression 20k steps into a run. Regenerate only
   deliberately, with `make_step_config_golden.py`.
3. **The rewired `Trainer`** — that it really runs on a `TrainerSession`
   (attach mode pulled the session's world into the model's host params, the
   session's `RunState` is what the model reports) and that it serves the
   native viewer.

Plus an opt-in end-to-end run (`SS_TEST_DATASET`, with
`SS_TEST_IMAGE_DIR` / `SS_TEST_DOWNSCALE` for pre-downscaled academic
sets) that trains through `TrainerSession` with the refine window pulled
inside the run — `refine_stop_num_iter` counts back from the *end*, so a short
run with the defaults never densifies.

```bash
pytest tests/python/test_trainer_parity.py -q
SS_TEST_DATASET=/path/to/mipnerf360/garden \
  SS_TEST_IMAGE_DIR=images_4 SS_TEST_DOWNSCALE=4 \
  pytest tests/python/test_trainer_parity.py -q
```

Why the golden is the config and not a loss curve: comparing loss curves
between the two drivers could never be exact, because they seeded the splats
from different RNGs (`std::mt19937` vs torch). A 2000-step run on Mip-NeRF 360
`garden` put them within 0.4% on final `rgb_loss` with an identical
splat-count trajectory — but two *identical* Python runs differed by the same
~1 dB PSNR, so that only ever bounded the drift at "below seeding noise".
Config equality was the exact statement, and it is what got frozen.

### Regenerating a golden

Both goldens (`step_config_golden.json`, `dataparser_golden.json`) encode a
proof that no longer has a second implementation to re-derive it from. So:

- Regenerate only when the change is *intended*, and describe it in the commit
  message. A regenerated golden with no explanation is indistinguishable from
  a silently broken schedule.
- Read the diff. The step-config golden splits step-invariant fields into
  `constant` and the scheduled ones into `per_step`, so a diff shows exactly
  which quantities moved.

## What to run before calling a change done

| change | gate |
|---|---|
| any kernel | CUDA build + Vulkan build + the relevant parity test on both |
| engine logic | both builds + `engine_render_parity` + `engine_train_step`-level check |
| config field | add the row in `src/config/TrainConfig.h`; check `spirula train --help` and the GUI's All Options editor |
| training-loop logic | `TrainerCore.cpp` — `build_step_config()` is the only place it lives |
| build system | every mode in [build.md](build.md) |
| anything | one short training run per backend on a public scene |

## Profiling

`SS_PROFILE=1` enables the env-gated per-stage timing breakdown
(H2D / D2H / D2D / memset / device / host). Header-only, works on both
backends — the right first tool when a backend is unexpectedly slow rather
than wrong.
