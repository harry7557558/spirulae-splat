# Architecture

## Layers

```
                ┌────────────────────────────────────────────────────────────┐
  front ends    │ ssplat train (CLI)    ssplat (ImGui GUI)   Python package  │
                │ src/app/cli/          src/app/gui/         spirulae_splat/ │
                └────────────┬────────────────┬──────────────────┬───────────┘
                             │                │                  │
                    ┌────────┴────────────────┴──────────────────┴┐ pybind11
                    │  TrainerCore (src/app/)                     │ src/bindings/
                    │  config → dataset → seeding → step loop     │   ext.cpp
                    │  Python enters here too: _C.TrainerSession  │   bind_data
                    └───────────────────────┬─────────────────────┘   bind_viewer
                                            │                         bind_trainer
                 ┌──────────────────────────┴────────────────────────────────┐
   engine        │  src/engine/*.cpp — torch-free, CUDA-free, process-global │
                 │  state in EngineState.h; src/data/ owns image I/O         │
                 └──────────────────────────┬────────────────────────────────┘
                                            │  launch declarations
                 ┌──────────────────────────┴────────────────────────────────┐
   backend seam  │  src/backend/api/*.h  +  BackendRuntime.h                 │
                 └───────────┬──────────────────────────────┬────────────────┘
                             │                              │
              ┌──────────────┴────────────┐    ┌────────────┴───────────────┐
   backends   │ CUDA: src/kernels/**/*.cu │    │ Vulkan: src/backend/vulkan │
              │  + src/instantiations/    │    │  launchers + shaders/      │
              └──────────────┬────────────┘    └────────────┬───────────────┘
                             └──────────┬───────────────────┘
                                        │
                    src/shaders/*.slang — shared device math
                    src/backend/common/SortScan.h — CUB / onesweep radix
```

The seam is the important part, and it is the best-documented part of the
repo: read `src/backend/README.md`, then `src/backend/vulkan/README.md`.

## Where responsibilities live

| responsibility | code |
|---|---|
| training config (source of truth) | Python dataclasses in `spirulae_splat/modules/` → codegen → `src/app/generated/cli_config.h` |
| dataset parsing (native) | `src/data/parsers/{Colmap,Nerfstudio,Metashape}Parser.cpp`, `DatasetCommon.cpp`, `DatasetParser.h` |
| dataset parsing (Python client) | `spirulae_splat/modules/native_dataparser.py` — an adapter, not a parser. The Python implementation is gone; `dataparser.py` is now just the config dataclass, `scripts/{colmap,metashape}_utils.py` keep a Python reader for preprocessing, and `camera_utils.py` is retained on no code path as the reference for the unported `orientation_method` / `center_method` (docs/notes/pose-normalization.md) |
| image cache / prefetch / warp | `src/data/DataManager.cpp` |
| training session orchestration | `src/app/TrainerCore.cpp` — single implementation, bound as `_C.TrainerSession`; `spirulae_splat/modules/trainer.py` keeps config construction, resume, eval and profiling |
| the actual training step | `Engine*.cpp`, entered via `engine_train_step_managed` |
| kernels | `src/kernels/**/*.cu` (CUDA) + `src/backend/vulkan/shaders/*.slang` + `src/backend/vulkan/kernels/*.cpp` |
| web viewer client | `src/app/webviewer/viewer.html` — single source, embedded into the engine library at build time; every front end serves the same bytes |
| web viewer server | `src/app/webviewer/{Viewer,HttpServer,RenderWorker}.cpp` — single implementation, bound as `_C.WebViewer`; `ss_trainer.py` and `ss_viewer.py` drive it |
| standalone WASM viewer | `viewer/` — independent, but compiles the C++ parsers in place |

Three subsystems used to be implemented twice, once in Python and once in
C++. All three are now single implementations in C++, reached from Python
through a binding: dataset parsing (§4.1) via
`spirulae_splat.modules.native_dataparser`, the viewer server (§4.2) via
`_C.WebViewer`, and the training driver (§4.3) via `_C.TrainerSession` /
`spirulae_splat.modules.native_trainer`. The Python implementations were
deleted after a parity gate proved each pair agreed; the gates survive as
golden-value regression tests on the C++ side ([testing.md](testing.md)
§4-6). **Do not add a fourth** — new functionality goes in C++ with a
binding.

What is left in Python is what has no C++ counterpart: the config
dataclasses (the source of truth codegen reads), checkpoint resume and its
layout adaptation, the eval metrics (LPIPS and the SSIM variants are torch
models), and the image loading that eval needs.

## The engine

`Engine.h` is hand-maintained (not generated) and groups its API as:
lifecycle (`engine_reset`), setup (`set_data_3dgs`, `set_camera_params`,
`set_training_data[_warped]`), forward, loss+backward, optimizer step,
bilagrid (RGB/depth/normal), PPISP, background blending, color space,
densification, fused train steps, DataManager-driven train step, debug
render, state queries, and viewer blit.

Two properties to keep in mind:

- **It is a process-global singleton.** `EngineState` holds the world splats,
  camera table, forward cache, GT data, gradients, optimizer moments, bilagrid
  / PPISP / background / color-space state and the viewer state. Call
  `engine_reset()` between runs that swap datasets — `ss_benchmark`'s
  per-scene subprocess isolation exists for exactly this reason.
- **It is torch-free and CUDA-free.** `Engine*.cpp` are `.cpp`, not `.cu`;
  they reach the device only through the generated launch declarations and
  `backend/api/BackendRuntime.h`. Keep it that way — it is what lets the same
  engine drive the Vulkan backend.

## Ordering rules worth knowing

Image-space corrections compose in a fixed order:

```
render → background blending → bilagrid → PPISP → loss
```

`Engine.h` states this at the declaration sites (`applied BEFORE
bilagrid/PPISP`, `applied AFTER bilagrid`). Changing the order changes results
silently, so it is asserted only by parity tests.

## Primitives

3DGS, Mip (anti-aliased 3DGS), and 3DGUT are compile-time *types*, not runtime
branches: `Primitive.cuh`, `Primitive3DGS.cuh`, `Primitive3DGUT.cuh`,
`PrimitiveBase3DGS.cuh`. Camera model, SH degree and distortion are value
parameters. This is what makes `src/instantiations/` (111 generated instantiation TUs)
necessary on the CUDA side, and what maps to Slang specialization constants on
the Vulkan side.
