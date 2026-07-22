# Architecture

## Layers

```
                 ┌──────────────────────────────────────────────────────────┐
  front ends     │ ssplat-train (CLI)    ssplat-gui (ImGui)   Python package │
                 │ src/app/cli/          src/app/gui/         spirulae_splat/│
                 └───────────┬────────────────┬──────────────────┬──────────┘
                             │                │                  │
                    ┌────────┴────────────────┴───┐              │ pybind11
                    │  TrainerCore (src/app/)     │              │ src/bindings/
                    │  config → dataset → seeding │              │   ext.cpp
                    │  → step loop                │              │
                    └──────────────┬──────────────┘              │
                                   │                             │
                 ┌─────────────────┴─────────────────────────────┴──────────┐
   engine        │  src/engine/*.cpp — torch-free, CUDA-free, process-global │
                 │  state in EngineState.h; src/data/ owns image I/O         │
                 └──────────────────────────┬───────────────────────────────┘
                                            │  launch declarations
                 ┌──────────────────────────┴───────────────────────────────┐
   backend seam  │  src/backend/api/*.h  +  BackendRuntime.h                 │
                 └───────────┬──────────────────────────────┬───────────────┘
                             │                              │
              ┌──────────────┴───────────┐     ┌────────────┴──────────────┐
   backends   │ CUDA: src/kernels/**/*.cu│     │ Vulkan: src/backend/vulkan│
              │  + src/instantiations/   │     │  launchers + shaders/     │
              └──────────────┬───────────┘     └────────────┬──────────────┘
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
| dataset parsing (Python) | `spirulae_splat/modules/{dataparser,colmap_utils,metashape_utils,camera_utils}.py` — *duplicate; being retired* |
| image cache / prefetch / warp | `src/data/DataManager.cpp` |
| training session orchestration | `src/app/TrainerCore.cpp` and `spirulae_splat/modules/trainer.py` — *duplicate; being retired* |
| the actual training step | `Engine*.cpp`, entered via `engine_train_step_managed` |
| kernels | `src/kernels/**/*.cu` (CUDA) + `src/backend/vulkan/shaders/*.slang` + `src/backend/vulkan/kernels/*.cpp` |
| web viewer client | `spirulae_splat/viewer/viewer.html` — single source; the C++ viewer embeds it at build time |
| web viewer server | `src/app/webviewer/{Viewer,HttpServer,RenderWorker}.cpp` and `spirulae_splat/viewer/*.py` — *duplicate; being retired* |
| standalone WASM viewer | `viewer/` — independent, but compiles the C++ parsers in place |

Three subsystems are currently implemented twice, once in Python and once in
C++. That is being collapsed onto the C++ side; see
[restructure-proposal.md](restructure-proposal.md) §4. **Do not add a fourth.**

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
