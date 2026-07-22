# Repository restructure proposal

Status: phases 0-4 applied (see §8); §4 and phases 5-7 still proposal.

Goal: make the tree navigable for humans and agents, remove duplicated
Python/C++ subsystems, split oversized files, and land durable documentation
that survives moving between machines. Hard constraints: **CUDA and Vulkan
backends both keep working**, and **Python/PyTorch keeps working**.

---

## 1. What's wrong today (measured)

### 1.1 Depth/breadth imbalance

Every piece of native code lives at `spirulae_splat/splat/cuda/csrc/` — five
path segments before the first real file, and each segment is misleading:

| segment | why it's wrong |
|---|---|
| `splat/` | vestige of the vendored gsplat drop. `spirulae_splat/splat/README.md` is literally gsplat's README, which reads as if it documents this directory. |
| `cuda/` | contains the Slang shaders and the **Vulkan** backend. |
| `csrc/` | contains the engine, the dataset parsers, the CLI, the native GUI, the web-viewer HTTP server, and the backend abstraction — not just "C sources". |

Then that one directory holds **144 tracked files flat** (excluding
`generated/`, `external/`, `ins/`). By prefix:

```
Bilagrid*  43    Engine*    20    Rasterization* 14    Projection* 13
Mesh*       9    Primitive*  4    + ~30 ungrouped singles
```

Meanwhile `backend/`, `app/`, `app/gui/` are properly nested. So the tree is
simultaneously too deep (the prefix) and too flat (the leaf).

### 1.2 Oversized files

Non-generated, non-vendored files over ~700 lines:

```
3668  csrc/PixelWise.cu                    1274  backend/vulkan/kernels/Bilagrid.cpp
2655  csrc/Optimizer.cu                    1182  slang/vulkan/visualizer.slang
2471  csrc/Densify.cu                      1172  csrc/Meshing.cu
2412  slang/harmonics.slang                1144  csrc/FusedProjectionBwdOptim_kernel.cuh
1980  csrc/Delaunay3D.cpp                  1142  slang/ppisp.slang
1957  csrc/DataManager.cpp                 1131  csrc/PerPixelLoss.cu
1806  csrc/Visualizer.cu                   1118  csrc/FusedSSIM.cu
1778  csrc/MeshingHost.cpp                 1097  backend/vulkan/VulkanRuntime.cpp
1478  app/gui/GuiApp.cpp                   1086  slang/vulkan/densify.slang
1377  csrc/Tensor.h                        1055  slang/projection_utils.slang
1329  csrc/MeshUV.cpp                       996  csrc/SplatTileIntersector.cu
1318  viewer/src/viewer.cpp                 995  csrc/Common.cuh
```

Plus `CMakeLists.txt` at 723 lines covering four build modes,
`spirulae_splat/modules/model.py` at 2141, `modules/trainer.py` at 1189.

### 1.3 Duplicated Python / C++ subsystems

Three overlapping stacks, each maintained twice:

| subsystem | Python | C++ | shared today |
|---|---|---|---|
| dataset parsing | `modules/dataparser.py` (903), `colmap_utils.py` (646), `metashape_utils.py` (397), `camera_utils.py` (648), `dataset.py` (313), `datamanager.py` (128) ≈ **3.0k lines** | `app/ColmapParser.cpp` (605), `NerfstudioParser.cpp` (571), `MetashapeParser.cpp` (602), `DatasetCommon.cpp` (429), `DatasetParser.h` (293) ≈ **2.5k lines** | nothing |
| web viewer server | `viewer/server.py` (100), `http_server.py` (262), `render_worker.py` (176), `annotation.py` (199) ≈ **0.7k** | `app/Viewer.cpp` (224), `HttpServer.cpp` (250), `RenderWorker.cpp` (583) ≈ **1.1k** | `viewer.html` only (CMake embeds it — good precedent) |
| training driver | `modules/trainer.py` + `model.py::engine_train_step_managed` | `app/TrainerCore.{h,cpp}` (914) | the engine underneath |

The third one is the most fragile: `app/TrainerCore.h`'s own header comment
documents itself as a *line-by-line port* of the Python managed path, with a
mapping table in `app/README.md`. That is a drift bomb — the comment is
already doing the job a shared implementation should do. The C++ WASM viewer
already proves the right pattern: `viewer/CMakeLists.txt` compiles
`csrc/app/*Parser.cpp` **in place** rather than reimplementing.

Encouragingly, the hard part is already done: `modules/trainer.py` and
`model.py` drive the C++ engine (`engine_train_step_managed`,
`engine_setup_data_manager`, …) — PyTorch is no longer the compute path. What
remains duplicated is dataset parsing, viewer serving, and the step-loop
orchestration.

### 1.4 Four different things named "viewer"

```
viewer/                        standalone WebGL2 + WASM viewer (independent)
spirulae_splat/viewer/         Python HTTP server + viewer.html
csrc/app/Viewer.cpp            C++ HTTP server serving the same viewer.html
csrc/app/gui/                  native ImGui desktop GUI
spirulae_splat/viewer_legacy/  stale, untracked, only __pycache__ survives
```

### 1.5 Two independent build systems

`setup.py` (torch `CUDAExtension`, globs `csrc/*.cu`, `csrc/*.cpp`, `ins/*.cu`)
and `CMakeLists.txt` (723 lines, backends × torch × CLI × GUI × tests) each
carry their own source lists, flags, and codegen invocation. They must be kept
in sync by hand.

### 1.6 Smaller papercuts

- `spirulae_splat/generate_{headers,cli_config,kernel_instantiation,backend_api,vulkan_stubs}.py` — five codegen tools sitting in the Python package root next to runtime modules.
- `csrc/tests/test_delaunay3d.py` — a Python test inside the C++ tree; separate from root `tests/` (5 Python tests) and `csrc/backend/tests/` (17 C++ parity tests) and `csrc/backend/vulkan/tests/` (3).
- `csrc/stb_image_impl.cpp`, `stb_image_write_impl.cpp` at csrc root instead of `external/`.
- `spirulae_splat/scripts/` contains only `.gitignore` + `models/.gitkeep` — a runtime cache directory shaped like a source directory.
- `spirulae_splat/{perf_timer,cuda_resource_viewer}.py` loose at package root.
- `.gitmodules` is a fully commented-out stale file (already deleted in the working tree — keep it deleted).
- Untracked `csrc/BilagridBackwardSelection.md` — a design note that belongs in `docs/`.

### 1.7 Private information currently in the public tree

```
scripts/batch_process_data.bash:10        /mnt/d/gs/data/vocab_tree_...
spirulae_splat/modules/edge_detector.py:40-47   /mnt/d/gs/data/{360_v2/bicycle_4,
                                          Atrium_Godiva, adr/20260221-queens_park,
                                          adr/20251109-elevator-insta360-...,
                                          zipnerf/nyc}
tests/test_ppisp.py:72                    /mnt/d/plot.png   (commented)
tests/test_rasterization.py:367           /mnt/d/plot.png   (commented)
```

`360_v2` and `zipnerf` are standard academic sets and are fine; the rest name
private captures and a local drive layout. `ss_benchmark.py` is clean (paths
are CLI arguments; only academic scene names appear).

---

## 2. Proposed target layout

```
spirulae-splat/
├── AGENTS.md                  # canonical agent/human orientation doc
├── CLAUDE.md                  # -> AGENTS.md (symlink, or 1-line include)
├── README.md                  # user-facing: install + quick start only
├── CMakeLists.txt             # thin: options + add_subdirectory + include(cmake/*)
├── cmake/
│   ├── Slang.cmake            # slangc discovery/fetch + SPIR-V embed
│   ├── BackendCuda.cmake
│   ├── BackendVulkan.cmake
│   ├── TorchExtension.cmake
│   ├── Apps.cmake             # ssplat-train / ssplat-gui / ssplat-mesh
│   └── Tests.cmake
├── docs/
│   ├── architecture.md        # layer diagram + data flow + engine state model
│   ├── build.md               # the build matrix, one table
│   ├── backends.md            # api seam, CUDA vs Vulkan, coverage table
│   ├── adding-a-kernel.md     # the both-backends checklist
│   ├── codegen.md             # what is generated from what; invariants
│   ├── datasets.md            # supported dataset layouts
│   ├── testing.md             # parity tests, python tests, how to run
│   └── notes/                 # design notes (bilagrid bwd selection, etc.)
├── src/                       # ALL native code (was spirulae_splat/splat/cuda/csrc)
│   ├── core/                  # Tensor.h, Camera.h, CameraModel.h, Common.cuh,
│   │                          #   Interpolation.cuh, PoolSlots.h, GradQuant.cuh,
│   │                          #   CheckpointIO.h, NonShQuantState.h
│   ├── primitives/            # Primitive.cuh, Primitive3DGS/3DGUT, PrimitiveBase3DGS
│   ├── kernels/
│   │   ├── projection/        # ProjectionFwd/Bwd/PackedFwd/QuantGrad (+_kernel.cuh)
│   │   ├── raster/            # Rasterization{Fwd,Bwd,Eval3D*,MomentsFwd}
│   │   ├── tile/              # IntersectTile, SplatTileIntersector
│   │   ├── pixelwise/         # PixelWise.* (split, see §3.1)
│   │   ├── bilagrid/          # the 43 Bilagrid* files
│   │   ├── optim/             # Optimizer.*, FusedProjectionBwdOptim.*
│   │   ├── densify/           # Densify.*, Quantile.cu
│   │   ├── loss/              # PerPixelLoss, PerSplatLoss, FusedSSIM, ColorShiftReg
│   │   ├── background/        # BackgroundSphericalHarmonics
│   │   ├── ppisp/             # PpispInit, PPISP parts of PixelWise
│   │   └── visualize/         # Visualizer.*
│   ├── mesh/                  # Meshing*, MeshUV, MeshExport, Delaunay3D
│   ├── engine/                # Engine*.cpp/.h
│   ├── data/                  # DataManager.*, parsers/{Colmap,Nerfstudio,Metashape,
│   │                          #   Common}, DatasetParser.h, FastFloat.h, Json.h, Xml.h
│   ├── backend/               # unchanged: api/ common/ cuda/ vulkan/
│   ├── shaders/               # was splat/cuda/slang/  (+ shaders/vulkan/)
│   ├── app/
│   │   ├── cli/               # main.cpp, mesh_main.cpp
│   │   ├── gui/               # unchanged
│   │   ├── webviewer/         # Viewer.cpp, HttpServer.cpp, RenderWorker.cpp,
│   │   │                      #   viewer.html  (single source, embedded)
│   │   └── TrainerCore.{h,cpp}
│   ├── bindings/              # ext.cpp split per subsystem (see §3.2)
│   ├── generated/             # unchanged, .gitattributes-marked
│   ├── instantiations/        # was splat/cuda/ins/
│   └── external/              # + stb_image_impl.cpp, stb_image_write_impl.cpp
├── spirulae_splat/            # thin Python package (see §4)
│   ├── __init__.py
│   ├── _native.py             # the one place that imports the extension
│   ├── cli/                   # ss_trainer, ss_benchmark, ss_viewer, ss_meshing
│   ├── config/                # the dataclasses that codegen the C++ CLI config
│   ├── training/              # trainer.py, model.py (split), optimizer.py, resume*
│   ├── data/                  # thin adapters over the C++ parsers
│   ├── metrics/               # lpips.py, lpips_models/, debug_image.py
│   └── utils/                 # perf_timer, verbose, cuda_resource_viewer, _profile
├── tools/
│   └── codegen/               # generate_{headers,cli_config,kernel_instantiation,
│                              #   backend_api,vulkan_stubs}.py
├── tests/
│   ├── python/                # existing tests/*.py + csrc/tests/test_delaunay3d.py
│   └── native/                # was csrc/backend/tests + backend/vulkan/tests
│                              #   + csrc/tests/delaunay3d_bench.cpp
├── scripts/                   # unchanged (data preprocessing utilities)
└── viewer/                    # standalone WASM/WebGL2 viewer (name fixed by
                               #   the published GitHub Pages URL)
```

Two things this deliberately does **not** change: the `backend/{api,cuda,
vulkan,common}` seam (it's already the best-organized part of the repo, and
its README is accurate), and the codegen contracts.

### Cost of the move

Local `#include "..."` directives across the native tree: **589 total**, with
a long tail (top: `generated/set_namespace.cuh` ×26, `BilagridConfig.cuh` ×25,
`backend/api/BackendTypes.h` ×20, `Common.cuh` ×19). Path-string references to
`splat/cuda/csrc` in tracked files: **65**, in 12 files (`CMakeLists.txt`, the
five codegen scripts, `.gitattributes`, `viewer/CMakeLists.txt`,
`viewer/src/dataset_bridge.cpp`, two `_kernel.cuh` comments,
`modules/resume_codecs.py`, `harmonics.slang`).

So the mechanical cost is one scripted commit. Recommendation: **rewrite
includes to be path-qualified** (`#include "core/Common.cuh"`) rather than
adding every subdirectory to `target_include_directories` — the latter keeps
the diff at zero but leaves the tree just as unnavigable to a reader, which is
the whole point of the exercise.

If the full move is judged too disruptive, the **fallback** is to do §2's
`src/` subdivision *inside* `spirulae_splat/splat/cuda/csrc/` and skip the
rename. That captures ~80% of the navigability win for ~40% of the churn.

---

## 3. File splits

### 3.1 The split mechanism — DONE

Three facts made this cheap:

1. The oversized files are already divided by `// ====` banner sections.
2. `generate_headers.py` can feed one `.cuh` from many `.cu` files.
   (It originally inferred them from a `Name.` filename prefix; that was
   replaced with an explicit `HEADER_SOURCES` map so a source file can be
   named after **what it does** rather than after the header it declares into,
   and so a rename fails loudly instead of silently dropping declarations.)
3. Both `CMakeLists.txt` and `setup.py` glob `*.cu`, so new files need no
   build-file change — *once* the CMake glob has `CONFIGURE_DEPENDS`, which it
   did not and now does. Without it a new source is silently not compiled and
   surfaces only as an undefined reference at link time.

Naming rules adopted: parts are named for their function; `<Name>_kernel.cuh`
stays reserved for device bodies consumed by
`generate_kernel_instantiation.py`; a family's shared preamble goes in
`<Family>Common.cuh`.

**`PixelWise.cu` (3668) → 7 files + 2 headers**, all declaring into the
unchanged `PixelWise.cuh`:

| new file | contents | lines |
|---|---|---|
| `ImageConvert.cu` | type conversion; rendered → expected depth | 330 |
| `ImageColorOps.cu` | background blending (plain + noise), log map, overexposure reg | 411 |
| `DepthGeometry.cu` | depth → points / normal, depth-normal loss, ray ↔ linear depth | 706 |
| `ImageDistort.cu` | distort / undistort | 109 |
| `ImageWarp.cu` | wide ↔ pinhole warps, incl. byte-fused | 524 |
| `GtDepthNormalWarp.cu` | GT depth/normal wide → pinhole warps | 767 |
| `Ppisp.cu` | PPISP | 670 |
| `PixelWiseCommon.cuh` | slang namespaces, host view helpers | 86 |
| `BilinearSample.cuh` | device samplers shared by the distort/warp parts | 139 |

**`Optimizer.cu` (2655) → 8 files + 1 header**: `TensorSetZero.cu`,
`AdamOptim.cu`, `NewtonOptim.cu`, `ScaleAgnosticMeanOptim.cu`,
`FusedGeometryOptim.cu`, `FusedAppearanceOptim.cu`,
`TrustRegion3DGS2Optim.cu`, `ColorOptim.cu`, `OptimizerCommon.cuh`. The common
header now also `#include`s `Optimizer.cuh`, which the original `.cu` did not —
so launcher signatures are checked against their declarations.

**`Densify.cu` (2471) → 5 files + 3 headers**: `DensifySampling.cu`,
`DensifyScoring.cu`, `Relocation.cu`, `McmcRelocation.cu`,
`DensifySplitFilter.cu`, plus `DensifyCommon.cuh`, `DensifyQuantCopy.cuh` (the
quantized optim-state copy/zero helpers shared by both relocation paths) and
`DensifyInternal.cuh` (two sampling primitives that were `static` but used from
three TUs — now external-linkage with a declared interface).

Verified: `PixelWise.cuh`, `Optimizer.cuh` and `Densify.cuh` keep byte-identical
declaration sets (43 / 26 / 20); all five generators produce no drift; both
backends build; 16/16 parity tests pass.

### 3.2 Remaining splits (not yet done)

| file | split into |
|---|---|
| `Visualizer.cu` 1806 | overlay / fields / blit |
| `DataManager.cpp` 1957 | cache / prefetch / warp / lifecycle |
| `MeshingHost.cpp` 1778 | extract / simplify / export |
| `Tensor.h` 1377 | `core/Tensor.h` (view types) + `core/DeviceTensor.h` + `core/Pool.h` |
| `Common.cuh` 995 | `core/Math.cuh` + `core/Launch.cuh` + `core/Common.cuh` |
| `ext.cpp` 559 (144 `m.def`s) | `bindings/bind_{kernels,engine,optim,mesh,data,viewer}.cpp` + `ext.cpp` calling `register_*(m)` |
| `app/gui/GuiApp.cpp` 1478 | panel-per-file, following the existing `*Panel.cpp` convention |
| `backend/vulkan/kernels/Bilagrid.cpp` 1274 | `Bilagrid{Sample,Optim,Tv}.cpp` |
| `CMakeLists.txt` 723 | `cmake/*.cmake` per §2 |
| `modules/model.py` 2141 | `training/{config,splats,forward,losses,engine_step}.py` |

`Delaunay3D.cpp` (1980) is a self-contained algorithm with a narrow header —
leave it, or mark `linguist-vendored`. Generated files (`generated/*.cuh`,
`ins/*.cu`) stay untouched per your note.

### 3.3 Bilagrid: 43 flat files

Beyond moving them into `kernels/bilagrid/`, note the shape:
six samplers (`Sample`, `UniformSample`, `DepthSample`, `NormalSample`,
`LoglinearSample`, `PpispUniformSample`) × `{Fwd, BwdV1, BwdV2}` `_kernel.cuh`
= 18 near-parallel files, plus the V1/V2 selector machinery
(`BilagridBwdSelector.h`, `BilagridBwdTiming.h`). Suggest:

```
kernels/bilagrid/
├── README.md              # the V1/V2 selector story (absorb BilagridBackwardSelection.md)
├── Config.cuh  Reader.cuh  AxisAngle.cuh  Utils.{cu,cuh}
├── selector/              # BwdSelector.h, BwdTiming.h
├── samplers/{affine,uniform,depth,normal,loglinear,ppisp_uniform}/{Fwd,BwdV1,BwdV2}.cuh
├── ppisp/                 # PpispMath{,Bwd}.cuh, PpispSample{Fwd,Bwd}.*
├── TvLoss{Fwd,Bwd}.cu  FusedAdam.cu  Init.cu
└── Bindings.{h,cpp}
```

---

## 4. Collapsing the Python/C++ duplication

Principle: **C++ is the single implementation; Python is a client.** This is
already the direction of travel (README calls the Python trainer "legacy…may
be deprecated"; the engine is torch-free; the WASM viewer compiles the C++
parsers in place). Three steps, each independently shippable and testable.

### 4.1 Dataset parsing → C++ (removes ~3.0k lines of Python)

Expose the existing parsers through pybind:

```cpp
// src/bindings/bind_data.cpp
py::class_<ssplat::DatasetSpec>(m, "DatasetSpec")   // cameras, image paths,
    .def_readonly("cameras", ...)                    // intrinsics, distortion,
    .def_readonly("points", ...)                     // c2w, train/eval split,
    ...;                                             // seed point cloud
m.def("parse_dataset", &ssplat::parse_dataset,
      py::arg("path"), py::arg("config"));           // auto-detects layout
```

Python side becomes a ~150-line adapter in `spirulae_splat/data/` that maps
`DatasetSpec` onto whatever structures `trainer.py` still wants. Delete
`dataparser.py`, `colmap_utils.py`, `metashape_utils.py`, and the parsing half
of `camera_utils.py` / `dataset.py`.

Watch-outs:
- `dataparser.py`, `metashape_utils.py`, `ss_meshing.py` import `nerfstudio`.
  Check whether the nerfstudio dependency survives this — if only the parsers
  needed it, dropping it removes a heavyweight dep from the Python install.
- Mask handling, depth/normal prior loading, and downscale-factor selection
  live on the Python side today; confirm parity in `DatasetCommon.cpp` before
  deleting, or port the gaps first.
- Keep `scripts/process_data_*.py` on the Python side — those are
  preprocessing tools, not the training data path.

**Verification gate:** a golden test that parses each of COLMAP / Nerfstudio /
Metashape fixtures through both paths and asserts identical camera matrices,
intrinsics, split assignment, and point counts. Run it *before* deleting the
Python implementation, keep it as `tests/python/test_dataparser_parity.py`
until the deletion commit.

### 4.2 Viewer server → C++ (removes ~0.7k lines of Python, kills the 4-way name clash)

`viewer.html` is already single-source (CMake embeds it into
`app_generated/viewer_html.h`). Do the same for the server:

```python
from spirulae_splat import _native
viewer = _native.WebViewer(port=7007, render_fn=my_render_callback)
viewer.start()
```

`app/Viewer.cpp` + `HttpServer.cpp` + `RenderWorker.cpp` move to
`src/app/webviewer/` and gain a pybind wrapper. Delete
`spirulae_splat/viewer/{server,http_server,render_worker}.py`. `annotation.py`
either ports to C++ or stays as a Python-only add-on layered on the C++
server's route table.

`ss_viewer.py` (the standalone `.ply` viewer) then becomes a thin CLI: load
PLY → `set_data_3dgs` → `WebViewer`. Longer term it can simply exec the native
binary.

Renames to end the ambiguity (`viewer/` stays — see phase 6);
`spirulae_splat/viewer/` deleted; `csrc/app/Viewer.cpp` → `app/webviewer/`;
`app/gui/` stays as the desktop GUI. Delete the stale
`spirulae_splat/viewer_legacy/`.

### 4.3 Training loop → shared `TrainerCore` (removes the documented port)

Highest drift risk, so do it last. `TrainerCore.h` already declares the phase
split (`check_config` → `load_dataset` → `setup_engine` → `train`) with
pause/stop atomics and a per-step callback — that is exactly a Python-friendly
API. Expose:

```python
sess = _native.TrainerSession(config)     # config from the same dataclasses
sess.check_config(); sess.load_dataset(); sess.setup_engine()
sess.train(on_step=lambda st: ...)        # GIL released inside
```

`modules/trainer.py` keeps: config construction, output-dir/JSON conventions,
eval metrics (LPIPS is a torch model), checkpoint resume/adapt/codecs, and
profiling. `model.py::engine_train_step_managed` and `_build_*_config` are
deleted in favour of the C++ ones — note those are the functions
`TrainerCore.h`'s comment says it mirrors, so this removes the mirror rather
than the original.

The config dataclasses stay the single source of truth
(`generate_cli_config.py` already generates the C++ side from them) — extend
that generator to also emit the pybind `TrainerConfig` struct so the three
representations can't diverge.

**Verification gate:** train a short run on a public scene through both paths
before/after and compare per-step loss to within float noise; the existing
`backend/tests/engine/engine_train_parity.cpp` is the model for this.

---

## 5. Documentation plan

### 5.1 `AGENTS.md` (root, canonical) + `CLAUDE.md` symlink

One screen of orientation plus links. Target contents:

- **What this is** — 3DGS trainer; native C++/CUDA+Vulkan engine; Python is a client.
- **Repo map** — the §2 table, one line per directory.
- **Build** — the four modes and the exact commands. Critically: *use
  `build_develop.bash` (Linux) / `build_develop.bat` (Windows), not `pip
  install`*, for development.
- **Codegen invariants** — what `generate_headers.py` (the
  `/*[AutoHeaderGeneratorExport]*/` marker and the `Name.*.cu` collection
  rule), `generate_kernel_instantiation.py`, `generate_cli_config.py`,
  `generate_backend_api.py`, `generate_vulkan_stubs.py` each own; never
  hand-edit below the `AUTO HEADER GENERATOR` splitter line; `.gitattributes`
  marks generated/vendored trees.
- **The two-backend rule** — every kernel change needs a CUDA impl, a Slang
  impl, and a parity test; link `docs/adding-a-kernel.md`.
- **Testing** — how to run native parity tests and Python tests.
- **Conventions** — file naming, `*_kernel.cuh` split, banner sections,
  namespace layout.
- **Gotchas worth writing down once** (currently only in per-machine agent
  memory):
  - never call `fs::remove_all` in a torch-linked executable;
  - quantized gradient codecs must decode code 0 to exactly `0.0`;
  - `refine_stop_num_iter` counts back from the end, so short runs never densify;
  - pass `--no-keep-viewer-alive` in scripted Python runs or training stalls at exit;
  - `SSPLAT_PROFILE=1` enables the backend timing breakdown.

Keep it under ~200 lines; push detail into `docs/`.

### 5.2 `docs/` set

Per §2. Notable moves: absorb `csrc/backend/README.md` +
`csrc/backend/vulkan/README.md` into `docs/backends.md` (leave short
pointer-READMEs in place), move the untracked
`csrc/BilagridBackwardSelection.md` into `docs/notes/`, and move the
gsplat-derivation content out of `spirulae_splat/splat/README.md` into a
root `THIRD_PARTY.md` / `NOTICE` so it stops looking like this repo's docs.

Also trim `README.md` (currently 25.8 KB, install + full CLI reference +
features): keep features/install/quick-start, move the exhaustive flag
reference to `docs/cli.md`.

### 5.3 Keep a per-directory `README.md`

One short paragraph in each `src/*` subdir saying what lives there and what
the invariant is. `backend/README.md` is a good model — accurate, has a
diagram, explains *why* the seam exists.

### 5.4 Privacy scrub (do this regardless of the restructure)

- `spirulae_splat/modules/edge_detector.py:40-47` — replace the hardcoded
  `/mnt/d/gs/data/...` demo paths (incl. private capture names) with an
  `argv[1]` argument; drop the commented alternatives.
- `spirulae_splat/modules/edge_detector.py:47` — `/mnt/d/temp.png` output →
  argument or `./`.
- `scripts/batch_process_data.bash:10` — `vocab_tree_path` → env var with a
  documented default, e.g. `${SSPLAT_VOCAB_TREE:?set to your vocab tree .bin}`.
- `tests/test_ppisp.py:72`, `tests/test_rasterization.py:367` — drop the
  commented `/mnt/d/plot.png` lines.
- Add a CI/pre-commit grep for `/mnt/`, `/home/<user>`, `C:\Users` so it
  doesn't regress.
- Note: these strings are in git *history* too. Rewriting history is probably
  not worth it for a handful of directory names — but worth a deliberate
  decision rather than an oversight.
- Standard academic dataset names (Mip-NeRF 360 / `360_v2`, ZipNeRF) stay.

---

## 6. Suggested sequencing

Each phase ends at a green gate: **CUDA build + Vulkan build + native parity
tests + `pytest tests/` all pass**, and one short training run on a public
scene per backend.

| # | phase | risk | reversible |
|---|---|---|---|
| 0 | Privacy scrub; delete `viewer_legacy/`, stale `.gitmodules`; move `BilagridBackwardSelection.md` → `docs/notes/` | none | yes |
| 1 | Write `AGENTS.md` + `docs/` skeleton **against the current layout** | none | yes |
| 2 | File splits inside existing dirs (§3) — pure `git mv`/text moves, no include churn thanks to the dotted-name codegen rule | low | yes |
| 3 | Split `CMakeLists.txt` into `cmake/*.cmake`; make `setup.py` consume a shared source-list (or delegate to CMake via scikit-build-core) | medium — touches all four build modes | yes |
| 4 | The directory move (§2) + path-qualified includes, as one scripted commit; update `AGENTS.md`/`docs/` paths in the same commit | medium (large diff, low semantic risk) | yes |
| 5 | §4.1 dataset parsing unification (parity test first, deletion second) | medium | yes |
| 6 | §4.2 viewer unification + the four-way rename | low | yes |
| 7 | §4.3 TrainerCore unification | high — do alone, with the loss-curve gate | yes |

Phases 0–2 are worth doing even if you stop there. Phase 4 is the one that
actually fixes navigability. Phase 7 is the one that stops the two trainers
from drifting, and is the one I'd want a real regression run behind.

---

## 7. Questions — answered

1. **Is the Python/PyTorch path a first-class target or a compatibility
   shim?** *Answered: a compatibility shim.* The project started as a
   Nerfstudio/gsplat fork and is being steered toward a C++ codebase, with the
   Python/PyTorch dependency on the way out. There are current users on it, so
   **breakage must be paced, not rushed** — deprecate with notice, keep the
   path working until users have moved. This makes §4.3 the eventual
   destination rather than an urgent fix, and means each of §4.1/4.2/4.3
   should land the C++ binding *first* and delete the Python implementation in
   a later, separately announced commit.
2. **Does the `nerfstudio` dependency survive §4.1?** *Answered: the engine
   and native apps already build and run without nerfstudio/gsplat.* The
   remaining imports are in `modules/{dataparser,metashape_utils}.py`,
   `ss_meshing.py` and three `scripts/`. So §4.1 removes the dependency from
   the *training* path; `scripts/` may keep it as a preprocessing-only dep.
3. **`src/` at repo root vs keeping native code inside `spirulae_splat/`.**
   Open. Given (1), root `src/` is the better end state — the Python package
   should look like a client, not the host.
4. **History rewrite for the leaked local paths** — open. The working tree is
   clean as of phase 0; only history remains.
5. **`Delaunay3D.cpp` / `MeshUV.cpp` provenance** — open. Affects whether they
   get `linguist-vendored` and are exempt from the split rules.

## 8. Progress

- **Phase 0 — done.** Private paths scrubbed (`edge_detector.py` now takes
  argv, `batch_process_data.bash` reads `SSPLAT_VOCAB_TREE`, commented
  `/mnt/d/plot.png` lines dropped); `check_private_paths.sh` added as
  the regression guard (now `tools/check_private_paths.sh`); stale
  `.gitmodules` and `spirulae_splat/viewer_legacy/` removed; `BilagridBackwardSelection.md` → `docs/notes/`.
- **Phase 1 — done.** `AGENTS.md` + `CLAUDE.md` pointer; `docs/` with
  `architecture`, `build`, `backends`, `codegen`, `datasets`, `testing`,
  `README` index, `notes/`. Written against the *current* layout, so they are
  correct today and get path updates in phase 4.
- **Phase 2 — done.** File splits inside existing directories: `PixelWise.cu`
  (3668) → 7 function-named `.cu` + `PixelWiseCommon.cuh` + `BilinearSample.cuh`;
  `Optimizer.cu` (2655) → 8 `.cu` + `OptimizerCommon.cuh`; `Densify.cu` (2471)
  → 5 `.cu` + 3 `.cuh`. `generate_headers.py` grew an explicit `HEADER_SOURCES`
  map so a source file's name is free of the header it feeds. Declaration sets
  of `PixelWise.cuh` / `Optimizer.cuh` / `Densify.cuh` verified unchanged
  (43 / 26 / 20).
- **Phase 3 — done.** `CMakeLists.txt` (723 lines) → a 30-line running order
  plus `cmake/`: `SsplatOptions`, `SsplatSources`, `SsplatBackendCuda`,
  `SsplatBackendVulkan`, `SsplatSlang`, `SsplatEmbed`, `SsplatApps`. The
  backend `if/else` is now `include(SsplatBackendCuda|Vulkan)`, each leaving
  `SSPLAT_WITH_TORCH` + `SSPLAT_APP_LIBS` for the shared app targets; the
  duplicated per-app boilerplate (libpython, `-static-libstdc++`, `BUILD_RPATH`,
  include dirs) collapsed into `ssplat_configure_app()`, and the two
  copy-pasted hex-embed blocks into `ssplat_embed_file()`.

  Source lists are no longer duplicated: `cmake/sources.txt` holds the glob
  patterns and is read by both `cmake/SsplatSources.cmake` and `setup.py`'s new
  `get_sources()`. Plain text rather than CMake so the pip build needs no CMake
  — full scikit-build-core delegation was rejected as too disruptive for users
  still on the pip path (§7.1).

  Verified by diffing the generated `build.ninja` old-vs-new for three configs
  (torch CUDA / no-torch CUDA+GUI+tests / Vulkan): identical target sets,
  **zero** flag, define, or link-library differences; the only deltas are
  source ordering and one extra harmless `-I<builddir>` on `ssplat-mesh`.
  Both backends rebuild clean and the parity suite is 17/17.

- **Phase 4 — done (native half).** The directory move. All native code is now
  under a root `src/`, subdivided per §2: `core/ primitives/ kernels/{projection,
  raster,tile,pixelwise,ppisp,bilagrid,optim,densify,loss,background,visualize}/
  engine/ data/{,parsers/} mesh/ backend/ shaders/ app/{cli,gui,webviewer}/
  bindings/ generated/ instantiations/ external/`. `spirulae_splat/splat/cuda/`
  no longer holds native code. Also moved: the five codegen tools →
  `tools/codegen/`, `tests/*.py` + the stray `csrc/tests/test_delaunay3d.py` →
  `tests/python/`, `csrc/tests/delaunay3d_bench.cpp` → `tests/native/`.

  **`src/` is now the include root** and every local include is path-qualified
  against it (`#include "core/Common.cuh"`) — 716 quoted plus 145 angle-bracket
  directives rewritten mechanically, none left relative. A file's include list
  now names the subsystems it depends on.

  Build/codegen plumbing updated in the same pass: `cmake/sources.txt` gained
  `[cuda]` / `[portable]` sections (the Vulkan build's old flat `csrc/*.cpp`
  glob no longer expresses "the torch-free subset"), `SpirvShaders.cmake` is
  parameterised on the shaders dir, all five generators emit and read
  src-relative paths, `.gitattributes`, `viewer/CMakeLists.txt` (which compiles
  the parsers in place) and `build_develop.{bash,bat}` follow.
  `generate_kernel_instantiation.py` now slugs its output filenames from
  include *basenames*, so path-qualifying did not rename all 111 generated TUs.

  Verified: every moved file's content diff is include lines only (checked
  blob-by-blob against HEAD, 52 residual lines, all include directives or the
  one comment naming a header); CUDA build, Vulkan build and Torch-extension
  build all green; parity 17/17 plus the three Vulkan smoke tests; all five
  generators re-run with zero drift.

  **Deliberately deferred:** §2's reorganisation *inside* the Python package
  (`spirulae_splat/{cli,config,training,data,metrics,utils}`). It changes
  public import paths (`spirulae_splat.splat.cuda`,
  `spirulae_splat.modules.*`) that `scripts/` and external users depend on,
  which per §7.1 should be paced. So `spirulae_splat/splat/` still exists,
  holding only Python.

- **Phase 4b — backend symmetry, partial (done).** Reviewed whether the CUDA
  and Vulkan halves should mirror each other. Conclusion: make each backend
  *self-contained*, not make the two *look* like peers.

  Done: `src/shaders/vulkan/` (39 files) → `src/backend/vulkan/shaders/`, and
  `SpirvShaders.cmake` + `spirv_tool.cpp` with it. The Vulkan backend was the
  only subsystem split across two top-level directories (host code under
  `backend/`, device code under `shaders/`); it is now one directory —
  runtime + kernels + shaders + its own SPIR-V build. `src/shaders/` now means
  exactly what it says: the 8 Slang files compiled *twice*, to
  `src/generated/*.cuh` and to SPIR-V. The Vulkan shaders reach the shared math
  by src-relative path (`#include "shaders/densify.slang"`, with `-I<src>`),
  matching the C++ convention and disambiguating `densify.slang`, which exists
  in both directories.

  Rejected: mirroring the CUDA corpus under `src/backend/cuda/kernels/`. The
  `<Name>.cuh` declaration headers are the *backend-neutral contract* (Vulkan
  TUs include them; they parse with no CUDA toolkit) and they are generated
  from the `<Name>.cu` sitting next to them. A full mirror forces a choice
  between putting the neutral contract inside `backend/cuda/` (so Vulkan
  includes headers out of the CUDA backend — incoherent) and separating each
  `.cuh` from the `.cu` it is generated from (breaking the repo's most
  load-bearing codegen invariant). Symmetry is not worth either. The two sides
  are also not peers in fact: CUDA is the reference implementation the contract
  is *extracted from*, and a symmetric tree would assert a peerhood the codegen
  does not have.

  Two bugs surfaced and were fixed on the way:
  - `backend/api/*.h` forwarders were emitting `#include "../../<Name>.cuh"`,
    stale since the phase-4 move. Nothing includes them yet, so nothing caught
    it. They are now src-relative and each one is compile-checked to parse
    under `-DSSPLAT_BACKEND_VULKAN` with no CUDA toolkit.
  - `build_develop.bash` always exited 0: the trailing `libcsrc.so` move
    masked the build's status, so a failed build reported green. Now
    propagated.

- **Phase 5 — done (binding + gate; deletion deferred).** §4.1 dataset parsing.
  `src/bindings/bind_data.cpp` exposes `DatasetParserConfig`, `ParsedDataset`
  and `parse_dataset` / `parse_{colmap,nerfstudio,metashape}_dataset` through
  pybind (numpy arrays, GIL released during the parse).
  `spirulae_splat/modules/native_dataparser.py` is the ~170-line Python client;
  it imports neither torch nor nerfstudio.

  The parsers and `external/miniz.c` moved from the three app targets'
  source lists into the engine library (`cmake/sources.txt`), so `ssplat-train`,
  `ssplat-mesh`, `ssplat-gui`, the Vulkan `csrc_portable` and the Python
  extension now share one build of them instead of four.

  **Gate:** `tests/python/test_dataparser_parity.py` — 4 formats (COLMAP text +
  binary, Nerfstudio, Metashape) x 4 config variants, comparing frame set and
  order, `camera_to_worlds`, intrinsics, distortion, width/height, camera
  model, seed cloud, `train_frame_scale`, `train_to_normalized` and the
  validation split against `modules/dataparser.py`. 18 passed. Fixtures are
  synthesised from a fixed seed (`tests/python/dataset_fixtures.py`) so no
  local dataset is referenced; `SSPLAT_TEST_DATASET` opts into running the
  same comparison on a real one. A cross-format check asserts the four
  fixtures describe the same scene, so the comparison cannot agree on garbage.

  The gate found a **real bug in the native splitter**: `eval_mode="fraction"`
  and `validation_fraction` used `llround` where numpy's
  `np.linspace(0, n-1, k, dtype=int)` *truncates*, so the native trainer
  silently trained on a different frame subset than the Python one for most
  (n, k) pairs — e.g. n=7, k=5 gives [0,1,3,4,6] vs [0,2,3,4,6]. Fixed in
  `DatasetCommon.cpp` (`linspace_indices`).

  **Not done, deliberately:** deleting `dataparser.py` / `colmap_utils.py` /
  `metashape_utils.py` and rewiring `trainer.py`. Per §7.1 that is a separate,
  separately-announced commit; this one only makes it possible and proves it
  safe. The nerfstudio dependency therefore also still stands.

- **Phase 6 — done (binding + renames; deletion deferred).** §4.2 viewer
  server.

  `src/bindings/bind_viewer.cpp` exposes `WebViewer` (the real `ViewerServer`),
  `ViewerRenderConfig`, `PostSplitCameras` and `bake_post_split`. Python drives
  the same HTTP server, render worker and `viewer.html` the CLI and GUI use.

  The render worker never calls into Python: its hooks read atomics plus a
  mutex-guarded progress string that Python *pushes* (`set_step`,
  `set_progress_json`), instead of invoking Python callbacks from a non-GIL
  thread while holding the engine mutex — which would deadlock against a
  training loop that holds the GIL and wants that mutex. The one shared lock is
  reached through `viewer.engine_lock()`, a context manager that releases the
  GIL while blocking.

  `bake_post_split` also lands the cubemap split (fisheye/equisolid -> 5 faces,
  equirect -> 6) that `trainer.py::_setup_cpp_data_manager` had reimplemented.

  **viewer.html is now genuinely single-source**: moved to
  `src/app/webviewer/viewer.html`, embedded into the engine library (so the
  binding serves it too, not just the app binaries) and read from there by the
  legacy Python server. The four meanings are now `viewer/` (standalone WASM),
  `src/app/webviewer/` (the server + the one client), `src/app/gui/` (desktop
  GUI) and the legacy `spirulae_splat/viewer/` — the ambiguity is resolved by
  each one owning a distinct, descriptive path, not by renaming the root
  `viewer/`.

  **`viewer/` -> `web/` was tried and reverted.** The directory name is the
  published GitHub Pages URL (`.../spirulae-splat/viewer/`, linked from the
  README), so renaming it breaks external backlinks. A public URL outranks an
  internal tidy-up; `viewer/` keeps its name.

  **Fixed a hang this phase owns.** `ssplat-train` never exited when the viewer
  was enabled and `--keep-viewer-alive 0`: `HttpServer::stop()` closed the
  listening socket to break `accept()`, which POSIX leaves undefined and Linux
  simply ignores — the accept never returned and the join blocked forever. The
  accept loop now `select()`s with a 100 ms timeout and checks `_running`, and
  `stop()` closes the socket after the join. Verified: 0.63 s clean exit,
  endpoints unaffected.

  **Gate:** `tests/python/test_webviewer.py` (4 tests) — engine setup through
  the bindings, all four endpoints, `engine_lock()`, and a `stop()`-returns +
  port-released assertion that would have hung before the fix.

  **Not done, deliberately:** deleting `spirulae_splat/viewer/{server,
  http_server,render_worker}.py`, porting `annotation.py`, and rewiring
  `ss_viewer.py` / `trainer.py` onto `_C.WebViewer`. Same §7.1 pacing as
  phase 5.
