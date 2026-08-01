# AGENTS.md — orientation for humans and coding agents

Read this first. Detail lives in `docs/`; this file is the map and the rules.

## What this project is

A 3D Gaussian Splatting trainer. It began as a Nerfstudio/gsplat fork and is
being steered toward a **standalone C++ codebase**. The compute engine is
torch-free C++ with two interchangeable backends (CUDA and Vulkan compute via
Slang). PyTorch is no longer on the compute path — the Python package is a
*client* of the engine, kept working for existing users while the dependency
is phased out. **Both backends must keep working, and the Python path must
keep working**, on every change.

Direction of travel, so you don't push the wrong way:

- New functionality goes in C++ first. Python gets a binding, not a port.
- Duplicated Python/C++ subsystems have been collapsed onto the C++ side.
  Don't add a second implementation of anything. Dataset parsing, the viewer
  server and the training driver each exist **once**, in C++, with a Python
  binding: `native_dataparser.py`, `_C.WebViewer`, `_C.TrainerSession`. The
  Python implementations were deleted once a parity gate proved each pair
  agreed; those gates now hold golden values (see
  [docs/testing.md](docs/testing.md) §4-6). Per-step training logic belongs in
  `TrainerCore::build_step_config`, not in `model.py`.
- `nerfstudio` / `gsplat` are **not** dependencies of anything here any more,
  engine or Python. Don't reintroduce them.

## Repo map

```
CMakeLists.txt              running order only; the logic is in cmake/
cmake/                      build modules (options, backends, apps) + sources.txt,
                              the source list BOTH build systems read
setup.py / pyproject.toml   the pip/torch-extension build (own flags; see docs/build.md)
build_develop.bash/.bat     the dev build entry points — USE THESE, not pip install
AGENTS.md  CLAUDE.md        this file (CLAUDE.md is a pointer to it)
docs/                       architecture, build, backends, codegen, testing, notes
src/                        ALL native code (see below)
tools/codegen/              the five codegen tools (see "Codegen" below)
tests/python/               Python tests (may be stale; see docs/testing.md)
tests/native/               standalone native benchmarks
scripts/                    dataset preprocessing CLI tools (Python, standalone)
viewer/                     standalone WebGL2 + WASM viewer, independent of training
                              (published as the online viewer -- the GitHub Pages
                               URL depends on this directory name)
                              (build-time exception: compiles src/data/parsers/*.cpp
                               in place)
spirulae_splat/             Python package — a *client* of the engine, on its way out
├── ss_{trainer,benchmark,viewer,meshing}.py   console-script entry points
├── modules/                what has no C++ counterpart: config dataclasses
│                             (the config source of truth), checkpoint resume,
│                             eval metrics (torch LPIPS/SSIM) + the image
│                             loading eval needs. native_{dataparser,trainer}.py
│                             are the adapters onto the C++ implementations.
│                             camera_utils.py is on NO code path -- it is the
│                             reference for the unported orientation_method /
│                             center_method (docs/notes/pose-normalization.md)
└── splat/cuda/*.py         the extension import + lazy function wrappers
```

`src/` is the include root: every local include is path-qualified relative to
it (`#include "core/Common.cuh"`, `#include "backend/api/BackendTypes.h"`), so
a file's include lines tell you which subsystems it depends on.

```
src/
├── core/                   Tensor.h, Camera.h, Common.cuh, GradQuant.cuh, …
│                             the types and device helpers everything uses
├── primitives/             Primitive*.cuh — 3DGS / Mip / 3DGUT traits
│                             (compile-time types, not runtime branches)
├── kernels/                one directory per family; each holds the launchers
│   │                         (<Name>.cu), the GENERATED declaration header
│   │                         (<Name>.cuh) and the device body (<Name>_kernel.cuh)
│   ├── projection/  raster/  tile/
│   ├── pixelwise/  ppisp/  bilagrid/     (bilagrid is 42 files, see docs/notes/)
│   ├── optim/  densify/  loss/  background/  visualize/
├── engine/                 Engine*.cpp/.h — the torch-free training engine
│                             (process-global singleton)
├── data/                   DataManager (image cache / prefetch / warp) and
│   └── parsers/              COLMAP / Nerfstudio / Metashape readers
├── mesh/                   meshing pipeline, Delaunay3D, UV, export
├── sfm/                    structure from motion -- Vulkan-only, self-contained
│                             (images in, a COLMAP sparse/ model out)
│                             -- READ src/sfm/README.md
├── nn/                     the reusable GPU inference layer: its own Vulkan
│                             runtime (nn/vk/), tensor + ops + Slang kernels,
│                             host image I/O. Knows nothing about any model.
│                             -- READ src/nn/README.md
├── sam/                    SAM 2 / SAM 3 segmentation, on top of nn/
│                             -- READ src/sam/README.md
├── video/                  container demux + VK_KHR_video_decode_*, on top of
│                             nn/. PATENT-GATED: compiled only with
│                             SSPLAT_ENABLE_PATENTED=ON -- READ src/video/README.md
├── backend/                the backend seam — READ backend/README.md
│   ├── api/                backend-neutral launch declarations (GENERATED forwarders)
│   ├── cuda/  common/      CUDA runtime shim, SortScan, Profiler
│   ├── vulkan/             the whole Vulkan backend: runtime + kernels/ +
│   │                         shaders/ (its own entry points and SPIR-V build)
│   │                         — READ backend/vulkan/README.md
│   └── tests/              native cross-backend parity tests
├── shaders/                Slang device math SHARED by both backends — compiled
│                             twice, to src/generated/*.cuh and to SPIR-V
├── app/                    the ONE application, `ssplat` -- every tool below in
│   │                         one executable, dispatched on argv[1]
│   ├── Tools.h  Main.cpp   the subcommand table and the only main()
│   ├── cli/                main.cpp (`ssplat train`), mesh_main.cpp (mesh),
│   │                         sfm_main.cpp (sfm), sam_main.cpp (sam)
│   ├── FrameExtract.{h,cpp}  video -> sharp (optionally masked) frames, shared
│   │                         by `ssplat sam` and the GUI
│   ├── gui/                Dear ImGui desktop app (`ssplat` with no arguments)
│   ├── webviewer/          HTTP server + render worker + viewer.html (the ONE
│   │                         browser client: embedded into the engine library,
│   │                         so CLI, GUI and _C.WebViewer serve the same bytes)
│   └── TrainerCore.{h,cpp} the ONE training driver: config -> dataset ->
│                             seeding -> step loop. CLI, GUI and Python all
│                             drive this; it lives in the engine library
│                             (cmake/sources.txt), not in the app targets.
├── bindings/               the pybind11 module
│   ├── ext.cpp             the bulk of it (144 m.def's)
│   ├── bind_data.cpp       native dataset parsers
│   ├── bind_viewer.cpp     native web-viewer server + post-split bake
│   └── bind_trainer.cpp    SsplatConfig + TrainerSession
├── generated/  app/generated/  instantiations/   GENERATED — do not hand-edit
└── external/               vendored (miniz, stb, npy)
```

## Building

Always use the dev scripts; they run codegen first and pick a sane job count.

```bash
# Linux
bash build_develop.bash -DSSPLAT_BUILD_CLI=ON -DSSPLAT_BUILD_GUI=ON -DSSPLAT_BACKEND=cuda
bash build_develop.bash -DSSPLAT_BUILD_CLI=ON -DSSPLAT_BACKEND=vulkan   # separate build dir advised
# Windows (cmd)
build_develop.bat -DSSPLAT_BUILD_CLI=ON -DSSPLAT_BACKEND=vulkan
```

**Do not use `pip install -e .` for development builds on Linux** — use
`build_develop.bash`. It produces `spirulae_splat/csrc.so` and the symlink
`ssplat`'s `$ORIGIN` lookup expects.

Everything builds into **one executable**, `build/ssplat`: no arguments opens
the GUI, `ssplat sfm|train|sam|mesh` are the command-line tools, and a symlink
named `ssplat-sfm` runs that tool directly (`src/app/Tools.h`). The GUI runs
reconstruction by re-running itself as a child process, so there is no sibling
binary to keep next to it. `-DSSPLAT_SEPARATE_TOOLS=ON` also builds the old
per-tool executables.

Backends build into different trees; keep them separate (`-B build_cuda`,
`-B build`) so you can test both without reconfiguring. Options:
`SSPLAT_BACKEND` (`cuda`|`vulkan`), `SSPLAT_BUILD_CLI`, `SSPLAT_BUILD_GUI`,
`SSPLAT_NO_TORCH`, `SSPLAT_BUILD_BACKEND_TESTS`, `SSPLAT_DEBUG_SYMBOLS`,
`SSPLAT_BUILD_SFM`, `SSPLAT_BUILD_SAM`, `SSPLAT_ENABLE_PATENTED`,
`SSPLAT_SEPARATE_TOOLS`.
Full matrix and per-platform notes: `docs/build.md`.

**`SSPLAT_ENABLE_PATENTED` is OFF by default and should stay that way in
anything you commit.** It gates `src/video/` -- the H.264 / H.265 / AV1
bitstream parsers and the VK_KHR_video_decode_* driver -- which is the only
patent-encumbered code in the tree. With it off, everything that wanted it
shells out to ffmpeg instead; no feature disappears, a subprocess appears. See
the comment on the option in `cmake/SsplatOptions.cmake` before changing it.

## Codegen — the invariants that bite

Generated trees are marked in `.gitattributes` (`src/generated/`,
`src/app/generated/`, `src/instantiations/`) and are **committed**, so a fresh checkout
builds with no Python at all. Five generators, all run from the repo root:

| generator | reads | writes |
|---|---|---|
| `tools/codegen/generate_headers.py` | `/*[AutoHeaderGeneratorExport]*/` markers in `src/kernels/**/*.cu` | the declaration section of the matching `<Name>.cuh` |
| `tools/codegen/generate_kernel_instantiation.py` | kernel decls in `src/kernels/**/*_kernel.cuh` | `src/instantiations/*.cu` |
| `tools/codegen/generate_cli_config.py` | the Python config dataclasses (`ast`-parsed, no torch import) | `src/app/generated/cli_config.h` |
| `tools/codegen/generate_backend_api.py` | per-kernel `.cuh` headers | `src/backend/api/*.h` forwarders |
| `tools/codegen/generate_vulkan_stubs.py` | link-probes the Vulkan build | throwing stubs for unported kernels |

Rules:

1. **Never hand-edit below the `AUTO HEADER GENERATOR — DO NOT EDIT` splitter
   line** in a `.cuh`. Everything above it is yours; everything below is
   regenerated from the `.cu`.
2. To export a launch function to the engine, put
   `/*[AutoHeaderGeneratorExport]*/` immediately above its definition and
   rerun `tools/codegen/generate_headers.py`.
3. **A header can be fed by several `.cu` files**, listed explicitly in
   `HEADER_SOURCES` in `tools/codegen/generate_headers.py`. This is the sanctioned way to
   split a large `.cu`: name each part after **what it does**, not after the
   header it feeds (e.g. `ImageWarp.cu`, `DepthGeometry.cu` → `PixelWise.cuh`),
   and add it to the list. Both build systems glob via `cmake/sources.txt`, so
   no build file changes are needed. A listed file that doesn't exist is a hard error,
   so a rename can't silently drop declarations.
4. The Python config dataclasses are the **single source of truth** for the
   training config. Adding a field there makes it appear in the native CLI,
   the GUI's "All Options" editor, and `--help` automatically after codegen.
   A new field that collides across groups must be listed in `RENAMES`.
5. `.cuh` declaration sections must stay CUDA-include-free — they have to
   parse under `-DSSPLAT_BACKEND_VULKAN` without the CUDA toolkit.

## The Vulkan-only subsystems

`src/sfm/`, `src/nn/`, `src/sam/` and `src/video/` are **not** part of the
two-backend rule below. They are Vulkan + Slang only, carry their own Vulkan
context, share nothing with the training engine, and are absent from a CUDA
build by default (`SSPLAT_BUILD_SFM` / `SSPLAT_BUILD_SAM` default OFF there).
Nothing in them goes through `cmake/sources.txt`, `setup.py` or the pybind
module; the pip build never sees them.

The layering runs one way and must keep doing so:

```
app/gui, app/cli ──► sam ──► nn ──► nn/vk
                 └──► video ──┘
                 └──► sfm  (its own vk/, independent of nn/)
```

`nn/` is the piece meant to be reused. A learned feature detector for SfM, or a
depth/geometry model, goes on top of it unchanged -- so model-specific
constants, weights formats and pipeline policy stay in `sam/` (or the next
`src/<model>/`), never in `nn/`.

## The two-backend rule

Every kernel-level change needs **three** things, or the Vulkan build breaks
or silently diverges:

1. the CUDA implementation (`src/kernels/<family>/<Kernel>.cu` + `_kernel.cuh`),
2. the Slang implementation (`cuda/slang/vulkan/*.slang`) and its launcher
   (`src/backend/vulkan/kernels/*.cpp`),
3. a parity test in `src/backend/tests/` that runs both and compares.

If the Vulkan side isn't ready, `tools/codegen/generate_vulkan_stubs.py` emits a throwing
stub so the portable engine still links — that's a deliberate TODO marker, not
a finished state. Coverage is tracked in `src/backend/vulkan/README.md`,
which is the authoritative and unusually detailed document on the Vulkan
design (device baseline, capability variants, memory model, atomics). Read it
before touching anything under `backend/vulkan/`.

## Testing

```bash
# native parity tests (CUDA build)
bash build_develop.bash -DSSPLAT_BUILD_BACKEND_TESTS=ON && ./build/<test_name>
# the Vulkan build produces the same test binaries unconditionally
pytest tests/python/                      # Python tests — currently unmaintained, may fail
```

Each `src/backend/tests/*.cpp` becomes an executable of the same name.
`backend/tests/engine/*` drive the real engine end to end. Details and the
CUDA-vs-Vulkan reference-dump workflow: `docs/testing.md`.

## Conventions

- Files are named after **what they do**, not after the file they were split
  out of, and not after the header they declare into.
- `<Name>_kernel.cuh` is reserved: it means "the device body that
  `tools/codegen/generate_kernel_instantiation.py` instantiates". Don't use the `_kernel`
  suffix for an ordinary shared-helper header — give it a descriptive name
  (e.g. `BilinearSample.cuh`).
- A `<Family>Common.cuh` holds the shared preamble when a family spans several
  translation units.
- Section banners inside long files:
  ```
  // ================
  // Section Name
  // ================
  ```
  These are load-bearing for the split workflow — keep them meaningful.
- Slang device math is shared by both backends; don't fork it per backend.
- A `.cpp` (as opposed to `.cu`) under `src/` means "portable, compiles for Vulkan
  builds too". Keep the engine layer in `.cpp` and CUDA-free.

## Gotchas worth knowing before you hit them

- **`fs::remove_all` in a torch-linked executable** — never. libtorch's `nftw`
  interposition makes it misbehave; use an explicit recursive delete.
- **Quantized gradient codecs must decode code 0 to exactly `0.0`.** Anything
  else and Adam amplifies the pseudo-gradient into visible floaters.
- **Python test/scripted runs need `--no-keep-viewer-alive`**, or training
  hangs at exit waiting on the viewer.
- **A kernel with one slot per image must fold its grid.** CUDA caps
  `gridDim.y/z` at 65535 and Vulkan caps *every* dispatch dimension at 65535,
  so a per-image axis put on y/z dies somewhere between 8k and 40k images.
  Fold it across two dimensions (`vkk::fold_1d`, or `bilagrid_tv_grid` in
  `kernels/bilagrid/BilagridConfig.cuh` for the 3D case) and pass the fold
  factor to the kernel. Related: the bilagrid samplers index cells with
  int32, which `engine_init_bilagrid_*` enforces up front.
- **`SSPLAT_PROFILE=1`** enables the per-stage backend timing breakdown
  (H2D / D2H / D2D / memset / device / host), header-only, both backends.
- **The engine is a process-global singleton.** Call `engine_reset()` between
  runs that swap datasets, or the new run inherits the old splats, camera
  table, optimizer moments and color-space matrices.
- **`build/_deps` on WSL drvfs** may need `git config safe.directory` entries
  for the FetchContent'd GLFW/imgui.
- **A static library whose only content is a static initializer is not
  linked.** The generated SPIR-V blob tables are archive members that nothing
  references, so registration is an explicit call at each library's entry point
  (`NN_ENSURE_EMBEDDED_MODULES`, `src/nn/vk/EmbeddedSpirv.h`), not a
  namespace-scope constructor. A new shader directory needs its own
  declare/ensure pair or the kernels come back "no shader module".
- **Three Vulkan devices can be live in one process** -- the engine's, the SfM
  module's and the inference layer's. The GUI sequences them (mask preview,
  then reconstruction as a *child process*, then training) so no two overlap;
  do not casually make two of them concurrent. Converging is
  `docs/notes/sfm-port-plan.md` phase 6.
- **Segmentation weights are never committed or bundled.** They are Meta's,
  under Meta's licences, and SAM 3's is not GPLv3-compatible. They are fetched
  at run time after the user has seen the terms -- `src/app/gui/ModelCache.cpp`
  is where that policy lives, and it is the only place that should grow one.
- **The inference layer's VRAM pool is process-wide and grow-only**, so
  destroying a `sam::Session` frees nothing by itself and a 2 GB checkpoint
  stays resident until the process exits. `Session::unload()` (called by the
  destructor) releases it; anything that keys a pool slot per instance must
  namespace the key per owner, as `Tracker`'s memory bank does -- two trackers
  numbering slots from zero write over each other.
- **A GUI worker that clears a `busy` flag at the end of its function will
  strand it.** Every early `return set_error(...)` skips the line, and the next
  request is refused forever. Use a scope guard (`SegmentPanel::start_job`).

## Do not commit

Local dataset paths, private capture names, personal directory layouts
(`/mnt/...`, `/home/<user>/...`, `C:\Users\...`). Standard academic dataset
names (e.g. Mip-NeRF 360, ZipNeRF) are fine. Run
`bash tools/check_private_paths.sh` before pushing.
