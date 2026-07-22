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
- Duplicated Python/C++ subsystems are being collapsed onto the C++ side.
  Don't add a second implementation of anything.
- `nerfstudio` / `gsplat` are **not** runtime dependencies of the engine or
  the native apps. A few Python modules and `scripts/` still import
  nerfstudio; that is legacy surface, not something to extend.

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
                              (build-time exception: compiles src/data/parsers/*.cpp
                               in place)
spirulae_splat/             Python package — a *client* of the engine, on its way out
├── ss_{trainer,benchmark,viewer,meshing}.py   console-script entry points
├── modules/                config dataclasses, training driver, eval metrics, resume
├── viewer/                 Python HTTP viewer server + viewer.html (html is shared
│                             with the C++ viewer, which embeds it at build time)
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
├── backend/                the backend seam — READ backend/README.md
│   ├── api/                backend-neutral launch declarations (GENERATED forwarders)
│   ├── cuda/  common/      CUDA runtime shim, SortScan, Profiler
│   ├── vulkan/             the whole Vulkan backend: runtime + kernels/ +
│   │                         shaders/ (its own entry points and SPIR-V build)
│   │                         — READ backend/vulkan/README.md
│   └── tests/              native cross-backend parity tests
├── shaders/                Slang device math SHARED by both backends — compiled
│                             twice, to src/generated/*.cuh and to SPIR-V
├── app/                    the native applications
│   ├── cli/                main.cpp (ssplat-train), mesh_main.cpp (ssplat-mesh)
│   ├── gui/                Dear ImGui desktop app (ssplat-gui)
│   ├── webviewer/          HTTP server + render worker for the embedded viewer
│   └── TrainerCore.{h,cpp} the CLI/GUI training loop
├── bindings/ext.cpp        the pybind11 module (144 m.def's)
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
`ssplat-train`'s `$ORIGIN` lookup expects.

Backends build into different trees; keep them separate (`-B build_cuda`,
`-B build`) so you can test both without reconfiguring. Options:
`SSPLAT_BACKEND` (`cuda`|`vulkan`), `SSPLAT_BUILD_CLI`, `SSPLAT_BUILD_GUI`,
`SSPLAT_NO_TORCH`, `SSPLAT_BUILD_BACKEND_TESTS`, `SSPLAT_DEBUG_SYMBOLS`.
Full matrix and per-platform notes: `docs/build.md`.

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
- **`SSPLAT_PROFILE=1`** enables the per-stage backend timing breakdown
  (H2D / D2H / D2D / memset / device / host), header-only, both backends.
- **The engine is a process-global singleton.** Call `engine_reset()` between
  runs that swap datasets, or the new run inherits the old splats, camera
  table, optimizer moments and color-space matrices.
- **`build/_deps` on WSL drvfs** may need `git config safe.directory` entries
  for the FetchContent'd GLFW/imgui.

## Do not commit

Local dataset paths, private capture names, personal directory layouts
(`/mnt/...`, `/home/<user>/...`, `C:\Users\...`). Standard academic dataset
names (e.g. Mip-NeRF 360, ZipNeRF) are fine. Run
`bash scripts/check_private_paths.sh` before pushing.
