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
CMakeLists.txt              all native builds (723 lines, 4 modes — see docs/build.md)
setup.py / pyproject.toml   the pip/torch-extension build (independent of CMake)
build_develop.bash/.bat     the dev build entry points — USE THESE, not pip install
AGENTS.md  CLAUDE.md        this file (CLAUDE.md is a pointer to it)
docs/                       architecture, build, backends, codegen, testing, notes
scripts/                    dataset preprocessing CLI tools (Python, standalone)
tests/                      Python tests (may be stale; see docs/testing.md)
viewer/                     standalone WebGL2 + WASM viewer, independent of training
                              (build-time exception: compiles csrc/app/*Parser.cpp in place)
spirulae_splat/             Python package
├── ss_{trainer,benchmark,viewer,meshing}.py   console-script entry points
├── generate_*.py           the five codegen tools (see "Codegen" below)
├── modules/                config dataclasses, training driver, eval metrics, resume
├── viewer/                 Python HTTP viewer server + viewer.html (html is shared
│                             with the C++ viewer, which embeds it at build time)
└── splat/cuda/
    ├── _backend.py _wrapper*.py   the extension import + lazy function wrappers
    ├── slang/              Slang device math, shared by CUDA and Vulkan
    │   └── vulkan/         Vulkan-only compute entry points
    ├── ins/                GENERATED kernel instantiations (111 files)
    └── csrc/               ALL native code (see below)
```

`spirulae_splat/splat/cuda/csrc/` is where nearly everything lives. The path
is vestigial — `splat/` is a gsplat leftover, `cuda/` also holds the Slang
shaders and the Vulkan backend, and `csrc/` holds the engine, the dataset
parsers, the CLI, the GUI and the web viewer. Inside it:

```
csrc/
├── Engine*.cpp/.h          the torch-free training engine (process-global singleton)
├── <Kernel>.cu             kernel launchers + __global__ kernels
├── <Kernel>.cuh            GENERATED declaration section (see "Codegen")
├── <Kernel>_kernel.cuh     the device-side kernel body, split out for reuse
├── Bilagrid*               43 files — bilateral grid; see docs/notes/
├── Primitive*.cuh          3DGS / Mip / 3DGUT primitive traits
├── DataManager.{cpp,h}     image cache / prefetch / warp pipeline
├── Tensor.h Camera.h Common.cuh   core types and device helpers
├── Mesh*  Delaunay3D.*     meshing pipeline
├── ext.cpp                 the pybind11 module (144 m.def's)
├── app/                    ssplat-train CLI, dataset parsers, web viewer, gui/
├── backend/                the backend seam — READ backend/README.md
│   ├── api/                backend-neutral launch declarations (GENERATED forwarders)
│   ├── cuda/  common/      CUDA runtime shim, SortScan, Profiler
│   └── vulkan/             Vulkan runtime + kernels/ — READ backend/vulkan/README.md
├── generated/  external/   do not hand-edit / vendored
└── tests/ backend/tests/   native parity tests
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

Generated trees are marked in `.gitattributes` (`csrc/generated/`,
`csrc/app/generated/`, `cuda/ins/`) and are **committed**, so a fresh checkout
builds with no Python at all. Five generators, all run from the repo root:

| generator | reads | writes |
|---|---|---|
| `generate_headers.py` | `/*[AutoHeaderGeneratorExport]*/` markers in `csrc/*.cu` | the declaration section of `csrc/<Name>.cuh` |
| `generate_kernel_instantiation.py` | kernel decls in `csrc/*.cuh` | `cuda/ins/*.cu` |
| `generate_cli_config.py` | the Python config dataclasses (`ast`-parsed, no torch import) | `csrc/app/generated/cli_config.h` |
| `generate_backend_api.py` | per-kernel `.cuh` headers | `csrc/backend/api/*.h` forwarders |
| `generate_vulkan_stubs.py` | link-probes the Vulkan build | throwing stubs for unported kernels |

Rules:

1. **Never hand-edit below the `AUTO HEADER GENERATOR — DO NOT EDIT` splitter
   line** in a `.cuh`. Everything above it is yours; everything below is
   regenerated from the `.cu`.
2. To export a launch function to the engine, put
   `/*[AutoHeaderGeneratorExport]*/` immediately above its definition and
   rerun `generate_headers.py`.
3. **A header can be fed by several `.cu` files**, listed explicitly in
   `HEADER_SOURCES` in `generate_headers.py`. This is the sanctioned way to
   split a large `.cu`: name each part after **what it does**, not after the
   header it feeds (e.g. `ImageWarp.cu`, `DepthGeometry.cu` → `PixelWise.cuh`),
   and add it to the list. Both CMake and `setup.py` glob `*.cu`, so no build
   file changes are needed. A listed file that doesn't exist is a hard error,
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

1. the CUDA implementation (`csrc/<Kernel>.cu` + `_kernel.cuh`),
2. the Slang implementation (`cuda/slang/vulkan/*.slang`) and its launcher
   (`csrc/backend/vulkan/kernels/*.cpp`),
3. a parity test in `csrc/backend/tests/` that runs both and compares.

If the Vulkan side isn't ready, `generate_vulkan_stubs.py` emits a throwing
stub so the portable engine still links — that's a deliberate TODO marker, not
a finished state. Coverage is tracked in `csrc/backend/vulkan/README.md`,
which is the authoritative and unusually detailed document on the Vulkan
design (device baseline, capability variants, memory model, atomics). Read it
before touching anything under `backend/vulkan/`.

## Testing

```bash
# native parity tests (CUDA build)
bash build_develop.bash -DSSPLAT_BUILD_BACKEND_TESTS=ON && ./build/<test_name>
# the Vulkan build produces the same test binaries unconditionally
pytest tests/                      # Python tests — currently unmaintained, may fail
```

Each `csrc/backend/tests/*.cpp` becomes an executable of the same name.
`backend/tests/engine/*` drive the real engine end to end. Details and the
CUDA-vs-Vulkan reference-dump workflow: `docs/testing.md`.

## Conventions

- Files are named after **what they do**, not after the file they were split
  out of, and not after the header they declare into.
- `<Name>_kernel.cuh` is reserved: it means "the device body that
  `generate_kernel_instantiation.py` instantiates". Don't use the `_kernel`
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
- `csrc/*.cpp` (as opposed to `.cu`) means "portable, compiles for Vulkan
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
