# Building

There are **two independent build systems**, by design and by history:

- `CMakeLists.txt` — every native artifact (engine, both backends, CLI, GUI,
  meshing tool, parity tests, and optionally the Python extension).
- `setup.py` / `pyproject.toml` — the pip path, using torch's `CUDAExtension`.
  It globs its own source list and carries its own flags. **It does not read
  `CMakeLists.txt`**, so changes that touch compile flags or add source
  directories must be applied to both.

For development, always use CMake via the dev scripts.

## Dev builds

```bash
# Linux
bash build_develop.bash [cmake args...]
# Windows (plain cmd prompt; locates VS via vswhere, calls vcvars64)
build_develop.bat [cmake args...]
```

Both run codegen first (skipped gracefully if `python3` is missing — the
generated files are committed) and then configure + build into `build/`.
`build_develop.bash` additionally caps the job count by available RAM
(~750 MB/job) and, for Torch builds, moves `build/libcsrc.so` to
`spirulae_splat/csrc.so` leaving a symlink behind so `ssplat-train`'s
`$ORIGIN` lookup keeps resolving.

> Do not use `pip install -e .` for Linux development builds.

## The matrix

| goal | command |
|---|---|
| CUDA CLI + GUI | `bash build_develop.bash -DSSPLAT_BUILD_CLI=ON -DSSPLAT_BUILD_GUI=ON -DSSPLAT_BACKEND=cuda` |
| Vulkan CLI + GUI | `bash build_develop.bash -DSSPLAT_BUILD_CLI=ON -DSSPLAT_BUILD_GUI=ON -DSSPLAT_BACKEND=vulkan` |
| CUDA, no Python/Torch | `bash build_develop.bash -DSSPLAT_NO_TORCH=ON` |
| Python extension (pip) | `pip install -e . --no-build-isolation` |
| parity tests | add `-DSSPLAT_BUILD_BACKEND_TESTS=ON` (Vulkan builds them unconditionally) |

Keep the two backends in **separate build directories** so you can test both
without reconfiguring, e.g. `-B build_cuda` and `-B build`.

## Options

| option | default | effect |
|---|---|---|
| `SSPLAT_BACKEND` | `cuda` | `cuda` \| `vulkan`. `vulkan` builds the portable engine layer + `backend/vulkan/` **without the CUDA toolkit**, and forces `SSPLAT_BUILD_CLI=ON`. |
| `SSPLAT_BUILD_CLI` | `OFF` | `ssplat-train` |
| `SSPLAT_BUILD_GUI` | `OFF` | `ssplat-gui`; FetchContent's GLFW 3.4 + Dear ImGui v1.92.8 (needs network once) |
| `SSPLAT_NO_TORCH` | `OFF` | skip Torch/Python even if present; `csrc` becomes a STATIC lib and the exe is self-contained |
| `SSPLAT_BUILD_BACKEND_TESTS` | `OFF` | build `backend/tests/*` (CUDA branch; Vulkan always builds them) |
| `SSPLAT_DEBUG_SYMBOLS` | `OFF` | host `-g`, CUDA cubin lineinfo, `slangc -g2`. Bloats binaries substantially — profiling/debugging only. |
| `SSPLAT_SLANGC` | *(empty)* | path to a `slangc` to use; empty means find on PATH and fetch the pinned release on miss/mismatch |

## Targets

`ssplat-train` (CLI trainer), `ssplat-gui` (native GUI), `ssplat-mesh`
(meshing), `csrc` (engine; SHARED with Torch, STATIC without),
`csrc_portable` (Vulkan branch's torch-free engine objects),
`ssplat_backend_vulkan`, plus one executable per file in
`backend/tests/`, `backend/tests/engine/`, and `backend/vulkan/tests/`.

## Backend specifics

**CUDA.** Needs a recent CUDA toolkit (+ a compatible MSVC on Windows). With
`SSPLAT_NO_TORCH=ON`, architectures come from
`nvidia-smi --query-gpu=compute_cap`; override with `-DTORCH_CUDA_ARCH_LIST`.
The libpython link and the static-libstdc++ / `nftw` interposition workarounds
exist *only* because of libtorch and are skipped in no-torch builds.

**Vulkan.** Needs the Vulkan SDK. No Python and no CUDA toolkit required.
Slang is pinned to a specific version (`SSPLAT_SLANG_VERSION` in
`CMakeLists.txt`); if the `slangc` on PATH doesn't match, CMake fetches the
pinned release. SPIR-V blobs are **never committed** — they are compiled at
build time (one `slangc` edge per blob, see `slang/SpirvShaders.cmake`) and
embedded into the binary. On an offline machine, transfer a matching `slangc`
and point `-DSSPLAT_SLANGC=` at it.

**Windows.** `build_develop.bat` always calls `vcvars64` even when `cl` is
already on PATH (an ambient `cl`/`INCLUDE` may reference an uninstalled SDK),
picks the newest installed CUDA toolkit unless `CUDA_PATH` is set, falls back
to the VS-bundled CMake/Ninja, and configures with `-DSSPLAT_NO_TORCH=ON` —
the torch extension build is not supported on Windows, and a broken torch
install aborts configure from inside `TorchConfig.cmake` (which `QUIET` cannot
suppress). Pass a trailing `-DSSPLAT_NO_TORCH=OFF` to try anyway.

## Build-time cost

The long poles are the biggest `.cu` translation units and the CUDA template
instantiation set in `cuda/ins/` (111 files). Splitting an oversized `.cu`
using the `Name.Part.cu` convention (see [codegen.md](codegen.md)) both
improves readability and parallelizes the build — no build-file edits needed,
since CMake and `setup.py` both glob.
