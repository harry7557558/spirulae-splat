# Building

There are **two independent build systems**, by design and by history:

- `CMakeLists.txt` + `cmake/` — every native artifact (engine, both backends,
  CLI, GUI, meshing tool, parity tests, and optionally the Python extension).
- `setup.py` / `pyproject.toml` — the pip path, using torch's `CUDAExtension`.
  It carries its own compile/link flags, deliberately: it targets a released
  wheel (its own `WITH_SYMBOLS` / `LINE_INFO` env gates), CMake targets local
  development. **Flag changes must be applied to both.**

They share exactly one thing: the source list, `cmake/sources.txt`. It is a
plain list of glob patterns relative to the repo root, read line-by-line by
`cmake/SsplatSources.cmake` and by `setup.py`'s `get_sources()` — plain text so
the pip build never needs CMake installed. Adding a source directory means
editing that file once.

Two files that belong to `src/app/` are in that list on purpose:
`app/webviewer/*.cpp` and `app/TrainerCore.cpp`. The Python extension binds
both (`bindings/bind_viewer.cpp`, `bindings/bind_trainer.cpp`), so they live
in the engine library and every app target — plus the extension — shares one
build of them, rather than each app compiling its own copy. The consequence
for `setup.py`: it has to reproduce the `viewer.html` byte-array embed that
CMake does (`embed_viewer_html()`), because `Viewer.cpp` includes the
generated header.

For development, always use CMake via the dev scripts.

## CMake layout

`CMakeLists.txt` is just the running order; the logic lives in `cmake/`:

| file | what it does |
|---|---|
| `SsplatOptions.cmake` | options (`SSPLAT_BUILD_CLI/GUI`, `SSPLAT_NO_TORCH`, `SSPLAT_DEBUG_SYMBOLS`, …), backend selection, tree paths (`SSPLAT_ROOT`, `SSPLAT_CSRC`, …) |
| `SsplatSources.cmake` | `ssplat_collect_sources()` — expands `sources.txt` with `CONFIGURE_DEPENDS` |
| `SsplatBackendCuda.cmake` | Torch probe, CUDA arch detection, flags, the `csrc` library, CUDA-side parity tools |
| `SsplatBackendVulkan.cmake` | portable engine object lib, slangc + SPIR-V embed, `ssplat_backend_vulkan`, the Vulkan-side tests |
| `SsplatSlang.cmake` | `ssplat_find_slangc()` — pinned version, PATH lookup, fetch on miss |
| `SsplatEmbed.cmake` | `ssplat_embed_file()` — bake a file into a byte-array header |
| `SsplatApps.cmake` | `ssplat-train`, `ssplat-mesh`, `ssplat-gui` (backend-agnostic) |

Exactly one backend module runs. It leaves behind `SSPLAT_WITH_TORCH` and
`SSPLAT_APP_LIBS`, which is the whole contract `SsplatApps.cmake` depends on.

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
| Vulkan GUI, everything on | `bash build_develop.bash -DSSPLAT_BUILD_CLI=ON -DSSPLAT_BUILD_GUI=ON -DSSPLAT_BACKEND=vulkan -DSSPLAT_ENABLE_PATENTED=ON` |

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
| `SSPLAT_BUILD_SFM` | `ON` for `vulkan`, `OFF` for `cuda` | `ssplat_sfm` + `ssplat-sfm` + `sfm_*_test`. Vulkan-only; a CUDA build can opt in if the Vulkan SDK is present. |
| `SSPLAT_BUILD_SAM` | `ON` for `vulkan`, `OFF` for `cuda` | `ssplat_nn` + `ssplat_sam` + `ssplat-sam` + `nn_ops_test` / `sam_pipeline_test`, and the GUI's in-process masking. Same rule as SfM. |
| `SSPLAT_ENABLE_PATENTED` | `OFF` | `ssplat_video` — container demux + `VK_KHR_video_decode_*`. **Read the note below before turning it on.** |
| `SSPLAT_SFM_REALS` / `SSPLAT_SFM_LOSSES` | all | trim the bundle-adjustment shader variant matrix while iterating (`src/sfm/README.md`) |

### `SSPLAT_ENABLE_PATENTED`

Off by default, and deliberately. This repository is GPLv3, and `src/video/` —
the H.264 / H.265 / AV1 bitstream parsers and the Vulkan Video driver — is the
one part of it carrying third-party patent exposure (H.264/H.265 via MPEG LA
and Access Advance, AV1 via the claims asserted against AOMedia). With it off,
that directory is neither compiled nor linked, and everything that wanted it
falls back to an external **ffmpeg**: `ssplat-sam extract` and `ssplat-sam
video` say so and exit, and the GUI extracts frames with ffmpeg and tells the
user why.

Turning it on buys in-process GPU decoding: roughly 15× faster frame
extraction (a 127-second 1080p30 clip in ten seconds rather than minutes),
masking that rides along on the same device pass, and no ffmpeg to install.
Nothing else in the build changes.

If you distribute binaries, decide for your jurisdiction and your users before
shipping one built with it on.

## Targets

`ssplat-train` (CLI trainer), `ssplat-gui` (native GUI), `ssplat-mesh`
(meshing), `ssplat-sfm` (structure from motion), `ssplat-sam` (segmentation,
tracking and frame extraction), `csrc` (engine; SHARED with Torch, STATIC
without), `csrc_portable` (Vulkan branch's torch-free engine objects),
`ssplat_backend_vulkan`, `ssplat_sfm`, `ssplat_nn`, `ssplat_sam`,
`ssplat_video`, plus one executable per file in `backend/tests/`,
`backend/tests/engine/`, `backend/vulkan/tests/`, `sfm/tests/`, `nn/tests/`
and `sam/tests/`.

`ssplat-gui` looks for `ssplat-sfm` **next to itself** (not on PATH), so a
release has to ship them in one directory. That is also the only thing the GUI
runs out of process; masking, decoding and training are all in-process.

### Known: `SSPLAT_SFM_REALS=df` does not compile on Windows

slangc 2026.12.0.1 hits an internal error on all three `ba_df_*` variants under
Windows (`error[E99998]: Slang compilation aborted due to internal error`),
deterministically; `float` and `double` are fine, and all three build on Linux.
The double-float real is not the default (`--ba-real` defaults to `double`), so
the workaround costs nothing on that platform:

```bat
build_develop.bat "-DSSPLAT_SFM_REALS=float;double"
```

`ssplat-sfm` then reports "variant not built into this binary" if something
asks for `df`. Not yet reduced to a minimal repro or filed upstream.

### Known: `sfm_mask_test` fails on Windows

Its mask-discovery case writes both `e.jpg.jpeg` and `e.jpg.JPEG` and asserts
which one is resolved. On a case-insensitive filesystem those are the same
file, so the assertion cannot hold — and a real Windows dataset cannot contain
both spellings either, so this is the test's assumption failing, not the
lookup. Every other `sfm_*_test`, `nn_ops_test` and `sam_pipeline_test` passes
on Windows.

## Backend specifics

**CUDA.** Needs a recent CUDA toolkit (+ a compatible MSVC on Windows). With
`SSPLAT_NO_TORCH=ON`, architectures come from
`nvidia-smi --query-gpu=compute_cap`; override with `-DTORCH_CUDA_ARCH_LIST`.
The libpython link and the static-libstdc++ / `nftw` interposition workarounds
exist *only* because of libtorch and are skipped in no-torch builds.

**Vulkan.** Needs the Vulkan SDK. No Python and no CUDA toolkit required.
Slang is pinned to a specific version (`SSPLAT_SLANG_VERSION` in
`cmake/SsplatSlang.cmake`); if the `slangc` on PATH doesn't match, CMake fetches the
pinned release. SPIR-V blobs are **never committed** — they are compiled at
build time (one `slangc` edge per blob, see `src/backend/vulkan/shaders/SpirvShaders.cmake`) and
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
instantiation set in `src/instantiations/` (111 files). Splitting an oversized `.cu`
into function-named parts (see [codegen.md](codegen.md)) both improves
readability and parallelizes the build — no build-file edits needed, since
both builds glob `cmake/sources.txt`. Only `HEADER_SOURCES` in
`generate_headers.py` has to learn the new file.
