# Building

`CMakeLists.txt` + `cmake/` builds every artifact: the engine, both
backends, the CLI, the GUI, the meshing tool and the parity tests. There is
no second build system and nothing to keep in sync with one.

The source list lives apart from it, in `cmake/sources.txt`: a plain list of
glob patterns relative to the repo root, read line-by-line by
`cmake/SsSources.cmake`. Adding a source directory means editing that file
once.

Two files that belong to `src/app/` are in that list on purpose:
`app/webviewer/*.cpp` and `app/TrainerCore.cpp`. They live in the engine
library so every app target shares one build of them rather than each
compiling its own copy.

Always build through the dev scripts.

## CMake layout

`CMakeLists.txt` is just the running order; the logic lives in `cmake/`:

| file | what it does |
|---|---|
| `SsOptions.cmake` | options (`SS_BUILD_CLI/GUI`, `SS_DEBUG_SYMBOLS`, …), backend selection, tree paths (`SS_ROOT`, `SS_CSRC`, …) |
| `SsSources.cmake` | `ss_collect_sources()` — expands `sources.txt` with `CONFIGURE_DEPENDS` |
| `SsBackendCuda.cmake` | Torch probe, CUDA arch detection, flags, the `csrc` library, CUDA-side parity tools |
| `SsBackendVulkan.cmake` | portable engine object lib, slangc + SPIR-V embed, `ss_backend_vulkan`, the Vulkan-side tests |
| `SsSlang.cmake` | `ss_find_slangc()` — pinned version, PATH lookup, fetch on miss |
| `SsEmbed.cmake` | `ss_embed_file()` — bake a file into a byte-array header |
| `SsApps.cmake` | the `spirula` executable — every tool the build has, in one binary (backend-agnostic) |

Exactly one backend module runs. It leaves behind `SS_WITH_TORCH` and
`SS_APP_LIBS`, which is the whole contract `SsApps.cmake` depends on.

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
(~750 MB/job) and runs `tools/check_ss_prefix.sh`.

## The matrix

| goal | command |
|---|---|
| CUDA CLI + GUI | `bash build_develop.bash -DSS_BUILD_CLI=ON -DSS_BUILD_GUI=ON -DSS_BACKEND=cuda` |
| Vulkan CLI + GUI | `bash build_develop.bash -DSS_BUILD_CLI=ON -DSS_BUILD_GUI=ON -DSS_BACKEND=vulkan` |
| parity tests | add `-DSS_BUILD_BACKEND_TESTS=ON` (Vulkan builds them unconditionally) |
| Vulkan GUI, everything on | `bash build_develop.bash -DSS_BUILD_CLI=ON -DSS_BUILD_GUI=ON -DSS_BACKEND=vulkan -DSS_ENABLE_PATENTED=ON` |

Keep the two backends in **separate build directories** so you can test both
without reconfiguring, e.g. `-B build_cuda` and `-B build`.

## Options

| option | default | effect |
|---|---|---|
| `SS_BACKEND` | `cuda` | `cuda` \| `vulkan`. `vulkan` builds the portable engine layer + `backend/vulkan/` **without the CUDA toolkit**, and forces `SS_BUILD_CLI=ON`. |
| `SS_BUILD_CLI` | `OFF` | the command-line tools (`spirula train`, `spirula mesh`) |
| `SS_BUILD_GUI` | `OFF` | the graphical application (`spirula` with no arguments); FetchContent's GLFW 3.4 + Dear ImGui v1.92.8 (needs network once) |
| `SS_SEPARATE_TOOLS` | `OFF` | *also* build `spirula-sfm` and `spirula-sam` standalone — same code, but neither links the engine (24 MB vs the combined 61 MB) |
| `SS_BUILD_BACKEND_TESTS` | `OFF` | build `backend/tests/*` (CUDA branch; Vulkan always builds them) |
| `SS_DEBUG_SYMBOLS` | `OFF` | host `-g`, CUDA cubin lineinfo, `slangc -g2`. Bloats binaries substantially — profiling/debugging only. |
| `SS_SLANGC` | *(empty)* | path to a `slangc` to use; empty means find on PATH and fetch the pinned release on miss/mismatch |
| `SS_BUILD_SFM` | `ON` for `vulkan`, `OFF` for `cuda` | `ss_sfm` + `spirula sfm` + `sfm_*_test`. Vulkan-only; a CUDA build can opt in if the Vulkan SDK is present. |
| `SS_BUILD_SAM` | `ON` for `vulkan`, `OFF` for `cuda` | `ss_nn` + `ss_sam` + `spirula sam` + `nn_ops_test` / `sam_pipeline_test`, and the GUI's in-process masking. Same rule as SfM. |
| `SS_ENABLE_PATENTED` | `OFF` | `ss_video` — container demux + `VK_KHR_video_decode_*`. **Read the note below before turning it on.** |
| `SS_SFM_REALS` / `SS_SFM_LOSSES` | all | trim the bundle-adjustment shader variant matrix while iterating (`src/sfm/README.md`) |

### `SS_ENABLE_PATENTED`

Off by default, and deliberately. This repository is GPLv3, and `src/video/` —
the H.264 / H.265 / AV1 bitstream parsers and the Vulkan Video driver — is the
one part of it carrying third-party patent exposure (H.264/H.265 via MPEG LA
and Access Advance, AV1 via the claims asserted against AOMedia). With it off,
that directory is neither compiled nor linked, and everything that wanted it
falls back to an external **ffmpeg**: `spirula sam extract` and `spirula sam
video` say so and exit, and the GUI extracts frames with ffmpeg and tells the
user why.

Turning it on buys in-process GPU decoding: roughly 15× faster frame
extraction (a 127-second 1080p30 clip in ten seconds rather than minutes),
masking that rides along on the same device pass, and no ffmpeg to install.
Nothing else in the build changes.

If you distribute binaries, decide for your jurisdiction and your users before
shipping one built with it on.

## Targets

`spirula` — one executable holding every tool this build has, dispatching on
its first argument (`src/app/Tools.h`):

| | |
|---|---|
| `spirula` | the graphical application |
| `spirula <file-or-folder>` | the same, opening what was named |
| `spirula train` | the CLI trainer |
| `spirula sfm` | structure from motion |
| `spirula sam` | segmentation, tracking and frame extraction |
| `spirula mesh` | meshing (CUDA backend only) |

A copy or symlink named `spirula-sfm` runs that tool directly, so scripts
written against the separate executables keep working; `SS_SEPARATE_TOOLS`
builds those executables for real.

Libraries: `csrc` (engine; SHARED with Torch, STATIC without), `csrc_portable`
(Vulkan branch's torch-free engine objects), `ss_backend_vulkan`,
`ss_sfm`, `ss_nn`, `ss_sam`, `ss_video`, plus one executable
per file in `backend/tests/`, `backend/tests/engine/`, `backend/vulkan/tests/`,
`sfm/tests/`, `nn/tests/` and `sam/tests/`.

The GUI runs reconstruction by re-running **itself** (`spirula sfm auto ...`) as
a child process -- see `src/app/gui/SfmRunner.h` for why that stage and only
that stage is out of process. Masking, decoding and training are in-process.

### Known: `SS_SFM_REALS=df` does not compile on Windows

slangc 2026.12.0.1 hits an internal error on all three `ba_df_*` variants under
Windows (`error[E99998]: Slang compilation aborted due to internal error`),
deterministically; `float` and `double` are fine, and all three build on Linux.
The double-float real is not the default (`--ba-real` defaults to `double`), so
the workaround costs nothing on that platform:

```bat
build_develop.bat "-DSS_SFM_REALS=float;double"
```

`spirula sfm` then reports "variant not built into this binary" if something
asks for `df`. Not yet reduced to a minimal repro or filed upstream.

### Known: `sfm_mask_test` fails on Windows

Its mask-discovery case writes both `e.jpg.jpeg` and `e.jpg.JPEG` and asserts
which one is resolved. On a case-insensitive filesystem those are the same
file, so the assertion cannot hold — and a real Windows dataset cannot contain
both spellings either, so this is the test's assumption failing, not the
lookup. Every other `sfm_*_test`, `nn_ops_test` and `sam_pipeline_test` passes
on Windows.

## Backend specifics

**CUDA.** Needs a recent CUDA toolkit (+ a compatible MSVC on Windows).
Architectures come from `nvidia-smi --query-gpu=compute_cap`; override with
`-DCMAKE_CUDA_ARCHITECTURES`.

**Vulkan.** Needs the Vulkan SDK. No Python and no CUDA toolkit required.
Slang is pinned to a specific version (`SS_SLANG_VERSION` in
`cmake/SsSlang.cmake`); if the `slangc` on PATH doesn't match, CMake fetches the
pinned release. SPIR-V blobs are **never committed** — they are compiled at
build time (one `slangc` edge per blob, see `src/backend/vulkan/shaders/SpirvShaders.cmake`) and
embedded into the binary. On an offline machine, transfer a matching `slangc`
and point `-DSS_SLANGC=` at it.

**Windows.** `build_develop.bat` always calls `vcvars64` even when `cl` is
already on PATH (an ambient `cl`/`INCLUDE` may reference an uninstalled SDK),
picks the newest installed CUDA toolkit unless `CUDA_PATH` is set, falls back
to the VS-bundled CMake/Ninja —
the torch extension build is not supported on Windows, and a broken torch
install aborts configure from inside `TorchConfig.cmake` (which `QUIET` cannot
suppress).

## Build-time cost

The long poles are the biggest `.cu` translation units and the CUDA template
instantiation set in `src/instantiations/` (111 files). Splitting an oversized `.cu`
into function-named parts (see [codegen.md](codegen.md)) both improves
readability and parallelizes the build — no build-file edits needed, since
both builds glob `cmake/sources.txt`. Only `HEADER_SOURCES` in
`generate_headers.py` has to learn the new file.
