# Porting the segmentation stack into this repo

Record of moving the standalone `ssam` tree (a Vulkan + Slang port of SAM 2 /
SAM 3) into `spirulae-splat`, and of what the GUI now does with it.

Landed 2026-08-01. What is left is at the bottom and in
[`src/sam/README.md`](../../src/sam/README.md).

---

## 1. What it replaces

The GUI's masking used to be a subprocess: an external Python running the
embedded `scripts/mask.py` against `lang-segment-anything`, which needs a CUDA
PyTorch install. That is a hard ask for the audience the GUI is for, and it
does not work at all on a machine without CUDA.

Frame extraction used to be a subprocess too — `ffmpeg` at an oversampled rate
into a temporary folder, then a sharpest-frame pass over the JPEGs
(`FrameSelect`).

Both now have an in-process implementation on the same Vulkan device as
everything else, with the subprocess kept as a fallback rather than deleted.
Neither is a requirement any more; both are still reachable, on purpose:

| stage | built in | fallback | why keep the fallback |
|---|---|---|---|
| frames | `VK_KHR_video_decode_*` (`src/video/`) | ffmpeg | patent gate (off by default); codecs and profiles a driver will not decode; Intel/llvmpipe expose no video queue at all |
| masks | SAM 2 / SAM 3 (`src/sam/`) | `python scripts/mask.py` | CUDA builds, where `SSPLAT_BUILD_SAM` is off; a user who already has lang-sam tuned |

## 2. The three-way split

The incoming tree was one library with `src/{core,vk,nn,model,pipeline,video}`.
It came in as three, because only one third of it is about SAM:

```
src/nn/      core/ vk/ io/ Tensor Ops Op*.cpp shaders/     namespace nn, nn::vk
src/sam/     Sam.h model/ pipeline/ Masking shaders/detr    namespace sam
src/video/   demux + bitstream parsers + Vulkan Video       namespace video
```

`nn/` is the piece meant to be reused — a learned feature detector for SfM
(replacing SIFT), or a depth/geometry model, sits on it unchanged. That is why
model-specific things (weight formats, prompt policy, the DETR kernels) are all
one level up.

`sam/Common.h` and `video/Common.h` carry the names those subsystems borrow
from `nn` (`vk::`, `Error`, `fail`, `parallel_for`, the log macros). Without
them ~400 call sites would each grow an `nn::` that says nothing a reader does
not already know from the directory.

### What changed on the way in

- `ssam::` → `nn::` / `sam::` / `video::`; `SSAM_LOG_*` → `NN_LOG_*`;
  `SSAM_CHECK` → `NN_CHECK`.
- Runtime knobs joined the repository's namespace: `SSAM_LOG` →
  `SSPLAT_NN_LOG`, `SSAM_DEBUG_SYNC` → `SSPLAT_NN_DEBUG_SYNC`. Device selection
  and validation deliberately took the *engine's* existing names —
  `SSPLAT_VK_DEVICE`, `SSPLAT_VK_VALIDATION`, `SSPLAT_PROFILE` — so one knob
  moves both subsystems, which is what a user expects and what the eventual
  one-device convergence needs.
- The vendored `stb_image` was dropped for `src/external/`'s (same version),
  and image I/O stopped instantiating it — the repository has one impl TU.
- The two apps became `ssplat-sam` subcommands (`segment`, `track`, `video`,
  `extract`, `devices`).
- `ssam`'s `tools/spirv_embed.cpp` was dropped for `spirv_tool`'s new
  `embed --nn <tag>` mode, so the repository still has one host tool for SPIR-V.

### The one non-mechanical change

The shader blob registry became **additive**. Three libraries each contribute a
generated table, and `nn::vk::Pipelines` rescans the tail whenever the registry
grows, so one pipeline cache serves all of them.

Registration is an explicit call at each library's entry point
(`NN_ENSURE_EMBEDDED_MODULES`), not a static initializer, for an unglamorous
reason worth remembering: these are static archives, and **an object file that
nothing references is not linked at all**. A namespace-scope constructor in the
generated blob TU silently never runs, and every kernel comes back "no shader
module". Adding a fourth shader directory needs its own declare/ensure pair.

## 3. What the GUI does now

`Screen::Colmap` became `Screen::NewDataset`, and the shared half of
`ColmapRunner` became `DatasetPrep` — frames, sharpest-frame selection, `.insv`
track splitting, masking, and the resume rules. Both dataset paths run it, so
how frames are chosen cannot differ between them.

- **Engine selector** — "Built-in (GPU)" vs "COLMAP (installed separately)",
  shown only when both are actually available. A Vulkan user with no COLMAP
  never learns COLMAP exists; a CUDA user is never offered a back end this
  build does not have.
- **Mask preview** (`SegmentPanel`) — the prompt run on one real frame of the
  real input, through the same `sam::Masker` the dataset run uses, with red
  marking what would be dropped. Clicks work too, which is the only prompt a
  SAM 2 checkpoint understands. The model is dropped when the panel closes, so
  a 2 GB backbone is not sitting in VRAM during reconstruction.
- **Model cache** (`ModelCache`) — checkpoints are fetched at run time into the
  cache directory, never bundled, after a one-time per-family consent dialog.
  SAM 2.1 is Apache-2.0 and says so; SAM 3 is under Meta's own licence, which
  is *not* open source and not GPLv3-compatible, so it gets a tick box. The
  wording is three sentences on purpose — a wall of legal text is read by
  nobody, which is the outcome the requirement exists to avoid.
- **Advanced** — an "extra flags" field passed to `ssplat-sfm auto` verbatim,
  and explicit "use ffmpeg" / "use the Python script" overrides, so an expert
  is never boxed in by what the panel chose to surface.

### Why SfM runs as a child process

`SfmRunner` spawns the sibling `ssplat-sfm` binary rather than calling the
library. This is a deliberate current state, not leftover COLMAP shape:

1. `src/sfm/` is still a CLI at heart — ~270 `printf` sites and no cancellation
   token (`sfm-port-plan.md` phase 3). In-process it could neither be stopped
   nor reported on.
2. Global bundle adjustment on a large model and a live trainer must not share
   a VRAM budget. A child process gives that separation for free: every byte it
   held is gone when it exits.
3. It keeps **one** Vulkan device live in the GUI process instead of two — the
   port plan's own §10 risk.

The user still installs nothing; it is our binary, shipped alongside. When
phase 3 lands, only `SfmRunner::run`'s body changes.

The cost is that progress has to be read out of the child's stdout
(`SfmRunner::note_progress`), which is grubby and will break if those lines
change format. They are the ones every stage prints on every run, so it is a
known, contained bet.

## 4. Verified

On this machine (RTX 5070 Laptop, Intel RPL-S, llvmpipe), Vulkan build:

- `nn_ops_test` (30 checks) and `sam_pipeline_test` pass on all three devices.
- `ssplat-sam segment` with a text prompt and with a click, on a real SAM 3 and
  a real SAM 2.1 checkpoint.
- `ssplat-sam extract` end to end, with and without masking.
- GUI: photos → built-in SfM → dataset → trainer (25/25 images registered,
  26.4k points, opens and previews).
- GUI: video → GPU decode → masked frames → dataset.
- GUI on Intel: the video panel reports why GPU decoding is unavailable and
  extracts with ffmpeg instead; the run completes.

## 5. Not done yet

1. **One Vulkan device.** Three contexts can exist in one process. The GUI
   sequences them; that is a schedule, not a guarantee. Same work as
   `sfm-port-plan.md` phase 6, and it should be done once for both.
2. **In-process SfM** — phase 3 of the SfM port plan, which is what removes the
   stdout parsing above.
3. **Interactive refinement in the preview.** `Tracker::refineInstance` exists;
   the panel re-runs from scratch on every click instead of adding conditioning
   memory.
4. **Masking a photo folder is per-image, not tracked.** Consecutive photos may
   be anywhere in the scene, so a memory bank does not apply — but for a folder
   that *is* an ordered walk-around, tracking would be both faster and steadier.
   Nothing detects which one it has.
5. **Checksums for downloaded checkpoints.** The size floor catches a truncated
   file; it does not catch a corrupted one.
6. **`scripts/mask.py` and `scripts/extract_frames.py`** are still the
   standalone Python tools with their own users, now duplicated in kind by
   `ssplat-sam`. Revisit once the native path has run on enough captures to be
   the obvious default.
