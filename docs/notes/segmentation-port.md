# Porting the segmentation stack into this repo

Record of moving the standalone `ssam` tree (a Vulkan + Slang port of SAM 2 /
SAM 3) into Spirula Studio, and of what the GUI now does with it.

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
| masks | SAM 2 / SAM 3 (`src/sam/`) | `python scripts/mask.py` | CUDA builds, where `SS_BUILD_SAM` is off; a user who already has lang-sam tuned |

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
  `SS_NN_LOG`, `SSAM_DEBUG_SYNC` → `SS_NN_DEBUG_SYNC`. Device selection
  and validation deliberately took the *engine's* existing names —
  `SS_VK_DEVICE`, `SS_VK_VALIDATION`, `SS_PROFILE` — so one knob
  moves both subsystems, which is what a user expects and what the eventual
  one-device convergence needs.
- The vendored `stb_image` was dropped for `src/external/`'s (same version),
  and image I/O stopped instantiating it — the repository has one impl TU.
- The two apps became `spirula sam` subcommands (`segment`, `track`, `video`,
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
  SAM 2 checkpoint understands: one object per thing, scrub to another frame to
  correct one, and the clicks are carried into the run rather than being a
  preview toy. The model is dropped when the panel closes, so a 2 GB backbone
  is not sitting in VRAM during reconstruction.
- **Model cache** (`ModelCache`) — checkpoints are fetched at run time into the
  cache directory, never bundled, after a one-time per-family consent dialog.
  SAM 2.1 is Apache-2.0 and says so; SAM 3 is under Meta's own licence, which
  is *not* open source and not GPLv3-compatible, so it gets a tick box. The
  wording is three sentences on purpose — a wall of legal text is read by
  nobody, which is the outcome the requirement exists to avoid.
- **Advanced** — an "extra flags" field passed to `spirula sfm auto` verbatim,
  and explicit "use ffmpeg" / "use the Python script" overrides, so an expert
  is never boxed in by what the panel chose to surface.

### Why SfM runs as a child process

`SfmRunner` re-runs **this executable** as `spirula sfm auto ...` rather than
calling the library. This is a deliberate current state, not leftover COLMAP
shape:

1. `src/sfm/` is still a CLI at heart — ~270 `printf` sites and no cancellation
   token (`sfm-port-plan.md` phase 3). In-process it could neither be stopped
   nor reported on.
2. Global bundle adjustment on a large model and a live trainer must not share
   a VRAM budget. A child process gives that separation for free: every byte it
   held is gone when it exits.
3. It keeps **one** Vulkan device live in the GUI process instead of two — the
   port plan's own §10 risk.

The user still installs nothing, and since the merge to a single executable
there is nothing to find either: the child is the same file, invoked with a
different first argument (`src/app/Tools.h`). When phase 3 lands, only
`SfmRunner::run`'s body changes.

The cost is that progress has to be read out of the child's stdout
(`SfmRunner::note_progress`), which is grubby and will break if those lines
change format. They are the ones every stage prints on every run, so it is a
known, contained bet.

## 4. Verified

On this machine (RTX 5070 Laptop, Intel RPL-S, llvmpipe), Vulkan build:

- `nn_ops_test` (30 checks) and `sam_pipeline_test` pass on all three devices.
- `spirula sam segment` with a text prompt and with a click, on a real SAM 3 and
  a real SAM 2.1 checkpoint.
- `spirula sam extract` end to end, with and without masking.
- GUI: photos → built-in SfM → dataset → trainer (25/25 images registered,
  26.4k points, opens and previews).
- GUI: video → GPU decode → masked frames → dataset.
- GUI on Intel: the video panel reports why GPU decoding is unavailable and
  extracts with ffmpeg instead; the run completes.

## 5. Follow-up, 2026-08-01

Five things came back from using it.

**The preview panel could deadlock.** Its worker cleared `_busy` at the end of
the function while every failure left through `return set_error(...)`, so the
first early exit stranded the flag at true and `start_job()` refused to run
again — a panel that never showed anything, whatever you typed. An empty prompt
was enough to trigger it, which is the state the panel *opens* in. Now an RAII
guard clears it, and the frame is published as soon as it is decoded, before
the model is consulted: you see the picture first, and can click on it, which
is the only prompt a SAM 2 checkpoint takes.

**The polarity now labels the prompts.** "Keep only what I name" turns "What to
remove" into "What to keep" and "...but keep" into "...but remove", in both the
dataset screen and the preview, with the radio moved above them since it
decides what they mean. `scripts/mask.py` grew the matching `--keep_prompted`:
the external path had been silently ignoring the polarity and always removing
what the prompt named.

**VRAM was never given back** — see `src/sam/README.md`, "VRAM lifetime". A
closed preview panel left 2.4 GB resident for the life of the GUI; measured at
2418 MiB before and 11 MiB after.

**Two trackers shared one memory bank.** `SlotAllocator` was per-tracker while
`VramPool` is process-wide, so `"people; cars"` — one tracker per phrase — had
both numbering their slots from zero and writing over each other. Keys now
carry the tracker id.

**Masking is not slower in the GUI than in the CLI**, which was the report.
Measured, `extract`-with-masking and `track` over the same frames come out
within 3% of each other; the cost is one SAM 3 backbone pass per frame and
decode is 1% of it. What was missing was a way to know that: the log said
"20 frames written" every tenth frame and nothing about the rate. It now
reports a rate and an estimate, anchored at the first frame so the checkpoint
upload does not skew it.

## 6. Second follow-up, 2026-08-01: clicked objects, and where the time is

**Clicks are prompts for the run now, not just for the preview.** They were
collected by `SegmentPanel` and thrown away when it closed; `PrepJob` carried
only text, so a SAM 2 checkpoint — which has no text tower — could be selected
in the GUI and then failed the whole dataset run with "this checkpoint has no
text encoder". A click now survives to `MaskOptions::seeds`, through both the
extraction path and the folder-masking path.

**One prompt was holding every click.** Every point went into a single
`VisualPrompt`, so pointing at two things asked SAM for one object covering
both and got a mask fitting neither. Clicks now carry an object id and a frame;
`src/sam/README.md` has the model. The CLI spells them `--object` and
`--at-frame`, and `sam_pipeline_test` covers arrival-on-the-right-frame,
refine-does-not-duplicate, and the seeding frame's own mask.

**Re-prompting reloaded the checkpoint.** `Masker::init` is how the preview
changes a prompt, and it called `Session::loadModel` every time, which re-read
and re-uploaded the file — seconds of stall per click on a 700 MB checkpoint.
`loadModel` now returns immediately when the request names the model already
resident.

**A third of a masked frame was CPU with the GPU idle**: ~75 ms to encode a
1080p mask PNG through stb's deflate and ~27 ms to decode the next JPEG,
against ~275 ms of model. Both are off the critical path now
(`src/app/WriterPool.h`, plus a one-frame read-ahead in the CLI); end-to-end
went 384 → 279 ms/frame on Hiera-T.

**The remaining gap to PyTorch is arithmetic throughput, not algorithm.** The
report was that `scripts/SAM2-GUI` tracks at ~10 fps because SAM's memory makes
video cheap, and that this treats frames as independent images. The memory bank
is there and has been since the port; what is missing is bf16 tensor cores.
Measured, profiled and written up in `src/sam/README.md` under "Speed",
including a negative result on the one cheap hypothesis (more flash-decoding
splits: flat). `VK_KHR_cooperative_matrix` is the 3–5× that is left, it is
available on the hardware and in the toolchain, and it is a project of its own.

## 7. Third follow-up, 2026-08-01: tensor cores

`VK_KHR_cooperative_matrix`, which the section above called the 3–5x that was
left. It is in, behind a probe, and it is worth ~1.4x end to end rather than
3–5x. The measurements and the two wrong turns are in `src/nn/README.md` under
"Cooperative matrix"; the short version is that tensor cores are ~5x the fp32
pipe here, so the staging path that was free next to slow arithmetic became the
whole cost next to fast arithmetic, and the first two tilings measured *slower*
than the fp32 kernel they replace. What fixed the GEMM was not tiling at all:
the weight matrix read column-major with stride K already *is* the fragment the
hardware wants, in the fp16 the checkpoint holds, so it never goes through
shared memory.

Attention took only half. `O = O * corr + P @ V` rescales its accumulator by a
per-query factor on every key tile, and the KHR extension does not expose which
element of a fragment a lane holds, so that rescale cannot be written. `Q @ K^T`
has no such problem — the softmax already reads the scores out of shared memory
by (row, column) — so that half moved and the rest did not.

Safety: the kernels live in their own SPIR-V modules, because capabilities are
per module and `vkCreateShaderModule` may reject a whole blob; `Pipelines`
creates the module lazily, so a device without the extension never sees those
words. `SS_NN_COOPMAT=0` forces the fp32 path, and `nn_ops_test` checks
*both* paths in one process against references rounded the way each path rounds,
so the tolerances did not have to be widened. Masks moved by 3 pixels in 13 M.

Also here: SAM 2.1 Base+ and Small joined the catalog. Four sizes across ~2x in
speed is a choice a user can actually make, and two of them are the ones most
captures want.

## 8. Not done yet

1. **One Vulkan device.** Three contexts can exist in one process. The GUI
   sequences them; that is a schedule, not a guarantee. Same work as
   `sfm-port-plan.md` phase 6, and it should be done once for both.
2. **In-process SfM** — phase 3 of the SfM port plan, which is what removes the
   stdout parsing above.
3. **The other half of attention.** `P @ V` is still scalar, for the reason in
   section 7. `VK_NV_cooperative_matrix2` exposes the per-element mapping and
   would close it; a portable alternative is a two-pass softmax that finds the
   row maximum before accumulating, at 1.5x the arithmetic. Neither has been
   costed against the ~35% of a Hiera-T frame it would touch.
4. **Cooperative matrix on a 64-wide subgroup is untested.** The kernels are
   written for it and the launcher passes the subgroup count, but there is no
   AMD or Intel Arc part here to run it on. A wrong answer would show up in
   `nn_ops_test` immediately; a slow one would not show up at all.
5. **The preview segments a still; the run tracks.** A click shows what it
   selects on the frame it was made on, which is honest but is not what the run
   will produce three hundred frames later. Propagating a few frames in the
   panel would close that, at a few seconds per attempt.
6. **A photo folder with clicks is tracked, without asking.** Clicks force
   video mode because a click means nothing on another frame without a memory
   bank to carry it — right for an ordered walk-around, wrong for an unordered
   collection. Nothing detects which one it has.
7. **Checksums for downloaded checkpoints.** The size floor catches a truncated
   file; it does not catch a corrupted one.
8. **`scripts/mask.py` and `scripts/extract_frames.py`** are still the
   standalone Python tools with their own users, now duplicated in kind by
   `spirula sam`. Revisit once the native path has run on enough captures to be
   the obvious default.
