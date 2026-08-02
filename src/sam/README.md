# Segmentation (`src/sam/`)

SAM 2 / SAM 2.1 / SAM 3 on the inference layer (`src/nn/`): text-prompted
detection, point and box segmentation, and video tracking with a memory bank.

It exists to replace the Python/PyTorch subprocess the GUI used to shell out to
for masking (`scripts/mask.py` with lang-segment-anything, which pulls in a
CUDA PyTorch). What lands in a dataset's `masks/` is the same idea and the same
defaults; what changed is that it runs on the same Vulkan device as everything
else, needs nothing installed, and can be tried on one frame before committing
to a capture.

Imported 2026-08-01 from a separate development tree (a Vulkan + Slang port of
[sam3.cpp](https://github.com/PABannier/sam3.cpp), itself a ggml port of Meta's
SAM 3). The split into `src/nn/`, `src/sam/` and `src/video/` happened on the
way in; that tree is now read-only and this is upstream.

```
                            text "yellow school bus"          SAM 3 only
                                      │              ┌ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ┐
           ┌ ViT, 32 blk ┐            │
   image ──┤             ├─► neck ────┼──► fusion enc ──► DETR dec ──► seg head ──► masks
           └ Hiera, 4 st ┘            │        ▲       └ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ┘
                                      │    text enc + exemplar enc
                                      │
                                      └──► SAM prompt enc ──► SAM mask dec ──► mask
                                                 ▲                  │
                                         memory attention ◄── memory encoder   (video)
```

## Rules

- **Vulkan only**, like `src/sfm/`. There is no CUDA path here and there will
  not be one; the module is built by default only for `SSPLAT_BACKEND=vulkan`
  (`SSPLAT_BUILD_SAM`), and a CUDA build keeps the `scripts/mask.py`
  subprocess.
- **Nothing model-agnostic belongs here.** Tensors, ops, the Vulkan runtime and
  host image I/O are `src/nn/`'s. If you find yourself adding a general kernel
  under `sam/`, it belongs one directory down — that separation is what lets a
  learned feature detector or a depth model reuse the same layer.
- **No heavy dependencies.** Vulkan, Slang, C++17 and the repository's vendored
  `stb_image`. No ggml, no onnxruntime, no PyTorch, no OpenCV.
- **No exceptions across the public API.** `sam/Sam.h` returns `false` / an
  empty result and leaves a description in `Session::lastError()`.
- **Weights are never bundled.** They are Meta's, under Meta's licences, and
  SAM 3's is not compatible with this repository's GPLv3. They are fetched at
  run time into a cache directory after the user has seen the terms — see
  `src/app/gui/ModelCache.{h,cpp}`, which is where that policy lives.

## Layout

```
Sam.h                  the one public header: prompts, Session, Tracker
Common.h               names borrowed from nn/ (see the comment in it)
Masking.{h,cpp}        prompt -> mask POLICY, shared by the CLI and the GUI
MaskIo.cpp             mask + overlay PNG writers
model/                 weights and forward passes for both families
├── Hparams.h          both header formats, and the derived Hiera block table
├── Weights            the .ggml loader (f32/f16/q4_0/q4_1/q8_0)
├── SamModel           the loaded model plus tables worth computing once
├── ImageEncoder       SAM 3's ViT + SimpleFPN necks; the shared entry point
├── HieraEncoder       SAM 2's Hiera trunk + FpnNeck
├── TextEncoder / Tokenizer / Detector / SamDecoder / MemoryModule
pipeline/              Session (encode + PCS + PVS), Tracker, PostProcess
shaders/detr.slang     the SAM-specific kernels; the rest are nn/shaders/
tests/                 sam_pipeline_test -- the whole library over a synthetic
                         checkpoint
```

## The two families

| | SAM 3 | SAM 2 / 2.1 |
|---|---|---|
| backbone | Perception Encoder ViT, 32 blocks at width 1024, stride 14 | Hiera, 4 stages at widths 96–768 (tiny) or 144–1152 (large), stride 16 |
| neck | two SimpleFPNs (detector + tracker), transposed convs | FpnNeck: 1×1 convs, one nearest-neighbour top-down add |
| prompts | text, exemplar boxes, clicks | clicks and boxes only |
| tracker head | identical | identical |
| input | 1008×1008 | 1024×1024, ImageNet mean/std |
| licence | Meta's SAM 3 licence (not open source) | Apache-2.0 |

A SAM 2 checkpoint is visual-only: `segmentConcept` refuses it with a message
saying so, and tracking is seeded with `addInstance`. That is why the GUI's
mask panel offers clicks as well as text — with a SAM 2 model they are the only
prompt there is.

Not ported: EdgeTAM, which sam3.cpp also carries. The backbone seam is one
function (`model::encode_image`), so adding it is bounded work; the loader
rejects unknown files by magic rather than mis-reading them.

## VRAM lifetime

The pool under `src/nn/` is **process-wide and grow-only**, so a destroyed
`Session` frees nothing by itself: a released checkpoint would sit in VRAM
until the process exited. `Session::unload()` releases the weights, the derived
tables and the current frame's features, and the destructor calls it. That
matters because everything that segments in this repository goes on to do
something else on the same GPU -- a reconstruction, a training run -- and 2 GB
is most of a laptop card.

A `Tracker`'s memory bank is its own: slot keys carry the tracker's id in their
high 16 bits, because one masking prompt of `"people; cars"` builds one tracker
per phrase and the pool is shared. `~Tracker` frees its own slots and only its
own.

The GUI takes this one step further and calls `nn::shutdown()` when a dataset
job ends (`DatasetPrep::run`), which destroys the device outright; the mask
preview only unloads the model, since reopening the panel is common and
rebuilding pipelines for the ~50 MB the device holds is a bad trade.

## The mask policy

`Masking.h` is the semantic layer, and it is deliberately the *only* copy: the
CLI, the GUI's dataset run and the GUI's preview all go through it, so what a
user sees in the preview is what gets written.

- several positive phrases, semicolon separated, unioned;
- negative phrases carved back out — a region matching one is KEPT even when it
  also matches a positive phrase;
- the output says what to **keep**, so by default the prompted objects come out
  black and everything else white. `keep_prompted` flips it, for a capture
  where the prompt names the subject rather than a distractor;
- a longest-side cap on what the model sees, with the mask returned at the
  source resolution;
- **clicked objects** (`MaskOptions::seeds`), described below.

The text half matches `scripts/mask.py`, defaults included, so a dataset masked
either way is the same dataset. Clicks have no counterpart there —
lang-segment-anything takes words and nothing else — and the GUI says so rather
than dropping them.

### Clicked objects

A `SeedPrompt` is one object's clicks on one frame, and both halves of that are
load-bearing.

**Objects are separate.** SAM segments a single thing per prompt: one instance
given a click on the dog and a click on the bicycle returns a mask that fits
neither. Each object id becomes its own tracked instance with its own memory
bank, and the masks are unioned at the end — with `resolve_overlaps` deciding
the boundary where two of them claim the same pixel.

**A click belongs to the frame it was drawn on.** Pointing at (900, 500) says
nothing about frame 400 of a moving capture. The first seed for an object
starts the track (`Tracker::addInstance`), and a later seed for the same object
is a correction at the frame where it drifted (`refineInstance`, which replaces
the conditioning memory so the fix sticks). That is SAM 2's conditioning-frame
model, and it is what the reference GUIs expose.

Two details that are easy to get wrong:

- The mask for the seeding frame comes back from `addInstance`/`refineInstance`
  themselves, through their `mask_out` argument. The propagation pass for that
  frame ran *before* the instance existed, so it is not in that `Result`, and
  running the pass again to collect it would advance the frame counter twice
  and write a second memory slot for one frame.
- A seed lands on the first frame at or after the one it names, not on that
  frame exactly. The frame it was drawn on may be one the run never sees: the
  extractor keeps the sharpest frame of each window and drops the rest.

## Using it

```bash
ssplat sam devices
ssplat sam segment --model sam3-q4_0.ggml --image street.jpg --text "school bus" --out out/
ssplat sam segment --model sam3-q4_0.ggml --image cat.jpg --point 315,250 --out out/
ssplat sam track   --model sam3-q4_0.ggml --frames frames/ --out masks/ --text "person; car"
ssplat sam extract clip.mp4 --skip 30 --model sam3-q4_0.ggml --text "person"

# Two clicked objects, the first corrected at frame 90. Works on a SAM 2
# checkpoint, which has no text tower and no other way to be prompted.
ssplat sam track --model sam2.1_hiera_tiny_f16.ggml --frames frames/ --out masks/ \
    --point 640,360 --at-frame 90 --point 700,300 \
    --object --point 120,500
```

```cpp
#include "sam/Sam.h"

sam::Session session;
session.loadModel({.model_path = "sam3-q4_0.ggml"});
session.encodeImage(nn::load_image("street.jpg"));   // once per frame
sam::Result r = session.segmentConcept({.text = "yellow school bus"});
```

Runtime knobs: `SSPLAT_NN_LOG=0..3`, `SSPLAT_VK_DEVICE=<index|name>`,
`SSPLAT_PROFILE=1`, `SSPLAT_VK_VALIDATION=1`, `SSPLAT_NN_DEBUG_SYNC=1`.

## Checkpoints

The format is byte-compatible with sam3.cpp, so its conversion script and every
published model work unchanged. Quantized files (`q4_0`, `q4_1`, `q8_0`) are
dequantized to fp16 during upload: they shrink the *file*, not VRAM.

Measured on the released SAM 3 checkpoint at its native 1008×1008
(`ssplat sam ... --vram` prints it):

```
  weights            weights        1658.5 MiB   (7 chunks)
  neck.pe/det/trk    features        322.8 MiB
  arena              workspace       512.0 MiB (peak 316.6 MiB)
  everything else                     20.3 MiB
  TOTAL                             2513.6 MiB
```

Weights are chunked at 256 MiB rather than allocated as one blob: NVIDIA's
Windows driver refuses a single 1 GiB device-local allocation on an 8 GiB
laptop part with 7 GiB free.

**A Windows footgun.** Every WDDM allocation reserves system *commit*, so
`VK_ERROR_OUT_OF_DEVICE_MEMORY` can mean "the pagefile is small and the commit
charge is exhausted" rather than anything about VRAM — a machine with a 2 GB
pagefile and 15 GB of RAM mostly in use will refuse GPU allocations at around
2.7 GiB no matter which GPU you pick. The error message checks for this and
says so.

## Speed

End-to-end `ssplat sam track` on an RTX 5070 Laptop, one instance through
1920×1080 frames on disk — decode, model and mask PNG included:

| | ms/frame | model only |
|---|---|---|
| SAM 3, q4_0, text prompt | 1471 | 1444 |
| SAM 2.1 Hiera-L, f16 | 604 | 602 |
| SAM 2.1 Hiera-T, f16 | 279 | 273 |
| SAM 2.1 Hiera-T, `--memory-frames 4` | 225 | 223 |
| SAM 2.1 Hiera-T, `--memory-frames 1` | 190 | 188 |

Everything outside "model only" — reading the JPEG, encoding the mask PNG — is
overlapped with the GPU (`app/WriterPool.h` and a one-frame read-ahead). It was
worth doing: those two are ~105 ms a frame of pure CPU, a third of a Hiera-T
frame, and the loop used to sit through them with the device idle.

### Where the model time goes, and why it is what it is

`--profile` with `SSPLAT_NN_LOG=2` breaks it down by kernel. For Hiera-T, 72%
of the GPU time is `flash_attn` and 15% is `gemm_nt_big`; for SAM 3 it is the
other way round, 59% GEMM and 34% attention.

**SAM 2 is memory-attention-bound.** 4096 queries against `num_maskmem × 4096`
memory tokens over a **single** 256-wide head, four layers deep: 480 GFLOP per
frame once the bank is full, about half of the whole frame. `--memory-frames N`
is linear in it, which is why the table above moves the way it does.

**SAM 3 is backbone-bound.** 32 blocks over 5184 tokens at width 1024 is ~90%
of its FLOPs — 830 ms/frame in `gemm_nt_big` alone.

Both are running at the rate those kernels measure at in isolation: ~5 TFLOP/s
for the GEMM, ~3 for the attention, against ~21 TFLOP/s of fp32 peak on this
part. PyTorch's SAM 2 video predictor does ~10 fps with Hiera-L on comparable
hardware (`scripts/SAM2-GUI`) where this does 1.7. **That gap is bf16 tensor
cores, and nothing else.** Confirming that, and what it costs to close:

- Raising the flash-decoding split target from 512 to 1024 and 2048 workgroups
  measured 275 and 279 ms against 279 (that is, flat), so the attention kernel
  is not short of parallelism at this shape. `src/nn/README.md` predicted this;
  it is now measured on a second architecture.
- Nothing algorithmic is missing. The memory bank is there, it stays in VRAM
  across frames, and the K/V projections it feeds are 0.2% of the frame — there
  is no per-frame recomputation left to remove.
- `VK_KHR_cooperative_matrix` is present on this device and Slang emits it, so
  an fp16 tensor-core path for `gemm_nt_big` and `flash_attn` is available and
  is the 3–5× that is left. It is a real project: two kernels, a capability
  query, a fallback for the devices without it (llvmpipe has none), and
  numerical validation against `nn_ops_test`. It has not been attempted.

The levers a user has today are the checkpoint (**SAM 2.1 Tiny is 5× faster
than SAM 3** and is the right default for tracking a clicked object through a
video), `--memory-frames`, the number of frames, and `--mask-mode image`, which
drops the tracker head for 20%. That last one is **not** the default: the
memory bank is what carries an instance through a frame the detector misses,
and losing an object for one frame of a reconstruction costs more than the 20%.
`--img-size` is the biggest lever of all on SAM 2 — 512 instead of 1024 takes
Hiera-T to 64 ms of model time, since memory attention is quadratic in the
token grid — and is not available for SAM 3, whose rotary tables ship sized for
a 72×72 grid.

Two things that are **not** true, both checked because they were reported:

- Masking during extraction is not slower than masking frames already on disk.
  Over the same frames with SAM 3 q4_0: `track` 1520 ms, `extract --mask-mode
  video` 1542, `extract --mask-mode image` 1240. Decoding the video is 1% of a
  frame.
- Tracking a video is not slower than segmenting the frames as stills *because*
  of the tracking. It is slower by the cost of the tracker head, and it buys
  temporal consistency; the memory bank is a fidelity feature here, not a
  speed one. What makes video segmentation fast in the reference
  implementations is that they can skip the detector between frames
  (`--detect-every`), which this does too.

Which is why the GUI reports a rate and an estimate per frame rather than a
bare counter — at a minute per forty frames, "how long" is the only question a
progress line can usefully answer.

Where the kernel time went, and what did *not* work, is
[`src/nn/README.md`](../nn/README.md#what-the-kernels-are-actually-bound-by).

## How it differs from the ggml original

Faithful to the reference numerics; what changed is the execution model.

- **Activations never leave the device.** sam3.cpp routes every inter-stage
  tensor through host `std::vector<float>` because ggml's graph allocator would
  otherwise overwrite live buffers. Here only the final masks, boxes and scores
  come back, as thresholded bytes rather than floats.
- **No permutes.** Feature maps are channel-last and token sequences are
  `[N, heads*dim]`, so a fused qkv projection lands where attention indexes.
- **Fused epilogues** — bias, activation and the residual add ride along.
- **Flash attention.** The score matrix is never materialized.
- **Masks are decoded only for survivors** of the score threshold and NMS.
- **The memory bank stays in VRAM** across frames.
- GELU is the erf form everywhere (the ggml port uses tanh in the memory
  encoder; the PyTorch source is `nn.GELU()` in both places).

## Tests

```bash
./build/nn_ops_test        # every GPU kernel vs an independent CPU reference
./build/sam_pipeline_test  # the whole library over a synthetic checkpoint
```

`sam_pipeline_test` writes real, format-correct checkpoints at a shrunken
configuration — one `sam3` and one `sam2`, including ggml's reversed `ne[]`
order exactly as the conversion scripts write it — and runs concept
segmentation, visual segmentation and video tracking over them. It proves the
plumbing: every weight name resolves, every shape lines up, the arena's scoping
is sound, the memory bank and association logic run. It does not prove
numerical fidelity, since the weights are noise; that is `nn_ops_test`'s job.

Neither substitutes for running a released checkpoint, which is how the
remaining shape and grid-size assumptions were shaken out. `ssplat sam segment`
on a real model, and the GUI's mask preview, are part of the manual checklist
before a change lands.

## Not done yet

- **One Vulkan device.** `nn::vk::Context`, `sfm::VkCtx` and the engine's
  `backend::vk::Context` are three separate devices in one process. The GUI
  sequences them (preview, then reconstruction as a child process, then
  training) so no two are live at once, but that is a schedule, not a
  guarantee. Converging is the same work as `docs/notes/sfm-port-plan.md`
  phase 6 and should be done once, for both.
- **fp16 activations + cooperative matrix**, which is the only remaining lever
  on SAM 2's memory attention. It changes the numerics and needs device
  features outside the Vulkan 1.2 baseline.
- **EdgeTAM**, per above.
- **Interactive refinement in the GUI preview.** `Tracker::refineInstance`
  exists and the CLI can reach it; the preview panel currently re-runs from
  scratch on every click instead of adding conditioning memory.
