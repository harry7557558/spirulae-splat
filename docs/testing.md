# Testing

The native parity tests are the suite. The Python one this file used to
describe is gone -- see §3 for what it covered and what now does not.

## 1. Native cross-backend parity tests (the important ones)

`src/backend/tests/*.cpp` — currently 18 tools covering projection (fwd, bwd,
quant-grad), rasterization bwd, tile intersect, warp, FPBO, optimizer (general
+ geometry), densify, per-pixel train, PPISP, bilagrid, multi-scale loss
(`mask_loss_semantics` is self-checking rather than dump-then-compare: it pins
what an image mask means in the loss, in both mask modes and with none),
meshing (activation, LBVH, occupancy/bisection/color, moment raster, the
per-camera samplers and the visibility cull), plus
`backend/tests/engine/` which drives the *real* engine end to end
(render parity, train parity, and two self-checking tools rather than
dump-then-compare: `engine_reset_state` trains the same scene twice across an
`engine_reset()` and the two must land in the same place, and
`gt_bilagrid_sentinels` pins the GT bilateral grids' invariants -- one depth
scalar per camera, and the no-GT sentinels passing through untouched).
`backend/vulkan/tests/` adds 3 Vulkan-only smoke tests (runtime, pipeline,
sort/scan).

**The same source builds under both backends.** Each file only touches
`backend::`, the generated launch declarations, and `Tensor.h`. The workflow
is dump-then-compare:

```bash
# on the CUDA machine
./build_cuda/projection_parity dump ref.bin
# on the target machine / device
./build/projection_parity compare ref.bin
```

Inputs are deterministic, and comparison is tolerance-based — fast-math
exp/sqrt chains legitimately differ across compilers, and borderline-cull
flips change whole rows, so a small allowance for those is built in.

`msloss_parity` splits its reference into two channels. **Tight**: per-pixel
gradients (deterministic given the raw-loss sums, which enter them only through
smooth reduce math), the densification loss map in every mode, equal-shape
`v_ref_depth` / `v_ref_normal` scatters (one atomic per cell), and the quantile
outputs. **Loose**: `LossValues`, the SSIM display scalar, and scaled-GT
scatters, which accumulate atomically in a backend-specific order. Each config
runs twice and compares the second return, because the scalars come back
through a one-iteration-behind async readout on both backends. One expected
mismatch survives: the CUDA SSIM scalar sums over TILE-GRID positions, so an
image whose dims are not a multiple of the tile picks up zero-padded
out-of-image contributions that differ between the 24- and 16-wide tiles
(`ssim_cs` cfg). It is display-only; gradients and loss values are unaffected.

Several tools also take a `*_DUMP_GOT` environment variable
(`FPBO_DUMP_GOT`, `PPISP_DUMP_GOT`, `MSLOSS_DUMP_GOT`, `PWTRAIN_DUMP_GOT`,
`BILAGRID_DUMP_GOT`, `DENSIFY_DUMP_GOT`) to write the *actual* values
alongside the reference, which is how you diff a mismatch numerically instead
of guessing.

### Building them

```bash
# CUDA branch: opt-in
bash build_develop.bash -B build_cuda -DSS_BACKEND=cuda -DSS_BUILD_BACKEND_TESTS=ON
# Vulkan branch: built unconditionally
bash build_develop.bash -DSS_BACKEND=vulkan
```

Each `.cpp` becomes an executable of the same base name in the build dir.

### Cross-machine / cross-vendor runs

The comparison target is often a different machine (e.g. an AMD GPU box), and
often offline. The pattern that works:

1. Transfer a matching `slangc` to the target and point `-DSS_SLANGC=` at
   it — SPIR-V is compiled at build time and never committed, so the target
   needs a compiler, and the version is pinned.
2. Dump references on the CUDA host.
3. Copy the `.bin` files over and run `compare` on the target.

Keep reference dumps out of git (`parity_refs/` is gitignored).

### macOS / MoltenVK

All 17 parity tools pass against a CUDA reference on Apple silicon, at
essentially the Linux numbers (`engine_render_parity`'s blit channel: 0.157%
of bytes on macOS against 0.152% on Linux, cap 0.2%).

`engine_render_parity` used to fail here on that channel, and the failure was
worth more than its number: the viewer's grid and frustum lines came out
*fragmented on macOS only*. It was not antialiasing, which is what this
section claimed for a while. `vis_blit`'s BVH descent read the popped node
inside a two-iteration child loop, and SPIRV-Cross re-materialized that
threadgroup read once per child instead of keeping it — so the second child
was read after the first child's push had overwritten the slot, and the
descent walked into the wrong subtree. The `[ForceUnroll]` on both child
loops is what keeps the pop a value; see the first MoltenVK rule in
`src/backend/vulkan/README.md`.

Three tools -- `msloss_parity`, `optimgeo_parity`, `meshing_parity` -- pass
only because `VulkanContext::init()` turns MoltenVK's default Metal fast-math
off. `SS_VK_FAST_MATH=1` puts it back, and they fail again; that is the knob
to reach for when measuring what the setting costs.

Speed is a separate question from parity, and macOS answers it differently. Two
kernel choices that cost nothing elsewhere cost an order of magnitude here, and
both are measured at run time rather than assumed: the GEMM tiling
(`OpGemm.cpp`, `SS_NN_GEMM_KERNEL` pins it) and the matcher's dot product
(`integerDotProduct4x8BitPackedUnsignedAccelerated`, `SS_SFM_NO_DOT4` pins it).
With those, an M2 runs SAM 3 image encoding at 12x an RTX 5070 (7.4 s vs 0.64 s,
was 140x) and brute-force matching at 22x (43 vs 1.9 ms per 8192x8192 pair, was
64x); the matcher's remainder is DP4A, which Apple has no instruction for. The
third is "slangc `[unroll]`" in `src/backend/vulkan/README.md`.

## 2. GUI / viewer checks

The web viewer can be driven headlessly over the Chrome DevTools Protocol.
Headless defaults to SwiftShader; to exercise a real GPU, run against a real
display (`DISPLAY=:0`).

A scripted run that serves the viewer needs **`--keep-viewer-alive 0`**, or
the process hangs at exit waiting on it.

## 3. What is gone

The Python suite this document used to describe -- `tests/python/`, the
dataparser and step-config goldens, the trainer and web-viewer gates -- was
deleted with the Python trainer it compared against. There is nothing to run
and nothing to regenerate.

Its job has not gone away, though, and nothing covers it today:

- **The dataset parsers** had a golden over 4 formats x 4 config variants x 2
  splits, checking the frame set, poses, intrinsics, distortion, seed cloud
  and train-frame scalars. A native replacement would generate its fixtures
  from a fixed seed, as that one did, so it needs no dataset on the machine.
  Two of its checks needed no golden at all and are worth rebuilding first:
  the train and eval splits must partition the frames (a bug dropping frames
  from *both* sides leaves each side self-consistent), and every fixture
  format must describe the same scene.
- **`build_step_config()`** had a golden over 8 config variants x 4 run states
  x 20 steps straddling every warmup and decay boundary. Drift in a ported LR
  schedule used to fail there; now it shows up as a quality regression 20k
  steps into a run.

Both are `src/backend/tests/`-shaped work: deterministic input, committed
expectation, one executable. Neither exists yet.

## What to run before calling a change done

| change | gate |
|---|---|
| any kernel | CUDA build + Vulkan build + the relevant parity test on both |
| engine logic | both builds + `engine_render_parity` + `engine_train_step`-level check |
| config field | add the row in `src/config/TrainConfig.h`; check `spirula train --help` and the GUI's All Options editor |
| training-loop logic | `TrainerCore.cpp` — `build_step_config()` is the only place it lives |
| build system | every mode in [build.md](build.md) |
| a comment you wrote | `python3 tools/check_comment_length.py` — the build runs it anyway ([lints](build.md#lints)) |
| anything | one short training run per backend on a public scene |

## Profiling

`SS_PROFILE=1` enables the env-gated per-stage timing breakdown
(H2D / D2H / D2D / memset / device / host). Header-only, works on both
backends — the right first tool when a backend is unexpectedly slow rather
than wrong.
