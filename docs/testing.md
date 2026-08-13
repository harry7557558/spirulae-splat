# Testing

The native parity tests are the suite. The Python one this file used to
describe is gone -- see §3 for what it covered and what now does not.

## 1. Native cross-backend parity tests (the important ones)

`src/backend/tests/*.cpp` — currently 18 tools covering projection (fwd, bwd,
quant-grad), rasterization bwd, tile intersect, warp, FPBO, optimizer (general
+ geometry), densify, per-pixel train, PPISP, bilagrid, multi-scale loss,
meshing (activation, LBVH, occupancy/bisection/color, moment raster, the
per-camera samplers and the visibility cull), plus
`backend/tests/engine/` which drives the *real* engine end to end
(render parity, train parity). `backend/vulkan/tests/` adds 3 Vulkan-only
smoke tests (runtime, pipeline, sort/scan).

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
| anything | one short training run per backend on a public scene |

## Profiling

`SS_PROFILE=1` enables the env-gated per-stage timing breakdown
(H2D / D2H / D2D / memset / device / host). Header-only, works on both
backends — the right first tool when a backend is unexpectedly slow rather
than wrong.
