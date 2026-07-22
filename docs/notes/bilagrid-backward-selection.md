# Bilagrid backward: adaptive implementation selection (design + plan)

Status: **PROPOSED** (not yet implemented). This document is the plan for
optimizing the bilateral-grid backward, which profiling flags as a consistent
top cost across NVIDIA, both AMD drivers (amdvlk/RADV), Intel iGPU, and
llvmpipe (garden `images_4`, 10 iters). The hot entry is the grid-gradient
kernel `bilagrid_ppisp_bwd_v1_grid` (and its per-family siblings).

## 1. The two backward strategies (and why neither wins everywhere)

The grid gradient `v_bilagrid` is a many-pixels-to-one-cell reduction. There
are two opposite ways to compute it:

- **v1 — gather (thread per grid cell).** Each thread owns one output cell,
  loops over exactly the pixels whose trilinear footprint touches it,
  accumulates the cell's per-channel gradient in registers, and writes it out
  **once** (block-reduce + one atomic per cell). Atomic traffic is O(cells).
  Cost is *redundant work*: neighboring cells' pixel footprints overlap, so
  each pixel's transform-backward math runs ~4x, and small grids force the
  `mult_x/mult_y` tiling that splits a cell across blocks.
  This is the codebase's current always-on path (`*_bwd_v1_grid` + the
  cheap thread-per-pixel `*_bwd_v1_rgb` for the image gradient).

- **v2 — scatter (thread per image pixel).** Each pixel visits its 8 trilinear
  corners x C channels and `atomicAdd`s directly into the grid, computing the
  image gradient in the same fused pass. No redundant work, one kernel, but
  atomic traffic is O(pixels x 8 x C) and collides hard when many pixels map to
  the same cell.
  This exists today only as the **disabled** `#if 0` CUDA
  `BilagridUniformSampleBwdV2_kernel.cuh` (affine) and has **no** Vulkan port.

**Crossover:** v2 wins at low resolution / few pixels per cell (few atomic
collisions, single fused kernel); v1 wins at high resolution / many pixels per
cell (avoids the atomic storm). The exact crossover depends on resolution,
grid shape `(L,gH,gW)`, family (channel count), and the device's atomic
throughput + L2 behavior + driver. No single choice is best across the
hardware matrix — hence a *selection* problem, not a *rewrite* problem.

## 2. Why not the two existing approaches verbatim

- **fused-bilagrid's offline calibration** (`calibration.py`): sweeps ~5740
  timed launches over resolutions x grid shapes x 41 launch presets, fits a
  log-linear cost model per hardware, and ships a checked-in results file
  keyed to the GPU it was tuned on. This is minutes of warmup per machine and
  a per-GPU artifact to maintain — exactly the warmup cost we want to avoid.

- **vksplat's Thompson-sampling scheduler** (`scheduler.{h,cpp}`): a Gaussian
  bandit over kernel variants with measured GPU time as the cost, made elastic
  by exponential forgetting (`adapt_tau`), warmup-noise annealing
  (`warmup_tau`), forced arm init, and forced re-exploration. Correct shape,
  but **context-free**: one flat arm array, no resolution keying. For bilagrid
  the optimum is essentially a *function of resolution*, so a context-free
  bandit would average over a run's resolutions and settle on a compromise (or
  thrash when resolution varies per iteration).

## 3. Proposed design: online **contextual** selection, no offline phase

Combine the two: reuse vksplat's per-arm statistics + annealing math, but make
the bandit **contextual**, keyed on the low-dimensional context that actually
determines the timing.

### Context key
`(family, H, W, L, gH, gW)`. In a run this collapses to a handful of distinct
keys (grid shape is fixed per family; resolution is single-valued per iteration
in this engine — `camera.{height,width}` for RGB, GT depth/normal resolution
for those). Because the reward is nearly a deterministic function of the key
plus small noise (the user's observation), each key's ranking is stable and
converges in a *few* samples.

### Arms (curated + tunable)
Arms are `(impl, target_tile_size)` pairs listed in one editable table
(`default_arms()` in `BilagridBwdSelector.h`). The current set is **5 arms**:
`v1_gather` at `target_tile_size ∈ {2, 4, 8, 16}` (a geometric spread —
`target_tile_size` feeds `mult_x/mult_y`, so small = aggressive footprint split
/ more parallelism + atomics, large = collapses to the no-split path / minimal
atomics + max serial work; the values separate further at higher resolution)
plus `v2_scatter` (thread-per-pixel, `target_tile_size` unused). Editing the
table to add/drop arms after profiling on real hardware is the intended tuning
knob. Exploration stays cheap: with K arms and R distinct resolutions, per-run
exploration is ~min_count·K·R backward passes, each still a *valid* gradient —
a few dozen mildly-slower steps amortized into training, **no separate
calibration phase, no warmup wall-clock the user has to wait through**. This is
the central win over fused-bilagrid.

### Selector (new component, backend-neutral)
`BilagridBwdSelector` — a small header-only class, reused by all five families
and BOTH backends (it is pure host bookkeeping, independent of CUDA/Slang):

```
map<ContextKey, vector<ArmStat>> table;   // ArmStat = vksplat's Stat: count, sum, sum2, num_no_update
ArmId  sample(key);                        // Gaussian TS over per-arm mean cost, keyed by `key`
void   update(key, arm, gpu_seconds);      // exponential forgetting + warmup-noise decay (per key)
```

- Port `Stat` (running mean/var with forgetting), `sample()` (draw one Gaussian
  per arm, pick min cost), forced init (`count < 2`), warmup-noise anneal, and
  the `max_no_update` staleness guard directly from vksplat — but index all of
  it by `key`. New keys start fresh; each pays its own tiny warmup.
- Because reward ~ deterministic per key, the default knobs should lean toward
  **explore-then-commit** (fast exploitation): a handful of samples per arm per
  key, then lock the winner with only occasional re-checks (the `max_no_update`
  probe). The full elastic machinery is retained but tuned aggressive; it
  degrades gracefully to greedy. (Optional later: seed a new key's priors from
  the nearest already-measured resolution to skip its warmup entirely.)
- **Overrides for reproducibility / parity / debugging:**
  `SSPLAT_BILAGRID_BWD = auto | v1 | v2 | v1:<tile> | #<arm-index>` pins the
  arm (`v1` = first v1 arm, `v1:8` = the tile-8 v1 arm, `#3` = raw index); a
  fixed RNG seed (as vksplat's `seed=42`) makes `auto` deterministic given
  identical timings. Reward unit is GPU **milliseconds**
  (`backend::event_elapsed_ms`).

### Measurement (async, no added stall)
Per-dispatch GPU time via `backend::Event` timestamp queries around the
backward dispatch, read back **one iteration behind** (the codebase already
uses this pattern for the AsyncReadout loss/SSIM display), attributed to the
`(key, arm)` stashed with the event. Mirrors vksplat's "measure after the
fence, then `update()`" loop. Up to three timers per iteration (RGB / depth /
normal are separate dispatches with separate keys). Zero extra synchronization
on the hot path: harvesting happens at the top of the backward hook, so every
pending event is >= 1 iteration old and its `event_synchronize` is effectively
free.

Implemented as `BwdTimingRing` (`BilagridBwdTiming.h`): one reused
start/end event pair per channel (RGB/depth/normal), `begin()`/`end()` bracket
the dispatch, `harvest()` completes prior-iteration measurements into a
`vector<BwdMeasurement>`. It brackets the WHOLE family backward (v1's grid+rgb
kernels, or v2's single fused kernel) so v1 vs v2 are compared on total cost.
**Vulkan caveat for Phase 5:** `event_record` flushes the stream, so always-on
timing adds 2 submits per timed channel per iteration on Vulkan; if that shows
up in profiling, time only every Nth iteration (the reward is near-stationary
per key, so sparse sampling is fine) rather than every step.

## 4. Correctness / stability notes

- v1 and v2 differ at the ULP level from atomic summation order (v1 already
  uses atomics on its fast path; cross-backend parity is already
  tolerance-based). Gradients feed Adam, which is robust to ULP noise.
  Switching arms mid-run injects only that noise; slow forgetting
  (`adapt_tau ~ 1000`) prevents oscillation.
- `bilagrid_parity` must test **each arm pinned** (`SSPLAT_BILAGRID_BWD`), so
  both strategies stay bit-parity-checked against their CUDA references.
- The grid-fold tail rule and null-`v_in` skip (see the Vulkan README slice 6)
  apply to the v2 kernels too.

## 5. Maintainability refactor (do this alongside, not after)

Today each family's v1 grid kernel is copy-pasted (5x CUDA headers x
uniform/patched, mirrored into 5 Slang files) differing only in channel count
and the transform-backward math. The v2 revival is the moment to collapse this:

- A single **family descriptor** `{ num_channels, needs_image_grad }` drives
  generic v1 *and* v2 kernels (channel count as a template arg / Slang
  spec-constant; `needs_image_grad` as a flag/spec-axis). affine/ppisp/loglinear
  set `needs_image_grad = true`; **depth/normal set it false** — their v2 then
  scatters *only* the grid gradient and skips the image-grad accumulation
  entirely (even cheaper than for RGB, and mirrors the existing `v_in=nullptr`
  skip on v1).
- The per-family files shrink to just the transform-backward math + param
  packing; the tiling/writeback/selection skeleton lives once.
- `target_tile_size` (currently hardcoded `5` at every call in
  `_engine_bilagrid_backward_hook`) becomes selector-owned, so its tuning is
  data-driven rather than a magic constant.

## 6. Implementation phases (each independently landable + parity-tested)

1. **Selector component** — `BilagridBwdSelector` (contextual Gaussian TS
   ported from vksplat) + `SSPLAT_BILAGRID_BWD` override + unit test with
   synthetic timings. No kernel changes. **DONE** — header-only
   `csrc/BilagridBwdSelector.h`; host unit test
   `csrc/backend/tests/bilagrid_selector_test.cpp` (auto-globbed; also builds
   standalone) verifies forced-init coverage, contextual separation (two
   resolutions → two optima), elasticity across a mid-run reward flip, and all
   override forms. ~97% tail exploitation on the synthetic model.
2. **Async timing plumbing** — per-dispatch `backend::Event` timing, one-iter
   readback, `(key, arm)` attribution at the backward hook. **DONE** —
   `csrc/BilagridBwdTiming.h` (`BwdTimingRing`) wired into
   `_engine_bilagrid_backward_hook` (harvest at hook entry; begin/end around the
   RGB/depth/normal dispatches). Gated by `SSPLAT_BILAGRID_PROFILE=1` (logs
   per-family GPU ms; zero overhead when unset). Validated on RTX 4080S / CUDA:
   20-step garden run reports PPISP RGB backward steady at ~0.94 ms (grid+rgb),
   matching the nsys profile; one-iteration-behind harvest, no hot-path stall.
3. **Revive + generalize v2 (CUDA)** — un-`#if 0` affine v2, refactor to the
   `{channels, needs_image_grad}` generic form, instantiate for all five
   families. Parity per arm. **PPISP DONE** (bottleneck + default-preset RGB
   family + hardest nonlinear case): `BilagridPpispUniformSampleBwdV2_kernel.cuh`
   (thread-per-pixel; reuses the v1 `_rgb` transform-backward verbatim and adds
   the corner scatter; templated on `NEEDS_IMAGE_GRAD`), launcher
   `bilagrid_ppisp_uniform_sample_backward_v2`, declared in `BilagridBindings.h`.
   Standalone CUDA equivalence check (scratchpad `ppisp_v2_equiv.cu`, links
   `libcsrc.a`): grid-grad v1-vs-v2 rel **2.2e-7** (atomic-order rounding),
   image-grad **bit-identical**, and v1 tile=2 vs tile=8 **bit-identical**
   (confirms the tile arms are output-invariant). REMAINING: affine (revive the
   existing v2), loglinear (same shape as ppisp), depth (`scalars` param) and
   normal (axis-angle; both `NEEDS_IMAGE_GRAD=false` -> scatter only, no image
   path) — same mechanical `_rgb`-body + scatter pattern.
4. **Port v2 to Vulkan (Slang)** — `bilagrid_*_bwd_v2` thread-per-pixel scatter
   via `atomic_add_f32`, `needs_image_grad` spec axis. Parity vs CUDA per arm.
   **PPISP DONE**: `bilagrid_ppisp.bilagrid_ppisp_bwd_v2` + `BpV2Params` +
   Vulkan launcher `bilagrid_ppisp_uniform_sample_backward_v2` (Bilagrid.cpp).
   Trains correctly on NVIDIA/AMD Vulkan (identical SSIM to v1).
5. **Wire the selector** into `_engine_bilagrid_backward_hook`: replace the
   three hardcoded `*_backward_v1(..., target_tile_size=5, ...)` calls with
   `sample(key)` -> dispatch chosen arm -> record timing -> `update`. **DONE**
   for PPISP RGB (affine/loglinear keep fixed v1 until their v2 lands);
   file-static selector + `BgSelectorDumper` (prints the per-arm table at exit
   under SSPLAT_BILAGRID_PROFILE). `off` restores the fixed-v1 baseline.

## Benchmark results (garden, 2026-07, PPISP RGB backward, per-arm GPU ms)

The optimal v1 tile is BOTH resolution- and hardware-dependent, and the fixed
`tile=5` default is catastrophic at low resolution (it collapses to `mult=1` ->
one thread serially scans a cell's whole ~thousands-of-pixels footprint). v2's
atomic scatter is catastrophic on AMD (correctly never chosen). Best arm and
gain vs the old fixed tile=5:

| resolution        | RTX4080S (CUDA)     | RX7800XT amdvlk     | RX7800XT RADV       |
|-------------------|---------------------|---------------------|---------------------|
| 648x420           | tile=3  0.31 (10.7x)| tile=3  0.67 (5.2x) | tile=3  0.67 (5.5x) |
| 1296x840          | tile=5  0.94 (1.0x) | tile=5  1.78 (1.0x) | tile=5  1.83 (1.0x) |
| 2593x1680         | tile=5  3.65 (1.0x) | tile=8  6.47 (1.16x)| tile=8  6.61 (1.12x)|

- **Optimum is often a non-power-of-2** (tile=3 wins at low res everywhere).
- **Hardware-dependent**: 4080S wants tile=5 at high res, AMD wants tile=8 --
  so a single checked-in calibration would be wrong somewhere; the online
  per-(hardware,resolution) selector is the point.
- v2 (scatter) never won here (26/81/384 ms on amdvlk) but is retained cheaply
  (see tuning) for regimes with few pixels-per-cell.
- **End-to-end**: RX7800XT amdvlk, 648x420, 500 iters: baseline (fixed tile=5)
  4.1 s vs auto 2.8 s = **1.46x whole-step speedup** -- the kernel saving
  (3.5 -> 0.67 ms/iter) dominates the timing overhead.

### Critical bug found + fixed (v2 grid-grad scatter target)
v1 READS grid values from the `g_id = grid_indices[ni]` slot but WRITES the
grad into a PER-BATCH buffer indexed by `ni` (`out_idx = ni*L*H*W*C`, size
[C_batch,...]; a later step reduces per-batch -> per-grid). The first v2 draft
scattered to the `g_id` base -> with batch size 1 (`ni=0`, `g_id` up to #cams)
that is a wild out-of-bounds atomic -> "operation not supported on global/shared
address space" a few steps in, and garbage gradients. Fix: read from
`grid_base=g_id*...`, scatter to `out_base=ni*...`. The standalone equivalence
test missed it (used `grid_indices=nullptr` so `g_id==ni`); it now runs a
non-identity `grid_indices` with `N_batch < N_grids`.

### Selector tuning (from the benchmark)
Catastrophic arms (v2 on AMD, ~10-60x the best) were being re-probed forever by
the staleness rule (~3 ms/iter waste at 2593x1680). Added a **competitive
re-probe gate**: staleness only re-probes arms whose trusted mean is within
`reprobe_factor` (8x) of the best; clearly-inferior arms are measured once
during forced init and then left alone. `max_no_update` 128 -> 256. Result: v2
drops to n=2 (forced-init only) and never re-probes; competitive arms still
re-probe so a genuine ranking change is still caught.

### Hardware variance: v2 is device-dependent, and a device heuristic is a trap
PPISP v2 vs the best v1 tile (garden 1296x840): NVIDIA/AMD discrete v2 LOSES;
**Intel iGPU v2 LOSES (~2.25x)** but **llvmpipe (CPU) v2 WINS (~2.2x)** -- v1's
gather leans on fast shared memory + cheap redundant ALU (great on real GPUs,
bad on a CPU with no shared memory), while v2's scatter leans on atomics (fine
on a serial CPU, contended on a GPU). A "discrete vs not" arm gate is WRONG:
Intel iGPU is non-discrete yet v2 loses there, and Apple Silicon can't be
predicted (untestable). The online selector already handles this by
MEASUREMENT -- it picked v2 on llvmpipe and v1 on Intel/discrete -- so v2 is
kept as a PPISP arm and the choice is left to the data, not a hardware rule.

Two device-agnostic tunings support keeping v2 everywhere cheaply:
- **`max_no_update` 256 -> 1500** (re-probe far less often). The reward is
  near-stationary per (device, resolution), so re-probing is mostly drift
  insurance. This matters where a losing arm is still "competitive" (Intel v2 is
  2.25x the best -> within `reprobe_factor` -> would be re-probed): at 1500 an
  Intel v2 re-probe (~3 s) amortizes to ~2 ms/iter.
- **Discard the FIRST measurement of each arm.** The first dispatch of a kernel
  variant pays a one-time pipeline-compile / kernel-load cost (Intel iGPU: first
  v2 ~6900 ms vs ~2000 ms steady) that would otherwise inflate the mean and can
  mislead the choice. `update()` drops the first sample per (key, arm) and
  forced-init keeps sampling until it has `min_count` *warm* samples.
6. **Consolidation refactor** — collapse the 5x v1 copy-paste onto the generic
   kernels; unify launchers; extend `bilagrid_parity` to sweep arms; optionally
   add a second `target_tile_size` arm. **PARTIAL:**
   - **Per-type arm sets** (`ppisp_arms()` / `affine_arms()` in the selector;
     the engine has one selector per family via `_bg_selector_for`). PPISP is
     **tile-only** `{2,3,4,5,6,8}` (v2 never won -> dropped, saving its warmup);
     affine is `{2,3,4,5,6,8}+v2` (v2 can beat v1 for the cheap LINEAR transform).
     Tiles 4 and 6 added, 13 dropped (never won). This immediately paid off:
     **affine @1296x840 on AMD converges to tile=4** (0.78 ms) -- beats tile=3
     (0.81) and the old fixed tile=5 (1.20 ms, 1.5x); no prior set had tile=4.
   - **Affine v2 wired** (`bilagrid_uniform_sample_backward_v2` on both backends)
     with the same `grid_indices` read=g_id/write=ni fix as ppisp (the existing
     affine v2 had the latent OOB bug + a params-layout with no grid_indices; the
     Vulkan public launcher didn't pass it). Equivalence test (12-ch): grid-grad
     rel 3.1e-7, image-grad 1.5e-7; trains on NVIDIA+AMD Vulkan and CUDA.
   - **loglinear/depth/normal DONE**: all five families are now selector-driven
     (per-family selectors via `_bg_selector_for`). loglinear got a full v2
     scatter kernel (`BilagridLoglinearSampleBwdV2_kernel.cuh` + Slang
     `bilagrid_loglinear_bwd_v2`; exp-diagonal math, image-grad on), verified
     (grid 3.6e-7 / img 1.0e-7) but v2 LOSES for loglinear (RGB nonlinear, like
     ppisp: 8 ms NV / 173 ms AMD vs tile ~0.6-0.9 ms), so it is DROPPED from the
     default arm set (kernel kept; re-add `{V2Scatter,-1}` to `loglinear_arms()`).
     depth (2-ch) + normal (3-ch)
     got real scatter-only v2 kernels (CUDA `BilagridDepth/NormalSampleBwdV2_kernel.cuh`
     + Slang `bilagrid_depth_bwd_v2` / `bilagrid_normal_bwd_v2`;
     `needs_image_grad=false` so no input-grad/gz tail; normal uses the
     hand-written Rodrigues `aa_rotate_bwd`). Equivalence (grid-grad, with
     grid_indices): depth rel 5.2e-7, normal rel 4.4e-7; tile-invariance
     bit-identical. Verified in-engine on the **kitti** dataset (depth/normal GT)
     across CUDA + NVIDIA/AMD Vulkan: depth converges to tile=4, normal to
     tile=3; v2 loses badly for these TINY geometry grids (4x8x8 -> ~2800 px/cell
     atomic contention, 4-6 ms vs 0.09-0.2 ms), so v2 is DROPPED from the depth/
     normal default arm sets (kernels/launchers kept; re-add `{V2Scatter,-1}` to
     `depth_arms()`/`normal_arms()` to re-enable). Net: v2 is a default arm only
     for affine (linear, cheap-per-pixel); it loses for ppisp/loglinear (RGB
     nonlinear) and depth/normal (tiny grids).
   - **Pre-existing hang found + fixed**: enabling depth on Vulkan exposed a
     latent bug -- the SEPARATE `launch_depth_bwd_v1` (Bilagrid.cpp) dispatched
     the `bilagrid_depth_bwd_v1_depth` input-grad kernel UNCONDITIONALLY, so with
     the engine's `v_depth=nullptr` it wrote N*h*w floats through a null device
     address -> device-lost, never-returning semaphore wait. (The shared
     `launch_family_bwd_v1` guards this for ppisp/loglinear/normal; the depth
     launcher was missing it, latent because depth bilagrid had never run on
     Vulkan.) Fixed with an `if (v_depth != nullptr)` guard. Localized via
     `SSPLAT_VK_DEBUG_SYNC=1` (the last un-`ok`'d dispatch named the kernel).
   - **Arm-set tuning validated**: the added tile=4 wins for affine@1296x840 and
     depth on AMD; tile=6 wins for ppisp@1296x840 on AMD -- both previously
     unreachable optima.
   - **DEFERRED** (noted): the full 5-family generic-kernel refactor, and adding
     v2 arms to the committed cross-backend `bilagrid_parity`.

## 7. Files this touches

- `csrc/BilagridUniformSampleBwdV2_kernel.cuh` (+ per-family v2 headers) —
  revive/generalize.
- `csrc/Bilagrid*UniformSample.cu`, `csrc/EngineBilagrid.cpp` (backward hook) —
  launch + selection wiring.
- `csrc/BilagridConfig.cuh` — shared block/atomic helpers.
- `slang/vulkan/bilagrid_{common,affine,ppisp,loglinear,depth,normal}.slang` —
  add v2 entries.
- `csrc/backend/vulkan/kernels/Bilagrid.cpp` — v2 dispatch in `launch_family_*`.
- `csrc/backend/tests/bilagrid_parity.cpp` — per-arm parity sweep.
- **New:** `csrc/BilagridBwdSelector.h` (the selector).
```
