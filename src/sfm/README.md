# Structure from Motion (`src/sfm/`)

A standalone SfM pipeline — images in, a COLMAP `sparse/` model out — with a
GPU compute backend (Vulkan + Slang) and no heavy dependencies. It exists to
replace the `colmap` subprocess: `ssplat-sfm` is the CLI, and the native GUI
will drive the same library in-process instead of shelling out
(`docs/notes/sfm-port-plan.md` phase 5).

Imported 2026-07-27 from a separate development tree; that tree is now
read-only and this is upstream. What the port has not yet done is at the bottom
of this file and in the port plan.

## Rules

- **Vulkan only.** There is no CUDA path here and there will not be one. The
  module carries its own Vulkan context (`vk/VkContext.h`) and its own embedded
  SPIR-V, and shares nothing with the training engine — including, for now, the
  device. Built by default only for `SSPLAT_BACKEND=vulkan`.
- **No heavy dependencies.** Vulkan, Slang, C++17, and the repository's
  vendored `stb_image`. No Ceres, no Eigen on the hot path, no OpenCV, no
  SQLite, no PyTorch.
- **Compute on the GPU, control flow on the host.** Slang kernels do the
  per-pixel / per-feature / per-observation work; the host owns graph
  structure, RANSAC bookkeeping and the mapper's decisions. RANSAC in
  particular stays on the host deliberately: thousands of tiny minimal solves
  with data-dependent control flow, and not the bottleneck.
- **Port algorithms, not code.** COLMAP is the reference for behaviour and
  parameter defaults; the implementation is ours, against our own data
  structures. Numerics that had to match are called out in the comments by
  decision number (`docs/notes/sfm-design.md`).
- **COLMAP's on-disk formats are the interchange format**, and every stage is a
  subcommand reading and writing files. Any one stage can be swapped for
  COLMAP's equivalent to bisect a failure.

## Stage graph

```
images/ ──► extract ──► features/ ─┐
                                   ├─► match ──► matches.bin ──► map ──► sparse/0..N
            (pairing: exhaustive / │             (two-view          │
             sequential / GPU      │              verification)     │
             pre-selection) ───────┘                                ▼
                                                    merge / manage rounds:
                                                    merge on shared poses (Sim3)
                                                    ► audit ► grow ► prune
                                                    ► reseed ► split a fold
                                                                    │
                                                                    ▼
                                                            sparse/0 ► training
```

`ssplat-sfm auto` runs all of it from two knobs, `--quality` and `--data-type`.

The `map` box is one incremental reconstruction of the whole capture by default.
`--mapper bottom-up` replaces it with the opposite schedule (`map/Bottomup.h`,
D57): the verified view graph is cut by normalized cut into *atoms* of a few
dozen images (`map/Partition.h`), each is reconstructed by the same `Mapper`
restricted to it, and the atoms are then merged upwards a level at a time —
every model absorbing at most one other per level, with one bundle adjustment
across all of them, intrinsics shared per camera group, between levels. It is a
*schedule*, not a different algorithm: every model is still built by
`Mapper::run`'s own rules and joined by the same Sim(3) merger the manage rounds
use, which is what makes a regression in it impossible to confuse with a
regression in the geometry. What it changes is where the cost and the risk sit —
an atom is too small to get its own focal wrong, and too small for a whole-model
pass to be expensive. `Bottomup.h` carries the numbers.

## Layout

```
SfmConfig.h  the one option surface: the aggregate of every stage's options
  .cpp         plus the descriptor table the CLI, --help and the GUI all read
vk/          VkContext.h   the Vulkan compute context every stage builds on
             EmbeddedSpirv.h  lookup for the SPIR-V compiled into the binary
core/        types shared by every stage, no Vulkan:
               Camera / CameraSetup   camera models, and which images share one
               Pose                   Rigid3, Sim3, angle-axis conversions
               Image / ImageLoader    decode, grayscale, the batch decode pool
               Exif                   focal prior + camera identity from headers
               Features / Matches     the on-disk feature and match formats
               Mask                   keypoint masking, sampled in uv
               Model                  Reconstruction + COLMAP binary IO
feature/     Sift (GPU), Matcher (GPU), Pairing, PairSelection (GPU),
               Verification (host worker pool)
geometry/    Essential, Fundamental, Homography, P3P, AbsolutePose,
               Triangulation, TwoView, LinAlg
optim/       Ransac   LO-RANSAC with MSAC scoring
ba/          Problem (model registry + problem layout), Solver (LM, dense
               Cholesky / implicit-Schur PCG), README.md
map/         Mapper, Bundle, CorrespondenceGraph, Merge, Manager, Profile,
               Partition (view-graph normalized cut), Bottomup (atom
               reconstruction + merge)
shaders/     common/ (Real, df, dmath, linalg, camera, loss), ba/, sift/, match/
tests/       one executable per file
```

`src/` is the include root, so every include is path-qualified:
`#include "sfm/core/Camera.h"`. Note that `sfm/core/Camera.h` and the engine's
`core/Camera.h` are different types with different jobs; the paths keep them
apart.

## Building

```bash
bash build_develop.bash -DSSPLAT_BACKEND=vulkan -DSSPLAT_BUILD_CLI=ON
```

`SSPLAT_BUILD_SFM` defaults ON for the Vulkan backend and OFF for CUDA (where
it can still be turned on if the Vulkan SDK is present). Shaders compile at
build time and are embedded, so nothing needs to sit next to the binary.

The bundle-adjustment kernels are compiled once per (Real, Loss) pair — nine
blobs at the default. Trim the matrix while iterating:

```bash
cmake -B build -DSSPLAT_SFM_REALS=df -DSSPLAT_SFM_LOSSES=trivial
```

Both are cached, so a trimmed value sticks until the full list is passed again
(`-DSSPLAT_SFM_REALS='float;double;df'`) or the build tree is wiped. Asking for
a variant that was trimmed out is a clear runtime error, not a crash.

## Running

```bash
ssplat-sfm auto IMAGES/ -o WORKSPACE/          # images -> sparse model
ssplat-sfm auto -o ws/                         # ./images + ./masks, all defaults
ssplat-sfm auto IMAGES/ -o ws/ --data-type video --quality medium
ssplat-sfm auto IMAGES/ -o ws/ --masks MASKS/  # drop keypoints on masked pixels
ssplat-sfm auto IMAGES/ -o ws/ --camera-model opencv-fisheye

ssplat-sfm extract IMAGES/ -o feats/
ssplat-sfm match   feats/ -o matches.bin
ssplat-sfm map     matches.bin feats/ -o sparse/ --images IMAGES/
ssplat-sfm merge   sparse/ -o merged/
ssplat-sfm ba      problem.txt --real df       # solver benchmark on a BAL problem
```

`ssplat-sfm --help` lists the commands, `ssplat-sfm <command> --help` (or
`ssplat-sfm help <command>`) prints that command's usage, its options with
their defaults and worked examples, and `ssplat-sfm --version` prints the
package version. A usage error names the flag, says what was wrong with it and
points at `--help`; it always exits 1, because `auto` spends exit codes 2 and 3
on *the reconstruction* being absent or partial.

Environment: `SSPLAT_SFM_MAP_PROF=1` prints a mapper stage breakdown,
`SSPLAT_SFM_DUMP_SG` / `SSPLAT_SFM_CMP_STEP` are BA solver debug hooks
(`ba/README.md`).

## Options

`SfmConfig.h` holds the pipeline's whole option surface: the stage option
structs unchanged (`SiftOptions`, `MatchOptions`, `PairSelectionOptions`,
`TwoViewOptions`, `CameraSetupOptions`, `MapperOptions`, `ManagerOptions`,
`MergeOptions`), the pipeline-level knobs that span stages, and
`SFM_CONFIG_FIELDS` — a descriptor table with one row per flag:

```c
F(member, name, cmds, tier, group, lo, hi, choices, help)
```

The table is the single source of truth for the CLI parser, for `--help`, and
(port plan phase 5) for the GUI's options editor; each is one macro expansion
over it in `SfmConfig.cpp`, so a new knob is one row and never three edits.
`cmds` is which subcommands accept the flag — a name may repeat across commands
with disjoint masks, which is how `--max-error` is the verification tolerance
for `auto`/`match`/`map` and the *alignment* tolerance for `merge`. `tier` is
`Basic` (what the GUI shows unfolded: quality, data type, camera model, camera
sharing, focal, masks, pair mode) or `Advanced`; `Alias` marks a second
spelling. A `bool` row is a switch and gets both `--name` and `--no-name`, and
`--help` prints whichever direction changes the default.

Three things stay hand-parsed in `src/app/cli/sfm_main.cpp`, because they do
not name one scalar field: `--camera-model PREFIX=MODEL` and `--focal PREFIX=F`
(which also feed the per-group override list), `--no-manage` (four fields at
once), and the flags that pick what a command *does* rather than how —
`-o/--output`, `map --audit`, `auto --no-masks`. The CLI offers a token to the
hand-parsed cases first, so those names always win.

Pipeline knobs are fanned out into the stage structs by `SfmConfig::finalize()`
and nowhere else, so the CLI and the GUI cannot disagree about what
`--max-error` (one tolerance, two struct fields — D47), `--device`, `--quiet`
or the camera settings mean.

`--quality low|medium|high|extreme` sets the working resolution, the feature cap
and the pair-selection breadth; `--data-type individual|video|internet` sets the
pairing mode, the seed angle and how cameras are grouped. Both are applied
*before* the table's overrides, so naming a flag explicitly always wins over the
preset, and `auto` reports every field a preset moved:

```
  preset    : --max-image-size 3200 -> 1000
  preset    : --max-features 8192 -> 2048
```

Everything else has a default that a beginner should not have to touch.

`--mapper flat|bottom-up|auto` picks the schedule (see the stage graph);
`--bup-atom-size` and `--bup-overlap` size the atoms and the overlap the
merge aligns on. `--merge-tracks` and `--rank-by-visibility` are on by default
and exist to be turned off when attributing a change.

Two things about the surface changed when it was unified, both deliberate:
`auto --no-merge` now disables *merging* only, which is what it already meant on
`map`; skipping the whole merge / grow / prune / reseed loop is `--no-manage`,
which now exists on both commands. And `auto` accepts every advanced flag
`extract`, `match` and `map` accept, since it runs those stages — a run can be
tuned without decomposing it into three commands, and the GUI's editor has one
command's worth of fields to show.

## Tests

Each `tests/*.cpp` builds to an executable of the same name that prints
PASS/FAIL and returns 0/1 — the same convention as `src/backend/tests/`.

| binary | covers | needs a GPU |
|---|---|---|
| `sfm_sift_test` | GPU SIFT, matcher, batch decode, camera-model and format round trips | yes |
| `sfm_map_test` | synthetic reconstruction end to end, incl. manage/audit/split | yes |
| `sfm_cholesky_test` | dense GPU Cholesky vs a CPU reference | yes |
| `sfm_geometry_test` | F, H, E, P3P, triangulation, RANSAC, SVD/eigen kernels | no |
| `sfm_merge_test` | Sim(3) algebra, model alignment, track splicing, fold detection | no |
| `sfm_mask_test` | mask uv sampling, decode, file discovery | no |

End to end, the check that matters is a reconstruction on a public dataset
scored against the reference that ships with it: `tools/sfm/eval_poses.py` reads
a COLMAP model (binary or text), a Nerfstudio `transforms.json` or a Metashape
`.xml` (`tools/sfm/metashape_ref.py`) and reports registration rate, an
alignment-free relative-pose AUC and Sim(3)-aligned absolute errors.

Two things about that metric are worth knowing before reading a number:

- **A pair touching an unregistered image counts as a 180 degree failure**, so
  the AUC is capped by `(registered / total)^2`. On a model whose registered
  poses are all good, AUC *is* the registration rate; do not read it as
  accuracy until coverage is accounted for.
- **The references mostly held the principal point at the image centre**, so a
  model that refines it (D51) is scored against one that absorbed the offset
  into its rotations. That shows up as a uniform relative-rotation error of
  about `dx / f` radians and it is not necessarily the model being worse.

## Not done yet

Inherited from the source tree — none of it is a regression introduced by the
port. Ordered by what blocks the most.

**Scale, past roughly 500-1000 images**

1. **Local BA.** Only global BA runs today. The measurement that decides its
   shape — host dense LM versus a persistent-kernel GPU BA on a 10-30 camera
   problem — was never taken. Take it first. Note what a 1000-image profile
   actually says before assuming this is the win: the mapper's bundle
   adjustments are *mostly small already*, because the model grows from two
   images, and its cost sits in the last few full-size ones. Those are the ones
   local BA does not replace. What made the difference at that size was the
   linear solver's dense/CG crossover (`--ba-solver`) and bounding track length
   (`kMergeMaxTrack`), both of which act on exactly those passes — and, once
   the CG path stopped requiring per-image intrinsics groups, the crossover
   started applying to ordinary captures at all (`ba/README.md`).
2. **Shared GPU primitives**: radix sort, prefix scan, segmented reduction, a
   descriptor-set cache, record-once/replay command buffers. The solver
   re-records a command buffer per LM iteration, which is fine at global-BA
   scale and fatal at local-BA scale. Also the clean fix for the SIFT
   extractor leaking buffers when they grow between images.
3. **Gauge fixing and constant-parameter masks**, per-observation weights, and
   mid-solve outlier down-weighting. LM damping regularizes the gauge today.
4. **Vocabulary tree / global-descriptor index** past ~3k images. Less urgent
   than it was: pair selection is now two-stage (a symmetric mini-vs-mini
   shortlist over every pair, then the reliable asymmetric score on the
   shortlist), which cut the quadratic term 5-6x — 59 s to 10 s on 1194 images —
   while keeping 98-99.5% of the same pair set. The term is still quadratic.

**Quality**

5. ~~Visibility-pyramid next-image scoring~~ — done (`--rank-by-visibility`,
   on by default): the next image is ranked by how its visible structure
   *spreads* over the frame, COLMAP's MIN_UNCERTAINTY default, and an image
   that already failed sorts behind every untried one.
6. ~~Track merging~~ — done (`--merge-tracks`, on by default):
   `Mapper::mergeTracks` fuses two 3D points a correspondence says are the same
   feature, subject to an all-inliers reprojection test, one observation per
   image, the union still subtending the minimum triangulation angle, and a cap
   on the merged track's length — the reduced camera system has an entry per
   image *pair* on a track, so unbounded fusion buys a twentieth observation of
   an already-pinned point and pays for it in the solver.
7. **Misregistration on large unordered sets** — images placed in the wrong part
   of the scene. Diagnosed, then parked in favour of throughput work; the
   visibility ranking of item 5 helped and did not close it.
7b. **The fold split's veto is calibrated on two points.** A real fold's cut
   severs 0.00% of the model's co-visibility (`sfm_merge_test`); the one sound
   model in an 80-dataset corpus that the conflicts talked into a cut severed
   1.30%, and taking it cost 568 images and 65 points of AUC. The veto is now
   0.5%, between them but nearer the fold. It fires on one dataset in eighty,
   so a third data point is worth having before trusting the number.
8. **Nister 5-point** (calibrated init) and **EPnP**.
9. **Automatic camera-model detection.** A fisheye capture run with the default
   rectilinear model reconstructs badly and nothing detects it. The focal
   bootstrap's peripheral inlier curve is a usable signal.

**Integration**

10. **Undistortion stage.** Never written. The dataset parser takes distortion
    parameters, so confirm this is wanted before building it.
11. **Equirectangular end to end.** The mapper writes `EQUIRECTANGULAR`
    (model 17); `ColmapParser` stops at model 10 and the renderer has no
    spherical camera, so such a model cannot currently be trained on.
12. **Faster decode.** A scaled JPEG decode straight to the working resolution
    would cut the CPU time. (The *peak* half of this is done: the decoder
    resamples out of stb's RGB buffer instead of building a full-resolution
    float image, so a concurrent decode costs 3 B per source pixel rather than
    7, and the pool's byte budget now comes from the machine's RAM instead of a
    fixed 1 GiB. Extraction is GPU-bound again on 21 MP inputs.)
13. **Verification, fewer model fits.** A pair with no real geometry still runs
    every RANSAC trial for both F and H. Do *not* re-attempt SPRT for this: it
    was measured and rejected (residual evaluation is a few percent of RANSAC's
    cost). The win is in not proposing hopeless pairs. Three cheaper things
    landed since: an exact early bail in the scoring loop (a model that cannot
    reach the incumbent's inlier count stops being scored), local optimization
    *inside* the trial loop rather than only after it (a better incumbent means
    fewer trials for the same confidence), and a homography residual with no
    transcendental or square root at all.
14. **Matcher register-blocking**, if it is still bandwidth-bound. Measured:
    the win was in the workgroup *width*, not in registers. A pair's train
    descriptors are streamed through groupshared once per workgroup, so the
    traffic is `ceil(nQuery / TQ) * nTrain`; `match_rows` (the pair-selection
    path, which has no column side to constrain its query packing) now runs
    256 wide instead of 64 and reads a scoring pair's train set twice instead
    of eight times. Register-blocking on top of that would halve the
    groupshared reads again at roughly half the occupancy -- untested.

**Unstarted**

15. Learned frontend (ALIKED / SuperPoint + LightGlue) behind the existing
    extractor and matcher interfaces.
16. A **global** (GLOMAP-style) mapper. The **bottom-up** one exists
    (`--mapper bottom-up`, `map/Partition.h` + `map/Bottomup.h`); what it has
    not got is parallel atom reconstruction, which needs a second `rec_` per
    worker and one shared `VkContext`.
17. Parity benchmarking on ETH3D / IMC.

**Deliberately out of scope**, so they are not silently skipped: rig
constraints in the mapper (a rig is used as a ground-truth-free *diagnostic*,
not a constraint), GPS / geo-registration, MVS / dense reconstruction,
incremental database updates, and relating two models that share no images.
