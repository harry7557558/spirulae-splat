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

## Layout

```
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
map/         Mapper, Bundle, CorrespondenceGraph, Merge, Manager, Profile
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

`--help` on any subcommand lists its flags. `--quality low|medium|high|extreme`
sets the working resolution, the feature cap and the pair-selection breadth;
`--data-type individual|video|internet` sets the pairing mode and how cameras
are grouped. Everything else has a default that a beginner should not have to
touch.

Environment: `SSPLAT_SFM_MAP_PROF=1` prints a mapper stage breakdown,
`SSPLAT_SFM_DUMP_SG` / `SSPLAT_SFM_CMP_STEP` are BA solver debug hooks
(`ba/README.md`).

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
scored against the COLMAP ground truth that ships with it
(`tools/sfm/eval_poses.py`, once the tooling lands in phase 7).

## Not done yet

Inherited from the source tree — none of it is a regression introduced by the
port. Ordered by what blocks the most.

**Scale, past roughly 500-1000 images**

1. **Local BA.** Only global BA runs today. The measurement that decides its
   shape — host dense LM versus a persistent-kernel GPU BA on a 10-30 camera
   problem — was never taken. Take it first.
2. **Shared GPU primitives**: radix sort, prefix scan, segmented reduction, a
   descriptor-set cache, record-once/replay command buffers. The solver
   re-records a command buffer per LM iteration, which is fine at global-BA
   scale and fatal at local-BA scale. Also the clean fix for the SIFT
   extractor leaking buffers when they grow between images.
3. **Gauge fixing and constant-parameter masks**, per-observation weights, and
   mid-solve outlier down-weighting. LM damping regularizes the gauge today.
4. **Vocabulary tree / global-descriptor index** past ~3k images.

**Quality**

5. **Visibility-pyramid next-image scoring** — a raw correspondence count today.
6. **Track merging** in the triangulator. Continuation, creation and completion
   exist; fusing two 3D points bridged by a correspondence does not.
7. **Misregistration on large unordered sets** — images placed in the wrong part
   of the scene. Diagnosed, then parked in favour of throughput work; re-examine
   with item 5.
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
    would cut both CPU time and the peak that sizes the decode pool.
13. **Verification, fewer model fits.** A pair with no real geometry still runs
    every RANSAC trial for both F and H. Do *not* re-attempt SPRT for this: it
    was measured and rejected (residual evaluation is a few percent of RANSAC's
    cost). The win is in not proposing hopeless pairs.
14. **Matcher register-blocking**, if it is still bandwidth-bound. Measure first.

**Unstarted**

15. Learned frontend (ALIKED / SuperPoint + LightGlue) behind the existing
    extractor and matcher interfaces.
16. Global (GLOMAP-style) and hierarchical mappers.
17. Parity benchmarking on ETH3D / IMC.

**Deliberately out of scope**, so they are not silently skipped: rig
constraints in the mapper (a rig is used as a ground-truth-free *diagnostic*,
not a constraint), GPS / geo-registration, MVS / dense reconstruction,
incremental database updates, and relating two models that share no images.
