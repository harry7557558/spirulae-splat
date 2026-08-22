# Datasets

## Supported layouts

Three formats, auto-detected in this order: **Nerfstudio**, then **COLMAP**,
then **Metashape**.

The native parsers (`src/data/parsers/*Parser.cpp`) are the one
implementation, shared by the CLI trainer, the GUI and the WASM viewer.

| format | inputs | parser |
|---|---|---|
| COLMAP | `cameras`/`images`/`points3D` in `.bin` or `.txt` | `ColmapParser.cpp` |
| Nerfstudio | `transforms.json` + a PLY point cloud | `NerfstudioParser.cpp` (PLY reader lives here) |
| Metashape | camera-export `.xml` + `.ply`, optionally a `.psx` project for filename disambiguation | `MetashapeParser.cpp` (XML via `app/Xml.h`, zips via `external/miniz`) |

Default subdirectory names: `images/`, `masks/`, `depths/`, `normals/`.
The COLMAP reconstruction directory is auto-detected over
`{sparse/0, colmap/sparse/0, sparse, colmap, .}` unless `recon_dir` is set.

## Camera models

Perspective (with full radial / tangential / thin-prism distortion),
equidistant and equisolid fisheye (including >180° FOV as produced by 360
cameras), and equirectangular/spherical. See `src/core/CameraModel.h` — it is
plain C++17 with no CUDA dependency, which is why the WASM viewer can reuse
it.

COLMAP writes 18 camera models and `ColmapParser.cpp` accepts every one.
`EQUIRECTANGULAR` (id 17) is the spherical one: its params are `(w, h)`
rather than a calibration, because the image *is* the calibration. It reaches
the same `CameraModelType::EQUIRECTANGULAR` as a Metashape `spherical` sensor,
with the same convention — +Z forward at the image centre, azimuth wrapping at
the left/right edge — so the two formats describe an identical camera and
`bake_post_split` treats them identically. Note the engine's canonical panorama
intrinsics assume a 2:1 (360°×180°) image; the parser warns when one is not.

## Lens distortion tiers

Distortion is a separate axis from the camera model, and a COMPILE-TIME one:
a template argument in CUDA, a `kDistortion` specialization constant on Vulkan.
The four tiers, cheapest first, are `None`, `OpenCV` (`k1 k2 p1 p2`),
`ThinPrism` (`k1 k2 k3 k4 p1 p2 sx1 sy1`, COLMAP's `THIN_PRISM_FISHEYE`) and
`Rational` (`k1..k6 p1 p2`, COLMAP's `FULL_OPENCV`, where `k4..k6` divide).
A slot index does NOT mean the same thing across tiers. `core/CameraModel.h`
is the one definition; `shaders/projection_utils.slang` is the one
implementation, generic over an `ICameraDistortion`.

The parser picks the CHEAPEST tier that represents the source camera exactly,
so a PINHOLE dataset costs no distortion registers at all and a `FULL_OPENCV`
camera whose `k4..k6` are zero demotes to `ThinPrism`
(`camera_distortion_demote`). Only eleven (model, tier) pairs are compiled —
no COLMAP fisheye model is rational, and EQUIRECTANGULAR carries no distortion.
`camera_distortion_is_compiled()` is that list, and it must stay in step with
`kCameraVariants` in `tools/codegen/generate_kernel_instantiation.py` and the
export list in `shaders/primitive_3dgs.slang`.

A `transforms.json` is ambiguous about what `k4` means -- OpenCV's first
rational DENOMINATOR term, or Kannala-Brandt's / Metashape's fourth RADIAL one.
An explicit `camera_distortion` key settles it (`MetashapeParser` always writes
one); failing that, a fisheye camera is Kannala-Brandt, `k5`/`k6` mean rational
on their own, and the mere PRESENCE of `b1`/`b2`/`sx1`/`sy1` -- keys a rational
camera never carries -- makes `k4` radial. That last rule is what a
Metashape-converted `transforms.json` needs: reading its `k4` as a denominator
is a different lens, not a small error.

`FOV`, `SIMPLE_DIVISION`, `DIVISION`, `EUCM` and `RAD_TAN_THIN_PRISM_FISHEYE`
have no exact tier, and neither does a Metashape sensor skew (`b2`), which is
an off-diagonal pixel term where every tier's pixel map is diagonal. They are
fitted onto a (model, `ThinPrism`) pair by near-minimax regression
(`data/DistortionFit.h`).

`fit_camera_auto` chooses the camera model from the source's MEASURED field of
view, not from its name: a COLMAP `FOV` lens is perspective at omega 0.3 and a
180-degree fisheye at omega 0.87, and forcing the second onto a pinhole target
puts `tan(85 deg) = 11.4` into a degree-8 polynomial and fails. It then walks a
coefficient ladder (all eight, no thin prism, no `k3`/`k4`, ..., none) until the
fitted distortion is invertible everywhere sampled. **It never fails**: a
dataset that took hours to reconstruct must not refuse to load over a lens
model, so the worst case is a plain fisheye and a warning, not an exception.

A fit closer than `dsfit::kExactFitPx` (0.1 px) is left alone -- the fitted
camera already reproduces the source to better than bilinear resampling can
resolve, so re-distorting would only cost VRAM and blur. In practice that
covers `FOV`, both division models and `EUCM`; `RAD_TAN_THIN_PRISM_FISHEYE` and
a real `b2` do not reach it. When it is not reached, the source model goes into
`ParsedDataset::redistort`, which makes `bake_post_split` set `any_warp` even at
K = 1 so the images route through the warp path's staging and get resampled.

The resampling reads the TRUE source projection (`shaders/camera_source.slang`,
the one place those models live on device; `data/SourceCamera.h` is the host
mirror, and the two must agree exactly) rather than the fit -- going through the
fit would be a no-op. Source model ids 0..17 are COLMAP's own `CameraModelId`
values with COLMAP's parameter array verbatim; ours start at 1000 so COLMAP can
keep appending to its enum.

A source model is sampled only where its image still grows outward as the ray
tilts off axis. Past that it has folded -- at the lens border for a polynomial
fisheye, well inside the frame for the division and unified models -- and two
directions share a pixel, so a warped face would be ringed by a mirrored copy
of the image instead of ending at the lens. Those rays are dropped, and the
synthesized FOV mask drops them by the same test. The fit domain stops earlier
still, where the radial rate falls below a quarter of its on-axis value: at the
fold itself no fitted distortion is invertible, so the coefficient ladder would
degrade a camera the tier otherwise reproduces to a fraction of a pixel.

Two paths, both gathering per destination pixel:

- **K = 1** (`kernels/pixelwise/ImageRedistort.cu`): destination pixel ->
  ray through the fitted camera -> source pixel. The fit leaves the pose alone,
  so a destination pixel and its source pixel are the SAME ray: depth (linear or
  ray) and camera-frame normals transfer unchanged and only the sampling
  coordinate moves. None of GtDepthNormalWarp.cu's point-space handling applies.
- **K > 1** (warp_to_pinhole): the face ray projects STRAIGHT through
  the source camera, so the fitted camera is never materialized and the two
  passes cost one kernel and no intermediate image. `RayToPixel<D, kFromSource>`
  is the seam; `kFromSource` is a template argument (a specialization constant
  on Vulkan) so an ordinary dataset pays neither the branch nor the 16 registers
  the source parameters occupy.

## Two-stage parse

1. **`parse_dataset`** → `ParsedDataset`: per-**input** cameras in the raw
   `train_frame="points"` frame — poses and points exactly as stored. The
   normalized-frame similarity is computed only to obtain the
   `train_frame_scale` scalar.
2. **`bake_post_split`** (`data/PostSplit.cpp`) → the post-split arrays
   `engine_setup_data_manager` consumes. This is either an identity (K=1)
   pass-through, or the split `camhost::plan_split_faces` plans for a wide
   camera when `warp_to_pinhole` is enabled (`warp_spherical_to_pinhole` for a
   panorama).

## The split

A wide camera is rendered as pinhole faces, one per frame of a fixed table:
five around the optical axis for a fisheye (front, +x, +y, -x, -y), the six
cube faces for a panorama. Every face has the focal a 90-degree face of
`ceil(sqrt(W*H/K))` pixels would have -- the density of an uncropped split --
but is **cropped to the rays the lens holds**: the planner rasterizes each
frame's visibility (a valid projection, by the GPU warp's own fold test, that
lands inside the image), takes the bounding box, and covers every box with one
common tile, chosen to minimize the pixel count plus a fixed cost per face. A
frame that holds only a corner sliver (under 4% of its 90 degrees) is dropped.

The fixed cost is what decides the outcome. A face pays a projection of every
splat and a sort whatever its size: measured at 0.5 S^2 pixel-equivalents with
115k splats and 1.2 S^2 with 219k on an RTX 5070 laptop (S the side of a
90-degree face), and growing with the splat count; the planner charges 0.75
S^2. So a 180-degree fisheye, or a cropped one like a 108 x 162 degree frame,
trains on six half-height faces -- the front cut in two and four side bands --
at 60-66% of the pixels five square faces took (measured 17% faster at 115k
splats), while a 200-degree lens whose side bands keep 80% of each face stays
at five square faces: cropping there would save 6% of the pixels for 20% more
faces. Past that, the lever is skipping fully-masked raster tiles, not the
layout. All faces of one camera share a size, so a batch stays one tensor and
the fused projection path still applies.

Each post camera carries its own intrinsics (an off-centre principal point is
what a crop is) and its frame's rotation, and the GT warp samples through
exactly that table, so the face that is rendered is the face that was warped.

### One face size, or one per lens

`--warp-face-fit` decides whether the faces of one camera share a size.

- **`uniform`** is the plan above: one tile covers every crop, so a batch
  renders all of a camera's faces in one pass -- which is what the fused
  projection-backward optimizer (`--use-fused-proj-bwd-optim`) needs.
- **`per-face`** gives each frame's crop its own face, rounded up to a multiple
  of 32 px so near-equal faces still share a size, and renders **one pass per
  distinct size**, accumulating gradient across the passes. Every pass is
  weighted by its face count, so a face carries exactly the weight it would in
  one batch -- the loss, the gradient and the densification score are the
  uniform plan's, only the pixels differ.
- **`auto`** (default) picks `per-face` only when the run already renders in
  several passes -- when the fused optimizer is off anyway, either because it
  was turned off or because `--split-batch` is active on a step of more than
  one image -- and the crops save at least 10% of the pixels. Otherwise, or
  when the saving would have to pay for the fused optimizer itself, it stays
  `uniform`.

Splitting is per input image, so a face pass renders that image's faces of one
size; the passes of a step accumulate into the same gradient buffers and one
optimizer step follows.

Measured, 3000 steps on an RTX 5070 laptop, gradient already accumulating
across passes in both columns (`--use-fused-proj-bwd-optim 0 --split-batch 1`):

| capture | fit | faces | passes | time | peak VRAM |
|---|---|---|---|---|---|
| 1000x1500 fisheye, 108x162 deg | uniform | 6 x 548x298 | 1 | 40.4 s | 804 MiB |
| | per-face | 5, 548 wide | 3 | 41.1 s | 684 MiB |
| 960x960 fisheye, ~200 deg | uniform | 5 x 430x430 | 1 | 60.1 s | 978 MiB |
| | per-face | 5, 430 wide | 2 | 60.0 s | 942 MiB |

So where the passes are already paid for, per-face costs nothing in time and
returns 4-15% of the VRAM. Where they are not, it costs the fused optimizer:
the same capture with it available runs 37.2 s / 722 MiB uniform against
40.9 s / 690 MiB per-face -- 10% slower to save 4% of the memory, which is why
`auto` leaves it alone there. The remaining waste in both fits is the masked
part of a face, which the rasterizer still runs; skipping fully-masked tiles
would take it without splitting anything.
Rasterization still runs over a face's masked-out pixels; a lens boundary that
is not a rectangle leaves some, which the synthesized FOV mask excludes from
the loss. A lens whose visible rays fit one face is not split at all.

`spirula geometry --check` verifies the planner: every visible ray of each
test camera must land in a face, and the faces must not exceed the uncropped
pixel count.

That normalized-frame similarity is where `orientation_method` and
`center_method` act — and only `up` / `poses` is implemented natively;
anything else is approximated with a warning. Since `train_frame="points"`
leaves splats in the raw frame, the choice moves `train_frame_scale` and the
viewer's default camera, not the training coordinates. The unported methods
(`pca`, `vertical`, `gsplat`, `focus`) still have a working Python reference:
[notes/pose-normalization.md](notes/pose-normalization.md).

## Train/eval split

`eval_mode` selects the strategy:

- `all` (default) — every image trains.
- `fraction` — linspace-spread `ceil(N * train_split_fraction)` train images.
- `interval` — index `% eval_interval == 0` is eval, rest train.
- `filename` — basename contains `train` / `eval`.

`validation_fraction` additionally holds out a linspace-spread slice for
validation. Frames whose camera position exceeds `outlier_threshold` MADs from
the geometric median of all camera positions are rejected (default: off).

`require_image_files=false` keeps frames whose image file is missing — the
standalone viewer uses this to load camera poses from a dataset shipped
without pixels.

## Downscaling

Stored intrinsics can be divided by a fixed factor (the Mip-NeRF 360
`images_2` / `images_4` convention). Note: the auto-detect mode that probes
the first image's actual resolution is implemented on the Python side only;
the native parser currently requires an explicit factor.

## Preprocessing tools

`reference/scripts/` holds standalone Python utilities that produce these
layouts: frame extraction with blur skipping, COLMAP/GLOMAP driving,
Metashape conversion, downscaling, undistortion, masking (including a SAM2
GUI), monocular depth/normal prediction, and raw conversion. See
`reference/scripts/README.md`. These are *preprocessing*, separate from the
training data path, and stay on the Python side.

`reference/scripts/batch_process_data.bash` needs a COLMAP vocabulary tree; set
`SS_VOCAB_TREE` to its path.

## Benchmarking

`reference/python/benchmark.py` drives multi-scene runs over standard
academic sets (Mip-NeRF 360 `360_v2`, ZipNeRF) by calling `spirula train`.
Dataset roots are passed as arguments — no paths are hardcoded. **Each scene runs in its own subprocess**: the engine is
a process-global singleton, and running several scenes in one process leaks
state between them and silently degrades metrics.

## Do not hardcode local paths

Dataset locations belong in arguments or environment variables. See the
"Do not commit" section of [`../AGENTS.md`](../AGENTS.md) and
`tools/check_private_paths.sh`.
