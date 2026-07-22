# Datasets

## Supported layouts

Three formats, auto-detected in this order: **Nerfstudio**, then **COLMAP**,
then **Metashape**. Both the native parsers (`src/data/parsers/*Parser.cpp`) and the
Python parsers (`spirulae_splat/modules/dataparser.py`) probe in the same
order.

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

## Two-stage parse

1. **`parse_dataset`** → `ParsedDataset`: per-**input** cameras in the raw
   `train_frame="points"` frame — poses and points exactly as stored. The
   normalized-frame similarity is computed only to obtain the
   `train_frame_scale` scalar.
2. **`bake_post_split`** → the post-split arrays `engine_setup_data_manager`
   consumes. This is either an identity (K=1) pass-through, or the cubemap
   split when `warp_to_pinhole` is enabled: **5 faces** for fisheye/equisolid,
   **6 faces** for equirectangular.

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

`scripts/` holds standalone Python utilities that produce these layouts:
frame extraction with blur skipping, COLMAP/GLOMAP driving, Metashape
conversion, downscaling, undistortion, masking (including a SAM2 GUI),
monocular depth/normal prediction, and raw conversion. See
`scripts/README.md`. These are *preprocessing*, separate from the training
data path, and stay on the Python side.

`scripts/batch_process_data.bash` needs a COLMAP vocabulary tree; set
`SSPLAT_VOCAB_TREE` to its path.

## Benchmarking

`spirulae-benchmark` drives multi-scene runs over standard academic sets
(Mip-NeRF 360 `360_v2`, ZipNeRF). Dataset roots are passed as arguments — no
paths are hardcoded. **Each scene runs in its own subprocess**: the engine is
a process-global singleton, and running several scenes in one process leaks
state between them and silently degrades metrics.

## Do not hardcode local paths

Dataset locations belong in arguments or environment variables. See the
"Do not commit" section of [`../AGENTS.md`](../AGENTS.md) and
`scripts/check_private_paths.sh`.
