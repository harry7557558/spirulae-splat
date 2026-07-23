# Standalone CLI Trainer (`ssplat-train`) + Native GUI (`ssplat-gui`)

Working notes for the no-Python trainer under `src/app/`. Written 2026-07-10,
GUI added 2026-07-12; covers code structure, verification status, and
remaining ports. Long-term context: Phase 1 (standalone C++ trainer) and
Phase 2 (Dear ImGui/GLFW GUI) of the GUI-app plan are done; Phase 3 is
Windows/Linux packaging; kernels stay in Slang for the eventual
Vulkan/cross-vendor port.

## Build & run

```bash
cmake -G Ninja -B build -DSSPLAT_BUILD_CLI=ON && cmake --build build --target ssplat-train
./build/ssplat-train [<preset>] --data <colmap_dataset_dir> [--flag value ...]
```

- Presets = tyro subcommands: `3dgs` (default), `360-camera`, `in-the-wild`,
  `linear-color`, `synthetic`, `meshing`, `academic-baseline`.
- Flag conventions: flattened names (`--sh-degree`, not `--model.sh-degree`);
  `-`/`_` interchangeable; bools take a value (`--warp-to-pinhole 1`);
  `--key=value` works (arity-1 only); `none` clears optionals; tuples take N
  values (`--bilagrid-shape 8 8 4`). `--help` shows preset-resolved defaults.

### Building without Torch/Python

Torch + Python are only needed for the Python extension (`ext.cpp`); the
engine and app code are torch-free. CMake falls back automatically when
`import torch` fails or python3 is missing; `-DSSPLAT_NO_TORCH=ON` forces it:

```bash
cmake -G Ninja -B build_notorch -DSSPLAT_NO_TORCH=ON && cmake --build build_notorch
```

In this mode `csrc` is a STATIC lib (the engine API has no dllexport
annotations, and this yields a self-contained exe — deps are just libcudart +
system libs), `SSPLAT_BUILD_CLI` is forced ON, CUDA archs come from
`nvidia-smi --query-gpu=compute_cap` (override with `-DTORCH_CUDA_ARCH_LIST`),
and the libpython link + static-libstdc++/nftw interposition workarounds are
skipped (they exist only because of libtorch). Generated headers
(`src/generated/`, `app/generated/`, `src/instantiations/`) are committed, so a fresh
checkout builds with no Python at all. With Torch present the extension build
is unchanged (shared `libcsrc`, same flags as before).

Windows (VS2022 + CUDA toolkit): run `build_develop.bat` from any cmd prompt —
it locates VS via vswhere and ALWAYS calls vcvars64 (an ambient cl/INCLUDE may
reference an uninstalled SDK), picks the newest installed CUDA toolkit
(ambient CUDA_PATH may pin one too old for the MSVC in use), falls back to
VS-bundled cmake/ninja, and configures with `-DSSPLAT_NO_TORCH=ON` (a broken
torch install would abort configure from inside TorchConfig.cmake — QUIET
can't suppress errors raised inside a found package config). Manual
equivalent from a vcvars64 shell:
`cmake -G Ninja -B build -DSSPLAT_NO_TORCH=ON -DCMAKE_BUILD_TYPE=Release`.

## Native GUI (`ssplat-gui`, Phase 2)

Dear ImGui + GLFW + OpenGL 3.2-core desktop app ("Spirulae Splat").
Off by default; enable with:

```bash
cmake -G Ninja -B build -DSSPLAT_BUILD_GUI=ON   # composes with -DSSPLAT_NO_TORCH=ON
cmake --build build --target ssplat-gui
```

GLFW 3.4 and imgui v1.92.8 are pinned and fetched at configure time via
FetchContent (network needed once; on WSL drvfs mounts git may need
`safe.directory` entries for `build/_deps/{glfw,imgui}-src`). Only system dep
is OpenGL + X11/Wayland dev packages on Linux; Windows (MSVC) and macOS
(GL 3.2 forward-compat) code paths are in place — macOS still needs the
CUDA→Vulkan port to actually run the engine.

Design: novice path is Home → "Open a Dataset" (or "Create Dataset from
Photos/Video", which drives external `colmap`/`ffmpeg` CLIs with live log +
cancel) → preset dropdown + a curated Basic Options list → Start Training →
live native viewport. Advanced path: "All Options" editor **generated from
the `SSPLAT_CONFIG_FIELDS` X-macro** — all config fields, grouped by
sub-config, searchable, Python docstrings as tooltips, modified-from-preset
highlighting with right-click reset; new Python config fields appear
automatically after codegen. Web viewer can be additionally served from the
GUI (Basic Options) for remote monitoring.

| File | Role |
|---|---|
| `gui/GuiMain.cpp` | GLFW window + GL 3.2 core context + ImGui bootstrap, dark style, DPI scale, frame loop, close-confirm flow. **Drag-and-drop** (glfwSetDropCallback → `GuiApp::handle_drop`): a dropped path is auto-detected as an SfM dataset folder (transforms.json / sparse/ / colmap/ marker, or a Metashape camera .xml + point-cloud .ply pair → open), a photo folder (contains images → COLMAP screen), or a video file (extension → COLMAP video screen, .insv preset applies); files from inside a dataset (transforms.json, .db/.bin/.txt/.xml) open their parent. Video/photo drops are ignored while training (datasets go through the stop-confirm flow). |
| `gui/GuiApp.h/.cpp` | Screens (Home / COLMAP / Train), layout, wiring; recents + tool paths persisted to `~/.config/spirulae-splat/gui.conf` (`%APPDATA%` on Windows). Session-destroying navigation (Home / open-dataset / quit during training) goes through a stop-and-save confirm modal with a deferred pending-action; output folder defaults to `<dataset>/outputs` with a Browse button + resolved-run-path preview. Editing any dataset-parsing option (dataparser group via the generated `parse_settings_equal`, plus warp/load/scale flags) marks the dataset dirty and auto-reloads it once the edited widget loses focus. NOTE: `open_dataset`/`add_recent` take `std::string` **by value** — callers pass `_recents` elements and `add_recent` mutates that vector (a const& dangles; this was a real bug that corrupted the path to ""). |
| `gui/ConfigUI.h/.cpp` | The X-macro-generated options editor. Zero per-field special cases (only `data` is hidden, managed by the dataset picker). |
| `gui/ViewportPanel.h/.cpp` | Native viewport with two backends behind the browser-identical NavCamera navigation: **Preview** = GL point cloud + frusta as soon as the dataset parses (PreviewRenderer), **Engine** = RenderWorker once training starts, with all four web-viewer camera models (Pinhole / Fisheye-equidistant / Fisheye-equisolid / Equirectangular; `fovToIntrinsics` + per-model FOV ranges ported from viewer.html). Buffer picker, camera-frusta overlay with a live **frustum-size slider** (`ViewRequest::cam_size_scale`), render-scale + live-refresh throttle (0.15 s). Initial framing: seed-point centroid target, median camera distance, camera-centroid direction. **Double-click centering** (viewer.html `recenterAt`): pan laterally so the 3D point under the cursor sits on the optical axis and make it the orbit pivot (rotate toward it when behind the camera plane, >180° models); preview mode picks the nearest displayed point along the cursor ray (3% angular cone, `PreviewRenderer::pick_point`), engine mode attaches `pick_px/py` to the next render and gets the point back in the ViewResult — depth-channel readback, no extra VRAM or render pass; all four display camera models via `viewer_pixel_ray`. |
| `gui/NavCamera.h/.cpp` | 1:1 port of viewer.html's `cam`/`quat`/`Nav`: quaternion camera, four modes (Turntable / Trackball / First Person / Free Fly), same sensitivities and mappings for **mouse** (LMB orbit-or-look, RMB/MMB/Shift pan, wheel dolly), **keyboard** (WASD/arrows, E/Q up-down or Fly-roll, active while the pointer is over the viewport), **gamepad** (GLFW gamepad API: left stick move, right stick look, triggers up-down/roll -- including the browser quirk that triggers only translate while the left stick is deflected), and **touch** via the OS's pointer/gesture emulation (single finger = orbit; system pinch/pan gestures arrive as wheel; GLFW exposes no raw multitouch). Keep in sync with viewer.html's Nav. The initial pose replicates `cam.reset()` verbatim (target = client-frame origin = the CAMERA-POSE center via center_method="poses", pos=[0,0,1], orbit(0,-250)) -- verified pixel-equivalent against a `/render` fetch from ssplat-train's web viewer at the client's default c2w. `SSPLAT_NAV_DEBUG=1` logs per-frame mouse-drag nav decisions (button/pan/target/pos) to stderr. |
| `gui/PreviewRenderer.h/.cpp` + `gui/GlLoader.h/.cpp` | Offscreen-FBO GL renderer for the dataset preview (sparse points, vertex-colored + per-input camera frusta as center+offset lines scaled by a uniform). Frustum wireframes are built from each camera's TRUE model + distortion via `frustum_template` — a C++ port of `viewer/js/dataset.js` `frustumTemplate`/`generateRay`/`undistort` (itself a port of `fill_frustum_segments_kernel` / `projection_utils.cuh`): pinhole = classic border pyramid, wide models (fisheye/equisolid/equirect) = dome/globe with dimmed interior gridlines, Newton undistortion with the `is_valid_distortion` bail-out + bisect-toward-center for out-of-domain border pixels. (Before this, frusta were tan-based pinhole pyramids — a ~200° THIN_PRISM_FISHEYE dataset drew as giant crisscrossing triangles.) The display-side vertex shader implements the same camera models as the engine viewer (pinhole / fisheye-equidistant / equisolid / equirectangular; `u_s` = fx/(W/2), fy/(H/2), linear view-distance depth), so the preview -> engine transition at training start keeps pose AND intrinsics with no jump (ViewportPanel `maybe_frame` re-frames only when the dataset identity changes; navigation done in the preview carries into training, verified pixel-aligned at fisheye 333°). The projection function is shared with the fragment shader, which **re-projects the interpolated view-space position and discards fragments whose reprojection error exceeds 5% of the viewport** — the web viewer LINE_FS seam fix, killing the full-width streaks that segments crossing the equirect ±180° seam / fisheye backward point / pinhole near plane used to rasterize as. That check alone still keeps ~5%-viewport stubs at both ends of a wrapped segment (the error vanishes at the endpoints), which stack into ladder artifacts at the equirect edges — so line vertices also carry `a_delta` (this-endpoint − other-endpoint, pre-scale) and the vertex shader kills the WHOLE segment when its chord crosses the equirect seam half-plane (x=0, z>0 view space); same mechanism as the web viewer's `aDelta`/`vKill`. Optional **grid + axes** (viewport "grid" checkbox): power-of-10 ground grid + colored positive half-axes generated in the **training (saved-splat) frame** — lines mark round coordinates of the exported model — then mapped into the preview's normalized frame; the cell decade adapts to the nav distance and the finite line patch recenters on the orbit target (lattice-snapped, with hysteresis), cell edges subdivided 4x so they curve correctly under the non-pinhole projections; depth-tested like all preview geometry (`ensure_grid`). GlLoader = minimal namespaced (`glx::`) proc loader via glfwGetProcAddress; no GLEW/GLAD. |
| `gui/TrainRunner.h/.cpp` | TrainerSession on a worker thread: phases Idle→Loading→Ready→Preparing→Training→Done, pause/stop, metric history for plots, optional ViewerServer. Lifetime rule: viewport detaches before the session is replaced. |
| `gui/ColmapRunner.h/.cpp` | images/video → COLMAP dataset. Requires **COLMAP >= 4.x** (version probed from `colmap help`; 3.x-era `use_gpu` flags are gone — flags follow `scripts/run_colmap.bash`): feature_extractor (**SIFT or ALIKED** via `FeatureExtraction.type`; camera model = any the parser supports; single camera / per-folder / per-image; optional **initial focal length** — `init_focal_factor` × width probed via `stbi_info`, composed into `ImageReader.camera_params` with centered principal point + zero distortion, or a raw `camera_params` string) → **explicitly chosen** exhaustive / sequential / vocab-tree matcher (no "auto"; the GUI presets sequential for video, exhaustive for photos). Sequential supports overlap + **quadratic overlap** (on by default) + **loop closure** (`SequentialMatching.loop_detection` via the vocab tree — SIFT only; the tree is auto-found near the workspace/cache or downloaded via curl into `~/.cache/spirulae-splat`). **LightGlue** matching (`FeatureMatching.type` `SIFT_LIGHTGLUE`/`ALIKED_LIGHTGLUE`, default for ALIKED) → mapper (`ba_use_gpu` — **forced OFF for fisheye models**, COLMAP's GPU BA doesn't support them; `Mapper.ba_refine_extra_params 0` for perspective models by default — distortion held fixed for stability and recovered in the final BA, per run_colmap.bash's advice; `min_num_matches`) → best-effort **model_merger** when the mapper splits (kept only if the merged model registers more images; written as the next `sparse/<N>` — on a real X5 capture this fused 86+39 partials into 116/118 frames) → optional bundle_adjuster refinement **on the largest/merged model** (was hardcoded `sparse/0`), with a **verify-and-revert guard**: mean reprojection error before/after via `model_analyzer`, refinement discarded when worse/non-finite (releasing pp + 8 thin-prism coefficients on a ~200° fisheye reliably diverges to ~1e150 px — pp additionally stays fixed for fisheye). `.insv` preset: THIN_PRISM_FISHEYE, per-folder cameras, focal factor 0.269 (Insta360 X5: fx=fy≈0.269·width), **exhaustive** matcher (the two lens tracks are concatenated, so sequential misses cross-lens pairs: 116/118 exhaustive vs 68/118 sequential+loop on the same capture). **Resume** (`ColmapJob::resume`, GUI checkbox appears when the workspace holds a previous run): extracted frames are kept per track, mask.py resumes on its own, COLMAP skips features/matches already in database.db, and existing `sparse/<N>` models skip the mapper (the mapper only writes models on completion, so partials are never trusted); with resume off a non-empty workspace is refused rather than mixed into. A measured resume replay of a full .insv run took 39 s vs 246 s. **Repetitive-scene knobs** (Advanced → "Repetitive scenes", with an Off/Low/Medium/High preset combo; all 0 = default): `SiftMatching.max_ratio`, `TwoViewGeometry.min_num_inliers`, `Mapper.abs_pose_min_num_inliers` / `abs_pose_min_inlier_ratio` / `abs_pose_max_error` — stricter matching + registration so similar-looking rooms don't weld together. The GUI's default workspace auto-suffixes `_2`, `_3`, ... instead of pointing at an existing non-empty folder. Optional **AI masking**: the CMake-embedded `scripts/mask.py` runs via external Python (lang-segment-anything / SAM-3), masks feed `ImageReader.mask_path` + the trainer's masks dir; missing-package output is detected and surfaced with install instructions. Photo-folder inputs are indexed **in place**, recursively (no copy) — the absolute image dir is handed to the GUI in-memory for the immediate open; no marker file is written, so on later re-opens set `data.image_dir` in the dataparser options (video datasets use the default `images/`). Multi-track `.insv` videos split into `images/cam<N>/` + per-folder cameras. The mapper may emit several partial models — the trainers auto-pick the largest (below). |
| `gui/FrameSelect.h/.cpp` | Blur-aware video frame selection (multithreaded C++ port of `scripts/extract_frames.py`): ffmpeg extracts candidates at (target fps x sharpness window) with `-nostdin` (**without it ffmpeg can hang reading stdin** — the "stuck at ffmpeg" bug), then the sharpest per window is kept (variance of 3x3 Laplacian on mean-subtracted 512^2 gray, scored across all cores). |
| `gui/FileDialog.h/.cpp` | Built-in ImGui file/folder browser (no native-dialog dep; works over WSLg/remote X). |
| `gui/Subprocess.h/.cpp` | Cross-platform subprocess with merged stdout/stderr line streaming + kill-on-cancel (fork/execvp + process group on POSIX, CreateProcess + pipe on Windows). Bare `\r` ends a log line (ffmpeg/curl progress). |

The frustum-size control also reaches the web viewer: viewer.html gained a
"Camera Size" slider sending `camera_size_scale` per /render; the render
worker applies it via the new `engine_viewer_set_camera_size()` (Engine.h /
Visualizer.cu), which updates the scale without re-running
`engine_viewer_init` (that would wipe the thumbnail cache).

Verified 2026-07-12 (WSLg, RTX 5070, COLMAP 4.2, driven via synthetic X
input): garden → dataset preview (points + frusta + size slider) → train
400-600 steps (8-24 ms/step incl. live viewport) → orbit/zoom/depth-buffer/
frusta-with-thumbnails → checkpoint saved → Train Again; Home-during-training
confirm modal; video (60-frame test clip) → ffmpeg candidates → sharpness
selection → COLMAP sequential → mapper → BA refine → "Open in Trainer" →
trained; `mask.py` invocation verified against a real lang-sam install
(24 masks). CLI re-verified unchanged after the TrainerCore/RenderWorker
extraction (see below).

## Files

| File | Role |
|---|---|
| `main.cpp` | CLI parsing (`--help`, flag table), stdout progress printing, web-viewer wiring. The engine plumbing lives in TrainerCore (below). |
| `TrainerCore.h/.cpp` | **Shared trainer session** (extracted from main.cpp 2026-07-12, code moved verbatim): color-space resolution, `seed_splats`, `build_step_config`/`build_loss_weights`, `scheduled_lr`, `save_config_json`, and `TrainerSession` = check_config / load_dataset / setup_engine / train(callbacks) / save_checkpoint / progress_json + the engine mutex + pause/stop/render-pending atomics + `make_viewer_config/hooks`. Used by both `ssplat-train` and `ssplat-gui`. |
| `RenderWorker.h/.cpp` | **Shared interactive render path** (extracted from Viewer.cpp 2026-07-12): latest-wins worker thread, c2w remap + engine render + display transforms + `engine_blit_view`, returns RGB8; `viewer_upload_cameras()` + `viewer_upload_grid()` (axes/grid overlay axis-aligned in the **engine/saved-splat frame**, drawn when `ViewRequest::show_grid`; `grid_dist` + `grid_target` — nav distance and orbit target in the client's normalized frame, remapped like the c2w — drive the engine's zoom-adaptive cell decade and lattice-snapped patch recentering per render). Buffer keys: distortion buffers are offered/computed **only when a distortion regularizer is configured** — they are full-resolution never-freed pool allocations, so zero reg weights = zero extra VRAM (previously `distortion_reg_on` forced the distortion channels on EVERY viewer render). `sh` and `refinement_score` debug renders are ported from model.py (`engine_debug_forward` with a **max_num_splats-sized** DC override — cur-sized throws a tensor-size mismatch; `engine_copy_accum_buffer` col0 for the score). **Pick** (`ViewRequest::pick_px/py` -> `ViewResult::pick_hit/pick_point`): 3D point under a pixel from the already-downloaded ray-depth channel (alpha-unpremultiplied like the blit kernel, alpha > 0.1 to reject background), CV pixel ray via `viewer_pixel_ray` (host generate_ray port, 4 display models) through the remapped c2w, result mapped back to the client's normalized frame. Viewer.cpp adds HTTP+JPEG on top; the GUI viewport uploads to a GL texture. |
| `DatasetParser.h` | Public dataset structs: `DatasetParserConfig`, `ParsedDataset` (per-INPUT cameras), `PostSplitCameras` + `bake_post_split()` (warp expansion), parse fns, shared `dsparse::` helpers. |
| `ColmapParser.cpp` | COLMAP binary + text reader + format auto-detect dispatcher (`parse_dataset`: transforms.json → nerfstudio, else try COLMAP, else Metashape — Python's probe order). When `colmap_recon_dir` is not set, models under `sparse/*` and `colmap/sparse/*` are enumerated and the one with the most registered images wins (sparse/0 is NOT necessarily the largest; count read cheaply from the images.bin header / images.txt comment). Python's `dataparser.py _parse_colmap_data` does the same. |
| `NerfstudioParser.cpp` | transforms.json reader + self-contained PLY point reader (ascii + binary_little_endian). Split into `parse_nerfstudio_dataset` (reads the file) and `parse_nerfstudio_meta` (consumes a transforms-shaped `JsonValue`; the Metashape front-end feeds this). |
| `MetashapeParser.cpp` | Metashape camera-export `.xml` + `.ply` reader (port of `_parser_metashape_data` + `metashape_utils.py`): sensor intrinsics (`calibration[@class!='initial']`, p1/p2 swap, b1/b2÷f), component transforms, OpenCV→OpenGL flip, camera→image matching (photo-path suffix via the optional `.psx` project's zipped camera table, else label substring). Builds a transforms-shaped meta → `parse_nerfstudio_meta`. |
| `DatasetCommon.cpp` | Shared bakes: normalized-frame scale, eval/val splits, aux-file discovery, geometric-median outlier filter, **`bake_post_split`** (fisheye 5-face / equirect 6-face cubemap expansion; axes MUST match `DataManager.cpp` `kAxesFisheye5`/`kAxesEquirect6`). |
| `Json.h` | Minimal dependency-free JSON parser (handles Python `Infinity`/`NaN`). |
| `Xml.h` | Minimal dependency-free XML parser (ElementTree subset: `attr`/`find`/`findall`/recursive `iter`; comments/PI/CDATA/DOCTYPE skipped, standard entities). Metashape nests `<camera>` in `<group>` and `<camera_ids>` in `<partition>` — the recursive `iter` is load-bearing. |
| `Knn.h` | Exact kd-tree kNN (multi-threaded queries) for `seed_splats` scale init; replaced the hash-grid approx that degenerated to O(N²) on SfM outliers. |
| `HttpServer.h/.cpp` | Minimal HTTP/1.0 GET server (POSIX sockets; winsock shim compiles but untested). Serial request handling — parity with Python's non-threading `HTTPServer`. |
| `Viewer.h/.cpp` | Web-viewer server port (viewer/server.py + http_server.py + render_worker.py + annotation.py): latest-wins render worker, `get_outputs` viewer subset, `engine_blit_view` GPU annotation/colormap, stb JPEG encode. Serves the **unchanged** `viewer.html` (embedded at configure time via CMake hex; `SSPLAT_VIEWER_HTML=<path>` env overrides for dev). `/pick?px=&py=&<camera params>` returns the 3D point under a pixel as JSON for viewer.html's double-click centering (the Python server has no /pick; the client treats non-OK responses as a no-op). |
| `generated/cli_config.h` | AUTO-GENERATED — do not edit. `SsplatConfig` struct (all 189 config fields, defaults baked), `SSPLAT_CONFIG_FIELDS(X)` X-macro flag table, `ssplat_apply_preset()`. |
| `../../../../generate_cli_config.py` | The generator. AST-parses the Python config dataclasses (no torch import → works on fresh checkout). Run by `build_develop.bash`. |
| `../external/` | All vendored third-party code (marked `linguist-vendored` in `.gitattributes` along with the generated dirs): `stb_image.h`/`stb_image_write.h` (images), `npy.hpp` (checkpoints), `miniz.c/.h` (zip reading for the Metashape `.psx` camera table; compiled into `ssplat-train` only). |

Debug: `SSPLAT_DUMP_CAMERAS=<path> ssplat-train ...` dumps parsed + post-split
camera arrays as JSON and exits before engine setup — used to diff against
the Python dataparser/trainer algebra (see verification notes below).

## Mesh extraction (`ssplat-mesh`)

C++ twin of `spirulae_splat/ss_meshing.py` (same engine code path:
`Meshing.h` / `MeshingHost.cpp` / `MeshUV.cpp` / `MeshExport.cpp`). Built by
the same `SSPLAT_BUILD_CLI` block; `mesh_main.cpp` + the dataset parsers, no
HTTP/viewer.

```bash
./build/ssplat-mesh <ckpt> [--data <dir>] [--format ply,obj,gltf,glb] \
    [--color none|vertex|texture] [--texture-size 0] [--flag value ...]
```

- `<ckpt>` = run dir (config.json + `step-*.ckpt/`), a `*.ckpt` dir, or a
  `splat.ply` directly. `--data` defaults to config.json's `data`; the
  dataparser settings (image_dir, recon dir, metashape paths, numeric
  rescale) and `model.relative_scale` are honored from config.json.
- Color modes: `none`, `vertex` (per-vertex RGB), `texture` (LSCM UV atlas +
  baked image). Unsupported pairs error out up front: PLY+texture, OBJ+vertex.
- With `--color texture`, a format token may carry the texture encoding:
  `glb+png` (default), `glb+jpg` (JPEG q95), `glb+jpeg75` (JPEG q75) — works
  for obj/gltf/glb; JPEG is part of core glTF (`image/jpeg`).
- Output: one file per `--format` next to the checkpoint's splat.ply
  (`mesh.ply` / `.obj`+`.mtl`+tex / `.gltf`+`.bin`(+tex) / `.glb`), all
  emitted dependency-free (PNG/JPEG via vendored stb_image_write). `.glb`
  embeds the texture; validated clean with gltf_validator.
- glTF compatibility: plain PBR material (no KHR_materials_unlit — partially
  supporting viewers rendered it solid white), uint16 indices when they fit,
  float VEC3 COLOR_0 — matching the Khronos sample-model conventions.
- `splat.ply` reader expects float32 binary-little-endian properties (what
  both the Python trainer and `EngineCheckpoint.cpp` write).

## Config codegen (source of truth = Python dataclasses)

`tools/codegen/generate_cli_config.py` parses `TrainerConfig` (+ preset
subclasses) in `modules/trainer.py` and the nested
`SpirulaeSplatDataParserConfig` / `SpirulaeSplatDataManagerConfig` /
`SpirulaeSplatModelConfig` / `OptimizerConfig`. Field docstrings become help
text; preset `default_factory` lambdas become `ssplat_apply_preset` branches.

- **Collisions**: flattening drops group prefixes; the generator hard-errors
  on un-listed collisions. Resolve in its `RENAMES` dict. Currently only
  `datamanager.split_batch` → `--dm-split-batch` (the model one keeps the
  plain name; datamanager's is the legacy Python-path OOM workaround, no-op
  on the managed path).
- **Type mapping**: `Optional[int/float]` → `std::optional`;
  `Optional[str/Path]` → `std::string` with `""` = None; `Literal[str...]` →
  string + validated choices; `Literal[True,False,None]` → `optional<bool>`;
  `Tuple[...]` → `std::array`; `Union[bool,int]`
  (`rescale_camera_to_fit`) via `TYPE_OVERRIDES` → float (0=off, -1=auto,
  >0=factor).
- **Verified** (2026-07-10): all 189 defaults + all 7 presets match
  instantiated Python configs; a full `in-the-wild` command line with mixed
  flag styles resolves identically to `tyro.cli` (189/189 fields).
  Verification scripts were session-scratch; the approach: flatten a tyro
  config, diff against the CLI's `config.json` dump.

## Python → C++ port mapping (all in `TrainerCore.cpp` unless noted)

As of phase 7 (docs/restructure-proposal.md §4.3) this table documents
*history*, not a live duplication: `TrainerCore.cpp` is compiled into the
engine library and bound to Python as `_C.TrainerSession`
(`bindings/bind_trainer.cpp`), so the C++ column is the single
implementation and the Python column is a client. The Python functions
listed here still exist — users are on that path and breakage is paced, not
rushed — and `tests/python/test_trainer_parity.py` asserts the two produce an
identical `EngineStepConfig` for every step of every preset. Change one side
without the other and that gate fails.

| Python | C++ |
|---|---|
| `optimizer.py:54 get_scheduled_lr` | `scheduled_lr()` — pass `nullopt` final when the `*_lr_final` attr doesn't exist in Python (getattr default = constant LR) |
| `model.py:1483 _build_loss_weights` | `build_loss_weights()` — 19 `LossWeightIndex` slots, YUV group normalized to `1-ssim_lambda` |
| `core.py:280 _build_optim_config` | in `build_step_config()` — incl. scale-agnostic mean, `scale_reg/alpha`, quantization level→bits (core.py:95), eps_tr schedule |
| `core.py:321 _build_densify_config` | in `build_step_config()` — `noise_lr_scalar = use_revised ? 1 : alpha`, `max_world_size * alpha` |
| `model.py:1887 engine_train_step_managed` per-step args | in `build_step_config()` — bilagrid/PPISP LRs (AdaGrad switches to constant `*_adagrad_lr`), background noise ramp / SH LRs, color-shift EMA beta |
| `model.py:546 populate_modules` (seeding) | `seed_splats()` — min_init/cap resolution, 4-NN log scales (exact kd-tree, `Knn.h`), logit opacity, DC from SfM color incl. gamut conversion |
| `model.py:823 _maybe_init_bilagrid` | static init before the loop (dataset modalities known up front); bilagrid_shape (X,Y,W) → engine (L=W,H=Y,W=X) |
| `model.py:505-543` background + color space | `resolve_color()` + init calls; gamut matrices ported from `_wrapper_per_pixel.py:111` |
| `trainer.py:207 _setup_cpp_data_manager` (incl. warp) | DataManager setup; batch-size policy from `datamanager.py:91` (non-stochastic); K/post arrays from `bake_post_split` |
| `trainer.py:257-414` warp expansion | `DatasetCommon.cpp bake_post_split()` — K (fisheye/equisolid=5, equirect=6 or direct), cubemap-axis c2w expansion, canonical PINHOLE intrins at `out_shape=ceil(sqrt(HW/K))`, viewmat bake over POST arrays, per-INPUT intrins for the wide-warp kernel, synthetic-FOV-mask `load_masks` enable (trainer.py:452-461), direct-equirect canonical intrins + depth/normal guard |
| `dataparser.py` COLMAP branch | `ColmapParser.cpp` — w2c→c2w flip (dataparser.py:644), `train_frame_scale` over ALL frames pre-split, eval_mode splits, aux mask/depth/normal discovery, outlier filter |
| `dataparser.py _parse_nerfstudio_data` | `NerfstudioParser.cpp` — frame/meta intrinsics fallback, DISTORTION_KEYS, .webp/missing-file skip, filename sort, geometric-median outlier filter (dataparser.py:138-173), equirect canonical intrins, `applied_transform` inverse on poses+points (train_frame="points", dataparser.py:501-531), frame-specified mask/depth/normal paths with directory-probe fallback |
| `trainer.py:741 train` / `:939 save_checkpoint` | loop + `save_checkpoint` lambda (step-%09d.ckpt dirs, latest-only pruning) |
| `viewer/*` + `trainer.py render/_render/get_progress/toggle_pause` | `Viewer.cpp` + main.cpp hooks — engine mutex shared with the train loop, `render_pending` fairness yield (trainer.py:146-153), pause gate, `/progress` JSON, c2w normalized→train remap via `ParsedDataset.train_to_normalized` (verified vs `dpo['train_to_normalized_transform']`), `camera_size` kNN heuristic (annotation.py:32) |

## Verified working

- All presets train and converge (banana COLMAP set — 3dgs: PSNR 26.7 @ 500
  steps); linear-color exercises color-space init + trust region; meshing
  exercises 3dgut + median losses; academic-baseline exercises interval eval
  split + non-FPBO + fp32.
- **360-camera preset works end-to-end** (SharkWipf_SampleDataset, 400×3840²
  nerfstudio fisheye + masks): 5-face warp (400→2000 post cameras), warped
  GT upload, synthetic FOV masks, bilagrid/PPISP at n_post slots. Behavior
  parity vs Python trainer at 300 steps: PSNR 17.7 vs 17.1, SSIM 0.775 vs
  0.737, same ~72 ms/step (deltas are RNG-level: seeding/batch order).
- Numeric verification (2026-07-10): C++ `SSPLAT_DUMP_CAMERAS` dump vs
  Python — nerfstudio c2w/points/order/train_frame_scale AND the full warp
  expansion (K, offsets, 2000×viewmats/intrins/dist, input intrins) all
  match to float32 precision.
- **Metashape parser** (2026-07-10, dye_alley 732-frame dual-fisheye rig):
  all 732 cameras match the Python metashape path to float32 precision
  (c2w, intrins, dist incl. p1/p2 swap + b1/b2 scaling, points, order,
  train_frame_scale); `.psx` camera-table disambiguation exercised (labels
  alone are ambiguous across `cam1`/`cam2`); trains end-to-end and
  auto-detect falls through COLMAP→Metashape correctly. Deliberate
  deviation: with multiple components we train on the **largest** (Python's
  reversed-sort quirk picks the smallest).
- Checkpoint cadence + `save_only_latest_checkpoint` pruning + final save.
- **Web viewer** (2026-07-10, SharkWipf run): browser client unchanged; `/`,
  `/render` (rgb / depth / alpha / depth_normal / distortion buffers, JPEG,
  training-camera frusta + thumbnails via `engine_blit_view`), `/buffers`,
  `/progress`, `/pause-toggle` all exercised live during training —
  step latency stayed ~72 ms/step with concurrent renders, pause froze and
  resumed the loop, `keep_viewer_alive` serves after training completes.
  Headless-friendly: `ssh -L 7007:localhost:7007 <box>`.

## ⚠️ Gotchas

- **libtorch interposes `std::filesystem`**: any exe linked against `csrc`
  binds `std::filesystem::remove_all` to libtorch.so's ABI-incompatible copy
  → segfault (jump to null). `-static-libstdc++` does NOT fix it (torch libs
  precede the archive on the link line). Checkpoint pruning uses POSIX
  `nftw` instead. Other fs calls (create_directories, directory_iterator,
  exists) currently bind compatibly but carry the same risk until no-torch.
- `libpython` is linked only because csrc's pybind layer leaves CPython
  symbols undefined; drops out with no-torch.
- Engine auto-resolves `split_batch`+FPBO conflicts via
  `max_input_batch_size` (prints a warning) — pass both through as Python does.
- **Pre-existing**: the process can dump core during exit teardown after a
  completed run (observed on both pre- and post-refactor `ssplat-train`
  2026-07-12; likely a DataManager/engine thread racing static destruction).
  Harmless — all output is already on disk — but worth fixing before
  packaging (Phase 3).

## TODOs (rough priority)

1. **Eval pass + metrics** — iterate `next_val_batch` / render train views,
   PSNR/SSIM from engine buffers; then `validation_fraction` early-stop
   (model config `overfit_score_*`, `early_stop_*` fields are parsed but
   unused).
2. **Resume** — `engine_load_checkpoint` after skeleton setup
   (trainer.py:132-137); config.json round-trip (CLI dump is close to but
   not tyro-compatible; decide on a shared format).
3. ~~ImGui native viewport (Phase 2)~~ DONE 2026-07-12 (`ssplat-gui`, see
   above). Still open within it: CUDA-GL interop upload (currently D2H +
   glTexImage2D — fine at viewport sizes), debug-only `sh` /
   `refinement_score` buffers (engine_debug_forward) not ported, COLMAP
   progress bar is stage-based only, no mesh-export UI (use `ssplat-mesh`).
4. Seeding fidelity: jitter repeated seed points toward a neighbor
   (model.py:550-556) instead of exact duplication; `suppress_initial_scales`
   (model.py:584/700s). (Exact kNN: DONE, `Knn.h`.)
5. `rescale_camera_to_fit` auto-detect (probe image resolution,
   dataparser.py:363-372); COLMAP **text** format fallback. (Metashape
   parser: DONE, `MetashapeParser.cpp`.)
6. Non-default orientation/center methods (`pca`/`vertical`/`gsplat`/`focus`)
   — currently approximated as `up`/`poses` with a warning (only affects
   `train_frame_scale`).
7. Windows: MSVC+nvcc build DONE (2026-07-10, VS2022 + CUDA 12.8 on an
   RTX 3090: no-torch static build links `ssplat-train.exe` and trains
   mipnerf360/garden; only source fix needed was an MSVC branch for a GCC
   atomic builtin in `MeshingHost.cpp`). Remaining: `cudart_static`, CI,
   installer (Phase 3).

## Unsupported-by-design (guarded with clear errors)

`--resume`, `--use-bvh`, `--use-camera-optimizer`,
`--deblur-training-images`, `--optimizer-offload`, `--save-eval-images`,
`--cache-images gpu`,
`--train-frame` ≠ points, `--rescale-camera-to-fit` auto mode,
direct-equirect (`--warp-spherical-to-pinhole 0`) with depth/normal
supervision. `--num-downscales` warns and is ignored (Python-data-path
feature).
