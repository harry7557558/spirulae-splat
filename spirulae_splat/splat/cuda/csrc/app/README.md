# Standalone CLI Trainer (`ssplat-train`)

Working notes for the no-Python trainer under `csrc/app/`. Written 2026-07-10;
covers code structure, verification status, and remaining ports. Long-term
context: this is Phase 1 of GUI-app plan (standalone C++ trainer → Dear
ImGui/GLFW GUI in the same binary → Windows/Linux packaging; kernels stay in
Slang for the eventual Vulkan/cross-vendor port).

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
- Links `libcsrc` (which still links libtorch + libpython on this branch);
  the app code itself has no torch dependency and should survive the
  `260524-no-torch` migration unchanged except for the CMake link line.

## Files

| File | Role |
|---|---|
| `main.cpp` | CLI parsing, engine setup, train loop, checkpointing. Ports of the Python managed path (see mapping below). |
| `DatasetParser.h` | Public dataset structs: `DatasetParserConfig`, `ParsedDataset` (per-INPUT cameras), `PostSplitCameras` + `bake_post_split()` (warp expansion), parse fns, shared `dsparse::` helpers. |
| `ColmapParser.cpp` | COLMAP **binary** reader (text TODO) + format auto-detect dispatcher (`parse_dataset`: transforms.json → nerfstudio, else COLMAP). |
| `NerfstudioParser.cpp` | transforms.json reader + self-contained PLY point reader (ascii + binary_little_endian). |
| `DatasetCommon.cpp` | Shared bakes: normalized-frame scale, eval/val splits, aux-file discovery, geometric-median outlier filter, **`bake_post_split`** (fisheye 5-face / equirect 6-face cubemap expansion; axes MUST match `DataManager.cpp` `kAxesFisheye5`/`kAxesEquirect6`). |
| `Json.h` | Minimal dependency-free JSON parser (handles Python `Infinity`/`NaN`). |
| `HttpServer.h/.cpp` | Minimal HTTP/1.0 GET server (POSIX sockets; winsock shim compiles but untested). Serial request handling — parity with Python's non-threading `HTTPServer`. |
| `Viewer.h/.cpp` | Web-viewer server port (viewer/server.py + http_server.py + render_worker.py + annotation.py): latest-wins render worker, `get_outputs` viewer subset, `engine_blit_view` GPU annotation/colormap, stb JPEG encode. Serves the **unchanged** `viewer.html` (embedded at configure time via CMake hex; `SSPLAT_VIEWER_HTML=<path>` env overrides for dev). |
| `generated/cli_config.h` | AUTO-GENERATED — do not edit. `SsplatConfig` struct (all 189 config fields, defaults baked), `SSPLAT_CONFIG_FIELDS(X)` X-macro flag table, `ssplat_apply_preset()`. |
| `../../../../generate_cli_config.py` | The generator. AST-parses the Python config dataclasses (no torch import → works on fresh checkout). Run by `build_develop.bash`. |
| `../stb_image_write.h` | vendored (public domain, v1.16) for viewer JPEG encoding. |

Debug: `SSPLAT_DUMP_CAMERAS=<path> ssplat-train ...` dumps parsed + post-split
camera arrays as JSON and exits before engine setup — used to diff against
the Python dataparser/trainer algebra (see verification notes below).

## Config codegen (source of truth = Python dataclasses)

`spirulae_splat/generate_cli_config.py` parses `TrainerConfig` (+ preset
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

## Python → C++ port mapping (all in main.cpp unless noted)

| Python | C++ |
|---|---|
| `optimizer.py:54 get_scheduled_lr` | `scheduled_lr()` — pass `nullopt` final when the `*_lr_final` attr doesn't exist in Python (getattr default = constant LR) |
| `model.py:1483 _build_loss_weights` | `build_loss_weights()` — 19 `LossWeightIndex` slots, YUV group normalized to `1-ssim_lambda` |
| `core.py:280 _build_optim_config` | in `build_step_config()` — incl. scale-agnostic mean, `scale_reg/alpha`, quantization level→bits (core.py:95), eps_tr schedule |
| `core.py:321 _build_densify_config` | in `build_step_config()` — `noise_lr_scalar = use_revised ? 1 : alpha`, `max_world_size * alpha` |
| `model.py:1887 engine_train_step_managed` per-step args | in `build_step_config()` — bilagrid/PPISP LRs (AdaGrad switches to constant `*_adagrad_lr`), background noise ramp / SH LRs, color-shift EMA beta |
| `model.py:546 populate_modules` (seeding) | `seed_splats()` — min_init/cap resolution, 4-NN log scales (spatial-hash approx), logit opacity, DC from SfM color incl. gamut conversion |
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

## TODOs (rough priority)

1. **Eval pass + metrics** — iterate `next_val_batch` / render train views,
   PSNR/SSIM from engine buffers; then `validation_fraction` early-stop
   (model config `overfit_score_*`, `early_stop_*` fields are parsed but
   unused).
2. **Resume** — `engine_load_checkpoint` after skeleton setup
   (trainer.py:132-137); config.json round-trip (CLI dump is close to but
   not tyro-compatible; decide on a shared format).
3. **ImGui native viewport (Phase 2)** — reuse `Viewer.cpp`'s render path
   (everything up to the JPEG encode); swap HTTP+JPEG for CUDA-GL interop.
   The web viewer stays as the headless/cloud fallback. Debug-only `sh` /
   `refinement_score` buffers (engine_debug_forward) not ported.
4. Seeding fidelity: jitter repeated seed points toward a neighbor
   (model.py:550-556) instead of exact duplication; `suppress_initial_scales`
   (model.py:584/700s); exact kNN if quality differs.
5. `rescale_camera_to_fit` auto-detect (probe image resolution,
   dataparser.py:363-372); COLMAP **text** format fallback; Metashape
   parser (XML; the JSON/PLY parsers in Json.h / NerfstudioParser.cpp are
   already dependency-free building blocks).
6. Non-default orientation/center methods (`pca`/`vertical`/`gsplat`/`focus`)
   — currently approximated as `up`/`poses` with a warning (only affects
   `train_frame_scale`).
7. Windows: MSVC+nvcc CI build, `cudart_static`, then installer (Phase 3).

## Unsupported-by-design (guarded with clear errors)

`--resume`, `--use-bvh`, `--use-camera-optimizer`,
`--deblur-training-images`, `--optimizer-offload`, `--save-eval-images`,
`--data-format metashape`, `--cache-images gpu`,
`--train-frame` ≠ points, `--rescale-camera-to-fit` auto mode,
direct-equirect (`--warp-spherical-to-pinhole 0`) with depth/normal
supervision. `--num-downscales` warns and is ignored (Python-data-path
feature).
