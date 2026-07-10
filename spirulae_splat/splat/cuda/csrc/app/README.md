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
| `ColmapParser.h/.cpp` | Standalone COLMAP **binary** reader + dataset bake into engine-native arrays (`ParsedDataset`). |
| `generated/cli_config.h` | AUTO-GENERATED — do not edit. `SsplatConfig` struct (all 189 config fields, defaults baked), `SSPLAT_CONFIG_FIELDS(X)` X-macro flag table, `ssplat_apply_preset()`. |
| `../../../../generate_cli_config.py` | The generator. AST-parses the Python config dataclasses (no torch import → works on fresh checkout). Run by `build_develop.bash`. |

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
| `trainer.py:207 _setup_cpp_data_manager` (non-warp) | DataManager setup; batch-size policy from `datamanager.py:91` (non-stochastic) |
| `dataparser.py` COLMAP branch | `ColmapParser.cpp` — w2c→c2w flip (dataparser.py:644), viewmat bake (trainer.py:403), `train_frame_scale` from orient-up/center/auto-scale (scalar only; `train_frame="points"` leaves poses raw), eval_mode splits, aux mask/depth/normal discovery |
| `trainer.py:741 train` / `:939 save_checkpoint` | loop + `save_checkpoint` lambda (step-%09d.ckpt dirs, latest-only pruning) |

## Verified working (banana dataset smoke runs)

- All presets train and converge (3dgs: PSNR 26.7 @ 500 steps); linear-color
  exercises color-space init + trust region; meshing exercises 3dgut + median
  losses; academic-baseline exercises interval eval split + non-FPBO + fp32.
- `360-camera` fails **cleanly** ("warp-to-pinhole not supported") — its
  defining feature is unported; better loud than mistrained.
- Checkpoint cadence + `save_only_latest_checkpoint` pruning + final save.

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

1. **warp_to_pinhole port** — unblocks the `360-camera` preset. Port
   `trainer.py:257-401`: per-camera K (fisheye/equisolid=5, equirect=6),
   post-split camera expansion (cubemap axes must match DataManager.cpp's
   `kAxesFisheye5`/`kAxesEquirect6`), pass K_per_camera/post_offsets/input
   intrins to `engine_setup_data_manager`. Also direct-equirect (K=1)
   canonical intrins (trainer.py:361-369).
2. **Eval pass + metrics** — iterate `next_val_batch` / render train views,
   PSNR/SSIM from engine buffers; then `validation_fraction` early-stop
   (model config `overfit_score_*`, `early_stop_*` fields are parsed but
   unused).
3. **Resume** — `engine_load_checkpoint` after skeleton setup
   (trainer.py:132-137); config.json round-trip (CLI dump is close to but
   not tyro-compatible; decide on a shared format).
4. **Viewer hook** — `engine_viewer_init` + render thread
   (`forward_3dgs`/`engine_blit_view`); this becomes the ImGui viewport in
   Phase 2.
5. Seeding fidelity: jitter repeated seed points toward a neighbor
   (model.py:550-556) instead of exact duplication; `suppress_initial_scales`
   (model.py:584/700s); exact kNN if quality differs.
6. `rescale_camera_to_fit` auto-detect (probe image resolution,
   dataparser.py:363-372); COLMAP **text** format fallback; nerfstudio
   `transforms.json` + Metashape parsers (need a JSON/XML dep decision).
7. Non-default orientation/center methods (`pca`/`vertical`/`gsplat`/`focus`)
   — currently approximated as `up`/`poses` with a warning (only affects
   `train_frame_scale`).
8. Windows: MSVC+nvcc CI build, `cudart_static`, then installer (Phase 3).

## Unsupported-by-design (guarded with clear errors)

`--resume`, `--use-bvh`, `--use-camera-optimizer`,
`--deblur-training-images`, `--optimizer-offload`, `--save-eval-images`,
`--data-format nerfstudio|metashape`, `--cache-images gpu`,
`--train-frame` ≠ points, `--rescale-camera-to-fit` auto mode.
`--num-downscales` warns and is ignored (Python-data-path feature).
