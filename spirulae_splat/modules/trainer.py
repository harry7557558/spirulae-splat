from dataclasses import dataclass, field, asdict
from typing import Type, List, Tuple, Dict, Optional
from pathlib import Path

import threading
import numpy as np
import torch
import time
from tqdm import tqdm
from copy import deepcopy
import cv2

from spirulae_splat.modules.camera import Cameras
from spirulae_splat.modules.model import SpirulaeSplatModelConfig, SpirulaeSplatModel
from spirulae_splat.modules.datamanager import SpirulaeSplatDataManagerConfig, SpirulaeSplatDataManager

from spirulae_splat.modules.dataparser import SpirulaeSplatDataParserConfig
from spirulae_splat.modules.dataset import SpirulaeSplatDataset
from spirulae_splat.modules.optimizer import OptimizerConfig

from spirulae_splat.modules.native_dataparser import (
    parse_dataset as parse_dataset_native,
    to_dataparser_outputs,
    to_native_parser_config,
)
from spirulae_splat.modules.native_trainer import make_session
from spirulae_splat.modules._profile import PROFILE_TRAIN_STEP



@dataclass
class TrainerConfig:
    """Generic method that works well for most datasets."""

    data: Path
    """Path to dataset. Can be a Nerfstudio or a COLMAP dataset."""

    resume: Optional[Path] = None
    """Resume training from a checkpoint. Pass a run output dir (latest
        step-*.ckpt is used) or a specific step-*.ckpt dir. The run's config.json
        supplies the architecture/model/data config (run-control flags like
        num_iterations / save cadence / viewer are still taken from the CLI);
        requires the checkpoint to have been written with save_full_checkpoint."""

    output_dir_prefix: Path = Path("outputs")
    """Prefix to output directory"""

    output_dir_name: Optional[Path] = None
    """Output directory name relative to output_dir_prefix.
        If not specified, will set a generic combining current timestamp and dataset name."""

    steps_per_save: int = 2000
    """Save checkpoint every this number of steps.
        If -1, save only at the end. If zero, never save (used in benchmark)."""
    save_only_latest_checkpoint: bool = True
    """Whether to save only last checkpoint"""
    save_full_checkpoint: bool = False
    """If True, each checkpoint's `state.tar` additionally includes the Resume
        slots -- world raw parameters (at max_num_splats) plus all optimizer
        state -- making the checkpoint sufficient to resume training, not just
        for inference. If False, only the Always slots (appearance/inference
        params) are saved alongside `splat.ply`."""
    save_eval_images: bool = False
    """Whether to save eval images at end of training"""

    num_iterations: int = 30000
    """Number of training iterations"""

    viewer_port: int = 7007
    """Port used by the web viewer"""

    disable_viewer: bool = False
    """If True, ss_trainer skips starting the viewer thread. Used by
        ss_benchmark so each scene runs without competing for the viewer port."""

    keep_viewer_alive: bool = True
    """If True, ss_trainer keeps the process (and thus the viewer) running
        after training + eval finish, so the result can still be inspected in
        the browser. Press Ctrl-C to exit. Ignored when disable_viewer=True."""

    dataparser: SpirulaeSplatDataParserConfig = field(default_factory=SpirulaeSplatDataParserConfig)
    """Specifies configurations for data parsing"""

    datamanager: SpirulaeSplatDataManagerConfig = field(default_factory=SpirulaeSplatDataManagerConfig)
    """Specifies configurations for data management during training"""

    model: SpirulaeSplatModelConfig = field(default_factory=SpirulaeSplatModelConfig)
    """Specifies configurations for main model, losses, and densification"""

    # optimizer: dict = field(default_factory=lambda: _DEFAULT_OPTIMIZERS_WITH_SCALE_SCHEDULER)
    optimizer: OptimizerConfig = field(default_factory=OptimizerConfig)
    """Specifies configurations for optimization"""


class Trainer:
    """Drives a native TrainerSession and adds what stays in Python.

    The session (src/app/TrainerCore.{h,cpp}, bound as `_C.TrainerSession`) owns
    config -> dataset parse -> splat seeding -> engine + DataManager setup ->
    the per-step config and the step itself. It is the same code `ssplat train`
    and the GUI run, so there is no second implementation to drift.

    What this class keeps, and why each has to stay here:
      * the output-dir / config.json conventions -- `ss_trainer.py --resume`
        reads the Python dataclass dump back;
      * checkpoint resume, including layout adaptation (resume_adapt);
      * the eval pass -- LPIPS and the SSIM variants are torch models;
      * the training loop itself, because resume needs a non-zero start step
        and the profiling / debug-dump probes hang off it;
      * `SpirulaeSplatModel`, in attach mode, purely as the eval / metrics
        renderer over the world the session seeded.
    """

    def __init__(
        self,
        config: TrainerConfig
    ):
        # Wipe any C++ EngineState left over from a prior Trainer in this
        # process (e.g. ss_benchmark looping over scenes). Without this the
        # next scene reuses the previous scene's world splats / camera table /
        # bilagrid / PPISP / background / color-space state and silently
        # produces wrong metrics and broken renders.
        from spirulae_splat.modules.core import engine_reset
        engine_reset()

        self.config = config

        # ---- the native session ------------------------------------------
        # Python owns the output dir and config.json, so the session is told
        # not to write its own (different, less complete) dump.
        self.session = make_session(config)
        self.session.write_config_json = False
        self.output_dir = self._setup_output_dir()
        self.session.out_dir_override = str(self.output_dir)
        # setup_engine() logs the output directory itself; don't print it twice.

        self.session.load_dataset()

        # A Python view of the SAME parse -- no second parse for train. Feeds
        # the model (cameras + seed cloud) and dataparser_transforms.json.
        self.dataparser_outputs_train = to_dataparser_outputs(
            self.session.dataset,
            depth_unit_scale_factor=config.dataparser.depth_unit_scale_factor,
            train_frame=config.dataparser.train_frame,
        )
        self._train_frame_scale = float(self.session.dataset.train_frame_scale)
        self._train_to_normalized_transform = \
            self.dataparser_outputs_train['train_to_normalized_transform']

        # Don't overwrite the source run's config.json when resuming.
        if self.config.resume is None:
            self._save_config_json()

        # Seeds the world, installs the DataManager, initializes background /
        # color space / bilagrid / PPISP.
        self.session.setup_engine()

        # ---- the Python-side model, attached to that world ----------------
        self.model = SpirulaeSplatModel(
            self.config,
            self.dataparser_outputs_train['metadata'],
            self.dataparser_outputs_train['cameras'],
            attach_engine=True,
            num_train_data=int(self.session.post.n_post),
        )
        self.model.core.train_frame_scale = self._train_frame_scale
        # Mirror the session's setup decisions onto the model, so resume's
        # eager-init logic and the eval renderer agree about what exists.
        st = self.session.run_state
        self.model._bilagrid_rgb_init    = st.bilagrid_rgb_init
        self.model._bilagrid_depth_init  = st.bilagrid_depth_init
        self.model._bilagrid_normal_init = st.bilagrid_normal_init
        self.model._ppisp_init           = st.ppisp_init

        # The engine skeleton is established -- exactly the precondition
        # engine_load_checkpoint needs. Restore if resuming, else start at 0.
        self._start_step = 0
        if self.config.resume is not None:
            self._resume_from_checkpoint(self.config.resume)

        # Eval state, built on demand by eval().
        self.dataset_eval = None
        self.datamanager = None
        self._viewer = None

        # Progress tracking. The step counter, pause flag, latency window and
        # /progress body all live on the session (shared with the viewer's
        # render worker), so there is nothing to mirror here.
        self.start_time = None
        self._training_time = None

    # ---- viewer ----------------------------------------------------------

    def start_viewer(self, host: str = "0.0.0.0", port: Optional[int] = None):
        """Serve the native web viewer for this session.

        Same HTTP server, render worker and viewer.html the CLI and GUI use
        (src/app/webviewer/). It reads the step counter, pause flag, progress
        JSON and engine mutex straight off the session, so the training loop
        below pushes nothing to it.
        """
        from spirulae_splat.splat.cuda import _C
        if self._viewer is not None:
            return self._viewer
        port = self.config.viewer_port if port is None else port
        self._viewer = _C.WebViewer()
        self._viewer.start_for_session(self.session, host, port)
        print(f"Viewer at http://{host}:{port}/")
        return self._viewer

    def stop_viewer(self):
        if self._viewer is not None:
            self._viewer.stop()
            self._viewer = None

    def toggle_pause(self) -> bool:
        self.session.paused = not self.session.paused
        return self.session.paused

    def get_progress(self):
        import json as _json
        return _json.loads(self.session.progress_json())

    def _setup_output_dir(self):
        if self.config.output_dir_name is not None:
            output_dir = self.config.output_dir_prefix / self.config.output_dir_name
        else:
            import datetime
            timestamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
            dataset_name = self.config.data.stem
            output_dir = self.config.output_dir_prefix / f"{dataset_name}_{timestamp}"
        output_dir.mkdir(parents=True, exist_ok=True)
        return output_dir

    def _save_config_json(self):
        import json

        # config
        config_dict = asdict(
            self.config, dict_factory=lambda data: {
                key: (str(value) if isinstance(value, Path)
                      else value.__name__ if isinstance(value, Type)
                      else value)
                for key, value in data
            })
        with open(self.output_dir / "config.json", "w") as f:
            json.dump(config_dict, f, indent=4)

        # dataparser transform
        dataparser_transform_dict = {
            'transform': self.dataparser_outputs_train['dataparser_transform'][:3, :].tolist(),
            'scale': float(self.dataparser_outputs_train['dataparser_scale']),
        }
        with open(self.output_dir / "dataparser_transforms.json", "w") as f:
            json.dump(dataparser_transform_dict, f, indent=4)

    # ---- training loop ---------------------------------------------------

    def _train_step(self, step: int):
        """One step through the session, plus the Python-only probes.

        The step config (loss weights, LR schedules, densify / bilagrid /
        PPISP / background args) is built in C++ by
        TrainerCore::build_step_config -- there is no Python copy to keep in
        sync. tests/python/test_trainer_parity.py is the gate on that.
        """
        if PROFILE_TRAIN_STEP:
            _prof = getattr(self, "_step_prof", None)
            if _prof is None:
                self._step_prof = _prof = {
                    "n": 0, "warmup": 10, "engine": 0, "gpu_drain": 0, "verbose": 0}
            from time import perf_counter_ns as _t
            t0 = _t()

        self.model.step = step
        with self.session.engine_lock():
            losses = self.session.train_step(step)

        if PROFILE_TRAIN_STEP:
            t1 = _t()
            # GPU drain probe: how much GPU work is still queued after the
            # engine call returns. With async D->H the host runs ahead, so
            # iter N+1's host work overlaps iter N's GPU tail.
            torch.cuda.synchronize()
            t2 = _t()

        self._record_metrics(step, losses)
        self.model.verbose()

        self._maybe_debug_dump(step)

        if PROFILE_TRAIN_STEP:
            t3 = _t()
            if _prof["warmup"] > 0:
                _prof["warmup"] -= 1
            else:
                _prof["n"] += 1
                _prof["engine"]    += (t1 - t0)
                _prof["gpu_drain"] += (t2 - t1)
                _prof["verbose"]   += (t3 - t2)
            if _prof["n"] >= 100:
                n = _prof["n"]
                print(f"\n[PROF n={n}] "
                      f"engine={_prof['engine']/n/1e6:.3f}ms "
                      f"gpu_drain={_prof['gpu_drain']/n/1e6:.3f}ms "
                      f"verbose={_prof['verbose']/n/1e6:.3f}ms",
                      flush=True)
                for k in ("engine", "gpu_drain", "verbose"):
                    _prof[k] = 0
                _prof["n"] = 0

    def _record_metrics(self, step: int, losses: Dict[str, float]):
        """Feed the session's loss dict to the verbose printer.

        `cur_num_splats` is read back from the engine rather than accumulated
        from `num_added`: the session's densify step is the only writer, and
        the counter is a host-side field, so this is both cheaper and immune
        to a missed update.
        """
        from spirulae_splat.splat.cuda import _C
        losses = dict(losses)
        for key in ("num_added", "cur_num_splats", "max_num_splats"):
            losses.pop(key, None)
        for key, value in losses.items():
            self.model.training_verboser.add_metric(key, value)

        self.model.core.cur_num_splats = int(_C.engine_get_cur_num_splats())
        sh_degree_to_use = step // max(self.config.model.sh_degree_warmup_every, 1)
        self.model.core.sh_degree_to_use = sh_degree_to_use
        v = self.model.training_verboser
        v.add_metric("num_splats",     self.model.core.cur_num_splats, last_only=True)
        v.add_metric("max_num_splats", self.model.core.max_num_splats, last_only=True)
        v.add_metric("num_sh",         min(sh_degree_to_use, self.config.model.sh_degree),
                     last_only=True)
        v.add_metric("max_num_sh",     self.config.model.sh_degree, last_only=True)

    @staticmethod
    def _maybe_debug_dump(step: int):
        """Dump GT + render at specific steps, for visual debugging.

        `SS_DEBUG_DUMP_STEPS=0,1,10,100 SS_DEBUG_DUMP_DIR=debug ss-train ...`
        """
        import os
        spec = os.environ.get("SS_DEBUG_DUMP_STEPS", "")
        if not spec:
            return
        try:
            steps = {int(s) for s in spec.split(",") if s.strip()}
        except ValueError:
            return
        if step not in steps:
            return
        from spirulae_splat.modules import debug_image as _dbg
        prefix = f"{os.environ.get('SS_DEBUG_DUMP_DIR', 'debug')}/step{step:06d}"
        for k, a in _dbg.dump_engine_step(prefix).items():
            print(f"[debug_image] {prefix}_{k}: {_dbg.buffer_stats(a)}", flush=True)

    # ---- eval ------------------------------------------------------------

    def _setup_eval_data(self):
        """Parse the eval split and build the Python loader for it.

        Training never touches this path -- the engine's own DataManager feeds
        the training loop. Eval stays in Python because the metrics (LPIPS,
        the SSIM variants) are torch models.
        """
        if self.datamanager is not None:
            return
        ds = parse_dataset_native(
            self.config.data,
            to_native_parser_config(self.config.dataparser, split="eval"))
        dpo = to_dataparser_outputs(
            ds,
            depth_unit_scale_factor=self.config.dataparser.depth_unit_scale_factor,
            train_frame=self.config.dataparser.train_frame)
        self.dataparser_outputs_eval = dpo
        self.dataset_eval = SpirulaeSplatDataset(dpo)

        config = deepcopy(self.config.datamanager)
        config.max_batch_per_epoch = 9**9
        config.load_depths = False
        config.load_normals = False
        config.split_batch = False
        config.patch_batch_size = None
        config.deblur_training_images = False
        config.cache_images = "disk"
        self.datamanager = SpirulaeSplatDataManager(config, device="cuda", eval=True)
        self.datamanager.train_dataset = self.dataset_eval

    def _get_eval_metrics_dict(self):
        inputs = self.datamanager.next_train(0, None)  # type: List[Tuple[Cameras, Dict]]
        assert not isinstance(inputs, tuple)

        for i, (camera, batch) in enumerate(inputs):
            assert i == 0
            model_outputs = self.model.get_outputs(camera)
            metrics_dict, img_dict = self.model.get_image_metrics_and_images(model_outputs, batch)

        return metrics_dict, img_dict

    def get_eval_metrics_dict(self, *args):
        with self.session.engine_lock():
            self.model.eval()
            return self._get_eval_metrics_dict(*args)

    @torch.no_grad()
    def train(self):
        # Engine VRAM peak: pool capacities grow monotonically (high-water
        # mark), so engine_get_pool_breakdown + engine_get_scratch_bytes after
        # the loop equals the training-time peak.
        from spirulae_splat.splat.cuda import _C
        train_wall_start = time.perf_counter()
        self.start_time = time.time()

        sess = self.session
        for step in range(self._start_step, self.config.num_iterations):
            # Pause gate and render fairness, both driven by the session's
            # atomics -- the viewer's render worker sets render_pending before
            # it contends for the engine mutex, and yielding here is what stops
            # the loop's immediate re-acquire from starving it.
            while sess.paused and not sess.stop_requested:
                time.sleep(0.05)
            if sess.stop_requested:
                break
            while sess.render_pending:
                time.sleep(0.0005)

            if step > self._start_step and self.config.steps_per_save > 0 \
                    and step % self.config.steps_per_save == 0:
                self.save_checkpoint(step)

            self._train_step(step)
            sess.cur_step = step + 1

        if self.config.steps_per_save != 0:
            self.save_checkpoint(sess.cur_step)
            print(f"Checkpoint saved to: {self.output_dir.absolute()}")
        self._training_time = time.perf_counter() - train_wall_start
        pool_cap = sum(cap for _, _, cap in _C.engine_get_pool_breakdown())
        self._engine_vram_bytes = pool_cap + _C.engine_get_scratch_bytes()
        print()

    def _train_with_profiling(self):
        def trace_handler(prof: torch.profiler.profile):
            prof.export_chrome_trace(str(self.output_dir / "memprof.json.gz"))
            prof.export_memory_timeline(str(self.output_dir / "memprof.html"), device="cuda:0")
        self.config.num_iterations = 100
        with torch.profiler.profile(
            activities=[
                torch.profiler.ProfilerActivity.CPU,
                torch.profiler.ProfilerActivity.CUDA,
            ],
            schedule=torch.profiler.schedule(wait=0, warmup=0, active=self.config.num_iterations, repeat=1),
            record_shapes=True,
            profile_memory=True,
            with_stack=True,
            on_trace_ready=trace_handler,
        ) as prof:
            for step in range(self.config.num_iterations):
                prof.step()
                with torch.profiler.record_function(str(step)):
                    self._train_step(step)

    @torch.no_grad()
    def eval(self):
        # eval_mode="all" trains on every frame, so there is no held-out set.
        if self.config.dataparser.eval_mode == "all":
            return
        self._setup_eval_data()
        if len(self.dataset_eval.cameras) == 0:
            return

        metrics = {}
        images = {}
        for i in tqdm(range(len(self.dataset_eval.cameras)), desc="Eval", unit="step"):
            metric_dict, image_dict = self.get_eval_metrics_dict()
            for key, value in metric_dict.items():
                if key not in metrics:
                    metrics[key] = []
                metrics[key].append(value)
            if self.config.save_eval_images:
                for key, value in image_dict.items():
                    if key not in images:
                        images[key] = []
                    images[key].append(value)
        for key, value in [*metrics.items()]:
            value = sum(value) / len(value)
            metrics['avg_'+key] = value
            print(f"{key}: {value}")

        if self.config.save_eval_images:
            for key, images in images.items():
                for idx, image in enumerate(images):
                    path = self.output_dir / f"eval-{key}-{idx:05d}.png"
                    image = (torch.clip(image, 0, 1) * 255).to(torch.uint8).cpu().numpy()
                    cv2.imwrite(str(path), cv2.cvtColor(image, cv2.COLOR_RGB2BGR))

        # Per-run training stats: subprocess-launched benchmarks read these
        # back out of metrics.json since they can't observe the in-process
        # trainer state.
        if hasattr(self, "_training_time"):
            metrics["training_time"] = self._training_time
        if hasattr(self, "_engine_vram_bytes"):
            metrics["engine_vram"] = self._engine_vram_bytes / 1024 ** 2

        import json
        with open(self.output_dir / "metrics.json", "w") as f:
            json.dump(metrics, f, indent=4)

    def _resume_from_checkpoint(self, load_dir) -> None:
        """Restore engine + training state from a full checkpoint. Called from
        __init__ after the engine skeleton is built (world seeded at
        max_num_splats, data manager installed).

        Supports resuming with a DIFFERENT layout than the checkpoint: fewer
        Gaussians (new smaller cap_max), a different SH degree, and with/without
        bilagrid/PPISP. When the layouts differ, the checkpoint's buffers are
        adapted on the CPU (resume_adapt: splat reduction, SH resample, quant
        re-encode, bilagrid resample) to the target layout before loading -- no
        extra VRAM. Appearance channels that both the checkpoint AND the current
        config want are eagerly re-initialized here (at the config's dims) so
        their pool slots exist as restore targets and the per-step lazy init
        won't re-zero them."""
        import shutil, tempfile
        from spirulae_splat.modules.resume import (
            resolve_checkpoint, read_state_json, check_resumable)
        from spirulae_splat.modules.resume_adapt import adapt_checkpoint

        run_dir, ckpt_dir = resolve_checkpoint(Path(load_dir))
        check_resumable(ckpt_dir)          # not-full checkpoint -> clear error
        state = read_state_json(ckpt_dir)

        model = self.model
        cfg = self.config.model
        opt = self.config.optimizer
        n_img = model.num_train_data

        # --- target layout: what the engine ACTUALLY holds ------------------
        # Taken from the session's RunState (mirrored onto the model in
        # __init__), not re-derived from the config. Deriving it from config
        # alone got this wrong: the depth / normal grids also require the
        # dataset to *have* depth / normal maps, which setup_engine() checks
        # and this did not. On any dataset without normals -- Mip-NeRF 360,
        # say -- the target claimed a normal bilagrid the checkpoint never
        # had, so every resume took the full CPU adaptation path (unpack,
        # resample, repack) to "fix" a difference that did not exist.
        def _lhw(shape):                       # (grid_X, grid_Y, grid_W) -> (L,H,W)
            X, Y, W_g = shape
            return (W_g, Y, X)
        target = {
            "max_num_splats": int(model.core.max_num_splats),
            "num_sh": int(model.gauss_params["features_sh"].shape[1]),
            "num_images": int(n_img),
            "bilagrid": {
                "rgb":    _lhw(cfg.bilagrid_shape) if model._bilagrid_rgb_init else None,
                "depth":  _lhw(cfg.bilagrid_shape_geometry) if model._bilagrid_depth_init else None,
                "normal": _lhw(cfg.bilagrid_shape_geometry) if model._bilagrid_normal_init else None,
            },
            "ppisp": bool(model._ppisp_init),
        }

        # --- adapt (CPU) if the layout differs, then load -------------------
        tmp = tempfile.mkdtemp(prefix="ss_resume_")
        try:
            adapted = adapt_checkpoint(ckpt_dir, target, tmp)
            load_from = adapted if adapted is not None else ckpt_dir
            step = model.core.engine_load_checkpoint(load_from)
            eff_state = read_state_json(load_from)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

        # Refresh host gauss_params from the (restored) device world + step.
        model.core.engine_sync_splats_to_host()
        model.core.cur_num_splats = int(eff_state.get("cur_num_splats",
                                                      model.core.cur_num_splats))
        model.step = step
        self._start_step = step
        self.current_step = step
        adapted_note = " (layout-adapted)" if adapted is not None else ""
        print(f"Resumed from {ckpt_dir} at step {step}{adapted_note} "
              f"(cur_num_splats={model.core.cur_num_splats})")

    def save_checkpoint(self, step: int) -> None:
        # Checkpoints are now a directory (PLY + npy on the C++ side);
        # torch.save(model.state_dict()) no longer captures the engine buffers.
        ckpt_path: Path = self.output_dir / f"step-{step:09d}.ckpt"
        self.model.core.engine_save_checkpoint(
            ckpt_path, step,
            full_dump=self.config.save_full_checkpoint,
        )
        if self.config.save_only_latest_checkpoint:
            import shutil
            for f in self.output_dir.glob("step-*.ckpt"):
                if f != ckpt_path:
                    if f.is_dir():
                        shutil.rmtree(f, ignore_errors=True)
                    else:
                        f.unlink()

    def print_vram_breakdown(self):
        from spirulae_splat.splat.cuda import _C

        name_map = {}
        def add_tensor(key, value):
            if isinstance(value, list) or isinstance(value, tuple):
                for i, v in enumerate(value):
                    add_tensor(f"{key}[{i}]", v)
                return
            if not isinstance(value, torch.Tensor):
                return
            if value.data_ptr() not in name_map:
                name_map[value.data_ptr()] = []
            name_map[value.data_ptr()].append(key)
        for key, value in self.model.core.__dict__.items():
            add_tensor(f"core.{key}", value)
        for key, value in self.model.__dict__.items():
            add_tensor(f"model.{key}", value)

        import gc
        ptr_map = set()
        breakdown = []
        for obj in gc.get_objects():
            try:
                if isinstance(obj, torch.Tensor):
                    if obj.is_cuda:
                        if obj.data_ptr() in ptr_map:
                            continue
                        ptr_map.add(obj.data_ptr())
                        breakdown.append((obj.nbytes, obj.shape, str(obj.dtype), str(obj.device), obj.data_ptr()))
            except:
                pass
        breakdown.sort(reverse=True)
        torch_total = 0
        print("=== PyTorch tensors ===")
        for nbytes, shape, dtype, device, ptr in breakdown:
            print(f"  {nbytes/1024**2:8.2f} MiB  {str(shape):24s} {dtype:16s}  {' '.join(name_map.get(ptr, []))}")
            torch_total += nbytes

        print(f"\n=== C++ DevicePool ===")
        pool_breakdown = _C.engine_get_pool_breakdown()
        pool_breakdown.sort(key=lambda x: -x[2])
        pool_total_used = 0
        pool_total_cap = 0
        for key, used, cap in pool_breakdown:
            print(f"  {cap/1024**2:8.2f} MiB  (used {used/1024**2:.2f})  {key}")
            pool_total_used += used
            pool_total_cap += cap
        scratch_bytes = _C.engine_get_scratch_bytes()
        if scratch_bytes > 0:
            print(f"  {scratch_bytes/1024**2:8.2f} MiB  DeviceScratch")
            pool_total_cap += scratch_bytes

        print(f"\n=== Summary ===")
        print(f"  PyTorch tensors:          {torch_total / 1024**2:8.2f} MiB")
        print(f"  C++ pool (capacity):      {pool_total_cap / 1024**2:8.2f} MiB")
        print(f"  torch.cuda.mem_allocated: {torch.cuda.memory_allocated() / 1024**2:8.2f} MiB")
        print(f"  torch.cuda.mem_reserved:  {torch.cuda.memory_reserved() / 1024**2:8.2f} MiB")
        free, total = torch.cuda.mem_get_info()
        print(f"  GPU memory in use:        {(total - free) / 1024**2:8.2f} MiB")
        print()



@dataclass
class TrainerConfig360Camera(TrainerConfig):
    """Preset for training on original distorted images captured by 360 cameras (e.g. Insta360, DJI Osmo). Recommended if your dataset contains fisheye images with a circle visible."""
    datamanager: SpirulaeSplatDataManagerConfig = field(default_factory=lambda: SpirulaeSplatDataManagerConfig(
        warp_to_pinhole=True,
        mask_boundary_offset=-0.025,
    ))
    model: SpirulaeSplatModelConfig = field(default_factory=lambda: SpirulaeSplatModelConfig(
        primitive="mip",
        # densify_score_mode="median",
        # densify_loss_map_mode="robust_edge_aware",
        # densify_robust_edge_aware_quantile=0.99,
        long_axis_split_opacity_k=(0.5, 0.6, 15000),
        input_depth_is_ray_depth=True,  # generated by predict_geometry.py
    ))


@dataclass
class TrainerConfigInTheWild(TrainerConfig):
    """Preset for datasets consisting of internet images, with extreme lighting variation, with un-masked outliers, and/or shot with long focal lengths."""
    dataparser: SpirulaeSplatDataParserConfig = field(default_factory=lambda: SpirulaeSplatDataParserConfig(
        center_method="focus",
        outlier_threshold=10.0,
    ))
    datamanager: SpirulaeSplatDataManagerConfig = field(default_factory=lambda: SpirulaeSplatDataManagerConfig(
        load_depths=True,
        load_normals=True,
        mask_boundary_offset=-0.025,
    ))
    model: SpirulaeSplatModelConfig = field(default_factory=lambda: SpirulaeSplatModelConfig(
        densify_score_mode="median",
        densify_loss_map_mode="robust_edge_aware",
        densify_robust_edge_aware_quantile=0.75,
        # densify_score_blend_world_grad=0.1,
        ssim_lambda=0.1,
        rgb_distortion_reg=0.1,
        depth_distortion_reg=0.01,
        sh_degree_warmup_every=0,
        long_axis_split_opacity_k=(0.5, 0.6, 30000),
        noise_lr=10.0,
        noise_lr_final=0.1,
        erank_reg=0.1,
    ))
    optimizer: OptimizerConfig = field(default_factory=lambda: OptimizerConfig(
        means_lr=5.0e-5,
        means_lr_final=1.0e-7,
    ))


@dataclass
class TrainerConfigLinear(TrainerConfig):
    """Preset for training splats in linear color spaces (e.g. ACEScg)."""
    model: SpirulaeSplatModelConfig = field(default_factory=lambda: SpirulaeSplatModelConfig(
        splat_color_gamut="ACEScg",  # configurable
        splat_color_is_linear=True,
        image_color_gamut="Rec.2020",  # configurable
        image_color_is_linear=False,  # configurable
        background_mode="noise",
    ))


@dataclass
class TrainerConfigSynthetic(TrainerConfig):
    """Preset for training splats on synthetic datasets rendered with constant exposure."""
    model: SpirulaeSplatModelConfig = field(default_factory=lambda: SpirulaeSplatModelConfig(
        min_init_fraction=0.1,
        use_bilateral_grid=False,
        use_ppisp=False,
        use_bilateral_grid_for_geometry=False,
        long_axis_split_opacity_k=(0.5, 0.6, 25000),
    ))


@dataclass
class TrainerConfigMeshing(TrainerConfig):
    """Preset for training splats for meshing. Use `spirulae-meshing` to convert trained splats to mesh."""
    model: SpirulaeSplatModelConfig = field(default_factory=lambda: SpirulaeSplatModelConfig(
        primitive="3dgut",
        sh_degree=0,
        sh_reg=10.0,
        overexposure_reg=10.0,
        background_mode="noise",
        depth_distortion_reg=0.01,
        normal_distortion_reg=0.01,
        mean_median_depth_weight=0.01,
        median_depth_normal_reg_weight=0.01,
        normal_supervision_weight=0.01,
        median_normal_supervision_weight=0.01,
        median_render_normal_reg_weight=0.01,
        erank_reg=0.01,
        erank_reg_s3=0.01,
    ))


@dataclass
class TrainerConfigAcademicBaseline(TrainerConfig):
    """Preset that replicates 3DGS MCMC as faithful as possible."""
    dataparser: SpirulaeSplatDataParserConfig = field(default_factory=lambda: SpirulaeSplatDataParserConfig(
        eval_mode="interval",
        eval_interval=8,
        center_method="gsplat",
        orientation_method="gsplat",
    ))
    datamanager: SpirulaeSplatDataManagerConfig = field(default_factory=lambda: SpirulaeSplatDataManagerConfig(
        max_batch_per_epoch=9**9,
        load_depths=False,
        load_normals=False,
    ))
    model: SpirulaeSplatModelConfig = field(default_factory=lambda: SpirulaeSplatModelConfig(
        primitive="3dgs",
        relative_scale=1.0,
        use_bilateral_grid=False,
        use_bilateral_grid_for_geometry=False,
        use_ppisp=False,
        use_revised_densification=False,
        densify_loss_map_mode="none",
        use_long_axis_split=False,
        use_fused_proj_bwd_optim=False,
        quantization_level=0,
        max_screen_size=float('inf'),
        max_world_size=float('inf'),
        suppress_initial_scales=False,
        scale_init=0.1,
        opacity_init=0.5,
        depth_distortion_reg=0.0,
        normal_reg_weight=0.0,
        alpha_reg_weight=0.0,
        alpha_loss_weight=0.0,
        alpha_loss_weight_under=0.0,
        erank_reg=0.0,
        erank_reg_s3=0.0,
        quat_norm_reg=0.0,
        sh_reg=0.0,
        normal_supervision_weight=0.0,
        opacity_reg=0.01,
        scale_reg=0.01,
    ))
    # optimizer: dict = field(default_factory=lambda: _DEFAULT_OPTIMIZERS)
    optimizer: OptimizerConfig = field(default_factory=lambda: OptimizerConfig(
        max_steps=30000,
        use_scale_agnostic_mean=False,
        use_per_splat_bias_correction=False,
        means_lr=1.6e-4,
        means_lr_final=1.6e-6,
        scales_lr=0.005,
        scales_lr_final=None,
        quats_lr=0.001,
        opacities_lr=0.05,
        features_dc_lr=0.0025,
        features_sh_lr=0.0025 / 20,
    ))
