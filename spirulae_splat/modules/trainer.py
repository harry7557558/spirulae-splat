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

from spirulae_splat.modules.dataparser import SpirulaeSplatDataparser, SpirulaeSplatDataParserConfig
from spirulae_splat.modules.dataset import SpirulaeSplatDataset
from spirulae_splat.modules.optimizer import OptimizerConfig

from spirulae_splat.viewer.annotation import annotate_train_cameras
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
        self.dataparser = SpirulaeSplatDataparser(config.dataparser, config.data)
        self.dataparser_outputs_train, self.dataparser_outputs_eval = self.dataparser.parse()
        self.dataset_train = SpirulaeSplatDataset(self.dataparser_outputs_train)
        self.dataset_eval = SpirulaeSplatDataset(self.dataparser_outputs_eval)

        self._train_frame_scale = float(
            self.dataparser_outputs_train.get('train_frame_scale', 1.0))
        self._train_to_normalized_transform = self.dataparser_outputs_train.get(
            'train_to_normalized_transform', torch.eye(4))

        self.datamanager = SpirulaeSplatDataManager(self.config.datamanager, device="cuda")
        self.datamanager.train_dataset = self.dataset_train

        self.model = SpirulaeSplatModel(
            self.config,
            self.dataparser_outputs_train['metadata'],
            self.dataparser_outputs_train['cameras']
        )
        # Push the frame scale into the engine wrapper so it can rescale
        # means_lr / scale_reg / max_world_size / noise_lr per step.
        self.model.core.train_frame_scale = self._train_frame_scale

        # When the C++ DataManager is enabled, install it on the engine once
        # here. Per-step training will then go through
        # ``model.engine_train_step_managed`` and bypass the Python data path.
        self._setup_cpp_data_manager()

        # The engine skeleton (world seeded at max_num_splats + data manager) is
        # now established -- exactly the precondition engine_load_checkpoint
        # needs. Restore from a checkpoint if resuming; otherwise start at step 0.
        self._start_step = 0
        if self.config.resume is not None:
            self._resume_from_checkpoint(self.config.resume)

        self.output_dir = self._setup_output_dir()
        print(f"Output directory: {self.output_dir.absolute()}")

        # Don't overwrite the source run's config.json when resuming.
        if self.config.resume is None:
            self._save_config_json()

        self.lock = threading.Lock()
        # Lock-fairness signal: the render thread sets this before contending
        # for `self.lock`. The training thread checks it at each iteration's
        # acquire site and yields the OS scheduler enough times for the render
        # thread to win the contended acquire. Without this, train_step's
        # tight re-acquire after release would starve render for many
        # iterations (observed: 7-30 s waits on a 200 ms/step run).
        self._render_pending = threading.Event()

        # Progress tracking
        self.current_step = 0
        self.start_time = None
        self.last_step_time = None
        self.step_latencies = []
        
        # Pause/resume control
        self._paused = False
        self._pause_lock = threading.Lock()
        self._pause_event = threading.Event()
        self._pause_event.set()  # Initially not paused (set = can proceed)

    def toggle_pause(self) -> bool:
        """Toggle pause state. Returns True if now paused, False if resumed."""
        with self._pause_lock:
            self._paused = not self._paused
            if self._paused:
                self._pause_event.clear()  # Paused: block training
            else:
                self._pause_event.set()  # Resumed: allow training
            return self._paused
    
    def _check_pause(self):
        """Wait if training is paused. Called during training loop."""
        self._pause_event.wait()

    def get_progress(self):
        if self.start_time is None:
            return {
                "step": 0,
                "total_steps": self.config.num_iterations,
                "elapsed_time": 0,
                "eta": None,
                "latency_ms": None,
                "paused": False,
            }
        elapsed = time.time() - self.start_time
        avg_latency = sum(self.step_latencies) / len(self.step_latencies) if self.step_latencies else None
        if avg_latency and self.current_step > 0:
            remaining_steps = self.config.num_iterations - self.current_step
            eta = remaining_steps * avg_latency
        else:
            eta = None
        return {
            "step": self.current_step,
            "total_steps": self.config.num_iterations,
            "elapsed_time": elapsed,
            "eta": eta,
            "latency_ms": avg_latency * 1000 if avg_latency else None,
            "paused": self._paused,
        }

    def _setup_cpp_data_manager(self):
        """Install the C++ DataManager on the engine from the parsed dataset.

        Bakes per-camera viewmats / intrins / dist_coeffs into the same flat
        layout the engine consumes — same R/T inverse + Y/Z flip + relative
        scale that ``model.engine_train_step`` would have applied per step,
        but done ONCE here so the C++ side can hand them to the engine
        directly each step.

        The managed path is currently restricted to plain, no-frills setups:
        no camera optimizer, no warp_to_pinhole, no equirectangular split,
        no split_batch, no BVH primitive.
        Anything else raises here so the user sees the limitation up front
        rather than silently quality-regressing.
        """
        import spirulae_splat.csrc as _C

        dm_cfg     = self.config.datamanager
        model_cfg  = self.config.model
        dpo        = self.dataparser_outputs_train

        # ---- Feature-combination validation ------------------------------
        # datamanager.split_batch (the Python data-path OOM workaround) is a
        # no-op here: the C++ DataManager builds a deterministic once-per-epoch
        # schedule and, on datasets with many small resolution groups, packs
        # sub-B remainders across resolutions into a single optimizer step that
        # the engine consumes via heterogeneous grad accumulation. Effective
        # batch size (and hence FPBO vs. grad-accumulator VRAM) is derived from
        # max_batch_per_epoch, so no extra flag is needed.
        if self.model.core.primitive not in ('3dgs', 'mip', '3dgut') or self.model.core.use_bvh:
            raise NotImplementedError("use_cpp_data_manager requires non-BVH 3dgs/mip/3dgut primitive")

        cameras = dpo['cameras']
        N = len(cameras)
        if N == 0:
            raise RuntimeError("use_cpp_data_manager: empty cameras")

        # Per-camera model enum int. Mixed models OK -- the C++ DataManager
        # partitions by (W, H, model).
        cam_types = [str(t).upper() for t in cameras.camera_type]
        camera_models = []
        for s in cam_types:
            e = _C.camera_model_from_name(s) if hasattr(_C, 'camera_model_from_name') \
                else getattr(_C.CameraModelType, s)
            camera_models.append(int(e))

        # Per-camera input shape.
        widths  = [int(w) for w in cameras.width.detach().cpu().flatten().tolist()]
        heights = [int(h) for h in cameras.height.detach().cpu().flatten().tolist()]

        # ---- Warp-to-pinhole split factor per camera --------------------
        # FISHEYE / EQUISOLID -> K=5 when warp_to_pinhole is on.
        # EQUIRECTANGULAR    -> K=6 always (regardless of flag).
        # everything else    -> K=1 (pass-through).
        from spirulae_splat.modules.camera import CameraType
        PINHOLE_VAL  = int(_C.camera_model_from_name("PINHOLE"))
        FISHEYE_VAL  = int(_C.camera_model_from_name("FISHEYE"))
        EQUISOLID_VAL = int(_C.camera_model_from_name("EQUISOLID"))
        EQUIRECT_VAL  = int(_C.camera_model_from_name("EQUIRECTANGULAR"))
        K_per_camera = []
        for cm in camera_models:
            if cm == EQUIRECT_VAL:
                K_per_camera.append(6)
            elif dm_cfg.warp_to_pinhole and cm in (FISHEYE_VAL, EQUISOLID_VAL):
                K_per_camera.append(5)
            else:
                K_per_camera.append(1)
        any_warp = any(k > 1 for k in K_per_camera)
        post_offsets = []
        acc = 0
        for k in K_per_camera:
            post_offsets.append(acc); acc += k
        n_post = acc

        # Warp path cannot supply depth / normal supervision (the warp kernels
        # operate only on RGB + binary mask). PPISP IS supported on warp: it's
        # render-side and runs on the POST-split pinhole sub-cameras, with one
        # parameter slot per sub-camera (n_post total).
        if any_warp:
            if dm_cfg.load_depths and dpo['metadata'].get('depth_filenames', None):
                raise NotImplementedError(
                    "use_cpp_data_manager warp (fisheye/equirectangular split) "
                    "does not support depth supervision yet.")
            if dm_cfg.load_normals and dpo['metadata'].get('normal_filenames', None):
                raise NotImplementedError(
                    "use_cpp_data_manager warp does not support normal supervision yet.")

        # ---- Per-INPUT c2w / intrins / dist_coeffs (raw, host) ----------
        c2w_all = cameras.camera_to_worlds.detach().cpu()
        if c2w_all.dim() == 2:
            c2w_all = c2w_all.unsqueeze(0)
        c2w_all = c2w_all.float()  # [N, 3, 4]

        input_intrins = cameras.intrins.detach().cpu().float()
        if input_intrins.dim() == 1:
            input_intrins = input_intrins.unsqueeze(0).repeat(N, 1)
        if cameras.distortion_params is None:
            input_dist_coeffs = torch.zeros(N, 10, dtype=torch.float32)
        else:
            input_dist_coeffs = cameras.distortion_params.detach().cpu().float()
            if input_dist_coeffs.dim() == 1:
                input_dist_coeffs = input_dist_coeffs.unsqueeze(0).repeat(N, 1)
            if input_dist_coeffs.shape[1] < 10:
                pad = torch.zeros(N, 10 - input_dist_coeffs.shape[1], dtype=torch.float32)
                input_dist_coeffs = torch.cat([input_dist_coeffs, pad], dim=1)

        # ---- Expand to POST-split c2w / intrins / dist_coeffs ------------
        # Same algebra as modules/resample.py warp_*_to_pinhole, but on cameras
        # only (no images): apply diag(1, -1, -1) flip on c2w R, rotate by the
        # cubemap axis, flip back; copy T. Output intrins are canonical
        # PINHOLE at fx=fy=cx=cy=out_shape/2; dist_coeffs are all zero.
        import math
        # 5- and 6-face cubemap axis tables (must match the C++ axes in
        # DataManager.cpp's kAxesFisheye5 / kAxesEquirect6).
        axes_5 = torch.tensor([
            [[1,0,0],[0,1,0],[0,0,1]],
            [[0,1,0],[0,0,1],[1,0,0]],
            [[-1,0,0],[0,0,1],[0,1,0]],
            [[0,-1,0],[0,0,1],[-1,0,0]],
            [[1,0,0],[0,0,1],[0,-1,0]],
        ], dtype=torch.float32)
        axes_6 = torch.tensor([
            [[1,0,0],[0,1,0],[0,0,1]],
            [[0,0,-1],[0,1,0],[1,0,0]],
            [[-1,0,0],[0,0,1],[0,1,0]],
            [[0,-1,0],[0,0,1],[-1,0,0]],
            [[1,0,0],[0,0,1],[0,-1,0]],
            [[1,0,0],[0,-1,0],[0,0,-1]],
        ], dtype=torch.float32)
        D = torch.tensor([1.0, -1.0, -1.0])  # diag flip applied per face

        post_c2w     = torch.zeros(n_post, 3, 4, dtype=torch.float32)
        post_intrins = torch.zeros(n_post, 4, dtype=torch.float32)
        post_dist    = torch.zeros(n_post, 10, dtype=torch.float32)
        post_widths        = [0] * n_post
        post_heights       = [0] * n_post
        post_camera_models = [0] * n_post   # int enum; viewer uses these
        for i in range(N):
            k_i = K_per_camera[i]
            off = post_offsets[i]
            ci  = c2w_all[i]  # [3, 4]
            if k_i == 1:
                post_c2w[off]     = ci
                post_intrins[off] = input_intrins[i]
                post_dist[off]    = input_dist_coeffs[i]
                post_widths[off]        = widths[i]
                post_heights[off]       = heights[i]
                post_camera_models[off] = camera_models[i]
                continue
            # Pick the right axes table.
            axes = axes_6 if k_i == 6 else axes_5
            out_shape = int(math.ceil(math.sqrt(heights[i] * widths[i] / k_i)))
            # Same algebra as modules/resample.py:
            #   R' = R_orig * diag(1, -1, -1)             # OpenGL -> OpenCV flip
            #   R_warp[k] = R' @ axes[k]^T                # cubemap rotation in OpenCV
            #   R_post[k] = R_warp[k] * diag(1, -1, -1)   # OpenCV -> OpenGL flip
            # The earlier `axes[k] @ R` form was wrong (left-multiply instead
            # of right-multiply by axes^T); for face 0 it happened to give
            # the right answer (axes[0] = I), so the bug only showed up for
            # the non-identity faces -- exactly the symptom the user saw
            # (a blurry / messy render after warp).
            R = ci[:3, :3] * D[None, :]                  # [3, 3]
            for k in range(k_i):
                R_new = (R @ axes[k].T) * D[None, :]     # [3, 3]
                post_c2w[off + k, :3, :3] = R_new
                post_c2w[off + k, :3, 3]  = ci[:3, 3]
                # Canonical pinhole intrins for the sub-image.
                s = float(out_shape) / 2.0
                post_intrins[off + k] = torch.tensor([s, s, s, s])
                # Dist coeffs already zero.
                post_widths[off + k]        = out_shape
                post_heights[off + k]       = out_shape
                # Each post-split face is rendered/displayed as PINHOLE.
                post_camera_models[off + k] = PINHOLE_VAL

        # ---- c2w -> viewmat (R/T inverse + Y/Z flip + relative_scale) ---
        # Matches model.engine_train_step's per-step formula.
        R_v = post_c2w[:, :3, :3].clone()
        T_v = post_c2w[:, :3, 3:4].clone()
        if model_cfg.relative_scale is not None:
            T_v = T_v * float(model_cfg.relative_scale)
        R_v = R_v * torch.tensor([[[1.0, -1.0, -1.0]]])
        R_inv = R_v.transpose(-1, -2)
        T_inv = -torch.bmm(R_inv, T_v)
        viewmats = torch.eye(4).unsqueeze(0).repeat(n_post, 1, 1).float()
        viewmats[:, :3, :3] = R_inv
        viewmats[:, :3, 3:4] = T_inv

        intrins     = post_intrins
        dist_coeffs = post_dist

        # ---- File-path lists --------------------------------------------
        def _strs(seq):
            return [str(p) for p in seq] if seq is not None else []

        image_filenames  = _strs(dpo['image_filenames'])
        mask_filenames   = _strs(dpo.get('mask_filenames', None))
        depth_filenames  = _strs(dpo['metadata'].get('depth_filenames',  None)) if dm_cfg.load_depths  else []
        normal_filenames = _strs(dpo['metadata'].get('normal_filenames', None)) if dm_cfg.load_normals else []

        # ---- Train / val split -------------------------------------------
        val_set = set(int(i) for i in self.dataset_train.val_indices)
        train_indices = [i for i in range(N) if i not in val_set]
        val_indices   = sorted(val_set)

        # Effective batch sizes — reuse the Python datamanager's policy so
        # behavior matches at iteration scale.
        train_bs = max(1, self.datamanager.train_batch_size(stochastic=False))
        val_bs   = max(1, self.datamanager.val_batch_size(stochastic=False)) if val_indices else 1

        # ---- Build the DataManagerConfig --------------------------------
        cmap = {
            "cpu":          _C.DataManagerCacheMode.CPU,
            "cpu-pageable": _C.DataManagerCacheMode.CPU,
            "disk":         _C.DataManagerCacheMode.DISK,
        }
        cache_key = dm_cfg.cache_images
        if cache_key not in cmap:
            raise NotImplementedError(
                f"use_cpp_data_manager: cache_images='{cache_key}' (only "
                "'cpu', 'cpu-pageable', or 'disk' are supported)")

        c_cfg = _C.DataManagerConfig()
        c_cfg.cache_mode       = cmap[cache_key]
        # Enable masks also when warp_to_pinhole is on AND the dataset has
        # any fisheye/equisolid input -- the C++ DataManager synthesizes
        # a 1x1 all-white placeholder per such image so the wide-warp mask
        # kernel produces a proper post-split FOV mask (1 inside lens, 0
        # outside). Without a mask the unseen regions of each cubemap
        # face stay black and the training loss bakes gray splats there.
        _needs_synth_mask = dm_cfg.warp_to_pinhole and any(
            cm in (FISHEYE_VAL, EQUISOLID_VAL) for cm in camera_models
        )
        c_cfg.load_masks       = (len(mask_filenames) > 0) or _needs_synth_mask
        c_cfg.load_depths      = (len(depth_filenames) > 0 and not any_warp)
        c_cfg.load_normals     = (len(normal_filenames) > 0 and not any_warp)
        c_cfg.train_batch_size = train_bs
        c_cfg.val_batch_size   = val_bs
        c_cfg.warp_to_pinhole  = bool(dm_cfg.warp_to_pinhole)
        c_cfg.mask_boundary_offset = float(dm_cfg.mask_boundary_offset)

        # Input intrins / dist_coeffs are needed by the wide warp kernel
        # (fisheye / equisolid). Pass them even when only equirectangular
        # warp is active -- the equirectangular kernel just ignores them.
        input_intrins_list     = input_intrins.contiguous().flatten().tolist() if any_warp else []
        input_dist_coeffs_list = input_dist_coeffs.contiguous().flatten().tolist() if any_warp else []

        _C.engine_setup_data_manager(
            c_cfg, camera_models,
            image_filenames, mask_filenames, depth_filenames, normal_filenames,
            widths, heights,
            K_per_camera, post_offsets,
            viewmats.contiguous().flatten().tolist(),
            intrins.contiguous().flatten().tolist(),
            dist_coeffs.contiguous().flatten().tolist(),
            input_intrins_list, input_dist_coeffs_list,
            train_indices, val_indices)

        # The bilagrid (+ PPISP, + camera-optimizer) tables are sized off
        # `self.model.num_train_data`. The model's pre-init uses a uniform
        # K=6 upper bound for the C++ datamanager path since it has no
        # per-camera K info at __init__ time; resolve it to the real
        # post-split count now. For non-warp datasets (n_post == N), this
        # shrinks the table back to N -- without this, bilagrid `n_grids`
        # is 6x too large and the TV-loss normalization (~1/N_grids) makes
        # the TV regularizer 6x weaker than intended.
        self.model.num_train_data = int(n_post)

        # Build a POST-split Cameras object for the viewer. For the warp
        # path each input image is exposed as K=5 (fisheye/equisolid) or
        # K=6 (equirectangular) PINHOLE sub-cameras at the cubemap face
        # poses -- exactly matching what the engine renders + caches
        # thumbnails for. For mixed datasets, the per-camera K table
        # already handles the interleave (K=1 pass-through for native
        # PINHOLE, K=5/6 split for the wide / equirect entries).
        _INT_TO_NAME = {
            PINHOLE_VAL:  "PINHOLE",
            FISHEYE_VAL:  "FISHEYE",
            EQUISOLID_VAL:"EQUISOLID",
            EQUIRECT_VAL: "EQUIRECTANGULAR",
        }
        post_camera_type_names = [_INT_TO_NAME[m] for m in post_camera_models]
        self._viewer_cameras = Cameras(
            intrins           = post_intrins,
            distortion_params = post_dist,
            height            = torch.tensor(post_heights, dtype=torch.int32),
            width             = torch.tensor(post_widths,  dtype=torch.int32),
            camera_to_worlds  = post_c2w,
            camera_type       = post_camera_type_names,
        )

        # Pass through to the model so engine_train_step_managed knows which
        # bilagrid lazy-init branches are actually live this run.
        self.model._dm_has_depths  = c_cfg.load_depths
        self.model._dm_has_normals = c_cfg.load_normals

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

    def _render(self, c2w, fx, fy, cx, cy, w, h, camera_model, buffer_key="rgb"):
        if self._train_frame_scale != 1.0:
            # Viewer client sends c2w in the legacy normalized frame; remap to
            # the actual training frame before rendering so navigation feel
            # (movement speed, axis-up) stays the same regardless of train_frame.
            # T is a similarity (scale * R | t); apply the full similarity to
            # the camera position but only the unit rotation to the c2w R block
            # — otherwise the scale would be baked into the camera basis
            # vectors and the camera would over-zoom by a factor of scale.
            T = self._train_to_normalized_transform.detach().cpu().numpy().astype(c2w.dtype)
            R_full = T[:3, :3]
            t_full = T[:3, 3]
            s = float(np.linalg.norm(R_full[:, 0]))
            R_unit = R_full / s
            R_in = c2w[:3, :3]
            t_in = c2w[:3, 3]
            c2w_new = np.empty((3, 4), dtype=c2w.dtype)
            c2w_new[:3, :3] = R_unit @ R_in
            c2w_new[:3, 3] = R_full @ t_in + t_full
            c2w = c2w_new
        camera = Cameras((fx, fy, cx, cy), [0.0]*10, h, w, torch.from_numpy(c2w), camera_model)
        camera = camera.to("cuda")
        # Pass want_keys=[buffer_key] so the model can skip the SH /
        # refinement_score debug renders + the depth_normal kernel when
        # the viewer didn't ask for them. Empty buffer_key ("") means
        # "render everything" -- used by the viewer's buffer-enumeration
        # path so the dropdown sees every supported key.
        want = None if not buffer_key else [buffer_key]
        outputs = self.model.get_outputs(camera, want_keys=want)

        # Stuff `None` placeholders into the output dict for every buffer the
        # viewer's dropdown should know about. This lets the buffer-enumeration
        # endpoint (and `last_keys` cache populated by the first hot render)
        # list every supported key even though `want_keys` only rendered the
        # currently selected one. The hot path stays fast because the actual
        # tensors are only produced for the requested key; when the user
        # switches to a previously-`None` slot, the next request renders it.
        for _k in ["rgb", "depth", "alpha"] + ["normal"] * 0 + ["depth_normal",
                   "depth_median", "normal_median",
                   "rgb_distortion", "depth_distortion",
                   "sh", "refinement_score"]:
            outputs.setdefault(_k, None)
        # Use the post-split Cameras when available (CPP DataManager path).
        # For the legacy Python DM path no warp split happens, so the
        # per-input cameras are already the right thing for the viewer.
        _viewer_cams = getattr(self, '_viewer_cameras', None)
        if _viewer_cams is None:
            _viewer_cams = self.dataset_train.cameras
        outputs['_post_processor'] = lambda tensor, **kwargs: annotate_train_cameras(
            buffer_key, tensor,
            outputs.get('depth', None), outputs['alpha'],
            camera, _viewer_cams,
            **kwargs
        )
        return outputs

    def render(self, c2w, fx, fy, cx, cy, w, h, camera_model,
               buffer_key="rgb", *, show_training_cameras: bool = False):
        # Invoking the post-processor + D->H copy inside this lock is what
        # makes viewer requests respond quickly during training: it ensures
        # the .cpu() Memcpy is queued on the default stream BEFORE the next
        # train_step queues iter N+1's kernels. Otherwise the Memcpy has to
        # wait for every queued training kernel, which can be ~1 s per iter
        # on large scenes -- the source of the multi-second viewer lag.
        # The "render desired" flag is driven by the RenderWorker's
        # on_submit/on_idle hooks -- the HTTP submit() flips it before we
        # ever get here, which is what lets train_step yield the lock at the
        # very next iteration boundary instead of after several missed ones.
        with self.lock:
            self.model.eval()
            outputs = self._render(c2w, fx, fy, cx, cy, w, h,
                                   camera_model, buffer_key)
            if buffer_key and outputs.get(buffer_key) is not None:
                pp = outputs.pop('_post_processor', None)
                if pp is not None:
                    annotated = pp(outputs[buffer_key],
                                   show_training_cameras=show_training_cameras)
                    outputs[buffer_key] = annotated.cpu().numpy()
            return outputs

    def _train_step(self, step: int):
        # ---- profiling probes (gated by PROFILE_TRAIN_STEP) ----
        if PROFILE_TRAIN_STEP:
            _prof = getattr(self, "_step_prof", None)
            if _prof is None:
                self._step_prof = _prof = {
                    "n": 0, "warmup": 10,
                    "step_cb": 0, "next_train": 0, "engine": 0,
                    "verbose": 0, "gpu_drain": 0,
                }
            from time import perf_counter_ns as _t
            t0 = _t()

        self.model.step_cb(step)
        if PROFILE_TRAIN_STEP: t1 = _t()

        # Managed path: engine pulls the next batch from its DataManager.
        # Python doesn't touch GT / camera tensors per step.
        if PROFILE_TRAIN_STEP: t2 = _t()
        self.model.engine_train_step_managed(step)
        # Optional: dump GT + render at specific steps for visual debugging.
        # Enable with e.g. `SS_DEBUG_DUMP_STEPS=0,1,10,100 SS_DEBUG_DUMP_DIR=debug ss-train ...`
        import os as _os
        _dump_steps = _os.environ.get("SS_DEBUG_DUMP_STEPS", "")
        if _dump_steps:
            try:
                _steps = {int(s) for s in _dump_steps.split(",") if s.strip()}
            except ValueError:
                _steps = set()
            if step in _steps:
                from spirulae_splat.modules import debug_image as _dbg
                _dir = _os.environ.get("SS_DEBUG_DUMP_DIR", "debug")
                prefix = f"{_dir}/step{step:06d}"
                arrs = _dbg.dump_engine_step(prefix)
                for k, a in arrs.items():
                    print(f"[debug_image] {prefix}_{k}: {_dbg.buffer_stats(a)}", flush=True)

        if PROFILE_TRAIN_STEP:
            t3 = _t()
            # GPU drain probe: how much GPU work is still queued after the
            # engine call returns. With option 3 (async D->H), the host gets
            # ahead and this becomes significant — meaning iter N+1's host
            # work overlaps with iter N's GPU tail.
            torch.cuda.synchronize()
            t3b = _t()

        self.model.verbose()

        if PROFILE_TRAIN_STEP:
            t4 = _t()
            if _prof["warmup"] > 0:
                _prof["warmup"] -= 1
            else:
                _prof["n"] += 1
                _prof["step_cb"]    += (t1 - t0)
                _prof["next_train"] += (t2 - t1)
                _prof["engine"]     += (t3 - t2)
                _prof["gpu_drain"]  += (t3b - t3)
                _prof["verbose"]    += (t4 - t3b)
            if _prof["n"] >= 100:
                n = _prof["n"]
                print(
                    f"\n[PROF n={n}] "
                    f"step_cb={_prof['step_cb']/n/1e6:.3f}ms "
                    f"next_train={_prof['next_train']/n/1e6:.3f}ms "
                    f"engine={_prof['engine']/n/1e6:.3f}ms "
                    f"gpu_drain={_prof['gpu_drain']/n/1e6:.3f}ms "
                    f"verbose={_prof['verbose']/n/1e6:.3f}ms "
                    f"sum={(t4-t0)/1e6:.3f}ms",
                    flush=True,
                )
                for k in ("step_cb", "next_train", "engine", "gpu_drain", "verbose"):
                    _prof[k] = 0
                _prof["n"] = 0

    def train_step(self, *args):
        # Yield to a pending render before contending for self.lock. A short
        # sleep is enough to let the render thread's blocked acquire win on
        # Linux (pthread_mutex wakes the longest-waiting thread when the
        # current owner has been idle); without it, the immediate re-acquire
        # at the top of the next iter beats the render thread to the lock.
        if self._render_pending.is_set():
            time.sleep(0.001)
        with self.lock:
            self.model.train()
            return self._train_step(*args)

    def _get_eval_metrics_dict(self):
        inputs = self.datamanager.next_train(0, None)  # type: List[Tuple[Cameras, Dict]]
        assert not isinstance(inputs, tuple)

        for i, (camera, batch) in enumerate(inputs):
            assert i == 0
            model_outputs = self.model.get_outputs(camera)
            metrics_dict, img_dict = self.model.get_image_metrics_and_images(model_outputs, batch)

        return metrics_dict, img_dict

    def get_eval_metrics_dict(self, *args):
        with self.lock:
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
        self.last_step_time = self.start_time
        for step in range(self._start_step, self.config.num_iterations):
            self._check_pause()  # Check if paused and wait if needed
            if step > self._start_step and self.config.steps_per_save > 0 and step % self.config.steps_per_save == 0:
                self.save_checkpoint(step)
            # if step % 100 == 50:
            #     self.print_vram_breakdown()
            #     exit(0)
            step_start = time.time()
            self.current_step = step + 1  # 1-based
            self.train_step(step)
            step_end = time.time()
            latency = step_end - step_start
            self.step_latencies.append(latency)
            if len(self.step_latencies) > 100:  # keep last 100
                self.step_latencies.pop(0)
            # pbar.update(1)
            avg_latency = sum(self.step_latencies) / len(self.step_latencies)
            elapsed = time.time() - self.start_time
            eta = (self.config.num_iterations - self.current_step) * avg_latency
        if self.config.steps_per_save != 0:
            self.save_checkpoint(self.config.num_iterations)
        self._training_time = time.perf_counter() - train_wall_start
        pool_cap = sum(cap for _, _, cap in _C.engine_get_pool_breakdown())
        self._engine_vram_bytes = pool_cap + _C.engine_get_scratch_bytes()
        print(f"Checkpoint saved to: {self.output_dir.absolute()}")
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
                    self.train_step(step)

    @torch.no_grad()
    def eval(self):
        if self.config.dataparser.eval_mode == "all" or len(self.dataset_eval.cameras) == 0:
            return

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
        from spirulae_splat.modules.resume import resolve_checkpoint, read_state_json
        from spirulae_splat.modules.resume_adapt import adapt_checkpoint

        run_dir, ckpt_dir = resolve_checkpoint(Path(load_dir))
        state = read_state_json(ckpt_dir)

        model = self.model
        cfg = self.config.model
        opt = self.config.optimizer
        n_img = model.num_train_data

        # --- target layout from the (possibly changed) config ---------------
        def _lhw(shape):                       # (grid_X, grid_Y, grid_W) -> (L,H,W)
            X, Y, W_g = shape
            return (W_g, Y, X)
        geo = cfg.use_bilateral_grid_for_geometry
        ppisp_lr = opt.ppisp_adagrad_lr if cfg.use_adagrad_ppisp_optim else opt.ppisp_lr
        target = {
            "max_num_splats": int(model.core.max_num_splats),
            "num_sh": int(model.gauss_params["features_sh"].shape[1]),
            "num_images": int(n_img),
            "bilagrid": {
                "rgb":    _lhw(cfg.bilagrid_shape) if (cfg.use_bilateral_grid and opt.bilagrid_lr > 0) else None,
                "depth":  _lhw(cfg.bilagrid_shape_geometry) if (geo and opt.bilagrid_depth_lr > 0 and cfg.depth_supervision_weight > 0) else None,
                "normal": _lhw(cfg.bilagrid_shape_geometry) if (geo and opt.bilagrid_normal_lr > 0 and cfg.normal_supervision_weight > 0) else None,
            },
            "ppisp": bool(cfg.use_ppisp and ppisp_lr > 0),
        }

        # --- eager-init appearance the checkpoint has AND the config wants ---
        def _ck_on(t):
            return bool(state.get(f"bilagrid_{t}", {}).get("enabled"))
        tbg = target["bilagrid"]
        if tbg["rgb"] is not None and _ck_on("rgb") and not model._bilagrid_rgb_init:
            model.core.engine_init_bilagrid(
                n_img, rgb_type=cfg.bilagrid_type, rgb_LHW=tbg["rgb"],
                use_adagrad=cfg.use_adagrad_bilagrid_optim)
            model._bilagrid_rgb_init = True
        if tbg["depth"] is not None and _ck_on("depth") and not model._bilagrid_depth_init:
            model.core.engine_init_bilagrid(
                n_img, depth_LHW=tbg["depth"], use_adagrad=cfg.use_adagrad_bilagrid_optim)
            model._bilagrid_depth_init = True
        if tbg["normal"] is not None and _ck_on("normal") and not model._bilagrid_normal_init:
            model.core.engine_init_bilagrid(
                n_img, normal_LHW=tbg["normal"], use_adagrad=cfg.use_adagrad_bilagrid_optim)
            model._bilagrid_normal_init = True
        if target["ppisp"] and bool(state.get("ppisp", {}).get("enabled")) and not model._ppisp_init:
            model.core.engine_init_ppisp(
                n_img, param_type=cfg.ppisp_param_type, use_adagrad=cfg.use_adagrad_ppisp_optim)
            model._ppisp_init = True

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
        load_depths=False,
        load_normals=False,  # TODO
        mask_boundary_offset=-0.025,
    ))
    model: SpirulaeSplatModelConfig = field(default_factory=lambda: SpirulaeSplatModelConfig(
        # densify_score_mode="median",
        densify_loss_map_mode="robust_edge_aware",
        densify_robust_edge_aware_quantile=0.99,
        long_axis_split_opacity_k=(0.5, 0.6, 5000),
    ))


@dataclass
class TrainerConfigInTheWild(TrainerConfig):
    """Preset for in-the-wild datasets, like datasets consisting of internet images, or datasets with extreme lighting variation and/or un-masked outliers."""
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
        ssim_lambda=0.1,
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
