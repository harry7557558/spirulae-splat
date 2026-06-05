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
from spirulae_splat.modules.optimizer import create_optimizers, FusedAdamOptimizerConfig, FusedNewtonOptimizerConfig, OptimizerConfig

from spirulae_splat.viewer.annotation import annotate_train_cameras
from spirulae_splat.modules._profile import PROFILE_TRAIN_STEP



# _DEFAULT_OPTIMIZERS = {
#     "_dummy": FusedAdamOptimizerConfig(lr=1.0, eps=0.0),
#     "means": FusedAdamOptimizerConfig(
#         lr=1.6e-4, eps=1e-15,
#         lr_final=1.6e-6, max_steps=30000
#     ),
#     "scales": FusedAdamOptimizerConfig(lr=0.005, eps=1e-15),
#     "quats": FusedAdamOptimizerConfig(lr=0.0005, eps=1e-15),
#     "features_dc": FusedAdamOptimizerConfig(lr=0.0025, eps=1e-15),
#     "features_sh": FusedAdamOptimizerConfig(
#         lr=0.0025 / 20, eps=1e-15,
#         tr=1.0e-6 / 20, tr_final=1.0e-8 / 20, max_steps=30000,
#     ),
#     "opacities": FusedAdamOptimizerConfig(lr=0.05, eps=1e-15),
#     "background_dc": FusedAdamOptimizerConfig(lr=0.0025, eps=1e-15),
#     "background_sh": FusedAdamOptimizerConfig(lr=0.0025 / 20, eps=1e-15),
#     "bilateral_grid": FusedAdamOptimizerConfig(
#         lr=2e-3, eps=1e-15,
#         lr_final=1e-4, max_steps=30000, warmup_steps=1000,
#     ),
#     "bilateral_grid_depth": FusedAdamOptimizerConfig(
#         lr=2e-3, eps=1e-15,
#         lr_final=1e-4, max_steps=30000, warmup_steps=2000,
#     ),
#     "bilateral_grid_normal": FusedAdamOptimizerConfig(
#         lr=5e-4, eps=1e-15,
#         lr_final=4e-5, max_steps=30000, warmup_steps=2000,
#     ),
#     "ppisp": FusedAdamOptimizerConfig(
#         lr=2e-3, eps=1e-15,
#         lr_final=2e-5, max_steps=30000, warmup_steps=500, lr_pre_warmup=2e-5
#     ),
#     "camera_opt": FusedAdamOptimizerConfig(
#         lr=1e-4, eps=1e-15,
#         lr_final=5e-7, max_steps=30000, warmup_steps=1000,
#     ),
# }

# _DEFAULT_OPTIMIZERS_WITH_SCALE_SCHEDULER = {**_DEFAULT_OPTIMIZERS}
# _DEFAULT_OPTIMIZERS_WITH_SCALE_SCHEDULER["scales"] = FusedAdamOptimizerConfig(
#     # https://arxiv.org/abs/2603.08661
#     lr=0.02, eps=1e-15,
#     lr_final=0.005, max_steps=10000, warmup_steps=1000, lr_pre_warmup=0.005
# )

# _TRIANGLE_OPTIMIZERS = {**_DEFAULT_OPTIMIZERS_WITH_SCALE_SCHEDULER}
# _TRIANGLE_OPTIMIZERS["means"] = FusedAdamOptimizerConfig(
#     lr=1.6e-4, eps=1e-15, lr_final=1.6e-6, max_steps=30000,
# )
# # _TRIANGLE_OPTIMIZERS["scales"] = FusedAdamOptimizerConfig(
# #     lr=0.005, eps=1e-15, lr_final=0.0002, max_steps=30000,
# # )
# # _TRIANGLE_OPTIMIZERS["quats"] = FusedAdamOptimizerConfig(
# #     lr=0.0005, eps=1e-15, lr_final=0.0001, max_steps=30000,
# # )
# _TRIANGLE_OPTIMIZERS["bilateral_grid"] = FusedAdamOptimizerConfig(
#     lr=5e-4, eps=1e-15,
#     lr_final=1e-6, max_steps=30000, warmup_steps=1000,
# )

# _SECOND_ORDER_POSITION_OPTIMIZERS = {**_DEFAULT_OPTIMIZERS_WITH_SCALE_SCHEDULER}
# _SECOND_ORDER_POSITION_OPTIMIZERS["means"] = FusedNewtonOptimizerConfig(
#     mode="mean", lr=1.0e-6, eps=1e-15,
#     lr_final=1.0e-8, max_steps=30000, #warmup_steps=1000,
# )

# _SECOND_ORDER_OPTIMIZERS = {**_SECOND_ORDER_POSITION_OPTIMIZERS}
# # _SECOND_ORDER_OPTIMIZERS["quats"] = FusedNewtonOptimizerConfig(
# #     mode="quat", lr=1.0e-6, eps=1e-15,
# #     lr_final=1.0e-8, max_steps=30000, #warmup_steps=1000,
# # )
# _SECOND_ORDER_OPTIMIZERS["scales"] = FusedNewtonOptimizerConfig(
#     mode="scale", lr=1.0e-6, eps=1e-15,
#     lr_final=1.0e-8, max_steps=30000, #warmup_steps=1000,
#     # mode="scale", lr=1.0e-5, eps=1e-15,
#     # lr_final=1.0e-7, max_steps=30000, #warmup_steps=1000,
# )
# # TODO: investigate whether this messes up MCMC densification
# # _SECOND_ORDER_OPTIMIZERS["opacities"] = FusedNewtonOptimizerConfig(
# #     mode="opacity", lr=1.0e-6, eps=1e-15,
# #     lr_final=1.0e-8, max_steps=30000, #warmup_steps=3000,
# # )
# # _SECOND_ORDER_OPTIMIZERS["features_dc"] = FusedNewtonOptimizerConfig(
# #     mode="color", lr=1.0e-6, eps=1e-15,
# #     lr_final=1.0e-8, max_steps=30000, #warmup_steps=1000,
# # )


@dataclass
class TrainerConfig:
    """Default 3DGS method"""

    data: Path
    """Path to dataset. Can be a Nerfstudio or a COLMAP dataset."""

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
    """If True, each checkpoint additionally dumps every world/bilagrid/PPISP
        parameter and Adam optimizer state as one .npy per buffer in a `full/`
        subfolder. Useful for offline inspection of training dynamics."""
    save_eval_images: bool = False
    """Whether to save eval images at end of training"""

    num_iterations: int = 30000
    """Number of training iterations"""

    viewer_port: int = 7007
    """Port used by the web viewer"""

    disable_viewer: bool = False
    """If True, ss_trainer skips starting the viewer thread. Used by
        ss_benchmark so each scene runs without competing for the viewer port."""

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
        # self.optimizers = create_optimizers(self.model, self.config.optimizer)

        # When the C++ DataManager is enabled, install it on the engine once
        # here. Per-step training will then go through
        # ``model.engine_train_step_managed`` and bypass the Python data path.
        self._use_cpp_dm = bool(self.config.datamanager.use_cpp_data_manager)
        if self._use_cpp_dm:
            self._setup_cpp_data_manager()

        self.output_dir = self._setup_output_dir()
        print(f"Output directory: {self.output_dir.absolute()}")

        self._save_config_json()

        self.lock = threading.Lock()

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
        if dm_cfg.split_batch:
            raise NotImplementedError("use_cpp_data_manager + split_batch")
        if getattr(model_cfg, 'use_camera_optimizer', False):
            raise NotImplementedError("use_cpp_data_manager + use_camera_optimizer")
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

        # Warp path cannot supply depth / normal / PPISP supervision.
        if any_warp:
            if dm_cfg.load_depths and dpo['metadata'].get('depth_filenames', None):
                raise NotImplementedError(
                    "use_cpp_data_manager warp (fisheye/equirectangular split) "
                    "does not support depth supervision yet.")
            if dm_cfg.load_normals and dpo['metadata'].get('normal_filenames', None):
                raise NotImplementedError(
                    "use_cpp_data_manager warp does not support normal supervision yet.")
            if getattr(model_cfg, 'use_ppisp', False):
                raise NotImplementedError(
                    "use_cpp_data_manager warp does not support PPISP yet.")

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
        for i in range(N):
            k_i = K_per_camera[i]
            off = post_offsets[i]
            ci  = c2w_all[i]  # [3, 4]
            if k_i == 1:
                post_c2w[off]     = ci
                post_intrins[off] = input_intrins[i]
                post_dist[off]    = input_dist_coeffs[i]
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
        c_cfg.load_masks       = (len(mask_filenames) > 0)
        c_cfg.load_depths      = (len(depth_filenames) > 0 and not any_warp)
        c_cfg.load_normals     = (len(normal_filenames) > 0 and not any_warp)
        c_cfg.train_batch_size = train_bs
        c_cfg.val_batch_size   = val_bs
        c_cfg.warp_to_pinhole  = bool(dm_cfg.warp_to_pinhole)

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

        # When the warp path is active, the bilagrid table is sized to the
        # POST-split camera count -- override the model's num_train_data
        # accordingly (the model's own pre-init guess uses a uniform K).
        if any_warp:
            self.model.num_train_data = int(n_post)

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

    def _render(self, c2w, fx, fy, cx, cy, w, h, camera_model):
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
        outputs = self.model.get_outputs(camera)
        outputs['_post_processor'] = lambda tensor, **kwargs: annotate_train_cameras(
            tensor, outputs['depth'], outputs['alpha'],
            camera, self.dataset_train.cameras, self.dataset_train.thumbnails,
            relative_scale=self.model.config.relative_scale,
            warp_to_pinhole=self.datamanager.config.warp_to_pinhole, **kwargs
        )
        return outputs

    def render(self, *args):
        with self.lock:
            self.model.eval()
            return self._render(*args)

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

        # for optim in self.optimizers.values():
        #     optim.zero_grad()
        # self.model.step_cb(self.optimizers, step)
        self.model.step_cb(step)
        if PROFILE_TRAIN_STEP: t1 = _t()

        if self._use_cpp_dm:
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
        else:
            inputs = self.datamanager.next_train(step, self.config.num_iterations)  # type: List[Tuple[Cameras, Dict]]
            if isinstance(inputs, tuple):
                train_inputs, val_inputs = inputs
                if len(train_inputs) > 1 or len(val_inputs) > 1:
                    raise NotImplementedError("Validation with split_batch is not supported")  # TODO
                train_inputs, val_inputs = train_inputs[0], val_inputs[0]
                inputs = [((train_inputs[0], val_inputs[0]), (train_inputs[1], val_inputs[1]))]
            if PROFILE_TRAIN_STEP: t2 = _t()

            if len(inputs) != 1:
                raise NotImplementedError()
            for i, (camera, batch) in enumerate(inputs):
                # Engine path: fused C++ training step (no render output H↔D copy)
                if (self.model.core.primitive in ['3dgs', 'mip', '3dgut']
                        and not self.model.core.use_bvh
                        and not isinstance(camera, tuple)):
                    self.model.engine_train_step(camera, batch)
                else:
                    model_outputs = self.model.get_outputs(camera)
                    loss_grad = self.model.get_loss_grad(model_outputs, batch, len(inputs))
                    self.model.backward(model_outputs, loss_grad)
                    self.model.optim_step()

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

        # for optim in self.optimizers.values():
        #     optim.step()
        # self.model.step_post_backward()

    def train_step(self, *args):
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
        for step in range(self.config.num_iterations):
            self._check_pause()  # Check if paused and wait if needed
            if step > 0 and self.config.steps_per_save > 0 and step % self.config.steps_per_save == 0:
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
                        shutil.rmtree(f)
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
class TrainerConfigSquaredPos(TrainerConfig):
    """Method with second-order optimizer for positions"""
    model: SpirulaeSplatModelConfig = field(default_factory=lambda: SpirulaeSplatModelConfig(
        compute_hessian_diagonal="position",
    ))
    # optimizer: dict = field(default_factory=lambda: _SECOND_ORDER_POSITION_OPTIMIZERS)  # TODO


@dataclass
class TrainerConfigSquared(TrainerConfig):
    """Method with second-order optimizer"""
    model: SpirulaeSplatModelConfig = field(default_factory=lambda: SpirulaeSplatModelConfig(
        compute_hessian_diagonal="all",
    ))
    # optimizer: dict = field(default_factory=lambda: _SECOND_ORDER_OPTIMIZERS)  # TODO


@dataclass
class TrainerConfigPatched(TrainerConfig):
    """Method with patched batching"""
    datamanager: SpirulaeSplatDataManagerConfig = field(default_factory=lambda: SpirulaeSplatDataManagerConfig(
        # patch_batch_size=-1,
        # patch_size=64,
        max_batch_per_epoch=800,
    ))
    model: SpirulaeSplatModelConfig = field(default_factory=lambda: SpirulaeSplatModelConfig(
        packed=True,
        use_bvh=True,
        use_camera_optimizer=False,
        use_bilateral_grid=False,
        use_bilateral_grid_for_geometry=False,  # TODO: slow
        primitive="mip", max_screen_size=float('inf'),  # TODO
    ))
    # optimizer: dict = field(default_factory=lambda: _DEFAULT_OPTIMIZERS_WITH_SCALE_SCHEDULER)  # TODO


_MODEL_PRESET_CONFINED = dict(
    background_mode="noise",
)
_MODEL_PRESET_OPEN = dict(
    background_mode="sh",
    background_sh_degree=4,
)
_MODEL_PRESET_3DGS2TR_POS = dict(
    compute_hessian_diagonal="position",
)
_MODEL_PRESET_3DGS2TR = dict(
    compute_hessian_diagonal="all",
)
_MODEL_PRESET_LOW_TEXTURE = dict(
    use_edge_aware_score=False,
    use_revised_densification=False,
    use_long_axis_split=False,
)
_MODEL_PRESET_RICH_TEXTURE = dict(
    use_edge_aware_score=True,
    use_revised_densification=True,
    use_long_axis_split=True,
    max_screen_size=0.2,  # default 0.3
)
_MODEL_PRESET_NO_COLOR_SHIFT = dict(
    use_bilateral_grid=True,
    bilagrid_shape=(8, 8, 4),
    bilagrid_type="loglinear",
)


@dataclass
class TrainerConfigConfinedLowTexture(TrainerConfig):
    """Preset for visually confined environments with large textureless surfaces."""
    model: SpirulaeSplatModelConfig = field(default_factory=lambda: SpirulaeSplatModelConfig(
        **_MODEL_PRESET_LOW_TEXTURE,
        **_MODEL_PRESET_CONFINED,
        **_MODEL_PRESET_NO_COLOR_SHIFT,
    ))
    # optimizer: dict = field(default_factory=lambda: _DEFAULT_OPTIMIZERS_WITH_SCALE_SCHEDULER)  # TODO

@dataclass
class TrainerConfigConfined(TrainerConfig):
    """Preset for visually confined environments with moderate to rich texture."""
    model: SpirulaeSplatModelConfig = field(default_factory=lambda: SpirulaeSplatModelConfig(
        **_MODEL_PRESET_CONFINED,
        **_MODEL_PRESET_3DGS2TR_POS,
        **_MODEL_PRESET_RICH_TEXTURE,
    ))
    # optimizer: dict = field(default_factory=lambda: _SECOND_ORDER_POSITION_OPTIMIZERS)  # TODO

@dataclass
class TrainerConfigConfinedSquared(TrainerConfig):
    """[Unstable] Preset for visually confined environments with moderate to rich texture."""
    model: SpirulaeSplatModelConfig = field(default_factory=lambda: SpirulaeSplatModelConfig(
        **_MODEL_PRESET_CONFINED,
        **_MODEL_PRESET_3DGS2TR,
        **_MODEL_PRESET_RICH_TEXTURE,
    ))
    # optimizer: dict = field(default_factory=lambda: _SECOND_ORDER_OPTIMIZERS)  # TODO

@dataclass
class TrainerConfigOpenLowTexture(TrainerConfig):
    """Preset for open environments large textureless surfaces, a sky box will be trained."""
    model: SpirulaeSplatModelConfig = field(default_factory=lambda: SpirulaeSplatModelConfig(
        **_MODEL_PRESET_LOW_TEXTURE,
        **_MODEL_PRESET_OPEN,
        **_MODEL_PRESET_NO_COLOR_SHIFT,
        # relative_scale=10.0,
    ))
    # optimizer: dict = field(default_factory=lambda: _DEFAULT_OPTIMIZERS_WITH_SCALE_SCHEDULER)  # TODO

@dataclass
class TrainerConfigOpen(TrainerConfig):
    """Preset for open environments, a sky box will be trained."""
    model: SpirulaeSplatModelConfig = field(default_factory=lambda: SpirulaeSplatModelConfig(
        **_MODEL_PRESET_OPEN,
        **_MODEL_PRESET_3DGS2TR_POS,
        **_MODEL_PRESET_RICH_TEXTURE,
    ))
    # optimizer: dict = field(default_factory=lambda: _SECOND_ORDER_POSITION_OPTIMIZERS)  # TODO

@dataclass
class TrainerConfigOpenSquared(TrainerConfig):
    """[Unstable] Preset for open environments, a sky box will be trained."""
    model: SpirulaeSplatModelConfig = field(default_factory=lambda: SpirulaeSplatModelConfig(
        **_MODEL_PRESET_OPEN,
        **_MODEL_PRESET_3DGS2TR,
        **_MODEL_PRESET_RICH_TEXTURE,
    ))
    # optimizer: dict = field(default_factory=lambda: _SECOND_ORDER_OPTIMIZERS)  # TODO

@dataclass
class TrainerConfigCenteredObject(TrainerConfig):
    """Preset for centered objects, can also be used for human avatar."""
    dataparser: SpirulaeSplatDataParserConfig = field(default_factory=lambda: SpirulaeSplatDataParserConfig(
        center_method="focus",
    ))
    model: SpirulaeSplatModelConfig = field(default_factory=lambda: SpirulaeSplatModelConfig(
        **_MODEL_PRESET_OPEN,
        **_MODEL_PRESET_3DGS2TR,
        **_MODEL_PRESET_RICH_TEXTURE,
        **_MODEL_PRESET_NO_COLOR_SHIFT,
        apply_loss_for_mask=True,
        cap_max=100000,
    ))
    # optimizer: dict = field(default_factory=lambda: _DEFAULT_OPTIMIZERS_WITH_SCALE_SCHEDULER)
    # optimizer: dict = field(default_factory=lambda: _SECOND_ORDER_POSITION_OPTIMIZERS)  # TODO

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
        use_edge_aware_score=False,
        use_loss_map=False,
        use_long_axis_split=False,
        use_fused_proj_bwd_optim=False,
        quantize_sh_optim=False,
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
