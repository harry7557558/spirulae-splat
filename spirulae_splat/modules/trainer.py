from dataclasses import dataclass, field, asdict
from typing import Type, List, Tuple, Dict, Optional
from pathlib import Path

import threading
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



_DEFAULT_OPTIMIZERS = {
    "_dummy": FusedAdamOptimizerConfig(lr=1.0, eps=0.0),
    "means": FusedAdamOptimizerConfig(
        lr=1.6e-4, eps=1e-15,
        lr_final=1.6e-6, max_steps=30000
    ),
    "scales": FusedAdamOptimizerConfig(lr=0.005, eps=1e-15),
    "quats": FusedAdamOptimizerConfig(lr=0.0005, eps=1e-15),
    "features_dc": FusedAdamOptimizerConfig(lr=0.0025, eps=1e-15),
    "features_sh": FusedAdamOptimizerConfig(
        lr=0.0025 / 20, eps=1e-15,
        tr=1.0e-6 / 20, tr_final=1.0e-8 / 20, max_steps=30000,
    ),
    "features_ch": FusedAdamOptimizerConfig(lr=0.0025 / 5, eps=1e-15),
    "sv_sites": FusedAdamOptimizerConfig(lr=0.01, eps=1e-15),
    "sv_colors": FusedAdamOptimizerConfig(lr=0.0005, eps=1e-15),
    "opacities": FusedAdamOptimizerConfig(lr=0.05, eps=1e-15),
    "densities": FusedAdamOptimizerConfig(
        lr=0.05, eps=1e-15,
        lr_final=0.0005, max_steps=30000,
    ),
    "background_color": FusedAdamOptimizerConfig(lr=0.0025, eps=1e-15),
    "background_sh": FusedAdamOptimizerConfig(lr=0.0025 / 5, eps=1e-15),
    "bilateral_grid": FusedAdamOptimizerConfig(
        lr=2e-3, eps=1e-15,
        lr_final=1e-4, max_steps=30000, warmup_steps=1000,
    ),
    "bilateral_grid_depth": FusedAdamOptimizerConfig(
        lr=2e-3, eps=1e-15,
        lr_final=1e-4, max_steps=30000, warmup_steps=2000,
    ),
    "bilateral_grid_normal": FusedAdamOptimizerConfig(
        lr=5e-4, eps=1e-15,
        lr_final=4e-5, max_steps=30000, warmup_steps=2000,
    ),
    "ppisp": FusedAdamOptimizerConfig(
        lr=2e-3, eps=1e-15,
        lr_final=2e-5, max_steps=30000, warmup_steps=500, lr_pre_warmup=2e-5
    ),
    "camera_opt": FusedAdamOptimizerConfig(
        lr=1e-4, eps=1e-15,
        lr_final=5e-7, max_steps=30000, warmup_steps=1000,
    ),
}

_DEFAULT_OPTIMIZERS_WITH_SCALE_SCHEDULER = {**_DEFAULT_OPTIMIZERS}
_DEFAULT_OPTIMIZERS_WITH_SCALE_SCHEDULER["scales"] = FusedAdamOptimizerConfig(
    # https://arxiv.org/abs/2603.08661
    lr=0.02, eps=1e-15,
    lr_final=0.005, max_steps=10000, warmup_steps=1000, lr_pre_warmup=0.005
)

_TRIANGLE_OPTIMIZERS = {**_DEFAULT_OPTIMIZERS_WITH_SCALE_SCHEDULER}
_TRIANGLE_OPTIMIZERS["means"] = FusedAdamOptimizerConfig(
    lr=1.6e-4, eps=1e-15, lr_final=1.6e-6, max_steps=30000,
)
# _TRIANGLE_OPTIMIZERS["scales"] = FusedAdamOptimizerConfig(
#     lr=0.005, eps=1e-15, lr_final=0.0002, max_steps=30000,
# )
# _TRIANGLE_OPTIMIZERS["quats"] = FusedAdamOptimizerConfig(
#     lr=0.0005, eps=1e-15, lr_final=0.0001, max_steps=30000,
# )
_TRIANGLE_OPTIMIZERS["bilateral_grid"] = FusedAdamOptimizerConfig(
    lr=5e-4, eps=1e-15,
    lr_final=1e-6, max_steps=30000, warmup_steps=1000,
)

_SECOND_ORDER_POSITION_OPTIMIZERS = {**_DEFAULT_OPTIMIZERS_WITH_SCALE_SCHEDULER}
_SECOND_ORDER_POSITION_OPTIMIZERS["means"] = FusedNewtonOptimizerConfig(
    mode="mean", lr=1.0e-6, eps=1e-15,
    lr_final=1.0e-8, max_steps=30000, #warmup_steps=1000,
)

_SECOND_ORDER_OPTIMIZERS = {**_SECOND_ORDER_POSITION_OPTIMIZERS}
# _SECOND_ORDER_OPTIMIZERS["quats"] = FusedNewtonOptimizerConfig(
#     mode="quat", lr=1.0e-6, eps=1e-15,
#     lr_final=1.0e-8, max_steps=30000, #warmup_steps=1000,
# )
_SECOND_ORDER_OPTIMIZERS["scales"] = FusedNewtonOptimizerConfig(
    mode="scale", lr=1.0e-6, eps=1e-15,
    lr_final=1.0e-8, max_steps=30000, #warmup_steps=1000,
    # mode="scale", lr=1.0e-5, eps=1e-15,
    # lr_final=1.0e-7, max_steps=30000, #warmup_steps=1000,
)
# TODO: investigate whether this messes up MCMC densification
# _SECOND_ORDER_OPTIMIZERS["opacities"] = FusedNewtonOptimizerConfig(
#     mode="opacity", lr=1.0e-6, eps=1e-15,
#     lr_final=1.0e-8, max_steps=30000, #warmup_steps=3000,
# )
# _SECOND_ORDER_OPTIMIZERS["features_dc"] = FusedNewtonOptimizerConfig(
#     mode="color", lr=1.0e-6, eps=1e-15,
#     lr_final=1.0e-8, max_steps=30000, #warmup_steps=1000,
# )


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
    save_eval_images: bool = False
    """Whether to save eval images at end of training"""

    num_iterations: int = 30000
    """Number of training iterations"""

    viewer_port: int = 7007
    """Port used by the web viewer"""

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
        self.config = config
        self.dataparser = SpirulaeSplatDataparser(config.dataparser, config.data)
        self.dataparser_outputs_train, self.dataparser_outputs_eval = self.dataparser.parse()
        self.dataset_train = SpirulaeSplatDataset(self.dataparser_outputs_train)
        self.dataset_eval = SpirulaeSplatDataset(self.dataparser_outputs_eval)

        self.datamanager = SpirulaeSplatDataManager(self.config.datamanager, device="cuda")
        self.datamanager.train_dataset = self.dataset_train

        self.model = SpirulaeSplatModel(
            self.config,
            self.dataparser_outputs_train['metadata'],
            self.dataparser_outputs_train['cameras']
        ).cuda()
        # self.optimizers = create_optimizers(self.model, self.config.optimizer)

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
            'scale': self.dataparser_outputs_train['dataparser_scale']
        }
        with open(self.output_dir / "dataparser_transforms.json", "w") as f:
            json.dump(dataparser_transform_dict, f, indent=4)

    def _render(self, c2w, fx, fy, cx, cy, w, h, camera_model):
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
        # for optim in self.optimizers.values():
        #     optim.zero_grad()
        # self.model.step_cb(self.optimizers, step)
        self.model.step_cb(step)

        inputs = self.datamanager.next_train(step, self.config.num_iterations)  # type: List[Tuple[Cameras, Dict]]
        if isinstance(inputs, tuple):
            train_inputs, val_inputs = inputs
            if len(train_inputs) > 1 or len(val_inputs) > 1:
                raise NotImplementedError("Validation with split_batch is not supported")  # TODO
            train_inputs, val_inputs = train_inputs[0], val_inputs[0]
            inputs = [((train_inputs[0], val_inputs[0]), (train_inputs[1], val_inputs[1]))]

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

        self.model.verbose()

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
        config.compute_visibility_masks = False
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

        import json
        with open(self.output_dir / "metrics.json", "w") as f:
            json.dump(metrics, f, indent=4)

    def save_checkpoint(self, step: int) -> None:
        ckpt_path: Path = self.output_dir / f"step-{step:09d}.ckpt"
        torch.save(
            {
                "step": step,
                "model": self.model.state_dict(),
                # "optimizers": {k: v.state_dict() for (k, v) in self.optimizers.items()},
            },
            ckpt_path,
        )
        # delete previous checkpoint
        if self.config.save_only_latest_checkpoint:
            for f in self.output_dir.glob("*.ckpt"):
                if f != ckpt_path:
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
        patch_batch_size=-1,
        patch_size=64,
        max_batch_per_epoch=800,
    ))
    model: SpirulaeSplatModelConfig = field(default_factory=lambda: SpirulaeSplatModelConfig(
        # packed=True,
        use_bvh=True,
        use_camera_optimizer=False,
        use_bilateral_grid=False,
        use_bilateral_grid_for_geometry=False,  # TODO: slow
        alpha_reg_weight=0.0,
        primitive="mip", max_screen_size=float('inf'),  # TODO
    ))
    # optimizer: dict = field(default_factory=lambda: _DEFAULT_OPTIMIZERS_WITH_SCALE_SCHEDULER)  # TODO

@dataclass
class TrainerConfigTriangle(TrainerConfig):
    """Method for triangle splatting"""
    datamanager: SpirulaeSplatDataManagerConfig = field(default_factory=lambda: SpirulaeSplatDataManagerConfig(
        compute_visibility_masks=True,
    ))
    model: SpirulaeSplatModelConfig = field(default_factory=lambda: SpirulaeSplatModelConfig(
        primitive="opaque_triangle",
        kernel_radius=0.5,
        sh_degree=0,
        bilagrid_shape=(16, 16, 8),
        refine_stop_num_iter=0,
        # alpha_reg_weight=0.0,
        scale_reg=0.001,
        # erank_reg=1.0,
        # supersampling=2,
        min_opacity=0.005,
        supervision_warmup=0,
        depth_supervision_weight=0.0,
        normal_supervision_weight=0.04,
        depth_distortion_reg=0.01,
        normal_distortion_reg=0.005,
        rgb_distortion_reg=0.01,
        ssim_lambda=0.4,
        preallocate_splat_tensors=False,  # TODO
    ))
    # optimizer: dict = field(default_factory=lambda: _TRIANGLE_OPTIMIZERS)  # TODO

@dataclass
class TrainerConfigTrianglePatched(TrainerConfig):
    """Method for triangle splatting with patched batching"""
    datamanager: SpirulaeSplatDataManagerConfig = field(default_factory=lambda: SpirulaeSplatDataManagerConfig(
        compute_visibility_masks=True,
        patch_batch_size=-1,
        patch_size=64,
        max_batch_per_epoch=800,
    ))
    model: SpirulaeSplatModelConfig = field(default_factory=lambda: SpirulaeSplatModelConfig(
        primitive="opaque_triangle",
        kernel_radius=0.5,
        sh_degree=0,
        bilagrid_shape=(16, 16, 8),
        refine_stop_num_iter=0,
        # alpha_reg_weight=0.0,
        scale_reg=0.001,
        # erank_reg=1.0,
        # supersampling=2,
        min_opacity=0.005,
        supervision_warmup=0,
        depth_supervision_weight=0.0,
        normal_supervision_weight=0.04,
        depth_distortion_reg=0.01,
        normal_distortion_reg=0.005,
        rgb_distortion_reg=0.01,
        ssim_lambda=0.4,
        reg_warmup_length=1000,
        preallocate_splat_tensors=False,  # TODO

        # packed=True,
        use_bvh=True,
        use_camera_optimizer=False,
        # use_bilateral_grid=False,
        # use_bilateral_grid_for_geometry=False,
        alpha_reg_weight=0.0,
        max_screen_size_clip_hardness=1.5,
    ))
    # optimizer: dict = field(default_factory=lambda: _TRIANGLE_OPTIMIZERS)  # TODO

@dataclass
class TrainerConfigVoxel(TrainerConfig):
    """[WIP] Does not work at this time, do not use"""
    # ^ This message will appear in `spirulae-train --help`
    model: SpirulaeSplatModelConfig = field(default_factory=lambda: SpirulaeSplatModelConfig(
        primitive="voxel",
        use_mcmc=False,
        sh_degree=0,
        background_color="black",
        train_background_color=False,
        # use_bilateral_grid=False,
        use_bilateral_grid_for_geometry=False,
        alpha_reg_weight=0.0,
        # alpha_loss_weight=0.0,
        # alpha_loss_weight_under=0.0,
        opacity_reg=0.0001,
        scale_reg=0.0,
        erank_reg=0.0,
        erank_reg_s3=0.0,
        depth_supervision_weight=0.0,
        supervision_warmup=0,
        depth_distortion_reg=0.01,
        reg_warmup_length=100,
        distortion_reg_warmup=100,
        preallocate_splat_tensors=False,  # TODO
    ))
    # optimizer: dict = field(default_factory=lambda: _DEFAULT_OPTIMIZERS)  # TODO


_MODEL_PRESET_CONFINED = dict(
    randomize_background=True,
    train_background_color=False,
)
_MODEL_PRESET_OPEN = dict(
    randomize_background=False,
    train_background_color=True,
    background_color="gray",
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
    relocate_heuristic_weight=0.0,
    use_long_axis_split=False,
)
_MODEL_PRESET_RICH_TEXTURE = dict(
    use_edge_aware_score=True,
    relocate_heuristic_weight=1.0,
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
        start_resolution=None,
        load_depths=False,
        load_normals=False,
    ))
    model: SpirulaeSplatModelConfig = field(default_factory=lambda: SpirulaeSplatModelConfig(
        primitive="3dgs",
        relative_scale=1.0,
        background_color="black",
        train_background_color=False,
        use_bilateral_grid=False,
        use_bilateral_grid_for_geometry=False,
        use_ppisp=False,
        relocate_heuristic_weight=0.0,
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
        randomize_background=0.0,
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
