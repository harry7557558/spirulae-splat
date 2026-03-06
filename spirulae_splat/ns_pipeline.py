"""
Nerfstudio Template Pipeline
"""

import typing
from dataclasses import dataclass, field
from typing import List, Literal, Optional, Tuple, Type, Dict

import torch
import torch.distributed as dist
from torch.cuda.amp.grad_scaler import GradScaler
from torch.nn.parallel import DistributedDataParallel as DDP

from spirulae_splat.ns_datamanager import SpirulaeDataManagerConfig
from spirulae_splat.ns_model import SpirulaeModel, SpirulaeModelConfig
from nerfstudio.data.datamanagers.base_datamanager import (
    DataManager,
    DataManagerConfig,
)
from nerfstudio.models.base_model import ModelConfig
from nerfstudio.pipelines.base_pipeline import (
    VanillaPipeline,
    VanillaPipelineConfig,
)
from nerfstudio.utils import profiler
from nerfstudio.cameras.cameras import Cameras, CameraType


@dataclass
class SpirulaePipelineConfig(VanillaPipelineConfig):
    """Configuration for pipeline instantiation"""

    _target: Type = field(default_factory=lambda: SpirulaePipeline)
    """target class to instantiate"""
    datamanager: DataManagerConfig = field(default_factory=SpirulaeDataManagerConfig)
    """specifies the datamanager config"""
    model: ModelConfig = field(default_factory=SpirulaeModelConfig)
    """specifies the model config"""


class SpirulaePipeline(VanillaPipeline):
    """Template Pipeline

    Args:
        config: the pipeline config used to instantiate class
    """

    def __init__(
        self,
        config: SpirulaePipelineConfig,
        device: str,
        test_mode: Literal["test", "val", "inference"] = "val",
        world_size: int = 1,
        local_rank: int = 0,
        grad_scaler: Optional[GradScaler] = None,
    ):
        super(VanillaPipeline, self).__init__()
        self.config = config
        self.test_mode = test_mode
        self.datamanager: DataManager = config.datamanager.setup(
            device=device, test_mode=test_mode, world_size=world_size, local_rank=local_rank
        )

        seed_pts = None
        if (
            hasattr(self.datamanager, "train_dataparser_outputs")
            and "points3D_xyz" in self.datamanager.train_dataparser_outputs.metadata
        ):
            pts = self.datamanager.train_dataparser_outputs.metadata["points3D_xyz"]
            pts_rgb = self.datamanager.train_dataparser_outputs.metadata["points3D_rgb"]
            seed_pts = (pts, pts_rgb)

        assert self.datamanager.train_dataset is not None, "Missing input dataset"
        self._model = config.model.setup(
            scene_box=self.datamanager.train_dataset.scene_box,
            num_train_data=len(self.datamanager.train_dataset),
            metadata=self.datamanager.train_dataset.metadata,
            device=device,
            grad_scaler=grad_scaler,
            seed_points=seed_pts,
        )
        self.model.to(device)

        self.world_size = world_size
        if world_size > 1:
            self._model = typing.cast(
                SpirulaeModel, DDP(self._model, device_ids=[local_rank], find_unused_parameters=True)
            )
            dist.barrier(device_ids=[local_rank])

    @profiler.time_function
    def get_train_loss_dict(self, step: int):
        """This function gets your training loss dict. This will be responsible for
        getting the next batch of data from the DataManager and interfacing with the
        Model class, feeding the data to the model's forward function.

        Args:
            step: current iteration step to update sampler if using DDP (distributed)
        """
        inputs = self.datamanager.next_train(step)  # type: List[Tuple[Cameras, Dict]]
        if isinstance(inputs, tuple):
            train_inputs, val_inputs = inputs
            if len(train_inputs) > 1 or len(val_inputs) > 1:
                raise NotImplementedError("Validation with split_batch is not supported")  # TODO
            train_inputs, val_inputs = train_inputs[0], val_inputs[0]
            inputs = [((train_inputs[0], val_inputs[0]), (train_inputs[1], val_inputs[1]))]

        for i, (camera, batch) in enumerate(inputs):
            torch.cuda.empty_cache()
            model_outputs = self._model(camera)
            torch.cuda.empty_cache()
            metrics_dict = self.model.get_metrics_dict(model_outputs, batch)
            is_not_last = (i != len(inputs) - 1)
            loss_dict = self.model.get_loss_dict(model_outputs, batch, metrics_dict, no_static_losses=is_not_last)
            torch.cuda.empty_cache()
            if is_not_last:
                torch.stack([
                    x for x in loss_dict.values() if isinstance(x, torch.Tensor)
                ]).sum().backward()
            torch.cuda.empty_cache()

        return model_outputs, loss_dict, metrics_dict

    def get_average_eval_image_metrics(
        self, step=None, output_path=None, get_std=False
    ):
        result = super().get_average_eval_image_metrics(
            step, output_path, get_std)
        def print_metric(key, get_std, decimals, name=None, sc=1.0):
            if key not in result:
                return
            if name is None:
                name = key
            name += ' '*(8-len(name))
            fmt = "{:."+str(decimals)+"f}"
            val = fmt.format(sc*result[key])
            if get_std:
                std = fmt.format(sc*result[key+'_std'])
                print(name, val, '±', std)
            else:
                print(name, val)
        print()
        print_metric('gaussian_count', False, 0, "splats")
        print('-'*8)
        print_metric('psnr', get_std, 2)
        print_metric('ssim', get_std, 3)
        print_metric('lpips', get_std, 3)
        print('-'*8)
        print_metric('num_rays_per_sec', get_std, 2, "Mrays/s", 1e-6)
        print_metric('fps', get_std, 1)
        print()
        return result
