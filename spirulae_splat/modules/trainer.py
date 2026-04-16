from dataclasses import dataclass, field
from typing import Type, List, Tuple, Dict
from pathlib import Path

import threading
import torch

from spirulae_splat.modules.camera import Cameras
from spirulae_splat.modules.model import SpirulaeSplatModelConfig, SpirulaeSplatModel
from spirulae_splat.modules.datamanager import SpirulaeSplatDataManagerConfig, SpirulaeSplatDataManager

from spirulae_splat.modules.dataparser import SpirulaeSplatDataparser, SpirualeSplatDataParserConfig
from spirulae_splat.modules.dataset import SpirulaeSplatDataset

from spirulae_splat.modules.optimizer import create_optimizers

from spirulae_splat.ns_config import _DEFAULT_OPTIMIZERS, _SECOND_ORDER_OPTIMIZERS


@dataclass
class SpirulaeSplatTrainerConfig:
    """🐚 Main training script for spirulae-splat"""

    data: Path
    """Path to dataset. Can be a Nerfstudio or a COLMAP dataset."""

    dataparser: SpirualeSplatDataParserConfig = field(default_factory=SpirualeSplatDataParserConfig)
    """Specifies configurations for data parsing"""

    datamanager: SpirulaeSplatDataManagerConfig = field(default_factory=SpirulaeSplatDataManagerConfig)
    """Specifies configurations for data management during training"""

    model: SpirulaeSplatModelConfig = field(default_factory=SpirulaeSplatModelConfig)
    """Specifies configurations for main model, losses, and densification"""

    optimizer: dict = field(default_factory=lambda: _SECOND_ORDER_OPTIMIZERS)
    """Specifies configurations for optimization"""


class SpirulaeSplatTrainer:

    def __init__(
        self,
        config: SpirulaeSplatTrainerConfig
    ):
        self.config = config
        self.dataparser = SpirulaeSplatDataparser(config.dataparser, config.data)
        self.dataparser_outputs = self.dataparser.parse()
        self.dataset = SpirulaeSplatDataset(self.dataparser_outputs)

        self.datamanager = SpirulaeSplatDataManager(self.config.datamanager, device="cuda")
        self.datamanager.train_dataset = self.dataset

        self.model = SpirulaeSplatModel(
            self.config.model,
            self.dataparser_outputs['metadata'],
            self.dataparser_outputs['cameras']
        ).cuda()

        self.lock = threading.Lock()

    def _render(self, c2w, fx, fy, cx, cy, w, h, camera_model):
        camera = Cameras((fx, fy, cx, cy), [0.0]*10, h, w, torch.from_numpy(c2w), camera_model)
        outputs = self.model.get_outputs(camera)
        # outputs['rgb'] = annotate_train_cameras(outputs['rgb'], outputs['depth'], outputs['alpha'], camera, self.model.cameras)
        return outputs

    def render(self, *args):
        with self.lock:
            self.model.eval()
            return self._render(*args)

    def _get_train_loss_dict(self, step: int):
        inputs = self.datamanager.next_train(step)  # type: List[Tuple[Cameras, Dict]]
        if isinstance(inputs, tuple):
            train_inputs, val_inputs = inputs
            if len(train_inputs) > 1 or len(val_inputs) > 1:
                raise NotImplementedError("Validation with split_batch is not supported")  # TODO
            train_inputs, val_inputs = train_inputs[0], val_inputs[0]
            inputs = [((train_inputs[0], val_inputs[0]), (train_inputs[1], val_inputs[1]))]

        for i, (camera, batch) in enumerate(inputs):
            # torch.cuda.empty_cache()
            model_outputs = self.model.get_outputs(camera)
            # torch.cuda.empty_cache()
            metrics_dict = self.model.get_metrics_dict(model_outputs, batch)
            is_not_last = (i != len(inputs) - 1)
            loss_dict = self.model.get_loss_dict(model_outputs, batch, metrics_dict, no_static_losses=is_not_last)
            # torch.cuda.empty_cache()
            if is_not_last:
                torch.stack([
                    x for x in loss_dict.values() if isinstance(x, torch.Tensor)
                ]).sum().backward()
            # torch.cuda.empty_cache()

        return model_outputs, loss_dict, metrics_dict

    def get_train_loss_dict(self, *args):
        with self.lock:
            self.model.train()
            return self._get_train_loss_dict(*args)

    def train(self):
        optims = create_optimizers(self.model, self.config.optimizer)
        for step in range(30000):
            for optim in optims.values():
                optim.zero_grad()
            self.model.step_cb(optims, step)
            model_outputs, loss_dict, metrics_dict = self.get_train_loss_dict(0)
            loss = torch.stack([x for x in loss_dict.values() if isinstance(x, torch.Tensor)]).sum()
            loss.backward()
            for optim in optims.values():
                optim.step()
            self.model.step_post_backward(step)
