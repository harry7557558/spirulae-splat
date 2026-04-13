from dataclasses import dataclass, field
from typing import Type, List, Tuple, Dict
from pathlib import Path

import torch

from spirulae_splat.modules.camera import Cameras
from spirulae_splat.modules.model import SpirulaeSplatModelConfig, SpirulaeSplatModel
from spirulae_splat.modules.datamanager import SpirulaeSplatDataManagerConfig, SpirulaeSplatDataManager

from spirulae_splat.modules.dataparser import SpirulaeSplatDataparser, SpirualeSplatDataParserConfig
from spirulae_splat.modules.dataset import SpirulaeSplatDataset


@dataclass
class SpirulaeSplatTrainerConfig:
    """🐚 Main training script for spirulae-splat"""

    data: Path
    """Path to dataset. Can be a Nerfstudio or a COLMAP dataset."""

    dataparser: SpirualeSplatDataParserConfig = field(default_factory=SpirualeSplatDataParserConfig)
    """Specifies the dataparser config"""

    datamanager: SpirulaeSplatDataManagerConfig = field(default_factory=SpirulaeSplatDataManagerConfig)
    """specifies the datamanager config"""

    model: SpirulaeSplatModelConfig = field(default_factory=SpirulaeSplatModelConfig)
    """specifies the model config"""


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



def entrypoint():
    import tyro
    config = tyro.cli(SpirulaeSplatTrainerConfig)
    trainer = SpirulaeSplatTrainer(config)

    for i in range(1000):
        model_outputs, loss_dict, metrics_dict = trainer.get_train_loss_dict(0)
        loss = torch.stack([x for x in loss_dict.values() if isinstance(x, torch.Tensor)]).sum()
        loss.backward()

if __name__ == "__main__":
    entrypoint()
