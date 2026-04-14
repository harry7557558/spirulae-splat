from dataclasses import dataclass, field
from typing import Type, List, Tuple, Dict
from pathlib import Path

import torch

from spirulae_splat.modules.camera import Cameras
from spirulae_splat.modules.model import SpirulaeSplatModelConfig, SpirulaeSplatModel
from spirulae_splat.modules.datamanager import SpirulaeSplatDataManagerConfig, SpirulaeSplatDataManager

from spirulae_splat.modules.dataparser import SpirulaeSplatDataparser, SpirualeSplatDataParserConfig
from spirulae_splat.modules.dataset import SpirulaeSplatDataset

from spirulae_splat.viewer.server import ViewerServer, SliderDef, DropdownDef


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

    def render(self, c2w, fx, fy, cx, cy, w, h, camera_model):
        camera = Cameras((fx, fy, cx, cy), [0.0]*10, h, w, torch.from_numpy(c2w), camera_model)
        self.model.eval()
        return self.model.get_outputs(camera)

    def get_train_loss_dict(self, step: int):
        inputs = self.datamanager.next_train(step)  # type: List[Tuple[Cameras, Dict]]
        if isinstance(inputs, tuple):
            train_inputs, val_inputs = inputs
            if len(train_inputs) > 1 or len(val_inputs) > 1:
                raise NotImplementedError("Validation with split_batch is not supported")  # TODO
            train_inputs, val_inputs = train_inputs[0], val_inputs[0]
            inputs = [((train_inputs[0], val_inputs[0]), (train_inputs[1], val_inputs[1]))]

        for i, (camera, batch) in enumerate(inputs):
            self.model.train()
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

    server = ViewerServer(
        render_fn=trainer.render,
        http_host="localhost",
        http_port=7007,
        ws_host="localhost",
        ws_port=8765,
        jpeg_quality=75,
        # Example extra controls:
        # extra_sliders=[
        #     SliderDef(id="brightness", label="Brightness", min=0.5, max=2.0, step=0.05, value=1.0, unit="×"),
        #     SliderDef(id="anim_speed", label="Anim Speed", min=0.0, max=5.0, step=0.1, value=1.0, unit="×"),
        # ],
        # extra_dropdowns=[
        #     DropdownDef(
        #         id="shading_mode",
        #         label="Shading",
        #         options=[
        #             {"value": "full", "label": "Full"},
        #             {"value": "diffuse", "label": "Diffuse"},
        #             {"value": "specular", "label": "Specular"},
        #         ],
        #         default="full",
        #     ),
        # ],
        open_browser=False,
    )

    server.start()
    server.wait()

    return
    for i in range(1000):
        model_outputs, loss_dict, metrics_dict = trainer.get_train_loss_dict(0)
        loss = torch.stack([x for x in loss_dict.values() if isinstance(x, torch.Tensor)]).sum()
        loss.backward()

if __name__ == "__main__":
    entrypoint()
