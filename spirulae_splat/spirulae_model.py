"""
Template Model File

Currently this subclasses the Nerfacto model. Consider subclassing from the base Model.
"""
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Type

from nerfstudio.models.splatfacto import SplatfactoModel, SplatfactoModelConfig  # for subclassing Nerfacto model
from nerfstudio.models.base_model import Model, ModelConfig  # for custom Model

from nerfstudio.cameras.rays import RayBundle, RaySamples


@dataclass
class SpirulaeModelConfig(SplatfactoModelConfig):
    """Template Model Configuration.

    Add your custom model config parameters here.
    """

    _target: Type = field(default_factory=lambda: SpirulaeModel)


class SpirulaeModel(SplatfactoModel):
    """Template Model."""

    config: SpirulaeModelConfig
