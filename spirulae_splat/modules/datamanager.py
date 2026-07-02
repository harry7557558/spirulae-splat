"""
Template DataManager
"""

from dataclasses import dataclass, field
from typing import Dict, Literal, List, Tuple, Union, Optional, Callable
from copy import deepcopy
import random
import math
import os

import torch
import numpy as np

from spirulae_splat.modules.camera import Cameras, CameraType
from spirulae_splat.modules.resample import warp_equirectangular_to_pinhole, warp_wide_to_pinhole

from spirulae_splat.modules.dataset import IndexedDatasetWrapper

from concurrent.futures import ThreadPoolExecutor
from torch.utils.data import DataLoader
from torch.utils.data._utils.collate import default_collate

from tqdm import tqdm


@dataclass
class SpirulaeSplatDataManagerConfig:
    """Configuration for Spirulae-Splat data manager."""

    max_batch_per_epoch: int = 800
    """Maximum number of batches per epoch, used for configuring batch size"""
    split_batch: bool = False
    """Whether to one large batch into many small batches to avoid OOM, at cost of slower training"""

    cache_images: Literal["cpu-pageable", "cpu", "gpu", "disk"] = "disk"
    """Image cache location. If "cpu", caches on cpu. If "gpu", caches on device.
        If "cpu-pageable", cache on cpu pageable memory (saves RAM but may cause error if spill to swap memory).
        If "disk", cache on disk (limited support). """

    load_depths: bool = True
    """Whether to load depth maps, if exist"""
    load_normals: bool = True
    """Whether to load normal maps, if exist"""

    mask_boundary_offset: float = 0
    """Signed boundary offset applied to binarized masks at decode time, as a
        fraction of sqrt(W*H) of the decoded mask. Positive dilates (grows)
        foreground, negative erodes (shrinks). Runs on CPU during data loading
        via separable Felzenszwalb-Huttenlocher squared-Euclidean DT (exact,
        O(N) per row + col)."""

    warp_to_pinhole: bool = False
    """Whether to split an image into 5 undistorted pinhole images.
        Can sometimes give better quality and compatibility for dataset captured by fisheye/360 cameras."""

    deblur_training_images: bool = False
    """Whether to use a custom trained deep learning model to deblur images before training"""



class SpirulaeSplatDataManager:
    """Data manager for spirulae-splat"""

    config: SpirulaeSplatDataManagerConfig

    def __init__(
        self,
        config: SpirulaeSplatDataManagerConfig,
        device: Union[torch.device, str] = "cpu",
        eval: bool = False,
        **kwargs,  # pylint: disable=unused-argument
    ):
        self.config = config
        self.device = device

        self.train_dataset = None
        self.eval_dataset = None

        self.eval = eval

    def train_batch_size(self, stochastic=True):
        n = (len(self.train_dataset) - len(self.train_dataset.val_indices)) / self.config.max_batch_per_epoch
        n = max(n, 1.0)
        if stochastic:
            return int(n) + int(random.random() < n-int(n))
        return int(n+0.5)

    def val_batch_size(self, stochastic=True):
        num_val = len(self.train_dataset.val_indices)
        num_train = len(self.train_dataset) - num_val
        n_train = max(num_train / self.config.max_batch_per_epoch, 1.0)
        n = n_train * num_val / num_train
        if stochastic:
            return int(n) + int(random.random() < n-int(n))
        return int(math.ceil(n))

    def next_train(self, step: int, indices=None) -> List[Tuple[Cameras, Dict]]:
        """Returns the next sample as a single-image batch.

        Training data is fed from the C++ engine; only the sequential eval
        path is ported here. Each call returns the next image in dataset
        order (wrapping around), matching the trainer's per-camera eval loop.
        """
        idx = getattr(self, "_next_train_idx", 0)
        self._next_train_idx = (idx + 1) % len(self.train_dataset)

        data = self.train_dataset.get_data(
            idx,
            load_depths=self.config.load_depths,
            load_normals=self.config.load_normals,
        )

        # Index with a 1-D tensor to keep the leading batch dimension the
        # model's get_outputs expects (camera.width[0], camera_type[0], ...).
        camera = self.train_dataset.cameras[torch.tensor([idx])].to(self.device)

        batch = {"image": data["image"].unsqueeze(0).to(self.device)}
        return [(camera, batch)]
