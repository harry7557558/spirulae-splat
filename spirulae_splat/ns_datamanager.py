"""
Template DataManager
"""

from dataclasses import dataclass, field
from typing import Dict, Literal, List, Tuple, Type, Union, Optional, TypeVar, Generic
from copy import deepcopy
import random
import math
import os
import hashlib
import re
import json

import torch
import numpy as np

from collections import deque
from functools import cached_property

from nerfstudio.cameras.rays import RayBundle
from nerfstudio.data.datamanagers.full_images_datamanager import (
    DataManagerConfig,
    FullImageDatamanager,
    FullImageDatamanagerConfig,
    _undistort_image,
    CONSOLE, track, TDataset
)
from nerfstudio.cameras.cameras import Cameras, CameraType

from spirulae_splat.ns_dataset import SpirulaeDataset

from spirulae_splat.splat._torch_impl import depth_map, depth_inv_map


from concurrent.futures import ThreadPoolExecutor
from torch.nn.parallel import DataParallel


@dataclass
class SpirulaeDataManagerConfig(FullImageDatamanagerConfig):
    """Template DataManager Config

    Add your custom datamanager config parameters here.
    """
    _target: Type = field(default_factory=lambda: SpirulaeDataManager)

    max_batch_per_epoch: int = 768
    """Maximum number of batches per epoch, used for configuring batch size"""

    patch_batch_size: Optional[int] = None  # 256
    """If set, batch patches instead of full images
        Make this a multiple of 16 (tile size)
        Affects training speed, optimal value depends on image focal length
        Too small may lead to suboptimal performance in SSIM, as well as increasing VRAM usage / floating point issues"""

    patch_size: int = 64
    """Patch size used in patch batching
        Make number of pixels (patch_batch_size * patch_size**2) consistent (e.g. 1M)
        """

    cache_images: Literal["cpu-pageable", "cpu", "gpu"] = "cpu-pageable"
    """Whether to cache images in memory. If "cpu", caches on cpu. If "gpu", caches on device."""
    cache_images_type: Literal["uint8", "float32"] = "uint8"
    """The image type returned from manager, caching images in uint8 saves memory"""


class SpirulaeDataManager(FullImageDatamanager):
    """Template DataManager

    Args:
        config: the DataManagerConfig used to instantiate class
    """

    config: SpirulaeDataManagerConfig

    def __init__(
        self,
        config: SpirulaeDataManagerConfig,
        device: Union[torch.device, str] = "cpu",
        test_mode: Literal["test", "val", "inference"] = "val",
        world_size: int = 1,
        local_rank: int = 0,
        **kwargs,  # pylint: disable=unused-argument
    ):
        super().__init__(
            config=config, device=device, test_mode=test_mode, world_size=world_size, local_rank=local_rank, **kwargs
        )

        self.train_cameras_splits = []  # type: list[list]
        self.train_unseen_cameras = []  # type: list[deque]

    def _load_images(
        self, split: Literal["train", "eval"], cache_images_device: Literal["cpu-pageable", "cpu", "gpu"]
    ) -> List[Dict[str, torch.Tensor]]:
        undistorted_images: List[Dict[str, torch.Tensor]] = []

        # Which dataset?
        if split == "train":
            dataset = self.train_dataset
        elif split == "eval":
            dataset = self.eval_dataset
        else:
            assert_never(split)

        def undistort_idx(idx: int) -> Dict[str, torch.Tensor]:
            data = dataset.get_data(idx, image_type=self.config.cache_images_type, _is_viewer=False)
            dataset.cameras.width[idx] = data["image"].shape[1]
            dataset.cameras.height[idx] = data["image"].shape[0]
            camera = dataset.cameras[idx].reshape(())
            assert data["image"].shape[1] == camera.width.item() and data["image"].shape[0] == camera.height.item(), (
                f'The size of image ({data["image"].shape[1]}, {data["image"].shape[0]}) loaded '
                f'does not match the camera parameters ({camera.width.item(), camera.height.item()})'
            )
            if camera.distortion_params is None or torch.all(camera.distortion_params == 0):
                return data
            K = camera.get_intrinsics_matrices().numpy()
            distortion_params = camera.distortion_params.numpy()

            if camera.camera_type.item() == CameraType.FISHEYE.value:
                # don't undistort
                mask = data["mask"] if "mask" in data else None
            else:
                image = data["image"].numpy()
                K, image, mask = _undistort_image(camera, distortion_params, data, image, K)
                data["image"] = torch.from_numpy(image)
            if mask is not None:
                data["mask"] = mask

            dataset.cameras.fx[idx] = float(K[0, 0])
            dataset.cameras.fy[idx] = float(K[1, 1])
            dataset.cameras.cx[idx] = float(K[0, 2])
            dataset.cameras.cy[idx] = float(K[1, 2])
            dataset.cameras.width[idx] = data["image"].shape[1]
            dataset.cameras.height[idx] = data["image"].shape[0]
            return data

        CONSOLE.log(f"Caching/undistorting {split} images")
        with ThreadPoolExecutor() as executor:
            undistorted_images = list(
                track(
                    executor.map(
                        undistort_idx,
                        range(len(dataset)),
                    ),
                    description=f"Caching/undistorting {split} images",
                    transient=True,
                    total=len(dataset),
                )
            )

        # Move to device.
        if cache_images_device == "gpu":
            for cache in undistorted_images:
                cache["image"] = cache["image"].to(self.device)
                if "mask" in cache:
                    cache["mask"] = cache["mask"].to(self.device)
                if "depth" in cache:
                    cache["depth"] = cache["depth"].to(self.device)
                if "normal" in cache:
                    cache["normal"] = cache["normal"].to(self.device)
                self.train_cameras = self.train_dataset.cameras.to(self.device)
        elif cache_images_device.startswith("cpu"):
            for cache in undistorted_images:
                if cache_images_device == "cpu":  # not cpu-pageable
                    cache["image"] = cache["image"].pin_memory()
                    if "mask" in cache:
                        cache["mask"] = cache["mask"].pin_memory()
                    if "depth" in cache:
                        cache["depth"] = cache["depth"].pin_memory()
                    if "normal" in cache:
                        cache["normal"] = cache["normal"].pin_memory()
                self.train_cameras = self.train_dataset.cameras
        else:
            assert False, "Invalid cache image device: " + cache_images_device

        return undistorted_images

    def get_tiles(self, batch_size: int):
        TILE_SIZE = 16
        assert self.config.patch_size > 0 and self.config.patch_size % TILE_SIZE == 0, "Patch size must be a multiple of tile size"

        # image_shapes = torch.stack([self.train_dataset.cameras.height, self.train_dataset.cameras.width], dim=-1)
        image_shapes = torch.tensor([im['image'].shape[:2] for im in self.cached_train])
        assert (image_shapes >= self.config.patch_size).all(), "Image shape must be at least patch size"
        effective_offsets = image_shapes - self.config.patch_size

        image_indices = torch.randint(0, len(self.cached_train), (batch_size,))
        offsets = (torch.rand([batch_size, 2]) * effective_offsets[image_indices] + 0.5).long()
        patches = []
        masks = []
        has_mask = False
        default_mask = torch.ones((self.config.patch_size, self.config.patch_size, 1), dtype=torch.bool)
        offset_slices = []
        for image_idx, (y0, x0) in zip(image_indices, offsets):
            y0, x0 = y0.item(), x0.item()
            y1 = y0 + self.config.patch_size
            x1 = x0 + self.config.patch_size
            image = self.cached_train[image_idx]
            patch = image['image'][y0:y1, x0:x1]
            mask = default_mask
            if 'mask' in image:
                assert image['mask'].shape[:2] == image['image'].shape[:2], \
                    f"image shape ({image['image'].shape[:2]}) and mask shape ({image['mask'].shape[:2]}) mismatch"
                mask = image['mask'][y0:y1, x0:x1]
                has_mask = True
            # TODO: depth, normal
            patches.append(patch)
            masks.append(mask)
            offset_slices.append((slice(y0, y1), slice(x0, x1)))
        batch = { 'image': torch.stack(patches).to(self.device) }
        if has_mask:
            batch['mask'] = torch.stack(masks).to(self.device)

        camera = self.train_dataset.cameras[image_indices]
        camera.height = self.config.patch_size * torch.ones_like(camera.height)
        camera.width = self.config.patch_size * torch.ones_like(camera.width)
        camera.cx = camera.cx - offsets[:, 1:2]
        camera.cy = camera.cy - offsets[:, 0:1]
        camera = camera.to(self.device)
        if camera.metadata is None:
            camera.metadata = {}
        camera.metadata["cam_idx"] = image_indices
        camera.metadata["slices"] = offset_slices

        return camera, batch

    def setup_batches(self):
        """Make separate lists for images with different shapes and camera models"""
        key_map = {}
        self.train_cameras_splits = []
        self.train_unseen_cameras = []
        for cam_idx in range(len(self.train_dataset)):
            cam = self.train_dataset.cameras[cam_idx]
            w, h = cam.width.item(), cam.height.item()
            model = cam.camera_type.item()
            key = (w, h, model)
            if key not in key_map:
                key_map[key] = len(self.train_cameras_splits)
                self.train_cameras_splits.append([])
                self.train_unseen_cameras.append(deque())
            key = key_map[key]
            self.train_cameras_splits[key].append(cam_idx)

        self.camera_split_weights = np.array([len(s) for s in self.train_cameras_splits]) / len(self.train_dataset)

    def next_train(self, step: int) -> Tuple[Cameras, Dict]:
        """Returns the next training batch

        Returns a Camera instead of raybundle"""

        if self.config.patch_batch_size is not None:
            return self.get_tiles(self.config.patch_batch_size)

        train_batch_size = (len(self.train_dataset) + self.config.max_batch_per_epoch - 1) \
            // self.config.max_batch_per_epoch
        # train_batch_size = 4

        if len(self.train_cameras_splits) == 0:
            self.setup_batches()

        image_indices = []
        split = np.random.choice(range(len(self.train_cameras_splits)), p=self.camera_split_weights)
        train_batch_size = min(train_batch_size, len(self.train_cameras_splits[split]))
        for i in range(train_batch_size):
            if len(self.train_unseen_cameras[split]) == 0:
                perm = self.train_cameras_splits[split][:]
                random.shuffle(perm)
                self.train_unseen_cameras[split].extend(perm)
            image_idx = self.train_unseen_cameras[split].popleft()
            image_indices.append(image_idx)
        image_indices = list(set(image_indices))  # TODO: better way for this

        batch = {}
        for image_idx in image_indices:
            img = self.cached_train[image_idx]
            for key, value in img.items():
                if key not in batch:
                    batch[key] = []
                if isinstance(value, torch.Tensor):
                    value = value.clone()
                batch[key].append(value)
        for key in batch:
            if isinstance(batch[key][0], torch.Tensor):
                batch[key] = torch.stack(batch[key]).to(self.device)

        camera = self.train_dataset.cameras[torch.tensor(image_indices)].to(self.device)
        if camera.metadata is None:
            camera.metadata = {}
        camera.metadata["cam_idx"] = image_indices
        return camera, batch

    @cached_property
    def dataset_type(self) -> Type[TDataset]:
        """Returns the dataset type passed as the generic argument"""
        return SpirulaeDataset

