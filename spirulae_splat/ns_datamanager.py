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
from spirulae_splat.modules.training_losses import SplatTrainingLosses
from spirulae_splat.splat.utils import interpolate_se3, random_c2w_on_unit_sphere, ls_camera_intersection

from concurrent.futures import ThreadPoolExecutor
from torch.nn.parallel import DataParallel

from tqdm import tqdm


@dataclass
class SpirulaeDataManagerConfig(FullImageDatamanagerConfig):
    """Template DataManager Config

    Add your custom datamanager config parameters here.
    """
    _target: Type = field(default_factory=lambda: SpirulaeDataManager)

    max_batch_per_epoch: int = 768
    """Maximum number of batches per epoch, used for configuring batch size"""

    patch_batch_size: Optional[int] = None  # 256
    """If not None, batch patches instead of full images
        Leave -1 to let program decide"""

    patch_size: int = 64
    """Patch size used in patch batching
        Make this a multiple of 16 (tile size)
        Affects training speed, optimal value depends on image focal length
        Too small may lead to suboptimal performance in SSIM, as well as increasing VRAM usage / floating point issues
        """

    cache_images: Literal["cpu-pageable", "cpu", "gpu"] = "cpu-pageable"
    """Whether to cache images in memory. If "cpu", caches on cpu. If "gpu", caches on device."""
    cache_images_type: Literal["uint8", "float32"] = "uint8"
    """The image type returned from manager, caching images in uint8 saves memory"""

    deblur_training_images: bool = False
    """Whether to use a custom trained deep learning model to deblur images before training"""

    compute_visibility_masks: bool = False


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
            if self.config.patch_batch_size is not None:
                for key in ['mask', 'depth', 'normal']:
                    if key not in data or data['image'].shape[:2] == data[key].shape[:2]:
                        continue
                    data[key] = torch.nn.functional.interpolate(
                        data[key].float().permute(2, 0, 1)[None],
                        size=data['image'].shape[:2], mode='bilinear', align_corners=False
                    )[0].permute(1, 2, 0).to(data[key].dtype)
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

        if self.config.deblur_training_images:
            from spirulae_splat.modules.enhancer import infer
            for idx, data in enumerate(tqdm(undistorted_images, "Deblurring images")):
                original_shape = data['image'].shape
                data['image'] = infer(data['image'][None])[0]
                H, W, _ = data['image'].shape
                for key in data.keys():
                    if key != 'image' and isinstance(data[key], torch.Tensor) \
                            and data[key].shape[:2] == original_shape[:2]:
                        data[key] = data[key][:H, :W]
                dataset.cameras.width[idx] = W
                dataset.cameras.height[idx] = H

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

        if not hasattr(self, 'image_shapes'):
            # self.image_shapes = torch.stack([self.train_dataset.cameras.height, self.train_dataset.cameras.width], dim=-1)
            self.image_shapes = torch.tensor([im['image'].shape[:2] for im in self.cached_train])
            assert (self.image_shapes >= self.config.patch_size).all(), "Image shape must be at least patch size"
            self.effective_offsets = self.image_shapes - self.config.patch_size
            
            self.pixels_per_image = sum([w*h for (w, h) in self.image_shapes]) / len(self.image_shapes)
            self.images_per_batch = max(len(self.cached_train) / self.config.max_batch_per_epoch, 1)

        if batch_size == -1:
            pixels_per_batch = self.pixels_per_image * self.images_per_batch
            pixels_per_patch = self.config.patch_size**2
            batch_size = max(int(pixels_per_batch // pixels_per_patch), 1)

        image_indices = torch.randint(0, len(self.cached_train), (batch_size,))
        offsets = (torch.rand([batch_size, 2]) * self.effective_offsets[image_indices] + 0.5).long()
        batch = { 'image': [] }
        patch_offsets = []
        for image_idx, (y0, x0) in zip(image_indices, offsets):
            y0, x0 = y0.item(), x0.item()
            y1 = y0 + self.config.patch_size
            x1 = x0 + self.config.patch_size
            image = self.cached_train[image_idx]
            patch = image['image'][y0:y1, x0:x1]
            batch['image'].append(patch)

            for key in ['depth', 'normal', 'mask']:
                if key not in image:
                    continue
                assert image[key].shape[:2] == image['image'].shape[:2], \
                    f"image shape ({image['image'].shape[:2]}) and {key} shape ({image[key].shape[:2]}) mismatch"
                patch = image[key][y0:y1, x0:x1]
                if len(patch.shape) == 2:
                    patch = patch.unsqueeze(-1)
                if key not in batch:
                    batch[key] = []
                batch[key].append(patch)

            patch_offsets.append((x0, y0))

        for key in batch:
            batch[key] = torch.stack(batch[key]).to(self.device)

        camera = self.train_dataset.cameras[image_indices]
        if camera.metadata is None:
            camera.metadata = {}
        camera.metadata['actual_height'] = camera.height.float().mean().item()
        camera.metadata['actual_width'] = camera.width.float().mean().item()
        camera.metadata['actual_images_per_batch'] = self.images_per_batch
        camera.height = self.config.patch_size * torch.ones_like(camera.height)
        camera.width = self.config.patch_size * torch.ones_like(camera.width)
        camera.cx = camera.cx - offsets[:, 1:2]
        camera.cy = camera.cy - offsets[:, 0:1]
        camera = camera.to(self.device)
        # camera.metadata["num_unique_cam_idx"] = len(set(*image_indices))
        camera.metadata["cam_idx"] = image_indices.to(self.device)
        camera.metadata["patch_offsets"] = torch.tensor(patch_offsets, dtype=torch.int32).to(self.device)

        if self.config.compute_visibility_masks:
            camera.metadata['visibility_masks'] = SplatTrainingLosses.get_visibility_masks(batch, self.device)

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

    def random_cameras(self, batch_size: int):
        """Generate random cameras by interpolating existing cameras"""
        num_cameras = len(self.train_dataset.cameras)
        device = self.train_dataset.cameras.device
        if False:
            # interpolate between training cameras
            i0 = torch.randint(0, num_cameras, (batch_size,), device=device)
            i1 = torch.randint(0, num_cameras, (batch_size,), device=device)
            camera_0 = self.train_dataset.cameras[i0]
            camera_1 = self.train_dataset.cameras[i1]
            t = torch.rand(batch_size, device=device)
            camera = Cameras(
                camera_to_worlds=interpolate_se3(camera_0.camera_to_worlds, camera_1.camera_to_worlds, t),
                fx=torch.lerp(camera_0.fx, camera_1.fx, t),
                fy=torch.lerp(camera_0.fy, camera_1.fy, t),
                cx=torch.lerp(camera_0.cx, camera_1.cx, t),
                cy=torch.lerp(camera_0.cy, camera_1.cy, t),
                width=int(torch.lerp(camera_0.width.float(), camera_1.width.float(), t).mean().item()+0.5),
                height=int(torch.lerp(camera_0.height.float(), camera_1.height.float(), t).mean().item()+0.5),
                camera_type=random.choice([camera_0.camera_type, camera_1.camera_type]),
            )
        else:
            # random camera-to-world on unit sphere, useful for object centric scenes
            if not hasattr(self, 'camera_center'):
                self.camera_center = ls_camera_intersection(self.train_dataset.cameras.camera_to_worlds).reshape(1, 3)
            camera_to_worlds = random_c2w_on_unit_sphere(batch_size, device=device)[:, :3, :]
            idx = torch.randint(0, num_cameras, (batch_size,), device=device)
            camera = self.train_dataset.cameras[idx]
            dist = torch.norm(camera.camera_to_worlds[:, :3, 3] - self.camera_center, dim=-1, keepdim=True)
            camera_to_worlds[:, :3, 3] *= dist
            camera_to_worlds[:, :3, 3] += self.camera_center
            camera = Cameras(
                camera_to_worlds=camera_to_worlds,
                fx=camera.fx, fy=camera.fy, cx=camera.cx, cy=camera.cy,
                width=camera.width, height=camera.height,
                camera_type=camera.camera_type, #distortion_params=camera.distortion_params
            )
        camera.metadata = {}
        return camera

    def next_train(self, step: int) -> Tuple[Cameras, Dict]:
        """Returns the next training batch

        Returns a Camera instead of raybundle"""

        if not hasattr(self, 'cached_train'):
            self.cached_train = self._load_images("train", cache_images_device=self.config.cache_images)

        if self.config.patch_batch_size is not None:
            return self.get_tiles(self.config.patch_batch_size)

        train_batch_size = (len(self.train_dataset) + self.config.max_batch_per_epoch - 1) \
            // self.config.max_batch_per_epoch

        # TODO
        if random.random() < (step - 10000) / (30000 - 10000) and False:
            return self.random_cameras(train_batch_size), {}

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
            # TODO: handle cases where different img has different set of auxiliary buffers
            img = self.cached_train[image_idx]
            for key, value in img.items():
                if key not in batch:
                    batch[key] = []
                if isinstance(value, torch.Tensor):
                    value = value.clone().cuda()
                batch[key].append(value)
            assert 'image' in img
            h, w, _ = batch['image'][-1].shape
            for key in ['depth', 'normal', 'mask']:
                if key not in batch:
                    continue
                if len(batch[key][-1].shape) == 2:
                    batch[key][-1] = batch[key][-1].unsqueeze(-1)
                batch[key][-1] = torch.nn.functional.interpolate(
                    batch[key][-1].float().permute(2, 0, 1)[None],
                    size=(h, w), mode='bilinear', align_corners=False
                )[0].permute(1, 2, 0).to(batch[key][-1].dtype)
        for key in batch:
            if isinstance(batch[key][0], torch.Tensor):
                batch[key] = torch.stack(batch[key]).to(self.device)

        camera = self.train_dataset.cameras[torch.tensor(image_indices)].to(self.device)
        if camera.metadata is None:
            camera.metadata = {}
        camera.metadata["cam_idx"] = image_indices

        if self.config.compute_visibility_masks:
            camera.metadata['visibility_masks'] = SplatTrainingLosses.get_visibility_masks(batch, self.device)

        return camera, batch

    @cached_property
    def dataset_type(self) -> Type[TDataset]:
        """Returns the dataset type passed as the generic argument"""
        return SpirulaeDataset

