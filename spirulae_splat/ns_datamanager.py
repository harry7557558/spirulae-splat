"""
Template DataManager
"""

from dataclasses import dataclass, field
from typing import Dict, Literal, List, Tuple, Type, Union, Optional, Callable
from copy import deepcopy
import random
import math
import os
import hashlib
import re
import json
import time

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

from spirulae_splat.ns_dataset import SpirulaeDataset, IndexedDatasetWrapper
from spirulae_splat.modules.training_losses import SplatTrainingLosses
from spirulae_splat.splat.utils import interpolate_se3, random_c2w_on_unit_sphere, ls_camera_intersection

from concurrent.futures import ThreadPoolExecutor
from torch.nn.parallel import DataParallel
from torch.utils.data import DataLoader
from torch.utils.data._utils.collate import default_collate

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

    cache_images: Literal["cpu-pageable", "cpu", "gpu", "disk"] = "cpu-pageable"
    """Image cache location. If "cpu", caches on cpu. If "gpu", caches on device.
        If "cpu-pageable", cache on cpu pageable memory (saves RAM but may cause error if spill to swap memory).
        If "disk", cache on disk (limited support). """
    cache_images_type: Literal["uint8", "float32"] = "uint8"
    """The image type returned from manager, caching images in uint8 saves memory"""

    deblur_training_images: bool = False
    """Whether to use a custom trained deep learning model to deblur images before training"""

    compute_visibility_masks: bool = False


class IndexGroup:
    def __init__(self, indices: List[int]):
        self.indices = indices[:]
        self.ptr = len(indices)

    def __len__(self):
        return len(self.indices)

    def __next__(self):
        if self.ptr >= len(self.indices):
            random.shuffle(self.indices)
            self.ptr = 0
        idx = self.indices[self.ptr]
        self.ptr += 1
        return idx


class IndexGroups:
    def __init__(self, indices: List[int], keys: List):
        groups = {}
        for idx, key in zip(indices, keys):
            if key not in groups:
                groups[key] = []
            groups[key].append(idx)
        self.groups = []  # type: List[IndexGroup]
        for key, idxs in groups.items():
            self.groups.append(IndexGroup(idxs))
        self.probs = np.array([float(len(g)) for g in self.groups])
        self.probs /= np.sum(self.probs)

    def random_idx(self) -> int:
        return np.random.choice(range(len(self.groups)), p=self.probs)


class IndexGroupsWithDataLoader(IndexGroups):
    def __init__(self, indices: List[int], keys: List, getitem: Callable, batch_size: int, parallel: bool):
        super().__init__(indices, keys)
        self.batch_size = batch_size
        self.getitem = getitem
        self.parallel = parallel

        if not self.parallel:
            return

        self.dataloaders = []
        max_num_workers = min(max(os.cpu_count()-2,1), batch_size)
        for group in self.groups:
            dataset = IndexedDatasetWrapper(getitem, [[i] for i in group.indices])
            dataloader = DataLoader(
                dataset,
                batch_size=min(batch_size, len(group)),
                num_workers=min(max_num_workers, len(group)),
                persistent_workers=True,
                shuffle=True
            )
            self.dataloaders.append([dataloader, iter(dataloader)])

    def get_batch_pytorch(self):
        dataloader = self.dataloaders[self.random_idx()]
        try:
            return next(dataloader[1])
        except StopIteration:
            dataloader[1] = iter(dataloader[0])
            return next(dataloader[1])

    def get_batch(self):
        if self.parallel:
            return self.get_batch_pytorch()
        group = self.groups[self.random_idx()]
        batch = []
        for i in range(min(self.batch_size, len(group))):
            batch.append(self.getitem(next(group)))
        return default_collate(batch)



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
        cache_images = config.cache_images
        super().__init__(
            config=config, device=device, test_mode=test_mode, world_size=world_size, local_rank=local_rank, **kwargs
        )
        self.config.cache_images = cache_images

        self.train_index_group_loader = None  # type: Optional[IndexGroupsWithDataLoader]

    def _undistort_idx(self, dataset, idx: int, return_idx=False) -> Dict[str, torch.Tensor]:
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
        if return_idx:
            return idx, data
        return data

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

    def _load_images(
        self, split: Literal["train", "eval"], cache_images_device: Literal["cpu-pageable", "cpu", "gpu", "disk"]
    ) -> List[Dict[str, torch.Tensor]]:
        undistorted_images: List[Dict[str, torch.Tensor]] = []

        # Which dataset?
        if split == "train":
            dataset = self.train_dataset
        elif split == "eval":
            dataset = self.eval_dataset
        else:
            assert_never(split)

        if cache_images_device == "disk":
            return []

        CONSOLE.log(f"Caching {split} images")
        with ThreadPoolExecutor() as executor:
            undistorted_images = list(
                track(
                    executor.map(
                        self._undistort_idx,
                        [dataset]*len(dataset),
                        range(len(dataset)),
                    ),
                    description=f"Caching {split} images",
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

    def get_train_image(self, idx: int):
        if self.config.cache_images == "disk":
            batch = self._undistort_idx(self.train_dataset, idx)
        else:
            batch = self.cached_train[idx]

        assert 'image' in batch
        for key in ['mask', 'depth', 'normal']:
            if key not in batch or batch['image'].shape[:2] == batch[key].shape[:2]:
                continue
            batch[key] = torch.nn.functional.interpolate(
                batch[key].float().permute(2, 0, 1)[None],
                size=batch['image'].shape[:2], mode='bilinear', align_corners=False
            )[0].permute(1, 2, 0).to(batch[key].dtype)

        camera = self.train_dataset.cameras[idx]
        if camera.metadata is None:
            camera.metadata = {}
        camera.metadata["cam_idx"] = idx

        if self.config.compute_visibility_masks:
            camera.metadata['visibility_masks'] = SplatTrainingLosses.get_visibility_masks(batch, self.device)

        camera_flattened = {}
        for key in ['camera_to_worlds', 'fx', 'fy', 'cx', 'cy', 'width', 'height', 'distortion_params', 'camera_type', 'times', 'metadata']:
            value = getattr(camera, key)
            if value is None:
                value = torch.empty((0,))
            camera_flattened[key] = value
        return camera_flattened, batch

    def get_tiles(self, batch_size: int):
        # TODO: multithread with pytorch data loader

        TILE_SIZE = 16
        assert self.config.patch_size > 0 and self.config.patch_size % TILE_SIZE == 0, "Patch size must be a multiple of tile size"

        # TODO: use IndexGroups with camera type support
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

    def setup_index_group_loaders(self):
        """Make separate lists for images with different shapes and camera models"""

        def get_key(idx):
            cam = self.train_dataset.cameras[idx]
            w, h, camera_type = cam.width.item(), cam.height.item(), cam.camera_type.item()
            additional_params = []
            for key in ['camera_to_worlds', 'fx', 'fy', 'cx', 'cy', 'distortion_params', 'times']:
                value = getattr(cam, key)  # type: Optional[torch.Tensor]
                additional_params.append(value if value is None else tuple(value.shape))
            return (w, h, camera_type, *additional_params)

        val_indices = self.train_dataset.val_indices
        train_indices = sorted(set(range(len(self.train_dataset))).difference(val_indices))
        self.train_index_group_loader = IndexGroupsWithDataLoader(
            train_indices, [get_key(idx) for idx in train_indices],
            self.get_train_image, self.train_batch_size(False), self.config.cache_images != "gpu"
        )
        self.val_index_group_loader = IndexGroupsWithDataLoader(
            val_indices, [get_key(idx) for idx in val_indices],
            self.get_train_image, self.val_batch_size(False), self.config.cache_images != "gpu"
        )

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

        if not hasattr(self, 'cached_train') and self.config.cache_images != "disk":
            self.cached_train = self._load_images("train", cache_images_device=self.config.cache_images)

        if self.config.patch_batch_size is not None:
            assert self.config.cache_images != "disk", "Disk caching not supported in patch batching mode"
            assert len(self.train_dataset.val_indices) == 0, "Validation is not supported in patch batching mode"
            return self.get_tiles(self.config.patch_batch_size)

        # TODO
        if random.random() < (step - 10000) / (30000 - 10000) and False:
            return self.random_cameras(self.train_batch_size()), {}

        if self.train_index_group_loader is None:
            self.setup_index_group_loaders()

        camera, batch = self.train_index_group_loader.get_batch()
        for key, value in camera.items():
            if isinstance(value, torch.Tensor) and value.numel() == 0:
                camera[key] = None
        camera = Cameras(**camera)

        camera = camera.to(self.device)
        for key, value in batch.items():
            if isinstance(value, torch.Tensor):
                batch[key] = value.to(self.device)

        val_batch_size = self.val_batch_size(True)
        if len(self.train_dataset.val_indices) > 0 and val_batch_size > 0:
            val_camera, val_batch = self.val_index_group_loader.get_batch()
            for key, value in val_camera.items():
                if isinstance(value, torch.Tensor):
                    if value.numel() == 0:
                        val_camera[key] = None
                    else:
                        val_camera[key] = value[:val_batch_size, ...].to(self.device)
            val_camera = Cameras(**val_camera)
            for key, value in val_batch.items():
                if isinstance(value, torch.Tensor):
                    val_batch[key] = value[:val_batch_size].to(self.device)
            return (camera, val_camera), (batch, val_batch)

        return camera, batch

    @cached_property
    def dataset_type(self) -> Type[TDataset]:
        """Returns the dataset type passed as the generic argument"""
        return SpirulaeDataset

