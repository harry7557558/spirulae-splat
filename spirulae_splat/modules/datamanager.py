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

    cache_images: Literal["cpu-pageable", "cpu", "gpu", "disk"] = "cpu-pageable"
    """Image cache location. If "cpu", caches on cpu. If "gpu", caches on device.
        If "cpu-pageable", cache on cpu pageable memory (saves RAM but may cause error if spill to swap memory).
        If "disk", cache on disk (limited support). """

    load_depths: bool = True
    """Whether to load depth maps, if exist"""
    load_normals: bool = True
    """Whether to load normal maps, if exist"""

    warp_to_pinhole: bool = False
    """Whether to split an image into 5 undistorted pinhole images.
        Can sometimes give better quality and compatibility for dataset captured by fisheye/360 cameras."""

    deblur_training_images: bool = False
    """Whether to use a custom trained deep learning model to deblur images before training"""

    use_cpp_data_manager: bool = False
    """If True, route training data loading + per-step camera/GT prep through the
        C++ DataManager (engine_setup_data_manager + engine_train_step_managed).
        The Python SpirulaeSplatDataManager is still constructed for the
        viewer / eval / render paths, but the training loop bypasses it. Maps
        ``cache_images`` -> CacheMode as: "cpu" / "cpu-pageable" -> CPU,
        "disk" -> DISK. ``warp_to_pinhole``, ``split_batch``,
        ``model.use_camera_optimizer``, supersampling > 1, and the BVH
        primitive are not supported on this path; the trainer raises a clear
        error when those are combined with the flag."""


class IndexGroup:
    def __init__(self, indices: List[int], eval: bool):
        self.indices = indices[:]
        self.eval = eval
        if not self.eval:
            random.shuffle(self.indices)
        self.ptr = len(indices)

    def __len__(self):
        return len(self.indices)

    def __next__(self):
        if self.ptr >= len(self.indices):
            if self.eval:
                random.shuffle(self.indices)
            self.ptr = 0
        idx = self.indices[self.ptr]
        self.ptr += 1
        return idx


class IndexGroups:
    def __init__(self, indices: List[int], keys: List, eval: bool):
        groups = {}
        for idx, key in zip(indices, keys):
            if key not in groups:
                groups[key] = []
            groups[key].append(idx)
        self.groups = []  # type: List[IndexGroup]
        for key, idxs in groups.items():
            self.groups.append(IndexGroup(idxs, eval))
        self.probs = np.array([float(len(g)) for g in self.groups])
        self.probs /= np.sum(self.probs)

    def random_idx(self) -> int:
        return np.random.choice(range(len(self.groups)), p=self.probs)


class IndexGroupsWithDataLoader(IndexGroups):
    def __init__(self, indices: List[int], keys: List, getitem: Callable, batch_size: int, parallel: bool, eval: bool):
        super().__init__(indices, keys, eval)
        self.batch_size = batch_size
        self.getitem = getitem
        self.parallel = parallel

        if not self.parallel:
            return

        self.dataloaders = []
        min_num_workers = 4
        max_num_workers = min(max(os.cpu_count()-2,1), max(batch_size, min_num_workers))
        while max_num_workers * len(self.groups) > 2*os.cpu_count():
            max_num_workers -= 1
            if max_num_workers in [0, 1]:
                self.parallel = False
                return  # TODO: concurrency
        for group in self.groups:
            dataset = IndexedDatasetWrapper(getitem, [[i] for i in group.indices])
            n_workers = min(max_num_workers, max(len(group), (min_num_workers+len(self.groups)-1)//len(self.groups)))
            dataloader = DataLoader(
                dataset,
                batch_size=min(batch_size, len(group)),
                num_workers=n_workers,
                persistent_workers=True,
                shuffle=(not eval),
                prefetch_factor=4 if n_workers > 0 else None,
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



class SpirulaeSplatDataManager:
    """Data manager for spirulae-splat"""

    config: SpirulaeSplatDataManagerConfig

    def __init__(
        self,
        config: SpirulaeSplatDataManagerConfig,
        device: Union[torch.device, str] = "cpu",
        test_mode: Literal["test", "val", "inference"] = "val",
        world_size: int = 1,
        local_rank: int = 0,
        eval: bool = False,
        **kwargs,  # pylint: disable=unused-argument
    ):
        self.config = config
        self.device = device

        self.train_dataset = None
        self.eval_dataset = None

        self.eval = eval

        self.train_index_group_loader = None  # type: Optional[IndexGroupsWithDataLoader]

    def _undistort_idx(self, dataset, idx: int, return_idx=False) -> Dict[str, torch.Tensor]:
        data = dataset.get_data(idx, load_depths=self.config.load_depths, load_normals=self.config.load_normals)
        dataset.cameras.width[idx] = data["image"].shape[1]
        dataset.cameras.height[idx] = data["image"].shape[0]
        camera = dataset.cameras[idx]
        assert data["image"].shape[1] == camera.width.item() and data["image"].shape[0] == camera.height.item(), (
            f'The size of image ({data["image"].shape[1]}, {data["image"].shape[0]}) loaded '
            f'does not match the camera parameters ({camera.width.item(), camera.height.item()})'
        )
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

        print(f"Caching {split} images")
        with ThreadPoolExecutor() as executor:
            with ThreadPoolExecutor() as executor:
                for result in tqdm(
                    executor.map(
                        self._undistort_idx,
                        [dataset]*len(dataset),
                        range(len(dataset)),
                    ),
                    f"Caching {split} images",
                    total=len(dataset),
                ):
                    undistorted_images.append(result)

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

        camera = self.train_dataset.cameras[idx]  # type: Cameras
        if camera.metadata is None:
            camera.metadata = {}
        camera.metadata["cam_idx"] = int(idx)

        camera_flattened = {}
        for key in ['camera_to_worlds', 'intrins', 'width', 'height', 'distortion_params', 'camera_type', 'metadata']:
            value = getattr(camera, key)
            if isinstance(value, torch.Tensor):
                value = value.cpu()
            if value is None:
                value = torch.empty((0,), device="cpu")
            camera_flattened[key] = value
        return camera_flattened, batch

    def setup_index_group_loaders(self):
        """Make separate lists for images with different shapes and camera models"""

        def get_key(idx):
            cam = self.train_dataset.cameras[idx]
            w, h, camera_type = cam.width.item(), cam.height.item(), cam.camera_type[0]
            additional_params = []
            for key in ['camera_to_worlds', 'intrins', 'distortion_params']:
                value = getattr(cam, key)  # type: Optional[torch.Tensor]
                additional_params.append(value if value is None else tuple(value.shape))
            return (w, h, camera_type, *additional_params)

        self.train_index_group_loader = IndexGroupsWithDataLoader(
            self.train_indices, [get_key(idx) for idx in self.train_indices],
            self.get_train_image, self.train_batch_size(False),
            #self.config.cache_images != "gpu"
            self.config.cache_images == "disk",
            self.eval,
        )
        self.val_index_group_loader = IndexGroupsWithDataLoader(
            self.val_indices, [get_key(idx) for idx in self.val_indices],
            self.get_train_image, self.val_batch_size(False),
            #self.config.cache_images != "gpu"
            self.config.cache_images == "disk",
            self.eval,
        )

    def next_train(self, step: int, max_steps: Optional[int]) -> List[Tuple[Cameras, Dict]]:
        """Returns the next training batch

        Returns a Camera instead of raybundle"""

        if not hasattr(self, 'cached_train') and self.config.cache_images != "disk":
            self.cached_train = self._load_images("train", cache_images_device=self.config.cache_images)

        if self.train_index_group_loader is None:
            self.val_indices = sorted(self.train_dataset.val_indices)
            self.train_indices = sorted(
                set(range(len(self.train_dataset))).difference(self.train_dataset.val_indices))
            self.setup_index_group_loaders()

        def pack_batch(camera, batch, max_batch_size=None, _no_split_batch=False) -> List[Tuple[Cameras, Dict]]:
            if self.config.split_batch and not _no_split_batch:
                n = len(batch['image'])
                results = []
                for i in range(n):
                    slice_value = lambda value: value[i:i+1] if (isinstance(value, torch.Tensor) or isinstance(value, list)) and len(value) == n else value
                    camera_i = {}
                    for key, value in camera.items():
                        camera_i[key] = slice_value(value)
                    camera_i['metadata'] = {}
                    for key, value in camera.get('metadata', '').items():
                        camera_i['metadata'][key] = slice_value(value)
                    batch_i = {}
                    for key, value in batch.items():
                        batch_i[key] = slice_value(value)
                    results.extend(pack_batch(camera_i, batch_i, max_batch_size, _no_split_batch=True))
                return results

            # Engine path: keep training data on CPU; C++ does H→D in pool buffers each step.
            for key, value in camera.items():
                if isinstance(value, torch.Tensor) and value.numel() == 0:
                    camera[key] = None
                elif isinstance(value, torch.Tensor) and max_batch_size is not None:
                    camera[key] = value[:max_batch_size]
                else:
                    camera[key] = value
            metadata = camera.pop('metadata', {})
            camera = Cameras(**camera)
            for key, value in metadata.items():
                if isinstance(value, torch.Tensor) and value.numel() == 0:
                    metadata[key] = None
            camera.metadata = metadata
            return [(camera, batch)]

        camera, batch = self.train_index_group_loader.get_batch()
        if camera['camera_type'][0] == CameraType.EQUIRECTANGULAR.value:
            camera['camera_to_worlds'], camera['intrins'], camera['distortion_params'], batch['image'] = \
                warp_equirectangular_to_pinhole(camera['camera_to_worlds'].cuda(), batch['image'].cuda())
            camera['height'], camera['width'] = batch['image'].shape[-3:-1]
            if 'cam_idx' in camera.get('metadata', {}):
                camera['metadata']['cam_idx'] = torch.stack(sum([[6*i+j for j in range(6)] for i in camera['metadata']['cam_idx']], []))
            camera['camera_type'] = [CameraType.PERSPECTIVE.value] * len(camera['intrins'])
        elif self.config.warp_to_pinhole:
            if 'mask' not in batch:
                batch['mask'] = torch.ones((*batch['image'].shape[:3], 1), dtype=torch.bool, device=batch['image'].device)
            for key in ['image', 'depth', 'normal', 'mask']:
                if key not in batch:
                    continue
                elif key in ['depth', 'normal']:  # TODO: add support
                    raise NotImplementedError("Depth and normal are not supported in warp to pinhole mode. Use --datamanager.no-load-depths and --datamanager.no-load-normals if needed.")
                if key in ['depth', 'mask'] and batch[key].ndim == 3:
                    batch[key] = batch[key].unsqueeze(-1)
                returns = warp_wide_to_pinhole(camera['camera_type'][0], camera['intrins'].cuda(), camera['distortion_params'].cuda(), camera['camera_to_worlds'].cuda(), batch[key].cuda())
                batch[key] = returns[-1]
            camera['camera_to_worlds'], camera['intrins'], camera['distortion_params'] = returns[:3]
            camera['height'], camera['width'] = batch['image'].shape[-3:-1]
            if 'cam_idx' in camera.get('metadata', {}):
                camera['metadata']['cam_idx'] = torch.stack(sum([[5*i+j for j in range(5)] for i in camera['metadata']['cam_idx']], []))
            camera['camera_type'] = [CameraType.PERSPECTIVE.value] * len(camera['intrins'])
        results = pack_batch(camera, batch)

        val_batch_size = self.val_batch_size(True)
        has_val = len(self.train_dataset.val_indices) > 0 and val_batch_size > 0
        if has_val:
            val_camera, val_batch = self.val_index_group_loader.get_batch()
            val_results = pack_batch(val_camera, val_batch, val_batch_size)

        if has_val:
            return results, val_results
        return results
