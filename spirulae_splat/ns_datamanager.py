"""
Template DataManager
"""

from dataclasses import dataclass, field
from typing import Dict, Literal, List, Tuple, Type, Union
from copy import deepcopy
import random

import torch

from nerfstudio.cameras.rays import RayBundle
from nerfstudio.data.datamanagers.full_images_datamanager import (
    FullImageDatamanager,
    FullImageDatamanagerConfig,
    _undistort_image,
    CONSOLE, track
)
from nerfstudio.cameras.cameras import Cameras, CameraType

from concurrent.futures import ThreadPoolExecutor


@dataclass
class SpirulaeDataManagerConfig(FullImageDatamanagerConfig):
    """Template DataManager Config

    Add your custom datamanager config parameters here.
    """
    _target: Type = field(default_factory=lambda: SpirulaeDataManager)

    cache_images: Literal["cpu", "gpu"] = "gpu"
    """Whether to cache images in memory. If "cpu", caches on cpu. If "gpu", caches on device."""
    cache_images_type: Literal["uint8", "float32"] = "float32"
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

    def _load_images(
        self, split: Literal["train", "eval"], cache_images_device: Literal["cpu", "gpu"]
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
            data = dataset.get_data(idx, image_type=self.config.cache_images_type)
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
            image = data["image"].numpy()

            K, image, mask = _undistort_image(camera, distortion_params, data, image, K)
            data["image"] = torch.from_numpy(image)
            if mask is not None:
                data["mask"] = mask

            dataset.cameras.fx[idx] = float(K[0, 0])
            dataset.cameras.fy[idx] = float(K[1, 1])
            dataset.cameras.cx[idx] = float(K[0, 2])
            dataset.cameras.cy[idx] = float(K[1, 2])
            dataset.cameras.width[idx] = image.shape[1]
            dataset.cameras.height[idx] = image.shape[0]
            return data

        CONSOLE.log(f"Caching / undistorting {split} images")
        with ThreadPoolExecutor() as executor:
            undistorted_images = list(
                track(
                    executor.map(
                        undistort_idx,
                        range(len(dataset)),
                    ),
                    description=f"Caching / undistorting {split} images",
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
                self.train_cameras = self.train_dataset.cameras.to(self.device)
        elif cache_images_device == "cpu":
            for cache in undistorted_images:
                cache["image"] = cache["image"].pin_memory()
                if "mask" in cache:
                    cache["mask"] = cache["mask"].pin_memory()
                self.train_cameras = self.train_dataset.cameras
        else:
            assert_never(cache_images_device)

        return undistorted_images

    def next_train(self, step: int) -> Tuple[Cameras, Dict]:
        """Returns the next training batch

        Returns a Camera instead of raybundle"""
        image_idx = self.train_unseen_cameras.pop(random.randint(0, len(self.train_unseen_cameras) - 1))
        # Make sure to re-populate the unseen cameras list if we have exhausted it
        if len(self.train_unseen_cameras) == 0:
            self.train_unseen_cameras = [i for i in range(len(self.train_dataset))]

        data = deepcopy(self.cached_train[image_idx])
        data["image"] = data["image"].to(self.device)

        assert len(self.train_dataset.cameras.shape) == 1, "Assumes single batch dimension"
        camera = self.train_dataset.cameras[image_idx : image_idx + 1].to(self.device)
        if camera.metadata is None:
            camera.metadata = {}
        camera.metadata["cam_idx"] = image_idx
        return camera, data
