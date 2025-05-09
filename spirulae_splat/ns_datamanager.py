"""
Template DataManager
"""

from dataclasses import dataclass, field
from typing import Dict, Literal, List, Tuple, Type, Union, Optional
from copy import deepcopy
import random
import os
import hashlib

import torch
import numpy as np

from nerfstudio.cameras.rays import RayBundle
from nerfstudio.data.datamanagers.full_images_datamanager import (
    DataManagerConfig,
    FullImageDatamanager,
    FullImageDatamanagerConfig,
    _undistort_image,
    CONSOLE, track
)
from nerfstudio.cameras.cameras import Cameras, CameraType

from spirulae_splat.splat._torch_impl import depth_map, depth_inv_map


from concurrent.futures import ThreadPoolExecutor


@dataclass
class SpirulaeDataManagerConfig(FullImageDatamanagerConfig):
    """Template DataManager Config

    Add your custom datamanager config parameters here.
    """
    _target: Type = field(default_factory=lambda: SpirulaeDataManager)

    max_batch_per_epoch: int = 768
    """Maximum number of batches per epoch, used for configuring batch size"""

    cache_images: Literal["cpu", "gpu"] = "gpu"
    """Whether to cache images in memory. If "cpu", caches on cpu. If "gpu", caches on device."""
    cache_images_type: Literal["uint8", "float32"] = "float32"
    """The image type returned from manager, caching images in uint8 saves memory"""

    depth_model: Literal[
        # disabled
        None,
        # MoGe - https://github.com/microsoft/moge
        "moge:Ruicheng/moge-vitl",
        # Depth Anything v2 - not recommended based on empirical results
        "depth_anything_v2_vits",
        "depth_anything_v2_vitb",
        "depth_anything_v2_vitl",
    ] = None  # "moge:Ruicheng/moge-vitl"



class DepthPredictor:

    MODELS = {
        "depth_anything_v2_vits": {
            "download_url": "https://huggingface.co/depth-anything/Depth-Anything-V2-Small/resolve/main/depth_anything_v2_vits.pth?download=true",
            "filename": "depth_anything_v2_vits.pth",
            "sky_threshold": 10.0,
        },
        "depth_anything_v2_vitb": {
            "download_url": "https://huggingface.co/depth-anything/Depth-Anything-V2-Base/resolve/main/depth_anything_v2_vitb.pth?download=true",
            "filename": "depth_anything_v2_vitb.pth",
            "sky_threshold": 10.0,
        },
        "depth_anything_v2_vitl": {
            "download_url": "https://huggingface.co/depth-anything/Depth-Anything-V2-Large/resolve/main/depth_anything_v2_vitl.pth?download=true",
            "filename": "depth_anything_v2_vitl.pth",
            "sky_threshold": 0.1,
        },
        "moge:Ruicheng/moge-vitl": {
            "sky_threshold": np.inf,
        }
    }

    cache_dir: Optional[str] = os.path.join("scripts", "depth_cache")

    max_depth_image_size: int = 1000

    def __init__(self, model_id: str):
        if model_id not in self.MODELS:
            raise ValueError(f"Model must be one of {self.MODELS.keys()}; Current: `{model_id}`")
        self.model_id = model_id
        self.sky_threshold = self.MODELS[model_id]['sky_threshold']

        self.load_model()
        CONSOLE.log(f"Depth model loaded: {model_id}")

    def __call__(self, image: torch.Tensor, camera: Cameras=None):

        depth = None
        depth_loaded = False

        # cache lookup
        if self.cache_dir is not None:
            # compute result file name (image hash ^ model ID hash)
            image_hash = image.flatten()
            image_hash = image_hash[:(image.numel()//16)*16].reshape(16, -1)
            image_hash = (image_hash.sum(dim=1) % 256).byte()
            image_hash = image_hash.cpu().numpy()
            model_hash = hashlib.md5(self.model_id.encode('utf-8')).digest()
            model_hash = np.frombuffer(model_hash, dtype=np.uint8)
            result_hash = np.bitwise_xor(model_hash, image_hash)
            result_hash = ''.join([f'{byte:02x}' for byte in result_hash])
            # locate cache file
            script_dir = os.path.dirname(os.path.realpath(__file__))
            cache_dir = os.path.join(script_dir, self.cache_dir, result_hash[:2])
            os.makedirs(cache_dir, exist_ok=True)
            cache_filename = os.path.join(cache_dir, f"{result_hash}.npz")
            if os.path.exists(cache_filename):
                depth = np.load(cache_filename)
                depth = depth.f.depth
                depth = torch.from_numpy(depth).float().to(image.device)
                depth_loaded = True

        # inference depth model
        if depth is None:
            with torch.inference_mode():
                depth = self.infer(image, camera)

        # save cache
        if self.cache_dir is not None and not depth_loaded:
            # np.savez_compressed(cache_filename, depth=depth.cpu().numpy())
            np.savez(cache_filename, depth=depth.half().cpu().numpy())  # 1.7MB -> 2.1MB, magnitudes faster loading

        if image.shape[:2] != depth.shape[:2]:
            h, w, _ = image.shape
            depth = torch.nn.functional.interpolate(
                depth.reshape(1, 1, *depth.shape[:2]),
                size=(h, w), mode='bilinear', align_corners=False
            ).reshape(h, w, 1)

        if np.isfinite(self.sky_threshold) and self.sky_threshold >= 0.0:
            depth[depth >= self.sky_threshold] = 0.0

        # return depth_map(depth)
        return depth

    def get_model_path(self, model: str):
        model = self.MODELS[model]
        download_url = model["download_url"]
        filename = model["filename"]

        script_dir = os.path.dirname(os.path.realpath(__file__))
        model_path = os.path.join(script_dir, "scripts", "models", filename)

        if not os.path.exists(model_path):
            CONSOLE.log(f"Downloading depth model {model} from {download_url}")
            import urllib.request
            urllib.request.urlretrieve(download_url, model_path)

        return model_path

    def load_model(self):

        if self.model_id.startswith("depth_anything_v2"):
            model_path = self.get_model_path(self.model_id)
            self.model = self.load_depth_anything_v2(model_path)
            self.infer = self.inference_depth_anything_v2

        elif self.model_id.startswith("moge:"):
            self.model = self.load_moge()
            self.infer = self.inference_moge

    def load_depth_anything_v2(self, model_path):
        from spirulae_splat.scripts.depth_anything_v2.dpt import DepthAnythingV2

        model_configs = {
            'vits': {'encoder': 'vits', 'features': 64, 'out_channels': [48, 96, 192, 384]},
            'vitb': {'encoder': 'vitb', 'features': 128, 'out_channels': [96, 192, 384, 768]},
            'vitl': {'encoder': 'vitl', 'features': 256, 'out_channels': [256, 512, 1024, 1024]},
            'vitg': {'encoder': 'vitg', 'features': 384, 'out_channels': [1536, 1536, 1536, 1536]}
        }
        encoder = self.model_id.split('_')[-1]

        model = DepthAnythingV2(**model_configs[encoder])
        model.load_state_dict(torch.load(model_path, map_location='cpu'))
        model = model.cuda().eval()

        return model

    def load_moge(self):
        try:
            from moge.model.v1 import MoGeModel
        except ImportError:
            raise ImportError("Import error, please install https://github.com/microsoft/moge")

        model = MoGeModel.from_pretrained(self.model_id.split(':')[-1]).cuda().eval()
        return model

    def inference_depth_anything_v2(self, image, camera=None):

        batch = image.permute(2, 0, 1).unsqueeze(0)
        batch = batch.float().cuda() / 255.0

        b, c, h, w = batch.shape
        sc = min(self.max_depth_image_size/np.sqrt(h*w), 1.0)  # save memory
        h, w = int(sc*h+0.5), int(sc*w+0.5)
        h1 = ((h-7)//14+1) * 14
        w1 = ((w-7)//14+1) * 14
        batch = torch.nn.functional.interpolate(
            batch, size=(h1, w1), mode='bilinear', align_corners=False)

        mean = torch.tensor([[[[0.485]], [[0.456]], [[0.406]]]]).float().cuda()
        std = torch.tensor([[[[0.229]], [[0.224]], [[0.225]]]]).float().cuda()
        batch = (batch - mean) / std

        # from time import perf_counter
        # torch.cuda.synchronize(); time0 = perf_counter()
        batch = self.model(batch).unsqueeze(1)
        # torch.cuda.synchronize(); time1 = perf_counter()
        # print(time1-time0)

        depth = batch[0].permute(1, 2, 0)
        depth = 1.0 / torch.clip(depth+0.1, 0.01)
        depth = depth.to(image.device)
        return depth

    def inference_moge(self, image, camera=None):

        image = image.permute(2, 0, 1)
        image = image.float().cuda() / 255.0

        fov_x = None
        if camera is not None:
            fov_x = 2.0 * torch.rad2deg(torch.arctan(camera.width / (2.0*camera.fx)))

        output = self.model.infer(image, fov_x)
        depth = output['depth']
        mask = output['mask']

        del output
        torch.cuda.empty_cache()

        depth[~mask] = 0.0
        return depth


class SpirulaeDataManager(FullImageDatamanager):
    """Template DataManager

    Args:
        config: the DataManagerConfig used to instantiate class
    """

    config: DataManagerConfig = field(default_factory=SpirulaeDataManagerConfig)

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
        self.num_train = 0
        self.num_eval = 0

    def _load_images(
        self, split: Literal["train", "eval"], cache_images_device: Literal["cpu", "gpu"]
    ) -> List[Dict[str, torch.Tensor]]:
        undistorted_images: List[Dict[str, torch.Tensor]] = []

        if self.config.depth_model is not None:
            self.depth_model = DepthPredictor(self.config.depth_model)

        # Which dataset?
        if split == "train":
            dataset = self.train_dataset
            self.num_train = len(dataset)
        elif split == "eval":
            dataset = self.eval_dataset
            self.num_eval = len(dataset)
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

            if camera.camera_type.item() == CameraType.FISHEYE.value:
                # don't undistort
                mask = None
            else:
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

        # predict depth
        def predict_depth(idx) -> Dict[str, torch.Tensor]:
            cache = undistorted_images[idx]
            # if idx != 0: return cache
            image = cache["image"]
            camera = dataset.cameras[idx]
            depth = self.depth_model(image, camera)
            cache["depth"] = depth
            return cache

        if self.config.depth_model is not None:
            CONSOLE.log(f"Predicting/loading depth for {split} images")
            undistorted_images = list(track(
                map(predict_depth, range(len(dataset))),
                description=f"Predicting depth for {split} images",
                transient=True,
                total=len(dataset),
            ))
            # for cache in undistorted_images:
            #     import matplotlib.pyplot as plt
            #     fig, (ax1, ax2) = plt.subplots(1, 2)
            #     ax1.imshow(cache['image'].cpu().numpy())
            #     ax2.imshow(depth_map(cache['depth']).cpu().numpy())
            #     plt.show()
            #     # break

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

        if self.config.depth_model is not None:
            del self.depth_model
            torch.cuda.empty_cache()

        return undistorted_images

    def next_train(self, step: int) -> Tuple[Cameras, Dict]:
        """Returns the next training batch

        Returns a Camera instead of raybundle"""

        train_batch_size = (self.num_train + self.config.max_batch_per_epoch - 1) \
            // self.config.max_batch_per_epoch
        if train_batch_size > 1 and step >= 0:
            cameras, batches = [], []
            for i in range(train_batch_size):
                camera, batch = self.next_train(-1)
                cameras.append(camera)
                batches.append(batch)
            return cameras, batches

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
