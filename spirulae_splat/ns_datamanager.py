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

    cache_images: Literal["cpu-pageable", "cpu", "gpu"] = "cpu-pageable"
    """Whether to cache images in memory. If "cpu", caches on cpu. If "gpu", caches on device."""
    cache_images_type: Literal["uint8", "float32"] = "uint8"
    """The image type returned from manager, caching images in uint8 saves memory"""

    # TODO: probably do this offline and load with data
    depth_model: Literal[
        # disabled
        None,
        # MoGe - https://github.com/microsoft/moge
        # fast one with good 3D points, with some memory leak
        "moge:Ruicheng/moge-vitl",
        # Depth Pro: https://github.com/apple/ml-depth-pro
        # fine details (e.g. wires, tree branches); slower, recommend >10GB VRAM
        "depth_pro",
        # VGGT: https://github.com/facebookresearch/vggt
        "vggt:facebook/VGGT-1B",
        # Depth Anything v2 - not recommended based on empirical results
        "depth_anything_v2_vits",
        "depth_anything_v2_vitb",
        "depth_anything_v2_vitl",
        # Composite model - Combine high resolution depth and low resolution mask to produce high resolution mask
        "depth_pro + moge:Ruicheng/moge-vitl",
    ] = None  # "moge:Ruicheng/moge-vitl"



class DepthPredictor(torch.nn.Module):

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
        },
        "depth_pro": {
            "sky_threshold": np.inf,  # inconsistent
        },
        "vggt:facebook/VGGT-1B": {
            "sky_threshold": np.inf,
        },
    }
    # TODO: https://github.com/AIVFI/Monocular-Depth-Estimation-Rankings-and-2D-to-3D-Video-Conversion-Rankings/blob/main/README.md

    cache_dir: Optional[str] = os.path.join("scripts", "depth_cache")

    def __init__(self, model_id: str, use_cache: bool=True, primary_device=None):
        super().__init__()

        if model_id not in self.MODELS:
            raise ValueError(f"Model must be one of {self.MODELS.keys()}; Current: `{model_id}`")
        self.model_id = model_id
        self.sky_threshold = self.MODELS[model_id]['sky_threshold']
        self.use_cache = use_cache

        # Find all available CUDA devices
        self.device_count = torch.cuda.device_count()
        self.device0 = "cuda:0"
        if primary_device is not None:
            self.device0 = primary_device

        self._load_model()

    @staticmethod
    def _get_hash_log_filename(model_id: str):
        disallowed = r'[<>:"/\\|?*\x00-\x1F\x7F\s]'
        model_id = re.sub(disallowed, '_', model_id)
        model_id = re.sub(r'_+', '_', model_id).lstrip('_')
        return f"depth_hashes_{model_id}.json"

    @staticmethod
    def _update_hash_log(log_filename: str, hash: str, idx: str):
        try:
            with open(log_filename, 'r') as fp:
                content = json.load(fp)
        except:
            content = {}
        idx = str(idx)
        if idx not in content or content[idx] != hash:
            content[idx] = hash
            with open(log_filename, 'w') as fp:
                json.dump(content, fp, indent=4)

    @staticmethod
    def _get_cache_filename(image: torch.Tensor, model_id: str) -> str:
        # compute result file name (image hash ^ model ID hash)
        image_hash = image.flatten()
        image_hash = image_hash[:(image.numel()//16)*16].reshape(16, -1)
        image_hash = (image_hash.sum(dim=1) % 256).byte()
        image_hash = image_hash.cpu().numpy()
        model_hash = hashlib.md5(model_id.encode('utf-8')).digest()
        model_hash = np.frombuffer(model_hash, dtype=np.uint8)
        result_hash = np.bitwise_xor(model_hash, image_hash)
        result_hash = ''.join([f'{byte:02x}' for byte in result_hash])
        # locate cache file
        script_dir = os.path.dirname(os.path.realpath(__file__))
        cache_dir = os.path.join(script_dir, DepthPredictor.cache_dir, result_hash[:2])
        os.makedirs(cache_dir, exist_ok=True)
        return os.path.join(cache_dir, f"{result_hash}.npz")

    def forward(self, image: torch.Tensor, camera: Cameras=None):
        image = image.to(self.device0)

        depth = None
        depth_loaded = False

        # cache lookup
        if self.use_cache and self.cache_dir is not None:
            cache_filename = self._get_cache_filename(image, self.model_id)
            if os.path.exists(cache_filename):
                depth = np.load(cache_filename)
                depth = depth.f.depth
                depth = torch.from_numpy(depth).float().to(image.device)
                depth_loaded = True

        # inference depth model
        if depth is None:
            with torch.inference_mode():
                depth = self.infer(image, camera)
                depth = depth.to(self.device0)
            torch.cuda.empty_cache()

        # save cache
        if self.use_cache and self.cache_dir is not None and not depth_loaded:
            np.savez(cache_filename, depth=depth.cpu().numpy())

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

    def _get_model_path(self, model: str):
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

    def _load_model(self):
        self.model = None

        if self.model_id.startswith("depth_anything_v2"):
            model_path = self._get_model_path(self.model_id)
            self.model = self._load_depth_anything_v2(model_path)
            self.infer = self._infer_depth_anything_v2

        elif self.model_id.startswith("depth_pro"):
            self.model = self._load_depth_pro()
            self.infer = self._infer_depth_pro

        elif self.model_id.startswith("moge:"):
            self.model = self._load_moge()
            self.infer = self._infer_moge

        elif self.model_id.startswith("vggt:"):
            self.model = self._load_vggt()
            self.infer = self._infer_vggt

        if self.model is not None:
            self.model = self.model.to(self.device0)
            CONSOLE.log(f"\nDepth model loaded: {self.model_id}, {self.model.__class__}")

    def _load_depth_anything_v2(self, model_path):
        from depth_anything_v2.dpt import DepthAnythingV2

        model_configs = {
            'vits': {'encoder': 'vits', 'features': 64, 'out_channels': [48, 96, 192, 384]},
            'vitb': {'encoder': 'vitb', 'features': 128, 'out_channels': [96, 192, 384, 768]},
            'vitl': {'encoder': 'vitl', 'features': 256, 'out_channels': [256, 512, 1024, 1024]},
            'vitg': {'encoder': 'vitg', 'features': 384, 'out_channels': [1536, 1536, 1536, 1536]}
        }
        encoder = self.model_id.split('_')[-1]

        model = DepthAnythingV2(**model_configs[encoder])
        model.load_state_dict(torch.load(model_path, map_location='cpu'))
        model = model.to(self.device0).eval()

        return model

    def _load_depth_pro(self):
        try:
            from depth_pro.utils import load_rgb
            from depth_pro import depth_pro
        except ImportError:
            raise ImportError("Import error, please install https://github.com/apple/ml-depth-pro")

        config = depth_pro.DEFAULT_MONODEPTH_CONFIG_DICT
        config.checkpoint_uri = os.path.join(os.path.dirname(depth_pro.__file__), "../../", config.checkpoint_uri)
        # use FP32 since FP16 gives noisy normal, need much more RAM
        model, _ = depth_pro.create_model_and_transforms(config, self.device0, torch.float32)
        return model.eval()

    def _load_moge(self):
        try:
            from moge.model.v1 import MoGeModel
        except ImportError:
            raise ImportError("Import error, please install https://github.com/microsoft/moge")
        model = MoGeModel.from_pretrained(self.model_id.split(':')[-1]).to(self.device0).eval()
        return model

    def _load_vggt(self):
        try:
            from vggt.models.vggt import VGGT
        except ImportError:
            raise ImportError("Import error, please install https://github.com/facebookresearch/vggt")
        model = VGGT.from_pretrained(self.model_id.split(':')[-1]).to(self.device0).eval()
        # model = model.half()
        return model

    def _infer_depth_anything_v2(self, image, camera=None):
        max_depth_image_size: int = 1000

        batch = image.permute(2, 0, 1).unsqueeze(0)
        batch = batch.float().to(self.device0) / 255.0

        b, c, h, w = batch.shape
        sc = min(max_depth_image_size/np.sqrt(h*w), 1.0)  # save memory
        h, w = int(sc*h+0.5), int(sc*w+0.5)
        h1 = ((h-7)//14+1) * 14
        w1 = ((w-7)//14+1) * 14
        batch = torch.nn.functional.interpolate(
            batch, size=(h1, w1), mode='bilinear', align_corners=False)

        mean = torch.tensor([[[[0.485]], [[0.456]], [[0.406]]]]).float().to(self.device0)
        std = torch.tensor([[[[0.229]], [[0.224]], [[0.225]]]]).float().to(self.device0)
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

    def _infer_depth_pro(self, image, camera=None):
        image = image.permute(2, 0, 1)
        image = image.float().to(self.device0) / 255.0
        image = 2.0*image-1.0

        f_px = None
        if camera is not None:
            f_px = camera.fx.to(image.device)

        prediction = self.model.infer(image, f_px=f_px)
        depth = prediction["depth"]

        del prediction
        return depth.unsqueeze(-1)

    def _infer_moge(self, image, camera=None):
        image = image.permute(2, 0, 1)
        image = image.float().to(self.device0) / 255.0

        fov_x = None
        if camera is not None:
            fov_x = 2.0 * torch.rad2deg(torch.arctan(camera.width / (2.0*camera.fx)))

        output = self.model.infer(image, fov_x)
        depth = output['depth']
        mask = output['mask']

        del output
        depth[~mask] = 0.0
        return depth

    def _infer_vggt(self, image, camera=None):
        # TODO: do this batched, the model wasn't trained on single image
        image = image.permute(2, 0, 1)
        image = image.to(self.device0) / 255.0

        target_size = 518
        block_size = 14

        height, width = image.shape[-2:]
        if width >= height:
            new_width = target_size
            new_height = round(height * (new_width / width) / block_size) * block_size
        else:
            # seems to work, not sure if this drops model performance
            new_height = target_size
            new_width = round(width * (new_height / height) / block_size) * block_size

        batch = torch.nn.functional.interpolate(
            image.unsqueeze(0), size=(new_height, new_width),
            mode='bilinear', align_corners=False)

        aggregated_tokens_list, ps_idx = self.model.aggregator(batch[None])
        depth, _ = self.model.depth_head(aggregated_tokens_list, batch[None], ps_idx)
        depth = depth.squeeze(0).squeeze(0)
        
        del aggregated_tokens_list
        del ps_idx
        return depth


class CompositeDepthPredictor(torch.nn.Module):

    MODELS = {
        frozenset({'depth_pro', 'moge:Ruicheng/moge-vitl'}): None,
    }

    _loaded_models: Dict[str, DepthPredictor] = {}

    def __init__(self, model_id: str, dataset_path: Optional[str] = None):
        super().__init__()
        self.model_id_str = model_id
        self.model_id = frozenset(model_id.split(' + '))
        if len(self.model_id) > 1 and self.model_id not in self.MODELS:
            raise ValueError(f"Model must be one of {self.MODELS.keys()}; Current: `{model_id}`")

        self.device_count = torch.cuda.device_count()
        self.device0 = "cuda:0"

        if len(self.model_id) == 1:
            self.infer = self._infer_single_model
        elif self.model_id == frozenset({'depth_pro', 'moge:Ruicheng/moge-vitl'}):
            self.infer = self._infer_depth_pro_moge

        self._hash_log_filename = os.path.join(
            dataset_path, DepthPredictor._get_hash_log_filename(model_id))

    def _init_models(self):
        device_idx = 0
        for mid in self.model_id:
            if mid in self._loaded_models:
                continue
            device_idx = (device_idx + 1) % self.device_count
            model = DepthPredictor(
                mid, use_cache=False,
                primary_device=f"cuda:{device_idx}"
            )
            if self.device_count > 0 and False:
                model = DataParallel(model, device_ids=range(self.device_count))
            self._loaded_models[mid] = model

    def _get_model(self, mid: str):
        if mid not in self._loaded_models:
            self._init_models()
        return self._loaded_models[mid]

    def forward(self, image: torch.Tensor, camera: Cameras=None, idx: Optional[int]=None):

        depth = None
        depth_loaded = False

        # cache lookup
        if DepthPredictor.cache_dir is not None:
            cache_filename = DepthPredictor._get_cache_filename(image, self.model_id_str)
            if idx is not None:
                DepthPredictor._update_hash_log(self._hash_log_filename, cache_filename, idx)
            if os.path.exists(cache_filename):
                depth = np.load(cache_filename)
                depth = depth.f.depth
                depth = torch.from_numpy(depth).float().to(image.device)
                depth_loaded = True

        # inference depth model
        if depth is None:
            with torch.inference_mode():
                depth = self.infer(image, camera)
                depth = depth.to(self.device0)

        # save cache
        if DepthPredictor.cache_dir is not None and not depth_loaded:
            # np.savez_compressed(cache_filename, depth=depth.cpu().numpy())
            # np.savez(cache_filename, depth=depth.half().cpu().numpy())
            np.savez(cache_filename, depth=depth.cpu().numpy())

        if image.shape[:2] != depth.shape[:2]:
            h, w, _ = image.shape
            depth = torch.nn.functional.interpolate(
                depth.reshape(1, 1, *depth.shape[:2]),
                size=(h, w), mode='bilinear', align_corners=False
            ).reshape(h, w, 1)

        return depth

    def _infer_single_model(self, image, camera=None):
        depth = self._get_model(self.model_id_str)(image, camera)
        return depth.float().to(self.device0)

    def _infer_depth_pro_moge(self, image, camera=None):
        depth_sharp = self._get_model('depth_pro')(image, camera).float().to(self.device0)
        depth_blurry = self._get_model('moge:Ruicheng/moge-vitl')(image, camera).float().to(self.device0)
        mask_blurry = depth_blurry > 0

        depth_sharp_downsampled = self._downscale_image(depth_sharp)
        z1 = torch.quantile(depth_sharp_downsampled, 0.98)
        depth_sharp_clipped = torch.clip(depth_sharp, max=z1)

        # fit mask; TODO: logistic regression
        dist_model = self._fit_distortion_model(depth_sharp_clipped, mask_blurry)
        mask = dist_model(depth_sharp_clipped) > 0.75

        # fit depth
        dist_model = self._fit_distortion_model(depth_sharp_clipped, depth_blurry, mask)
        depth_undistorted = dist_model(depth_sharp_clipped)

        depth = depth_undistorted
        depth[~mask] = 0.0
        return depth.to(image.device)

    @staticmethod
    def _downscale_image(im, size=384):
        im = im.float()
        h0, w0 = im.shape[:2]
        sc = size / min(h0, w0)
        if sc >= 1:
            return im.reshape((h0, w0))
        h1, w1 = int(sc*h0+1), int(sc*w0+1)
        lq = torch.nn.functional.interpolate(
            im.reshape(1, 1, *im.shape[:2]),
            size=(h1, w1), mode='bilinear', align_corners=False
        )
        return lq.reshape((h1, w1))

    @staticmethod
    def _generate_depth_embedding(x, _uv_cache={}):
        h, w = x.shape[:2]
        x = x.flatten()

        uv_degree = 1
        z_degree = 3

        if (h, w) not in _uv_cache:
            u = (torch.arange(w, dtype=torch.float32)+0.5)/w * 2.0 - 1.0
            v = (torch.arange(h, dtype=torch.float32)+0.5)/h * 2.0 - 1.0
            u = u[None, :].float().to(x.device).repeat((h, 1))
            v = v[:, None].float().to(x.device).repeat((1, w))
            u, v = u.flatten(), v.flatten()

            uv = []
            for i in range(uv_degree+1):
                for j in range(uv_degree+1):
                    uv.append(torch.cos(math.pi/2*i*u)*torch.cos(math.pi/2*j*v))
                    if j != 0:
                        uv.append(torch.cos(math.pi/2*i*u)*torch.sin(math.pi/2*j*v))
                    if i != 0:
                        uv.append(torch.sin(math.pi/2*i*u)*torch.cos(math.pi/2*j*v))
                    if i != 0 and j != 0:
                        uv.append(torch.sin(math.pi/2*i*u)*torch.sin(math.pi/2*j*v))

            uv = torch.stack(uv)
            _uv_cache[(h, w)] = uv

        else:
            uv = _uv_cache[(h, w)].to(x.device)

        A = []
        for k in range(z_degree+1):
            ed = x**k / math.factorial(k) * torch.exp(-x)
            A.append(ed * uv)
        return torch.concatenate(A)

    @staticmethod
    def _fit_distortion_model(x, y, weights=None):
        x = CompositeDepthPredictor._downscale_image(x)
        y = CompositeDepthPredictor._downscale_image(y)
        if weights is not None:
            weights = CompositeDepthPredictor._downscale_image(weights)

        if weights is None:
            x_mean, x_std = torch.mean(x).item(), torch.std(x).item()
            x = (x - x_mean) / x_std
        else:
            weights = weights + 1.0 / weights.numel()  # prevent degeneracy
            m_sum = weights.sum().item()
            x_mean = (x*weights).sum().item() / m_sum
            x = x - x_mean
            x_std = math.sqrt(((x**2)*weights).sum().item() / m_sum)
            x = x / x_std

        A = CompositeDepthPredictor._generate_depth_embedding(x)

        # linear least squares
        if weights is None:
            mat = A @ A.T
            vec = A @ y.flatten()
        else:
            mat = (weights.reshape(1,-1) * A) @ A.T
            vec = A @ (weights.flatten() * y.flatten())
        c = torch.linalg.solve(mat, vec)

        def pred(x):
            h, w = x.shape[:2]
            x = (x - x_mean) / x_std
            A = CompositeDepthPredictor._generate_depth_embedding(x)
            return (c @ A).reshape(h, w)

        return pred


class SpirulaeDataManager(FullImageDatamanager):
    """Template DataManager

    Args:
        config: the DataManagerConfig used to instantiate class
    """

    config: SpirulaeDataManagerConfig

    _train_call_count: int = 0

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

        self.train_unseen_cameras = deque()

        self.cached_train_images = []

    def _load_images(
        self, split: Literal["train", "eval"], cache_images_device: Literal["cpu-pageable", "cpu", "gpu"]
    ) -> List[Dict[str, torch.Tensor]]:
        undistorted_images: List[Dict[str, torch.Tensor]] = []

        if self.config.depth_model is not None and split == "train":
            SpirulaeDataManager._train_call_count += 1
            if SpirulaeDataManager._train_call_count > 1 or torch.cuda.device_count() > 1:
                message = (
                    "\nWARNING: Multi-GPU with depth prediction is currently not supported,\n"
                    "and can lead to poor performance and system instability if used.\n"
                    "If multi-GPU training is intended, use single GPU on the first time training on the dataset,\n"
                    "and abort it after depth prediction and rerun with multi-GPU.\n"
                    "If not, set CUDA_VISIBLE_DEVICES to the index of the intended GPU before training.\n"
                )
                CONSOLE.print('[bold yellow]'+message)

            if ' + ' in self.config.depth_model or True:
                self.depth_model = CompositeDepthPredictor(self.config.depth_model, dataset_path=str(self.config.data))
            else:
                self.depth_model = DepthPredictor(self.config.depth_model)

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
                mask = None
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

        # predict depth
        def predict_depth(idx) -> Dict[str, torch.Tensor]:
            cache = undistorted_images[idx]
            # if idx != 0: return cache
            image = cache["image"]
            camera = dataset.cameras[idx]
            try:
                depth = self.depth_model(image, camera, idx)
            except:
                import traceback
                traceback.print_exc()
            cache["depth"] = depth
            return cache

        if self.config.depth_model is not None and split == "train":
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
        elif cache_images_device.startswith("cpu"):
            for cache in undistorted_images:
                if cache_images_device == "cpu":  # not cpu-pageable
                    cache["image"] = cache["image"].pin_memory()
                    if "mask" in cache:
                        cache["mask"] = cache["mask"].pin_memory()
                self.train_cameras = self.train_dataset.cameras
        else:
            assert_never(cache_images_device)

        if self.config.depth_model is not None and split == "train":
            del self.depth_model
            torch.cuda.empty_cache()

        return undistorted_images

    def next_train(self, step: int) -> Tuple[Cameras, Dict]:
        """Returns the next training batch

        Returns a Camera instead of raybundle"""

        if len(self.cached_train_images) == 0:
            self.cached_train_images = self._load_images("train", cache_images_device=self.config.cache_images)

        train_batch_size = (len(self.train_dataset) + self.config.max_batch_per_epoch - 1) \
            // self.config.max_batch_per_epoch
        # train_batch_size = 4

        image_indices = []
        for i in range(train_batch_size):
            if len(self.train_unseen_cameras) == 0:
                perm = list(range(len(self.train_dataset)))
                random.shuffle(perm)
                self.train_unseen_cameras.extend(perm)
            image_idx = self.train_unseen_cameras.popleft()
            image_indices.append(image_idx)

        batch = {}
        for image_idx in image_indices:
            img = self.cached_train_images[image_idx]
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

