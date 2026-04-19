# Copyright 2022 the Regents of the University of California, Nerfstudio Team and contributors. All rights reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


import torch
from torch import Tensor
import numpy as np
from pathlib import Path
from typing import Callable, Optional, Literal, List, Tuple, Dict
from jaxtyping import Float, UInt8
import cv2
from PIL import Image

from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm

import math

from spirulae_splat.splat.utils import resize_image
from spirulae_splat.modules.camera import Cameras
from copy import deepcopy


def get_image_mask_tensor_from_path(filepath: Path, height: int, width: int) -> torch.Tensor:
    """
    Utility function to read a mask image from the given path and return a boolean tensor
    Modified to handle masks that are apparently grayscale but actually RGB
    """
    pil_mask = Image.open(filepath).convert("L")
    if pil_mask.size[0] != width or pil_mask.size[1] != height:
        pil_mask = pil_mask.resize((width, height), resample=Image.Resampling.NEAREST)
    mask_tensor = torch.from_numpy(np.array(pil_mask)).unsqueeze(-1).bool()
    # TODO: configurable all nonzero vs 128/0.5?
    return mask_tensor


def get_image_from_path(
    filename: Path,
    scale_factor: float,
):
    filename = Path(filename)
    assert filename.exists(), f"File `{filename}` does not exist"
    is_exr = filename.suffix.lower() == ".exr"
    is_dng = filename.suffix.lower() == ".dng"
    if is_exr:
        try:
            import OpenEXR
        except ImportError:
            print("OpenEXR not properly installed. Please install using `pip install OpenEXR`.")
            exit(0)
        exr_file = OpenEXR.File(str(filename))
        exr_channels = exr_file.channels()
        if 'RGBA' in exr_channels:
            image = exr_channels['RGBA'].pixels
        elif 'RGB' in exr_channels:
            image = exr_channels['RGB'].pixels
        else:
            raise ValueError("Unsupported EXR file. Make sure it contains channel `RGB` or `RGBA`.")
    elif is_dng:
        try:
            import rawpy
        except ImportError:
            print("rawpy not properly installed (required for DNG). Please install using `pip install rawpy`.")
            exit(0)
        with rawpy.imread(str(filename)) as raw:
            # matching https://github.com/GNOME/shotwell/blob/b7c5957b4400664f255a6a8c8a63daf2d72958c6/src/photos/GRaw.vala
            image = raw.postprocess(
                half_size=False,
                no_auto_bright=True,
                bright=1.0,
                auto_bright_thr=0.01,
                use_camera_wb=True,
                use_auto_wb=True,
                # output_color=rawpy.ColorSpace.ACES,  # ACES2065-1
                # gamma=(1.0, 1.0),
                # output_bps=16,
                output_color=rawpy.ColorSpace.sRGB,  # ACES2065-1
                output_bps=8,
                highlight_mode=rawpy.HighlightMode.Clip,
                demosaic_algorithm=rawpy.DemosaicAlgorithm.PPG,
            )
    else:
        image = cv2.imread(filename, cv2.IMREAD_UNCHANGED)
    # image = cv2.cvtColor(cv2.imread(image_filename), cv2.COLOR_BGR2RGB)
    if scale_factor != 1.0:
        width, height, _ = image.shape
        newsize = (int(width * scale_factor), int(height * scale_factor))
        image = cv2.resize(image, newsize, cv2.INTER_LINEAR)
    if len(image.shape) == 2:
        image = image[:, :, None]
    if image.shape[2] == 1:
        image = np.repeat(image, 3, axis=2)
    assert len(image.shape) == 3
    if not (is_exr or is_dng):
        image = cv2.cvtColor(image, cv2.COLOR_BGR2RGB if image.shape[2] == 3 else cv2.COLOR_BGRA2RGBA)
    if image.dtype == np.bool_:
        image = image.astype(np.uint8) * 255
    return image


def get_depth_image_from_path(
    filepath: Path,
    height: int,
    width: int,
    scale_factor: float,
    interpolation: int = cv2.INTER_NEAREST,
) -> torch.Tensor:
    """Loads, rescales and resizes depth images.
    Filepath points to a 16-bit or 32-bit depth image, or a numpy array `*.npy`.

    Args:
        filepath: Path to depth image.
        height: Target depth image height.
        width: Target depth image width.
        scale_factor: Factor by which to scale depth image.
        interpolation: Depth value interpolation for resizing.

    Returns:
        Depth image torch tensor with shape [height, width, 1].
    """
    filepath = Path(filepath)
    assert filepath.exists(), f"File `{filepath}` does not exist"
    if filepath.suffix == ".npy":
        image = np.load(filepath).astype(np.float32) * scale_factor
        # image = cv2.resize(image, (width, height), interpolation=interpolation)
    else:
        image = cv2.imread(str(filepath.absolute()), cv2.IMREAD_ANYDEPTH)
        image = image.astype(np.float32) * scale_factor
        # image = cv2.resize(image, (width, height), interpolation=interpolation)
    return torch.from_numpy(image[:, :, np.newaxis])


def get_normal_image_from_path(
    filepath: Path,
    height: int,
    width: int,
    interpolation: int = cv2.INTER_NEAREST,
) -> torch.Tensor:
    filepath = Path(filepath)
    assert filepath.exists(), f"File `{filepath}` does not exist"
    if filepath.suffix == ".npy":
        image = np.load(filepath).astype(np.float32)
        # image = cv2.resize(image, (width, height), interpolation=interpolation)
    else:
        image = cv2.cvtColor(cv2.imread(str(filepath.absolute())), cv2.COLOR_BGR2RGB)
        # image = (image.astype(np.float32) / 255.0) * 2.0 - 1.0
        # image = cv2.resize(image, (width, height), interpolation=interpolation)
    image = torch.from_numpy(image[:, :, :3])
    return image
    # return image / torch.norm(image, dim=-1, keepdim=True).clip(min=1e-12)



THUMBNAIL_RESOLUTION = 64


class SpirulaeSplatDataset:
    """Dataset that returns images.

    Args:
        dataparser_outputs: description of where and how to read input images.
        scale_factor: The scaling factor for the dataparser outputs
    """

    exclude_batch_keys_from_device: List[str] = ["image", "mask"]
    cameras: Cameras

    def __init__(self, dataparser_outputs: dict, scale_factor: float = 1.0, *args, **kwargs):
        self._dataparser_outputs = dataparser_outputs
        self.scale_factor = scale_factor
        self.metadata = deepcopy(dataparser_outputs['metadata'])
        self.cameras = deepcopy(dataparser_outputs['cameras'])
        self.cameras.rescale(scale_factor)
        self.mask_color = dataparser_outputs['metadata'].get("mask_color", None)
        self.val_indices = set(dataparser_outputs['metadata'].get("val_indices", []))

        self.thumbnails = torch.empty(
            (len(self.cameras), THUMBNAIL_RESOLUTION, THUMBNAIL_RESOLUTION, 4),
            dtype=torch.uint8
        )
        self.thumbnails[..., :3] = 128
        self.thumbnails[..., 3:] = 255
        self._thumbnail_loaded = [False] * len(self.cameras)
        self._mask_thumbnail_loaded = [False] * len(self.cameras)

    def __len__(self):
        return len(self._dataparser_outputs['image_filenames'])

    def get_numpy_image(self, image_idx: int):
        """Returns the image of shape (H, W, 3 or 4).

        Args:
            image_idx: The image index in the dataset.
        """
        image_filename = self._dataparser_outputs['image_filenames'][image_idx]
        image = get_image_from_path(image_filename, self.scale_factor)
        assert image.dtype in [np.uint8, np.uint16, np.float16, np.float32], f"Unsupported image dtype {image.dtype} for image {image_filename}"
        assert image.shape[2] in [3, 4], f"Image shape of {image.shape} is in correct."
        return image

    def get_image(self, image_idx: int) -> UInt8[Tensor, "image_height image_width num_channels"]:
        """Returns a 3 channel image in torch.Tensor, in original format (uint8, uint16, float16, float32).

        Args:
            image_idx: The image index in the dataset.
        """
        image = self.get_numpy_image(image_idx)

        # TODO: this doesn't work in PyTorch dataloader in disk caching mode
        if not self._thumbnail_loaded[image_idx]:
            thumbnail = image[:, :, :3]  # TODO: RGBA
            if thumbnail.dtype == np.uint16:
                thumbnail = thumbnail >> 16
            elif thumbnail.dtype in [np.float16, np.float32]:
                thumbnail = (np.clip(thumbnail.astype(np.float32) / 255.0, 0.0, 1.0) * 255).astype(np.uint8)
            thumbnail = cv2.resize(thumbnail, (THUMBNAIL_RESOLUTION, THUMBNAIL_RESOLUTION))
            self.thumbnails[image_idx, :, :, :3] = torch.from_numpy(thumbnail)
            self._thumbnail_loaded[image_idx] = True

        if image.dtype == np.uint16:
            try:
                # test if torch supports uint16, which was added in PyTorch 2.3
                torch.tensor([0], dtype=torch.uint16)
                image_dtype = torch.uint16
            except:
                raise NotImplementedError("16-bit images are not supported in this version of PyTorch. Please upgrade to PyTorch 2.3 or later, or convert your images to 8-bit.")
        image = torch.from_numpy(image)
        return image

    def get_data(self, image_idx: int, load_depths=True, load_normals=True) -> Dict:
        """Returns the ImageDataset data as a dictionary.

        Args:
            image_idx: The image index in the dataset.
            image_type: the type of images returned
        """
        image = self.get_image(image_idx)

        data = {"image_idx": image_idx, "image": image}
        if self._dataparser_outputs['mask_filenames'] is not None:
            mask_filepath = self._dataparser_outputs['mask_filenames'][image_idx]
            data["mask"] = get_image_mask_tensor_from_path(
                filepath=mask_filepath,
                width=data["image"].shape[1], height=data["image"].shape[0]
            )
            assert (
                data["mask"].shape[:2] == data["image"].shape[:2]
            ), f"Mask and image have different shapes. Got {data['mask'].shape[:2]} and {data['image'].shape[:2]}"
            if not self._mask_thumbnail_loaded[image_idx]:
                thumbnail = data["mask"].numpy().astype(np.uint8) * 255
                if thumbnail.ndim == 3:
                    thumbnail = thumbnail[:, :, 0]
                thumbnail = cv2.resize(thumbnail, (THUMBNAIL_RESOLUTION, THUMBNAIL_RESOLUTION))
                self.thumbnails[image_idx, :, :, 3] = torch.from_numpy(thumbnail)
                self._mask_thumbnail_loaded[image_idx] = True
        if load_depths:
            if self._dataparser_outputs['metadata'].get("depth_filenames", None) is not None:
                depth_filepath = self._dataparser_outputs['metadata']["depth_filenames"][image_idx]
                data["depth"] = get_depth_image_from_path(
                    filepath=depth_filepath, scale_factor=self.scale_factor,
                    width=data["image"].shape[1], height=data["image"].shape[0]
                )
                # assert data["depth"].shape[:2] == data["image"].shape[:2]
        if load_normals:
            if self._dataparser_outputs['metadata'].get("normal_filenames", None) is not None:
                normal_filepath = self._dataparser_outputs['metadata']["normal_filenames"][image_idx]
                data["normal"] = get_normal_image_from_path(
                    filepath=normal_filepath,
                    width=data["image"].shape[1], height=data["image"].shape[0]
                )
                # assert data["normal"].shape[:2] == data["image"].shape[:2]
        if self.mask_color:
            data["image"] = torch.where(
                data["mask"] == 1.0, data["image"], torch.ones_like(data["image"]) * torch.tensor(self.mask_color)
            )

        metadata = self.get_metadata(data)
        data.update(metadata)
        return data

    def get_metadata(self, data: Dict) -> Dict:
        """Method that can be used to process any additional metadata that may be part of the model inputs.

        Args:
            image_idx: The image index in the dataset.
        """
        del data
        return {}

    def __getitem__(self, image_idx: int) -> Dict:
        return self.get_data(image_idx)

    @property
    def image_filenames(self) -> List[Path]:
        """
        Returns image filenames for this dataset.
        The order of filenames is the same as in the Cameras object for easy mapping.
        """

        return self._dataparser_outputs['image_filenames']


class IndexedDatasetWrapper(torch.utils.data.Dataset):

    def __init__(
            self,
            getitem: Callable,
            getitem_args: List,
            *args, **kwargs
        ):
        super().__init__(*args, **kwargs)
        self.getitem = getitem
        self.getitem_args = getitem_args

    def __len__(self):
        return len(self.getitem_args)

    def __getitem__(self, idx):
        return self.getitem(*self.getitem_args[idx])
