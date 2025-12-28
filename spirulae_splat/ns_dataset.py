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


from nerfstudio.data.datasets.base_dataset import *
from typing import Callable, Optional
import cv2


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


def compute_overexposure_mask(img: torch.Tensor, image_type: Literal['uint8', 'float32']):
    """Compute over exposure mask, should suffice for most unedited real-world photos"""
    import cv2
    img = img.numpy()
    if image_type == "uint8":
        img = img.astype(np.float32) / 255.0
    gray = 0.2126 * img[...,0] + 0.7152 * img[...,1] + 0.0722 * img[...,2]
    mask = (gray > 0.97) | ((img[...,0] > 0.98) & (img[...,1] > 0.98) & (img[...,2] > 0.98))
    mask = cv2.GaussianBlur(mask.astype(np.float32), (5,5), 0)
    mask = (mask < 0.5).astype(np.bool_)
    return torch.from_numpy(mask).unsqueeze(-1)


from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm

import math

from spirulae_splat.splat.utils import resize_image



class SpirulaeDataset(InputDataset):
    """Dataset that returns images.

    Args:
        dataparser_outputs: description of where and how to read input images.
        scale_factor: The scaling factor for the dataparser outputs
    """

    exclude_batch_keys_from_device: List[str] = ["image", "mask"]
    cameras: Cameras

    def __init__(self, dataparser_outputs: DataparserOutputs, scale_factor: float = 1.0, *args, **kwargs):
        super().__init__(dataparser_outputs, scale_factor, *args, **kwargs)
        self._dataparser_outputs = dataparser_outputs
        self.scale_factor = scale_factor
        self.scene_box = deepcopy(dataparser_outputs.scene_box)
        self.metadata = deepcopy(dataparser_outputs.metadata)
        self.cameras = deepcopy(dataparser_outputs.cameras)
        self.cameras.rescale_output_resolution(scaling_factor=scale_factor)
        self.mask_color = dataparser_outputs.metadata.get("mask_color", None)
        self.mask_overexposure = dataparser_outputs.metadata.get("mask_overexposure", False)

    def __len__(self):
        return len(self._dataparser_outputs.image_filenames)

    def get_numpy_image(self, image_idx: int) -> npt.NDArray[np.uint8]:
        """Returns the image of shape (H, W, 3 or 4).

        Args:
            image_idx: The image index in the dataset.
        """
        image_filename = self._dataparser_outputs.image_filenames[image_idx]
        pil_image = Image.open(image_filename)
        if self.scale_factor != 1.0:
            width, height = pil_image.size
            newsize = (int(width * self.scale_factor), int(height * self.scale_factor))
            pil_image = pil_image.resize(newsize, resample=Image.Resampling.BILINEAR)
        image = np.array(pil_image, dtype="uint8")  # shape is (h, w) or (h, w, 3 or 4)
        if len(image.shape) == 2:
            image = image[:, :, None].repeat(3, axis=2)
        assert len(image.shape) == 3
        assert image.dtype == np.uint8
        assert image.shape[2] in [3, 4], f"Image shape of {image.shape} is in correct."
        return image

    def get_image_float32(self, image_idx: int) -> Float[Tensor, "image_height image_width num_channels"]:
        """Returns a 3 channel image in float32 torch.Tensor.

        Args:
            image_idx: The image index in the dataset.
        """
        image = torch.from_numpy(self.get_numpy_image(image_idx).astype("float32") / 255.0)
        if self._dataparser_outputs.alpha_color is not None and image.shape[-1] == 4:
            assert (self._dataparser_outputs.alpha_color >= 0).all() and (
                self._dataparser_outputs.alpha_color <= 1
            ).all(), "alpha color given is out of range between [0, 1]."
            image = image[:, :, :3] * image[:, :, -1:] + self._dataparser_outputs.alpha_color * (1.0 - image[:, :, -1:])
        return image

    def get_image_uint8(self, image_idx: int) -> UInt8[Tensor, "image_height image_width num_channels"]:
        """Returns a 3 channel image in uint8 torch.Tensor.

        Args:
            image_idx: The image index in the dataset.
        """
        image = torch.from_numpy(self.get_numpy_image(image_idx))
        if self._dataparser_outputs.alpha_color is not None and image.shape[-1] == 4:
            assert (self._dataparser_outputs.alpha_color >= 0).all() and (
                self._dataparser_outputs.alpha_color <= 1
            ).all(), "alpha color given is out of range between [0, 1]."
            image = image[:, :, :3] * (image[:, :, -1:] / 255.0) + 255.0 * self._dataparser_outputs.alpha_color * (
                1.0 - image[:, :, -1:] / 255.0
            )
            image = torch.clamp(image, min=0, max=255).to(torch.uint8)
        return image

    def get_data(self, image_idx: int, image_type: Literal["uint8", "float32"] = "float32", _is_viewer=True, _load_auxiliary=True, _viewer_thumbnail_cache={}) -> Dict:
        """Returns the ImageDataset data as a dictionary.

        Args:
            image_idx: The image index in the dataset.
            image_type: the type of images returned
        """
        # multithreaded image loading, speed up nerfstudio viewer.py without having to modify nerfstudio
        # data manager will explicitly pass _is_viewer=False
        if _is_viewer:
            if len(_viewer_thumbnail_cache) == 0:
                _viewer_thumbnail_cache_1 = []
                def load_data(idx):
                    data = self.get_data(idx, image_type, _is_viewer=False, _load_auxiliary=False)
                    if 'mask' in data:
                        # background = torch.ones_like(data['image']) * [0.125, 32][image_type == 'uint8']
                        background = torch.ones_like(data['image']) * torch.tensor([(0,0,1), (0,0,255)][image_type == 'uint8']).to(data['image'])
                        data['image'] = torch.where(data['mask'], data['image'], background)
                    data['image'] = resize_image(data['image'][None], 2**int(max(0.5*math.log2(data['image'].numel()/10000), 0.0)))[0]
                    return {
                        'image_idx': data['image_idx'],
                        'image': data['image'],
                    }
                with ThreadPoolExecutor() as executor:
                    for result in tqdm(
                        executor.map(load_data, range(len(self))),
                        total=len(self),
                        desc="Loading images"
                    ):
                        _viewer_thumbnail_cache_1.append(result)
                for i, im in enumerate(_viewer_thumbnail_cache_1):
                    _viewer_thumbnail_cache[i] = im
            if image_idx in _viewer_thumbnail_cache:
                # print('cache hit:', image_idx, image_type)
                data = _viewer_thumbnail_cache[image_idx]
                del _viewer_thumbnail_cache[image_idx]
                return data
            # print('cache miss:', image_idx, image_type)
        else:
            del _viewer_thumbnail_cache

        # regular loading
        if image_type == "float32":
            image = self.get_image_float32(image_idx)
        elif image_type == "uint8":
            image = self.get_image_uint8(image_idx)
        else:
            raise NotImplementedError(f"image_type (={image_type}) getter was not implemented, use uint8 or float32")

        data = {"image_idx": image_idx, "image": image}
        if self._dataparser_outputs.mask_filenames is not None:
            mask_filepath = self._dataparser_outputs.mask_filenames[image_idx]
            data["mask"] = get_image_mask_tensor_from_path(
                filepath=mask_filepath,
                width=data["image"].shape[1], height=data["image"].shape[0]
            )
            assert (
                data["mask"].shape[:2] == data["image"].shape[:2]
            ), f"Mask and image have different shapes. Got {data['mask'].shape[:2]} and {data['image'].shape[:2]}"
        if _load_auxiliary:
            if self._dataparser_outputs.metadata.get("depth_filenames", None) is not None:
                depth_filepath = self._dataparser_outputs.metadata["depth_filenames"][image_idx]
                data["depth"] = get_depth_image_from_path(
                    filepath=depth_filepath, scale_factor=self.scale_factor,
                    width=data["image"].shape[1], height=data["image"].shape[0]
                )
                # assert data["depth"].shape[:2] == data["image"].shape[:2]
            if self._dataparser_outputs.metadata.get("normal_filenames", None) is not None:
                normal_filepath = self._dataparser_outputs.metadata["normal_filenames"][image_idx]
                data["normal"] = get_normal_image_from_path(
                    filepath=normal_filepath,
                    width=data["image"].shape[1], height=data["image"].shape[0]
                )
                # assert data["normal"].shape[:2] == data["image"].shape[:2]
        if self.mask_color:
            data["image"] = torch.where(
                data["mask"] == 1.0, data["image"], torch.ones_like(data["image"]) * torch.tensor(self.mask_color)
            )
        if self.mask_overexposure:
            o_mask = compute_overexposure_mask(data["image"], image_type)
            if "mask" not in data:
                data['mask'] = o_mask
            else:
                data['mask'] = data['mask'] & o_mask

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

        return self._dataparser_outputs.image_filenames


class IndexedDatasetWrapper(Dataset):

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
