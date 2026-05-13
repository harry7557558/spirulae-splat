import numpy as np
import torch
from torch import Tensor
from enum import Enum
from typing import Union, Tuple, List, Optional, Any, Literal

class CameraType(Enum):
    ORTHOGRAPHIC = None
    PERSPECTIVE = "PINHOLE"
    EQUIDISTANT = "FISHEYE"
    EQUISOLID = "EQUISOLID"
    EQUIRECTANGULAR = "EQUIRECTANGULAR"
    CYLINDRICAL = None

class CameraDistortionType(Enum):
    UNDISTORTED = None
    SIMPLE = None  # k1
    RADIAL = None  # k1 k2
    OPENCV = None  # k1 k2 p1 p2
    OPENCV_FISHEYE = None  # k1 k2 k3 k4
    THIN_PRISM_FISHEYE = None  # k1 k2 k3 k4 p1 p2 sx1 sy1
    METASHAPE = None  # k1 k2 k3 k4 p1 p2 b1 b2, with swapped p1/p2 and normalized b1/b2
    FULL_METASHAPE = None  # above + 96 Fourier coefficients
    FULL_OPENCV = None
    RAD_TAN_THIN_PRISM_FISHEYE = None
    SIMPLE_DIVISION = None
    DIVISION = None

_COLMAP_CAMERA_MODEL_TO_TYPE = {
    "SIMPLE_PINHOLE": CameraType.PERSPECTIVE,
    "PINHOLE": CameraType.PERSPECTIVE,
    "SIMPLE_RADIAL": CameraType.PERSPECTIVE,
    "SIMPLE_RADIAL_FISHEYE": CameraType.EQUIDISTANT,
    "RADIAL": CameraType.PERSPECTIVE,
    "RADIAL_FISHEYE": CameraType.EQUIDISTANT,
    "OPENCV": CameraType.PERSPECTIVE,
    "OPENCV_FISHEYE": CameraType.EQUIDISTANT,
    "FULL_OPENCV": None,
    "FOV": None,
    "THIN_PRISM_FISHEYE": CameraType.EQUIDISTANT,
    "RAD_TAN_THIN_PRISM_FISHEYE": None,
    "SIMPLE_DIVISION": None,
    "DIVISION": None,
    "SIMPLE_FISHEYE": CameraType.EQUIDISTANT,
    "FISHEYE": CameraType.EQUIDISTANT,
    # additional ones
    "EQUISOLID": CameraType.EQUISOLID,
    "EQUIRECTANGULAR": CameraType.EQUIRECTANGULAR,
    "CYLINDRICAL": None,
}

def colmap_camera_model_to_type(camera_model: str):
    assert camera_model in _COLMAP_CAMERA_MODEL_TO_TYPE, f"Unknown camera model {camera_model}"
    value = _COLMAP_CAMERA_MODEL_TO_TYPE[camera_model]
    assert value is not None, f"Unsupported camera model {camera_model}"
    return value.value


class Cameras:

    def __init__(
        self,
        intrins: Union[Tensor, Tuple[Tensor, Tensor, Tensor, Tensor], Tuple[float, float, float, float]],
        distortion_params: Union[Tensor, Tuple[Tensor], Tuple[float]],
        height: Union[int, Tensor],
        width: Union[int, Tensor],
        camera_to_worlds: Tensor,
        camera_type: Union[str, List[str]],
        metadata: Optional[dict] = {}
    ):
        # camera to worlds
        if camera_to_worlds.ndim == 2:
            camera_to_worlds = camera_to_worlds[None]
        assert camera_to_worlds.ndim == 3 and \
            camera_to_worlds.shape[1] == 3 and \
            camera_to_worlds.shape[2] == 4, \
            f"camera_to_worlds must be (N, 3, 4), you have {camera_to_worlds.shape}"
        num_cameras = len(camera_to_worlds)

        # intrinsics
        if isinstance(intrins, Tensor):
            if intrins.ndim == 1:
                intrins = intrins[None].repeat(num_cameras, 1)
            assert intrins.ndim == 2 and intrins.shape[-1] == 4, f"Intrinsics must be (N, 4), got {intrins.shape}"
            assert intrins.shape[0] == num_cameras
        elif isinstance(intrins, tuple) or isinstance(intrins, list):
            assert len(intrins) == 4, "Intrinsics must be (N, 4)"
            intrins = [x.float().reshape(num_cameras) if isinstance(x, torch.Tensor)
                       else torch.full((), x).float().repeat(num_cameras)
                       for x in intrins]
            intrins = torch.stack(intrins, dim=-1)
        else:
            raise ValueError("Invalid intrinsics")

        # distortion parameters
        if isinstance(distortion_params, Tensor):
            if distortion_params.ndim == 1:
                distortion_params = distortion_params[None].repeat(num_cameras, 1)
            assert distortion_params.ndim == 2 and distortion_params.shape[0] == num_cameras, distortion_params.shape
        elif isinstance(distortion_params, tuple) or isinstance(distortion_params, list):
            distortion_params = [x.float().reshape(num_cameras) if isinstance(x, torch.Tensor)
                       else torch.full((), x).float().repeat(num_cameras)
                       for x in distortion_params]
            distortion_params = torch.stack(distortion_params, dim=-1)
        else:
            raise ValueError("Invalid distortion parameters")

        # height and width
        if isinstance(height, Tensor):
            assert height.ndim == 1 and len(height) == num_cameras
        else:
            height = torch.full((), height).repeat(num_cameras)
        if isinstance(width, Tensor):
            assert width.ndim == 1 and len(width) == num_cameras
            width = width.int()
        else:
            width = torch.full((), width).repeat(num_cameras)

        # camera type
        if (isinstance(camera_type, list) or isinstance(camera_type, tuple)) and len(camera_type) == 1:
            camera_type = camera_type[0]
        if isinstance(camera_type, str):
            camera_type = [camera_type] * num_cameras
        else:
            assert len(camera_type) == num_cameras

        self.intrins = intrins.float()
        self.distortion_params = distortion_params.float()
        self.width = width.int()  # type: Tensor
        self.height = height.int()
        self.camera_to_worlds = camera_to_worlds.float()
        self.camera_type = camera_type
        self.metadata = metadata

    def rescale(
        self,
        factor: float, 
        rounding_mode: Literal['floor', 'round', 'ceil'] = "floor"
    ):
        self.intrins = self.intrins * factor  # TODO: consider rounding

        assert rounding_mode in ['floor', 'round', 'ceil'], rounding_mode
        if rounding_mode == "floor":
            self.height = (self.height * factor).int()
            self.width = (self.width * factor).int()
        elif rounding_mode == "round":
            self.height = (self.height * factor + 0.5).int()
            self.width = (self.width * factor + 0.5).int()
        elif rounding_mode == "ceil":
            self.height = torch.ceil(self.height * factor).int()
            self.width = torch.ceil(self.width * factor).int()

    def __len__(self):
        return len(self.intrins) if self.intrins.ndim == 2 else 1

    def __getitem__(self, idx: Union[int, Tensor]) -> 'Cameras':
        squeeze = False
        if isinstance(idx, int):
            # idx = torch.full((), idx)
            idx = torch.full((1,), idx)
            squeeze = True
        assert idx.ndim in [0, 1], "key must be scalar or 1-dimensional"
        cam = Cameras(
            intrins=self.intrins[idx],
            distortion_params=self.distortion_params[idx],
            height=self.height[idx],
            width=self.width[idx],
            camera_to_worlds=self.camera_to_worlds[idx],
            camera_type=[self.camera_type[i] for i in idx.flatten().cpu().numpy().tolist()]
        )
        if squeeze:
            cam.intrins = cam.intrins.squeeze(0)
            cam.distortion_params = cam.distortion_params.squeeze(0)
            cam.height = cam.height.squeeze(0)
            cam.width = cam.width.squeeze(0)
            cam.camera_to_worlds = cam.camera_to_worlds.squeeze(0)
            cam.camera_type = cam.camera_type[0]
            idx = idx.item()
        metadata = {}
        for key, value in self.metadata.items():
            if value is None:
                metadata[key] = value
                continue
            if not isinstance(value, torch.Tensor):
                raise NotImplementedError()
            if value.shape[0] != len(self.intrins):
                raise NotImplementedError()
            metadata[key] = value[idx]
        cam.metadata = metadata
        return cam

    def to(self, device):
        self.intrins = self.intrins.to(device)
        self.distortion_params = self.distortion_params.to(device)
        self.height = self.height.to(device)
        self.width = self.width.to(device)
        self.camera_to_worlds = self.camera_to_worlds.to(device)
        for key, value in self.metadata.items():
            if isinstance(value, torch.Tensor):
                self.metadata[key] = value.to(device)
        return self
