import numpy as np
import yaml
import json
from typing import Literal

class Camera:
    def __init__(self, config_path: str, transform_path: str=None):
        with open(config_path, "r") as f:
            data = yaml.safe_load(f)

        self.w = data.get("width", 1280)
        self.h = data.get("height", 720)
        self.fx = data.get("fx", 568.0)
        self.fy = data.get("fy", 568.0)
        self.cx = data.get("cx", self.w / 2)
        self.cy = data.get("cy", self.h / 2)
        self.model = data.get("model", "OPENCV")  # type: Literal["OPENCV", "OPENCV_FISHEYE"]
        self.distortion = tuple(data.get("distortion", [0.0, 0.0, 0.0, 0.0]))

        if transform_path is None:
            self.scale = 1.0
            self.pose = np.eye(4)
        else:
            with open(transform_path, 'r') as f:
                transform = np.array(json.load(f))
            transform /= transform[3, 3]
            transform = np.linalg.inv(transform)
            self.scale = np.linalg.det(transform) ** (1/3)
            self.pose = transform / self.scale

    def scale_resolution(self, sc: float):
        self.w = int(self.w*sc+0.5)
        self.h = int(self.h*sc+0.5)
        self.fx *= sc
        self.fy *= sc
        self.cx *= sc
        self.cy *= sc

    def c2w_to_c2o(self, c2w):
        c2o = self.pose @ c2w
        c2o[:3, 3] *= self.scale
        return c2o

    @property
    def intrins(self):
        return (self.fx, self.fy, self.cx, self.cy)

