import numpy as np
import yaml
import json
from typing import Literal, Union

from spirulae_splat.splat._camera import _Camera

class Camera:
    def __init__(self, config_path: str, transform_path: str=None):
        if isinstance(config_path, str) and config_path.split('.')[-1] in ['yaml', 'yml']:
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

        elif isinstance(config_path, dict) or config_path.endswith("transforms.json"):
            config = config_path
            if isinstance(config, str):
                with open(config, 'r') as fp:
                    config = json.load(fp)
            if 'w' not in config:
                config = config['frames'][0]
            self.w = config['w']
            self.h = config['h']
            self.fx = config['fl_x']
            self.fy = config['fl_y']
            self.cx = config['cx']
            self.cy = config['cy']
            self.model = config['camera_model']
            if self.model == "OPENCV":
                self.distortion = tuple(config[k] for k in ['k1', 'k2', 'p1', 'p2'])
            elif self.model == "OPENCV_FISHEYE":
                self.distortion = tuple(config[k] for k in ['k1', 'k2', 'k3', 'k4'])

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

    def resize_rel(self, sc: float):
        self.w = int(self.w*sc+0.5)
        self.h = int(self.h*sc+0.5)
        self.fx *= sc
        self.fy *= sc
        self.cx *= sc
        self.cy *= sc

    def resize(self, w: int, h: int):
        self.fx *= w / self.w
        self.fy *= h / self.h
        self.cx *= w / self.w
        self.cy *= h / self.h
        self.w = int(w)
        self.h = int(h)

    def c2w_to_c2o(self, c2w):
        c2o = self.pose @ c2w
        c2o[:3, 3] *= self.scale
        return c2o

    @property
    def intrins(self):
        return (self.fx, self.fy, self.cx, self.cy)

    def _to_ssplat_camera(self, device="cuda"):
        return _Camera(self.h, self.w, self.model, self.intrins, self.distortion, device)
