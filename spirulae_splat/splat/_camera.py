"""Lightweight camera that doesn't break with nerfstudio update"""

import torch
from typing import Tuple, Literal, Dict, List, Optional
from jaxtyping import Float

import spirulae_splat.splat.cuda as _C


class _Camera:
    _undist_maps: Dict[Tuple, Float[torch.Tensor, "h w 2"]] = {}
    _undist_maps_list: List[Tuple[Tuple, int]] = []
    _undist_maps_cache_size: List[int] = [0]
    _undist_maps_max_cache_size: int = 1024*1024*2*128  # number of floats, 1GB

    def __init__(self,
                 height: int,
                 width: int,
                 model: Literal["", "OPENCV", "OPENCV_FISHEYE"],
                 intrins: Tuple[float, float, float, float],
                 dist_coeffs: Tuple[float, float, float, float]=(0.,0.,0.,0.,),
                 device="cuda"
            ):
        self.h = height
        self.w = width
        self.model = "" if model == "OPENCV" else model
        self.intrins = intrins
        self.dist_coeffs = dist_coeffs
        if all([x == 0 for x in dist_coeffs]):
            self.model = ""
        if len(dist_coeffs) > 4:
            dist_coeffs = [*dist_coeffs]
            if any([x != 0 for x in dist_coeffs[4:]]) and False:
                raise ValueError("Only support at most 4 distortion coefficients")
            self.dist_coeffs = tuple(dist_coeffs[:4])
        self.device = device

    def validate_model(self):
        assert self.model in ["", "OPENCV", "OPENCV_FISHEYE"], \
            f"Invalid camera model: {self.model}"

    def is_distorted(self):
        return self.model != ""

    def get_undist_map(self, always=False) -> Optional[torch.Tensor]:
        tup = (self.h, self.w, self.model, self.intrins, self.dist_coeffs, str(self.device))
        if tup in self._undist_maps:
            return self._undist_maps[tup]

        def update_cache(undist_map):
            nonlocal tup
            self._undist_maps[tup] = undist_map

            # limit cache size
            map_size = torch.numel(undist_map)
            self._undist_maps_list.append((tup, map_size))
            self._undist_maps_cache_size[0] += map_size
            while self._undist_maps_cache_size[0] > self._undist_maps_max_cache_size:
                if len(self._undist_maps_list) == 0:
                    break
                tup, size = self._undist_maps_list[0]
                del self._undist_maps_list[0]
                del self._undist_maps[tup]
                self._undist_maps_cache_size[0] -= size

            return undist_map

        # not distorted
        if self.model == "":

            if always:
                fx, fy, cx, cy = self.intrins
                x, y = torch.meshgrid(
                    torch.arange(self.w).to(self.device) + 0.5,
                    torch.arange(self.h).to(self.device) + 0.5,
                    indexing="xy",
                )  # [H, W]
                undist_map = torch.stack([(x - cx) / fx, (y - cy) / fy], dim=-1)
                return update_cache(undist_map)

            return None

        # fisheye
        with torch.cuda.device(self.device):
            undist_map = _C.render_undistortion_map(
                self.w, self.h, self.model,
                self.intrins, self.dist_coeffs
            )
        return update_cache(undist_map)
