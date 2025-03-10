"""Lightweight camera that doesn't break with nerfstudio update"""

import torch
from typing import Tuple, Literal, Dict, List, Optional
from jaxtyping import Float

import spirulae_splat.splat.cuda as _C


class _Camera:
    _undist_maps: Dict[Tuple, Float[torch.Tensor, "h w 2"]] = {}
    _undist_maps_list: List[Tuple[Tuple, int]] = []
    _undist_maps_cache_size: int = 0
    _undist_maps_max_cache_size: int = 1024*1024*2*128  # number of floats, 1GB

    BLOCK_WIDTH = 16

    def __init__(self,
                 height: int,
                 width: int,
                 model: Literal["", "OPENCV", "OPENCV_FISHEYE"],
                 intrins: Tuple[float, float, float, float],
                 dist_coeffs: Tuple[float, float, float, float]=(0.,0.,0.,0.,)
            ):
        self.h = height
        self.w = width
        self.model = model
        self.intrins = intrins
        self.dist_coeffs = dist_coeffs
        if all([x == 0 for x in dist_coeffs]):
            self.model = ""

    def validate_model(self):
        assert self.model in ["", "OPENCV", "OPENCV_FISHEYE"], \
            f"Invalid camera model: {self.model}"
        assert 1 < self.BLOCK_WIDTH <= 16, "block_width must be between 2 and 16"

    def is_distorted(self):
        return self.model != ""

    def get_undist_map(self) -> Optional[torch.Tensor]:
        tup = (self.h, self.w, self.model, self.intrins, self.dist_coeffs)
        if tup in self._undist_maps:
            return self._undist_maps[tup]

        if self.model == "":
            return None

        undist_map = _C.render_undistortion_map(
            self.w, self.h, self.model,
            self.intrins, self.dist_coeffs
        )
        self._undist_maps[tup] = undist_map
        return undist_map

        # TODO: limit cache size
