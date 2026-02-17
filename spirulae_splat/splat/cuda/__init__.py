from typing import Callable

try:
    from ._backend import _C
except:
    # possibly from setup.py
    _C = None

from ._wrapper import *

def _get_constant(name: str) -> int:
    return getattr(_C, name, -1)


def _make_lazy_cuda_func(name: str) -> Callable:

    def call_cuda(*args, **kwargs):
        global _C
        if _C == None:
            # should throw an error
            from ._backend import _C

        return getattr(_C, name)(*args, **kwargs)

    return call_cuda


# constants
BLOCK_WIDTH = _get_constant("BLOCK_WIDTH")
MAX_SORTED_SPLATS = _get_constant("MAX_SORTED_SPLATS")
SORTED_INDEX_INF = _get_constant("SORTED_INDEX_INF")

# background
render_background_sh_forward = _make_lazy_cuda_func("render_background_sh_forward")
render_background_sh_backward = _make_lazy_cuda_func("render_background_sh_backward")
render_undistortion_map = _make_lazy_cuda_func("render_undistortion_map")

# projection
map_gaussian_to_intersects = _make_lazy_cuda_func("map_gaussian_to_intersects")
get_tile_bin_edges = _make_lazy_cuda_func("get_tile_bin_edges")
compute_sh_forward = _make_lazy_cuda_func("compute_sh_forward")
compute_sh_backward = _make_lazy_cuda_func("compute_sh_backward")

# misc
blend_background_forward = _make_lazy_cuda_func("blend_background_forward")
blend_background_backward = _make_lazy_cuda_func("blend_background_backward")

# sparse voxel hash grid
svhash_create_initial_volume = _make_lazy_cuda_func("svhash_create_initial_volume")
svhash_get_voxels = _make_lazy_cuda_func("svhash_get_voxels")
svhash_split_voxels = _make_lazy_cuda_func("svhash_split_voxels")
