from typing import Callable

from ._backend import _C


def _get_constant(name: str) -> int:
    return getattr(_C, name)


def _make_lazy_cuda_func(name: str) -> Callable:
    def call_cuda(*args, **kwargs):
        return getattr(_C, name)(*args, **kwargs)

    return call_cuda


# constants
BLOCK_WIDTH = _get_constant("BLOCK_WIDTH")
MAX_SORTED_SPLATS = _get_constant("MAX_SORTED_SPLATS")
SORTED_INDEX_INF = _get_constant("SORTED_INDEX_INF")

# splat rasterization
rasterize_simple_forward = _make_lazy_cuda_func("rasterize_simple_forward")
rasterize_simple_backward = _make_lazy_cuda_func("rasterize_simple_backward")
rasterize_depth_forward = _make_lazy_cuda_func("rasterize_depth_forward")
rasterize_depth_backward = _make_lazy_cuda_func("rasterize_depth_backward")
rasterize_forward = _make_lazy_cuda_func("rasterize_forward")
rasterize_backward = _make_lazy_cuda_func("rasterize_backward")
rasterize_simplified_forward = _make_lazy_cuda_func("rasterize_simplified_forward")
rasterize_simplified_backward = _make_lazy_cuda_func("rasterize_simplified_backward")

# sorted splat rasterization
rasterize_indices = _make_lazy_cuda_func("rasterize_indices")
sort_per_pixel = _make_lazy_cuda_func("sort_per_pixel")
rasterize_simple_sorted_forward = _make_lazy_cuda_func("rasterize_simple_sorted_forward")
rasterize_simple_sorted_backward = _make_lazy_cuda_func("rasterize_simple_sorted_backward")
rasterize_depth_sorted_forward = _make_lazy_cuda_func("rasterize_depth_sorted_forward")
rasterize_depth_sorted_backward = _make_lazy_cuda_func("rasterize_depth_sorted_backward")
rasterize_sorted_forward = _make_lazy_cuda_func("rasterize_sorted_forward")
rasterize_sorted_backward = _make_lazy_cuda_func("rasterize_sorted_backward")
rasterize_simplified_sorted_forward = _make_lazy_cuda_func("rasterize_simplified_sorted_forward")
rasterize_simplified_sorted_backward = _make_lazy_cuda_func("rasterize_simplified_sorted_backward")

# background
render_background_sh_forward = _make_lazy_cuda_func("render_background_sh_forward")
render_background_sh_backward = _make_lazy_cuda_func("render_background_sh_backward")
render_undistortion_map = _make_lazy_cuda_func("render_undistortion_map")

# projection
map_gaussian_to_intersects = _make_lazy_cuda_func("map_gaussian_to_intersects")
get_tile_bin_edges = _make_lazy_cuda_func("get_tile_bin_edges")
project_gaussians_forward = _make_lazy_cuda_func("project_gaussians_forward")
project_gaussians_backward = _make_lazy_cuda_func("project_gaussians_backward")
compute_sh_forward = _make_lazy_cuda_func("compute_sh_forward")
compute_sh_backward = _make_lazy_cuda_func("compute_sh_backward")

# misc
compute_relocation = _make_lazy_cuda_func("compute_relocation")
compute_relocation_split = _make_lazy_cuda_func("compute_relocation_split")
