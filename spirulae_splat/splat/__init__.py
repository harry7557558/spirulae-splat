from .project_gaussians import project_gaussians
from .rasterize_simple import rasterize_gaussians_simple
from .rasterize_depth import rasterize_gaussians_depth
from .rasterize import rasterize_gaussians
from .rasterize_simplified import rasterize_gaussians_simplified
from .rasterize_indices import rasterize_gaussians_indices
from .rasterize_simple_sorted import rasterize_gaussians_simple_sorted
from .rasterize_depth_sorted import rasterize_gaussians_depth_sorted
from .rasterize_sorted import rasterize_gaussians_sorted
from .rasterize_simplified_sorted import rasterize_gaussians_simplified_sorted
from .utils import (
    bin_and_sort_gaussians,
    compute_cumulative_intersects,
)
from .sh import spherical_harmonics
from .cuda import (
    BLOCK_WIDTH,
)

from torch import Tensor

try:
    import ssplat_camera_utils

    from spirulae_splat.splat._camera import _Camera
    from typing import Optional

    def depth_to_normal(
        depths: Tensor, camera: _Camera, c2w: Optional[Tensor]=None, z_depth: bool = True, alpha: Optional[Tensor] = None,
        return_points = False
    ):
        assert c2w is None
        normals = ssplat_camera_utils.depth_to_normal(depths, camera.get_undist_map(always=True), z_depth)
        if alpha is not None:
            return normals * (alpha>0).float().reshape((*normals.shape[:2], 1)), alpha
        if return_points:
            points = ssplat_camera_utils.depth_to_points(depths, camera.get_undist_map(always=True), z_depth)
            return normals, points
        return normals

except ImportError:
    from .utils import depth_to_normal
    print("ssplat_camera_utils not found, fall back to PyTorch depth_to_normal")


__version__ = "0.1.0"

__all__ = [
    "__version__",
    "project_gaussians",
    "rasterize_gaussians_simple",
    "rasterize_gaussians_depth",
    "rasterize_gaussians",
    "rasterize_gaussians_simplified",
    "rasterize_gaussians_indices",
    "rasterize_gaussians_simple_sorted",
    "rasterize_gaussians_depth_sorted",
    "rasterize_gaussians_sorted",
    "rasterize_gaussians_simplified_sorted",
    "spherical_harmonics",
    # utils
    "bin_and_sort_gaussians",
    "compute_cumulative_intersects",
    "map_gaussian_to_intersects",
    "depth_to_normal",
    "BLOCK_WIDTH",
]
