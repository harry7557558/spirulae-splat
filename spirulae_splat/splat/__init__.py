from typing import Any
import torch
from .project_gaussians import project_gaussians
from .rasterize_simple import rasterize_gaussians_simple
from .rasterize_depth import rasterize_gaussians_depth
from .rasterize import rasterize_gaussians
from .rasterize_simplified import rasterize_gaussians_simplified
from .rasterize_indices import rasterize_gaussians_indices
from .rasterize_simple_sorted import rasterize_gaussians_simple_sorted
from .utils import (
    map_gaussian_to_intersects,
    bin_and_sort_gaussians,
    compute_cumulative_intersects,
    get_tile_bin_edges,
    depth_to_points,
    depth_to_normal,
)
from .sh import spherical_harmonics
import warnings

__version__ = "0.1.0"

__all__ = [
    "__version__",
    "project_gaussians",
    "rasterize_gaussians",
    "spherical_harmonics",
    # utils
    "bin_and_sort_gaussians",
    "compute_cumulative_intersects",
    "get_tile_bin_edges",
    "map_gaussian_to_intersects",
    # Function.apply() will be deprecated
    "ProjectGaussians",
    "RasterizeGaussians",
    "BinAndSortGaussians",
    "ComputeCumulativeIntersects",
    "ComputeCov2dBounds",
    "GetTileBinEdges",
    "MapGaussiansToIntersects",
    "SphericalHarmonics",
    "NDRasterizeGaussians",
]

# Define these for backwards compatibility


class MapGaussiansToIntersects(torch.autograd.Function):
    @staticmethod
    def forward(ctx, *args, **kwargs):
        warnings.warn(
            "MapGaussiansToIntersects is deprecated, use map_gaussian_to_intersects instead",
            DeprecationWarning,
        )
        return map_gaussian_to_intersects(*args, **kwargs)

    @staticmethod
    def backward(ctx: Any, *grad_outputs: Any) -> Any:
        raise NotImplementedError


class ComputeCumulativeIntersects(torch.autograd.Function):
    @staticmethod
    def forward(ctx, *args, **kwargs):
        warnings.warn(
            "ComputeCumulativeIntersects is deprecated, use compute_cumulative_intersects instead",
            DeprecationWarning,
        )
        return compute_cumulative_intersects(*args, **kwargs)

    @staticmethod
    def backward(ctx: Any, *grad_outputs: Any) -> Any:
        raise NotImplementedError


class GetTileBinEdges(torch.autograd.Function):
    @staticmethod
    def forward(ctx, *args, **kwargs):
        warnings.warn(
            "GetTileBinEdges is deprecated, use get_tile_bin_edges instead",
            DeprecationWarning,
        )
        return get_tile_bin_edges(*args, **kwargs)

    @staticmethod
    def backward(ctx: Any, *grad_outputs: Any) -> Any:
        raise NotImplementedError


class BinAndSortGaussians(torch.autograd.Function):
    @staticmethod
    def forward(ctx, *args, **kwargs):
        warnings.warn(
            "BinAndSortGaussians is deprecated, use bin_and_sort_gaussians instead",
            DeprecationWarning,
        )
        return bin_and_sort_gaussians(*args, **kwargs)

    @staticmethod
    def backward(ctx: Any, *grad_outputs: Any) -> Any:
        raise NotImplementedError


class ProjectGaussians(torch.autograd.Function):
    @staticmethod
    def forward(ctx, *args, **kwargs):
        warnings.warn(
            "ProjectGaussians is deprecated, use project_gaussians instead",
            DeprecationWarning,
        )
        return project_gaussians(*args, **kwargs)

    @staticmethod
    def backward(ctx: Any, *grad_outputs: Any) -> Any:
        raise NotImplementedError


class RasterizeGaussians(torch.autograd.Function):
    @staticmethod
    def forward(ctx, *args, **kwargs):
        warnings.warn(
            "RasterizeGaussians is deprecated, use rasterize_gaussians instead",
            DeprecationWarning,
        )
        return rasterize_gaussians(*args, **kwargs)

    @staticmethod
    def backward(ctx: Any, *grad_outputs: Any) -> Any:
        raise NotImplementedError


class NDRasterizeGaussians(torch.autograd.Function):
    @staticmethod
    def forward(ctx, *args, **kwargs):
        warnings.warn(
            "NDRasterizeGaussians is deprecated, use rasterize_gaussians instead",
            DeprecationWarning,
        )
        return rasterize_gaussians(*args, **kwargs)

    @staticmethod
    def backward(ctx: Any, *grad_outputs: Any) -> Any:
        raise NotImplementedError


class SphericalHarmonics(torch.autograd.Function):
    @staticmethod
    def forward(ctx, *args, **kwargs):
        warnings.warn(
            "SphericalHarmonics is deprecated, use spherical_harmonics instead",
            DeprecationWarning,
        )
        return spherical_harmonics(*args, **kwargs)

    @staticmethod
    def backward(ctx: Any, *grad_outputs: Any) -> Any:
        raise NotImplementedError
