"""Python bindings for custom Cuda functions"""

from typing import Optional

import torch
from jaxtyping import Float, Int
from typing import Tuple
from torch import Tensor
from torch.autograd import Function

import spirulae_splat.splat.cuda as _C

from .rasterize_simple import rasterize_preprocess

from spirulae_splat.perf_timer import PerfTimer
timer = PerfTimer("per_pixel_sort", ema_tau=5)


RETURN_IDX = False


def rasterize_gaussians_indices(
    positions: Float[Tensor, "*batch 3"],
    axes_u: Float[Tensor, "*batch 3"],
    axes_v: Float[Tensor, "*batch 3"],
    opacities: Float[Tensor, "*batch 1"],
    bounds: Int[Tensor, "*batch 4"],
    num_tiles_hit: Int[Tensor, "*batch 1"],
    intrins: Tuple[float, float, float, float],
    img_height: int,
    img_width: int,
    block_width: int,
) -> Tuple[Tensor, Tensor]:
    assert block_width > 1 and block_width <= 16, "block_width must be between 2 and 16"
    tile_bounds = (
        (img_width + block_width - 1) // block_width,
        (img_height + block_width - 1) // block_width,
        1,
    )
    block = (block_width, block_width, 1)
    img_size = (img_width, img_height, 1)

    timer.start()

    (
        num_intersects, gaussian_ids_sorted, tile_bins
    ) = rasterize_preprocess(
        positions, bounds, num_tiles_hit,
        tile_bounds, block_width
    )
    timer.mark("sort")  # ?us

    size_params = (tile_bounds, block, img_size)
    num_intersects, indices, depths = _C.rasterize_indices(
        *size_params, intrins,
        gaussian_ids_sorted, tile_bins,
        positions, axes_u, axes_v,
        opacities,
    )
    timer.mark("rasterize")  # ?us

    # _test_sort(size_params, num_intersects, indices, depths)

    if False:
        sorted_depths, argsort = torch.sort(depths, dim=-1)
        timer.mark("sort")
        indices = torch.gather(indices, dim=-1, index=argsort)
        timer.end("gather")
    else:
        # benchmark in one dataset:
        # - insertion: ~7ms
        # - quick: ~6ms
        # - heap: ~14ms
        # - random_quick: ~6ms
        # recursive ones (quick, random_quick) may cause stack overflow on some devices
        # i.e. "CUDA error: an illegal memory access was encountered"
        _C.sort_per_pixel(
            "quick",
            *size_params,
            num_intersects, indices, depths
        )
        timer.end("sort")

    return num_intersects, indices


def _test_sort(size_params, num_intersects, indices, depths):
    """Making this work requires changing `#if 0` at the end of `rasterize_indices_kernel()` to `1`
        Neglectable performance difference in practice"""

    # reference
    sorted_depths_ref, argsort = torch.sort(depths, dim=-1)
    sorted_indices_ref = torch.gather(indices, dim=-1, index=argsort)

    # custom sorting algorithm
    sorted_indices = indices.clone()
    sorted_depths = depths.clone()
    _C.sort_per_pixel(
        "random_quick",
        *size_params,
        num_intersects, sorted_indices, sorted_depths
    )

    torch.testing.assert_close(sorted_depths, sorted_depths_ref)
    # torch.testing.assert_close(sorted_indices, sorted_indices_ref)
