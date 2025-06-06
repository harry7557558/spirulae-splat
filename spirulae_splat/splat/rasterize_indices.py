"""Python bindings for custom Cuda functions"""

from typing import Optional

import torch
from jaxtyping import Float, Int
from typing import Tuple
from torch import Tensor
from torch.autograd import Function

import spirulae_splat.splat.cuda as _C
from spirulae_splat.splat._camera import _Camera

from .rasterize_simple import rasterize_preprocess
from .cuda import MAX_SORTED_SPLATS, SORTED_INDEX_INF

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
    camera: _Camera,
) -> Tuple[Tensor, Tensor]:

    timer.start()

    (
        num_intersects, gaussian_ids_sorted, tile_bins
    ) = rasterize_preprocess(
        positions, bounds, num_tiles_hit,
        camera.h, camera.w,
    )
    timer.mark("sort")  # ?us

    if num_intersects == 0:
        num_intersects = torch.zeros(
            (camera.h, camera.w), dtype=torch.int32, device=positions.device)
        indices = torch.empty(
            (camera.h, camera.w, MAX_SORTED_SPLATS),
            dtype=torch.int32, device=positions.device).fill_(SORTED_INDEX_INF)
        return num_intersects, indices

    num_intersects, indices, depths = _C.rasterize_indices(
        camera.h, camera.w,
        camera.model, camera.intrins, camera.get_undist_map(),
        gaussian_ids_sorted, tile_bins,
        positions, axes_u, axes_v,
        opacities,
    )
    timer.mark("rasterize")  # ?us

    num_intersects, indices, depths = [x.contiguous() for x in (num_intersects, indices, depths)]

    # _test_sort(camera.h, camera.w, num_intersects, indices, depths)

    if False:
        sorted_depths, argsort = torch.sort(depths, dim=-1)
        timer.mark("sort")
        indices = torch.gather(indices, dim=-1, index=argsort)
        timer.end("gather")
    else:
        # benchmark in one dataset (compiled without vs with using shared memory):
        # - insertion: ~7ms | ~2.7ms
        # - quick: ~6ms | stack overflow
        # - heap: ~14ms | ~3.3ms
        # - random_quick: ~6ms | stack overflow
        # - network: ? | slowest
        # recursive ones (quick, random_quick) may cause stack overflow depending on device
        # i.e. "CUDA error: an illegal memory access was encountered"
        _C.sort_per_pixel(
            "insertion",
            camera.h, camera.w,
            num_intersects, indices, depths
        )
        timer.end("sort")

    return num_intersects, indices


def _test_sort(h, w, num_intersects, indices, depths):
    """Making this work requires changing `#if 0` at the end of `rasterize_indices_kernel()` to `1`
        Neglectable performance difference in practice"""

    # reference
    sorted_depths_ref, argsort = torch.sort(depths, dim=-1)
    sorted_indices_ref = torch.gather(indices, dim=-1, index=argsort)

    # custom sorting algorithm
    sorted_indices = indices.clone()
    sorted_depths = depths.clone()
    _C.sort_per_pixel(
        "network",
        h, w, num_intersects, sorted_indices, sorted_depths
    )

    torch.testing.assert_close(sorted_depths, sorted_depths_ref)
    # torch.testing.assert_close(sorted_indices, sorted_indices_ref)
