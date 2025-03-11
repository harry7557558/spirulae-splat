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

from spirulae_splat.perf_timer import PerfTimer
timerf = PerfTimer("rasterize_depth_f", ema_tau=100)
timerb = PerfTimer("rasterize_depth_b", ema_tau=100)


DEBUG = False


def rasterize_gaussians_depth(
    positions: Float[Tensor, "*batch 3"],
    axes_u: Float[Tensor, "*batch 3"],
    axes_v: Float[Tensor, "*batch 3"],
    opacities: Float[Tensor, "*batch 1"],
    bounds: Int[Tensor, "*batch 4"],
    num_tiles_hit: Int[Tensor, "*batch 1"],
    camera: _Camera,
    depth_mode: str,
) -> Tuple[Tensor]:

    if not (num_tiles_hit > 0).any():
        return torch.zeros((camera.h, camera.w, 1)).float().to(positions)

    if positions.ndimension() != 2 or positions.size(1) != 3:
        raise ValueError("positions must have dimensions (N, 3)")

    if depth_mode not in ["mean", "median"]:
        raise ValueError("Depth mode must be `mean` or `median`.")

    return _RasterizeGaussiansDepth.apply(
        positions.contiguous(),
        axes_u.contiguous(),
        axes_v.contiguous(),
        opacities.contiguous(),
        bounds.contiguous(),
        num_tiles_hit.contiguous(),
        camera,
        depth_mode,
    )


class _RasterizeGaussiansDepth(Function):
    """Rasterizes 2D gaussians"""

    @staticmethod
    def forward(
        ctx,
        positions: Float[Tensor, "*batch 3"],
        axes_u: Float[Tensor, "*batch 3"],
        axes_v: Float[Tensor, "*batch 3"],
        opacities: Float[Tensor, "*batch 1"],
        bounds: Int[Tensor, "*batch 4"],
        num_tiles_hit: Int[Tensor, "*batch 1"],
        camera: _Camera,
        depth_mode: str,
    ) -> Tuple[Tensor]:
        timerf.start()
        device = positions.device

        (
            num_intersects, gaussian_ids_sorted, tile_bins
        ) = rasterize_preprocess(
            positions, bounds, num_tiles_hit,
            camera.h, camera.w,
        )
        timerf.mark("sort")  # 200us-350us

        if num_intersects < 1:
            out_depth = torch.ones(camera.h, camera.w, 1, device=device)
            gaussian_ids_sorted = torch.zeros(0, 1, device=device)
            tile_bins = torch.zeros(0, 2, device=device)
            final_idx = torch.zeros(camera.h, camera.w, device=device)
        else:
            final_idx, out_depth, out_visibility = _C.rasterize_depth_forward(
                depth_mode,
                camera.h, camera.w, camera.BLOCK_WIDTH,
                camera.model, camera.intrins, camera.get_undist_map(),
                gaussian_ids_sorted, tile_bins,
                positions, axes_u, axes_v,
                opacities,
            )
        timerf.mark("rasterize")  # 200us-600us

        ctx.camera = camera
        ctx.num_intersects = num_intersects
        ctx.depth_mode = depth_mode
        ctx.save_for_backward(
            gaussian_ids_sorted, tile_bins,
            positions, axes_u, axes_v,
            opacities,
            final_idx, out_depth, out_visibility,
        )
        timerf.end("save")  # ~10us -> 450us-950us

        if DEBUG:
            return out_depth, out_visibility, final_idx
        return out_depth

    @staticmethod
    def backward(ctx, v_out_depth, v_out_visibility=None, v_idx=None):

        camera = ctx.camera  # type: _Camera
        if camera.is_distorted():
            raise NotImplementedError("Unsupported distorted camera for backward")
        num_intersects = ctx.num_intersects

        (
            gaussian_ids_sorted, tile_bins,
            positions, axes_u, axes_v,
            opacities,
            final_idx, out_depth, out_visibility,
        ) = ctx.saved_tensors

        if num_intersects < 1:
            v_positions = torch.zeros_like(positions)
            v_positions_xy_abs = torch.zeros_like(positions)[..., :2]
            v_axes_u = torch.zeros_like(axes_u)
            v_axes_v = torch.zeros_like(axes_v)
            v_opacities = torch.zeros_like(opacities)

        else:
            timerb.start()

            backward_return = _C.rasterize_depth_backward(
                ctx.depth_mode,
                camera.h, camera.w, camera.BLOCK_WIDTH,
                camera.intrins,
                gaussian_ids_sorted, tile_bins,
                positions, axes_u, axes_v,
                opacities,
                final_idx, out_depth, out_visibility,
                v_out_depth,
            )
            timerb.mark("rasterize")  # 600us-2100us

            if DEBUG:
                clean = lambda x: torch.nan_to_num(x)
            else:
                clean = lambda x: torch.nan_to_num(torch.clip(x, -1., 1.))
            v_positions = clean(backward_return[0])
            v_positions_xy_abs = backward_return[1]
            v_axes_u = clean(backward_return[2])
            v_axes_v = clean(backward_return[3])
            v_opacities = clean(backward_return[4])
            timerb.mark("clean")  # 60us-80us

        # Abs grad for gaussian splitting criterion. See
        # - "AbsGS: Recovering Fine Details for 3D Gaussian Splatting"
        # - "EfficientGS: Streamlining Gaussian Splatting for Large-Scale High-Resolution Scene Representation"
        for key, value in [('grad', v_positions), ('absgrad', v_positions_xy_abs)][1:]:
            if not hasattr(positions, key) or getattr(positions, key) is None:
                setattr(positions, key, value)
            else:
                setattr(positions, key, getattr(positions, key)+value)

        if num_intersects >= 1:
            timerb.end("absgrad")  # ~10us -> 650us-2200us

        return (
            v_positions, v_axes_u, v_axes_v,
            v_opacities,
            None, None, None, None,
        )
