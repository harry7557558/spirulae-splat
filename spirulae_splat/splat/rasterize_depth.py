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
timerf = PerfTimer("rasterize_depth_f", ema_tau=100)
timerb = PerfTimer("rasterize_depth_b", ema_tau=100)


RETURN_IDX = False


def rasterize_gaussians_depth(
    positions: Float[Tensor, "*batch 3"],
    axes_u: Float[Tensor, "*batch 3"],
    axes_v: Float[Tensor, "*batch 3"],
    opacities: Float[Tensor, "*batch 1"],
    anisotropies: Float[Tensor, "*batch 2"],
    bounds: Int[Tensor, "*batch 4"],
    num_tiles_hit: Int[Tensor, "*batch 1"],
    intrins: Tuple[float, float, float, float],
    img_height: int,
    img_width: int,
    block_width: int,
    depth_mode: str,
) -> Tuple[Tensor]:
    assert block_width > 1 and block_width <= 16, "block_width must be between 2 and 16"

    if not (num_tiles_hit > 0).any():
        return torch.zeros((img_height, img_width, 1)).float().to(positions)

    if positions.ndimension() != 2 or positions.size(1) != 3:
        raise ValueError("positions must have dimensions (N, 3)")

    if isinstance(depth_mode, str):
        if depth_mode.startswith("mean"):
            depth_mode = 0
        elif depth_mode.startswith("median"):
            depth_mode = 1
    if depth_mode not in [0, 1]:
        raise ValueError("Depth mode must start with `mean` or `median`.")

    return _RasterizeGaussiansDepth.apply(
        positions.contiguous(),
        axes_u.contiguous(),
        axes_v.contiguous(),
        opacities.contiguous(),
        anisotropies.contiguous(),
        bounds.contiguous(),
        num_tiles_hit.contiguous(),
        intrins,
        img_height, img_width,
        block_width,
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
        anisotropies: Float[Tensor, "*batch 2"],
        bounds: Int[Tensor, "*batch 4"],
        num_tiles_hit: Int[Tensor, "*batch 1"],
        intrins: Tuple[float, float, float, float],
        img_height: int,
        img_width: int,
        block_width: int,
        depth_mode: int,
    ) -> Tuple[Tensor]:
        timerf.start()
        num_points = positions.size(0)
        tile_bounds = (
            (img_width + block_width - 1) // block_width,
            (img_height + block_width - 1) // block_width,
            1,
        )
        block = (block_width, block_width, 1)
        img_size = (img_width, img_height, 1)
        device = positions.device

        (
            num_intersects, gaussian_ids_sorted, tile_bins
        ) = rasterize_preprocess(
            positions, bounds, num_tiles_hit,
            tile_bounds, block_width
        )
        timerf.mark("sort")  # 200ns-350ns

        if num_intersects < 1:
            out_depth = torch.ones(img_height, img_width, 1, device=device)
            gaussian_ids_sorted = torch.zeros(0, 1, device=device)
            tile_bins = torch.zeros(0, 2, device=device)
            final_idx = torch.zeros(img_height, img_width, device=device)
        else:
            final_idx, out_depth, out_visibility = _C.rasterize_depth_forward(
                depth_mode,
                tile_bounds, block, img_size,
                *intrins,
                gaussian_ids_sorted, tile_bins,
                positions, axes_u, axes_v,
                opacities, anisotropies,
            )
        timerf.mark("rasterize")  # 200ns-600ns

        ctx.img_width = img_width
        ctx.img_height = img_height
        ctx.num_intersects = num_intersects
        ctx.block_width = block_width
        ctx.depth_mode = depth_mode
        ctx.intrins = intrins
        ctx.save_for_backward(
            gaussian_ids_sorted, tile_bins,
            positions, axes_u, axes_v,
            opacities, anisotropies,
            final_idx, out_depth, out_visibility,
        )
        timerf.end("save")  # ~10ns -> 450ns-950ns

        if RETURN_IDX:
            return out_depth, out_visibility, final_idx
        return out_depth

    @staticmethod
    def backward(ctx, v_out_depth, v_out_visibility=None, v_idx=None):

        img_height = ctx.img_height
        img_width = ctx.img_width
        num_intersects = ctx.num_intersects
        intrins = ctx.intrins

        (
            gaussian_ids_sorted, tile_bins,
            positions, axes_u, axes_v,
            opacities, anisotropies,
            final_idx, out_depth, out_visibility,
        ) = ctx.saved_tensors

        if num_intersects < 1:
            v_positions = torch.zeros_like(positions)
            v_positions_xy_abs = torch.zeros_like(positions)[..., :2]
            v_axes_u = torch.zeros_like(axes_u)
            v_axes_v = torch.zeros_like(axes_v)
            v_opacities = torch.zeros_like(opacities)
            v_anisotropies = torch.zeros_like(anisotropies)

        else:
            timerb.start()

            backward_return = _C.rasterize_depth_backward(
                ctx.depth_mode,
                img_height, img_width, ctx.block_width,
                *intrins,
                gaussian_ids_sorted, tile_bins,
                positions, axes_u, axes_v,
                opacities, anisotropies,
                final_idx, out_depth, out_visibility,
                v_out_depth,
            )
            timerb.mark("rasterize")  # 600ns-2100ns

            clean = lambda x: torch.nan_to_num(torch.clip(x, -1., 1.))
            v_positions = clean(backward_return[0])
            v_positions_xy_abs = backward_return[1]
            v_axes_u = clean(backward_return[2])
            v_axes_v = clean(backward_return[3])
            v_opacities = clean(backward_return[4])
            v_anisotropies = clean(backward_return[5])
            timerb.mark("clean")  # 60ns-80ns

        # Abs grad for gaussian splitting criterion. See
        # - "AbsGS: Recovering Fine Details for 3D Gaussian Splatting"
        # - "EfficientGS: Streamlining Gaussian Splatting for Large-Scale High-Resolution Scene Representation"
        for key, value in [('grad', v_positions), ('absgrad', v_positions_xy_abs)][1:]:
            if not hasattr(positions, key) or getattr(positions, key) is None:
                setattr(positions, key, value)
            else:
                setattr(positions, key, getattr(positions, key)+value)

        if num_intersects >= 1:
            timerb.end("absgrad")  # ~10ns -> 650ns-2200ns

        return (
            v_positions, v_axes_u, v_axes_v,
            v_opacities, v_anisotropies,
            None, None, None, None, None, None, None
        )
