"""Python bindings for custom Cuda functions"""

from typing import Optional

import torch
from jaxtyping import Float, Int
from typing import Tuple
from torch import Tensor
from torch.autograd import Function

import spirulae_splat.splat.cuda as _C

from .utils import bin_and_sort_gaussians, compute_cumulative_intersects


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
) -> Tuple[Tensor]:
    """
    TODO
    """
    assert block_width > 1 and block_width <= 16, "block_width must be between 2 and 16"

    if positions.ndimension() != 2 or positions.size(1) != 3:
        raise ValueError("positions must have dimensions (N, 3)")

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
    ) -> Tuple[Tensor]:
        num_points = positions.size(0)
        tile_bounds = (
            (img_width + block_width - 1) // block_width,
            (img_height + block_width - 1) // block_width,
            1,
        )
        block = (block_width, block_width, 1)
        img_size = (img_width, img_height, 1)
        device = positions.device

        num_intersects, cum_tiles_hit = compute_cumulative_intersects(num_tiles_hit)

        if num_intersects < 1:
            out_depth = torch.ones(img_height, img_width, 1, device=device)
            gaussian_ids_sorted = torch.zeros(0, 1, device=device)
            tile_bins = torch.zeros(0, 2, device=device)
            final_idx = torch.zeros(img_height, img_width, device=device)
        else:
            (
                isect_ids_unsorted,
                gaussian_ids_unsorted,
                isect_ids_sorted,
                gaussian_ids_sorted,
                tile_bins,
            ) = bin_and_sort_gaussians(
                num_points,
                num_intersects,
                positions,
                bounds,
                cum_tiles_hit,
                tile_bounds,
                block_width,
            )

            final_idx, out_depth, out_visibility = _C.rasterize_depth_forward(
                tile_bounds, block, img_size,
                *intrins,
                gaussian_ids_sorted, tile_bins,
                positions, axes_u, axes_v,
                opacities, anisotropies,
            )

        ctx.img_width = img_width
        ctx.img_height = img_height
        ctx.num_intersects = num_intersects
        ctx.block_width = block_width
        ctx.save_for_backward(
            torch.tensor(intrins),
            gaussian_ids_sorted, tile_bins,
            positions, axes_u, axes_v,
            opacities, anisotropies,
            final_idx, out_depth, out_visibility,
        )

        if RETURN_IDX:
            return out_depth, out_visibility, final_idx
        return out_depth

    @staticmethod
    def backward(ctx, v_out_depth, v_out_visibility=None, v_idx=None):

        img_height = ctx.img_height
        img_width = ctx.img_width
        num_intersects = ctx.num_intersects

        (
            intrins,
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
            backward_return = _C.rasterize_depth_backward(
                img_height, img_width, ctx.block_width,
                *intrins,
                gaussian_ids_sorted, tile_bins,
                positions, axes_u, axes_v,
                opacities, anisotropies,
                final_idx, out_depth, out_visibility,
                v_out_depth,
            )

            clean = lambda x: torch.nan_to_num(torch.clip(x, -1., 1.))
            (
                v_positions, v_positions_xy_abs,
                v_axes_u, v_axes_v,
                v_opacities, v_anisotropies
            ) = [clean(v) for v in backward_return]
            v_positions_xy_abs = backward_return[1]

        # Abs grad for gaussian splitting criterion. See
        # - "AbsGS: Recovering Fine Details for 3D Gaussian Splatting"
        # - "EfficientGS: Streamlining Gaussian Splatting for Large-Scale High-Resolution Scene Representation"
        if not hasattr(positions, 'absgrad') or positions.absgrad is None:
            positions.absgrad = v_positions_xy_abs
        else:
            positions.absgrad += v_positions_xy_abs

        return (
            v_positions, v_axes_u, v_axes_v,
            v_opacities, v_anisotropies,
            None, None, None, None, None, None,
        )
