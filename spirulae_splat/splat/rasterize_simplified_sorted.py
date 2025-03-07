"""Python bindings for custom Cuda functions"""

from typing import Optional

import torch
from jaxtyping import Float, Int
from typing import Tuple
from torch import Tensor
from torch.autograd import Function

import spirulae_splat.splat.cuda as _C


def rasterize_gaussians_simplified_sorted(
    positions: Float[Tensor, "*batch 3"],
    axes_u: Float[Tensor, "*batch 3"],
    axes_v: Float[Tensor, "*batch 3"],
    colors: Float[Tensor, "*batch channels"],
    opacities: Float[Tensor, "*batch 1"],
    num_intersects: Int[Tensor, "h w"],
    sorted_indices: Int[Tensor, "h w MAX_SORTED_INDICES"],
    intrins: Tuple[float, float, float, float],
    img_height: int,
    img_width: int,
    block_width: int,
) -> Tuple[Tensor, Tensor, Tensor, Tensor, Tensor]:
    assert block_width > 1 and block_width <= 16, "block_width must be between 2 and 16"

    if colors.dtype == torch.uint8:
        # make sure colors are float [0,1]
        colors = colors.float() / 255

    if positions.ndimension() != 2 or positions.size(1) != 3:
        raise ValueError("xys must have dimensions (N, 3)")

    if colors.ndimension() != 2:
        raise ValueError("colors must have dimensions (N, D)")

    return _RasterizeGaussiansSimplifiedSorted.apply(
        positions.contiguous(),
        axes_u.contiguous(),
        axes_v.contiguous(),
        colors.contiguous(),
        opacities.contiguous(),
        num_intersects.contiguous(),
        sorted_indices.contiguous(),
        intrins,
        img_height, img_width,
        block_width,
    )


class _RasterizeGaussiansSimplifiedSorted(Function):
    """Rasterizes 2D gaussians"""

    @staticmethod
    def forward(
        ctx,
        positions: Float[Tensor, "*batch 3"],
        axes_u: Float[Tensor, "*batch 3"],
        axes_v: Float[Tensor, "*batch 3"],
        colors: Float[Tensor, "*batch channels"],
        opacities: Float[Tensor, "*batch 1"],
        num_intersects: Int[Tensor, "h w"],
        sorted_indices: Int[Tensor, "h w MAX_SORTED_INDICES"],
        intrins: Tuple[float, float, float, float],
        img_height: int,
        img_width: int,
        block_width: int,
    ) -> Tensor:

        (
            out_alpha,
            out_img, out_depth, out_normal, out_reg_depth
        ) = _C.rasterize_simplified_sorted_forward(
            img_height, img_width, block_width,
            intrins,
            sorted_indices,
            positions, axes_u, axes_v,
            colors, opacities,
        )

        ctx.img_width = img_width
        ctx.img_height = img_height
        ctx.num_intersects = num_intersects
        ctx.block_width = block_width
        ctx.intrins = intrins
        ctx.save_for_backward(
            num_intersects, sorted_indices,
            positions, axes_u, axes_v,
            colors, opacities,
            out_alpha, out_depth,
        )

        return (
            out_img, out_alpha,
            out_depth, out_normal, out_reg_depth
        )

    @staticmethod
    def backward(
        ctx,
        v_out_img,
        v_out_alpha,
        v_out_depth,
        v_out_normal,
        v_out_depth_reg
        ):

        img_height = ctx.img_height
        img_width = ctx.img_width
        num_intersects = ctx.num_intersects
        intrins = ctx.intrins

        (
            num_intersects, sorted_indices,
            positions, axes_u, axes_v,
            colors, opacities,
            out_alpha, out_depth,
        ) = ctx.saved_tensors

        assert colors.shape[-1] == 3
        backward_return = _C.rasterize_simplified_sorted_backward(
            img_height, img_width, ctx.block_width,
            intrins,
            num_intersects, sorted_indices,
            positions, axes_u, axes_v,
            colors, opacities,
            out_alpha, out_depth,
            *[v.contiguous() for v in [
                v_out_alpha,
                v_out_img,
                v_out_depth,
                v_out_normal,
                v_out_depth_reg,
            ]]
        )

        clean = lambda x: torch.nan_to_num(torch.clip(x, -1., 1.))
        v_positions = clean(backward_return[0])
        v_positions_xy_abs = backward_return[1]
        v_axes_u = clean(backward_return[2])
        v_axes_v = clean(backward_return[3])
        v_colors = clean(backward_return[4])
        v_opacities = clean(backward_return[5])

        # Abs grad for gaussian splitting criterion. See
        # - "AbsGS: Recovering Fine Details for 3D Gaussian Splatting"
        # - "EfficientGS: Streamlining Gaussian Splatting for Large-Scale High-Resolution Scene Representation"
        for key, value in [('grad', v_positions), ('absgrad', v_positions_xy_abs)][1:]:
            if not hasattr(positions, key) or getattr(positions, key) is None:
                setattr(positions, key, value)
            else:
                setattr(positions, key, getattr(positions, key)+value)

        return (
            v_positions, v_axes_u, v_axes_v,
            v_colors, v_opacities,
            None, None, None, None, None, None,
        )
