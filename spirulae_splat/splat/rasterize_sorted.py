"""Python bindings for custom Cuda functions"""

from typing import Optional

import torch
from jaxtyping import Float, Int
from typing import Tuple
from torch import Tensor
from torch.autograd import Function

import spirulae_splat.splat.cuda as _C



def rasterize_gaussians_sorted(
    positions: Float[Tensor, "*batch 3"],
    axes_u: Float[Tensor, "*batch 3"],
    axes_v: Float[Tensor, "*batch 3"],
    colors: Float[Tensor, "*batch channels"],
    ch_degree_r: int,
    ch_degree_r_to_use: int,
    ch_degree_phi: int,
    ch_degree_phi_to_use: int,
    ch_coeffs: Float[Tensor, "*batch dim_ch channels"],
    opacities: Float[Tensor, "*batch 1"],
    depth_ref: Float[Tensor, "h w 1"],
    # background: Optional[Float[Tensor, "channels"]],
    depth_reg_pairwise_factor: float,
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

    return _RasterizeGaussiansSorted.apply(
        positions.contiguous(),
        axes_u.contiguous(),
        axes_v.contiguous(),
        colors.contiguous(),
        ch_degree_r, ch_degree_r_to_use,
        ch_degree_phi, ch_degree_phi_to_use,
        ch_coeffs.contiguous(),
        opacities.contiguous(),
        depth_ref.contiguous(),
        # background.contiguous(),
        depth_reg_pairwise_factor,
        num_intersects.contiguous(),
        sorted_indices.contiguous(),
        intrins,
        img_height,
        img_width,
        block_width,
    )


class _RasterizeGaussiansSorted(Function):
    """Rasterizes 2D gaussians"""

    @staticmethod
    def forward(
        ctx,
        positions: Float[Tensor, "*batch 3"],
        axes_u: Float[Tensor, "*batch 3"],
        axes_v: Float[Tensor, "*batch 3"],
        colors: Float[Tensor, "*batch channels"],
        ch_degree_r: int,
        ch_degree_r_to_use: int,
        ch_degree_phi: int,
        ch_degree_phi_to_use: int,
        ch_coeffs: Float[Tensor, "*batch dim_ch channels"],
        opacities: Float[Tensor, "*batch 1"],
        depth_ref: Float[Tensor, "h w 1"],
        # background: Optional[Float[Tensor, "channels"]],
        depth_reg_pairwise_factor: float,
        num_intersects: Int[Tensor, "h w"],
        sorted_indices: Int[Tensor, "h w MAX_SORTED_INDICES"],
        intrins: Tuple[float, float, float, float],
        img_height: int,
        img_width: int,
        block_width: int,
    ) -> Tensor:

        (
            out_alpha,
            out_img, out_depth, out_normal,
            out_reg_depth,
        ) = _C.rasterize_sorted_forward(
            img_height, img_width, block_width,
            intrins,
            depth_reg_pairwise_factor,
            sorted_indices,
            positions, axes_u, axes_v,
            colors,
            ch_degree_r, ch_degree_r_to_use,
            ch_degree_phi, ch_degree_phi_to_use,
            ch_coeffs,
            opacities, #background,
            depth_ref,
        )

        ctx.img_width = img_width
        ctx.img_height = img_height
        ctx.num_intersects = num_intersects
        ctx.block_width = block_width
        ctx.depth_reg_pairwise_factor = depth_reg_pairwise_factor
        ctx.intrins = intrins
        ctx.ch_degrees = (ch_degree_r, ch_degree_r_to_use,
                          ch_degree_phi, ch_degree_phi_to_use)
        ctx.save_for_backward(
            num_intersects, sorted_indices,
            positions, axes_u, axes_v,
            colors, ch_coeffs, opacities, #background,
            depth_ref,
            out_alpha, out_depth, out_normal, out_reg_depth,
        )

        return (
            out_img, out_alpha,
            out_depth, out_normal, out_reg_depth,
        )

    @staticmethod
    def backward(
        ctx,
        v_out_img,
        v_out_alpha,
        v_out_depth,
        v_out_normal,
        v_out_reg_depth,
        ):

        img_height = ctx.img_height
        img_width = ctx.img_width
        num_intersects = ctx.num_intersects
        intrins = ctx.intrins
        ch_degrees = ctx.ch_degrees

        (
            num_intersects, sorted_indices,
            positions, axes_u, axes_v,
            colors, ch_coeffs, opacities, #background,
            depth_ref,
            out_alpha, out_depth, out_normal, out_reg_depth,
        ) = ctx.saved_tensors

        assert colors.shape[-1] == 3
        backward_return = _C.rasterize_sorted_backward(
            img_height, img_width, ctx.block_width,
            intrins, *ch_degrees,
            ctx.depth_reg_pairwise_factor,
            num_intersects, sorted_indices,
            positions, axes_u, axes_v,
            colors, ch_coeffs, opacities, #background,
            depth_ref,
            out_alpha, out_depth,
            v_out_alpha, v_out_img, v_out_depth, v_out_normal,
            v_out_reg_depth.contiguous(),
        )

        clean = lambda x: torch.nan_to_num(torch.clip(x, -1., 1.))
        v_positions = clean(backward_return[0])
        v_positions_xy_abs = backward_return[1]
        v_axes_u = clean(backward_return[2])
        v_axes_v = clean(backward_return[3])
        v_colors = clean(backward_return[4])
        v_ch_coeffs = clean(backward_return[5])
        # v_ch_coeffs_abs = backward_return[6]
        v_opacities = clean(backward_return[6])
        # v_background = backward_return[8]
        v_depth_ref = clean(backward_return[7])

        # Abs grad for gaussian splitting criterion. See
        # - "AbsGS: Recovering Fine Details for 3D Gaussian Splatting"
        # - "EfficientGS: Streamlining Gaussian Splatting for Large-Scale High-Resolution Scene Representation"
        for key, value in [('grad', v_positions), ('absgrad', v_positions_xy_abs)][1:]:
            if not hasattr(positions, key) or getattr(positions, key) is None:
                setattr(positions, key, value)
            else:
                setattr(positions, key, getattr(positions, key)+value)
        # ch_coeffs.absgrad = v_ch_coeffs_abs

        return (
            v_positions, v_axes_u, v_axes_v,
            v_colors, *([None]*4), v_ch_coeffs, v_opacities,
            v_depth_ref,
            # v_background,
            None,
            None, None, None, None, None, None,
        )
