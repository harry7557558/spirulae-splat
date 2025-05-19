"""Python bindings for custom Cuda functions"""

from typing import Optional

import torch
from jaxtyping import Float, Int
from typing import Tuple
from torch import Tensor
from torch.autograd import Function

import spirulae_splat.splat.cuda as _C
from spirulae_splat.splat._camera import _Camera



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
    camera: _Camera,
) -> Tuple[Tensor, Tensor, Tensor, Tensor, Tensor]:

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
        camera,
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
        camera: _Camera,
    ) -> Tensor:

        (
            out_alpha,
            out_img, out_depth, out_normal,
            out_reg_depth,
        ) = _C.rasterize_sorted_forward(
            camera.h, camera.w, camera.intrins,
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

        ctx.camera = camera
        ctx.num_intersects = num_intersects
        ctx.depth_reg_pairwise_factor = depth_reg_pairwise_factor
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

        camera = ctx.camera  # type: _Camera
        num_intersects = ctx.num_intersects
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
            camera.h, camera.w, camera.intrins,
            *ch_degrees,
            ctx.depth_reg_pairwise_factor,
            num_intersects, sorted_indices,
            positions, axes_u, axes_v,
            colors, ch_coeffs, opacities, #background,
            depth_ref,
            out_alpha, out_depth,
            v_out_alpha, v_out_img, v_out_depth, v_out_normal,
            v_out_reg_depth.contiguous(),
        )

        def clean(x, h=1.0):
            return torch.nan_to_num(torch.clip(x, -h, h))
        v_positions = clean(backward_return[0])
        v_positions_xy_abs = clean(backward_return[1], 10.0)
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
            None, None, None, None,
        )
