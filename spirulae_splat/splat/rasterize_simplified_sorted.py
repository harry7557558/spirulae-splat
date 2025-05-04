"""Python bindings for custom Cuda functions"""

from typing import Optional

import torch
from jaxtyping import Float, Int
from typing import Tuple
from torch import Tensor
from torch.autograd import Function

import spirulae_splat.splat.cuda as _C
from spirulae_splat.splat._camera import _Camera


def rasterize_gaussians_simplified_sorted(
    positions: Float[Tensor, "*batch 3"],
    axes_u: Float[Tensor, "*batch 3"],
    axes_v: Float[Tensor, "*batch 3"],
    colors: Float[Tensor, "*batch channels"],
    opacities: Float[Tensor, "*batch 1"],
    num_intersects: Int[Tensor, "h w"],
    sorted_indices: Int[Tensor, "h w MAX_SORTED_INDICES"],
    camera: _Camera,
    intersect_count_reg_start: int
) -> Tuple[Tensor, Tensor, Tensor, Tensor, Tensor, Tensor]:

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
        camera,
        intersect_count_reg_start
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
        camera: _Camera,
        intersect_count_reg_start: int
    ) -> Tensor:

        (
            out_alpha,
            out_img, out_depth, out_normal,
            out_reg_depth, out_reg_intersect_count
        ) = _C.rasterize_simplified_sorted_forward(
            camera.h, camera.w,
            camera.model, camera.intrins, camera.get_undist_map(),
            sorted_indices,
            positions, axes_u, axes_v,
            colors, opacities,
            intersect_count_reg_start
        )

        ctx.camera = camera
        ctx.intersect_count_reg_start = intersect_count_reg_start
        ctx.num_intersects = num_intersects
        ctx.save_for_backward(
            num_intersects, sorted_indices,
            positions, axes_u, axes_v,
            colors, opacities,
            out_alpha, out_depth,
        )

        return (
            out_img, out_alpha,
            out_depth, out_normal,
            out_reg_depth, out_reg_intersect_count
        )

    @staticmethod
    def backward(
        ctx,
        v_out_img,
        v_out_alpha,
        v_out_depth,
        v_out_normal,
        v_out_depth_reg,
        v_out_intersect_count_reg,
        ):

        camera = ctx.camera  # type: _Camera
        intersect_count_reg_start = ctx.intersect_count_reg_start
        num_intersects = ctx.num_intersects

        (
            num_intersects, sorted_indices,
            positions, axes_u, axes_v,
            colors, opacities,
            out_alpha, out_depth,
        ) = ctx.saved_tensors

        assert colors.shape[-1] == 3
        backward_return = _C.rasterize_simplified_sorted_backward(
            camera.h, camera.w,
            camera.model, camera.intrins, camera.get_undist_map(),
            num_intersects, sorted_indices,
            positions, axes_u, axes_v,
            colors, opacities,
            intersect_count_reg_start,
            out_alpha, out_depth,
            *[v.contiguous() for v in [
                v_out_alpha,
                v_out_img,
                v_out_depth,
                v_out_normal,
                v_out_depth_reg,
                v_out_intersect_count_reg,
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
            None, None, None, None
        )
