"""Python bindings for 3D gaussian projection"""

from typing import Optional, Tuple

import torch
from jaxtyping import Float
from torch import Tensor
from torch.autograd import Function

import spirulae_splat.splat.cuda as _C
from spirulae_splat.splat._camera import _Camera


def project_gaussians(
    means3d: Float[Tensor, "*batch 3"],
    scales: Float[Tensor, "*batch 2"],
    quats: Float[Tensor, "*batch 4"],
    viewmat: Float[Tensor, "4 4"],
    camera: _Camera,
    clip_thresh: float = 0.01,
) -> Tuple[Tensor, Tensor, Tensor, Tensor, Tensor, Tensor]:
    assert (quats.norm(dim=-1) - 1 < 1e-4).all(), "quats must be normalized"
    camera.validate_model()

    return _ProjectGaussians.apply(
        means3d.contiguous(),
        scales.contiguous(),
        quats.contiguous(),
        viewmat.contiguous(),
        camera,
        clip_thresh,
    )


class _ProjectGaussians(Function):
    """Project 3D gaussians to 2D."""

    @staticmethod
    def forward(
        ctx,
        means3d: Float[Tensor, "*batch 3"],
        scales: Float[Tensor, "*batch 2"],
        quats: Float[Tensor, "*batch 4"],
        viewmat: Float[Tensor, "4 4"],
        camera: _Camera,
        clip_thresh: float = 0.01,
    ):
        num_points = means3d.shape[-2]
        if num_points < 1 or means3d.shape[-1] != 3:
            raise ValueError(f"Invalid shape for means3d: {means3d.shape}")

        (
            bounds, num_tiles_hit,
            positions, axes_u, axes_v,
        ) = _C.project_gaussians_forward(
            num_points,
            means3d,
            scales,
            quats,
            viewmat,
            camera.model, camera.intrins, camera.dist_coeffs,
            camera.h, camera.w,
            clip_thresh,
        )

        # Save non-tensors.
        ctx.img_height = camera.h
        ctx.img_width = camera.w
        ctx.num_points = num_points
        ctx.camera = camera

        # Save tensors.
        ctx.save_for_backward(
            means3d, scales, quats,
            viewmat,
            bounds, num_tiles_hit,
            positions, axes_u, axes_v,
        )

        return (positions, axes_u, axes_v, bounds, num_tiles_hit)

    @staticmethod
    def backward(
        ctx,
        v_positions, v_axes_u, v_axes_v,
        v_bounds, v_num_tiles_hit,
    ):
        camera = ctx.camera  # type: _Camera

        (
            means3d, scales, quats,
            viewmat,
            bounds, num_tiles_hit,
            positions, axes_u, axes_v,
        ) = ctx.saved_tensors

        backward_return = _C.project_gaussians_backward(
            ctx.num_points,
            means3d, scales, quats,
            viewmat,
            camera.intrins,
            num_tiles_hit,
            v_positions, v_axes_u, v_axes_v
        )

        clean = lambda x, h: torch.nan_to_num(torch.clip(x, -h, h))
        (v_means3d, v_scales, v_quats) = [clean(v, 10.) for v in backward_return[:3]]
        v_viewmat = clean(backward_return[3], 20.)  # 4x4
        v_viewmat = v_viewmat[:viewmat.shape[0], :viewmat.shape[1]]  # 4x4 or 3x4

        # Return a gradient for each input.
        return (
            # means3d: Float[Tensor, "*batch 3"],
            v_means3d,
            # scales: Float[Tensor, "*batch 2"],
            v_scales,
            # quats: Float[Tensor, "*batch 4"],
            v_quats,
            # viewmat: Float[Tensor, "4 4"],
            v_viewmat,
            # camera: _Camera
            None,
            # clip_thresh,
            None,
        )
