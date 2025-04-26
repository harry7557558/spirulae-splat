"""Python bindings for custom Cuda functions"""

from typing import Optional

import torch
from jaxtyping import Float, Int
from typing import Tuple
from torch import Tensor
from torch.autograd import Function

import spirulae_splat.splat.cuda as _C
from spirulae_splat.splat._camera import _Camera


DEBUG = False


def rasterize_gaussians_depth_sorted(
    positions: Float[Tensor, "*batch 3"],
    axes_u: Float[Tensor, "*batch 3"],
    axes_v: Float[Tensor, "*batch 3"],
    opacities: Float[Tensor, "*batch 1"],
    sorted_indices: Int[Tensor, "h w MAX_SORTED_INDICES"],
    camera: _Camera,
    depth_mode: str,
) -> Tensor:
    if positions.ndimension() != 2 or positions.size(1) != 3:
        raise ValueError("positions must have dimensions (N, 3)")

    if depth_mode not in ["mean", "median"]:
        raise ValueError("Depth mode be `mean` or `median`.")

    return _RasterizeGaussiansDepthSorted.apply(
        positions.contiguous(),
        axes_u.contiguous(),
        axes_v.contiguous(),
        opacities.contiguous(),
        sorted_indices.contiguous(),
        camera,
        depth_mode,
    )



class _RasterizeGaussiansDepthSorted(Function):
    """Rasterizes 2D gaussians"""

    @staticmethod
    def forward(
        ctx,
        positions: Float[Tensor, "*batch 3"],
        axes_u: Float[Tensor, "*batch 3"],
        axes_v: Float[Tensor, "*batch 3"],
        opacities: Float[Tensor, "*batch 1"],
        sorted_indices: Int[Tensor, "h w MAX_SORTED_INDICES"],
        camera: _Camera,
        depth_mode: str,
    ) -> Tuple[Tensor, Tensor]:

        final_idx, out_depth, out_visibility = _C.rasterize_depth_sorted_forward(
            depth_mode,
            camera.h, camera.w,
            camera.model, camera.intrins, camera.get_undist_map(),
            sorted_indices,
            positions, axes_u, axes_v,
            opacities,
        )

        ctx.camera = camera
        ctx.depth_mode = depth_mode
        ctx.save_for_backward(
            final_idx, sorted_indices,
            positions, axes_u, axes_v, opacities,
            out_depth, out_visibility,
        )

        if DEBUG:
            return out_depth, out_visibility
        return out_depth

    @staticmethod
    def backward(ctx, v_out_depth, v_out_visibility=None):

        camera = ctx.camera  # type: _Camera
        if camera.is_distorted():
            raise NotImplementedError("Unsupported distorted camera for backward")

        (
            final_idx, sorted_indices,
            positions, axes_u, axes_v, opacities,
            out_depth, out_visibility,
        ) = ctx.saved_tensors

        backward_return = _C.rasterize_depth_sorted_backward(
            ctx.depth_mode,
            camera.h, camera.w,
            camera.intrins,
            final_idx, sorted_indices,
            positions, axes_u, axes_v, opacities,
            out_depth, out_visibility,
            v_out_depth,
        )

        if DEBUG:
            clean = lambda x: torch.nan_to_num(x)
        else:
            clean = lambda x: torch.nan_to_num(torch.clip(x, -1., 1.))
        v_positions = clean(backward_return[0])
        v_positions_xy_abs = backward_return[1]
        v_axes_u = clean(backward_return[2])
        v_axes_v = clean(backward_return[3])
        v_opacities = clean(backward_return[4])

        # Abs grad for gaussian splitting criterion. See
        # - "AbsGS: Recovering Fine Details for 3D Gaussian Splatting"
        # - "EfficientGS: Streamlining Gaussian Splatting for Large-Scale High-Resolution Scene Representation"
        for key, value in [('grad', v_positions), ('absgrad', v_positions_xy_abs)][1:]:
            if not hasattr(positions, key) or getattr(positions, key) is None:
                setattr(positions, key, value)
            else:
                setattr(positions, key, getattr(positions, key)+value)

        return (
            v_positions, v_axes_u, v_axes_v, v_opacities,
            None, None, None,
        )
