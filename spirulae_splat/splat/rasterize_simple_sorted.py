"""Python bindings for custom Cuda functions"""

from typing import Optional

import torch
from jaxtyping import Float, Int
from typing import Tuple
from torch import Tensor
from torch.autograd import Function

import spirulae_splat.splat.cuda as _C
from spirulae_splat.splat._camera import _Camera


def rasterize_gaussians_simple_sorted(
    positions: Float[Tensor, "*batch 3"],
    axes_u: Float[Tensor, "*batch 3"],
    axes_v: Float[Tensor, "*batch 3"],
    colors: Float[Tensor, "*batch channels"],
    opacities: Float[Tensor, "*batch 1"],
    num_intersects: Int[Tensor, "h w"],
    sorted_indices: Int[Tensor, "h w MAX_SORTED_INDICES"],
    camera: _Camera,
    background: Optional[Float[Tensor, "channels"]] = None
) -> Tuple[Tensor, Tensor]:
    if colors.dtype == torch.uint8:
        # make sure colors are float [0,1]
        colors = colors.float() / 255

    if background is not None:
        assert (
            background.shape[0] == colors.shape[-1]
        ), f"incorrect shape of background color tensor, expected shape {colors.shape[-1]}"
        background = background.cpu()
    else:
        background = torch.zeros(
            colors.shape[-1], dtype=torch.float32, device="cpu"
        )

    if positions.ndimension() != 2 or positions.size(1) != 3:
        raise ValueError("positions must have dimensions (N, 3)")

    if colors.ndimension() != 2:
        raise ValueError("colors must have dimensions (N, D)")

    return _RasterizeGaussiansSimpleSorted.apply(
        positions.contiguous(),
        axes_u.contiguous(),
        axes_v.contiguous(),
        colors.contiguous(),
        opacities.contiguous(),
        num_intersects.contiguous(),
        sorted_indices.contiguous(),
        camera,
        background.contiguous(),
    )



class _RasterizeGaussiansSimpleSorted(Function):
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
        background: Float[Tensor, "channels"],
    ) -> Tuple[Tensor, Tensor]:

        out_img, out_alpha = _C.rasterize_simple_sorted_forward(
            camera.h, camera.w, camera.model,
            camera.intrins, camera.get_undist_map(),
            sorted_indices,
            positions, axes_u, axes_v,
            colors, opacities,
            background,
        )

        ctx.camera = camera
        ctx.intrins = camera.intrins
        ctx.save_for_backward(
            num_intersects, sorted_indices,
            positions, axes_u, axes_v,
            colors, opacities, background,
            out_alpha,
        )

        return out_img, out_alpha

    @staticmethod
    def backward(ctx, v_out_img, v_out_alpha):

        camera = ctx.camera  # type: _Camera
        if camera.is_distorted():
            raise NotImplementedError("Unsupported distorted camera for backward")
        intrins = ctx.intrins

        (
            num_intersects, sorted_indices,
            positions, axes_u, axes_v,
            colors, opacities, background,
            out_alpha,
        ) = ctx.saved_tensors

        backward_return = _C.rasterize_simple_sorted_backward(
            camera.h, camera.w,
            intrins,
            num_intersects, sorted_indices,
            positions, axes_u, axes_v,
            colors, opacities, background,
            out_alpha,
            v_out_img, v_out_alpha,
        )

        def clean(x, h=1.0):
            return torch.nan_to_num(torch.clip(x, -h, h))
        (
            v_positions, v_positions_xy_abs,
            v_axes_u, v_axes_v,
            v_colors, v_opacities,
        ) = [clean(v) for v in backward_return]
        v_positions_xy_abs = clean(backward_return[1], 10.0)

        v_background = None
        if background.requires_grad:
            v_background = torch.matmul(
                v_out_img.float().view(-1, 3).t(),
                (1.0-out_alpha).float().view(-1, 1)
            ).squeeze()

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
            None, None, None,
            v_background.cpu()
        )
