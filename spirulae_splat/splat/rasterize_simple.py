"""Python bindings for custom Cuda functions"""

from typing import Optional

import torch
from jaxtyping import Float, Int
from typing import Tuple
from torch import Tensor
from torch.autograd import Function

import spirulae_splat.splat.cuda as _C
from spirulae_splat.splat._camera import _Camera

from .utils import bin_and_sort_gaussians, compute_cumulative_intersects


RETURN_IDX = False


def rasterize_gaussians_simple(
    positions: Float[Tensor, "*batch 3"],
    axes_u: Float[Tensor, "*batch 3"],
    axes_v: Float[Tensor, "*batch 3"],
    colors: Float[Tensor, "*batch channels"],
    opacities: Float[Tensor, "*batch 1"],
    bounds: Int[Tensor, "*batch 4"],
    num_tiles_hit: Int[Tensor, "*batch 1"],
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
    else:
        background = torch.zeros(
            colors.shape[-1], dtype=torch.float32, device='cpu'
        )

    if not (num_tiles_hit > 0).any():
        shape = (camera.h, camera.w)
        device = positions.device
        return (
            background.reshape((1, 1, 3)).repeat((*shape, 1)),
            torch.zeros(shape).float().to(device)
        )

    if positions.ndimension() != 2 or positions.size(1) != 3:
        raise ValueError("positions must have dimensions (N, 3)")

    if colors.ndimension() != 2:
        raise ValueError("colors must have dimensions (N, D)")

    return _RasterizeGaussiansSimple.apply(
        positions.contiguous(),
        axes_u.contiguous(),
        axes_v.contiguous(),
        colors.contiguous(),
        opacities.contiguous(),
        bounds.contiguous(),
        num_tiles_hit.contiguous(),
        camera,
        background.contiguous(),
    )


def rasterize_preprocess(
        positions: Float[Tensor, "*batch 3"],
        bounds: Int[Tensor, "*batch 4"],
        num_tiles_hit: Int[Tensor, "*batch 1"],
        img_height: int,
        img_width: int,
    ):

    num_points = positions.size(0)

    num_intersects, cum_tiles_hit = compute_cumulative_intersects(num_tiles_hit)

    if num_intersects < 1:
        return num_intersects, None, None

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
        img_height, img_width,
        _Camera.BLOCK_WIDTH,
    )

    return num_intersects, gaussian_ids_sorted, tile_bins


class _RasterizeGaussiansSimple(Function):
    """Rasterizes 2D gaussians"""

    @staticmethod
    def forward(
        ctx,
        positions: Float[Tensor, "*batch 3"],
        axes_u: Float[Tensor, "*batch 3"],
        axes_v: Float[Tensor, "*batch 3"],
        colors: Float[Tensor, "*batch channels"],
        opacities: Float[Tensor, "*batch 1"],
        bounds: Int[Tensor, "*batch 4"],
        num_tiles_hit: Int[Tensor, "*batch 1"],
        camera: _Camera,
        background: Float[Tensor, "channels"],
    ) -> Tuple[Tensor, Tensor]:
        device = positions.device

        (
            num_intersects, gaussian_ids_sorted, tile_bins
        ) = rasterize_preprocess(
            positions, bounds, num_tiles_hit,
            camera.h, camera.w,
        )

        if num_intersects < 1:
            out_img = (
                torch.ones(camera.h, camera.w, colors.shape[-1], device=device)
                * background.to(device)
            )
            gaussian_ids_sorted = torch.zeros(0, 1, device=device)
            tile_bins = torch.zeros(0, 2, device=device)
            final_idx = torch.zeros(camera.h, camera.w, device=device)
            out_alpha = torch.ones(camera.h, camera.w, device=device)
        else:
            final_idx, out_img, out_alpha = _C.rasterize_simple_forward(
                camera.h, camera.w, camera.BLOCK_WIDTH,
                camera.model, camera.intrins, camera.get_undist_map(),
                gaussian_ids_sorted, tile_bins,
                positions, axes_u, axes_v,
                colors, opacities,
                background.cpu(),
            )

        ctx.num_intersects = num_intersects
        ctx.camera = camera
        ctx.save_for_backward(
            gaussian_ids_sorted, tile_bins,
            positions, axes_u, axes_v,
            colors, opacities, background,
            final_idx, out_alpha,
        )

        if RETURN_IDX:
            return out_img, out_alpha, final_idx
        return out_img, out_alpha

    @staticmethod
    def backward(ctx, v_out_img, v_out_alpha, v_idx=None):

        num_intersects = ctx.num_intersects
        camera = ctx.camera  # type: _Camera
        if camera.is_distorted():
            raise NotImplementedError("Unsupported distorted camera for backward")

        (
            gaussian_ids_sorted, tile_bins,
            positions, axes_u, axes_v,
            colors, opacities, background,
            final_idx, out_alpha,
        ) = ctx.saved_tensors

        if num_intersects < 1:
            v_positions = torch.zeros_like(positions)
            v_positions_xy_abs = torch.zeros_like(positions)[..., :2]
            v_axes_u = torch.zeros_like(axes_u)
            v_axes_v = torch.zeros_like(axes_v)
            v_colors = torch.zeros_like(colors)
            v_opacities = torch.zeros_like(opacities)

        else:
            backward_return = _C.rasterize_simple_backward(
                camera.h, camera.w, camera.BLOCK_WIDTH,
                camera.intrins,
                gaussian_ids_sorted, tile_bins,
                positions, axes_u, axes_v,
                colors, opacities, background.cpu(),
                final_idx, out_alpha,
                v_out_img, v_out_alpha,
            )

            clean = lambda x: torch.nan_to_num(torch.clip(x, -1., 1.))
            (
                v_positions, v_positions_xy_abs,
                v_axes_u, v_axes_v,
                v_colors, v_opacities,
            ) = [clean(v) for v in backward_return]
            v_positions_xy_abs = backward_return[1]

        v_background = None
        if background.requires_grad:
            v_background = torch.matmul(
                v_out_img.float().view(-1, 3).t(),
                (1.0-out_alpha).float().view(-1, 1)
            ).squeeze().to(background.device)

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
            v_background,
        )
