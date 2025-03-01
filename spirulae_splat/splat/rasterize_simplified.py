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
timerf = PerfTimer("rasterize_sp_f", ema_tau=100)
timerb = PerfTimer("rasterize_sp_b", ema_tau=100)


RETURN_IDX = False


def rasterize_gaussians_simplified(
    positions: Float[Tensor, "*batch 3"],
    axes_u: Float[Tensor, "*batch 3"],
    axes_v: Float[Tensor, "*batch 3"],
    colors: Float[Tensor, "*batch channels"],
    opacities: Float[Tensor, "*batch 1"],
    anisotropies: Float[Tensor, "*batch 2"],
    bounds: Int[Tensor, "*batch 4"],
    num_tiles_hit: Int[Tensor, "*batch 1"],
    intrins: Tuple[float, float, float, float],
    img_height: int,
    img_width: int,
    block_width: int,
) -> Tuple[Tensor, Tensor, Tensor, Tensor, Tensor]:
    assert block_width > 1 and block_width <= 16, "block_width must be between 2 and 16"

    if colors.dtype == torch.uint8:
        # make sure colors are float [0,1]
        colors = colors.float() / 255

    if not (num_tiles_hit > 0).any():
        shape = (img_height, img_width)
        device = positions.device
        # out_img = background.reshape((1, 1, 3)).repeat((*shape, 1))
        out_img = torch.zeros((*shape, 3)).float().to(device)
        out_alpha = torch.zeros((*shape, 1)).float().to(device)
        depth_im = torch.zeros((*shape, 1)).float().to(device)
        normal_im = torch.zeros((*shape, 3)).float().to(device)
        reg_depth = torch.zeros((*shape, 1)).float().to(device)
        return (out_img, out_alpha, depth_im, normal_im, reg_depth)

    if positions.ndimension() != 2 or positions.size(1) != 3:
        raise ValueError("xys must have dimensions (N, 3)")

    if colors.ndimension() != 2:
        raise ValueError("colors must have dimensions (N, D)")

    return _RasterizeGaussiansSimplified.apply(
        positions.contiguous(),
        axes_u.contiguous(),
        axes_v.contiguous(),
        colors.contiguous(),
        opacities.contiguous(),
        anisotropies.contiguous(),
        bounds.contiguous(),
        num_tiles_hit.contiguous(),
        intrins,
        img_height,
        img_width,
        block_width,
    )


class _RasterizeGaussiansSimplified(Function):
    """Rasterizes 2D gaussians"""

    @staticmethod
    def forward(
        ctx,
        positions: Float[Tensor, "*batch 3"],
        axes_u: Float[Tensor, "*batch 3"],
        axes_v: Float[Tensor, "*batch 3"],
        colors: Float[Tensor, "*batch channels"],
        opacities: Float[Tensor, "*batch 1"],
        anisotropies: Float[Tensor, "*batch 2"],
        bounds: Int[Tensor, "*batch 4"],
        num_tiles_hit: Int[Tensor, "*batch 1"],
        intrins: Tuple[float, float, float, float],
        img_height: int,
        img_width: int,
        block_width: int,
    ) -> Tensor:
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
        timerf.mark("sort")  # us

        if num_intersects < 1:
            shape = (img_height, img_width)
            gaussian_ids_sorted = torch.zeros(0, 1, device=device)
            tile_bins = torch.zeros(0, 2, device=device)
            final_idx = torch.zeros(img_height, img_width, device=device)
            out_img = torch.zeros((*shape, 3)).float().to(device)
            out_alpha = torch.zeros((*shape, 1)).float().to(device)
            out_depth = torch.zeros((*shape, 1)).float().to(device)
            out_normal = torch.zeros((*shape, 3)).float().to(device)
            out_reg_depth = torch.zeros((*shape, 1)).float().to(device)
        else:
            assert colors.shape == torch.Size([num_points, 3])
            assert opacities.shape == torch.Size([num_points, 1])

            (
                final_idx, out_alpha,
                out_img, out_depth, out_normal, out_reg_depth
            ) = _C.rasterize_simplified_forward(
                tile_bounds, block, img_size,
                intrins,
                gaussian_ids_sorted, tile_bins,
                positions, axes_u, axes_v,
                colors, opacities, anisotropies,
            )
        timerf.mark("rasterize")  # us

        ctx.img_width = img_width
        ctx.img_height = img_height
        ctx.num_intersects = num_intersects
        ctx.block_width = block_width
        ctx.intrins = intrins
        ctx.save_for_backward(
            gaussian_ids_sorted, tile_bins,
            positions, axes_u, axes_v,
            colors, opacities, anisotropies,
            final_idx, out_alpha, out_depth,
        )
        timerf.end("save")  # us

        output = (
            out_img, out_alpha,
            out_depth, out_normal, out_reg_depth
        )
        if RETURN_IDX:
            return (*output, final_idx)
        return output

    @staticmethod
    def backward(
        ctx,
        v_out_img,
        v_out_alpha,
        v_out_depth,
        v_out_normal,
        v_out_depth_reg,
        v_idx = None
        ):

        img_height = ctx.img_height
        img_width = ctx.img_width
        num_intersects = ctx.num_intersects
        intrins = ctx.intrins

        (
            gaussian_ids_sorted, tile_bins,
            positions, axes_u, axes_v,
            colors, opacities, anisotropies,
            final_idx, out_alpha, out_depth,
        ) = ctx.saved_tensors

        if num_intersects < 1:
            v_positions = torch.zeros_like(positions)
            v_positions_xy_abs = torch.zeros_like(positions)[..., :2]
            v_axes_u = torch.zeros_like(axes_u)
            v_axes_v = torch.zeros_like(axes_v)
            v_colors = torch.zeros_like(colors)
            v_opacities = torch.zeros_like(opacities)
            v_anisotropies = torch.zeros_like(anisotropies)

        else:
            timerb.start()

            assert colors.shape[-1] == 3
            backward_return = _C.rasterize_simplified_backward(
                img_height, img_width, ctx.block_width,
                intrins,
                gaussian_ids_sorted, tile_bins,
                positions, axes_u, axes_v,
                colors, opacities, anisotropies,
                final_idx, out_alpha, out_depth,
                *[v.contiguous() for v in [
                    v_out_alpha,
                    v_out_img,
                    v_out_depth,
                    v_out_normal,
                    v_out_depth_reg,
                ]]
            )
            timerb.mark("rasterize")  # us

            clean = lambda x: torch.nan_to_num(torch.clip(x, -1., 1.))
            v_positions = clean(backward_return[0])
            v_positions_xy_abs = backward_return[1]
            v_axes_u = clean(backward_return[2])
            v_axes_v = clean(backward_return[3])
            v_colors = clean(backward_return[4])
            v_opacities = clean(backward_return[5])
            v_anisotropies = clean(backward_return[6])
            timerb.mark("clean")  # us

        # Abs grad for gaussian splitting criterion. See
        # - "AbsGS: Recovering Fine Details for 3D Gaussian Splatting"
        # - "EfficientGS: Streamlining Gaussian Splatting for Large-Scale High-Resolution Scene Representation"
        for key, value in [('grad', v_positions), ('absgrad', v_positions_xy_abs)][1:]:
            if not hasattr(positions, key) or getattr(positions, key) is None:
                setattr(positions, key, value)
            else:
                setattr(positions, key, getattr(positions, key)+value)

        if num_intersects >= 1:
            timerb.end("absgrad")  # us

        return (
            v_positions, v_axes_u, v_axes_v,
            v_colors, v_opacities, v_anisotropies,
            None, None, None, None, None, None,
        )
