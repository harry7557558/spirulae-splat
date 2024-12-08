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


def rasterize_gaussians(
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
    anisotropies: Float[Tensor, "*batch 2"],
    depth_grads: Float[Tensor, "*batch 2"],
    depth_normal_ref: Float[Tensor, "*batch 2"],
    bounds: Int[Tensor, "*batch 4"],
    num_tiles_hit: Int[Tensor, "*batch 1"],
    intrins: Tuple[float, float, float, float],
    img_height: int,
    img_width: int,
    block_width: int,
    background: Optional[Float[Tensor, "channels"]] = None
) -> Tuple[Tensor, Tensor, Tensor, Tensor, Tensor]:
    """
    TODO
    """
    assert block_width > 1 and block_width <= 16, "block_width must be between 2 and 16"

    if colors.dtype == torch.uint8:
        # make sure colors are float [0,1]
        colors = colors.float() / 255

    if background is not None:
        assert (
            background.shape[0] == colors.shape[-1]
        ), f"incorrect shape of background color tensor, expected shape {colors.shape[-1]}"
    else:
        background = torch.zeros(
            colors.shape[-1], dtype=torch.float32, device=colors.device
        )

    if not (num_tiles_hit > 0).any():
        shape = (img_height, img_width)
        device = positions.device
        out_img = background.reshape((1, 1, 3)).repeat((*shape, 1))
        out_depth_grad = torch.zeros((*shape, 4)).float().to(device)
        out_reg_depth = torch.zeros(shape).float().to(device)
        out_reg_normal = torch.zeros(shape).float().to(device)
        out_alpha = torch.zeros(shape).float().to(device)
        return (
            out_img, out_depth_grad,
            out_reg_depth, out_reg_normal,
            out_alpha
        )

    if positions.ndimension() != 2 or positions.size(1) != 3:
        raise ValueError("xys must have dimensions (N, 3)")

    if colors.ndimension() != 2:
        raise ValueError("colors must have dimensions (N, D)")

    return _RasterizeGaussians.apply(
        positions.contiguous(),
        axes_u.contiguous(),
        axes_v.contiguous(),
        colors.contiguous(),
        ch_degree_r, ch_degree_r_to_use,
        ch_degree_phi, ch_degree_phi_to_use,
        ch_coeffs.contiguous(),
        opacities.contiguous(),
        anisotropies.contiguous(),
        depth_grads.contiguous(),
        depth_normal_ref.contiguous(),
        bounds.contiguous(),
        num_tiles_hit.contiguous(),
        intrins,
        img_height,
        img_width,
        block_width,
        background.contiguous()
    )


class _RasterizeGaussians(Function):
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
        anisotropies: Float[Tensor, "*batch 2"],
        depth_grads: Float[Tensor, "*batch 2"],
        depth_normal_ref: Float[Tensor, "*batch 2"],
        bounds: Int[Tensor, "*batch 4"],
        num_tiles_hit: Int[Tensor, "*batch 1"],
        intrins: Tuple[float, float, float, float],
        img_height: int,
        img_width: int,
        block_width: int,
        background: Optional[Float[Tensor, "channels"]] = None
    ) -> Tensor:
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
            out_img = (
                torch.ones(img_height, img_width, colors.shape[-1], device=device)
                * background
            )
            gaussian_ids_sorted = torch.zeros(0, 1, device=device)
            tile_bins = torch.zeros(0, 2, device=device)
            final_idx = torch.zeros(img_height, img_width, device=device)
            out_alpha = torch.zeros(img_height, img_width, device=device)
            out_reg_depth = torch.zeros(img_height, img_width, device=device)
            out_reg_normal = torch.zeros(img_height, img_width, device=device)
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
            assert colors.shape == torch.Size([num_points, 3])
            assert opacities.shape == torch.Size([num_points, 1])
            assert depth_grads.shape == torch.Size([num_points, 2])
            assert depth_normal_ref.shape == torch.Size([img_height, img_width, 3])

            (
                final_idx, out_alpha,
                out_img, out_depth_grad,
                out_reg_depth, out_reg_normal,
            ) = _C.rasterize_forward(
                tile_bounds, block, img_size,
                *intrins,
                gaussian_ids_sorted, tile_bins,
                positions, axes_u, axes_v,
                colors,
                ch_degree_r, ch_degree_r_to_use,
                ch_degree_phi, ch_degree_phi_to_use,
                ch_coeffs,
                opacities, anisotropies, background,
                depth_grads, depth_normal_ref,
            )

        ctx.img_width = img_width
        ctx.img_height = img_height
        ctx.num_intersects = num_intersects
        ctx.block_width = block_width
        ctx.save_for_backward(
            torch.tensor(intrins),
            torch.tensor([ch_degree_r, ch_degree_r_to_use,
                          ch_degree_phi, ch_degree_phi_to_use]),
            gaussian_ids_sorted, tile_bins,
            positions, axes_u, axes_v,
            colors, ch_coeffs, opacities, anisotropies, background,
            depth_grads, depth_normal_ref,
            final_idx, out_alpha, out_depth_grad,
        )

        output = (
            out_img, out_depth_grad,
            out_reg_depth, out_reg_normal,
            out_alpha
        )
        if RETURN_IDX:
            return (*output, final_idx)
        return output

    @staticmethod
    def backward(
        ctx,
        v_out_img,
        v_out_depth_grad,
        v_out_reg_depth,
        v_out_reg_normal,
        v_out_alpha,
        v_idx = None
        ):

        img_height = ctx.img_height
        img_width = ctx.img_width
        num_intersects = ctx.num_intersects

        (
            intrins, ch_degrees,
            gaussian_ids_sorted, tile_bins,
            positions, axes_u, axes_v,
            colors, ch_coeffs, opacities, anisotropies, background,
            depth_grads, depth_normal_ref,
            final_idx, out_alpha, out_depth_grad,
        ) = ctx.saved_tensors

        if num_intersects < 1:
            v_positions = torch.zeros_like(positions)
            v_positions_xy_abs = torch.zeros_like(positions)[..., :2]
            v_axes_u = torch.zeros_like(axes_u)
            v_axes_v = torch.zeros_like(axes_v)
            v_colors = torch.zeros_like(colors)
            v_ch_coeffs = torch.zeros_like(ch_coeffs)
            v_opacities = torch.zeros_like(opacities)
            v_anisotropies = torch.zeros_like(anisotropies)
            v_depth_grads = torch.zeros_like(depth_grads)
            v_depth_normal_ref = torch.zeros_like(depth_normal_ref)

        else:
            assert colors.shape[-1] == 3
            backward_return = _C.rasterize_backward(
                img_height, img_width, ctx.block_width,
                *intrins, *ch_degrees,
                gaussian_ids_sorted, tile_bins,
                positions, axes_u, axes_v,
                colors, ch_coeffs, opacities, anisotropies, background,
                depth_grads, depth_normal_ref,
                final_idx, out_alpha, out_depth_grad,
                v_out_alpha, v_out_img, v_out_depth_grad,
                v_out_reg_depth.contiguous(),
                v_out_reg_normal.contiguous(),
            )

            clean = lambda x: torch.nan_to_num(torch.clip(x, -1., 1.))
            (
                v_positions, v_positions_xy_abs,
                v_axes_u, v_axes_v,
                v_colors, v_ch_coeffs, v_ch_coeffs_abs,
                v_opacities, v_anisotropies,
                v_depth_grads, v_depth_normal_ref
            ) = [clean(v) for v in backward_return]
            v_positions_xy_abs = backward_return[1]
            v_ch_coeffs_abs = backward_return[6]

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
        ch_coeffs.absgrad = v_ch_coeffs_abs

        return (
            v_positions, v_axes_u, v_axes_v,
            v_colors, *([None]*4), v_ch_coeffs, v_opacities, v_anisotropies,
            v_depth_grads, v_depth_normal_ref,
            None, None, None, None, None, None,
            v_background,
        )
