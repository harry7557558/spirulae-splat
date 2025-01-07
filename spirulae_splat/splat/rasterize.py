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
timerf = PerfTimer("rasterize_f", ema_tau=100)
timerb = PerfTimer("rasterize_b", ema_tau=100)


RETURN_IDX = False


def rasterize_gaussians(
    positions: Float[Tensor, "*batch 3"],
    axes_u: Float[Tensor, "*batch 3"],
    axes_v: Float[Tensor, "*batch 3"],
    normals: Float[Tensor, "*batch 3"],
    colors: Float[Tensor, "*batch channels"],
    ch_degree_r: int,
    ch_degree_r_to_use: int,
    ch_degree_phi: int,
    ch_degree_phi_to_use: int,
    ch_coeffs: Float[Tensor, "*batch dim_ch channels"],
    opacities: Float[Tensor, "*batch 1"],
    anisotropies: Float[Tensor, "*batch 2"],
    depth_ref: Float[Tensor, "h w 1"],
    # background: Optional[Float[Tensor, "channels"]],
    depth_reg_pairwise_factor: float,
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

    # if background is not None:
    #     assert (
    #         background.shape[0] == colors.shape[-1]
    #     ), f"incorrect shape of background color tensor, expected shape {colors.shape[-1]}"
    # else:
    #     background = torch.zeros(
    #         colors.shape[-1], dtype=torch.float32, device=colors.device
    #     )

    if not (num_tiles_hit > 0).any():
        shape = (img_height, img_width)
        device = positions.device
        # out_img = background.reshape((1, 1, 3)).repeat((*shape, 1))
        out_img = torch.zeros((*shape, 3)).float().to(device)
        out_depth = torch.zeros((*shape, 2)).float().to(device)
        out_normal = torch.zeros((*shape, 3)).float().to(device)
        out_reg_depth = torch.zeros((*shape, 1)).float().to(device)
        out_alpha = torch.zeros((*shape, 1)).float().to(device)
        return (
            out_img, out_alpha,
            out_depth, out_normal, out_reg_depth,
        )

    if positions.ndimension() != 2 or positions.size(1) != 3:
        raise ValueError("xys must have dimensions (N, 3)")

    if colors.ndimension() != 2:
        raise ValueError("colors must have dimensions (N, D)")

    return _RasterizeGaussians.apply(
        positions.contiguous(),
        axes_u.contiguous(),
        axes_v.contiguous(),
        normals.contiguous(),
        colors.contiguous(),
        ch_degree_r, ch_degree_r_to_use,
        ch_degree_phi, ch_degree_phi_to_use,
        ch_coeffs.contiguous(),
        opacities.contiguous(),
        anisotropies.contiguous(),
        depth_ref.contiguous(),
        # background.contiguous(),
        depth_reg_pairwise_factor,
        bounds.contiguous(),
        num_tiles_hit.contiguous(),
        intrins,
        img_height,
        img_width,
        block_width,
    )


class _RasterizeGaussians(Function):
    """Rasterizes 2D gaussians"""

    @staticmethod
    def forward(
        ctx,
        positions: Float[Tensor, "*batch 3"],
        axes_u: Float[Tensor, "*batch 3"],
        axes_v: Float[Tensor, "*batch 3"],
        normals: Float[Tensor, "*batch 3"],
        colors: Float[Tensor, "*batch channels"],
        ch_degree_r: int,
        ch_degree_r_to_use: int,
        ch_degree_phi: int,
        ch_degree_phi_to_use: int,
        ch_coeffs: Float[Tensor, "*batch dim_ch channels"],
        opacities: Float[Tensor, "*batch 1"],
        anisotropies: Float[Tensor, "*batch 2"],
        depth_ref: Float[Tensor, "h w 1"],
        # background: Optional[Float[Tensor, "channels"]],
        depth_reg_pairwise_factor: float,
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

        # TODO: reuse this from previous render
        (
            num_intersects, gaussian_ids_sorted, tile_bins
        ) = rasterize_preprocess(
            positions, bounds, num_tiles_hit,
            tile_bounds, block_width
        )
        timerf.mark("sort")  # 200us-350us

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
            assert depth_ref.shape == torch.Size([img_height, img_width, 1])

            (
                final_idx, out_alpha,
                out_img, out_depth, out_normal,
                out_reg_depth,
            ) = _C.rasterize_forward(
                tile_bounds, block, img_size,
                *intrins,
                depth_reg_pairwise_factor,
                gaussian_ids_sorted, tile_bins,
                positions, axes_u, axes_v, normals,
                colors,
                ch_degree_r, ch_degree_r_to_use,
                ch_degree_phi, ch_degree_phi_to_use,
                ch_coeffs,
                opacities, anisotropies, #background,
                depth_ref,
            )
        timerf.mark("rasterize")  # 250us-600us

        ctx.img_width = img_width
        ctx.img_height = img_height
        ctx.num_intersects = num_intersects
        ctx.block_width = block_width
        ctx.depth_reg_pairwise_factor = depth_reg_pairwise_factor
        ctx.intrins = intrins
        ctx.ch_degrees = (ch_degree_r, ch_degree_r_to_use,
                          ch_degree_phi, ch_degree_phi_to_use)
        ctx.save_for_backward(
            gaussian_ids_sorted, tile_bins,
            positions, axes_u, axes_v, normals,
            colors, ch_coeffs, opacities, anisotropies, #background,
            depth_ref,
            final_idx, out_alpha, out_depth, out_normal, out_reg_depth,
        )
        timerf.end("save")  # ~10us -> 450us-950us

        output = (
            out_img, out_alpha,
            out_depth, out_normal, out_reg_depth,
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
        v_out_reg_depth,
        v_idx = None
        ):

        img_height = ctx.img_height
        img_width = ctx.img_width
        num_intersects = ctx.num_intersects
        intrins = ctx.intrins
        ch_degrees = ctx.ch_degrees

        (
            gaussian_ids_sorted, tile_bins,
            positions, axes_u, axes_v, normals,
            colors, ch_coeffs, opacities, anisotropies, #background,
            depth_ref,
            final_idx, out_alpha, out_depth, out_normal, out_reg_depth,
        ) = ctx.saved_tensors

        if num_intersects < 1:
            v_positions = torch.zeros_like(positions)
            v_positions_xy_abs = torch.zeros_like(positions)[..., :2]
            v_axes_u = torch.zeros_like(axes_u)
            v_axes_v = torch.zeros_like(axes_v)
            v_normals = torch.zeros_like(normals)
            v_colors = torch.zeros_like(colors)
            v_ch_coeffs = torch.zeros_like(ch_coeffs)
            v_opacities = torch.zeros_like(opacities)
            v_anisotropies = torch.zeros_like(anisotropies)
            v_depth_ref = torch.zeros_like(depth_ref)

        else:
            timerb.start()

            assert colors.shape[-1] == 3
            backward_return = _C.rasterize_backward(
                img_height, img_width, ctx.block_width,
                *intrins, *ch_degrees,
                ctx.depth_reg_pairwise_factor,
                gaussian_ids_sorted, tile_bins,
                positions, axes_u, axes_v, normals,
                colors, ch_coeffs, opacities, anisotropies, #background,
                depth_ref,
                final_idx, out_alpha, out_depth,
                v_out_alpha, v_out_img, v_out_depth, v_out_normal,
                v_out_reg_depth.contiguous(),
            )
            timerb.mark("rasterize")  # 600us-2600us

            clean = lambda x: torch.nan_to_num(torch.clip(x, -1., 1.))
            v_positions = clean(backward_return[0])
            v_positions_xy_abs = backward_return[1]
            v_axes_u = clean(backward_return[2])
            v_axes_v = clean(backward_return[3])
            v_normals = clean(backward_return[4])
            v_colors = clean(backward_return[5])
            v_ch_coeffs = clean(backward_return[6])
            # v_ch_coeffs_abs = backward_return[7]
            v_opacities = clean(backward_return[7])
            v_anisotropies = clean(backward_return[8])
            # v_background = backward_return[8]
            v_depth_ref = clean(backward_return[9])
            timerb.mark("clean")  # 150us-200us

        # Abs grad for gaussian splitting criterion. See
        # - "AbsGS: Recovering Fine Details for 3D Gaussian Splatting"
        # - "EfficientGS: Streamlining Gaussian Splatting for Large-Scale High-Resolution Scene Representation"
        for key, value in [('grad', v_positions), ('absgrad', v_positions_xy_abs)][1:]:
            if not hasattr(positions, key) or getattr(positions, key) is None:
                setattr(positions, key, value)
            else:
                setattr(positions, key, getattr(positions, key)+value)
        # ch_coeffs.absgrad = v_ch_coeffs_abs

        if num_intersects >= 1:
            timerb.end("absgrad")  # ~10us -> 900us-2800us

        return (
            v_positions, v_axes_u, v_axes_v, v_normals,
            v_colors, *([None]*4), v_ch_coeffs, v_opacities, v_anisotropies,
            v_depth_ref,
            # v_background,
            None,
            None, None, None, None, None, None,
        )
