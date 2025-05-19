"""Python bindings for custom Cuda functions"""

from typing import Optional

import torch
from jaxtyping import Float, Int
from typing import Tuple
from torch import Tensor
from torch.autograd import Function

import spirulae_splat.splat.cuda as _C
from spirulae_splat.splat._camera import _Camera

from .rasterize_simple import rasterize_preprocess

from spirulae_splat.perf_timer import PerfTimer
timerf = PerfTimer("rasterize_f", ema_tau=100)
timerb = PerfTimer("rasterize_b", ema_tau=100)


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
    depth_ref: Float[Tensor, "h w 1"],
    # background: Optional[Float[Tensor, "channels"]],
    depth_reg_pairwise_factor: float,
    bounds: Int[Tensor, "*batch 4"],
    num_tiles_hit: Int[Tensor, "*batch 1"],
    camera: _Camera
) -> Tuple[Tensor, Tensor, Tensor, Tensor, Tensor]:

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
        shape = (camera.h, camera.w)
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
        colors.contiguous(),
        ch_degree_r, ch_degree_r_to_use,
        ch_degree_phi, ch_degree_phi_to_use,
        ch_coeffs.contiguous(),
        opacities.contiguous(),
        depth_ref.contiguous(),
        # background.contiguous(),
        depth_reg_pairwise_factor,
        bounds.contiguous(),
        num_tiles_hit.contiguous(),
        camera
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
        depth_ref: Float[Tensor, "h w 1"],
        # background: Optional[Float[Tensor, "channels"]],
        depth_reg_pairwise_factor: float,
        bounds: Int[Tensor, "*batch 4"],
        num_tiles_hit: Int[Tensor, "*batch 1"],
        camera: _Camera
    ) -> Tensor:
        timerf.start()
        num_points = positions.size(0)
        device = positions.device

        # TODO: reuse this from previous render
        (
            num_intersects, gaussian_ids_sorted, tile_bins
        ) = rasterize_preprocess(
            positions, bounds, num_tiles_hit,
            camera.h, camera.w,
        )
        timerf.mark("sort")  # 200us-350us

        if num_intersects < 1:
            shape = (camera.h, camera.w)
            gaussian_ids_sorted = torch.zeros(0, 1, device=device)
            tile_bins = torch.zeros(0, 2, device=device)
            final_idx = torch.zeros(camera.h, camera.w, device=device)
            out_img = torch.zeros((*shape, 3)).float().to(device)
            out_alpha = torch.zeros((*shape, 1)).float().to(device)
            out_depth = torch.zeros((*shape, 1)).float().to(device)
            out_normal = torch.zeros((*shape, 3)).float().to(device)
            out_reg_depth = torch.zeros((*shape, 1)).float().to(device)
        else:
            assert colors.shape == torch.Size([num_points, 3])
            assert opacities.shape == torch.Size([num_points, 1])
            assert depth_ref.shape == torch.Size([camera.h, camera.w, 1])

            (
                final_idx, out_alpha,
                out_img, out_depth, out_normal,
                out_reg_depth,
            ) = _C.rasterize_forward(
                camera.h, camera.w, camera.intrins,
                depth_reg_pairwise_factor,
                gaussian_ids_sorted, tile_bins,
                positions, axes_u, axes_v,
                colors,
                ch_degree_r, ch_degree_r_to_use,
                ch_degree_phi, ch_degree_phi_to_use,
                ch_coeffs,
                opacities, #background,
                depth_ref,
            )
        timerf.mark("rasterize")  # 250us-600us

        ctx.camera = camera
        ctx.num_intersects = num_intersects
        ctx.depth_reg_pairwise_factor = depth_reg_pairwise_factor
        ctx.ch_degrees = (ch_degree_r, ch_degree_r_to_use,
                          ch_degree_phi, ch_degree_phi_to_use)
        ctx.save_for_backward(
            gaussian_ids_sorted, tile_bins,
            positions, axes_u, axes_v,
            colors, ch_coeffs, opacities, #background,
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

        camera = ctx.camera  # type: _Camera
        num_intersects = ctx.num_intersects
        ch_degrees = ctx.ch_degrees

        (
            gaussian_ids_sorted, tile_bins,
            positions, axes_u, axes_v,
            colors, ch_coeffs, opacities, #background,
            depth_ref,
            final_idx, out_alpha, out_depth, out_normal, out_reg_depth,
        ) = ctx.saved_tensors

        if num_intersects < 1:
            v_positions = torch.zeros_like(positions)
            v_positions_xy_abs = torch.zeros_like(positions)[..., :2]
            v_axes_u = torch.zeros_like(axes_u)
            v_axes_v = torch.zeros_like(axes_v)
            v_colors = torch.zeros_like(colors)
            v_ch_coeffs = torch.zeros_like(ch_coeffs)
            v_opacities = torch.zeros_like(opacities)
            v_depth_ref = torch.zeros_like(depth_ref)

        else:
            timerb.start()

            assert colors.shape[-1] == 3
            backward_return = _C.rasterize_backward(
                camera.h, camera.w, camera.intrins,
                *ch_degrees,
                ctx.depth_reg_pairwise_factor,
                gaussian_ids_sorted, tile_bins,
                positions, axes_u, axes_v,
                colors, ch_coeffs, opacities, #background,
                depth_ref,
                final_idx, out_alpha, out_depth,
                v_out_alpha, v_out_img, v_out_depth, v_out_normal,
                v_out_reg_depth.contiguous(),
            )
            timerb.mark("rasterize")  # 600us-2600us

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
            v_positions, v_axes_u, v_axes_v,
            v_colors, *([None]*4), v_ch_coeffs, v_opacities,
            v_depth_ref,
            # v_background,
            None, None, None, None,
        )
