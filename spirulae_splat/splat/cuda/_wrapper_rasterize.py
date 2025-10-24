# Modified from https://github.com/nerfstudio-project/gsplat/blob/2323de5905d5e90e035f792fe65bad0fedd413e7/gsplat/cuda/_wrapper.py

import math
import warnings
from dataclasses import dataclass
from enum import Enum
from typing import Any, Callable, Optional, Tuple

import torch
from torch import Tensor
from typing_extensions import Literal

import gsplat.cuda._wrapper

def _make_lazy_cuda_obj(name: str) -> Any:
    # pylint: disable=import-outside-toplevel
    from ._backend import _C

    obj = _C
    for name_split in name.split("."):
        obj = getattr(_C, name_split)
    return obj


def _make_lazy_cuda_func(name: str) -> Callable:
    from ._backend import _C

    def call_cuda(*args, **kwargs):
        return getattr(_C, name)(*args, **kwargs)

    return call_cuda



def rasterize_to_pixels(
    primitive: Literal["3dgs", "mip", "opaque_triangle"],
    splats: tuple[Tensor],  # means, quats, scales, opacities
    # colors: Tensor,  # [..., N, channels] or [nnz, channels]
    image_width: int,
    image_height: int,
    tile_size: int,
    isect_offsets: Tensor,  # [..., tile_height, tile_width]
    flatten_ids: Tensor,  # [n_isects]
    backgrounds: Optional[Tensor] = None,  # [..., channels]
    masks: Optional[Tensor] = None,  # [..., tile_height, tile_width]
    packed: bool = False,
    absgrad: bool = False,
) -> Tuple[Tensor, Tensor]:
    # TODO
    # image_dims = colors.shape[:-2]
    # channels = colors.shape[-1]
    # device = colors.device
    # if packed:
    #     nnz = colors.size(0)
    #     assert colors.shape[0] == nnz, colors.shape
    #     # assert means2d.shape == (nnz, 2), means2d.shape
    #     # assert conics.shape == (nnz, 3), conics.shape
    #     # assert opacities.shape == (nnz,), opacities.shape
    # else:
    #     N = colors.size(-2)
    #     assert colors.shape == image_dims + (N, channels), colors.shape
    #     # assert means2d.shape == image_dims + (N, 2), means2d.shape
    #     # assert conics.shape == image_dims + (N, 3), conics.shape
    #     # assert opacities.shape == image_dims + (N,), opacities.shape
    # if backgrounds is not None:
    #     assert backgrounds.shape == image_dims + (channels,), backgrounds.shape
    #     backgrounds = backgrounds.contiguous()
    # if masks is not None:
    #     assert masks.shape == isect_offsets.shape, masks.shape
    #     masks = masks.contiguous()

    # # Pad the channels to the nearest supported number if necessary
    # if channels > 513 or channels == 0:
    #     # TODO: maybe worth to support zero channels?
    #     raise ValueError(f"Unsupported number of color channels: {channels}")
    # if channels not in (1, 2, 3, 4, 5, 6, 7):
    #     padded_channels = (1 << (channels - 1).bit_length()) - channels
    #     colors = torch.cat(
    #         [
    #             colors,
    #             torch.zeros(*colors.shape[:-1], padded_channels, device=device),
    #         ],
    #         dim=-1,
    #     )
    #     if backgrounds is not None:
    #         backgrounds = torch.cat(
    #             [
    #                 backgrounds,
    #                 torch.zeros(
    #                     *backgrounds.shape[:-1], padded_channels, device=device
    #                 ),
    #             ],
    #             dim=-1,
    #         )
    # else:
    #     padded_channels = 0

    assert backgrounds is None, "TODO"

    tile_height, tile_width = isect_offsets.shape[-2:]
    assert (
        tile_height * tile_size >= image_height
    ), f"Assert Failed: {tile_height} * {tile_size} >= {image_height}"
    assert (
        tile_width * tile_size >= image_width
    ), f"Assert Failed: {tile_width} * {tile_size} >= {image_width}"

    if primitive in ["3dgs", "mip"]:
        _RasterizeToPixels = _RasterizeToPixels3DGS
    elif primitive in ["opaque_triangle"]:
        _RasterizeToPixels = _RasterizeToPixelsOpaqueTriangle
    render_outputs = _RasterizeToPixels.apply(
        *[(x.contiguous() if x is not None else None) for x in splats],
        # colors.contiguous(),
        backgrounds,
        masks,
        image_width,
        image_height,
        tile_size,
        isect_offsets.contiguous(),
        flatten_ids.contiguous(),
        absgrad,
    )

    # if padded_channels > 0:
    #     render_colors = render_colors[..., :-padded_channels]
    return tuple(render_outputs[:-1]), render_outputs[-1]



class _RasterizeToPixels3DGS(torch.autograd.Function):
    """Rasterize gaussians"""

    @staticmethod
    def forward(
        ctx,
        means2d: Tensor,  # [..., N, 2] or [nnz, 2]
        depths: Tensor,  # [..., N, 1] or [nnz, 1]
        conics: Tensor,  # [..., N, 3] or [nnz, 3]
        opacities: Tensor,  # [..., N] or [nnz]
        colors: Tensor,  # [..., N, channels] or [nnz, channels]
        backgrounds: Tensor,  # [..., channels], Optional
        masks: Tensor,  # [..., tile_height, tile_width], Optional
        width: int,
        height: int,
        tile_size: int,
        isect_offsets: Tensor,  # [..., tile_height, tile_width]
        flatten_ids: Tensor,  # [n_isects]
        absgrad: bool,
    ) -> Tuple[Tensor, Tensor]:
        (render_rgbs, render_depths), render_Ts, last_ids = _make_lazy_cuda_func(
            "rasterization_3dgs_forward"
        )(
            (means2d, depths, conics, opacities, colors), backgrounds, masks,
            width, height, tile_size, isect_offsets, flatten_ids,
        )

        ctx.save_for_backward(
            means2d, depths, conics, colors, opacities, backgrounds, masks,
            isect_offsets, flatten_ids, render_Ts, last_ids,
        )
        ctx.width = width
        ctx.height = height
        ctx.tile_size = tile_size
        ctx.absgrad = absgrad

        render_alphas = 1.0 - render_Ts
        return render_rgbs, render_depths, render_alphas

    @staticmethod
    def backward(
        ctx,
        v_render_rgbs: Tensor,  # [..., H, W, 3]
        v_render_depths: Tensor,  # [..., H, W, 3]
        v_render_alphas: Tensor,  # [..., H, W, 1]
    ):
        (
            means2d, depths, conics, colors, opacities, backgrounds, masks,
            isect_offsets, flatten_ids, render_Ts, last_ids,
        ) = ctx.saved_tensors
        width = ctx.width
        height = ctx.height
        tile_size = ctx.tile_size
        absgrad = ctx.absgrad

        (
            (v_means2d, v_depths, v_conics, v_opacities, v_colors),
            v_means2d_abs,
        ) = _make_lazy_cuda_func("rasterization_3dgs_backward")(
            (means2d, depths, conics, opacities, colors), backgrounds, masks,
            width, height, tile_size, isect_offsets, flatten_ids, render_Ts, last_ids,
            (v_render_rgbs.contiguous(), v_render_depths.contiguous()),
            v_render_alphas.contiguous(), absgrad,
        )
        if absgrad:
            means2d.absgrad = v_means2d_abs

        v_backgrounds = None
        if ctx.needs_input_grad[5]:
            v_backgrounds = (torch.cat([v_render_rgbs, v_render_depths], dim=-1) * \
                             render_Ts.float()).sum(dim=(-3, -2))

        return (
            v_means2d, v_depths, v_conics, v_opacities, v_colors, v_backgrounds,
            *([None]*(len(ctx.needs_input_grad)-6))
        )


class _RasterizeToPixelsOpaqueTriangle(torch.autograd.Function):
    """Rasterize opaque triangles"""

    @staticmethod
    def forward(
        ctx,
        means2d: Tensor,  # [..., N, 3, 2] or [nnz, 3, 2]
        hardness: Tensor,  # [..., N] or [nnz]
        colors: Tensor,  # [..., N, channels] or [nnz, channels]
        backgrounds: Tensor,  # [..., channels], Optional
        masks: Tensor,  # [..., tile_height, tile_width], Optional
        width: int,
        height: int,
        tile_size: int,
        isect_offsets: Tensor,  # [..., tile_height, tile_width]
        flatten_ids: Tensor,  # [n_isects]
        absgrad: bool,
    ) -> Tuple[Tensor, Tensor]:
        render_colors, render_Ts, last_ids = _make_lazy_cuda_func(
            "rasterization_opaque_triangle_forward"
        )(
            (means2d, hardness), colors, backgrounds, masks,
            width, height, tile_size, isect_offsets, flatten_ids,
        )

        ctx.save_for_backward(
            means2d, hardness, colors, backgrounds, masks,
            isect_offsets, flatten_ids, render_Ts, last_ids,
        )
        ctx.width = width
        ctx.height = height
        ctx.tile_size = tile_size
        ctx.absgrad = absgrad

        render_alphas = 1.0 - render_Ts
        return render_colors, render_alphas

    @staticmethod
    def backward(
        ctx,
        v_render_colors: Tensor,  # [..., H, W, 3]
        v_render_alphas: Tensor,  # [..., H, W, 1]
    ):
        (
            means2d, hardness, colors, backgrounds, masks,
            isect_offsets, flatten_ids, render_Ts, last_ids,
        ) = ctx.saved_tensors
        width = ctx.width
        height = ctx.height
        tile_size = ctx.tile_size
        absgrad = ctx.absgrad

        (
            (v_means2d, v_hardness),
            v_colors, v_means2d_abs,
        ) = _make_lazy_cuda_func("rasterization_opaque_triangle_backward")(
            (means2d, hardness), colors, backgrounds, masks,
            width, height, tile_size, isect_offsets, flatten_ids, render_Ts, last_ids,
            v_render_colors.contiguous(), v_render_alphas.contiguous(), absgrad,
        )
        if absgrad:
            means2d.absgrad = v_means2d_abs
            # means2d.absgrad = v_means2d.mean(-2)

        v_backgrounds = None
        if ctx.needs_input_grad[3]:
            v_backgrounds = (v_render_colors * render_Ts.float()).sum(dim=(-3, -2))

        return (
            v_means2d, v_hardness, v_colors, v_backgrounds,
            *([None]*(len(ctx.needs_input_grad)-4))
        )
