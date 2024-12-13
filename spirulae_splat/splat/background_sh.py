"""Python bindings for custom CUDA functions"""

import torch
from jaxtyping import Float, Int
from typing import Tuple
from torch import Tensor
from torch.autograd import Function

import spirulae_splat.splat.cuda as _C


def render_background_sh(
    w: int,
    h: int,
    intrins: Tuple[float, float, float, float],
    rotation: Float[Tensor, "3 3"],
    sh_degree: int,
    sh_coeffs: Float[Tensor, "sh_degree**2 3"],
    block_width: int
) -> Float[Tensor, "h w 3"]:
    assert block_width > 1 and block_width <= 16, "block_width must be between 2 and 16"

    return _RenderBackgroundSH.apply(
        w, h, intrins,
        rotation.contiguous(),
        sh_degree, sh_coeffs.contiguous(),
        block_width
    )


class _RenderBackgroundSH(Function):

    @staticmethod
    def forward(
        ctx,
        w: int,
        h: int,
        intrins: Tuple[float, float, float, float],
        rotation: Float[Tensor, "3 3"],
        sh_degree: int,
        sh_coeffs: Float[Tensor, "sh_degree**2 3"],
        block_width: int
    ) -> Tuple[Tensor]:

        out_color = _C.render_background_sh_forward(
            w, h, block_width, *intrins,
            rotation, sh_degree, sh_coeffs
        )

        ctx.w = w
        ctx.h = h
        ctx.block_width = block_width
        ctx.intrins = intrins
        ctx.sh_degree = sh_degree
        ctx.save_for_backward(rotation, sh_coeffs, out_color)

        return out_color

    @staticmethod
    def backward(ctx, v_out_color):

        w, h = ctx.w, ctx.h
        block_width = ctx.block_width
        intrins = ctx.intrins
        sh_degree = ctx.sh_degree
        rotation, sh_coeffs, out_color = ctx.saved_tensors

        v_rotation, v_sh_coeffs = _C.render_background_sh_backward(
            w, h, block_width, *intrins,
            rotation, sh_degree, sh_coeffs,
            out_color, v_out_color.contiguous()
        )

        return (
            None, None, None,
            v_rotation,
            None, v_sh_coeffs,
            None
        )
