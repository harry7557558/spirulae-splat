"""Python bindings for custom CUDA functions"""

import torch
from jaxtyping import Float, Int
from typing import Tuple, Literal
from torch import Tensor
from torch.autograd import Function

import spirulae_splat.splat.cuda as _C


def render_background_sh(
    width: int,
    height: int,
    camera_model: Literal['pinhole', 'fisheye'],
    Ks: Float[Tensor, "3 3"],
    rotation: Float[Tensor, "3 3"],
    sh_degree: int,
    sh_coeffs: Float[Tensor, "(sh_degree+1)**2 3"],
) -> Float[Tensor, "h w 3"]:

    return _RenderBackgroundSH.apply(
        width, height, camera_model,
        Ks.contiguous(), rotation.contiguous(),
        sh_degree, sh_coeffs.contiguous(),
    )


class _RenderBackgroundSH(Function):

    @staticmethod
    def forward(
        ctx,
        width: int,
        height: int,
        camera_model: Literal['pinhole', 'fisheye'],
        Ks: Float[Tensor, "3 3"],
        rotation: Float[Tensor, "3 3"],
        sh_degree: int,
        sh_coeffs: Float[Tensor, "(sh_degree+1)**2 3"],
    ) -> Tuple[Tensor]:
        
        out_color = _C.render_background_sh_forward(
            width, height, camera_model,
            Ks, rotation,
            sh_degree+1, sh_coeffs
        )

        ctx.meta = (width, height, camera_model, sh_degree)
        ctx.save_for_backward(Ks, rotation, sh_coeffs, out_color)

        return out_color

    @staticmethod
    def backward(ctx, v_out_color):

        (width, height, camera_model, sh_degree) = ctx.meta
        Ks, rotation, sh_coeffs, out_color = ctx.saved_tensors

        v_rotation, v_sh_coeffs = _C.render_background_sh_backward(
            width, height, camera_model,
            Ks, rotation,
            sh_degree+1, sh_coeffs,
            out_color, v_out_color.contiguous()
        )

        clean = lambda x, h: torch.nan_to_num(torch.clip(x, -h, h))
        v_rotation = clean(v_rotation, 40.)
        v_sh_coeffs = clean(v_sh_coeffs, 40.)

        return (
            None, None, None,
            None, v_rotation,
            None, v_sh_coeffs,
        )
