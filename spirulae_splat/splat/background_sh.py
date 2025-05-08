"""Python bindings for custom CUDA functions"""

import torch
from jaxtyping import Float, Int
from typing import Tuple
from torch import Tensor
from torch.autograd import Function

import spirulae_splat.splat.cuda as _C
from spirulae_splat.splat._camera import _Camera


def render_background_sh(
    camera: _Camera,
    rotation: Float[Tensor, "3 3"],
    sh_degree: int,
    sh_coeffs: Float[Tensor, "(sh_degree+1)**2 3"],
) -> Float[Tensor, "h w 3"]:

    return _RenderBackgroundSH.apply(
        camera,
        rotation.contiguous(),
        sh_degree, sh_coeffs.contiguous(),
    )


class _RenderBackgroundSH(Function):

    @staticmethod
    def forward(
        ctx,
        camera: _Camera,
        rotation: Float[Tensor, "3 3"],
        sh_degree: int,
        sh_coeffs: Float[Tensor, "(sh_degree+1)**2 3"],
    ) -> Tuple[Tensor]:

        out_color = _C.render_background_sh_forward(
            camera.w, camera.h,
            camera.model, camera.intrins, camera.get_undist_map(),
            rotation, sh_degree+1, sh_coeffs
        )

        ctx.camera = camera
        ctx.sh_degree = sh_degree
        ctx.save_for_backward(rotation, sh_coeffs, out_color)

        return out_color

    @staticmethod
    def backward(ctx, v_out_color):

        camera = ctx.camera  # type: _Camera

        sh_degree = ctx.sh_degree
        rotation, sh_coeffs, out_color = ctx.saved_tensors

        v_rotation, v_sh_coeffs = _C.render_background_sh_backward(
            camera.w, camera.h,
            camera.model, camera.intrins, camera.get_undist_map(),
            rotation, sh_degree+1, sh_coeffs,
            out_color, v_out_color.contiguous()
        )

        return (
            None,
            v_rotation,
            None, v_sh_coeffs,
        )
