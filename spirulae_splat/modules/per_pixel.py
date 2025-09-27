import torch

from spirulae_splat.splat.cuda import _C


def blend_background(rgb, alpha, background):
    return _BlendBackground.apply(rgb, alpha, background)


class _BlendBackground(torch.autograd.Function):

    @staticmethod
    def forward(
        ctx,
        rgb, alpha, background
    ):
        out_rgb = _C.blend_background_forward(rgb, alpha, background)

        ctx.save_for_backward(rgb, alpha, background)

        return out_rgb

    @staticmethod
    def backward(ctx, v_out_rgb):

        rgb, alpha, background = ctx.saved_tensors

        return _C.blend_background_backward(
            rgb, alpha, background,
            v_out_rgb
        )
