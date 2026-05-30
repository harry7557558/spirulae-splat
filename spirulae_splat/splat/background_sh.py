"""Python bindings for the C++ background-SH (skybox) renderer.

The C++ side fills caller-allocated device tensors via TorchTensorView; we
allocate the output buffers here on the matching device and convert them to
the (data_ptr, element_size, shape) tuple form the engine expects.

The C++ entry points take a world->camera viewmat (row-major 4x4), per-camera
intrinsics, and optional distortion coeffs + cam_indices. This thin Python
wrapper accepts a camera-to-world rotation (the legacy API) and builds a
viewmat from it (rotation^T as the upper-left 3x3, zero translation), which is
sufficient when the caller only needs the skybox in a fixed world frame.
"""

import torch
from jaxtyping import Float
from typing import Literal
from torch import Tensor
from torch.autograd import Function

import spirulae_splat.splat.cuda as _C


def _tv(t: Tensor):
    """Convert a contiguous torch.Tensor to TorchTensorView tuple."""
    assert t.is_contiguous(), f"Tensor must be contiguous, got strides {t.stride()}"
    return (t.data_ptr(), t.element_size(), list(t.shape))


def _null_tv():
    return (0, 4, [0])


def _build_viewmat_from_rotation_c2w(rotation_c2w: Float[Tensor, "B 3 3"]) -> Tensor:
    """Build a [B, 4, 4] row-major world->camera viewmat from a camera-to-world
    rotation. Only the rotation block is meaningful for sky rendering.

    The engine's viewmat convention bakes a Y/Z-flip (diag(1, -1, -1)) into the
    upper-left 3x3 before transposing, to convert NeRF-style camera axes (Y up,
    Z behind) to the slang/gsplat convention (Y down, Z forward) the kernel
    expects. We apply the same here so a caller passing a vanilla c2w rotation
    gets the right skybox orientation.
    """
    B = rotation_c2w.shape[0]
    flip = torch.tensor([[1.0, -1.0, -1.0]], device=rotation_c2w.device,
                        dtype=rotation_c2w.dtype)
    R = rotation_c2w * flip                       # NeRF -> gsplat axes
    R_inv = R.transpose(-1, -2)                    # camera-to-world -> world-to-camera
    vm = torch.zeros((B, 4, 4), dtype=rotation_c2w.dtype, device=rotation_c2w.device)
    vm[:, :3, :3] = R_inv
    vm[:, 3, 3] = 1.0
    return vm.contiguous()


def render_background_sh(
    width: int,
    height: int,
    camera_model: Literal['pinhole', 'fisheye'],
    intrins: Float[Tensor, "B 4"],
    rotation: Float[Tensor, "B 3 3"],
    sh_degree: int,
    sh_coeffs: Float[Tensor, "(sh_degree+1)**2 3"],
) -> Float[Tensor, "B H W 3"]:

    return _RenderBackgroundSH.apply(
        width, height, camera_model,
        intrins.contiguous(), rotation.contiguous(),
        sh_degree, sh_coeffs.contiguous(),
    )


class _RenderBackgroundSH(Function):

    @staticmethod
    def forward(
        ctx,
        width: int,
        height: int,
        camera_model: Literal['pinhole', 'fisheye'],
        intrins: Float[Tensor, "B 4"],
        rotation: Float[Tensor, "B 3 3"],
        sh_degree: int,
        sh_coeffs: Float[Tensor, "(sh_degree+1)**2 3"],
    ) -> Tensor:
        b = intrins.numel() // 4
        out_color = torch.empty(
            (b, height, width, 3),
            dtype=sh_coeffs.dtype, device=sh_coeffs.device,
        ).contiguous()

        viewmat = _build_viewmat_from_rotation_c2w(rotation)
        _C.render_background_sh_forward(
            int(width), int(height), str(camera_model).upper(),
            int(sh_degree),
            _tv(viewmat), _tv(intrins), _null_tv(),
            _tv(sh_coeffs),
            _tv(out_color),
        )

        ctx.meta = (width, height, camera_model, sh_degree)
        ctx.save_for_backward(intrins, viewmat, sh_coeffs, out_color)
        return out_color

    @staticmethod
    def backward(ctx, v_out_color):
        width, height, camera_model, sh_degree = ctx.meta
        intrins, viewmat, sh_coeffs, out_color = ctx.saved_tensors

        v_sh_coeffs = torch.zeros_like(sh_coeffs).contiguous()
        v_out_color = v_out_color.contiguous()

        _C.render_background_sh_backward(
            int(width), int(height), str(camera_model).upper(),
            int(sh_degree),
            _tv(viewmat), _tv(intrins), _null_tv(),
            _tv(sh_coeffs),
            _tv(out_color), _tv(v_out_color),
            _tv(v_sh_coeffs),
        )

        v_sh_coeffs = torch.nan_to_num(torch.clip(v_sh_coeffs, -40.0, 40.0))
        return (
            None, None, None,
            None, None,            # intrins, rotation: no grad
            None, v_sh_coeffs,
        )
