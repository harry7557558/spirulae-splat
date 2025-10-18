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


def spherical_harmonics(
    degrees_to_use: int,
    dirs: Tensor,  # [..., 3]
    coeffs_dc: Tensor,  # [..., K, 3]
    coeffs_sh: Tensor  # [..., K, 3]
) -> Tensor:
    """Computes spherical harmonics.

    Args:
        degrees_to_use: The degree to be used.
        dirs: Directions. [..., 3]
        coeffs: Coefficients. [..., K, 3]
        masks: Optional boolen masks to skip some computation. [...,] Default: None.

    Returns:
        Spherical harmonics. [..., 3]
    """
    assert (degrees_to_use + 1) ** 2 - 1 <= coeffs_sh.shape[-2], coeffs_sh.shape
    batch_dims = dirs.shape[:-1]
    assert dirs.shape == batch_dims + (3,), dirs.shape
    assert coeffs_dc.shape == batch_dims + (3,), coeffs_dc.shape
    assert (
        (len(coeffs_sh.shape) == len(batch_dims) + 2)
        and coeffs_sh.shape[:-2] == batch_dims
        and coeffs_sh.shape[-1] == 3
    ), coeffs_sh.shape
    sh_degree = {
        0: 0,
        3: 1,
        8: 2,
        15: 3,
        24: 4
    }[coeffs_sh.shape[-2]]
    return _SphericalHarmonics.apply(
        sh_degree, degrees_to_use, dirs.contiguous(), coeffs_dc.contiguous(), coeffs_sh.contiguous()
    )



def fully_fused_projection(
    primitive: Literal["3dgs", "mip", "opaque_triangle"],
    splats: tuple[Tensor],  # means, quats, scales, opacities
    viewmats: Tensor,  # [..., C, 4, 4]
    Ks: Tensor,  # [..., C, 3, 3]
    width: int,
    height: int,
    eps2d: float = 0.3,
    near_plane: float = 0.01,
    far_plane: float = 1e10,
    radius_clip: float = 0.0,
    packed: bool = False,
    sparse_grad: bool = False,
    camera_model: Literal["pinhole", "ortho", "fisheye", "ftheta"] = "pinhole",
) -> Tuple[Tensor, Tensor, Tensor, Tensor, Tensor]:
    if primitive in ["3dgs", "mip"]:
        means, quats, scales, opacities = splats
        batch_dims = means.shape[:-2]
        N = means.shape[-2]
        assert means.shape == batch_dims + (N, 3), means.shape
        assert quats.shape == batch_dims + (N, 4), quats.shape
        assert scales.shape == batch_dims + (N, 3), scales.shape
        assert opacities.shape == batch_dims + (N,), opacities.shape
    C = viewmats.shape[-3]
    assert viewmats.shape == batch_dims + (C, 4, 4), viewmats.shape
    assert Ks.shape == batch_dims + (C, 3, 3), Ks.shape
    if sparse_grad:
        assert packed, "sparse_grad is only supported when packed is True"
        assert batch_dims == (), "sparse_grad does not support batch dimensions"

    assert (
        camera_model != "ftheta"
    ), "ftheta camera is only supported via UT, please set with_ut=True in the rasterization()"

    viewmats = viewmats.contiguous()
    Ks = Ks.contiguous()
    if packed:
        raise NotImplementedError("Packed not supported for fully_fused_projection without use_bvh")
    else:
        proj_return = _FullyFusedProjection3DGS.apply(
            *[x.contiguous() for x in splats],
            viewmats,
            Ks,
            width,
            height,
            eps2d,
            near_plane,
            far_plane,
            radius_clip,
            primitive == "mip",
            camera_model,
        )
        return (
            proj_return[0],  # radii
            proj_return[1],  # depth
            tuple(proj_return[2:])
        )


def fully_fused_projection_hetero(
    means: Tensor,  # [..., N, 3]
    quats: Tensor,  # [..., N, 4]
    scales: Tensor,  # [..., N, 3]
    viewmats: Tensor,  # [..., C, 4, 4]
    Ks: Tensor,  # [..., C, 3, 3]
    width: int,
    height: int,
    intersection_count_map: Tensor,  # [C+1]
    intersection_splat_id: Tensor,  # [nnz]
    eps2d: float = 0.3,
    near_plane: float = 0.01,
    far_plane: float = 1e10,
    radius_clip: float = 0.0,
    sparse_grad: bool = False,
    calc_compensations: bool = False,
    camera_model: Literal["pinhole", "ortho", "fisheye", "ftheta"] = "pinhole",
    opacities: Optional[Tensor] = None,  # [..., N] or None
) -> Tuple[Tensor, Tensor, Tensor, Tensor, Tensor]:
    batch_dims = means.shape[:-2]
    N = means.shape[-2]
    C = viewmats.shape[-3]
    assert means.shape == batch_dims + (N, 3), means.shape
    assert viewmats.shape == batch_dims + (C, 4, 4), viewmats.shape
    assert Ks.shape == batch_dims + (C, 3, 3), Ks.shape
    means = means.contiguous()
    assert quats.shape == batch_dims + (N, 4), quats.shape
    assert scales.shape == batch_dims + (N, 3), scales.shape
    quats = quats.contiguous()
    scales = scales.contiguous()
    if sparse_grad:
        assert batch_dims == (), "sparse_grad does not support batch dimensions"
    if opacities is not None:
        assert opacities.shape == batch_dims + (N,), opacities.shape
        opacities = opacities.contiguous()

    assert (
        camera_model != "ftheta"
    ), "ftheta camera is only supported via UT, please set with_ut=True in the rasterization()"

    viewmats = viewmats.contiguous()
    Ks = Ks.contiguous()
    return _FullyFusedProjectionHetero.apply(
        means,
        quats,
        scales,
        opacities,
        viewmats,
        Ks,
        width,
        height,
        eps2d,
        near_plane,
        far_plane,
        radius_clip,
        sparse_grad,
        calc_compensations,
        camera_model,
        intersection_count_map,
        intersection_splat_id,
    )


def rasterize_to_pixels(
    primitive: Literal["3dgs", "mip", "opaque_triangle"],
    splats: tuple[Tensor],  # means, quats, scales, opacities
    colors: Tensor,  # [..., N, channels] or [nnz, channels]
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
    image_dims = colors.shape[:-2]
    channels = colors.shape[-1]
    device = colors.device
    if packed:
        nnz = colors.size(0)
        assert colors.shape[0] == nnz, colors.shape
        # assert means2d.shape == (nnz, 2), means2d.shape
        # assert conics.shape == (nnz, 3), conics.shape
        # assert opacities.shape == (nnz,), opacities.shape
    else:
        N = colors.size(-2)
        assert colors.shape == image_dims + (N, channels), colors.shape
        # assert means2d.shape == image_dims + (N, 2), means2d.shape
        # assert conics.shape == image_dims + (N, 3), conics.shape
        # assert opacities.shape == image_dims + (N,), opacities.shape
    if backgrounds is not None:
        assert backgrounds.shape == image_dims + (channels,), backgrounds.shape
        backgrounds = backgrounds.contiguous()
    if masks is not None:
        assert masks.shape == isect_offsets.shape, masks.shape
        masks = masks.contiguous()

    # Pad the channels to the nearest supported number if necessary
    if channels > 513 or channels == 0:
        # TODO: maybe worth to support zero channels?
        raise ValueError(f"Unsupported number of color channels: {channels}")
    if channels not in (1, 2, 3, 4, 5):
        padded_channels = (1 << (channels - 1).bit_length()) - channels
        colors = torch.cat(
            [
                colors,
                torch.zeros(*colors.shape[:-1], padded_channels, device=device),
            ],
            dim=-1,
        )
        if backgrounds is not None:
            backgrounds = torch.cat(
                [
                    backgrounds,
                    torch.zeros(
                        *backgrounds.shape[:-1], padded_channels, device=device
                    ),
                ],
                dim=-1,
            )
    else:
        padded_channels = 0

    tile_height, tile_width = isect_offsets.shape[-2:]
    assert (
        tile_height * tile_size >= image_height
    ), f"Assert Failed: {tile_height} * {tile_size} >= {image_height}"
    assert (
        tile_width * tile_size >= image_width
    ), f"Assert Failed: {tile_width} * {tile_size} >= {image_width}"

    render_colors, render_alphas = _RasterizeToPixels3DGS.apply(
        *[x.contiguous() for x in splats],
        colors.contiguous(),
        backgrounds,
        masks,
        image_width,
        image_height,
        tile_size,
        isect_offsets.contiguous(),
        flatten_ids.contiguous(),
        absgrad,
    )

    if padded_channels > 0:
        render_colors = render_colors[..., :-padded_channels]
    return render_colors, render_alphas



class _SphericalHarmonics(torch.autograd.Function):
    """Spherical Harmonics"""

    @staticmethod
    def forward(
        ctx, sh_degree: int, degrees_to_use: int, dirs: Tensor, coeffs_dc: Tensor, coeffs_sh: Tensor
    ) -> Tensor:
        colors = _make_lazy_cuda_func("compute_sh_forward")(
            sh_degree, degrees_to_use,
            dirs, coeffs_dc, coeffs_sh
        )
        ctx.save_for_backward(dirs, coeffs_sh, colors)
        ctx.sh_degree = (sh_degree, degrees_to_use)
        return colors

    @staticmethod
    def backward(ctx, v_colors: Tensor):
        dirs, coeffs_sh, colors = ctx.saved_tensors
        (sh_degree, degrees_to_use) = ctx.sh_degree
        compute_v_dirs = ctx.needs_input_grad[1]  # TODO
        v_coeffs_dc, v_coeffs_sh, v_dirs = _make_lazy_cuda_func("compute_sh_backward")(
            sh_degree, degrees_to_use,
            dirs, coeffs_sh, colors, v_colors.contiguous()
        )
        return None, None, v_dirs, v_coeffs_dc, v_coeffs_sh


class _FullyFusedProjection3DGS(torch.autograd.Function):
    """Projects Gaussians to 2D."""

    @staticmethod
    def forward(
        ctx,
        means: Tensor,  # [..., N, 3]
        quats: Tensor,  # [..., N, 4]
        scales: Tensor,  # [..., N, 3]
        opacities: Tensor,  # [..., N]
        viewmats: Tensor,  # [..., C, 4, 4]
        Ks: Tensor,  # [..., C, 3, 3]
        width: int,
        height: int,
        eps2d: float,
        near_plane: float,
        far_plane: float,
        radius_clip: float,
        is_antialiased: bool,
        camera_model: Literal["pinhole", "ortho", "fisheye", "ftheta"] = "pinhole",
    ) -> Tuple[Tensor, Tensor, Tensor, Tensor, Tensor]:
        assert (
            camera_model != "ftheta"
        ), "ftheta camera is only supported via UT, please set with_ut=True in the rasterization()"

        camera_model_type = gsplat.cuda._wrapper._make_lazy_cuda_obj(
            f"CameraModelType.{camera_model.upper()}"
        )

        # "covars" and {"quats", "scales"} are mutually exclusive
        radii, depths, (means2d, conics, proj_opacities) = _make_lazy_cuda_func(
            "projection_ewa_3dgs_forward"
        )(
            is_antialiased,
            (means, quats, scales, opacities),
            viewmats,
            Ks,
            width,
            height,
            eps2d,
            near_plane,
            far_plane,
            radius_clip,
            camera_model_type,
        )
        ctx.save_for_backward(
            means, quats, scales, opacities, viewmats, Ks, radii
        )
        ctx.is_antialiased = is_antialiased
        ctx.width = width
        ctx.height = height
        ctx.eps2d = eps2d
        ctx.camera_model_type = camera_model_type

        return radii, depths, means2d, conics, proj_opacities

    @staticmethod
    def backward(ctx, v_radii, v_depths, v_means2d, v_conics, v_opacities):
        (
            means,
            quats,
            scales,
            opacities,
            viewmats,
            Ks,
            radii,
        ) = ctx.saved_tensors
        (v_means, v_quats, v_scales, v_opacities), v_viewmats = _make_lazy_cuda_func(
            "projection_ewa_3dgs_backward"
        )(
            ctx.is_antialiased,
            (means, quats, scales, opacities),
            viewmats,
            Ks,
            ctx.width,
            ctx.height,
            ctx.eps2d,
            ctx.camera_model_type,
            radii,
            v_depths.contiguous(),
            (v_means2d.contiguous(), v_conics.contiguous(), v_opacities.contiguous()),
            ctx.needs_input_grad[4],  # viewmats_requires_grad
        )
        if not ctx.needs_input_grad[4]:
            v_viewmats = None
        return (
            v_means,
            v_quats,
            v_scales,
            v_opacities,
            v_viewmats,
            *([None]*(len(ctx.needs_input_grad)-5))
        )


class _FullyFusedProjectionHetero(torch.autograd.Function):
    """Projects Gaussians to 2D. Return packed tensors."""

    @staticmethod
    def forward(
        ctx,
        means: Tensor,  # [..., N, 3]
        quats: Tensor,  # [..., N, 4]
        scales: Tensor,  # [..., N, 3]
        opacities: Tensor,  # [..., N]
        viewmats: Tensor,  # [..., C, 4, 4]
        Ks: Tensor,  # [..., C, 3, 3]
        width: int,
        height: int,
        eps2d: float,
        near_plane: float,
        far_plane: float,
        radius_clip: float,
        sparse_grad: bool,
        calc_compensations: bool,
        camera_model: Literal["pinhole", "ortho", "fisheye", "ftheta"],
        intersection_count_map: Tensor,  # [C+1]
        intersection_splat_id: Tensor,  # [nnz]
    ) -> Tuple[Tensor, Tensor, Tensor, Tensor]:
        assert (
            camera_model != "ftheta"
        ), "ftheta camera is only supported via UT, please set with_ut=True in the rasterization()"

        camera_model_type = gsplat.cuda._wrapper._make_lazy_cuda_obj(
            f"CameraModelType.{camera_model.upper()}"
        )

        nnz = len(intersection_splat_id)

        (
            camera_ids,
            gaussian_ids,
            radii,
            means2d,
            depths,
            conics,
            compensations,
        ) = _make_lazy_cuda_func("projection_ewa_3dgs_hetero_forward")(
            means,
            quats,
            scales,
            opacities,
            viewmats,
            Ks,
            width,
            height,
            eps2d,
            near_plane,
            far_plane,
            radius_clip,
            calc_compensations,
            camera_model_type,
            intersection_count_map,
            intersection_splat_id
        )
        if not calc_compensations:
            compensations = None
        ctx.save_for_backward(
            camera_ids,
            gaussian_ids,
            means,
            quats,
            scales,
            viewmats,
            Ks,
            conics,
            compensations,
        )
        ctx.width = width
        ctx.height = height
        ctx.eps2d = eps2d
        ctx.sparse_grad = sparse_grad
        ctx.camera_model_type = camera_model_type

        return (
            camera_ids,
            gaussian_ids,
            radii,
            means2d,
            depths,
            conics,
            compensations,
        )

    @staticmethod
    def backward(
        ctx,
        v_camera_ids,
        v_gaussian_ids,
        v_radii,
        v_means2d,
        v_depths,
        v_conics,
        v_compensations,
    ):
        (
            camera_ids,
            gaussian_ids,
            means,
            quats,
            scales,
            viewmats,
            Ks,
            conics,
            compensations,
        ) = ctx.saved_tensors
        width = ctx.width
        height = ctx.height
        eps2d = ctx.eps2d
        sparse_grad = ctx.sparse_grad
        camera_model_type = ctx.camera_model_type

        if v_compensations is not None:
            v_compensations = v_compensations.contiguous()
        v_means, v_quats, v_scales, v_viewmats = _make_lazy_cuda_func(
            "projection_ewa_3dgs_hetero_backward"
        )(
            means,
            quats,
            scales,
            viewmats,
            Ks,
            width,
            height,
            eps2d,
            camera_model_type,
            camera_ids,
            gaussian_ids,
            conics,
            compensations,
            v_means2d.contiguous(),
            v_depths.contiguous(),
            v_conics.contiguous(),
            v_compensations,
            ctx.needs_input_grad[4],  # viewmats_requires_grad
            sparse_grad,
        )

        if sparse_grad:
            batch_dims = means.shape[:-2]
            B = math.prod(batch_dims)
            N = means.shape[-2]
        if not ctx.needs_input_grad[0]:
            v_means = None
        else:
            if sparse_grad:
                # TODO: gaussian_ids is duplicated so not ideal.
                # An idea is to directly set the attribute (e.g., .sparse_grad) of
                # the tensor but this requires the tensor to be leaf node only. And
                # a customized optimizer would be needed in this case.
                v_means = torch.sparse_coo_tensor(
                    indices=gaussian_ids[None],
                    values=v_means,  # [nnz, 3]
                    size=means.shape,
                    is_coalesced=len(viewmats) == 1,
                )
        if not ctx.needs_input_grad[1]:
            v_quats = None
        else:
            if sparse_grad:
                v_quats = torch.sparse_coo_tensor(
                    indices=gaussian_ids[None],
                    values=v_quats,  # [nnz, 4]
                    size=quats.shape,
                    is_coalesced=len(viewmats) == 1,
                )
        if not ctx.needs_input_grad[2]:
            v_scales = None
        else:
            if sparse_grad:
                v_scales = torch.sparse_coo_tensor(
                    indices=gaussian_ids[None],
                    values=v_scales,  # [nnz, 3]
                    size=scales.shape,
                    is_coalesced=len(viewmats) == 1,
                )
        if not ctx.needs_input_grad[4]:
            v_viewmats = None

        return (
            v_means,
            v_quats,
            v_scales,
            None,
            v_viewmats,
            *([None]*(len(ctx.needs_input_grad)-5))
        )


class _RasterizeToPixels3DGS(torch.autograd.Function):
    """Rasterize gaussians"""

    @staticmethod
    def forward(
        ctx,
        means2d: Tensor,  # [..., N, 2] or [nnz, 2]
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
        render_colors, render_Ts, last_ids = _make_lazy_cuda_func(
            "rasterization_3dgs_forward"
        )(
            (means2d, conics, opacities),
            colors,
            backgrounds,
            masks,
            width,
            height,
            tile_size,
            isect_offsets,
            flatten_ids,
        )

        ctx.save_for_backward(
            means2d,
            conics,
            colors,
            opacities,
            backgrounds,
            masks,
            isect_offsets,
            flatten_ids,
            render_Ts,
            last_ids,
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
            means2d,
            conics,
            colors,
            opacities,
            backgrounds,
            masks,
            isect_offsets,
            flatten_ids,
            render_Ts,
            last_ids,
        ) = ctx.saved_tensors
        width = ctx.width
        height = ctx.height
        tile_size = ctx.tile_size
        absgrad = ctx.absgrad

        # print(last_ids)
        # print(render_Ts)
        # print(tile_size, flatten_ids.shape)

        (
            (v_means2d, v_conics, v_opacities),
            v_colors,
            v_means2d_abs,
        ) = _make_lazy_cuda_func("rasterization_3dgs_backward")(
            (means2d, conics, opacities),
            colors,
            backgrounds,
            masks,
            width,
            height,
            tile_size,
            isect_offsets,
            flatten_ids,
            render_Ts,
            last_ids,
            v_render_colors.contiguous(),
            v_render_alphas.contiguous(),
            absgrad,
        )

        if absgrad:
            means2d.absgrad = v_means2d_abs

        if ctx.needs_input_grad[4]:
            v_backgrounds = (v_render_colors * render_Ts.float()).sum(
                dim=(-3, -2)
            )
        else:
            v_backgrounds = None

        return (
            v_means2d,
            v_conics,
            v_opacities,
            v_colors,
            v_backgrounds,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
        )


intersect_splat_tile = _make_lazy_cuda_func("intersect_splat_tile")
