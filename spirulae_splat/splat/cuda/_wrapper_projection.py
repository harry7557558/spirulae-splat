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



def fully_fused_projection(
    primitive: Literal["3dgs", "mip", "opaque_triangle"],
    splats: tuple[Tensor],  # means, quats, scales, opacities
    viewmats: Tensor,  # [..., C, 4, 4]
    Ks: Tensor,  # [..., C, 3, 3]
    width: int,
    height: int,
    near_plane: float = 0.01,
    far_plane: float = 1e10,
    packed: bool = False,
    sparse_grad: bool = False,
    camera_model: Literal["pinhole", "ortho", "fisheye", "ftheta"] = "pinhole",
    radial_coeffs: Optional[Tensor] = None,
    tangential_coeffs: Optional[Tensor] = None,
    thin_prism_coeffs: Optional[Tensor] = None,
) -> Tuple[Tensor, Tensor, Optional[Tensor], Tuple[Tensor]]:
    # if primitive in ["3dgs", "mip"]:
    if primitive in ["3dgs", "mip", "opaque_triangle"]:
        means, quats, scales, opacities, features_dc, features_sh = splats
        batch_dims = means.shape[:-2]
        N = means.shape[-2]
        assert means.shape == batch_dims + (N, 3), means.shape
        assert quats.shape == batch_dims + (N, 4), quats.shape
        assert scales.shape == batch_dims + (N, 3), scales.shape
        if primitive in ["3dgs", "mip"]:
            assert opacities.shape == batch_dims + (N,), opacities.shape
        else:
            assert opacities.shape == batch_dims + (N, 2), opacities.shape
        assert features_dc.shape == batch_dims + (N, 3), features_dc.shape
    elif primitive in ["opaque_triangle"]:
        means, hardness = splats
        batch_dims = means.shape[:-3]
        N = means.shape[-3]
        assert means.shape == batch_dims + (N, 3, 3), means.shape
        assert hardness.shape == batch_dims + (N,), hardness.shape
    C = viewmats.shape[-3]
    assert viewmats.shape == batch_dims + (C, 4, 4), viewmats.shape
    assert Ks.shape == batch_dims + (C, 3, 3), Ks.shape
    if sparse_grad:
        assert packed, "sparse_grad is only supported when packed is True"
        assert batch_dims == (), "sparse_grad does not support batch dimensions"

    if camera_model != "fisheye":
        assert radial_coeffs is None and tangential_coeffs is None and thin_prism_coeffs is None, \
            "Camera distortion is only supported for fisheye model"
    if radial_coeffs is not None:
        assert radial_coeffs.shape == batch_dims + (C, 4), radial_coeffs.shape
        radial_coeffs = radial_coeffs.contiguous().to(viewmats)
    if tangential_coeffs is not None:
        assert tangential_coeffs.shape == batch_dims + (C, 2), tangential_coeffs.shape
        tangential_coeffs = tangential_coeffs.contiguous().to(viewmats)
    if thin_prism_coeffs is not None:
        assert thin_prism_coeffs.shape == batch_dims + (C, 2), thin_prism_coeffs.shape
        thin_prism_coeffs = thin_prism_coeffs.contiguous().to(viewmats)

    assert (
        camera_model != "ftheta"
    ), "ftheta camera is only supported via UT, please set with_ut=True in the rasterization()"

    if packed:
        raise NotImplementedError("Packed not supported for fully_fused_projection without use_bvh")
    else:
        additional_args = []
        if primitive in ["3dgs", "mip"]:
            _RasterizeToPixels = _FullyFusedProjection3DGS
            additional_args = [primitive == "mip"]
        elif primitive in ["opaque_triangle"]:
            _RasterizeToPixels = _FullyFusedProjectionOpaqueTriangle
        proj_return = _RasterizeToPixels.apply(
            *[x.contiguous() for x in [*splats] + [viewmats, Ks]],
            width,
            height,
            near_plane,
            far_plane,
            camera_model,
            (radial_coeffs, tangential_coeffs, thin_prism_coeffs),
            *additional_args
        )
        return (
            proj_return[0],  # aabb
            proj_return[2],  # depth
            tuple(proj_return[1:])
        )

class _FullyFusedProjection3DGS(torch.autograd.Function):
    """Projects Gaussians to 2D."""

    @staticmethod
    def forward(
        ctx,
        means: Tensor,  # [..., N, 3]
        quats: Tensor,  # [..., N, 4]
        scales: Tensor,  # [..., N, 3]
        opacities: Tensor,  # [..., N]
        features_dc: Tensor,  # [..., N, 3]
        features_sh: Tensor,  # [..., N, x, 3]
        viewmats: Tensor,  # [..., C, 4, 4]
        Ks: Tensor,  # [..., C, 3, 3]
        width: int,
        height: int,
        near_plane: float,
        far_plane: float,
        camera_model: Literal["pinhole", "ortho", "fisheye", "ftheta"],
        dist_coeffs: Tuple[Tensor, Tensor, Tensor],
        is_antialiased: bool,
    ) -> Tuple[Tensor, Tensor, Tensor, Tensor, Tensor]:
        assert (
            camera_model != "ftheta"
        ), "ftheta camera is only supported via UT, please set with_ut=True in the rasterization()"

        camera_model_type = gsplat.cuda._wrapper._make_lazy_cuda_obj(
            f"CameraModelType.{camera_model.upper()}"
        )

        aabb, (means2d, depths, conics, proj_opacities, proj_rgbs) = _make_lazy_cuda_func(
            "projection_ewa_3dgs_forward"
        )(
            is_antialiased,
            (means, quats, scales, opacities, features_dc, features_sh),
            viewmats, Ks, width, height,
            near_plane, far_plane,
            camera_model_type, dist_coeffs
        )
        ctx.save_for_backward(means, quats, scales, opacities, features_dc, features_sh, viewmats, Ks, aabb)
        ctx.is_antialiased = is_antialiased
        ctx.width = width
        ctx.height = height
        ctx.camera_model_type = camera_model_type
        ctx.dist_coeffs = dist_coeffs

        return aabb, means2d, depths, conics, proj_opacities, proj_rgbs

    @staticmethod
    def backward(ctx, v_aabb, v_means2d, v_depths, v_conics, v_opacities, v_proj_rgbs):
        means, quats, scales, opacities, features_dc, features_sh, viewmats, Ks, aabb = ctx.saved_tensors
        (v_means, v_quats, v_scales, v_opacities, v_features_dc, v_features_sh), v_viewmats = _make_lazy_cuda_func(
            "projection_ewa_3dgs_backward"
        )(
            ctx.is_antialiased,
            (means, quats, scales, opacities, features_dc, features_sh),
            viewmats, Ks, ctx.width, ctx.height, ctx.camera_model_type, ctx.dist_coeffs, aabb,
            tuple([x.contiguous() for x in (v_means2d, v_depths, v_conics, v_opacities, v_proj_rgbs)]),
            ctx.needs_input_grad[6],  # viewmats_requires_grad
        )
        if not ctx.needs_input_grad[6]:
            v_viewmats = None
        return (
            v_means, v_quats, v_scales, v_opacities, v_features_dc, v_features_sh, v_viewmats,
            *([None]*(len(ctx.needs_input_grad)-7))
        )

class _FullyFusedProjectionOpaqueTriangle(torch.autograd.Function):
    """Projects Gaussians to 2D."""

    @staticmethod
    def forward(
        ctx,
        means: Tensor,  # [..., N, 3]
        quats: Tensor,  # [..., N, 4]
        scales: Tensor,  # [..., N, 3]
        hardness: Tensor,  # [..., N]
        features_dc: Tensor,  # [..., N, 3]
        features_sh: Tensor,  # [..., N, x, 3]
        viewmats: Tensor,  # [..., C, 4, 4]
        Ks: Tensor,  # [..., C, 3, 3]
        width: int,
        height: int,
        near_plane: float,
        far_plane: float,
        camera_model: Literal["pinhole", "ortho", "fisheye", "ftheta"],
        dist_coeffs: Tuple[Tensor, Tensor, Tensor],
    ) -> Tuple[Tensor, Tensor, Tensor, Tensor, Tensor]:
        assert (
            camera_model != "ftheta"
        ), "ftheta camera is only supported via UT, please set with_ut=True in the rasterization()"

        camera_model_type = gsplat.cuda._wrapper._make_lazy_cuda_obj(
            f"CameraModelType.{camera_model.upper()}"
        )

        aabb, (means2d, depths, proj_hardness, rgbs, normals) = _make_lazy_cuda_func(
            "projection_opaque_triangle_forward"
        )(
            (means, quats, scales, hardness, features_dc, features_sh),
            # (means, hardness),
            viewmats, Ks, width, height,
            near_plane, far_plane,
            camera_model_type, dist_coeffs
        )
        ctx.save_for_backward(means, quats, scales, hardness, features_dc, features_sh, viewmats, Ks, aabb)
        # ctx.save_for_backward(means, hardness, viewmats, Ks, aabb)
        ctx.width = width
        ctx.height = height
        ctx.camera_model_type = camera_model_type
        ctx.dist_coeffs = dist_coeffs

        return aabb, means2d, depths, proj_hardness, rgbs, normals

    @staticmethod
    def backward(ctx, v_aabb, v_means2d, v_depths, v_proj_hardness, v_rgbs, v_normals):
        means, quats, scales, hardness, features_dc, features_sh, viewmats, Ks, aabb = ctx.saved_tensors
        (v_means, v_quats, v_scales, v_hardness, v_features_dc, v_features_sh), v_viewmats = _make_lazy_cuda_func(
        # means, hardness, viewmats, Ks, aabb = ctx.saved_tensors
        # (v_means, v_hardness), v_viewmats = _make_lazy_cuda_func(
            "projection_opaque_triangle_backward"
        )(
            (means, quats, scales, hardness, features_dc, features_sh),
            # (means, hardness),
            viewmats, Ks, ctx.width, ctx.height, ctx.camera_model_type, ctx.dist_coeffs, aabb,
            tuple([x.contiguous() for x in (v_means2d, v_depths, v_proj_hardness, v_rgbs, v_normals)]),
            ctx.needs_input_grad[6],  # viewmats_requires_grad
        )
        if not ctx.needs_input_grad[6]:
            v_viewmats = None
        return (
            v_means, v_quats, v_scales, v_hardness, v_features_dc, v_features_sh, v_viewmats,
            # v_means, v_hardness, v_viewmats,
            *([None]*(len(ctx.needs_input_grad)-7))
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


intersect_splat_tile = _make_lazy_cuda_func("intersect_splat_tile")
