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
    eval3d: bool,
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
    dist_coeffs: Optional[Tensor] = None,
) -> Tuple[Tensor, Tensor, Optional[Tensor], Tuple[Tensor]]:
    # if primitive in ["3dgs", "mip"]:
    if primitive in ["3dgs", "mip", "opaque_triangle"]:
        if primitive in ["3dgs", "mip"]:
            means, quats, scales, opacities, features_dc, features_sh = splats
        else:
            means, quats, scales, opacities, features_dc, features_sh, features_ch = splats
        batch_dims = means.shape[:-2]
        N = means.shape[-2]
        assert means.shape == batch_dims + (N, 3), means.shape
        assert quats.shape == batch_dims + (N, 4), quats.shape
        assert scales.shape == batch_dims + (N, 3), scales.shape
        if primitive in ["3dgs", "mip"]:
            assert opacities.shape == batch_dims + (N,), opacities.shape
        else:
            assert opacities.shape == batch_dims + (N, 2), opacities.shape
            assert features_ch.shape == batch_dims + (N, 2, 3), features_ch.shape
        assert features_dc.shape == batch_dims + (N, 3), features_dc.shape
    C = viewmats.shape[-3]
    assert viewmats.shape == batch_dims + (C, 4, 4), viewmats.shape
    assert Ks.shape == batch_dims + (C, 3, 3), Ks.shape
    if sparse_grad:
        assert packed, "sparse_grad is only supported when packed is True"
        assert batch_dims == (), "sparse_grad does not support batch dimensions"

    if dist_coeffs is not None:
        assert dist_coeffs.shape == batch_dims + (C, 10), dist_coeffs.shape
        dist_coeffs = dist_coeffs.contiguous().to(viewmats)

    assert (
        camera_model != "ftheta"
    ), "ftheta camera is only supported via UT, please set with_ut=True in the rasterization()"

    if packed:
        raise NotImplementedError("Packed not supported for fully_fused_projection without use_bvh")
    else:
        additional_args = []
        if not eval3d:
            if primitive in ["3dgs", "mip"]:
                _FullyFusedProjection = _FullyFusedProjection3DGS
                additional_args = [primitive == "mip"]
            elif primitive in ["opaque_triangle"]:
                _FullyFusedProjection = _FullyFusedProjectionOpaqueTriangle
        else:
            if primitive in ["3dgs", "mip"]:
                _FullyFusedProjection = _FullyFusedProjection3DGSEval3D
            elif primitive in ["opaque_triangle"]:
                _FullyFusedProjection = _FullyFusedProjectionOpaqueTriangleEval3D
        in_splats = [x.contiguous() for x in splats]
        proj_return = _FullyFusedProjection.apply(
            *in_splats,
            viewmats.contiguous(),
            Ks.contiguous(),
            width,
            height,
            near_plane,
            far_plane,
            camera_model,
            dist_coeffs.contiguous() if dist_coeffs is not None else None,
            *additional_args
        )
        if not eval3d:
            depths = proj_return[2]
            out_splats = tuple(proj_return[1:])
            if primitive == "opaque_triangle":
                depths = depths.mean(-1)
        else:
            depths = proj_return[1]
            out_splats = (*in_splats, *proj_return[1:])
        return (
            proj_return[0],  # aabb
            depths,  # depth
            out_splats
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
        dist_coeffs: Optional[Tensor],
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
        features_ch: Tensor,  # [..., N, 2, 3]
        viewmats: Tensor,  # [..., C, 4, 4]
        Ks: Tensor,  # [..., C, 3, 3]
        width: int,
        height: int,
        near_plane: float,
        far_plane: float,
        camera_model: Literal["pinhole", "ortho", "fisheye", "ftheta"],
        dist_coeffs: Optional[Tensor],
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
            (means, quats, scales, hardness, features_dc, features_sh, features_ch),
            # (means, hardness),
            viewmats, Ks, width, height,
            near_plane, far_plane,
            camera_model_type, dist_coeffs
        )
        ctx.save_for_backward(means, quats, scales, hardness, features_dc, features_sh, features_ch, viewmats, Ks, aabb)
        # ctx.save_for_backward(means, hardness, viewmats, Ks, aabb)
        ctx.width = width
        ctx.height = height
        ctx.camera_model_type = camera_model_type
        ctx.dist_coeffs = dist_coeffs

        return aabb, means2d, depths, proj_hardness, rgbs, normals

    @staticmethod
    def backward(ctx, v_aabb, v_means2d, v_depths, v_proj_hardness, v_rgbs, v_normals):
        means, quats, scales, hardness, features_dc, features_sh, features_ch, viewmats, Ks, aabb = ctx.saved_tensors
        (v_means, v_quats, v_scales, v_hardness, v_features_dc, v_features_sh, v_features_ch), v_viewmats = _make_lazy_cuda_func(
        # means, hardness, viewmats, Ks, aabb = ctx.saved_tensors
        # (v_means, v_hardness), v_viewmats = _make_lazy_cuda_func(
            "projection_opaque_triangle_backward"
        )(
            (means, quats, scales, hardness, features_dc, features_sh, features_ch),
            # (means, hardness),
            viewmats, Ks, ctx.width, ctx.height, ctx.camera_model_type, ctx.dist_coeffs, aabb,
            tuple([x.contiguous() for x in (v_means2d, v_depths, v_proj_hardness, v_rgbs, v_normals)]),
            ctx.needs_input_grad[7],  # viewmats_requires_grad
        )
        if not ctx.needs_input_grad[7]:
            v_viewmats = None
        return (
            v_means, v_quats, v_scales, v_hardness, v_features_dc, v_features_sh, v_features_ch, v_viewmats,
            # v_means, v_hardness, v_viewmats,
            *([None]*(len(ctx.needs_input_grad)-8))
        )

class _FullyFusedProjection3DGSEval3D(torch.autograd.Function):

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
        dist_coeffs: Optional[Tensor],
    ) -> Tuple[Tensor, Tensor, Tensor, Tensor, Tensor]:
        assert (
            camera_model != "ftheta"
        ), "ftheta camera is only supported via UT, please set with_ut=True in the rasterization()"

        camera_model_type = gsplat.cuda._wrapper._make_lazy_cuda_obj(
            f"CameraModelType.{camera_model.upper()}"
        )

        aabb, (depths, proj_scales, proj_opacities, proj_rgbs) = _make_lazy_cuda_func(
            "projection_ewa_3dgs_eval3d_forward"
        )(
            (means, quats, scales, opacities, features_dc, features_sh),
            viewmats, Ks, width, height,
            near_plane, far_plane,
            camera_model_type, dist_coeffs
        )
        ctx.save_for_backward(means, quats, scales, opacities, features_dc, features_sh, viewmats, Ks, aabb)
        ctx.width = width
        ctx.height = height
        ctx.camera_model_type = camera_model_type
        ctx.dist_coeffs = dist_coeffs

        return aabb, depths, proj_scales, proj_opacities, proj_rgbs

    @staticmethod
    def backward(ctx, v_aabb, v_depths, v_proj_scales, v_proj_opacities, v_proj_rgbs):
        means, quats, scales, opacities, features_dc, features_sh, viewmats, Ks, aabb = ctx.saved_tensors
        (v_means, v_quats, v_scales, v_opacities, v_features_dc, v_features_sh), v_viewmats = _make_lazy_cuda_func(
            "projection_ewa_3dgs_eval3d_backward"
        )(
            (means, quats, scales, opacities, features_dc, features_sh),
            viewmats, Ks, ctx.width, ctx.height, ctx.camera_model_type, ctx.dist_coeffs, aabb,
            tuple([x.contiguous() for x in (v_depths, v_proj_scales, v_proj_opacities, v_proj_rgbs)]),
            ctx.needs_input_grad[6],  # viewmats_requires_grad
        )
        if not ctx.needs_input_grad[6]:
            v_viewmats = None
        else:
            raise NotImplementedError("viewmat grad")
        return (
            v_means, v_quats, v_scales, v_opacities, v_features_dc, v_features_sh, v_viewmats,
            *([None]*(len(ctx.needs_input_grad)-7))
        )

class _FullyFusedProjectionOpaqueTriangleEval3D(torch.autograd.Function):
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
        features_ch: Tensor,  # [..., N, 2, 3]
        viewmats: Tensor,  # [..., C, 4, 4]
        Ks: Tensor,  # [..., C, 3, 3]
        width: int,
        height: int,
        near_plane: float,
        far_plane: float,
        camera_model: Literal["pinhole", "ortho", "fisheye", "ftheta"],
        dist_coeffs: Optional[Tensor],
    ) -> Tuple[Tensor, Tensor, Tensor, Tensor, Tensor]:
        assert (
            camera_model != "ftheta"
        ), "ftheta camera is only supported via UT, please set with_ut=True in the rasterization()"

        camera_model_type = gsplat.cuda._wrapper._make_lazy_cuda_obj(
            f"CameraModelType.{camera_model.upper()}"
        )

        aabb, (depths, verts, rgbs, normals) = _make_lazy_cuda_func(
            "projection_opaque_triangle_eval3d_forward"
        )(
            (means, quats, scales, hardness, features_dc, features_sh, features_ch),
            viewmats, Ks, width, height,
            near_plane, far_plane,
            camera_model_type, dist_coeffs
        )
        ctx.save_for_backward(means, quats, scales, hardness, features_dc, features_sh, features_ch, viewmats, Ks, aabb)
        ctx.width = width
        ctx.height = height
        ctx.camera_model_type = camera_model_type
        ctx.dist_coeffs = dist_coeffs

        return aabb, depths, verts, rgbs, normals

    @staticmethod
    def backward(ctx, v_aabb, v_depths, v_verts, v_rgbs, v_normals):
        means, quats, scales, hardness, features_dc, features_sh, features_ch, viewmats, Ks, aabb = ctx.saved_tensors
        (v_means, v_quats, v_scales, v_hardness, v_features_dc, v_features_sh, v_features_ch), v_viewmats = _make_lazy_cuda_func(
            "projection_opaque_triangle_eval3d_backward"
        )(
            (means, quats, scales, hardness, features_dc, features_sh, features_ch),
            viewmats, Ks, ctx.width, ctx.height, ctx.camera_model_type, ctx.dist_coeffs, aabb,
            tuple([x.contiguous() for x in (v_depths, v_verts, v_rgbs, v_normals)]),
            ctx.needs_input_grad[7],  # viewmats_requires_grad
        )
        if not ctx.needs_input_grad[7]:
            v_viewmats = None
        return (
            v_means, v_quats, v_scales, v_hardness, v_features_dc, v_features_sh, v_features_ch, v_viewmats,
            *([None]*(len(ctx.needs_input_grad)-8))
        )



def fully_fused_projection_hetero(
    primitive: Literal["3dgs", "mip", "opaque_triangle"],
    splats: tuple[Tensor],  # means, quats, scales, opacities
    viewmats: Tensor,  # [..., C, 4, 4]
    Ks: Tensor,  # [..., C, 3, 3]
    image_width: int,
    image_height: int,
    tile_width: int,
    tile_height: int,
    intersection_count_map: Tensor,  # [C+1]
    intersection_splat_id: Tensor,  # [nnz]
    near_plane: float = 0.01,
    far_plane: float = 1e10,
    sparse_grad: bool = False,
    camera_model: Literal["pinhole", "ortho", "fisheye", "ftheta"] = "pinhole",
    dist_coeffs: Optional[Tensor] = None,
) -> Tuple[Tensor, Tensor, Tensor, Tensor, Tensor]:
    # if primitive in ["3dgs", "mip"]:
    if primitive in ["3dgs", "mip", "opaque_triangle"]:
        if primitive in ["3dgs", "mip"]:
            means, quats, scales, opacities, features_dc, features_sh = splats
        else:
            means, quats, scales, opacities, features_dc, features_sh, features_ch = splats
        batch_dims = means.shape[:-2]
        batch_dims = means.shape[:-2]
        N = means.shape[-2]
        C = viewmats.shape[-3]
        assert means.shape == batch_dims + (N, 3), means.shape
        assert quats.shape == batch_dims + (N, 4), quats.shape
        assert scales.shape == batch_dims + (N, 3), scales.shape
        if primitive in ["3dgs", "mip"]:
            assert opacities.shape == batch_dims + (N,), opacities.shape
        else:
            assert opacities.shape == batch_dims + (N, 2), opacities.shape
            assert features_ch.shape == batch_dims + (N, 2, 3), features_ch.shape
        assert features_dc.shape == batch_dims + (N, 3), features_dc.shape
    assert viewmats.shape == batch_dims + (C, 4, 4), viewmats.shape
    assert Ks.shape == batch_dims + (C, 3, 3), Ks.shape
    if sparse_grad:
        assert batch_dims == (), "sparse_grad does not support batch dimensions"

    additional_args = []
    if primitive in ["3dgs", "mip"]:
        _FullyFusedProjection = _FullyFusedProjection3DGSHetero
    elif primitive in ["opaque_triangle"]:
        _FullyFusedProjection = _FullyFusedProjectionOpaqueTriangleHetero
    in_splats = [x.contiguous() for x in splats]
    proj_return = _FullyFusedProjection.apply(
        *in_splats,
        viewmats.contiguous(), Ks.contiguous(),
        image_width, image_height, tile_width, tile_height,
        near_plane, far_plane,
        camera_model, dist_coeffs.contiguous() if dist_coeffs is not None else None,
        intersection_count_map.contiguous(), intersection_splat_id.contiguous(),
        *additional_args
    )
    splat_ids = proj_return[1]
    depths = proj_return[3]
    # out_splats = (*in_splats, *proj_return[3:])
    out_splats = (*[x[splat_ids] for x in in_splats], *proj_return[3:])
    return (
        proj_return[0],  # camera_ids
        splat_ids,  # splat_ids
        proj_return[2],  # aabb
        depths,  # depth
        out_splats
    )


class _FullyFusedProjection3DGSHetero(torch.autograd.Function):
    """Projects Gaussians to 2D. Return packed tensors."""

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
        image_width: int,
        image_height: int,
        tile_width: int,
        tile_height: int,
        near_plane: float,
        far_plane: float,
        camera_model: Literal["pinhole", "ortho", "fisheye", "ftheta"],
        dist_coeffs: Optional[Tensor],
        intersection_count_map: Tensor,  # [C+1]
        intersection_splat_id: Tensor,  # [nnz]
    ):
        assert (
            camera_model != "ftheta"
        ), "ftheta camera is only supported via UT, please set with_ut=True in the rasterization()"

        camera_model_type = gsplat.cuda._wrapper._make_lazy_cuda_obj(
            f"CameraModelType.{camera_model.upper()}"
        )

        (
            camera_ids, gaussian_ids, aabb,
            (depths, proj_scales, proj_opacities, proj_rgbs)
        ) = _make_lazy_cuda_func("projection_ewa_3dgs_hetero_forward")(
            (means, quats, scales, opacities, features_dc, features_sh),
            viewmats, Ks, image_width, image_height, tile_width, tile_height,
            near_plane, far_plane, camera_model_type, dist_coeffs,
            intersection_count_map, intersection_splat_id
        )
        ctx.save_for_backward(
            means, quats, scales, opacities, features_dc, features_sh,
            viewmats, Ks, dist_coeffs, aabb,
            camera_ids, gaussian_ids
        )
        ctx.camera = (image_width, image_height, tile_width, tile_height, camera_model_type)

        return camera_ids, gaussian_ids, aabb, depths, proj_scales, proj_opacities, proj_rgbs

    @staticmethod
    def backward(
        ctx,
        v_camera_ids,
        v_gaussian_ids,
        v_aabb,
        v_depths,
        v_proj_scales,
        v_proj_opacities,
        v_proj_rgbs,
    ):
        (
            means, quats, scales, opacities, features_dc, features_sh,
            viewmats, Ks, dist_coeffs, aabb,
            camera_ids, gaussian_ids
        ) = ctx.saved_tensors
        (image_width, image_height, tile_width, tile_height, camera_model_type) = ctx.camera
        sparse_grad = False

        (v_means, v_quats, v_scales, v_opacities, v_features_dc, v_features_sh), v_viewmats = _make_lazy_cuda_func(
            "projection_ewa_3dgs_hetero_backward"
        )(
            (means, quats, scales, opacities, features_dc, features_sh),
            viewmats, Ks, image_width, image_height, tile_width, tile_height, camera_model_type, dist_coeffs,
            camera_ids, gaussian_ids, aabb,
            tuple([x.contiguous() for x in (v_depths, v_proj_scales, v_proj_opacities, v_proj_rgbs)]),
            ctx.needs_input_grad[6], sparse_grad,
        )

        if sparse_grad:
            batch_dims = means.shape[:-2]
            B = math.prod(batch_dims)
            N = means.shape[-2]
        if ctx.needs_input_grad[0] and sparse_grad:
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

        return (
            v_means, v_quats, v_scales, v_opacities, v_features_dc, v_features_sh,
            v_viewmats,
            *([None]*(len(ctx.needs_input_grad)-7))
        )


class _FullyFusedProjectionOpaqueTriangleHetero(torch.autograd.Function):
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
        features_ch: Tensor,  # [..., N, 2, 3]
        viewmats: Tensor,  # [..., C, 4, 4]
        Ks: Tensor,  # [..., C, 3, 3]
        image_width: int,
        image_height: int,
        tile_width: int,
        tile_height: int,
        near_plane: float,
        far_plane: float,
        camera_model: Literal["pinhole", "ortho", "fisheye", "ftheta"],
        dist_coeffs: Optional[Tensor],
        intersection_count_map: Tensor,  # [C+1]
        intersection_splat_id: Tensor,  # [nnz]
    ) -> Tuple[Tensor, Tensor, Tensor, Tensor, Tensor]:
        assert (
            camera_model != "ftheta"
        ), "ftheta camera is only supported via UT, please set with_ut=True in the rasterization()"

        camera_model_type = gsplat.cuda._wrapper._make_lazy_cuda_obj(
            f"CameraModelType.{camera_model.upper()}"
        )

        (
            camera_ids, gaussian_ids, aabb,
            (depths, verts, rgbs, normals)
        ) = _make_lazy_cuda_func(
            "projection_opaque_triangle_hetero_forward"
        )(
            (means, quats, scales, hardness, features_dc, features_sh, features_ch),
            viewmats, Ks, image_width, image_height, tile_width, tile_height,
            near_plane, far_plane, camera_model_type, dist_coeffs,
            intersection_count_map, intersection_splat_id
        )
        ctx.save_for_backward(
            means, quats, scales, hardness, features_dc, features_sh, features_ch,
            viewmats, Ks, dist_coeffs, aabb, 
            camera_ids, gaussian_ids
        )
        ctx.camera = (image_width, image_height, tile_width, tile_height, camera_model_type)

        return camera_ids, gaussian_ids, aabb, depths, verts, rgbs, normals

    @staticmethod
    def backward(ctx, v_camera_ids, v_gaussian_ids, v_aabb, v_depths, v_verts, v_rgbs, v_normals):
        (
            means, quats, scales, hardness, features_dc, features_sh, features_ch,
            viewmats, Ks, dist_coeffs, aabb, 
            camera_ids, gaussian_ids
        ) = ctx.saved_tensors
        (image_width, image_height, tile_width, tile_height, camera_model_type) = ctx.camera
        sparse_grad = False

        (v_means, v_quats, v_scales, v_hardness, v_features_dc, v_features_sh, v_features_ch), v_viewmats = _make_lazy_cuda_func(
            "projection_opaque_triangle_hetero_backward"
        )(
            (means, quats, scales, hardness, features_dc, features_sh, features_ch),
            viewmats, Ks, image_width, image_height, tile_width, tile_height, camera_model_type, dist_coeffs,
            camera_ids, gaussian_ids, aabb,
            tuple([x.contiguous() for x in (v_depths, v_verts, v_rgbs, v_normals)]),
            ctx.needs_input_grad[7], sparse_grad,
        )
        if not ctx.needs_input_grad[7]:
            v_viewmats = None
        if sparse_grad:
            raise NotImplementedError()
        return (
            v_means, v_quats, v_scales, v_hardness, v_features_dc, v_features_sh, v_features_ch, v_viewmats,
            *([None]*(len(ctx.needs_input_grad)-8))
        )


@torch.no_grad()
def intersect_splat_tile(
    primitive: Literal["3dgs", "mip", "opaque_triangle"],
    splats: Tuple,
    width: int,
    height: int,
    viewmats: Tensor,
    Ks: Tensor,
    camera_model: Literal["pinhole", "ortho", "fisheye", "ftheta"] = "pinhole",
    dist_coeffs: Optional[Tensor] = None,
    relative_scale: float = 1.0
) -> Tuple[Tensor, Tensor]:
    camera_model = gsplat.cuda._wrapper._make_lazy_cuda_obj(
        f"CameraModelType.{camera_model.upper()}"
    )
    return _make_lazy_cuda_func(
        "intersect_splat_tile_" + {
            '3dgs': '3dgs',
            'mip': '3dgs',
            'opaque_triangle': 'opaque_triangle',
        }[primitive]
    )(
        splats,
        width,
        height,
        viewmats.contiguous(),
        Ks.contiguous(),
        camera_model,
        dist_coeffs.contiguous() if dist_coeffs is not None else None,
        relative_scale
    )
