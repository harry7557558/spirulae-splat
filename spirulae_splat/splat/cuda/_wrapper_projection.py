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

def add_gradient_component(tensor: Tensor, attr: str, grad: Tensor):
    grad = grad.view(tensor.shape)
    if hasattr(tensor, attr):
        setattr(tensor, attr, getattr(tensor, attr) + grad)
    else:
        setattr(tensor, attr, grad)


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
    primitive: Literal["3dgs", "mip", "3dgut", "3dgut_sv", "opaque_triangle"],
    splats: tuple[Tensor],  # means, quats, scales, opacities
    viewmats: Tensor,  # [..., C, 4, 4]
    intrins: Tensor,  # [..., C, 4]
    width: int,
    height: int,
    near_plane: float = 0.0,
    far_plane: float = 1e10,
    packed: bool = False,
    sparse_grad: bool = False,
    camera_model: Literal["pinhole", "fisheye"] = "pinhole",
    dist_coeffs: Optional[Tensor] = None,
    compute_hessian_diagonal: Literal[None, "position", "all"] = None,
) -> Tuple[Tensor, Tensor, Optional[Tensor], Tuple[Tensor]]:
    batch_dims = splats[0].shape[:-2]
    N = splats[0].shape[-2]
    C = viewmats.shape[-3]
    assert viewmats.shape == batch_dims + (C, 4, 4), viewmats.shape
    assert intrins.shape == batch_dims + (C, 4), intrins.shape
    if sparse_grad:
        assert packed, "sparse_grad is only supported when packed is True"
        assert batch_dims == (), "sparse_grad does not support batch dimensions"

    if dist_coeffs is not None:
        assert dist_coeffs.shape == batch_dims + (C, 10), dist_coeffs.shape
        dist_coeffs = dist_coeffs.contiguous().to(viewmats)

    if packed:
        raise NotImplementedError("Packed not supported for fully_fused_projection without use_bvh")
    else:

        splats = [x.contiguous() for x in splats]

        in_splats = splats[:]
        if primitive in ["3dgs", "mip", "3dgut", "3dgut_sv"]:
            _FullyFusedProjection = _FullyFusedProjection3DGS
            in_splats = [primitive] + in_splats
        elif primitive in ["opaque_triangle"]:
            _FullyFusedProjection = _FullyFusedProjectionOpaqueTriangle
        elif primitive in ["voxel"]:
            _FullyFusedProjection = _FullyFusedProjectionVoxel
        proj_returns = _FullyFusedProjection.apply(
            *in_splats,
            viewmats.contiguous(),
            intrins.contiguous(),
            width,
            height,
            near_plane,
            far_plane,
            camera_model,
            dist_coeffs.contiguous() if dist_coeffs is not None else None,
            *([compute_hessian_diagonal] if primitive in ["3dgs", "mip", "3dgut"] else [])
        )

        if primitive in ["3dgs", "mip"]:
            aabb, means2d, depths, conics, opacities, rgbs = proj_returns
            return aabb, depths, (means2d, depths, conics, opacities, rgbs)
        elif primitive in ['3dgut', '3dgut_sv']:
            means, quats, scales, opacities, features_dc, features_sh = splats
            aabb, depths, scales, opacities, rgbs = proj_returns
            return aabb, depths, (means, quats, depths, scales, opacities, rgbs)
        elif primitive in ['opaque_triangle']:
            means, quats, scales, hardness, features_dc, features_sh, features_ch = splats
            aabb, depths, verts, rgbs, normals = proj_returns
            return aabb, depths, (hardness, depths, verts, rgbs, normals)
        elif primitive in ['voxel']:
            pos_sizes, densities, features_dc, features_sh = splats
            aabb, depths, rgbs = proj_returns
            return aabb, depths, (pos_sizes, densities, rgbs)
        else:
            raise NotImplementedError()


class _FullyFusedProjection3DGS(torch.autograd.Function):
    """Projects Gaussians to 2D."""

    @staticmethod
    def forward(
        ctx,
        primitive: Literal["3dgs", "mip", "3dgut", "3dgut_sv"],
        means: Tensor,  # [..., N, 3]
        quats: Tensor,  # [..., N, 4]
        scales: Tensor,  # [..., N, 3]
        opacities: Tensor,  # [..., N]
        features_dc: Tensor,  # [..., N, 3]
        features_sh: Tensor,  # [..., N, x, 3]
        viewmats: Tensor,  # [..., C, 4, 4]
        intrins: Tensor,  # [..., C, 4]
        width: int,
        height: int,
        near_plane: float,
        far_plane: float,
        camera_model: Literal["pinhole", "fisheye"],
        dist_coeffs: Optional[Tensor],
        compute_hessian_diagonal: Literal[None, "position", "all"] = None,
    ) -> Tuple[Tensor, Tensor, Tensor, Tensor, Tensor]:

        camera_model_type = gsplat.cuda._wrapper._make_lazy_cuda_obj(
            f"CameraModelType.{camera_model.upper()}"
        )

        aabb, proj_returns = _make_lazy_cuda_func(
            f"projection_{primitive}_forward"
        )(
            (means, quats, scales, opacities, features_dc, features_sh),
            viewmats, intrins, width, height,
            near_plane, far_plane,
            camera_model_type, dist_coeffs
        )
        ctx.save_for_backward(means, quats, scales, opacities, features_dc, features_sh, viewmats, intrins, aabb)
        ctx.primitive = primitive
        ctx.width = width
        ctx.height = height
        ctx.camera_model_type = camera_model_type
        ctx.dist_coeffs = dist_coeffs
        ctx.compute_hessian_diagonal = compute_hessian_diagonal

        return aabb, *proj_returns

    @staticmethod
    def backward(ctx, v_aabb, *v_proj_returns):
        if ctx.compute_hessian_diagonal:
            assert all([hasattr(v, 'hess') for v in v_proj_returns])

        means, quats, scales, opacities, features_dc, features_sh, viewmats, intrins, aabb = ctx.saved_tensors
        if ctx.compute_hessian_diagonal is not None:
            if ctx.primitive not in ["3dgs", "mip", "3dgut"]:
                raise NotImplementedError()
            (
                (v_means, v_quats, v_scales, v_opacities, v_features_dc, v_features_sh),
                v_viewmats,
                vr_proj, h_proj
            ) = _make_lazy_cuda_func(
                f"projection_{ctx.primitive}_backward_with_hessian_diagonal" if ctx.compute_hessian_diagonal == "all"
                else f"projection_{ctx.primitive}_backward_with_position_hessian_diagonal"
            )(
                (means, quats, scales, opacities, features_dc, features_sh),
                viewmats, intrins, ctx.width, ctx.height, ctx.camera_model_type, ctx.dist_coeffs, aabb,
                [x.contiguous() for x in v_proj_returns],
                [x.gradr for x in v_proj_returns],
                [x.hess for x in v_proj_returns],
                ctx.needs_input_grad[7],  # viewmats_requires_grad
            )
            assert vr_proj is not None
            assert h_proj is not None
            if ctx.compute_hessian_diagonal == "all":  # all Gaussian parameters except SH
                # print("\n"*10)
                # print('v_proj_returns')
                # for (v, vr, h) in zip([None, None, *v_proj_returns], v_proj_returns[0].gradr_all, v_proj_returns[0].hess_all):
                #     # print(x.mean().item(), y.mean().item())
                #     if (h != 0).any():
                #         print(h.shape, (vr / h)[h != 0].median())
                #     else:
                #         print(h.shape, vr.mean(), h.mean())
                # print('vr_proj and h_proj')
                # for (v, vr, h) in zip(vr_proj, vr_proj[:-1], h_proj[:-1]):
                #     # print(x.mean().item(), y.mean().item())
                #     if (h != 0).any():
                #         print(h.shape, (vr / h)[h != 0].median())
                #     else:
                #         print(h.shape, vr.mean(), h.mean())
                vr_means, vr_quats, vr_scales, vr_opacities, vr_features_dc, vr_features_sh = vr_proj
                h_means, h_quats, h_scales, h_opacities, h_features_dc, h_features_sh = h_proj
                # means
                add_gradient_component(means, 'gradr', vr_means)
                add_gradient_component(means, 'hess', h_means)
                if ctx.primitive == "3dgut":
                    add_gradient_component(means, 'gradr', v_proj_returns[0].gradr_all[0])
                    add_gradient_component(means, 'hess', v_proj_returns[0].hess_all[0])
                means.scales = scales
                means.quats = quats
                means.opacities = opacities
                # quats
                add_gradient_component(quats, 'gradr', vr_quats)
                add_gradient_component(quats, 'hess', h_quats)
                if ctx.primitive == "3dgut":
                    add_gradient_component(quats, 'gradr', v_proj_returns[0].gradr_all[1])
                    add_gradient_component(quats, 'hess', v_proj_returns[0].hess_all[1])
                quats.scales = scales
                quats.opacities = opacities
                # scales
                add_gradient_component(scales, 'gradr', vr_scales)
                add_gradient_component(scales, 'hess', h_scales)
                scales.opacities = opacities
                # opacities
                add_gradient_component(opacities, 'gradr', vr_opacities)
                add_gradient_component(opacities, 'hess', h_opacities)
                # features_dc
                add_gradient_component(features_dc, 'gradr', vr_features_dc)
                add_gradient_component(features_dc, 'hess', h_features_dc)
                features_dc.opacities = opacities
                # features_sh
                assert vr_features_sh is None
                assert h_features_sh is None
                # print('total')
                # for x in (means, quats, scales, opacities, features_dc):
                #     if (x.hess != 0).any():
                #         print(x.shape, (x.gradr / x.hess)[x.hess != 0].median().item())
                #     else:
                #         print(x.shape)
                # print("\n"*10)
            else:  # means only
                # print(h_means_from_proj.view(v_means.shape).mean(), h_means.mean())
                # print(v_proj_returns[0].gradr_all[0].mean().item())
                # print(v_means.mean().item())
                means.gradr = vr_proj.view(v_means.shape)
                means.hess = h_proj.view(v_means.shape)
                if ctx.primitive == "3dgut":
                    means.gradr += v_proj_returns[0].gradr_all[0]
                    means.hess += v_proj_returns[0].hess_all[0]
                means.scales = scales
                means.quats = quats
                means.opacities = opacities
        else:
            (v_means, v_quats, v_scales, v_opacities, v_features_dc, v_features_sh), v_viewmats = _make_lazy_cuda_func(
                f"projection_{ctx.primitive}_backward"
            )(
                (means, quats, scales, opacities, features_dc, features_sh),
                viewmats, intrins, ctx.width, ctx.height, ctx.camera_model_type, ctx.dist_coeffs, aabb,
                [x.contiguous() for x in v_proj_returns],
                ctx.needs_input_grad[7],  # viewmats_requires_grad
            )
        if not ctx.needs_input_grad[7]:
            v_viewmats = None
        return (
            None, v_means, v_quats, v_scales, v_opacities, v_features_dc, v_features_sh, v_viewmats,
            *([None]*(len(ctx.needs_input_grad)-8))
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
        intrins: Tensor,  # [..., C, 4]
        width: int,
        height: int,
        near_plane: float,
        far_plane: float,
        camera_model: Literal["pinhole", "fisheye"],
        dist_coeffs: Optional[Tensor],
    ) -> Tuple[Tensor, Tensor, Tensor, Tensor, Tensor]:

        camera_model_type = gsplat.cuda._wrapper._make_lazy_cuda_obj(
            f"CameraModelType.{camera_model.upper()}"
        )

        aabb, proj_returns = _make_lazy_cuda_func(
            "projection_opaque_triangle_forward"
        )(
            (means, quats, scales, hardness, features_dc, features_sh, features_ch),
            # (means, hardness),
            viewmats, intrins, width, height,
            near_plane, far_plane,
            camera_model_type, dist_coeffs
        )
        ctx.save_for_backward(means, quats, scales, hardness, features_dc, features_sh, features_ch, viewmats, intrins, aabb)
        # ctx.save_for_backward(means, hardness, viewmats, intrins, aabb)
        ctx.width = width
        ctx.height = height
        ctx.camera_model_type = camera_model_type
        ctx.dist_coeffs = dist_coeffs

        return aabb, *proj_returns

    @staticmethod
    def backward(ctx, v_aabb, *v_proj_returns):
        means, quats, scales, hardness, features_dc, features_sh, features_ch, viewmats, intrins, aabb = ctx.saved_tensors
        (v_means, v_quats, v_scales, v_hardness, v_features_dc, v_features_sh, v_features_ch), v_viewmats = _make_lazy_cuda_func(
        # means, hardness, viewmats, intrins, aabb = ctx.saved_tensors
        # (v_means, v_hardness), v_viewmats = _make_lazy_cuda_func(
            "projection_opaque_triangle_backward"
        )(
            (means, quats, scales, hardness, features_dc, features_sh, features_ch),
            # (means, hardness),
            viewmats, intrins, ctx.width, ctx.height, ctx.camera_model_type, ctx.dist_coeffs, aabb,
            [x.contiguous() for x in v_proj_returns],
            ctx.needs_input_grad[7],  # viewmats_requires_grad
        )
        if not ctx.needs_input_grad[7]:
            v_viewmats = None
        return (
            v_means, v_quats, v_scales, v_hardness, v_features_dc, v_features_sh, v_features_ch, v_viewmats,
            # v_means, v_hardness, v_viewmats,
            *([None]*(len(ctx.needs_input_grad)-8))
        )


class _FullyFusedProjectionVoxel(torch.autograd.Function):

    @staticmethod
    def forward(
        ctx,
        pos_sizes: Tensor,  # [..., N, 4]
        densities: Tensor,  # [..., N, 8]
        features_dc: Tensor,  # [..., N, 3]
        features_sh: Tensor,  # [..., N, x, 3]
        viewmats: Tensor,  # [..., C, 4, 4]
        intrins: Tensor,  # [..., C, 4]
        width: int,
        height: int,
        near_plane: float,
        far_plane: float,
        camera_model: Literal["pinhole", "fisheye"],
        dist_coeffs: Optional[Tensor],
    ) -> Tuple[Tensor, Tensor, Tensor, Tensor, Tensor]:

        camera_model_type = gsplat.cuda._wrapper._make_lazy_cuda_obj(
            f"CameraModelType.{camera_model.upper()}"
        )

        aabb, (depths, proj_rgbs) = _make_lazy_cuda_func(
            "projection_voxel_forward"
        )(
            (pos_sizes, densities, features_dc, features_sh),
            viewmats, intrins, width, height,
            near_plane, far_plane,
            camera_model_type, dist_coeffs
        )
        ctx.save_for_backward(pos_sizes, densities, features_dc, features_sh, viewmats, intrins, aabb)
        ctx.width = width
        ctx.height = height
        ctx.camera_model_type = camera_model_type
        ctx.dist_coeffs = dist_coeffs

        return aabb, depths, proj_rgbs

    @staticmethod
    def backward(ctx, v_aabb, v_depths, v_proj_rgbs):
        pos_sizes, densities, features_dc, features_sh, viewmats, intrins, aabb = ctx.saved_tensors
        (v_pos_sizes, v_densities, v_features_dc, v_features_sh), v_viewmats = _make_lazy_cuda_func(
            "projection_voxel_backward"
        )(
            (pos_sizes, densities, features_dc, features_sh),
            viewmats, intrins, ctx.width, ctx.height, ctx.camera_model_type, ctx.dist_coeffs, aabb,
            (v_depths, v_proj_rgbs.contiguous()),
            ctx.needs_input_grad[4],  # viewmats_requires_grad
        )
        assert v_pos_sizes is None
        assert v_densities is None
        assert v_features_dc is not None
        assert v_features_sh is not None

        if not ctx.needs_input_grad[4]:
            v_viewmats = None
        return (
            v_pos_sizes, v_densities, v_features_dc, v_features_sh, v_viewmats,
            *([None]*(len(ctx.needs_input_grad)-5))
        )



def fully_fused_projection_hetero(
    primitive: Literal["3dgs", "mip", "3dgut", "opaque_triangle"],
    splats: tuple[Tensor],  # means, quats, scales, opacities
    viewmats: Tensor,  # [..., C, 4, 4]
    intrins: Tensor,  # [..., C, 4]
    image_width: int,
    image_height: int,
    tile_width: int,
    tile_height: int,
    intersection_count_map: Tensor,  # [C+1]
    intersection_splat_id: Tensor,  # [nnz]
    near_plane: float = 0.0,
    far_plane: float = 1e10,
    sparse_grad: bool = False,
    camera_model: Literal["pinhole", "ortho", "fisheye", "ftheta"] = "pinhole",
    dist_coeffs: Optional[Tensor] = None,
) -> Tuple[Tensor, Tensor, Tensor, Tensor, Tensor]:
    batch_dims = splats[0].shape[:-2]
    C = viewmats.shape[-3]
    assert viewmats.shape == batch_dims + (C, 4, 4), viewmats.shape
    assert intrins.shape == batch_dims + (C, 4), intrins.shape
    if sparse_grad:
        assert batch_dims == (), "sparse_grad does not support batch dimensions"

    additional_args = []
    if primitive in ["3dgs", "mip", "3dgut"]:
        _FullyFusedProjection = _FullyFusedProjection3DGSHetero
        in_splats = [primitive] + in_splats
    elif primitive in ["opaque_triangle"]:
        _FullyFusedProjection = _FullyFusedProjectionOpaqueTriangleHetero
    proj_returns = _FullyFusedProjection.apply(
        *in_splats,
        viewmats.contiguous(), intrins.contiguous(),
        image_width, image_height, tile_width, tile_height,
        near_plane, far_plane,
        camera_model, dist_coeffs.contiguous() if dist_coeffs is not None else None,
        intersection_count_map.contiguous(), intersection_splat_id.contiguous(),
        *additional_args
    )

    if primitive in ["3dgs", "mip"]:
        means2d, depths, conics, opacities, rgbs, absgrad = proj_returns[3:]
        return *proj_returns[:3], depths, (means2d, depths, conics, opacities, rgbs)
    elif primitive in ['3dgut']:
        means, quats, depths, scales, opacities, rgbs = proj_returns[3:]
        return *proj_returns[:3], depths, (means, quats, depths, scales, opacities, rgbs)
    elif primitive in ['opaque_triangle']:
        hardness, depths, verts, rgbs, normals = proj_returns[3:]
        return *proj_returns[:3], depths, (hardness, depths, verts, rgbs, normals)
    elif primitive in ['voxel']:
        pos_sizes, densities, features_dc, features_sh = splats
        depths, rgbs = proj_returns[3:]
        return *proj_returns[:3], depths, (pos_sizes, densities, rgbs)
    else:
        raise NotImplementedError()


class _FullyFusedProjection3DGSHetero(torch.autograd.Function):
    """Projects Gaussians to 2D. Return packed tensors."""

    @staticmethod
    def forward(
        ctx,
        primitive: Literal["3dgs", "mip", "3dgut"],
        means: Tensor,  # [..., N, 3]
        quats: Tensor,  # [..., N, 4]
        scales: Tensor,  # [..., N, 3]
        opacities: Tensor,  # [..., N]
        features_dc: Tensor,  # [..., N, 3]
        features_sh: Tensor,  # [..., N, x, 3]
        viewmats: Tensor,  # [..., C, 4, 4]
        intrins: Tensor,  # [..., C, 4]
        image_width: int,
        image_height: int,
        tile_width: int,
        tile_height: int,
        near_plane: float,
        far_plane: float,
        camera_model: Literal["pinhole", "fisheye"],
        dist_coeffs: Optional[Tensor],
        intersection_count_map: Tensor,  # [C+1]
        intersection_splat_id: Tensor,  # [nnz]
    ):
        camera_model_type = gsplat.cuda._wrapper._make_lazy_cuda_obj(
            f"CameraModelType.{camera_model.upper()}"
        )

        (
            camera_ids, gaussian_ids, aabb, proj_returns
        ) = _make_lazy_cuda_func(f"projection_{primitive}_hetero_forward")(
            (means, quats, scales, opacities, features_dc, features_sh),
            viewmats, intrins, image_width, image_height, tile_width, tile_height,
            near_plane, far_plane, camera_model_type, dist_coeffs,
            intersection_count_map, intersection_splat_id
        )
        ctx.save_for_backward(
            means, quats, scales, opacities, features_dc, features_sh,
            viewmats, intrins, dist_coeffs, aabb,
            camera_ids, gaussian_ids
        )
        ctx.primitive = primitive
        ctx.camera = (image_width, image_height, tile_width, tile_height, camera_model_type)

        return camera_ids, gaussian_ids, aabb, *proj_returns

    @staticmethod
    def backward(
        ctx,
        v_camera_ids,
        v_gaussian_ids,
        v_aabb,
        *v_proj_returns,
    ):
        (
            means, quats, scales, opacities, features_dc, features_sh,
            viewmats, intrins, dist_coeffs, aabb,
            camera_ids, gaussian_ids
        ) = ctx.saved_tensors
        (image_width, image_height, tile_width, tile_height, camera_model_type) = ctx.camera

        (v_means, v_quats, v_scales, v_opacities, v_features_dc, v_features_sh), v_viewmats = _make_lazy_cuda_func(
            f"projection_{ctx.primitive}_hetero_backward"
        )(
            (means, quats, scales, opacities, features_dc, features_sh),
            viewmats, intrins, image_width, image_height, tile_width, tile_height, camera_model_type, dist_coeffs,
            camera_ids, gaussian_ids, aabb,
            tuple([(x.contiguous() if isinstance(x, torch.Tensor) else None) for x in v_proj_returns]),
            ctx.needs_input_grad[7], False,
        )

        return (
            None,
            v_means, v_quats, v_scales, v_opacities, v_features_dc, v_features_sh,
            v_viewmats,
            *([None]*(len(ctx.needs_input_grad)-8))
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
        intrins: Tensor,  # [..., C, 4]
        image_width: int,
        image_height: int,
        tile_width: int,
        tile_height: int,
        near_plane: float,
        far_plane: float,
        camera_model: Literal["pinhole", "fisheye"],
        dist_coeffs: Optional[Tensor],
        intersection_count_map: Tensor,  # [C+1]
        intersection_splat_id: Tensor,  # [nnz]
    ) -> Tuple[Tensor, Tensor, Tensor, Tensor, Tensor]:

        camera_model_type = gsplat.cuda._wrapper._make_lazy_cuda_obj(
            f"CameraModelType.{camera_model.upper()}"
        )

        (
            camera_ids, gaussian_ids, aabb, proj_returns
        ) = _make_lazy_cuda_func(
            "projection_opaque_triangle_hetero_forward"
        )(
            (means, quats, scales, hardness, features_dc, features_sh, features_ch),
            viewmats, intrins, image_width, image_height, tile_width, tile_height,
            near_plane, far_plane, camera_model_type, dist_coeffs,
            intersection_count_map, intersection_splat_id
        )
        ctx.save_for_backward(
            means, quats, scales, hardness, features_dc, features_sh, features_ch,
            viewmats, intrins, dist_coeffs, aabb, 
            camera_ids, gaussian_ids
        )
        ctx.camera = (image_width, image_height, tile_width, tile_height, camera_model_type)

        return camera_ids, gaussian_ids, aabb, *proj_returns

    @staticmethod
    def backward(ctx, v_camera_ids, v_gaussian_ids, v_aabb, *v_proj_returns):
        (
            means, quats, scales, hardness, features_dc, features_sh, features_ch,
            viewmats, intrins, dist_coeffs, aabb,
            camera_ids, gaussian_ids
        ) = ctx.saved_tensors
        (image_width, image_height, tile_width, tile_height, camera_model_type) = ctx.camera
        sparse_grad = False

        (v_means, v_quats, v_scales, v_hardness, v_features_dc, v_features_sh, v_features_ch), v_viewmats = _make_lazy_cuda_func(
            "projection_opaque_triangle_hetero_backward"
        )(
            (means, quats, scales, hardness, features_dc, features_sh, features_ch),
            viewmats, intrins, image_width, image_height, tile_width, tile_height, camera_model_type, dist_coeffs,
            camera_ids, gaussian_ids, aabb,
            tuple([(x.contiguous() if isinstance(x, torch.Tensor) else None) for x in v_proj_returns]),
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
    primitive: Literal["3dgs", "mip", "3dgut", "opaque_triangle", "voxel"],
    splats: Tuple,
    width: int,
    height: int,
    viewmats: Tensor,
    intrins: Tensor,
    camera_model: Literal["pinhole", "fisheye"],
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
            '3dgut': '3dgs',
            'opaque_triangle': 'opaque_triangle',
            'voxel': 'voxel',
        }[primitive]
    )(
        splats,
        width,
        height,
        viewmats.contiguous(),
        intrins.contiguous(),
        camera_model,
        dist_coeffs.contiguous() if dist_coeffs is not None else None,
        relative_scale
    )
