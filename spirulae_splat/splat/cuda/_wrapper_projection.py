# Modified from https://github.com/nerfstudio-project/gsplat/blob/2323de5905d5e90e035f792fe65bad0fedd413e7/gsplat/cuda/_wrapper.py

import math
import warnings
from dataclasses import dataclass
from enum import Enum
from typing import Any, Callable, Optional, Tuple, List, Union

import torch
from torch import Tensor
from typing_extensions import Literal

# import gsplat.cuda._wrapper

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

def add_gradient_component(obj: Union[Tensor, dict], key: str, grad: Tensor):
    if isinstance(obj, torch.Tensor):
        assert False, "add_gradient_component called with a tensor instead of backward_info"
    if key in obj:
        obj[key] += grad
    else:
        obj[key] = grad


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


def scatter_add(ref_tensor: Tensor, tensor: Tensor, indices: Tensor):
    v_tensor = torch.zeros_like(ref_tensor)
    _make_lazy_cuda_func("inplace_scatter_add")(
        indices, tensor, v_tensor
    )
    return v_tensor


class _Index(torch.autograd.Function):
    @staticmethod
    def forward(
        ctx, tensor: Tensor, indices: Tensor
    ):
        ctx.save_for_backward(tensor, indices)
        # return tensor[indices]
        result = torch.empty(indices.shape[0], *tensor.shape[1:], device=tensor.device, dtype=tensor.dtype)
        _make_lazy_cuda_func("inplace_index")(
            indices, tensor, result
        )
        return result

    @staticmethod
    def backward(ctx, v_output: Tensor):
        tensor, indices = ctx.saved_tensors
        v_tensor = scatter_add(tensor, v_output, indices)
        return v_tensor, None


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
    backward_info: Optional[dict] = None,
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

    if packed and False:
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
            packed,
            *([compute_hessian_diagonal, backward_info] if primitive in ["3dgs", "mip", "3dgut"] else [])
        )

        if primitive in ["3dgs", "mip"]:
            aabb, means2d, depths, conics, opacities, rgbs = proj_returns
            return aabb, depths, (means2d, depths, conics, opacities, rgbs)
        elif primitive in ['3dgut', '3dgut_sv']:
            means, quats, scales, opacities, features_dc, features_sh = splats
            aabb, depths, scales, opacities, rgbs = proj_returns
            # if packed:
            #     means = _Index.apply(means, aabb.gaussian_ids)
            #     quats = _Index.apply(quats, aabb.gaussian_ids)
            backward_info['gaussian_ids'] = aabb.gaussian_ids
            return aabb, depths, (means, quats, depths, scales, opacities, rgbs)
        elif primitive in ['opaque_triangle']:
            means, quats, scales, hardness, features_dc, features_sh, features_ch = splats
            aabb, depths, verts, rgbs, normals = proj_returns
            if packed:
                hardness = _Index.apply(hardness, aabb.gaussian_ids)
            return aabb, depths, (hardness, depths, verts, rgbs, normals)
        elif primitive in ['voxel']:
            pos_sizes, densities, features_dc, features_sh = splats
            aabb, depths, rgbs = proj_returns
            if packed:
                pos_sizes = _Index.apply(pos_sizes, aabb.gaussian_ids)
                densities = _Index.apply(densities, aabb.gaussian_ids)
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
        packed: bool,
        compute_hessian_diagonal: Literal[None, "position", "all"] = None,
        backward_info: Optional[dict] = None,
    ) -> Tuple[Tensor, Tensor, Tensor, Tensor, Tensor]:

        camera_model_type = camera_model.upper()

        cuda_returns = _make_lazy_cuda_func(
            f"projection_{primitive}_packed_forward" if packed else
            f"projection_{primitive}_forward"
        )(
            (means, quats, scales, opacities, features_dc, features_sh),
            viewmats, intrins, width, height,
            near_plane, far_plane,
            camera_model_type, dist_coeffs
        )
        if packed:
            camera_ids, gaussian_ids, aabb, proj_returns = cuda_returns
            aabb.camera_ids = camera_ids
            aabb.gaussian_ids = gaussian_ids
        else:
            aabb, proj_returns = cuda_returns
            camera_ids, gaussian_ids = None, None

        ctx.save_for_backward(
            means, quats, scales, opacities, features_dc, features_sh,
            viewmats, intrins, camera_ids, gaussian_ids, aabb
        )
        ctx.primitive = primitive
        ctx.width = width
        ctx.height = height
        ctx.camera_model_type = camera_model_type
        ctx.dist_coeffs = dist_coeffs
        ctx.compute_hessian_diagonal = compute_hessian_diagonal
        ctx.backward_info = backward_info

        return aabb, *proj_returns

    @staticmethod
    def backward(ctx, v_aabb, *v_proj_returns):

        (
            means, quats, scales, opacities, features_dc, features_sh,
            viewmats, intrins, camera_ids, gaussian_ids, aabb
        ) = ctx.saved_tensors

        grad_accum_mask = [
            x.is_leaf and hasattr(x, 'grad') and x.grad is not None
            for x in (means, quats, scales, opacities, features_dc, features_sh)
        ]
        v_proj = (*[
            x.grad if is_accum else torch.zeros_like(x)
            for is_accum, x in zip(grad_accum_mask, (means, quats, scales, opacities, features_dc, features_sh))
        ],)
        if ctx.compute_hessian_diagonal == "all":
            # TODO: in place gradient accumulation
            assert ctx.backward_info is not None
            vr_proj = (*[
                torch.zeros_like(x)
                for x in (means, quats, scales, opacities, features_dc)
            ], None)
            h_proj = (*[
                torch.zeros_like(x)
                for x in (means, quats, scales, opacities, features_dc)
            ], None)
        elif ctx.compute_hessian_diagonal == "position":
            vr_proj = torch.zeros_like(means)
            h_proj = torch.zeros_like(means)
        v_viewmats = torch.zeros_like(viewmats) if ctx.needs_input_grad[7] else None

        if ctx.compute_hessian_diagonal is not None:
            if ctx.primitive not in ["3dgs", "mip", "3dgut"]:
                raise NotImplementedError()
            raster_return_keys = 'means2d depths conics proj_opacities colors'.split() \
                if ctx.primitive in ['3dgs', 'mip'] else 'depths proj_scales proj_opacities colors'.split()
            assert len(v_proj_returns) == len(raster_return_keys), ([x.shape for x in v_proj_returns], raster_return_keys)
            proj_return_keys = 'means quats scales opacities features_dc'.split() \
                if ctx.compute_hessian_diagonal == "all" else ['means']
            _make_lazy_cuda_func(
                f"projection_{ctx.primitive}_backward_with_hessian_diagonal" if ctx.compute_hessian_diagonal == "all"
                else f"projection_{ctx.primitive}_backward_with_position_hessian_diagonal"
            )(
                (means, quats, scales, opacities, features_dc, features_sh),
                viewmats, intrins, ctx.width, ctx.height, ctx.camera_model_type, ctx.dist_coeffs,
                camera_ids, gaussian_ids, aabb,
                [x.contiguous() for x in v_proj_returns],
                [ctx.backward_info[key+'.gradr'] for key in raster_return_keys],
                [ctx.backward_info[key+'.hess'] for key in raster_return_keys],
                v_proj, v_viewmats, vr_proj, h_proj,
            )
            for key in raster_return_keys:
                if key not in proj_return_keys:
                    del ctx.backward_info[key+'.gradr']
                    del ctx.backward_info[key+'.hess']
            for key in proj_return_keys:
                if 'proj_'+key+'.gradr' not in ctx.backward_info and 'proj_'+key+'.hess' not in ctx.backward_info:
                    continue
                assert 'proj_'+key+'.gradr' in ctx.backward_info
                assert 'proj_'+key+'.hess' in ctx.backward_info
                if gaussian_ids is not None:  # packed
                    assert ctx.primitive in ["3dgut"]
                    # ref_tensor = {"means": means, "quats": quats}[key]
                    # ctx.backward_info[key+'.gradr'] = scatter_add(ref_tensor, ctx.backward_info['proj_'+key+'.gradr'], gaussian_ids)
                    # ctx.backward_info[key+'.hess'] = scatter_add(ref_tensor, ctx.backward_info['proj_'+key+'.hess'], gaussian_ids)
                    ctx.backward_info[key+'.gradr'] = ctx.backward_info['proj_'+key+'.gradr']
                    ctx.backward_info[key+'.hess'] = ctx.backward_info['proj_'+key+'.hess']
                    del ctx.backward_info['proj_'+key+'.gradr']
                    del ctx.backward_info['proj_'+key+'.hess']
            if ctx.compute_hessian_diagonal == "position":
                vr_proj, h_proj = [vr_proj], [h_proj]
            for key, x, vr, h in zip(proj_return_keys, ctx.saved_tensors, vr_proj, h_proj):
                add_gradient_component(ctx.backward_info, key+'.gradr', vr)
                add_gradient_component(ctx.backward_info, key+'.hess', h)
                assert hasattr(x, 'optim_info')
                x.optim_info['gradr'] = ctx.backward_info[key+'.gradr']
                x.optim_info['hess'] = ctx.backward_info[key+'.hess']
                del ctx.backward_info[key+'.gradr']
                del ctx.backward_info[key+'.hess']
        else:
            _make_lazy_cuda_func(f"projection_{ctx.primitive}_backward")(
                (means, quats, scales, opacities, features_dc, features_sh),
                viewmats, intrins, ctx.width, ctx.height, ctx.camera_model_type, ctx.dist_coeffs,
                camera_ids, gaussian_ids, aabb,
                [x.contiguous() for x in v_proj_returns],
                v_proj, v_viewmats,
            )
        if not ctx.needs_input_grad[7]:
            v_viewmats = None

        return (
            None,
            *[None if is_accum else v for is_accum, v in zip(grad_accum_mask, v_proj)],
            v_viewmats,
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
        packed: bool,
        backward_info: Optional[dict] = None,
    ) -> Tuple[Tensor, Tensor, Tensor, Tensor, Tensor]:

        camera_model_type = camera_model.upper()

        cuda_returns = _make_lazy_cuda_func(
            f"projection_opaque_triangle_packed_forward" if packed else
            f"projection_opaque_triangle_forward"
        )(
            (means, quats, scales, hardness, features_dc, features_sh, features_ch),
            # (means, hardness),
            viewmats, intrins, width, height,
            near_plane, far_plane,
            camera_model_type, dist_coeffs
        )
        if packed:
            camera_ids, gaussian_ids, aabb, proj_returns = cuda_returns
            aabb.camera_ids = camera_ids
            aabb.gaussian_ids = gaussian_ids
        else:
            aabb, proj_returns = cuda_returns
            camera_ids, gaussian_ids = None, None

        ctx.save_for_backward(
            means, quats, scales, hardness, features_dc, features_sh, features_ch,
            viewmats, intrins, camera_ids, gaussian_ids, aabb
        )
        # ctx.save_for_backward(means, hardness, viewmats, intrins, aabb)
        ctx.width = width
        ctx.height = height
        ctx.camera_model_type = camera_model_type
        ctx.dist_coeffs = dist_coeffs
        ctx.backward_info = backward_info

        return aabb, *proj_returns

    @staticmethod
    def backward(ctx, v_aabb, *v_proj_returns):
        (
            means, quats, scales, hardness, features_dc, features_sh, features_ch,
            viewmats, intrins, camera_ids, gaussian_ids, aabb
        ) = ctx.saved_tensors

        params = (means, quats, scales, hardness, features_dc, features_sh, features_ch)
        grad_accum_mask = [
            x.is_leaf and hasattr(x, 'grad') and x.grad is not None
            for x in params
        ]
        v_proj = (*[
            x.grad if is_accum else torch.zeros_like(x)
            for is_accum, x in zip(grad_accum_mask, params)
        ],)
        v_viewmats = torch.zeros_like(viewmats) if ctx.needs_input_grad[7] else None

        _make_lazy_cuda_func("projection_opaque_triangle_backward")(
            params,
            # (means, hardness),
            viewmats, intrins, ctx.width, ctx.height, ctx.camera_model_type, ctx.dist_coeffs,
            camera_ids, gaussian_ids, aabb,
            [x.contiguous() for x in v_proj_returns],
            v_proj, v_viewmats,
        )
        if not ctx.needs_input_grad[7]:
            v_viewmats = None
        return (
            *[None if is_accum else v for is_accum, v in zip(grad_accum_mask, v_proj)],
            v_viewmats,
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
        packed: bool,
        backward_info: Optional[dict] = None,
    ) -> Tuple[Tensor, Tensor, Tensor, Tensor, Tensor]:

        camera_model_type = camera_model.upper()

        cuda_returns = _make_lazy_cuda_func(
            "projection_voxel_packed_forward" if packed else
            "projection_voxel_forward"
        )(
            (pos_sizes, densities, features_dc, features_sh),
            viewmats, intrins, width, height,
            near_plane, far_plane,
            camera_model_type, dist_coeffs
        )
        if packed:
            camera_ids, gaussian_ids, aabb, (depths, proj_rgbs) = cuda_returns
            aabb.camera_ids = camera_ids
            aabb.gaussian_ids = gaussian_ids
        else:
            aabb, (depths, proj_rgbs) = cuda_returns
            camera_ids, gaussian_ids = None, None

        ctx.save_for_backward(
            pos_sizes, densities, features_dc, features_sh,
            viewmats, intrins, camera_ids, gaussian_ids, aabb
        )
        ctx.width = width
        ctx.height = height
        ctx.camera_model_type = camera_model_type
        ctx.dist_coeffs = dist_coeffs
        ctx.backward_info = backward_info

        return aabb, depths, proj_rgbs

    @staticmethod
    def backward(ctx, v_aabb, v_depths, v_proj_rgbs):
        (
            pos_sizes, densities, features_dc, features_sh,
            viewmats, intrins, camera_ids, gaussian_ids, aabb
        ) = ctx.saved_tensors

        (v_pos_sizes, v_densities, v_features_dc, v_features_sh), v_viewmats = _make_lazy_cuda_func(
            "projection_voxel_backward"
        )(
            (pos_sizes, densities, features_dc, features_sh),
            viewmats, intrins, ctx.width, ctx.height, ctx.camera_model_type, ctx.dist_coeffs,
            camera_ids, gaussian_ids, aabb,
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
    in_splats: tuple[Tensor],  # means, quats, scales, opacities
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
    batch_dims = in_splats[0].shape[:-2]
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
        camera_model_type = camera_model.upper()

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

        camera_model_type = camera_model.upper()

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
    camera_model = camera_model.upper()
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
