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



def rasterize_to_pixels(
    primitive: Literal["3dgs", "mip", "3dgut", "opaque_triangle", "voxel"],
    splats: tuple[Tensor],  # means, quats, scales, opacities
    # colors: Tensor,  # [..., N, channels] or [nnz, channels]
    image_width: int,
    image_height: int,
    isect_offsets: Tensor,  # [..., tile_height, tile_width]
    flatten_ids: Tensor,  # [n_isects]
    backgrounds: Optional[Tensor] = None,  # [..., channels]
    masks: Optional[Tensor] = None,  # [..., image_height, image_width]
    packed: bool = False,
    absgrad: bool = False,
    output_distortion: bool = False,
    compute_hessian_diagonal: bool = False,
    **kwargs
) -> Tuple[Tensor, Tensor]:
    # TODO
    # image_dims = colors.shape[:-2]
    # channels = colors.shape[-1]
    # device = colors.device
    # if packed:
    #     nnz = colors.size(0)
    #     assert colors.shape[0] == nnz, colors.shape
    #     # assert means2d.shape == (nnz, 2), means2d.shape
    #     # assert conics.shape == (nnz, 3), conics.shape
    #     # assert opacities.shape == (nnz,), opacities.shape
    # else:
    #     N = colors.size(-2)
    #     assert colors.shape == image_dims + (N, channels), colors.shape
    #     # assert means2d.shape == image_dims + (N, 2), means2d.shape
    #     # assert conics.shape == image_dims + (N, 3), conics.shape
    #     # assert opacities.shape == image_dims + (N,), opacities.shape
    # if backgrounds is not None:
    #     assert backgrounds.shape == image_dims + (channels,), backgrounds.shape
    #     backgrounds = backgrounds.contiguous()
    # if masks is not None:
    #     assert masks.shape == isect_offsets.shape, masks.shape
    #     masks = masks.contiguous()

    assert backgrounds is None, "TODO"

    additional_args = [kwargs['viewmats'], kwargs['intrins'], kwargs['camera_model'], kwargs['dist_coeffs']]
    if primitive in ["3dgs", "mip"]:
        _RasterizeToPixels = _RasterizeToPixels3DGS
        additional_args = [primitive == "mip", compute_hessian_diagonal]
    if primitive in ["3dgut", "3dgut_sv"]:
        _RasterizeToPixels = _RasterizeToPixels3DGUT
        additional_args += [output_distortion, compute_hessian_diagonal]
    elif primitive in ["opaque_triangle"]:
        _RasterizeToPixels = _RasterizeToPixelsOpaqueTriangle
    elif primitive in ["voxel"]:
        _RasterizeToPixels = _RasterizeToPixelsVoxelEval3D

    render_outputs = _RasterizeToPixels.apply(
        *[(x.contiguous() if x is not None else None) for x in splats],
        backgrounds,
        masks,
        image_width,
        image_height,
        isect_offsets.contiguous(),
        flatten_ids.contiguous(),
        *additional_args
    )

    if primitive in ["3dgs", "mip"]:
        render_rgbs, render_depths, render_alphas = render_outputs
        return (render_rgbs, render_depths), render_alphas, {}
    elif primitive in ["3dgut", "3dgut_sv"]:
        (
            render_rgbs, render_depths, render_alphas,
            distortion_rgbs, distortion_depths
        ) = render_outputs
        if output_distortion:
            return (render_rgbs, render_depths), render_alphas, {
                'rgb_distortion': distortion_rgbs,
                'depth_distortion': distortion_depths,
            }
        return (render_rgbs, render_depths), render_alphas, {}
    elif primitive in ["opaque_triangle"]:
        (
            render_rgbs, render_depths, render_normals, render_alphas,
            distortion_rgbs, distortion_depths, distortion_normals, max_blending
        ) = render_outputs
        return (
            (render_rgbs, render_depths, render_normals), render_alphas,
            {
                'max_blending': max_blending,
                'rgb_distortion': distortion_rgbs,
                'depth_distortion': distortion_depths,
                'normal_distortion': distortion_normals,
            }
        )
    elif primitive in ["voxel"]:
        (
            render_rgbs, render_depths, render_alphas,
            distortion_rgbs, distortion_depths, max_blending
        ) = render_outputs
        return (
            (render_rgbs, render_depths), render_alphas,
            {
                'max_blending': max_blending,
                'rgb_distortion': distortion_rgbs,
                'depth_distortion': distortion_depths,
            }
        )
    else:
        raise NotImplementedError()



class _RasterizeToPixels3DGS(torch.autograd.Function):

    @staticmethod
    def forward(
        ctx,
        means2d: Tensor,  # [..., N, 2] or [nnz, 2]
        depths: Tensor,  # [..., N, 1] or [nnz, 1]
        conics: Tensor,  # [..., N, 3] or [nnz, 3]
        opacities: Tensor,  # [..., N] or [nnz]
        colors: Tensor,  # [..., N, channels] or [nnz, channels]
        backgrounds: Tensor,  # [..., channels], Optional
        masks: Tensor,  # [..., tile_height, tile_width], Optional
        width: int,
        height: int,
        isect_offsets: Tensor,  # [..., tile_height, tile_width]
        flatten_ids: Tensor,  # [n_isects]
        antialiased: bool,
        compute_hessian_diagonal: bool = False
    ) -> Tuple[Tensor, Tensor]:
        (render_rgbs, render_depths), render_Ts, last_ids = _make_lazy_cuda_func(
            f"rasterization_{['3dgs', 'mip'][antialiased]}_forward"
        )(
            (means2d, depths, conics, opacities, colors, None), backgrounds, masks,
            width, height, isect_offsets, flatten_ids,
        )

        ctx.save_for_backward(
            means2d, depths, conics, colors, opacities, backgrounds, masks,
            isect_offsets, flatten_ids, render_Ts, last_ids,
        )
        ctx.width = width
        ctx.height = height
        ctx.antialiased = antialiased
        ctx.compute_hessian_diagonal = compute_hessian_diagonal

        render_alphas = 1.0 - render_Ts
        return render_rgbs, render_depths, render_alphas

    @staticmethod
    def backward(
        ctx,
        v_render_rgbs: Tensor,  # [..., H, W, 3]
        v_render_depths: Tensor,  # [..., H, W, 1]
        v_render_alphas: Tensor,  # [..., H, W, 1]
    ):
        (
            means2d, depths, conics, colors, opacities, backgrounds, masks,
            isect_offsets, flatten_ids, render_Ts, last_ids,
        ) = ctx.saved_tensors
        width = ctx.width
        height = ctx.height
        if ctx.compute_hessian_diagonal:
            assert hasattr(v_render_rgbs, 'loss_map')

        if ctx.compute_hessian_diagonal:
            v_splats, vr_splats, h_splats = _make_lazy_cuda_func(
                f"rasterization_{['3dgs', 'mip'][ctx.antialiased]}_backward_with_hessian_diagonal"
            )(
                (means2d, depths, conics, opacities, colors, None), backgrounds, masks,
                width, height, isect_offsets, flatten_ids, render_Ts, last_ids,
                None, None,
                v_render_rgbs.loss_map if ctx.compute_hessian_diagonal else None,
                (v_render_rgbs.contiguous(), v_render_depths.contiguous()),
                v_render_alphas.contiguous(),
                None
            )
            v_means2d, v_depths, v_conics, v_opacities, v_colors, v_means2d_abs = v_splats
            for v, vr, h in zip(v_splats, vr_splats, h_splats):
                v.gradr = vr
                v.hess = h
        else:
            (
                v_means2d, v_depths, v_conics, v_opacities, v_colors, v_means2d_abs
            ) = _make_lazy_cuda_func(
                f"rasterization_{['3dgs', 'mip'][ctx.antialiased]}_backward"
            )(
                (means2d, depths, conics, opacities, colors, None), backgrounds, masks,
                width, height, isect_offsets, flatten_ids, render_Ts, last_ids,
                (v_render_rgbs.contiguous(), v_render_depths.contiguous()),
                v_render_alphas.contiguous(),
            )
        if v_means2d_abs is not None:
            means2d.absgrad = v_means2d_abs

        v_backgrounds = None
        if ctx.needs_input_grad[5]:
            v_backgrounds = (torch.cat([v_render_rgbs, v_render_depths], dim=-1) * \
                             render_Ts.float()).sum(dim=(-3, -2))

        return (
            v_means2d, v_depths, v_conics, v_opacities, v_colors, v_backgrounds,
            *([None]*(len(ctx.needs_input_grad)-6))
        )


class _RasterizeToPixels3DGUT(torch.autograd.Function):

    @staticmethod
    def forward(
        ctx,
        # gauss params
        means: Tensor,  # [..., N, 3]
        quats: Tensor,  # [..., N, 4]
        depths: Tensor,  # [..., N, 1] or [nnz, 1]
        proj_scales: Tensor,  # [..., N, 3]
        proj_opacities: Tensor,  # [..., N]
        colors: Tensor,  # [..., N, channels] or [nnz, channels]
        # rest
        backgrounds: Tensor,  # [..., channels], Optional
        masks: Tensor,  # [..., tile_height, tile_width], Optional
        width: int,
        height: int,
        isect_offsets: Tensor,  # [..., tile_height, tile_width]
        flatten_ids: Tensor,  # [n_isects]
        viewmats: Tensor,
        intrins: Tensor,
        camera_model: Literal["pinhole", "fisheye"],
        dist_coeffs: Optional[Tensor],
        output_distortion: bool = False,
        compute_hessian_diagonal: bool = False
    ) -> Tuple[Tensor, Tensor]:

        camera_model = gsplat.cuda._wrapper._make_lazy_cuda_obj(
            f"CameraModelType.{camera_model.upper()}"
        )

        (
            (render_rgbs, render_depths), render_Ts, last_ids, 
            render2_outputs, distortion_outputs, max_blending
        ) = _make_lazy_cuda_func(
            "rasterization_3dgut_forward"
        )(
            (means, quats, depths, proj_scales, proj_opacities, colors),
            viewmats, intrins, camera_model, dist_coeffs,
            backgrounds, masks,
            width, height, isect_offsets, flatten_ids,
            output_distortion
        )

        if output_distortion:
            render2_rgbs, render2_depths = render2_outputs
            distortion_rgbs, distortion_depths = distortion_outputs
        else:
            render2_rgbs, render2_depths, distortion_rgbs, distortion_depths = None, None, None, None
        assert max_blending is None

        ctx.save_for_backward(
            means, quats, depths, proj_scales, proj_opacities, colors,
            viewmats, intrins, dist_coeffs,
            backgrounds, masks,
            isect_offsets, flatten_ids, render_Ts, last_ids,
            *(
                [render_rgbs, render_depths, render2_rgbs, render2_depths]
                if output_distortion else [None]*4
            )
        )
        ctx.width = width
        ctx.height = height
        ctx.camera_model = camera_model
        ctx.output_distortion = output_distortion
        ctx.compute_hessian_diagonal = compute_hessian_diagonal

        render_alphas = 1.0 - render_Ts
        return render_rgbs, render_depths, render_alphas, distortion_rgbs, distortion_depths

    @staticmethod
    def backward(
        ctx,
        v_render_rgbs: Tensor,  # [..., H, W, 3]
        v_render_depths: Tensor,  # [..., H, W, 1]
        v_render_alphas: Tensor,  # [..., H, W, 1]
        v_distortion_rgbs: Tensor,  # [..., H, W, 3]
        v_distortion_depths: Tensor,  # [..., H, W, 1]
    ):
        (
            means, quats, depths, proj_scales, proj_opacities, colors,
            viewmats, intrins, dist_coeffs,
            backgrounds, masks,
            isect_offsets, flatten_ids, render_Ts, last_ids,
            render_rgbs, render_depths, render2_rgbs, render2_depths,
        ) = ctx.saved_tensors
        width = ctx.width
        height = ctx.height
        if ctx.compute_hessian_diagonal:
            assert hasattr(v_render_rgbs, 'loss_map')

        cuda_return = _make_lazy_cuda_func("rasterization_3dgut_backward" + "_with_hessian_diagonal"*ctx.compute_hessian_diagonal)(
            (means, quats, depths, proj_scales, proj_opacities, colors),
            viewmats, intrins, ctx.camera_model, dist_coeffs,
            backgrounds, masks,
            width, height, isect_offsets, flatten_ids, render_Ts, last_ids,
            (render_rgbs, render_depths) if ctx.output_distortion else None,
            (render2_rgbs, render2_depths) if ctx.output_distortion else None,
            v_render_rgbs.loss_map if ctx.compute_hessian_diagonal else None,
            (v_render_rgbs.contiguous(), v_render_depths.contiguous()),
            v_render_alphas.contiguous(),
            (v_distortion_rgbs.contiguous(), v_distortion_depths.contiguous()) if ctx.output_distortion else None,
            ctx.needs_input_grad[13]
        )

        h_splats = None
        if ctx.compute_hessian_diagonal:
            v_splats, v_viewmats, vr_splats, h_splats = cuda_return
        else:
            v_splats, v_viewmats = cuda_return
        (v_means, v_quats, v_depths, v_proj_scales, v_proj_opacities, v_colors) = v_splats
        # print(v_means.mean().item())

        if h_splats is not None:
            for v, vr, h in zip(v_splats, vr_splats, h_splats):
                add_gradient_component(v, 'gradr', vr)
                add_gradient_component(v, 'hess', h)
                v.gradr_all = vr_splats
                v.hess_all = h_splats  # so we can get all hess from projection backward pass

        v_backgrounds = None
        if ctx.needs_input_grad[6]:
            v_backgrounds = (torch.cat([v_render_rgbs, v_render_depths], dim=-1) * \
                             render_Ts.float()).sum(dim=(-3, -2))

        return (
            v_means, v_quats, v_depths, v_proj_scales, v_proj_opacities, v_colors, v_backgrounds,
            None, None, None, None, None, None, v_viewmats, None, None, None, None, None
        )


class _RasterizeToPixelsOpaqueTriangle(torch.autograd.Function):

    @staticmethod
    def forward(
        ctx,
        # gauss params
        hardness: Tensor,  # [..., N]
        # proj outputs
        depths: Tensor,  # [..., N, 1] or [nnz, 1]
        verts: Tensor,  # [..., N, 3, 3] or [nnz, 3, 3]
        rgbs: Tensor,  # [..., N, 3, 3] or [nnz, 3, 3]
        normals: Tensor,  # [..., N, 3] or [nnz, 3]
        # rest
        backgrounds: Optional[Tensor],  # [..., channels]
        max_blending_masks: Optional[Tensor],  # [..., image_height, image_width]
        width: int,
        height: int,
        isect_offsets: Tensor,  # [..., tile_height, tile_width]
        flatten_ids: Tensor,  # [n_isects]
        viewmats: Tensor,
        intrins: Tensor,
        camera_model: Literal["pinhole", "ortho", "fisheye", "ftheta"],
        dist_coeffs: Optional[Tensor],
    ) -> Tuple[Tensor, Tensor]:

        camera_model = gsplat.cuda._wrapper._make_lazy_cuda_obj(
            f"CameraModelType.{camera_model.upper()}"
        )

        # from time import perf_counter
        # torch.cuda.synchronize()
        # time0 = perf_counter()
        (
            (render_rgbs, render_depths, render_normals),
            render_Ts, last_ids,
            (render2_rgbs, render2_depths, render2_normals),
            (distortion_rgbs, distortion_depths, distortion_normals),
            max_blending
        ) = _make_lazy_cuda_func("rasterization_opaque_triangle_forward")(
            (hardness, depths, verts, rgbs, normals),
            viewmats, intrins, camera_model, dist_coeffs,
            backgrounds, max_blending_masks,
            width, height, isect_offsets, flatten_ids,
        )
        # torch.cuda.synchronize()
        # time1 = perf_counter()
        # print(f"fwd: {1e3*(time1-time0):.2f} ms")

        ctx.save_for_backward(
            hardness, depths, verts, rgbs, normals,
            viewmats, intrins, dist_coeffs,
            backgrounds,
            isect_offsets, flatten_ids, render_Ts, last_ids,
            render_rgbs, render_depths, render_normals,
            render2_rgbs, render2_depths, render2_normals
        )
        ctx.width = width
        ctx.height = height
        ctx.camera_model = camera_model

        render_alphas = 1.0 - render_Ts
        return (
            render_rgbs, render_depths, render_normals, render_alphas,
            distortion_rgbs, distortion_depths, distortion_normals, max_blending
        )

    @staticmethod
    def backward(
        ctx,
        v_render_rgbs: Tensor,  # [..., H, W, 3]
        v_render_depths: Tensor,  # [..., H, W, 1]
        v_render_normals: Tensor,  # [..., H, W, 3]
        v_render_alphas: Tensor,  # [..., H, W, 1]
        v_distortion_rgbs: Tensor,  # [..., H, W, 3]
        v_distortion_depths: Tensor,  # [..., H, W, 1]
        v_distortion_normals: Tensor,  # [..., H, W, 3]
        v_max_blending = None
    ):
        # assert v_max_blending is None, "max_blending does not support gradient"
        (
            hardness, depths, verts, rgbs, normals,
            viewmats, intrins, dist_coeffs,
            backgrounds,
            isect_offsets, flatten_ids, render_Ts, last_ids,
            render_rgbs, render_depths, render_normals,
            render2_rgbs, render2_depths, render2_normals
        ) = ctx.saved_tensors
        width = ctx.width
        height = ctx.height

        # from time import perf_counter
        # torch.cuda.synchronize()
        # time0 = perf_counter()
        (
            (v_hardness, v_depths, v_verts, v_rgbs, v_normals),
            v_viewmats,
        ) = _make_lazy_cuda_func("rasterization_opaque_triangle_backward")(
            (hardness, depths, verts, rgbs, normals),
            viewmats, intrins, ctx.camera_model, dist_coeffs,
            backgrounds,
            width, height, isect_offsets, flatten_ids, render_Ts, last_ids,
            (render_rgbs, render_depths, render_normals),
            (render2_rgbs, render2_depths, render2_normals),
            (v_render_rgbs.contiguous(), v_render_depths.contiguous(), v_render_normals.contiguous()),
            v_render_alphas.contiguous(),
            (v_distortion_rgbs.contiguous(), v_distortion_depths.contiguous(), v_distortion_normals.contiguous()),
            # ctx.needs_input_grad[12]
        )
        if ctx.needs_input_grad[12]:
            raise NotImplementedError("Camera optimizer not supported for opaque triangle")
        # torch.cuda.synchronize()
        # time1 = perf_counter()
        # print(f"bwd: {1e3*(time1-time0):.2f} ms")

        v_backgrounds = None
        if ctx.needs_input_grad[5]:
            v_backgrounds = (torch.cat([v_render_rgbs, v_render_depths, v_render_normals], dim=-1) * \
                             render_Ts.float()).sum(dim=(-3, -2))

        return (
            v_hardness, v_depths, v_verts, v_rgbs, v_normals, v_backgrounds,
            None, None, None, None, None, None, v_viewmats, None, None, None
        )

class _RasterizeToPixelsVoxelEval3D(torch.autograd.Function):
    """Rasterize gaussians"""

    @staticmethod
    def forward(
        ctx,
        # gauss params
        pos_sizes: Tensor,  # [..., N, 4]
        densities: Tensor,  # [..., N, 8]
        # proj outputs
        colors: Tensor,  # [..., N, channels] or [nnz, channels]
        # rest
        backgrounds: Tensor,  # [..., channels], Optional
        masks: Tensor,  # [..., tile_height, tile_width], Optional
        width: int,
        height: int,
        isect_offsets: Tensor,  # [..., tile_height, tile_width]
        flatten_ids: Tensor,  # [n_isects]
        viewmats: Tensor,
        intrins: Tensor,
        camera_model: Literal["pinhole", "ortho", "fisheye", "ftheta"],
        dist_coeffs: Optional[Tensor],
    ) -> Tuple[Tensor, Tensor]:

        camera_model = gsplat.cuda._wrapper._make_lazy_cuda_obj(
            f"CameraModelType.{camera_model.upper()}"
        )

        (
            (render_rgbs, render_depths), render_Ts, last_ids,
            (render2_rgbs, render2_depths),
            (distortion_rgbs, distortion_depths),
            max_blending
        ) = _make_lazy_cuda_func("rasterization_voxel_forward")(
            (pos_sizes, None, densities, colors),
            viewmats, intrins, camera_model, dist_coeffs,
            backgrounds, masks,
            width, height, isect_offsets, flatten_ids,
        )

        ctx.save_for_backward(
            pos_sizes, densities, colors,
            viewmats, intrins, dist_coeffs,
            backgrounds, masks,
            isect_offsets, flatten_ids, render_Ts, last_ids,
            render_rgbs, render_depths,
            render2_rgbs, render2_depths,
        )
        ctx.width = width
        ctx.height = height
        ctx.camera_model = camera_model

        render_alphas = 1.0 - render_Ts
        return (
            render_rgbs, render_depths, render_alphas,
            distortion_rgbs, distortion_depths, max_blending
        )

    @staticmethod
    def backward(
        ctx,
        v_render_rgbs: Tensor,  # [..., H, W, 3]
        v_render_depths: Tensor,  # [..., H, W, 1]
        v_render_alphas: Tensor,  # [..., H, W, 1]
        v_distortion_rgbs: Tensor,  # [..., H, W, 3]
        v_distortion_depths: Tensor,  # [..., H, W, 1]
        v_max_blending = None
    ):
        # assert v_max_blending is None, "max_blending does not support gradient"
        (
            pos_sizes, densities, colors,
            viewmats, intrins, dist_coeffs,
            backgrounds, masks,
            isect_offsets, flatten_ids, render_Ts, last_ids,
            render_rgbs, render_depths,
            render2_rgbs, render2_depths,
        ) = ctx.saved_tensors
        width = ctx.width
        height = ctx.height

        (
            (v_pos_sizes, v_depths, v_densities, v_colors),
            v_viewmats,
        ) = _make_lazy_cuda_func("rasterization_voxel_backward")(
            (pos_sizes, None, densities, colors),
            viewmats, intrins, ctx.camera_model, dist_coeffs,
            backgrounds, masks,
            width, height, isect_offsets, flatten_ids, render_Ts, last_ids,
            (render_rgbs, render_depths),
            (render2_rgbs, render2_depths),
            (v_render_rgbs.contiguous(), v_render_depths.contiguous()),
            v_render_alphas.contiguous(),
            (v_distortion_rgbs.contiguous(), v_distortion_depths.contiguous()),
            ctx.needs_input_grad[10]
        )
        assert v_pos_sizes is None
        assert v_depths is None
        assert v_densities is not None
        assert v_colors is not None

        v_backgrounds = None
        if ctx.needs_input_grad[3]:
            v_backgrounds = (torch.cat([v_render_rgbs, v_render_depths], dim=-1) * \
                             render_Ts.float()).sum(dim=(-3, -2))

        return (
            v_pos_sizes, v_densities,
            v_colors, v_backgrounds,
            None, None, None, None, None, None, v_viewmats, None, None, None
        )

