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



def rasterize_to_pixels(
    primitive: Literal["3dgs", "mip", "opaque_triangle"],
    eval3d: bool,
    splats: tuple[Tensor],  # means, quats, scales, opacities
    # colors: Tensor,  # [..., N, channels] or [nnz, channels]
    image_width: int,
    image_height: int,
    tile_size: int,
    isect_offsets: Tensor,  # [..., tile_height, tile_width]
    flatten_ids: Tensor,  # [n_isects]
    backgrounds: Optional[Tensor] = None,  # [..., channels]
    masks: Optional[Tensor] = None,  # [..., tile_height, tile_width]
    packed: bool = False,
    absgrad: bool = False,
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

    tile_height, tile_width = isect_offsets.shape[-2:]
    assert (
        tile_height * tile_size >= image_height
    ), f"Assert Failed: {tile_height} * {tile_size} >= {image_height}"
    assert (
        tile_width * tile_size >= image_width
    ), f"Assert Failed: {tile_width} * {tile_size} >= {image_width}"

    if not eval3d:
        if primitive in ["3dgs", "mip"]:
            _RasterizeToPixels = _RasterizeToPixels3DGS
        elif primitive in ["opaque_triangle"]:
            _RasterizeToPixels = _RasterizeToPixelsOpaqueTriangle
        kwargs = []
    else:
        assert not packed, "packed is not supported for eval3d"
        for key in ['viewmats', 'Ks']:
            assert key in kwargs, "Camera extrinsics and intrinsics must be provided for Eval3D"
        kwargs = [kwargs.get(key, None) for key in [
            'viewmats', 'Ks', 'camera_model',
            'dist_coeffs'
        ]]
        if primitive in ["3dgs", "mip"]:
            _RasterizeToPixels = _RasterizeToPixels3DGSEval3D
        elif primitive in ["opaque_triangle"]:
            _RasterizeToPixels = _RasterizeToPixelsOpaqueTriangleEval3D
    render_outputs = _RasterizeToPixels.apply(
        *[(x.contiguous() if x is not None else None) for x in splats],
        backgrounds,
        masks,
        image_width,
        image_height,
        tile_size,
        isect_offsets.contiguous(),
        flatten_ids.contiguous(),
        absgrad,
        *kwargs
    )

    if eval3d:  # (rgbd,), alpha, last_id, (distortion,)
        meta = {}
        if primitive in ["opaque_triangle"]:
            meta["max_blending"] = render_outputs[-1]
            render_outputs = render_outputs[:-1]
        n_dist = len(render_outputs) // 2
        for tensor, name in zip(render_outputs[-n_dist:], ['rgb', 'depth', 'normal']):
            if tensor is not None:
                meta[f"{name}_distortion"] = tensor
        return tuple(render_outputs[:n_dist]), render_outputs[n_dist], meta

    # (rgbd,), alpha, last_id
    return tuple(render_outputs[:-1]), render_outputs[-1], {}



class _RasterizeToPixels3DGS(torch.autograd.Function):
    """Rasterize gaussians"""

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
        tile_size: int,
        isect_offsets: Tensor,  # [..., tile_height, tile_width]
        flatten_ids: Tensor,  # [n_isects]
        absgrad: bool,
    ) -> Tuple[Tensor, Tensor]:
        (render_rgbs, render_depths), render_Ts, last_ids = _make_lazy_cuda_func(
            "rasterization_3dgs_forward"
        )(
            (means2d, depths, conics, opacities, colors), backgrounds, masks,
            width, height, tile_size, isect_offsets, flatten_ids,
        )

        ctx.save_for_backward(
            means2d, depths, conics, colors, opacities, backgrounds, masks,
            isect_offsets, flatten_ids, render_Ts, last_ids,
        )
        ctx.width = width
        ctx.height = height
        ctx.tile_size = tile_size
        ctx.absgrad = absgrad

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
        tile_size = ctx.tile_size
        absgrad = ctx.absgrad

        (
            (v_means2d, v_depths, v_conics, v_opacities, v_colors),
            v_means2d_abs,
        ) = _make_lazy_cuda_func("rasterization_3dgs_backward")(
            (means2d, depths, conics, opacities, colors), backgrounds, masks,
            width, height, tile_size, isect_offsets, flatten_ids, render_Ts, last_ids,
            (v_render_rgbs.contiguous(), v_render_depths.contiguous()),
            v_render_alphas.contiguous(), absgrad,
        )
        if absgrad:
            means2d.absgrad = v_means2d_abs

        v_backgrounds = None
        if ctx.needs_input_grad[5]:
            v_backgrounds = (torch.cat([v_render_rgbs, v_render_depths], dim=-1) * \
                             render_Ts.float()).sum(dim=(-3, -2))

        return (
            v_means2d, v_depths, v_conics, v_opacities, v_colors, v_backgrounds,
            *([None]*(len(ctx.needs_input_grad)-6))
        )


class _RasterizeToPixelsOpaqueTriangle(torch.autograd.Function):
    """Rasterize opaque triangles"""

    @staticmethod
    def forward(
        ctx,
        means2d: Tensor,  # [..., N, 3, 2] or [nnz, 3, 2]
        depths: Tensor,  # [..., N, 1] or [nnz, 1]
        hardness: Tensor,  # [..., N] or [nnz]
        colors: Tensor,  # [..., N, 3] or [nnz, 3]
        normals: Tensor,  # [..., N, 3] or [nnz, 3]
        backgrounds: Tensor,  # [..., channels], Optional
        masks: Tensor,  # [..., tile_height, tile_width], Optional
        width: int,
        height: int,
        tile_size: int,
        isect_offsets: Tensor,  # [..., tile_height, tile_width]
        flatten_ids: Tensor,  # [n_isects]
        absgrad: bool,
    ) -> Tuple[Tensor, Tensor]:
        (render_rgbs, render_depths, render_normals), render_Ts, last_ids = _make_lazy_cuda_func(
            "rasterization_opaque_triangle_forward"
        )(
            (means2d, depths, hardness, colors, normals), backgrounds, masks,
            width, height, tile_size, isect_offsets, flatten_ids,
        )

        ctx.save_for_backward(
            means2d, depths, hardness, colors, normals, backgrounds, masks,
            isect_offsets, flatten_ids, render_Ts, last_ids,
        )
        ctx.width = width
        ctx.height = height
        ctx.tile_size = tile_size
        ctx.absgrad = absgrad

        render_alphas = 1.0 - render_Ts
        return render_rgbs, render_depths, render_normals, render_alphas

    @staticmethod
    def backward(
        ctx,
        v_render_rgbs: Tensor,  # [..., H, W, 3]
        v_render_depths: Tensor,  # [..., H, W, 1]
        v_render_normals: Tensor,  # [..., H, W, 3]
        v_render_alphas: Tensor,  # [..., H, W, 1]
    ):
        (
            means2d, depths, hardness, colors, normals, backgrounds, masks,
            isect_offsets, flatten_ids, render_Ts, last_ids,
        ) = ctx.saved_tensors
        width = ctx.width
        height = ctx.height
        tile_size = ctx.tile_size
        absgrad = ctx.absgrad

        (
            (v_means2d, v_depths, v_hardness, v_rgbs, v_normals), v_means2d_abs,
        ) = _make_lazy_cuda_func("rasterization_opaque_triangle_backward")(
            (means2d, depths, hardness, colors, normals), backgrounds, masks,
            width, height, tile_size, isect_offsets, flatten_ids, render_Ts, last_ids,
            (v_render_rgbs.contiguous(), v_render_depths.contiguous(), v_render_normals.contiguous()),
            v_render_alphas.contiguous(), absgrad,
        )
        if absgrad:
            means2d.absgrad = v_means2d_abs
            # means2d.absgrad = v_means2d.mean(-2)

        v_backgrounds = None
        if ctx.needs_input_grad[5]:
            v_backgrounds = (torch.cat([v_render_rgbs, v_render_depths, v_render_normals], dim=-1) * \
                             render_Ts.float()).sum(dim=(-3, -2))

        return (
            v_means2d, v_depths, v_hardness, v_rgbs, v_normals, v_backgrounds,
            *([None]*(len(ctx.needs_input_grad)-6))
        )


class _RasterizeToPixels3DGSEval3D(torch.autograd.Function):
    """Rasterize gaussians"""

    @staticmethod
    def forward(
        ctx,
        # gauss params
        means: Tensor,  # [..., N, 3]
        quats: Tensor,  # [..., N, 4]
        scales: Tensor,  # [..., N, 3]
        opacities: Tensor,  # [..., N]
        features_dc: Tensor,  # [..., N, 3]
        features_sh: Tensor,  # [..., N, x, 3]
        # proj outputs
        depths: Tensor,  # [..., N, 1] or [nnz, 1]
        proj_scales: Tensor,  # [..., N, 3]
        proj_opacities: Tensor,  # [..., N]
        colors: Tensor,  # [..., N, channels] or [nnz, channels]
        # rest
        backgrounds: Tensor,  # [..., channels], Optional
        masks: Tensor,  # [..., tile_height, tile_width], Optional
        width: int,
        height: int,
        tile_size: int,
        isect_offsets: Tensor,  # [..., tile_height, tile_width]
        flatten_ids: Tensor,  # [n_isects]
        absgrad: bool,
        viewmats: Tensor,
        Ks: Tensor,
        camera_model: Literal["pinhole", "ortho", "fisheye", "ftheta"],
        dist_coeffs: Optional[Tensor],
    ) -> Tuple[Tensor, Tensor]:
        if absgrad:
            raise NotImplementedError("Absgrad not supported with Eval3D")

        camera_model = gsplat.cuda._wrapper._make_lazy_cuda_obj(
            f"CameraModelType.{camera_model.upper()}"
        )

        (render_rgbs, render_depths), render_Ts, last_ids, _1, _2 = _make_lazy_cuda_func(
            "rasterization_3dgs_eval3d_forward"
        )(
            (means, quats, depths, proj_scales, proj_opacities, colors),
            viewmats, Ks, camera_model, dist_coeffs,
            backgrounds, masks,
            width, height, tile_size, isect_offsets, flatten_ids,
        )

        ctx.save_for_backward(
            means, quats, depths, proj_scales, proj_opacities, colors,
            viewmats, Ks, dist_coeffs,
            backgrounds, masks,
            isect_offsets, flatten_ids, render_Ts, last_ids,
        )
        ctx.width = width
        ctx.height = height
        ctx.camera_model = camera_model
        ctx.tile_size = tile_size

        render_alphas = 1.0 - render_Ts
        return render_rgbs, render_depths, render_alphas, None, None

    @staticmethod
    def backward(
        ctx,
        v_render_rgbs: Tensor,  # [..., H, W, 3]
        v_render_depths: Tensor,  # [..., H, W, 1]
        v_render_alphas: Tensor,  # [..., H, W, 1]
        v_render_rgb_distortions: Optional[Tensor],
        v_render_depth_distortions: Optional[Tensor]
    ):
        (
            means, quats, depths, proj_scales, proj_opacities, colors,
            viewmats, Ks, dist_coeffs,
            backgrounds, masks,
            isect_offsets, flatten_ids, render_Ts, last_ids,
        ) = ctx.saved_tensors
        width = ctx.width
        height = ctx.height
        tile_size = ctx.tile_size

        (
            (v_means, v_quats, v_depths, v_proj_scales, v_proj_opacities, v_colors),
            v_viewmats,
        ) = _make_lazy_cuda_func("rasterization_3dgs_eval3d_backward")(
            (means, quats, depths, proj_scales, proj_opacities, colors),
            viewmats, Ks, ctx.camera_model, dist_coeffs,
            backgrounds, masks,
            width, height, tile_size, isect_offsets, flatten_ids,
            render_Ts, last_ids, None, None,
            (v_render_rgbs.contiguous(), v_render_depths.contiguous()),
            v_render_alphas.contiguous(), None
        )

        v_backgrounds = None
        if ctx.needs_input_grad[10]:
            v_backgrounds = (torch.cat([v_render_rgbs, v_render_depths], dim=-1) * \
                             render_Ts.float()).sum(dim=(-3, -2))

        return (
            v_means, v_quats, None, None, None, None,
            v_depths, v_proj_scales, v_proj_opacities, v_colors, v_backgrounds,
            *([None]*(len(ctx.needs_input_grad)-11))
        )


class _RasterizeToPixelsOpaqueTriangleEval3D(torch.autograd.Function):
    """Rasterize opaque triangles"""

    @staticmethod
    def forward(
        ctx,
        # gauss params
        means: Tensor,  # [..., N, 3]
        quats: Tensor,  # [..., N, 4]
        scales: Tensor,  # [..., N, 3]
        hardness: Tensor,  # [..., N]
        features_dc: Tensor,  # [..., N, 3]
        features_sh: Tensor,  # [..., N, x, 3]
        features_ch: Tensor,  # [..., N, x, 3]
        # proj outputs
        depths: Tensor,  # [..., N, 1] or [nnz, 1]
        verts: Tensor,  # [..., N, 3, 3] or [nnz, 3, 3]
        rgbs: Tensor,  # [..., N, 3, 3] or [nnz, 3, 3]
        normals: Tensor,  # [..., N, 3] or [nnz, 3]
        # rest
        backgrounds: Tensor,  # [..., channels], Optional
        masks: Tensor,  # [..., tile_height, tile_width], Optional
        width: int,
        height: int,
        tile_size: int,
        isect_offsets: Tensor,  # [..., tile_height, tile_width]
        flatten_ids: Tensor,  # [n_isects]
        absgrad: bool,
        viewmats: Tensor,
        Ks: Tensor,
        camera_model: Literal["pinhole", "ortho", "fisheye", "ftheta"],
        dist_coeffs: Optional[Tensor],
    ) -> Tuple[Tensor, Tensor]:
        if absgrad:
            raise NotImplementedError("Absgrad not supported with Eval3D")

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
        ) = _make_lazy_cuda_func("rasterization_opaque_triangle_eval3d_forward")(
            (hardness, depths, verts, rgbs, normals),
            viewmats, Ks, camera_model, dist_coeffs,
            backgrounds, masks,
            width, height, tile_size, isect_offsets, flatten_ids,
        )
        # torch.cuda.synchronize()
        # time1 = perf_counter()
        # print(f"fwd: {1e3*(time1-time0):.2f} ms")

        ctx.save_for_backward(
            hardness, depths, verts, rgbs, normals,
            viewmats, Ks, dist_coeffs,
            backgrounds, masks,
            isect_offsets, flatten_ids, render_Ts, last_ids,
            render_rgbs, render_depths, render_normals,
            render2_rgbs, render2_depths, render2_normals
        )
        ctx.width = width
        ctx.height = height
        ctx.camera_model = camera_model
        ctx.tile_size = tile_size

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
        (
            hardness, depths, verts, rgbs, normals,
            viewmats, Ks, dist_coeffs,
            backgrounds, masks,
            isect_offsets, flatten_ids, render_Ts, last_ids,
            render_rgbs, render_depths, render_normals,
            render2_rgbs, render2_depths, render2_normals
        ) = ctx.saved_tensors
        width = ctx.width
        height = ctx.height
        tile_size = ctx.tile_size

        # from time import perf_counter
        # torch.cuda.synchronize()
        # time0 = perf_counter()
        (
            (v_hardness, v_depths, v_verts, v_rgbs, v_normals),
            v_viewmats,
        ) = _make_lazy_cuda_func("rasterization_opaque_triangle_eval3d_backward")(
            (hardness, depths, verts, rgbs, normals),
            viewmats, Ks, ctx.camera_model, dist_coeffs,
            backgrounds, masks,
            width, height, tile_size, isect_offsets, flatten_ids, render_Ts, last_ids,
            (render_rgbs, render_depths, render_normals),
            (render2_rgbs, render2_depths, render2_normals),
            (v_render_rgbs.contiguous(), v_render_depths.contiguous(), v_render_normals.contiguous()),
            v_render_alphas.contiguous(),
            (v_distortion_rgbs.contiguous(), v_distortion_depths.contiguous(), v_distortion_normals.contiguous()),
        )
        # torch.cuda.synchronize()
        # time1 = perf_counter()
        # print(f"bwd: {1e3*(time1-time0):.2f} ms")

        v_backgrounds = None
        if ctx.needs_input_grad[11]:
            v_backgrounds = (torch.cat([v_render_rgbs, v_render_depths, v_render_normals], dim=-1) * \
                             render_Ts.float()).sum(dim=(-3, -2))

        return (
            None, None, None, v_hardness, None, None, None, v_depths, v_verts, v_rgbs, v_normals, v_backgrounds,
            *([None]*(len(ctx.needs_input_grad)-12))
        )
