import torch
from torch import Tensor

from typing import Any, Callable, Literal, Optional, Union, Tuple

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




def blend_background(rgb, alpha, background) -> Tensor:
    return _BlendBackground.apply(rgb, alpha, background)


class _BlendBackground(torch.autograd.Function):

    @staticmethod
    def forward(
        ctx,
        rgb, alpha, background
    ):
        out_rgb = _make_lazy_cuda_func("blend_background_forward")(rgb, alpha, background)

        ctx.save_for_backward(rgb, alpha, background)

        return out_rgb

    @staticmethod
    def backward(ctx, v_out_rgb):

        rgb, alpha, background = ctx.saved_tensors

        return _make_lazy_cuda_func("blend_background_backward")(
            rgb, alpha, background,
            v_out_rgb
        )




def depth_to_normal(
    depths: Tensor,  # [B, H, W, 1]
    camera_model: Literal["pinhole", "fisheye"],
    Ks: Union[Tensor, Tuple[float, float, float, float]],
    dist_coeffs: Optional[Tensor] = None,  # [..., C, 10]
    is_ray_depth: bool = True,
) -> Tensor:
    if isinstance(Ks, tuple):
        fx, fy, cx, cy = Ks
        Ks = torch.tensor([[[fx, 0, cx], [0, fy, cy], [0, 0, 1]]])
        Ks = Ks.to(depths).repeat(len(depths), 1, 1)
    return _DepthToNormal.apply(
        depths,
        camera_model,
        Ks.contiguous(),
        dist_coeffs.contiguous() if dist_coeffs is not None else None,
        is_ray_depth
    )


class _DepthToNormal(torch.autograd.Function):
    """Projects Gaussians to 2D."""

    @staticmethod
    def forward(
        ctx,
        depths: Tensor,  # [B, H, W, 1]
        camera_model: Literal["pinhole", "fisheye"],
        Ks: Tensor,
        dist_coeffs: Optional[Tensor],  # [..., C, 10]
        is_ray_depth: bool
    ) -> Tensor:

        camera_model_type = gsplat.cuda._wrapper._make_lazy_cuda_obj(
            f"CameraModelType.{camera_model.upper()}"
        )

        normals = _make_lazy_cuda_func("depth_to_normal_forward")(
            camera_model_type, Ks,
            dist_coeffs,
            is_ray_depth, depths
        )
        ctx.save_for_backward(depths, Ks, dist_coeffs)
        ctx.camera_model_type = camera_model_type
        ctx.is_ray_depth = is_ray_depth

        return normals

    @staticmethod
    def backward(ctx, v_normals):
        depths, Ks, dist_coeffs = ctx.saved_tensors
        v_depths = _make_lazy_cuda_func("depth_to_normal_backward")(
            ctx.camera_model_type, Ks,
            dist_coeffs,
            ctx.is_ray_depth, depths, v_normals
        )
        return (v_depths, *([None]*(len(ctx.needs_input_grad)-1)))



def ray_depth_to_linear_depth(
    depths: Tensor,  # [B, H, W, 1]
    camera_model: Literal["pinhole", "fisheye"],
    Ks: Union[Tensor, Tuple[float, float, float, float]],
    dist_coeffs: Optional[Tensor] = None,  # [..., C, 10]
) -> Tensor:
    if isinstance(Ks, tuple):
        fx, fy, cx, cy = Ks
        Ks = torch.tensor([[[fx, 0, cx], [0, fy, cy], [0, 0, 1]]])
        Ks = Ks.to(depths).repeat(len(depths), 1, 1)
    return _RayDepthToLinearDepth.apply(
        depths,
        camera_model,
        Ks.contiguous(),
        dist_coeffs.contiguous() if dist_coeffs is not None else None,
    )


class _RayDepthToLinearDepth(torch.autograd.Function):
    """Projects Gaussians to 2D."""

    @staticmethod
    def forward(
        ctx,
        depths: Tensor,  # [B, H, W, 1]
        camera_model: Literal["pinhole", "fisheye"],
        Ks: Tensor,
        dist_coeffs: Optional[Tensor],  # [..., C, 10]
    ) -> Tensor:

        camera_model_type = gsplat.cuda._wrapper._make_lazy_cuda_obj(
            f"CameraModelType.{camera_model.upper()}"
        )

        out_depths = _make_lazy_cuda_func("ray_depth_to_linear_depth_forward")(
            camera_model_type, Ks, dist_coeffs,
            depths
        )

        ctx.save_for_backward(Ks, dist_coeffs)
        ctx.camera_model_type = camera_model_type
        
        return out_depths

    @staticmethod
    def backward(ctx, v_out_depths):

        Ks, dist_coeffs = ctx.saved_tensors
    
        v_in_depths = _make_lazy_cuda_func("ray_depth_to_linear_depth_backward")(
            ctx.camera_model_type, Ks, dist_coeffs,
            v_out_depths
        )

        return (v_in_depths, *([None]*(len(ctx.needs_input_grad)-1)))


def distort_image(
    image: Tensor,  # [B, H, W, C]
    camera_model: Literal["pinhole", "fisheye"],
    Ks: Union[Tensor, Tuple[float, float, float, float]],
    dist_coeffs: Optional[Tensor] = None,  # [..., C, 10]
) -> Tensor:
    return _DistortOrUndistortImage.apply(
        False, image, camera_model, Ks, dist_coeffs
    )

def undistort_image(
    image: Tensor,  # [B, H, W, C]
    camera_model: Literal["pinhole", "fisheye"],
    Ks: Union[Tensor, Tuple[float, float, float, float]],
    dist_coeffs: Optional[Tensor] = None,  # [..., C, 10]
) -> Tensor:
    return _DistortOrUndistortImage.apply(
        True, image, camera_model, Ks, dist_coeffs
    )


class _DistortOrUndistortImage(torch.autograd.Function):
    """Projects Gaussians to 2D."""

    @staticmethod
    def forward(
        ctx,
        is_undistort: bool,
        image: Tensor,  # [B, H, W, C]
        camera_model: Literal["pinhole", "fisheye"],
        Ks: Tensor,
        dist_coeffs: Optional[Tensor],  # [..., C, 10]
    ) -> Tensor:

        if isinstance(Ks, tuple):
            fx, fy, cx, cy = Ks
            Ks = torch.tensor([[[fx, 0, cx], [0, fy, cy], [0, 0, 1]]])
            Ks = Ks.to(image).repeat(len(image), 1, 1)
        if isinstance(dist_coeffs, tuple):
            dist_coeffs = torch.tensor(dist_coeffs)[None].to(image).repeat(len(image), 1)

        camera_model_type = gsplat.cuda._wrapper._make_lazy_cuda_obj(
            f"CameraModelType.{camera_model.upper()}"
        )

        return _make_lazy_cuda_func("un"*int(is_undistort) + "distort_image")(
            camera_model_type, Ks.contiguous(),
            dist_coeffs.contiguous() if dist_coeffs is not None else None,
            image.contiguous()
        )

    @staticmethod
    def backward(ctx, v_out_image):
        raise NotImplementedError("undistort_image does not support backward")
