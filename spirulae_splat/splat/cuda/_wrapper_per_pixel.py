import torch
from torch import Tensor

from typing import Any, Callable, Literal, Optional, Union, Tuple
import functools

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


def linear_rgb_to_srgb(rgb: Tensor) -> Tensor:
    return _LinearRgbToSrgb.apply(rgb)

class _LinearRgbToSrgb(torch.autograd.Function):

    @staticmethod
    def forward(
        ctx,
        rgb
    ):
        out_rgb = _make_lazy_cuda_func("linear_rgb_to_srgb_forward")(rgb)

        ctx.save_for_backward(rgb)

        return out_rgb

    @staticmethod
    def backward(ctx, v_out_rgb):

        (rgb,) = ctx.saved_tensors

        return _make_lazy_cuda_func("linear_rgb_to_srgb_backward")(
            rgb, v_out_rgb
        ), None


@functools.cache
def get_color_transform_matrix(in_color_space: str, out_color_space: str = "Rec.709", device="cuda"):
    if out_color_space != "Rec.709":
        raise ValueError(f"Unsupported output color space {out_color_space}")
    # https://www.colour-science.org:8010/apps/rgb_colourspace_transformation_matrix
    if in_color_space == "ACES2065-1":
        return torch.tensor([
            [2.5247180476, -1.1325619434, -0.3921561044],
            [-0.2776344819, 1.3709123773, -0.0932778953],
            [-0.0165202369, -0.1479259606, 1.1644461975],
        ], device=device)
    if in_color_space == "ACEScg":
        return torch.tensor([
            [1.7072552160, -0.6200352595, -0.0872199564],
            [-0.1311566587, 1.1391010566, -0.0079443978],
            [-0.0245499075, -0.1248045805, 1.1493544880],
        ], device=device)
    if in_color_space == "Rec.2020":
        return torch.tensor([
            [1.6604910021, -0.5876411388, -0.0728498633],
            [-0.1245504745, 1.1328998971, -0.0083494226],
            [-0.0181507634, -0.1005788980, 1.1187296614],
        ], device=device)
    if in_color_space == "AdobeRGB":
        return torch.tensor([
            [1.3983671735, -0.3983451225, 0.0000054016],
            [-0.0000103176, 0.9999916496, -0.0000039459],
            [-0.0000003709, -0.0429269510, 1.0429319656],
        ], device=device)
    if in_color_space == "DCI-P3":
        return torch.tensor([
            [1.1548337042, -0.1451763523, -0.0096573518],
            [-0.0393300117, 1.0378282998, 0.0015017119],
            [-0.0184786235, -0.0689101110, 1.0873887345],
        ], device=device)
    raise ValueError(f"Unsupported input color space {in_color_space}")


def depth_to_normal(
    depths: Tensor,  # [B, H, W, 1]
    camera_model: Literal["pinhole", "fisheye"],
    intrins: Union[Tensor, Tuple[float, float, float, float]],
    dist_coeffs: Optional[Tensor] = None,  # [..., C, 10]
    is_ray_depth: bool = True,
    **kwargs
) -> Tensor:
    if isinstance(intrins, tuple):
        fx, fy, cx, cy = intrins
        intrins = torch.tensor([[fx, fy, cx, cy]])
        intrins = intrins.to(depths).repeat(len(depths), 1)
    return _DepthToNormal.apply(
        depths,
        camera_model,
        intrins.contiguous(),
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
        intrins: Tensor,
        dist_coeffs: Optional[Tensor],  # [..., C, 10]
        is_ray_depth: bool
    ) -> Tensor:

        camera_model_type = gsplat.cuda._wrapper._make_lazy_cuda_obj(
            f"CameraModelType.{camera_model.upper()}"
        )

        normals = _make_lazy_cuda_func("depth_to_normal_forward")(
            camera_model_type, intrins,
            dist_coeffs,
            is_ray_depth, depths
        )
        ctx.save_for_backward(depths, intrins, dist_coeffs)
        ctx.camera_model_type = camera_model_type
        ctx.is_ray_depth = is_ray_depth

        return normals

    @staticmethod
    def backward(ctx, v_normals):
        depths, intrins, dist_coeffs = ctx.saved_tensors
        v_depths = _make_lazy_cuda_func("depth_to_normal_backward")(
            ctx.camera_model_type, intrins,
            dist_coeffs,
            ctx.is_ray_depth, depths, v_normals
        )
        return (v_depths, *([None]*(len(ctx.needs_input_grad)-1)))



def ray_depth_to_linear_depth(
    depths: Tensor,  # [B, H, W, 1]
    camera_model: Literal["pinhole", "fisheye"],
    intrins: Union[Tensor, Tuple[float, float, float, float]],
    dist_coeffs: Optional[Tensor] = None,  # [..., C, 10]
) -> Tensor:
    if isinstance(intrins, tuple):
        fx, fy, cx, cy = intrins
        intrins = torch.tensor([[fx, fy, cx, cy]])
        intrins = intrins.to(depths).repeat(len(depths), 1)
    return _RayDepthToLinearDepth.apply(
        depths,
        camera_model,
        intrins.contiguous(),
        dist_coeffs.contiguous() if dist_coeffs is not None else None,
    )


class _RayDepthToLinearDepth(torch.autograd.Function):
    """Projects Gaussians to 2D."""

    @staticmethod
    def forward(
        ctx,
        depths: Tensor,  # [B, H, W, 1]
        camera_model: Literal["pinhole", "fisheye"],
        intrins: Tensor,
        dist_coeffs: Optional[Tensor],  # [..., C, 10]
    ) -> Tensor:

        camera_model_type = gsplat.cuda._wrapper._make_lazy_cuda_obj(
            f"CameraModelType.{camera_model.upper()}"
        )

        out_depths = _make_lazy_cuda_func("ray_depth_to_linear_depth_forward")(
            camera_model_type, intrins, dist_coeffs,
            depths
        )

        ctx.save_for_backward(intrins, dist_coeffs)
        ctx.camera_model_type = camera_model_type
        
        return out_depths

    @staticmethod
    def backward(ctx, v_out_depths):

        intrins, dist_coeffs = ctx.saved_tensors
    
        v_in_depths = _make_lazy_cuda_func("ray_depth_to_linear_depth_backward")(
            ctx.camera_model_type, intrins, dist_coeffs,
            v_out_depths
        )

        return (v_in_depths, *([None]*(len(ctx.needs_input_grad)-1)))


def distort_image(
    image: Tensor,  # [B, H, W, C]
    camera_model: Literal["pinhole", "fisheye"],
    intrins: Union[Tensor, Tuple[float, float, float, float]],
    dist_coeffs: Optional[Tensor] = None,  # [..., C, 10]
) -> Tensor:
    return _DistortOrUndistortImage.apply(
        False, image, camera_model, intrins, dist_coeffs
    )

def undistort_image(
    image: Tensor,  # [B, H, W, C]
    camera_model: Literal["pinhole", "fisheye"],
    intrins: Union[Tensor, Tuple[float, float, float, float]],
    dist_coeffs: Optional[Tensor] = None,  # [..., C, 10]
) -> Tensor:
    return _DistortOrUndistortImage.apply(
        True, image, camera_model, intrins, dist_coeffs
    )


class _DistortOrUndistortImage(torch.autograd.Function):
    """Projects Gaussians to 2D."""

    @staticmethod
    def forward(
        ctx,
        is_undistort: bool,
        image: Tensor,  # [B, H, W, C]
        camera_model: Literal["pinhole", "fisheye"],
        intrins: Tensor,
        dist_coeffs: Optional[Tensor],  # [..., C, 10]
    ) -> Tensor:

        if isinstance(intrins, tuple):
            fx, fy, cx, cy = intrins
            intrins = torch.tensor([[fx, fy, cx, cy]])
            intrins = intrins.to(image).repeat(len(image), 1)
        if isinstance(dist_coeffs, tuple):
            dist_coeffs = torch.tensor(dist_coeffs)[None].to(image).repeat(len(image), 1)

        camera_model_type = gsplat.cuda._wrapper._make_lazy_cuda_obj(
            f"CameraModelType.{camera_model.upper()}"
        )

        return _make_lazy_cuda_func("un"*int(is_undistort) + "distort_image")(
            camera_model_type, intrins.contiguous(),
            dist_coeffs.contiguous() if dist_coeffs is not None else None,
            image.contiguous()
        )

    @staticmethod
    def backward(ctx, v_out_image):
        raise NotImplementedError("undistort_image does not support backward")


def apply_ppisp(
    image: Tensor,  # [B, H, W, 3]
    ppisp_params: Tensor,  # [B, P]
    intrins: Union[Tensor, Tuple[float, float, float, float]], # [B, 4]
    actual_image_width: Optional[int] = None,
    actual_image_height: Optional[int] = None
) -> Tensor:
    if isinstance(intrins, tuple):
        fx, fy, cx, cy = intrins
        intrins = torch.tensor([[fx, fy, cx, cy]])
        intrins = intrins.to(image).repeat(len(image), 1)
    if actual_image_width is None:
        actual_image_width = image.shape[-2]
    if actual_image_height is None:
        actual_image_height = image.shape[-3]
    return _PPISP.apply(
        image,
        ppisp_params.contiguous(),
        intrins.contiguous(),
        actual_image_width,
        actual_image_height
    )

class _PPISP(torch.autograd.Function):
    """Applies Physically-Plausible Compensation and Control of Photometric Variations (PPISP) to images."""

    @staticmethod
    def forward(
        ctx,
        image: Tensor,  # [B, H, W, 3]
        ppisp_params: Tensor,  # [B, P]
        intrins: Tensor, # [B, 4]
        actual_image_width: int,
        actual_image_height: int
    ) -> Tensor:

        out_image = _make_lazy_cuda_func("ppisp_forward")(
            image, ppisp_params,
            intrins, actual_image_width, actual_image_height
        )

        ctx.save_for_backward(image, ppisp_params, intrins)
        ctx.actual_image_width = actual_image_width
        ctx.actual_image_height = actual_image_height

        return out_image

    @staticmethod
    def backward(ctx, v_out_image):

        image, ppisp_params, intrins = ctx.saved_tensors

        v_image, v_ppisp_params = _make_lazy_cuda_func("ppisp_backward")(
            image, ppisp_params,
            intrins, ctx.actual_image_width, ctx.actual_image_height,
            v_out_image
        )

        return v_image, v_ppisp_params, None, None, None
