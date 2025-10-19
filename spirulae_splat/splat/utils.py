"""Python bindings for binning and sorting gaussians"""

from typing import Tuple

import torch
import numpy as np
from typing import Optional
from jaxtyping import Float, Int
from torch import Tensor

import spirulae_splat.splat.cuda as _C
from spirulae_splat.splat._camera import _Camera


_TORCH_COMPILE_ARGS = {
    'dynamic': True,
    # 'force_parameter_static_shapes': False
}
torch._dynamo.config.force_nn_module_property_static_shapes = False
torch._dynamo.config.force_parameter_static_shapes = False
torch._dynamo.config.capture_scalar_outputs = True


def compute_cumulative_intersects(
    num_tiles_hit: Float[Tensor, "batch 1"]
) -> Tuple[int, Float[Tensor, "batch 1"]]:
    """Computes cumulative intersections of gaussians. This is useful for creating unique gaussian IDs and for sorting.

    Note:
        This function is not differentiable to any input.

    Args:
        num_tiles_hit (Tensor): number of intersected tiles per gaussian.

    Returns:
        A tuple of {int, Tensor}:

        - **num_intersects** (int): total number of tile intersections.
        - **cum_tiles_hit** (Tensor): a tensor of cumulated intersections (used for sorting).
    """
    cum_tiles_hit = torch.cumsum(num_tiles_hit, dim=0, dtype=torch.int32)
    num_intersects = cum_tiles_hit[-1].item()
    return num_intersects, cum_tiles_hit


def bin_and_sort_gaussians(
    num_points: int,
    num_intersects: int,
    positions: Float[Tensor, "batch 3"],
    bounds: Int[Tensor, "batch 4"],
    cum_tiles_hit: Float[Tensor, "batch 1"],
    img_height: int,
    img_width: int,
) -> Tuple[
    Float[Tensor, "num_intersects 1"],
    Float[Tensor, "num_intersects 1"],
    Float[Tensor, "num_intersects 1"],
    Float[Tensor, "num_intersects 1"],
    Float[Tensor, "num_intersects 2"],
]:
    """Mapping gaussians to sorted unique intersection IDs and tile bins used for fast rasterization.

    We return both sorted and unsorted versions of intersect IDs and gaussian IDs for testing purposes.

    Note:
        This function is not differentiable to any input.

    Args:
        num_points (int): number of gaussians.
        num_intersects (int): cumulative number of total gaussian intersections
        xys (Tensor): x,y locations of 2D gaussian projections.
        depths (Tensor): z depth of gaussians.
        radii (Tensor): radii of 2D gaussian projections.
        cum_tiles_hit (Tensor): list of cumulative tiles hit.
        tile_bounds (Tuple): tile dimensions as a len 3 tuple (tiles.x , tiles.y, 1).

    Returns:
        A tuple of {Tensor, Tensor, Tensor, Tensor, Tensor}:

        - **isect_ids_unsorted** (Tensor): unique IDs for each gaussian in the form (tile | depth id).
        - **gaussian_ids_unsorted** (Tensor): Tensor that maps isect_ids back to cum_tiles_hit. Useful for identifying gaussians.
        - **isect_ids_sorted** (Tensor): sorted unique IDs for each gaussian in the form (tile | depth id).
        - **gaussian_ids_sorted** (Tensor): sorted Tensor that maps isect_ids back to cum_tiles_hit. Useful for identifying gaussians.
        - **tile_bins** (Tensor): range of gaussians hit per tile.
    """
    isect_ids, gaussian_ids = _C.map_gaussian_to_intersects(
        num_points,
        num_intersects,
        positions.contiguous(),
        bounds.contiguous(),
        cum_tiles_hit.contiguous(),
        img_height, img_width
    )
    isect_ids_sorted, sorted_indices = torch.sort(isect_ids)
    gaussian_ids_sorted = torch.gather(gaussian_ids, 0, sorted_indices)
    tile_bins = _C.get_tile_bin_edges(
        num_intersects, isect_ids_sorted.contiguous(),
        img_height, img_width
    )
    return isect_ids, gaussian_ids, isect_ids_sorted, gaussian_ids_sorted, tile_bins


# @torch.compile(**_TORCH_COMPILE_ARGS)
def depth_to_points(
    depths: Tensor, camera: _Camera, c2w: Optional[Tensor]=None, z_depth: bool = True
) -> Tensor:
    """Convert depth maps to 3D points

    Args:
        depths: Depth maps [..., H, W, 1]
        camera: Camera intrinsics
        c2w: Camera-to-world transformation matrices [..., 4, 4]
        z_depth: Whether the depth is in z-depth (True) or ray depth (False)

    Returns:
        points: 3D points in the world coordinate system [..., H, W, 3]
    """
    assert depths.shape[-1] == 1, f"Invalid depth shape: {depths.shape}"
    if c2w is not None:
        assert c2w.shape[-2:] == (4, 4), f"Invalid viewmats shape: {c2w.shape}"

    # camera directions in camera coordinates
    undist_map = camera.get_undist_map(always=True)
    undist_map = undist_map.reshape(*([1]*(len(depths.shape)-3)), *depths.shape[-3:-1], 2)
    camera_dirs = torch.concatenate([
        undist_map, torch.ones_like(undist_map[..., :1])
    ], dim=-1)  # [..., H, W, 3]

    # ray directions in world coordinates
    if c2w is not None:
        directions = torch.einsum(
            "...ij,...hwj->...hwi", c2w[..., :3, :3], camera_dirs
        )  # [..., H, W, 3]
        origins = c2w[..., :3, -1]  # [..., 3]
    else:
        directions = camera_dirs

    if not z_depth:
        directions = torch.nn.functional.normalize(directions, dim=-1)

    if c2w is None:
        return depths * directions
    return origins[..., None, None, :] + depths * directions


# @torch.compile(**_TORCH_COMPILE_ARGS)
def depth_to_normal(
    depths: Tensor, camera: _Camera, c2w: Optional[Tensor]=None, z_depth: bool = True, alpha: Optional[Tensor] = None,
    return_points = False
) -> Tensor:
    """Convert depth maps to surface normals

    Args:
        depths: Depth maps [H, W, 1]
        camera: Camera intrinsics
        c2w: Camera-to-world transformation matrices [4, 4]
        z_depth: Whether the depth is in z-depth (True) or ray depth (False)

    Returns:
        normals: Surface normals in the world coordinate system [H, W, 3]
    """
    # TODO: write a fully fused CUDA kernel for this

    points = depth_to_points(depths, camera, c2w, z_depth=z_depth)  # [H, W, 3]
    dx = torch.cat(
        [points[2:, 1:-1, :] - points[:-2, 1:-1, :]], dim=-3
    )  # [H-2, W-2, 3]
    dy = torch.cat(
        [points[1:-1, 2:, :] - points[1:-1, :-2, :]], dim=-2
    )  # [H-2, W-2, 3]
    normals = torch.nn.functional.normalize(torch.cross(dx, dy, dim=-1), dim=-1)  # [H-2, W-2, 3]
    normals = torch.nn.functional.pad(normals, (0, 0, 1, 1, 1, 1), value=0.0)  # [H, W, 3]
    # normals = torch.nn.functional.pad(normals.permute(2, 0, 1), (1, 1, 1, 1), mode="replicate").permute(1, 2, 0)  # [H, W, 3]
    normals = normals.contiguous()

    # apply mask
    if alpha is not None:
        return normals * (alpha>0).float().reshape((*normals.shape[:2], 1)), alpha
        kernel = torch.tensor([[[
            [0, 1, 0], [1, 1, 1], [0, 1, 0]
        ]]], dtype=normals.dtype, device=alpha.device)
        conv_input = (alpha>0).float().reshape((1, 1, *alpha.shape[:2]))
        conv_result = torch.nn.functional.conv2d(conv_input, kernel, padding=1)
        alpha = (conv_result == 5).reshape((*normals.shape[:2], 1))
        return normals * alpha, alpha

    if return_points:
        return normals, points
    return normals


# @torch.compile(**_TORCH_COMPILE_ARGS)
def resize_image(image: torch.Tensor, d: int):
    """
    Downscale images using the same 'area' method in opencv

    :param image shape [B, H, W, C]
    :param d downscale factor (must be 2, 4, 8, etc.)

    return downscaled image in shape [B, H//d, W//d, C]
    """
    # weight = (1.0 / (d * d)) * torch.ones((1, 1, d, d), dtype=torch.float32, device=image.device)
    # return F.conv2d(image.float().permute(2, 0, 1)[:, None, ...], weight, stride=d).squeeze(1).permute(1, 2, 0).to(image)
    B, H, W, C = image.shape
    reshaped = image[:, :H//d*d, :W//d*d, :].view(B, H//d, d, W//d, d, C)
    blocks = reshaped.permute(0, 1, 3, 2, 4, 5).contiguous().view(B, H//d, W//d, d*d, C)
    return blocks.float().mean(dim=-2)
