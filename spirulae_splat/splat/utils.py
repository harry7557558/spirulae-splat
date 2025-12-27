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
    if d == 1:
        return image
    # weight = (1.0 / (d * d)) * torch.ones((1, 1, d, d), dtype=torch.float32, device=image.device)
    # return F.conv2d(image.float().permute(2, 0, 1)[:, None, ...], weight, stride=d).squeeze(1).permute(1, 2, 0).to(image)
    B, H, W, C = image.shape
    reshaped = image[:, :H//d*d, :W//d*d, :].view(B, H//d, d, W//d, d, C)
    blocks = reshaped.permute(0, 1, 3, 2, 4, 5).contiguous().view(B, H//d, W//d, d*d, C)
    return blocks.float().mean(dim=-2)


def rot_to_quat(R):
    """
    Convert rotation matrices to quaternions.
    R: (B,3,3)
    Returns q: (B,4) as (w,x,y,z)
    """
    B = R.shape[0]
    q = torch.zeros(B, 4, device=R.device, dtype=R.dtype)

    trace = R[:,0,0] + R[:,1,1] + R[:,2,2]
    pos = trace > 0

    # Case 1: trace positive → stable branch
    s = torch.sqrt(trace[pos] + 1.0) * 2
    q[pos, 0] = 0.25 * s
    q[pos, 1] = (R[pos, 2, 1] - R[pos, 1, 2]) / s
    q[pos, 2] = (R[pos, 0, 2] - R[pos, 2, 0]) / s
    q[pos, 3] = (R[pos, 1, 0] - R[pos, 0, 1]) / s

    # Case 2: diagonal branch
    neg = ~pos
    if neg.any():
        Rn = R[neg]
        diag = torch.stack([Rn[:,0,0],Rn[:,1,1],Rn[:,2,2]], dim=1)
        idx = diag.argmax(1)

        for k in range(3):
            mask = (idx == k)
            if not mask.any(): continue
            i,j,l = k, (k+1)%3, (k+2)%3
            Rm = Rn[mask]
            t = 1.0 + Rm[:,i,i] - Rm[:,j,j] - Rm[:,l,l]
            s = torch.sqrt(t + 1.0) * 2
            qc = q[neg][mask]
            qc[:,0] = (Rm[:,l,j] - Rm[:,j,l]) / s
            qc[:,1+i] = 0.25 * s
            qc[:,1+j] = (Rm[:,j,i] + Rm[:,i,j]) / s
            qc[:,1+l] = (Rm[:,l,i] + Rm[:,i,l]) / s
            q[neg][mask] = qc

    return torch.nn.functional.normalize(q, dim=1)


def quat_to_rot(q):
    """Quaternion → rotation matrix. q: (B,4)"""
    q = torch.nn.functional.normalize(q, dim=1)
    w, x, y, z = q[:,0], q[:,1], q[:,2], q[:,3]

    R = torch.stack([
        torch.stack([1-2*(y*y+z*z),   2*(x*y - z*w),   2*(x*z + y*w)], dim=1),
        torch.stack([2*(x*y + z*w),   1-2*(x*x+z*z),   2*(y*z - x*w)], dim=1),
        torch.stack([2*(x*z - y*w),   2*(y*z + x*w),   1-2*(x*x+y*y)], dim=1)
    ], dim=1)
    return R


def slerp(q0, q1, t):
    """
    Quaternion SLERP.
    q0, q1: (B,4)
    t: (B,)
    """
    q0 = torch.nn.functional.normalize(q0, dim=1)
    q1 = torch.nn.functional.normalize(q1, dim=1)

    dot = (q0 * q1).sum(1)
    flip = dot < 0
    q1 = q1.clone()
    q1[flip] *= -1
    dot = (q0 * q1).sum(1).clamp(-1, 1)

    theta = torch.acos(dot)
    sin_theta = torch.sin(theta)

    eps = 1e-6
    mask = sin_theta < eps

    # SLERP proper
    w1 = torch.sin((1 - t) * theta) / sin_theta
    w2 = torch.sin(t * theta) / sin_theta
    out = (q0 * w1.unsqueeze(1) + q1 * w2.unsqueeze(1))

    # LERP fallback (for near-parallel)
    out[mask] = q0[mask] * (1 - t[mask]).unsqueeze(1) + q1[mask] * t[mask].unsqueeze(1)
    return torch.nn.functional.normalize(out, dim=1)


def interpolate_se3(T0, T1, t):
    """
    Interpolate two SE3 transforms.
    T0,T1: (B,3,4) or (B,4,4)
    t: (B,) in [0,1]
    Returns: (B,4,4) valid SE3
    """
    if T0.shape[-2:] == (3,4):
        R0, t0 = T0[:,:3,:3], T0[:,:,3]
        R1, t1 = T1[:,:3,:3], T1[:,:,3]
    else:
        R0, t0 = T0[:,:3,:3], T0[:,:3,3]
        R1, t1 = T1[:,:3,:3], T1[:,:3,3]

    # Rotation: convert to quat → SLERP → convert back
    q0 = rot_to_quat(R0)
    q1 = rot_to_quat(R1)
    q = slerp(q0, q1, t)
    R = quat_to_rot(q)

    # Translation: simple lerp
    trans = (1 - t).unsqueeze(1) * t0 + t.unsqueeze(1) * t1

    # Assemble SE3
    B = t.shape[0]
    T = torch.zeros(B, 4, 4, device=T0.device, dtype=T0.dtype)
    T[:, :3, :3] = R
    T[:, :3, 3] = trans
    T[:, 3, 3] = 1.0
    return T


def random_c2w_on_unit_sphere(B, device="cpu", dtype=torch.float32):
    """
    Generate B random SE3 camera-to-world matrices.
    - camera center on unit sphere
    - z > 0 hemisphere
    - camera looks at origin (forward dir = -c)
    Returns: (B,4,4)
    """
    # Random points on S2 with positive z
    c = torch.randn(B, 3, device=device, dtype=dtype)
    # c[:,2].abs_()  # make z positive
    c = c / c.norm(dim=-1, keepdim=True)

    # Forward direction (camera looks at origin)
    # Convention: camera forward is +z axis in camera coords → use f = +c
    f = +c
    f = f / f.norm(dim=-1, keepdim=True)

    # Random world-up that isn't parallel to f
    up = torch.randn(B, 3, device=device, dtype=dtype)
    # Remove projection onto f (Gram-Schmidt)
    up = up - (up * f).sum(-1, keepdim=True) * f
    up = up / up.norm(dim=-1, keepdim=True)

    # Right = up × forward
    r = torch.cross(up, f, dim=1)
    r = r / r.norm(dim=-1, keepdim=True)

    # Recompute true orthonormal up = f × r
    u = torch.cross(f, r, dim=1)

    # Assemble rotation (R = [r u f])
    R = torch.stack([r, u, f], dim=2)   # (B,3,3)

    # Assemble final camera-to-world SE3
    T = torch.zeros(B, 4, 4, device=device, dtype=dtype)
    T[:, :3, :3] = R
    T[:, :3, 3] = c
    T[:, 3, 3] = 1.0
    return T


def ls_camera_intersection(T):
    """
    Compute the least-squares intersection point of camera view axes.
    
    T: (B,4,4) camera-to-world matrices
       - camera center = T[:, :3, 3]
       - camera view axis = T[:, :3, 2]   (world-space +z axis)
    
    Returns: (3,) tensor for the best point
    """
    B = T.shape[0]
    c = T[:, :3, 3]         # (B,3)
    d = T[:, :3, 2]         # (B,3), assumed unit or nearly unit
    d = d / d.norm(dim=1, keepdim=True)

    I = torch.eye(3, device=T.device, dtype=T.dtype)
    
    # A_i = I - d_i d_i^T
    outer = d.unsqueeze(2) @ d.unsqueeze(1)      # (B,3,3)
    A = I - outer                                 # (B,3,3)

    # Solve: (sum A_i) p = (sum A_i c_i)
    M = A.sum(dim=0)                              # (3,3)
    b = (A @ c.unsqueeze(2)).sum(dim=0).squeeze() # (3,)

    # Solve M p = b
    p = torch.linalg.solve(M, b)
    return p
