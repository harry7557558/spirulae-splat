"""Pure PyTorch implementations of various functions"""

import struct

import torch
import torch.nn.functional as F
from jaxtyping import Float, Int
from torch import Tensor
from typing import Tuple, Literal, Optional
from .cuda import BLOCK_WIDTH


class BesselJ0Function(torch.autograd.Function):
    @staticmethod
    def forward(ctx, x):
        ctx.save_for_backward(x)
        return torch.special.bessel_j0(x)
    @staticmethod
    def backward(ctx, grad_output):
        x, = ctx.saved_tensors
        # Derivative of J0 is -J1
        grad_x = -torch.special.bessel_j1(x) * grad_output
        return grad_x

class BesselJ1Function(torch.autograd.Function):
    @staticmethod
    def forward(ctx, x):
        ctx.save_for_backward(x)
        return torch.special.bessel_j1(x)
    @staticmethod
    def backward(ctx, grad_output):
        x, = ctx.saved_tensors
        # Derivative of J1 is (J0 - J1/x)
        grad_x = (torch.special.bessel_j0(x) - torch.special.bessel_j1(x) / x) * grad_output
        grad_x[x == 0] = 0  # Handle the division by zero
        return grad_x

def bessel_j0(x):
    return BesselJ0Function.apply(x)

def bessel_j1(x):
    return BesselJ1Function.apply(x)

def bessel_j(m, x):
    if m == 0:
        return bessel_j0(x)
    if m == 1:
        return bessel_j1(x)
    return 2*(m-1)/x * bessel_j(m-1, x) - bessel_j(m-2, x)



def quat_mult(q1, q2):
    w1, x1, y1, z1 = torch.unbind(q1, dim=-1)
    w2, x2, y2, z2 = torch.unbind(q2, dim=-1)
    w = w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2
    x = w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2
    y = w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2
    z = w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2
    return torch.stack([w, x, y, z]).T.contiguous()


def quat_to_rotmat(quat: Tensor) -> Tensor:
    assert quat.shape[-1] == 4, quat.shape
    w, x, y, z = torch.unbind(quat, dim=-1)
    mat = torch.stack(
        [
            1 - 2 * (y**2 + z**2),
            2 * (x * y - w * z),
            2 * (x * z + w * y),
            2 * (x * y + w * z),
            1 - 2 * (x**2 + z**2),
            2 * (y * z - w * x),
            2 * (x * z - w * y),
            2 * (y * z + w * x),
            1 - 2 * (x**2 + y**2),
        ],
        dim=-1,
    )
    return mat.reshape(quat.shape[:-1] + (3, 3))


def scale_rot_to_cov3d(scale: Tensor, quat: Tensor) -> Tensor:
    assert scale.shape[-1] == 2, scale.shape
    assert quat.shape[-1] == 4, quat.shape
    assert scale.shape[:-1] == quat.shape[:-1], (scale.shape, quat.shape)
    scale = torch.concat((scale, torch.zeros_like(scale[:,0:1])), axis=1)
    R = quat_to_rotmat(quat)  # (..., 3, 3)
    M = R * scale[..., None, :]  # (..., 3, 3)
    # TODO: save upper right because symmetric
    return M @ M.transpose(-1, -2)  # (..., 3, 3)


def quat_scale_to_covar_preci(
    quats: Tensor,  # [N, 4],
    scales: Tensor,  # [N, 3],
    compute_covar: bool = True,
    compute_preci: bool = True,
    triu: bool = False,
) -> Tuple[Optional[Tensor], Optional[Tensor]]:
    """PyTorch implementation of `gsplat.cuda._wrapper.quat_scale_to_covar_preci()`."""
    R = quat_to_rotmat(quats)  # (..., 3, 3)

    if compute_covar:
        M = R * scales[..., None, :]  # (..., 3, 3)
        covars = torch.bmm(M, M.transpose(-1, -2))  # (..., 3, 3)
        if triu:
            covars = covars.reshape(covars.shape[:-2] + (9,))  # (..., 9)
            covars = (
                covars[..., [0, 1, 2, 4, 5, 8]] + covars[..., [0, 3, 6, 4, 7, 8]]
            ) / 2.0  # (..., 6)
    if compute_preci:
        P = R * (1 / scales[..., None, :])  # (..., 3, 3)
        precis = torch.bmm(P, P.transpose(-1, -2))  # (..., 3, 3)
        if triu:
            precis = precis.reshape(precis.shape[:-2] + (9,))
            precis = (
                precis[..., [0, 1, 2, 4, 5, 8]] + precis[..., [0, 3, 6, 4, 7, 8]]
            ) / 2.0

    return covars if compute_covar else None, precis if compute_preci else None


def compute_compensation(cov2d_mat: Tensor):
    """
    params: cov2d matrix (*, 2, 2)
    returns: compensation factor as calculated in project_cov3d_ewa
    """
    det_denom = cov2d_mat[..., 0, 0] * cov2d_mat[..., 1, 1] - cov2d_mat[..., 0, 1] ** 2
    det_nomin = (cov2d_mat[..., 0, 0] - 0.3) * (cov2d_mat[..., 1, 1] - 0.3) - cov2d_mat[
        ..., 0, 1
    ] ** 2
    return torch.sqrt(torch.clamp(det_nomin / det_denom, min=0))


def project_ellipse_bound(T, V0, V1, fx, fy, cx, cy):

    V01 = torch.cross(V0, V1, dim=-1)
    V0T = torch.cross(T, V0, dim=-1)
    V1T = torch.cross(T, V1, dim=-1)

    A = V0T[:,0]**2 + V1T[:,0]**2 - V01[:,0]**2
    B = -V01[:,1]*V01[:,0] + V1T[:,1]*V1T[:,0] + V0T[:,1]*V0T[:,0]
    C = V0T[:,1]**2 + V1T[:,1]**2 - V01[:,1]**2
    D = 2.0 * (V0T[:,2]*V0T[:,0] + V1T[:,2]*V1T[:,0] - V01[:,2]*V01[:,0])
    E = 2.0 * (-V01[:,2]*V01[:,1] + V1T[:,2]*V1T[:,1] + V0T[:,2]*V0T[:,1])
    F = V0T[:,2]**2 + V1T[:,2]**2 - V01[:,2]**2

    valid = (B * B < A * C)

    # Translate to origin
    U = (C * D - B * E) / (2.0 * (B * B - A * C))
    V = (A * E - B * D) / (2.0 * (B * B - A * C))
    S = -(A * U**2 + 2.0 * B * U * V + C * V**2 + D * U + E * V + F)

    # Image transform
    U_T = fx * U + cx
    V_T = fy * V + cy
    A_T = A / (fx**2)
    B_T = B / (fx*fy)
    C_T = C / (fy**2)

    # Axis-aligned bounding box
    W_T = fx * torch.sqrt(C * S / (A * C - B * B))
    H_T = fy * torch.sqrt(A * S / (A * C - B * B))

    # Bounding circle
    L_T = 0.5 * (A_T + C_T - torch.sqrt((A_T - C_T)**2 + 4.0 * B_T**2))
    R_T = torch.sqrt(S / L_T)

    # Output
    center = torch.stack([U_T, V_T], dim=1)
    bound = torch.stack([W_T, H_T, R_T], dim=1)
    return center, bound, valid


def project_pix(fxfy, p_view, center, eps=1e-6):
    fx, fy = fxfy
    cx, cy = center

    rw = 1.0 / (p_view[..., 2] + 1e-6)
    p_proj = (p_view[..., 0] * rw, p_view[..., 1] * rw)
    u, v = (p_proj[0] * fx + cx, p_proj[1] * fy + cy)
    return torch.stack([u, v], dim=-1)


def clip_near_plane(p, viewmat, clip_thresh=0.01):
    R = viewmat[:3, :3]
    T = viewmat[:3, 3]
    p_view = torch.einsum("ij,nj->ni", R, p) + T[None]
    return p_view, p_view[..., 2] < clip_thresh


def get_tile_bbox(center, bound, tile_bounds):
    tile_size = torch.tensor(
        [BLOCK_WIDTH, BLOCK_WIDTH], dtype=torch.float32, device=center.device
    )
    tile_center = center / tile_size
    tile_radius = bound[..., :2] / tile_size

    top_left = (tile_center - tile_radius).to(torch.int32)
    bottom_right = (tile_center + tile_radius).to(torch.int32) + 1
    tile_min = torch.stack(
        [
            torch.clamp(top_left[..., 0], 0, tile_bounds[0]),
            torch.clamp(top_left[..., 1], 0, tile_bounds[1]),
        ],
        -1,
    )
    tile_max = torch.stack(
        [
            torch.clamp(bottom_right[..., 0], 0, tile_bounds[0]),
            torch.clamp(bottom_right[..., 1], 0, tile_bounds[1]),
        ],
        -1,
    )
    return tile_min, tile_max


def map_gaussian_to_intersects(
    num_points, xys, depths, radii, cum_tiles_hit, tile_bounds
):
    num_intersects = cum_tiles_hit[-1]
    isect_ids = torch.zeros(num_intersects, dtype=torch.int64, device=xys.device)
    gaussian_ids = torch.zeros(num_intersects, dtype=torch.int32, device=xys.device)

    for idx in range(num_points):
        if radii[idx] <= 0:
            break

        tile_min, tile_max = get_tile_bbox(xys[idx], radii[idx], tile_bounds)

        cur_idx = 0 if idx == 0 else cum_tiles_hit[idx - 1].item()

        # Get raw byte representation of the float value at the given index
        raw_bytes = struct.pack("f", depths[idx])

        # Interpret those bytes as an int32_t
        depth_id_n = struct.unpack("i", raw_bytes)[0]

        for i in range(tile_min[1], tile_max[1]):
            for j in range(tile_min[0], tile_max[0]):
                tile_id = i * tile_bounds[0] + j
                isect_ids[cur_idx] = (tile_id << 32) | depth_id_n
                gaussian_ids[cur_idx] = idx
                cur_idx += 1

    return isect_ids, gaussian_ids


def get_tile_bin_edges(num_intersects, isect_ids_sorted, tile_bounds):
    tile_bins = torch.zeros(
        (tile_bounds[0] * tile_bounds[1], 2),
        dtype=torch.int32,
        device=isect_ids_sorted.device,
    )

    for idx in range(num_intersects):
        cur_tile_idx = isect_ids_sorted[idx] >> 32

        if idx == 0:
            tile_bins[cur_tile_idx, 0] = 0
            continue

        if idx == num_intersects - 1:
            tile_bins[cur_tile_idx, 1] = num_intersects
            break

        prev_tile_idx = isect_ids_sorted[idx - 1] >> 32

        if cur_tile_idx != prev_tile_idx:
            tile_bins[prev_tile_idx, 1] = idx
            tile_bins[cur_tile_idx, 0] = idx

    return tile_bins


def visibility_kernel(r2):
    # return (1.0-r2) * (1.0-r2) * (r2 < 1.0)
    return (1.0-r2) * (r2 < 1.0)

def get_alpha(uv, opac) -> Tuple[Tensor, bool]:
    r2 = torch.norm(uv)**2
    vis = visibility_kernel(r2)
    t = 0.0 # torch.dot(uv, aniso)
    # t = torch.clamp(t, torch.zeros_like(t), torch.ones_like(t))
    m = t*t*(2.0*t-3.0)+1.0
    alpha = opac * vis * m
    return alpha, r2 >= 0.0 and alpha >= 1e-3


def get_intersection(position, axis_uv, pos_2d):
    pos_2d_3 = torch.concat((pos_2d, torch.ones_like(pos_2d[0:1])))
    A = torch.concat((*axis_uv, pos_2d_3)).reshape((3,3)).T
    uvt = -torch.linalg.inv(A) @ position
    uv = uvt[:2]
    t = -uvt[2:]
    return pos_2d_3*t, uv, torch.norm(uv) < 1.0


