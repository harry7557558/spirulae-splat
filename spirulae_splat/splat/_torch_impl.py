"""Pure PyTorch implementations of various functions"""

import struct

import torch
import torch.nn.functional as F
from jaxtyping import Float, Int
from torch import Tensor
from typing import Tuple, Literal, Optional
from spirulae_splat.splat.utils import (
    compute_cumulative_intersects, bin_and_sort_gaussians,
    _TORCH_COMPILE_ARGS
)
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


@torch.compile(**_TORCH_COMPILE_ARGS)
def depth_map(z):
    return torch.where(z>0.0, torch.log(z+1.0), z)

@torch.compile(**_TORCH_COMPILE_ARGS)
def depth_inv_map(z):
    return torch.where(z>0.0, torch.exp(z)-1.0, z)


def compute_sh_color(
    viewdirs: Float[Tensor, "*batch 3"],
    sh_coeffs: Float[Tensor, "*batch D C"],
    method: Literal["poly", "fast"] = "fast",
):
    """
    :param viewdirs (*, C)
    :param sh_coeffs (*, D, C) sh coefficients for each color channel
    return colors (*, C)
    """
    *dims, dim_sh, C = sh_coeffs.shape
    if method == "poly":
        bases = eval_sh_bases(dim_sh, viewdirs)  # (*, dim_sh)
    elif method == "fast":
        bases = eval_sh_bases_fast(dim_sh, viewdirs)  # (*, dim_sh)
    else:
        raise RuntimeError(f"Unknown mode: {method} for compute sh color.")
    return (bases[..., None] * sh_coeffs).sum(dim=-2)


"""
Taken from https://github.com/sxyu/svox2
"""

SH_C0 = 0.28209479177387814
SH_C1 = 0.4886025119029199
SH_C2 = [
    1.0925484305920792,
    -1.0925484305920792,
    0.31539156525252005,
    -1.0925484305920792,
    0.5462742152960396,
]
SH_C3 = [
    -0.5900435899266435,
    2.890611442640554,
    -0.4570457994644658,
    0.3731763325901154,
    -0.4570457994644658,
    1.445305721320277,
    -0.5900435899266435,
]
SH_C4 = [
    2.5033429417967046,
    -1.7701307697799304,
    0.9461746957575601,
    -0.6690465435572892,
    0.10578554691520431,
    -0.6690465435572892,
    0.47308734787878004,
    -1.7701307697799304,
    0.6258357354491761,
]

MAX_SH_BASIS = 10


def eval_sh_bases(basis_dim: int, dirs: torch.Tensor):
    """
    Evaluate spherical harmonics bases at unit directions,
    without taking linear combination.
    At each point, the final result may the be
    obtained through simple multiplication.

    :param basis_dim: int SH basis dim. Currently, 1-25 square numbers supported
    :param dirs: torch.Tensor (..., 3) unit directions

    :return: torch.Tensor (..., basis_dim)
    """
    result = torch.empty(
        (*dirs.shape[:-1], basis_dim), dtype=dirs.dtype, device=dirs.device
    )
    result[..., 0] = SH_C0
    if basis_dim > 1:
        x, y, z = dirs.unbind(-1)
        result[..., 1] = -SH_C1 * y
        result[..., 2] = SH_C1 * z
        result[..., 3] = -SH_C1 * x
        if basis_dim > 4:
            xx, yy, zz = x * x, y * y, z * z
            xy, yz, xz = x * y, y * z, x * z
            result[..., 4] = SH_C2[0] * xy
            result[..., 5] = SH_C2[1] * yz
            result[..., 6] = SH_C2[2] * (2.0 * zz - xx - yy)
            result[..., 7] = SH_C2[3] * xz
            result[..., 8] = SH_C2[4] * (xx - yy)

            if basis_dim > 9:
                result[..., 9] = SH_C3[0] * y * (3 * xx - yy)
                result[..., 10] = SH_C3[1] * xy * z
                result[..., 11] = SH_C3[2] * y * (4 * zz - xx - yy)
                result[..., 12] = SH_C3[3] * z * (2 * zz - 3 * xx - 3 * yy)
                result[..., 13] = SH_C3[4] * x * (4 * zz - xx - yy)
                result[..., 14] = SH_C3[5] * z * (xx - yy)
                result[..., 15] = SH_C3[6] * x * (xx - 3 * yy)

                if basis_dim > 16:
                    result[..., 16] = SH_C4[0] * xy * (xx - yy)
                    result[..., 17] = SH_C4[1] * yz * (3 * xx - yy)
                    result[..., 18] = SH_C4[2] * xy * (7 * zz - 1)
                    result[..., 19] = SH_C4[3] * yz * (7 * zz - 3)
                    result[..., 20] = SH_C4[4] * (zz * (35 * zz - 30) + 3)
                    result[..., 21] = SH_C4[5] * xz * (7 * zz - 3)
                    result[..., 22] = SH_C4[6] * (xx - yy) * (7 * zz - 1)
                    result[..., 23] = SH_C4[7] * xz * (xx - 3 * yy)
                    result[..., 24] = SH_C4[8] * (
                        xx * (xx - 3 * yy) - yy * (3 * xx - yy)
                    )
    return result


def eval_sh_bases_fast(basis_dim: int, dirs: torch.Tensor):
    """
    Evaluate spherical harmonics bases at unit direction for high orders
    using approach described by
    Efficient Spherical Harmonic Evaluation, Peter-Pike Sloan, JCGT 2013
    https://jcgt.org/published/0002/02/06/


    :param basis_dim: int SH basis dim. Currently, only 1-25 square numbers supported
    :param dirs: torch.Tensor (..., 3) unit directions

    :return: torch.Tensor (..., basis_dim)

    See reference C++ code in https://jcgt.org/published/0002/02/06/code.zip
    """
    result = torch.empty(
        (*dirs.shape[:-1], basis_dim), dtype=dirs.dtype, device=dirs.device
    )

    result[..., 0] = 0.2820947917738781

    if basis_dim <= 1:
        return result

    x, y, z = dirs.unbind(-1)

    fTmpA = -0.48860251190292
    result[..., 2] = 0.4886025119029199 * z
    result[..., 3] = fTmpA * x
    result[..., 1] = fTmpA * y

    if basis_dim <= 4:
        return result

    z2 = z * z
    fTmpB = -1.092548430592079 * z
    fTmpA = 0.5462742152960395
    fC1 = x * x - y * y
    fS1 = 2 * x * y
    result[..., 6] = 0.9461746957575601 * z2 - 0.3153915652525201
    result[..., 7] = fTmpB * x
    result[..., 5] = fTmpB * y
    result[..., 8] = fTmpA * fC1
    result[..., 4] = fTmpA * fS1

    if basis_dim <= 9:
        return result

    fTmpC = -2.285228997322329 * z2 + 0.4570457994644658
    fTmpB = 1.445305721320277 * z
    fTmpA = -0.5900435899266435
    fC2 = x * fC1 - y * fS1
    fS2 = x * fS1 + y * fC1
    result[..., 12] = z * (1.865881662950577 * z2 - 1.119528997770346)
    result[..., 13] = fTmpC * x
    result[..., 11] = fTmpC * y
    result[..., 14] = fTmpB * fC1
    result[..., 10] = fTmpB * fS1
    result[..., 15] = fTmpA * fC2
    result[..., 9] = fTmpA * fS2

    if basis_dim <= 16:
        return result

    fTmpD = z * (-4.683325804901025 * z2 + 2.007139630671868)
    fTmpC = 3.31161143515146 * z2 - 0.47308734787878
    fTmpB = -1.770130769779931 * z
    fTmpA = 0.6258357354491763
    fC3 = x * fC2 - y * fS2
    fS3 = x * fS2 + y * fC2
    result[..., 20] = (
        1.984313483298443 * z * result[..., 12] + -1.006230589874905 * result[..., 6]
    )
    result[..., 21] = fTmpD * x
    result[..., 19] = fTmpD * y
    result[..., 22] = fTmpC * fC1
    result[..., 18] = fTmpC * fS1
    result[..., 23] = fTmpB * fC2
    result[..., 17] = fTmpB * fS2
    result[..., 24] = fTmpA * fC3
    result[..., 16] = fTmpA * fS3
    return result


def ch_coeffs_to_color(
        degree_r, degree_r_to_use,
        degree_phi, degree_phi_to_use,
        coeffs, uv):
    r = torch.norm(uv, dim=-1)
    phi = torch.atan2(uv[...,1], uv[...,0])
    pi = torch.pi

    idx = 0
    color = torch.zeros(3, device=uv.device)
    for k in range(1, degree_r+1):
        if k > degree_r_to_use:
            idx += 2*degree_phi+1
            continue
        w = bessel_j(0, k*pi*r)
        color = color + w * coeffs[idx]
        idx += 1
        for m in range(1, degree_phi+1):
            if m > degree_phi_to_use:
                idx += 2
                continue
            wr = bessel_j(m, k*pi*r)
            wc = wr * torch.cos(m*phi)
            ws = wr * torch.sin(m*phi)
            color = color + wc * coeffs[idx]
            idx += 1
            color = color + ws * coeffs[idx]
            idx += 1

    assert idx == degree_r * (2*degree_phi+1)
    return color


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



def render_background_sh(
    w: int,
    h: int,
    intrins: Tuple[float, float, float, float],
    rotation: Float[Tensor, "3 3"],
    sh_degree: int,
    sh_coeffs: Float[Tensor, "sh_degree**2 3"],
) -> Float[Tensor, "h w 3"]:
    device = sh_coeffs.device
    fx, fy, cx, cy = intrins

    y, x = torch.meshgrid(torch.arange(h), torch.arange(w), indexing="ij")
    x = x.flatten().float().to(device) + 0.5
    y = y.flatten().float().to(device) + 0.5

    coord = torch.stack([(x - cx) / fx, -(y - cy) / fy, -torch.ones_like(x)], -1)
    if True:
        directions = torch.matmul(coord, rotation.T)
        norm = torch.linalg.norm(directions, dim=-1, keepdims=True)
        directions = directions / norm
    else:
        # no much magnitude difference, use above for denser gradient
        norm = torch.linalg.norm(coord, dim=-1, keepdims=True)
        directions = coord / norm @ rotation.T

    # TODO: deprecated in recent nerfstudio versions
    from nerfstudio.utils.spherical_harmonics import components_from_spherical_harmonics
    sh_components = components_from_spherical_harmonics(sh_degree, directions)  # [w*h, deg^2]
    bg_flat = torch.matmul(sh_components, sh_coeffs)  # [w*h, 3]
    bg_flat = torch.relu(bg_flat+0.5)

    return bg_flat.view(h, w, 3)

