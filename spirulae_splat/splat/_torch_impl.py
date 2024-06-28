"""Pure PyTorch implementations of various functions"""

import struct

import torch
import torch.nn.functional as F
from jaxtyping import Float, Int
from torch import Tensor
from typing import Tuple, Literal, Optional
from spirulae_splat.splat.utils import (
    compute_cumulative_intersects, bin_and_sort_gaussians
)


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


def get_tile_bbox(center, bound, tile_bounds, block_width):
    tile_size = torch.tensor(
        [block_width, block_width], dtype=torch.float32, device=center.device
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


def projected_depth_grad(p, R, fx, fy):
    n1 = R[:,:,2]
    invJ = torch.zeros((len(p), 3, 3), dtype=p.dtype, device=p.device)
    invJ[:,0,0] = p[:,2]/fx
    invJ[:,1,1] = p[:,2]/fy
    invJ[:,2,0] = p[:,0]/p[:,2]
    invJ[:,2,1] = p[:,1]/p[:,2]
    invJ[:,2,2] = 1.0
    n = torch.einsum('nij,nj->ni', invJ, n1)
    # return -n[:,:2] / n[:,2:]
    return -n[:,:2] * n[:,2:]


def project_gaussians(
    means3d, scales, quats,
    viewmat, intrins, img_size, block_width,
    clip_thresh=0.01,
):
    tile_bounds = (
        (img_size[0] + block_width - 1) // block_width,
        (img_size[1] + block_width - 1) // block_width,
        1,
    )
    fx, fy, cx, cy = intrins
    p_view, is_close = clip_near_plane(means3d, viewmat, clip_thresh)
    R1 = viewmat[:3, :3]
    R2 = quat_to_rotmat(quats)
    R = R1 @ R2
    V0 = scales[:,0,None] * R[:,:,0]
    V1 = scales[:,1,None] * R[:,:,1]

    center, bound, valid = project_ellipse_bound(p_view, V0, V1, *intrins)

    tile_min, tile_max = get_tile_bbox(center, bound, tile_bounds, block_width)
    assert (tile_max >= tile_min).all()
    tile_area = (tile_max[...,0]-tile_min[...,0]) * (tile_max[...,1]-tile_min[...,1])
    mask = (tile_area > 0) & (~is_close) & valid
    depth_grads = projected_depth_grad(p_view, R, fx, fy)

    num_tiles_hit = tile_area

    bounds = torch.concatenate((tile_min, tile_max), axis=-1)
    bounds = torch.where(~mask[..., None], 0, bounds.to(torch.int32))
    num_tiles_hit = torch.where(~mask, 0, num_tiles_hit)
    positions = torch.where(~mask[..., None], 0, p_view)
    axes_u = torch.where(~mask[..., None], 0, V0)
    axes_v = torch.where(~mask[..., None], 0, V1)
    depth_grads = torch.where(~mask[..., None], 0, depth_grads)

    return (
        positions, axes_u, axes_v,
        depth_grads,
        bounds, num_tiles_hit,
    )


def map_gaussian_to_intersects(
    num_points, xys, depths, radii, cum_tiles_hit, tile_bounds, block_width
):
    num_intersects = cum_tiles_hit[-1]
    isect_ids = torch.zeros(num_intersects, dtype=torch.int64, device=xys.device)
    gaussian_ids = torch.zeros(num_intersects, dtype=torch.int32, device=xys.device)

    for idx in range(num_points):
        if radii[idx] <= 0:
            break

        tile_min, tile_max = get_tile_bbox(
            xys[idx], radii[idx], tile_bounds, block_width
        )

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

def get_alpha(uv, opac, aniso) -> Tuple[Tensor, bool]:
    r2 = torch.norm(uv)**2
    vis = visibility_kernel(r2)
    t = torch.dot(uv, aniso)
    t = torch.clamp(t, torch.zeros_like(t), torch.ones_like(t))
    m = t*t*(2.0*t-3.0)+1.0
    alpha = opac * vis * m
    return alpha, r2 >= 0.0 and alpha >= 1e-4


def get_intersection(position, axis_uv, pos_2d):
    pos_2d_3 = torch.concat((pos_2d, torch.ones_like(pos_2d[0:1])))
    A = torch.concat((*axis_uv, pos_2d_3)).reshape((3,3)).T
    uvt = -torch.linalg.inv(A) @ position
    uv = uvt[:2]
    t = -uvt[2:]
    return pos_2d_3*t, uv, torch.norm(uv) < 1.0


def rasterize_gaussians_simple(
    positions: Float[Tensor, "*batch 3"],
    axes_u: Float[Tensor, "*batch 3"],
    axes_v: Float[Tensor, "*batch 3"],
    colors: Float[Tensor, "*batch channels"],
    opacities: Float[Tensor, "*batch 1"],
    anisotropies: Float[Tensor, "*batch 2"],
    bounds: Int[Tensor, "*batch 4"],
    num_tiles_hit: Int[Tensor, "*batch 1"],
    intrins: Tuple[float, float, float, float],
    img_height: int,
    img_width: int,
    block_width: int,
    background: Float[Tensor, "channels"] = None
):
    device = positions.device
    float32_param = { 'dtype': torch.float32, 'device': device }
    int32_param = { 'dtype': torch.int32, 'device': device }
    if background is None:
        background = torch.zeros(colors.shape[-1], **float32_param)

    num_points = positions.size(0)
    tile_bounds = (
        (img_width + block_width - 1) // block_width,
        (img_height + block_width - 1) // block_width,
        1,
    )
    block = (block_width, block_width, 1)
    img_size = (img_width, img_height, 1)

    num_intersects, cum_tiles_hit = compute_cumulative_intersects(num_tiles_hit)
    assert num_intersects >= 1

    (
        isect_ids_unsorted,
        gaussian_ids_unsorted,
        isect_ids_sorted,
        gaussian_ids_sorted,
        tile_bins,
    ) = bin_and_sort_gaussians(
        num_points, num_intersects,
        positions,
        bounds, cum_tiles_hit,
        tile_bounds, block_width,
    )

    final_idx = torch.zeros((img_size[1], img_size[0]), **int32_param)
    out_img = torch.zeros((img_size[1], img_size[0], 3), **float32_param)
    out_alpha = torch.zeros((img_size[1], img_size[0]), **float32_param)

    fx, fy, cx, cy = intrins

    for i in range(img_size[1]):
        for j in range(img_size[0]):
            tile_id = (i // block[0]) * tile_bounds[0] + (j // block[1])
            tile_bin_start = tile_bins[tile_id, 0]
            tile_bin_end = tile_bins[tile_id, 1]
            pos_screen = [j+0.5, i+0.5]
            pos_2d = [(pos_screen[0]-cx)/fx, (pos_screen[1]-cy)/fy]
            pos_2d = torch.tensor(pos_2d, **float32_param)

            T = 1.0
            cur_idx = 0
            pix_out = torch.zeros(3, **float32_param)

            for idx in range(tile_bin_start, tile_bin_end):
                gid = gaussian_ids_sorted[idx]
                pos = positions[gid]
                color  = colors[gid]
                opac = opacities[gid]
                aniso = anisotropies[gid]
                axis_uv = (axes_u[gid], axes_v[gid])

                poi, uv, valid = get_intersection(pos, axis_uv, pos_2d)
                if torch.linalg.norm(uv) > 1.0:
                    continue
                alpha, valid  = get_alpha(uv, opac, aniso)
                if  not valid:
                    continue
                
                next_T = T * (1. - alpha)

                vis = alpha * T
                pix_out += vis * color
                T = next_T
                cur_idx = idx
                if next_T <= 1e-3:
                    break

            final_idx[i, j] = cur_idx
            out_img[i, j] = pix_out + T * background
            out_alpha[i, j] = 1.0 - T

    return out_img, out_alpha, final_idx


def rasterize_gaussians_depth(
    positions: Float[Tensor, "*batch 3"],
    axes_u: Float[Tensor, "*batch 3"],
    axes_v: Float[Tensor, "*batch 3"],
    opacities: Float[Tensor, "*batch 1"],
    anisotropies: Float[Tensor, "*batch 2"],
    bounds: Int[Tensor, "*batch 4"],
    num_tiles_hit: Int[Tensor, "*batch 1"],
    intrins: Tuple[float, float, float, float],
    img_height: int,
    img_width: int,
    block_width: int,
):
    device = positions.device
    float32_param = { 'dtype': torch.float32, 'device': device }
    int32_param = { 'dtype': torch.int32, 'device': device }

    num_points = positions.size(0)
    tile_bounds = (
        (img_width + block_width - 1) // block_width,
        (img_height + block_width - 1) // block_width,
        1,
    )
    block = (block_width, block_width, 1)
    img_size = (img_width, img_height, 1)

    num_intersects, cum_tiles_hit = compute_cumulative_intersects(num_tiles_hit)
    assert num_intersects >= 1

    (
        isect_ids_unsorted,
        gaussian_ids_unsorted,
        isect_ids_sorted,
        gaussian_ids_sorted,
        tile_bins,
    ) = bin_and_sort_gaussians(
        num_points, num_intersects,
        positions,
        bounds, cum_tiles_hit,
        tile_bounds, block_width,
    )

    final_idx = torch.zeros((img_size[1], img_size[0]), **int32_param)
    out_depth = torch.zeros((img_size[1], img_size[0], 1), **float32_param)
    out_visibility = torch.zeros((img_size[1], img_size[0], 2), **float32_param)

    fx, fy, cx, cy = intrins

    for i in range(img_size[1]):
        for j in range(img_size[0]):
            tile_id = (i // block[0]) * tile_bounds[0] + (j // block[1])
            tile_bin_start = tile_bins[tile_id, 0]
            tile_bin_end = tile_bins[tile_id, 1]
            pos_screen = [j+0.5, i+0.5]
            pos_2d = [(pos_screen[0]-cx)/fx, (pos_screen[1]-cy)/fy]
            pos_2d = torch.tensor(pos_2d, **float32_param)

            T = 1.0
            interp = 1.0
            cur_idx = 0
            median_depth = 0.0

            for idx in range(tile_bin_start, tile_bin_end):
                gid = gaussian_ids_sorted[idx]
                pos = positions[gid]
                opac = opacities[gid]
                aniso = anisotropies[gid]
                axis_uv = (axes_u[gid], axes_v[gid])

                poi, uv, valid = get_intersection(pos, axis_uv, pos_2d)
                if torch.linalg.norm(uv) > 1.0:
                    continue
                alpha, valid  = get_alpha(uv, opac, aniso)
                if  not valid:
                    continue
                
                next_T = T * (1. - alpha)

                next_depth = poi[2]

                median_th = 0.5
                if next_T < median_th:
                    if T < 0.99999:
                        interp = (1.0-alpha)/alpha * (T-median_th)/median_th
                        median_depth = median_depth + (next_depth-median_depth)*interp
                    else:
                        median_depth = next_depth
                    T = next_T
                    cur_idx = idx
                    break

                median_depth = next_depth
                T = next_T
                cur_idx = idx

            final_idx[i, j] = cur_idx
            out_depth[i, j] = median_depth
            out_visibility[i, j, 0] = T
            out_visibility[i, j, 1] = interp

    return out_depth, out_visibility, final_idx


def rasterize_gaussians( 
    positions: Float[Tensor, "*batch 3"],
    axes_u: Float[Tensor, "*batch 3"],
    axes_v: Float[Tensor, "*batch 3"],
    colors: Float[Tensor, "*batch channels"],
    ch_degree_r: int,
    ch_degree_r_to_use: int,
    ch_degree_phi: int,
    ch_degree_phi_to_use: int,
    ch_coeffs: Float[Tensor, "*batch dim_ch 3"],
    opacities: Float[Tensor, "*batch 1"],
    anisotropies: Float[Tensor, "*batch 2"],
    depth_grads: Float[Tensor, "*batch 2"],
    depth_ref_im: Float[Tensor, "*batch 2"],
    bounds: Int[Tensor, "*batch 4"],
    num_tiles_hit: Int[Tensor, "*batch 1"],
    intrins: Tuple[float, float, float, float],
    img_height: int,
    img_width: int,
    block_width: int,
    background: Optional[Float[Tensor, "channels"]] = None
):
    device = positions.device
    float32_param = { 'dtype': torch.float32, 'device': device }
    int32_param = { 'dtype': torch.int32, 'device': device }
    if background is None:
        background = torch.zeros(colors.shape[-1], **float32_param)

    num_points = positions.size(0)
    tile_bounds = (
        (img_width + block_width - 1) // block_width,
        (img_height + block_width - 1) // block_width,
        1,
    )
    block = (block_width, block_width, 1)
    img_size = (img_width, img_height, 1)

    num_intersects, cum_tiles_hit = compute_cumulative_intersects(num_tiles_hit)
    assert num_intersects >= 1

    (
        isect_ids_unsorted,
        gaussian_ids_unsorted,
        isect_ids_sorted,
        gaussian_ids_sorted,
        tile_bins,
    ) = bin_and_sort_gaussians(
        num_points, num_intersects,
        positions,
        bounds, cum_tiles_hit,
        tile_bounds, block_width,
    )

    final_idx = torch.zeros((img_size[1], img_size[0]), **int32_param)
    out_alpha = torch.zeros((img_size[1], img_size[0]), **float32_param)
    out_img = torch.zeros((img_size[1], img_size[0], 3), **float32_param)
    out_depth_grad = torch.zeros((img_size[1], img_size[0], 4), **float32_param)
    out_reg_depth = torch.zeros((img_size[1], img_size[0]), **float32_param)
    out_reg_normal = torch.zeros((img_size[1], img_size[0]), **float32_param)
    float32_param['requires_grad'] = True

    fx, fy, cx, cy = intrins
    dim_ch = ch_degree_r * (2*ch_degree_phi+1)

    for i in range(img_size[1]):
        for j in range(img_size[0]):
            tile_id = (i // block[0]) * tile_bounds[0] + (j // block[1])
            tile_bin_start = tile_bins[tile_id, 0]
            tile_bin_end = tile_bins[tile_id, 1]
            pos_screen = [j+0.5, i+0.5]
            pos_2d = [(pos_screen[0]-cx)/fx, (pos_screen[1]-cy)/fy]
            pos_2d = torch.tensor(pos_2d, **float32_param)

            T = torch.ones(1, **float32_param)
            cur_idx = 0
            g_sum = torch.zeros(2, **float32_param)
            pix_out = torch.zeros(3, **float32_param)
            vis_sum = torch.zeros(1, **float32_param)
            depth_sum = torch.zeros(1, **float32_param)
            depth_squared_sum = torch.zeros(1, **float32_param)
            depth_normal_ref = depth_ref_im[i, j][:2]
            depth_ref = depth_ref_im[i, j][2]
            reg_depth = torch.zeros(1, **float32_param)
            reg_normal = torch.zeros(1, **float32_param)

            for idx in range(tile_bin_start, tile_bin_end):
                gid = gaussian_ids_sorted[idx]
                pos = positions[gid]
                opac = opacities[gid]
                aniso = anisotropies[gid]
                axis_uv = (axes_u[gid], axes_v[gid])

                poi, uv, valid = get_intersection(pos, axis_uv, pos_2d)
                if torch.linalg.norm(uv) > 1.0:
                    continue
                alpha, valid  = get_alpha(uv, opac, aniso)
                if  not valid:
                    continue

                color_0 = colors[gid]
                color = color_0
                if dim_ch > 0:
                    ch_color = ch_coeffs_to_color(
                        ch_degree_r, ch_degree_r_to_use,
                        ch_degree_phi, ch_degree_phi_to_use,
                        ch_coeffs[gid], uv
                    )
                    # ch_color = ch_coeffs[gid][0]
                    color = color_0 * torch.sigmoid(ch_color)

                next_T = T * (1. - alpha)

                vis = alpha * T
                # depth = pos[2]
                depth = poi[2]
                g_i = depth_grads[gid]
                n_i = g_i / (torch.norm(g_i)+1e-6)

                pix_out = pix_out + vis * color
                # reg_depth = reg_depth + vis*depth * vis_sum - vis * depth_sum
                # reg_depth = reg_depth + vis * (vis_sum*depth**2 + depth_squared_sum - 2.0*depth*depth_sum)
                reg_depth = reg_depth + vis * abs(depth-depth_ref)
                # reg_depth = reg_depth + vis * (depth-depth_ref)**2
                reg_normal = reg_normal + vis * (1.0 - torch.dot(n_i, depth_normal_ref))
                vis_sum = vis_sum + vis
                depth_sum = depth_sum + vis*depth
                depth_squared_sum = depth_squared_sum + vis*depth**2
                g_sum = g_sum + vis * g_i

                T = next_T
                cur_idx = idx
                if next_T <= 1e-3:
                    break

            final_idx[i, j] = cur_idx
            out_alpha[i, j] = 1.0-T
            out_img[i, j] = pix_out + T * background
            out_depth_grad[i, j] = torch.concat((g_sum, depth_sum, depth_squared_sum))
            out_reg_depth[i, j] = reg_depth
            out_reg_normal[i, j] = reg_normal

    return (
        out_img, out_depth_grad,
        out_reg_depth, out_reg_normal,
        out_alpha, final_idx
    )
