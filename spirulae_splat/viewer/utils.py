import torch
import torch.nn.functional as F
import numpy as np

from typing import Literal, Tuple, Optional

from spirulae_splat.splat._torch_impl import quat_to_rotmat


def generate_right_angle_rotation_matrices() -> np.ndarray:
    """List all "right angle" rotation matrices"""
    import itertools
    valid_matrices = []
    for matrix in itertools.product([-1, 0, 1], repeat=9):
        R = np.array(matrix).reshape(3, 3)
        if np.linalg.det(R) == 1 and (np.dot(R.T, R) == np.eye(3)).all():
            valid_matrices.append(R)
    return valid_matrices


def generate_sh_basis(dirs: torch.Tensor, mode: Literal['nerfstudio', 'gsplat']) -> torch.Tensor:
    x, y, z = torch.unbind(dirs, -1)
    one = torch.ones_like(x)
    xx, yy, zz = x*x, y*y, z*z

    if mode == "nerfstudio":
        return torch.stack([
            0.28209479177387814 * one,
            0.4886025119029199 * y,
            0.4886025119029199 * z,
            0.4886025119029199 * x,
            1.0925484305920792 * x * y,
            1.0925484305920792 * y * z,
            (0.9461746957575601 * zz - 0.31539156525251999),
            1.0925484305920792 * x * z,
            0.5462742152960396 * (xx - yy),
            0.5900435899266435 * y * (3.0 * xx - yy),
            2.890611442640554 * x * y * z,
            0.4570457994644658 * y * (5.0 * zz - 1.0),
            0.3731763325901154 * z * (5.0 * zz - 3.0),
            0.4570457994644658 * x * (5.0 * zz - 1.0),
            1.445305721320277 * z * (xx - yy),
            0.5900435899266435 * x * (xx - 3.0 * yy),
            2.5033429417967046 * x * y * (xx - yy),
            1.7701307697799304 * y * z * (3.0 * xx - yy),
            0.9461746957575601 * x * y * (7.0 * zz - 1.0),
            0.6690465435572892 * y * z * (7.0 * zz - 3.0),
            0.10578554691520431 * (35.0 * zz * zz - 30.0 * zz + 3.0),
            0.6690465435572892 * x * z * (7.0 * zz - 3.0),
            0.47308734787878004 * (xx - yy) * (7.0 * zz - 1.0),
            1.7701307697799304 * x * z * (xx - 3.0 * yy),
            0.6258357354491761 * (xx * (xx - 3.0 * yy) - yy * (3.0 * xx - yy)),
        ], dim=1)

    elif mode == "gsplat":
        from spirulae_splat.splat._torch_impl import eval_sh_bases
        return eval_sh_bases(25, dirs)


def rotate_sh_coeffs(coeffs: torch.Tensor, R: torch.Tensor,
                     mode: Literal['nerfstudio', 'gsplat']) -> torch.Tensor:
    if isinstance(R, np.ndarray):
        R = torch.from_numpy(R).to(coeffs)

    # fibonacci sample on a sphere
    M = coeffs.shape[1]
    N = int(1.2*M+1)  # >= M, more is slower but more numerically stable
    idx = torch.arange(N, dtype=torch.float32, device=coeffs.device) + 0.5
    phi = 2 * np.pi * idx / ((1 + np.sqrt(5)) / 2)
    cos_theta = 1.0 - 2.0 * idx / N
    sin_theta = torch.sqrt(torch.relu(1.0 - cos_theta * cos_theta))
    dirs = torch.stack([sin_theta * torch.cos(phi),
                        sin_theta * torch.sin(phi),
                        cos_theta], dim=1)  # (N, 3)

    # evaluate basis
    Y_src = generate_sh_basis(dirs, mode)[:, 1:M+1]  # (N, M)
    dirs_rot = dirs @ R.t()  # (N, 3)
    Y_rot = generate_sh_basis(dirs_rot, mode)[:, 1:M+1]  # (N, M)
    # coeffs: (batch, M, 3)

    # change of basis using direct linear transform
    f_vals = torch.matmul(Y_src, coeffs)  # (batch, M, 3)
    Y_rot_inv = torch.linalg.inv(Y_rot.T @ Y_rot) @ Y_rot.T  # (N, M)
    coeffs_rot = torch.matmul(Y_rot_inv, f_vals)  # (batch, N, 3)

    # numerical check
    if False:
        print(torch.linalg.eigvalsh(Y_rot.T @ Y_rot))
        f_vals_rot = torch.einsum('ni,bic->bnc', Y_rot, coeffs_rot)
        print(abs(f_vals - f_vals_rot).max().item())

    return coeffs_rot


def quat_to_rotmat(quat: torch.Tensor) -> torch.Tensor:
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

def rotmat_to_quat(R: torch.tensor):
    tr = R[..., 0, 0] + R[..., 1, 1] + R[..., 2, 2]
    qw = 0.5 * torch.sqrt(torch.clamp(tr + 1.0, min=0.0))
    qx = torch.where(
        qw > 1e-8,
        (R[..., 2, 1] - R[..., 1, 2]) / (4 * qw),
        torch.sqrt(torch.clamp(1 + R[..., 0, 0] - R[..., 1, 1] - R[..., 2, 2], min=0.0)) / 2,
    )
    qy = torch.where(
        qw > 1e-8,
        (R[..., 0, 2] - R[..., 2, 0]) / (4 * qw),
        (R[..., 0, 1] + R[..., 1, 0]) / (4 * qx + 1e-12),
    )
    qz = torch.where(
        qw > 1e-8,
        (R[..., 1, 0] - R[..., 0, 1]) / (4 * qw),
        (R[..., 0, 2] + R[..., 2, 0]) / (4 * qx + 1e-12),
    )

    quats = torch.stack([qw, qx, qy, qz], dim=-1)
    return F.normalize(quats, dim=-1)


def quat_scale_to_triangle_verts(
        quats: torch.Tensor,
        scales: torch.Tensor,
        means: Optional[torch.Tensor] = None,
        features_dc: Optional[torch.Tensor] = None,
        features_ch: Optional[torch.Tensor] = None
    ) -> torch.Tensor:
    quats = quats / torch.norm(quats, dim=-1, keepdim=True)
    M = quat_to_rotmat(quats)
    
    sx = torch.exp(scales[..., 0])
    sy = torch.exp(scales[..., 1])
    sz = scales[..., 2] - 0.5*(scales[..., 0]+scales[..., 1])
    zero = torch.zeros_like(sx)
    verts = torch.stack([
        sx, zero, zero,
        sx*(-0.5+sz), sy, zero,
        sx*(-0.5-sz), -sy, zero
    ], dim=-1).reshape(scales.shape[:-1] + (3, 3))

    verts = torch.bmm(verts, M.transpose(-1, -2))
    if means is not None:
        verts = verts + means[..., None, :]

    if features_dc is not None:
        features_dc = 0.28209479177387814 * features_dc + 0.5
        colors = [features_dc, features_dc, features_dc]
        if features_ch is not None:
            colors[0] = colors[0] + features_ch[..., 0, :]
            colors[1] = colors[1] + -0.5 * features_ch[..., 0, :] + 0.75**0.5 * features_ch[..., 1, :]
            colors[2] = colors[2] + -0.5 * features_ch[..., 0, :] - 0.75**0.5 * features_ch[..., 1, :]
        return verts, torch.stack(colors, dim=-2)

    return verts


def split_triangles(verts: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
    """
    Split each triangle into two along its longest edge,
    choosing the split point near the perpendicular projection of the opposite vertex
    (to avoid needle-like triangles).
    verts: (..., 3, 3)
    Returns: (tri1, tri2) each (..., 3, 3)
    """
    a, b, c = verts[..., 0, :], verts[..., 1, :], verts[..., 2, :]

    e01 = torch.norm(b - a, dim=-1)
    e12 = torch.norm(c - b, dim=-1)
    e20 = torch.norm(a - c, dim=-1)

    edges = torch.stack([e01, e12, e20], dim=-1)
    longest = torch.argmax(edges, dim=-1)

    # Prepare containers
    tri1 = verts.clone()
    tri2 = verts.clone()

    def split_edge(a, b, c):
        ab = b - a
        ab_len2 = (ab ** 2).sum(dim=-1, keepdim=True)
        # projection factor t of c onto line ab
        t = ((c - a) * ab).sum(dim=-1, keepdim=True) / (ab_len2 + 1e-9)
        t = t.clamp(0.25, 0.75)
        p = a + t * ab
        return p

    # edge (0,1), opposite c
    mask = (longest == 0)[..., None, None]
    p01 = split_edge(a, b, c)
    tri1 = torch.where(mask, torch.stack([a, p01, c], dim=-2), tri1)
    tri2 = torch.where(mask, torch.stack([b, p01, c], dim=-2), tri2)

    # edge (1,2), opposite a
    mask = (longest == 1)[..., None, None]
    p12 = split_edge(b, c, a)
    tri1 = torch.where(mask, torch.stack([b, p12, a], dim=-2), tri1)
    tri2 = torch.where(mask, torch.stack([c, p12, a], dim=-2), tri2)

    # edge (2,0), opposite b
    mask = (longest == 2)[..., None, None]
    p20 = split_edge(c, a, b)
    tri1 = torch.where(mask, torch.stack([c, p20, b], dim=-2), tri1)
    tri2 = torch.where(mask, torch.stack([a, p20, b], dim=-2), tri2)

    return tri1, tri2


def triangle_verts_to_quat_scale_mean(
        verts: torch.Tensor
    ) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """
    Invert `quat_scale_to_triangle_verts`.
    verts: (..., 3, 3)
    Returns quats (..., 4), scales (..., 3), means (..., 3)
    """
    means = verts.mean(dim=-2)
    centered = verts - means[..., None, :]

    # Get rotation from normal and edge directions
    e0 = centered[..., 1, :] - centered[..., 0, :]
    e1 = centered[..., 2, :] - centered[..., 0, :]
    n = torch.cross(e0, e1, dim=-1)
    n = n / torch.norm(n, dim=-1, keepdim=True)

    xdir = centered[..., 0, :] / torch.norm(centered[..., 0, :], dim=-1, keepdim=True)
    zdir = n
    ydir = torch.cross(zdir, xdir, dim=-1)
    ydir = ydir / torch.norm(ydir, dim=-1, keepdim=True)
    zdir = torch.cross(xdir, ydir, dim=-1)

    M = torch.stack([xdir, ydir, zdir], dim=-1)  # (..., 3, 3)
    quats = rotmat_to_quat(M)

    # Transform verts into local frame
    local = torch.matmul(centered, M)
    v0, v1, v2 = local[..., 0, :], local[..., 1, :], local[..., 2, :]

    sx = v0[..., 0]
    sy = 0.5 * (v1[..., 1] - v2[..., 1])
    sz = 0.5 * ((v1[..., 0] - v2[..., 0]) / (sx + 1e-9))
    s0 = torch.log(sx)
    s1 = torch.log(sy)
    s2 = sz + 0.5 * (s0 + s1)
    scales = torch.stack([s0, s1, s2], dim=-1)

    return quats, scales, means


if __name__ == "__main__":
    torch.random.manual_seed(42)
    N = 10
    quats = torch.randn(N, 4)
    scales = torch.randn(N, 3)
    means = torch.randn(N, 3)

    verts = quat_scale_to_triangle_verts(quats, scales, means)
    quats1, scales1, means1 = triangle_verts_to_quat_scale_mean(verts)

    quats = quats / torch.norm(quats, dim=-1, keepdim=True)
    print((quats*quats1).sum(-1))
    print(scales-scales1)
    print(means-means1)

    def print_trig(verts):
        print([tuple(x) for x in verts.cpu().numpy().tolist()])
    print_trig(verts[0])
    verts1, verts2 = split_triangles(verts)
    print_trig(verts1[0])
    print_trig(verts2[0])
