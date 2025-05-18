import torch
import torch.nn.functional as F
import numpy as np

from typing import Literal


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
