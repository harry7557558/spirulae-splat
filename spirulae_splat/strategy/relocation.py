from typing import Optional, Tuple

import math
import torch
from torch import Tensor
import spirulae_splat.splat.cuda as _C

N_MAX = 51
BINOMS = torch.zeros((N_MAX, N_MAX)).float().cuda()
for n in range(N_MAX):
    for k in range(n + 1):
        BINOMS[n, k] = math.comb(n, k)


def compute_relocation(
    opacities: Tensor,  # [N]
    scales: Tensor,  # [N, k]
    ratios: Tensor,  # [N]
) -> Tuple[Tensor, Tensor]:
    """Compute new Gaussians from a set of old Gaussians.
    This function interprets the Gaussians as samples from a likelihood distribution.
    It uses the old opacities and scales to compute the new opacities and scales.
    This is an implementation of the paper
    `3D Gaussian Splatting as Markov Chain Monte Carlo <https://arxiv.org/pdf/2404.09591>`_,
    Args:
        opacities: The opacities of the Gaussians. [N]
        scales: The scales of the Gaussians. [N, k]
        ratios: The relative frequencies for each of the Gaussians. [N]
    Returns:
        A tuple:
        **new_opacities**: The opacities of the new Gaussians. [N]
        **new_scales**: The scales of the Gaussians. [N, k]
    """
    N = opacities.shape[0]
    assert scales.shape[0] == N, scales.shape
    assert ratios.shape == (N,), ratios.shape
    opacities = opacities.contiguous()
    scales = scales.contiguous()
    ratios = ratios.int().contiguous()
    assert (ratios >= 2).all()

    new_opacities, new_scales = _C.compute_relocation(
        opacities, scales, ratios, BINOMS.to(opacities.device), N_MAX
    )
    return new_opacities, new_scales



def compute_relocation_split(
    positions: Tensor,  # [N, 3]
    quats: Tensor,  # [N, 4]
    opacities: Tensor,  # [N]
    scales: Tensor,  # [N, 2]
) -> Tuple[Tensor, Tensor, Tensor]:
    N = opacities.shape[0]
    assert positions.shape == (N, 3), positions.shape
    assert quats.shape == (N, 4), quats.shape
    assert scales.shape == (N, 2), scales.shape

    return _C.compute_relocation_split(
        positions.contiguous(),
        quats.contiguous(),
        opacities.contiguous(),
        scales.contiguous(),
    )
