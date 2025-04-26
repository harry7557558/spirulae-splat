import pytest
import torch
import time

from spirulae_splat.splat import _torch_impl
from spirulae_splat.splat.sh import num_sh_bases, spherical_harmonics


device = torch.device("cuda:0")


def _test_sh(method):
    torch.random.manual_seed(42)

    num_points = 1024
    degree = 4
    gt_colors = torch.ones(num_points, 3, device=device) * 0.5
    viewdirs = torch.randn(num_points, 3, device=device)
    viewdirs /= torch.linalg.norm(viewdirs, dim=-1, keepdim=True)
    sh_coeffs = torch.rand(
        num_points, num_sh_bases(degree), 3, device=device, requires_grad=True
    )
    sh_coeffs_0 = sh_coeffs[:, 0, :].contiguous().detach().requires_grad_(True)
    sh_coeffs_1 = sh_coeffs[:, 1:, :].contiguous().detach().requires_grad_(True)
    optim_check = torch.optim.Adam([sh_coeffs], lr=1e-2)
    optim = torch.optim.Adam([sh_coeffs_0, sh_coeffs_1], lr=1e-2)

    num_iters = 1000
    for _ in range(num_iters):
        optim_check.zero_grad()
        optim.zero_grad()

        # compute PyTorch's color and grad
        check_colors = _torch_impl.compute_sh_color(viewdirs, sh_coeffs, method)
        check_loss = torch.square(check_colors - gt_colors).mean()
        check_loss.backward()
        check_grad = sh_coeffs.grad.detach()

        # compute our colors and grads
        colors = spherical_harmonics(degree, viewdirs, sh_coeffs_0, sh_coeffs_1, method)
        loss = torch.square(colors - gt_colors).mean()
        loss.backward()
        grad = torch.cat((sh_coeffs_0.grad.unsqueeze(1), sh_coeffs_1.grad), dim=1).detach()

        torch.testing.assert_close(check_colors, colors)
        torch.testing.assert_close(check_grad, grad)

        optim_check.step()
        optim.step()

    # check final optimized color
    torch.testing.assert_close(check_colors, gt_colors)


def profile_sh(
    method, num_points: int = 1_000_000, degree: int = 4, n_iters: int = 1000
):
    torch.random.manual_seed(42)

    viewdirs = torch.randn(num_points, 3, device=device)
    viewdirs /= torch.linalg.norm(viewdirs, dim=-1, keepdim=True)
    sh_coeffs = torch.rand(
        num_points, num_sh_bases(degree), 3, device=device, requires_grad=True
    )
    sh_coeffs_0 = sh_coeffs[:, 0, :].contiguous().detach().requires_grad_(True)
    sh_coeffs_1 = sh_coeffs[:, 1:, :].contiguous().detach().requires_grad_(True)

    for _ in range(10):  # warmup
        colors = spherical_harmonics(degree, viewdirs, sh_coeffs_0, sh_coeffs_1, method)

    n_iters_fwd = n_iters
    torch.cuda.synchronize()
    tic = time.time()
    for _ in range(n_iters_fwd):
        _ = spherical_harmonics(degree, viewdirs, sh_coeffs_0, sh_coeffs_1, method)
    torch.cuda.synchronize()
    toc = time.time()
    ellipsed = (toc - tic) / n_iters_fwd * 1000  # ms
    print(f"[Fwd] Method: {method}, ellipsed: {ellipsed:.2f} ms")

    loss = spherical_harmonics(degree, viewdirs, sh_coeffs_0, sh_coeffs_1, method).sum()
    for _ in range(10):  # warmup
        loss.backward(retain_graph=True)

    n_iters_bwd = n_iters // 20
    torch.cuda.synchronize()
    tic = time.time()
    for _ in range(n_iters_bwd):
        loss.backward(retain_graph=True)
    torch.cuda.synchronize()
    toc = time.time()
    ellipsed = (toc - tic) / n_iters_bwd * 1000  # ms
    print(f"[Bwd] Method: {method}, ellipsed: {ellipsed:.2f} ms")


@pytest.mark.skipif(not torch.cuda.is_available(), reason="No CUDA device")
def test_sh_poly():
    _test_sh("poly")

@pytest.mark.skipif(not torch.cuda.is_available(), reason="No CUDA device")
def test_sh_fast():
    _test_sh("poly")

if __name__ == "__main__":
    test_sh_poly()
    test_sh_fast()
    profile_sh("poly")
    profile_sh("fast")
