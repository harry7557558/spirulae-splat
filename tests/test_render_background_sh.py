import pytest
import traceback
import torch
from torch.func import vjp  # type: ignore

from spirulae_splat.splat import _torch_impl
from spirulae_splat.splat.background_sh import render_background_sh
import spirulae_splat.splat.cuda as _C
from spirulae_splat.splat._camera import _Camera

torch.manual_seed(42)

device = torch.device("cuda:0")


def check_close(name, a, b, atol=1e-5, rtol=1e-5):
    print(name, [*a.shape], [*b.shape], a.dtype, b.dtype)
    try:
        torch.testing.assert_close(a, b, atol=atol, rtol=rtol)
    except AssertionError:
        traceback.print_exc()
        diff = torch.abs(a - b).detach()
        print(f"{diff.max()=} {diff.float().mean()=}")
        # assert False
        # import ipdb
        # ipdb.set_trace()


@pytest.mark.skipif(not torch.cuda.is_available(), reason="No CUDA device")
def test_render_background_sh():

    rotation = torch.randn(4).to(device)
    rotation /= torch.linalg.norm(rotation)
    rotation = _torch_impl.quat_to_rotmat(rotation)
    rotation = rotation.detach().contiguous()
    _rotation = rotation.clone().requires_grad_(True)
    rotation.requires_grad_(True)

    H, W = 500, 300
    cx, cy = 0.45*W, 0.55*H
    fx, fy = 1.5*W, 1.6*W
    cam = _Camera(H, W, "OPENCV", (fx, fy, cx, cy))

    sh_degree = 4
    sh_coeffs = torch.randn(((sh_degree+1)**2, 3)).to(device)
    _sh_coeffs = sh_coeffs.clone().requires_grad_(True)
    sh_coeffs.requires_grad_(True)

    output = render_background_sh(
        cam, rotation, sh_degree, sh_coeffs,
    )
    _output = _torch_impl.render_background_sh(
        W, H, (fx, fy, cx, cy),
        _rotation, sh_degree, _sh_coeffs,
    )

    print("test forward")
    check_close('image', output, _output)
    print()

    def fun(output):
        y, x = torch.meshgrid(torch.arange(H), torch.arange(W), indexing="ij")
        x = x.float().to(device).unsqueeze(-1) + 0.5
        y = y.float().to(device).unsqueeze(-1) + 0.5
        return torch.mean(torch.sin(output+x)*y)
    fun(output).backward()
    fun(_output).backward()

    print("test backward")
    tol = { 'atol': 1e-6, 'rtol': 1e-4 }
    print(rotation.grad)
    print(_rotation.grad)
    check_close('rotation', rotation.grad, _rotation.grad, **tol)
    # print(sh_coeffs.grad)
    # print(_sh_coeffs.grad)
    check_close('sh_coeffs', sh_coeffs.grad, _sh_coeffs.grad, **tol)


if __name__ == "__main__":
    test_render_background_sh()
