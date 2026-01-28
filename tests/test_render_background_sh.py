import pytest
import torch
from torch.func import vjp  # type: ignore

from spirulae_splat.splat import _torch_impl
from spirulae_splat.splat.background_sh import render_background_sh
import spirulae_splat.splat.cuda as _C
from spirulae_splat.splat._camera import _Camera

from utils import check_close, timeit

torch.manual_seed(42)

device = torch.device("cuda:0")


MODEL = "pinhole"  # pinhole or fisheye


@pytest.mark.skipif(not torch.cuda.is_available(), reason="No CUDA device")
def test_render_background_sh():

    B, H, W = 3, 500, 300
    cx, cy = 0.45*W, 0.55*H
    fx, fy = 1.5*W, 1.6*W

    intrins = torch.tensor([[fx, fy, cx, cy]]).repeat(B, 1).to(device)
    intrins *= torch.exp(0.2*torch.randn_like(intrins))

    rotation = torch.randn(B, 4).to(device)
    rotation /= torch.norm(rotation, dim=-1, keepdim=True)
    rotation = _torch_impl.quat_to_rotmat(rotation)
    rotation = rotation.detach().contiguous()
    _rotation = rotation.clone().requires_grad_(True)
    rotation.requires_grad_(True)

    sh_degree = 4
    sh_coeffs = torch.randn(((sh_degree+1)**2, 3)).to(device)
    _sh_coeffs = sh_coeffs.clone().requires_grad_(True)
    sh_coeffs.requires_grad_(True)

    output = render_background_sh(
        W, H, MODEL,
        intrins, rotation, sh_degree, sh_coeffs
    )
    _output = torch.stack([
        _torch_impl.render_background_sh(
            W, H,
            (intrins_i[0].item(), intrins_i[1].item(), intrins_i[2].item(), intrins_i[3].item()),
            rotation_i, sh_degree, _sh_coeffs,
        ) for (intrins_i, rotation_i) in zip(intrins, _rotation)
    ])

    print("test forward")
    check_close('image', output, _output)
    print()

    weight = torch.randn_like(output)
    def fun(output):
        return (weight*output).mean()
    fun(output).backward()
    fun(_output).backward()

    print("test backward")
    tol = { 'atol': 1e-6, 'rtol': 1e-4 }
    # print(rotation.grad)
    # print(_rotation.grad)
    check_close('rotation', rotation.grad, _rotation.grad, **tol)
    # print(sh_coeffs.grad)
    # print(_sh_coeffs.grad)
    check_close('sh_coeffs', sh_coeffs.grad, _sh_coeffs.grad, **tol)
    print()


@pytest.mark.skipif(not torch.cuda.is_available(), reason="No CUDA device")
def profile_render_background_sh():

    # B, H, W = 1, 1440, 1080
    B, H, W = 3, 1920, 1440
    cx, cy = 0.45*W, 0.55*H
    fx, fy = 1.5*W, 1.6*W

    intrins = torch.tensor([[fx, fy, cx, cy]]).repeat(B, 1).to(device)
    intrins *= torch.exp(0.2*torch.randn_like(intrins))

    rotation = torch.randn(B, 4).to(device)
    rotation /= torch.linalg.norm(rotation)
    rotation = _torch_impl.quat_to_rotmat(rotation)
    rotation = rotation.detach().contiguous()
    _rotation = rotation.clone().requires_grad_(True)
    rotation.requires_grad_(True)

    sh_degree = 4
    sh_coeffs = torch.randn(((sh_degree+1)**2, 3)).to(device)
    _sh_coeffs = sh_coeffs.clone().requires_grad_(True)
    sh_coeffs.requires_grad_(True)

    print("profile forward")

    output = render_background_sh(
        W, H, MODEL,
        intrins, rotation, sh_degree, sh_coeffs
    )[0]
    # _output = _torch_impl.render_background_sh(
    #     W, H, (fx, fy, cx, cy),
    #     _rotation, sh_degree, _sh_coeffs,
    # )

    # timeit(lambda: _torch_impl.render_background_sh(
    #     W, H, (fx, fy, cx, cy),
    #     _rotation, sh_degree, _sh_coeffs,
    # ), "forward torch")
    timeit(lambda: render_background_sh(
        W, H, MODEL,
        intrins[None], rotation[None], sh_degree, sh_coeffs
    )[0], "forward fused")

    print()

    print("profile backward")

    weights = torch.randn_like(output)
    loss= (weights * output).sum()
    # _loss = (weights * _output).sum()
    output.retain_grad()
    # _output.retain_grad()
    loss.backward(retain_graph=True)
    # _loss.backward(retain_graph=True)

    # timeit(lambda: _loss.backward(retain_graph=True), "backward torch")
    timeit(lambda: loss.backward(retain_graph=True), "backward fused", 20)

    print()


if __name__ == "__main__":

    test_render_background_sh()
    print()

    profile_render_background_sh()
