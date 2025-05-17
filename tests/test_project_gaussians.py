import pytest
import traceback
import torch
import torch.nn.functional as F

from spirulae_splat.splat import _torch_impl
from spirulae_splat.splat.project_gaussians import project_gaussians
import spirulae_splat.splat.cuda as _C
from spirulae_splat.splat._camera import _Camera

from test_render_background_sh import timeit

torch.manual_seed(42)

device = torch.device("cuda:0")


def projection_matrix(fx, fy, W, H, n=0.01, f=1000.0):
    return torch.tensor(
        [
            [2.0 * fx / W, 0.0, 0.0, 0.0],
            [0.0, 2.0 * fy / H, 0.0, 0.0],
            [0.0, 0.0, (f + n) / (f - n), -2 * f * n / (f - n)],
            [0.0, 0.0, 1.0, 0.0],
        ],
        device=device,
    )


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
def test_project_gaussians():
    num_points = 50

    means3d0 = torch.randn((num_points, 3), device=device, requires_grad=True)
    scales0 = torch.exp(torch.randn((num_points, 2), device=device))
    quats0 = torch.randn((num_points, 4), device=device)
    quats0 /= torch.linalg.norm(quats0, dim=-1, keepdim=True)

    H, W = 512, 512
    cx, cy = 0.45*W, 0.55*H
    fx, fy = 1.5*W, 1.6*W
    cam = _Camera(H, W, "OPENCV", (fx, fy, cx, cy))
    clip_thresh = 0.01
    viewmat = torch.tensor(
        [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 8.0],
            [0.0, 0.0, 0.0, 1.0],
        ],
        device=device,
    )
    viewmat[:3, :3] = _torch_impl.quat_to_rotmat(F.normalize(torch.randn(4), dim=0)).detach()

    params = (means3d0, scales0, quats0, viewmat)
    def decode_params(params):
        params = [*params]
        for i in range(len(params)):
            params[i] = torch.tensor(params[i].detach().cpu().numpy()).to(device)
            try:
                params[i] = params[i].requires_grad_(True)
            except RuntimeError:
                pass
        return params
    means3d, scales, quats, viewmat = decode_params(params)
    _means3d, _scales, _quats, _viewmat = decode_params(params)

    output = project_gaussians(
        means3d, scales, quats,
        viewmat, cam, clip_thresh,
    )
    _output = _torch_impl.project_gaussians(
        _means3d, _scales, _quats,
        _viewmat, (fx, fy, cx, cy), (W, H),
        clip_thresh,
    )

    print("test forward")
    check_close('positions', output[0], _output[0])
    check_close('axes_u', output[1], _output[1])
    check_close('axes_v', output[2], _output[2])
    check_close('bounds', output[3], _output[3])
    check_close('num_tiles_hit', output[4], _output[4])
    print()

    positions_w = torch.randn_like(output[0])
    axes_u_w = torch.randn_like(output[1])
    axes_v_w = torch.randn_like(output[2])

    def fun(output):
        positions, axes_u, axes_v, bounds, num_tiles_hit = output
        return (positions_w*positions).mean() + \
            (axes_u_w*axes_u).mean() + \
            (axes_v_w*axes_v).mean()
    fun(output).backward(retain_graph=True)
    fun(_output).backward(retain_graph=True)

    print("test backward")
    tol = { 'atol': 1e-6, 'rtol': 1e-5 }
    check_close('v_means3d', means3d.grad, _means3d.grad, **tol)
    check_close('v_scales', scales.grad, _scales.grad, **tol)
    check_close('v_quats', quats.grad, _quats.grad, **tol)
    check_close('v_viewmat', viewmat.grad, _viewmat.grad, **tol)
    # print(viewmat.grad)
    # print(_viewmat.grad)

    print()


@pytest.mark.skipif(not torch.cuda.is_available(), reason="No CUDA device")
def profile_project_gaussians():
    num_points = 2_000_000

    means3d0 = torch.randn((num_points, 3), device=device, requires_grad=True)
    scales0 = torch.exp(torch.randn((num_points, 2), device=device))
    quats0 = torch.randn((num_points, 4), device=device)
    quats0 /= torch.linalg.norm(quats0, dim=-1, keepdim=True)

    H, W = 512, 512
    cx, cy = 0.45*W, 0.55*H
    fx, fy = 1.5*W, 1.6*W
    cam = _Camera(H, W, "OPENCV", (fx, fy, cx, cy))
    clip_thresh = 0.01
    viewmat = torch.tensor(
        [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 8.0],
            [0.0, 0.0, 0.0, 1.0],
        ],
        device=device,
    )
    viewmat[:3, :3] = _torch_impl.quat_to_rotmat(F.normalize(torch.randn(4), dim=0)).detach()

    params = (means3d0, scales0, quats0, viewmat)
    def decode_params(params):
        params = [*params]
        for i in range(len(params)):
            params[i] = torch.tensor(params[i].detach().cpu().numpy()).to(device)
            try:
                params[i] = params[i].requires_grad_(True)
            except RuntimeError:
                pass
        return params
    means3d, scales, quats, viewmat = decode_params(params)
    _means3d, _scales, _quats, _viewmat = decode_params(params)

    forward = lambda: project_gaussians(
        means3d, scales, quats,
        viewmat, cam, clip_thresh,
    )
    _forward = lambda: _torch_impl.project_gaussians(
        _means3d, _scales, _quats,
        _viewmat, (fx, fy, cx, cy), (W, H),
        clip_thresh,
    )
    output = forward()
    _output = _forward()

    print("profile forward")
    timeit(_forward, "torch forward")
    timeit(forward, "fused forward")
    print()

    positions_w = torch.randn_like(output[0])
    axes_u_w = torch.randn_like(output[1])
    axes_v_w = torch.randn_like(output[2])

    def fun(output):
        positions, axes_u, axes_v, bounds, num_tiles_hit = output
        return (positions_w*positions).mean() + \
            (axes_u_w*axes_u).mean() + \
            (axes_v_w*axes_v).mean()
    fun(output).backward(retain_graph=True)
    fun(_output).backward(retain_graph=True)

    print("profile backward")
    timeit(lambda: fun(_output).backward(retain_graph=True), "torch backward")
    timeit(lambda: fun(output).backward(retain_graph=True), "fused backward")
    print()


if __name__ == "__main__":
    test_project_gaussians()
    profile_project_gaussians()
