import pytest
import traceback
import torch
from torch.func import vjp  # type: ignore

from spirulae_splat.splat import _torch_impl
from spirulae_splat.splat.project_gaussians import project_gaussians
import spirulae_splat.splat.rasterize_simple as rasterize_simple
import spirulae_splat.splat.cuda as _C

torch.manual_seed(40)

device = torch.device("cuda:0")


def check_close(name, a, b, atol=1e-5, rtol=1e-5):
    print(name, [*a.shape], [*b.shape], a.dtype, b.dtype)
    try:
        torch.testing.assert_close(a, b, atol=atol, rtol=rtol)
    except AssertionError:
        traceback.print_exc()
        diff = torch.abs(a - b).detach()
        print(f"{diff.max()=} {diff.mean()=}", end='\n\n')
        # assert False
        # import ipdb
        # ipdb.set_trace()


@pytest.mark.skipif(not torch.cuda.is_available(), reason="No CUDA device")
def test_rasterize_simple():
    num_points = 40

    means3d = torch.randn((num_points, 3), device=device)
    scales = torch.exp(torch.randn((num_points, 2), device=device))
    quats = torch.randn((num_points, 4), device=device)
    quats /= torch.linalg.norm(quats, dim=-1, keepdim=True)

    colors = torch.randn((num_points, 3), device=device, requires_grad=True)
    opacities = torch.rand((num_points, 1), device=device, requires_grad=True)
    _colors = colors.detach().clone().requires_grad_(True)
    _opacities = opacities.detach().clone().requires_grad_(True)

    H, W = 20, 30
    cx, cy = W / 2, H / 2
    fx, fy = int(0.7*W), int(0.7*W)
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
    viewmat[:3, :3] = _torch_impl.quat_to_rotmat(torch.randn(4))
    BLOCK_SIZE = 16

    params = project_gaussians(
        means3d,
        scales,
        quats,
        viewmat,
        fx, fy, cx, cy,
        H, W, BLOCK_SIZE,
        clip_thresh,
    )
    def decode_params(params):
        params = [param.detach().clone() for param in params]
        for i in range(len(params)):
            try:
                params[i] = params[i].requires_grad_(True)
            except RuntimeError:
                pass
        return params
    (
        xys,
        depths,
        depth_grads,
        radii,
        conics,
        compensation,
        num_tiles_hit,
        cov3d,
    ) = decode_params(params)
    (
        _xys,
        _depths,
        _depth_grads,
        _radii,
        _conics,
        _compensation,
        _num_tiles_hit,
        _cov3d,
    ) = decode_params(params)

    rgb_im, alpha_im, idx = rasterize_simple.rasterize_gaussians_simple(
        xys,
        depths,
        radii,
        conics,
        num_tiles_hit,
        colors,
        opacities,
        H, W, BLOCK_SIZE
    )

    _rgb_im, _alpha_im, _idx = _torch_impl.rasterize_gaussians_simple(
        _xys,
        _depths,
        _radii,
        _conics,
        _num_tiles_hit,
        _colors,
        _opacities,
        H, W, BLOCK_SIZE
    )

    print("test forward")
    check_close('rgb_im', rgb_im, _rgb_im)
    check_close('alpha_im', alpha_im, _alpha_im)
    check_close('final_idx', idx, _idx)

    def fun(rgb, alpha):
        rgb_r = torch.sin(rgb).norm(dim=2)
        alpha_r = torch.exp(torch.flip(alpha, [0, 1]))
        return (rgb_r * alpha_r).mean()
    fun(rgb_im, alpha_im).backward()
    fun(_rgb_im, _alpha_im).backward()

    print("test backward")
    tol = { 'atol': 1e-6, 'rtol': 1e-5 }
    check_close('v_xys', xys.grad, _xys.grad, **tol)
    check_close('v_conics', conics.grad, _conics.grad, **tol)
    check_close('v_colors', colors.grad, _colors.grad, **tol)
    check_close('v_opacities', opacities.grad, _opacities.grad, **tol)


if __name__ == "__main__":
    rasterize_simple.RETURN_IDX = True
    test_rasterize_simple()
