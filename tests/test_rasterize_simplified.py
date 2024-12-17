import pytest
import traceback
import torch
from torch.func import vjp  # type: ignore

from spirulae_splat.splat import _torch_impl
from spirulae_splat.splat.project_gaussians import project_gaussians
from spirulae_splat.splat import rasterize_simplified, rasterize, rasterize_depth
import spirulae_splat.splat.cuda as _C

torch.manual_seed(43)

device = torch.device("cuda:0")


def check_close(name, a, b, atol=1e-5, rtol=1e-5):
    print(name, [*a.shape], [*b.shape], a.dtype, b.dtype)
    try:
        torch.testing.assert_close(a, b, atol=atol, rtol=rtol)
    except AssertionError:
        traceback.print_exc()
        diff = (torch.abs(a - b) / torch.abs(b).mean()).detach()
        print(f"{diff.max()=} {diff.mean()=}", end='\n\n')


@pytest.mark.skipif(not torch.cuda.is_available(), reason="No CUDA device")
def test_rasterize_simplified():

    num_points, H, W = 40, 20, 30
    # num_points, H, W = 20, 10, 15
    cx, cy = 0.45*W, 0.55*H
    fx, fy = 0.7*W, 0.8*W
    intrins = (fx, fy, cx, cy)
    clip_thresh = 0.01
    viewmat = torch.tensor(
        [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 6.0],
            [0.0, 0.0, 0.0, 1.0],
        ],
        device=device,
    )
    viewmat[:3, :3] = _torch_impl.quat_to_rotmat(torch.randn(4))
    BLOCK_SIZE = 16

    means3d = 0.8*torch.randn((num_points, 3), device=device)
    scales = 0.8*torch.exp(torch.randn((num_points, 2), device=device))
    quats = torch.randn((num_points, 4), device=device)
    quats /= torch.linalg.norm(quats, dim=-1, keepdim=True)

    colors = torch.randn((num_points, 3), device=device, requires_grad=True)
    opacities = (0.995 * torch.rand((num_points, 1), device=device)).requires_grad_(True)
    anisotropies = torch.randn((num_points, 2), device=device, requires_grad=True)
    _colors = colors.detach().clone().requires_grad_(True)
    _opacities = opacities.detach().clone().requires_grad_(True)
    _anisotropies = anisotropies.detach().clone().requires_grad_(True)

    params = project_gaussians(
        means3d, scales, quats,
        viewmat, intrins,
        H, W, BLOCK_SIZE,
        clip_thresh,
    )
    def decode_params(params):
        params = [*params]
        for i in range(len(params)):
            params[i] = torch.tensor(params[i].detach().cpu().numpy()).to(device)
            try:
                params[i] = params[i].requires_grad_(True)
            except RuntimeError:
                pass
        return params
    (
        positions,
        axes_u,
        axes_v,
        depth_grads,
        bounds,
        num_tiles_hit,
    ) = decode_params(params)
    (
        _positions,
        _axes_u,
        _axes_v,
        _depth_grads,
        _bounds,
        _num_tiles_hit,
    ) = decode_params(params)

    output = rasterize_simplified.rasterize_gaussians_simplified(
        positions,
        axes_u,
        axes_v,
        colors,
        opacities,
        anisotropies,
        depth_grads,
        bounds,
        num_tiles_hit,
        intrins,
        H, W, BLOCK_SIZE,
    )

    output_r = rasterize.rasterize_gaussians(
        positions, axes_u, axes_v,
        colors,
        0, 0, 0, 0, torch.zeros(len(positions), 0, 3).to(positions),
        opacities, anisotropies, depth_grads,
        torch.zeros((H, W, 3)).float().to(device),
        1.0,
        bounds, num_tiles_hit,
        intrins, H, W, BLOCK_SIZE,
    )

    output_d = rasterize_depth.rasterize_gaussians_depth(
        positions, axes_u, axes_v,
        opacities, anisotropies,
        bounds, num_tiles_hit,
        intrins, H, W, BLOCK_SIZE,
        "mean"
    )

    print("test consistency")
    check_close('rgb', output[0], output_r[0])
    check_close('depth', output[2], output_r[1])
    check_close('reg_depth', output[3], output_r[2])
    check_close('alpha', output[1], output_r[4])
    check_close('idx', output[4], output_r[5])
    check_close('depth1', output[2][..., 2:3]/(output[1]+1e-6), output_d)
    print()

    _output = _torch_impl.rasterize_gaussians_simplified(
        _positions,
        _axes_u,
        _axes_v,
        _colors,
        _opacities,
        _anisotropies,
        _depth_grads,
        _bounds,
        _num_tiles_hit,
        intrins,
        H, W, BLOCK_SIZE,
    )

    print("test forward")
    check_close('out_img', output[0], _output[0])
    check_close('out_alpha', output[1], _output[1])
    check_close('out_depth', output[2], _output[2])
    check_close('out_depth_reg', output[3], _output[3])
    check_close('final_idx', output[4], _output[4])
    print()

    def fun(output):
        img, alpha, depth, depth_reg, idx = output
        img_r = torch.sin(img).norm(dim=2, keepdim=True) + 1
        depth_r = torch.flip(torch.sigmoid(depth.norm(dim=2, keepdim=True)), [0, 1])
        reg_r = torch.exp(-0.1*(depth_reg+1))
        alpha_r = torch.exp(torch.flip(alpha, [0, 1]))
        return ((img_r+1) * (alpha_r+1) + (depth_r+1) * (reg_r+1)).mean()
    fun(output).backward()
    fun(_output).backward()

    print("test backward")
    tol = { 'atol': 1e-6, 'rtol': 1e-5 }
    check_close('v_positions', positions.grad, _positions.grad, **tol)
    check_close('v_axes_u', axes_u.grad, _axes_u.grad, **tol)
    check_close('v_axes_v', axes_v.grad, _axes_v.grad, **tol)
    check_close('v_colors', colors.grad, _colors.grad, **tol)
    check_close('v_opacities', opacities.grad, _opacities.grad, **tol)
    check_close('v_anisotropies', anisotropies.grad, _anisotropies.grad, **tol)
    check_close('v_depth_grad', depth_grads.grad, _depth_grads.grad, **tol)

    assert (positions.absgrad > 0).any()
    assert (positions.absgrad >= abs(positions.grad)[:,:2]).all()


if __name__ == "__main__":
    # torch.autograd.set_detect_anomaly(True)
    rasterize.RETURN_IDX = True
    rasterize_simplified.RETURN_IDX = True
    test_rasterize_simplified()
