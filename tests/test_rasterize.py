import pytest
import traceback
import torch
from torch.func import vjp  # type: ignore

from spirulae_splat.splat import _torch_impl
from spirulae_splat.splat.project_gaussians import project_gaussians
import spirulae_splat.splat.rasterize_simple as rasterize_simple
import spirulae_splat.splat.rasterize as rasterize
import spirulae_splat.splat.cuda as _C

torch.manual_seed(43)

device = torch.device("cuda:0")


def check_close(name, a, b, atol=1e-5, rtol=1e-5):
    print(name, [*a.shape], [*b.shape], a.dtype, b.dtype)
    try:
        torch.testing.assert_close(a, b, atol=atol, rtol=rtol)
    except AssertionError:
        traceback.print_exc()
        diff = torch.abs(a - b).detach()
        print(f"{diff.max()=} {diff.mean()=}", end='\n\n')


@pytest.mark.skipif(not torch.cuda.is_available(), reason="No CUDA device")
def test_rasterize():

    num_points, H, W = 40, 20, 30
    # num_points, H, W = 20, 10, 15
    cx, cy = W / 2, H / 2
    fx, fy = int(0.7*W), int(0.7*W)
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
    background = torch.rand(3, device=device, requires_grad=True)
    _background = background.detach().clone().requires_grad_(True)

    colors = torch.randn((num_points, 3), device=device, requires_grad=True)
    opacities = (0.995 * torch.rand((num_points, 1), device=device)).requires_grad_(True)
    anisotropies = torch.randn((num_points, 2), device=device, requires_grad=True)
    _colors = colors.detach().clone().requires_grad_(True)
    _opacities = opacities.detach().clone().requires_grad_(True)
    _anisotropies = anisotropies.detach().clone().requires_grad_(True)
    depth_ref_im = torch.randn((H, W, 3), device=device, requires_grad=True)
    _depth_normal_ref = depth_ref_im.detach().clone().requires_grad_(True)

    ch_degree_r, ch_degree_r_to_use = 3, 2
    ch_degree_phi, ch_degree_phi_to_use = 2, 1
    dim_ch = ch_degree_r * (2*ch_degree_phi+1)
    ch_coeffs = torch.randn((num_points, dim_ch, 3), device=device).requires_grad_(True)
    _ch_coeffs = ch_coeffs.detach().clone().requires_grad_(True)

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

    output = rasterize.rasterize_gaussians(
        positions,
        axes_u,
        axes_v,
        colors,
        ch_degree_r, ch_degree_r_to_use,
        ch_degree_phi, ch_degree_phi_to_use,
        ch_coeffs,
        opacities,
        anisotropies,
        depth_grads,
        depth_ref_im,
        bounds,
        num_tiles_hit,
        intrins,
        H, W, BLOCK_SIZE,
        background
    )

    rgb_im, alpha_im, idx = rasterize_simple.rasterize_gaussians_simple(
        positions,
        axes_u,
        axes_v,
        colors,
        opacities,
        anisotropies,
        bounds,
        num_tiles_hit,
        intrins,
        H, W, BLOCK_SIZE,
        background.detach().clone()
    )

    print("test consistency")
    # check_close('rgb', output[0], rgb_im)
    check_close('alpha', output[4], alpha_im)
    check_close('idx', output[5], idx)
    print()

    _output = _torch_impl.rasterize_gaussians(
        _positions,
        _axes_u,
        _axes_v,
        _colors,
        ch_degree_r, ch_degree_r_to_use,
        ch_degree_phi, ch_degree_phi_to_use,
        _ch_coeffs,
        _opacities,
        _anisotropies,
        _depth_grads,
        _depth_normal_ref,
        _bounds,
        _num_tiles_hit,
        intrins,
        H, W, BLOCK_SIZE,
        _background
    )

    print("test forward")
    check_close('out_img', output[0], _output[0])
    check_close('out_depth_grad', output[1], _output[1])
    check_close('out_reg_depth', output[2], _output[2])
    check_close('out_reg_normal', output[3], _output[3])
    check_close('out_alpha', output[4], _output[4])
    check_close('final_idx', output[5], _output[5])
    print()

    def fun(output):
        img, depth_grad, reg_depth, reg_normal, alpha, idx = output
        # reg_normal *= 0.0
        # reg_depth *= 0.0
        # depth_grad = depth_grad * torch.tensor([[1,1,1,1]]).cuda().float()
        img_r = torch.sin(img).norm(dim=2) + 1
        depth_grad_r = torch.flip(torch.sigmoid(depth_grad.norm(dim=2)), [0, 1])
        reg_r = torch.exp(-0.1*(reg_depth+1)*(reg_normal+1))
        alpha_r = torch.exp(torch.flip(alpha, [0, 1]))
        return ((img_r+1) * (alpha_r+1) + (depth_grad_r+1) * (reg_r+1)).mean()
    fun(output).backward()
    fun(_output).backward()

    print("test backward")
    tol = { 'atol': 1e-6, 'rtol': 1e-5 }
    check_close('v_positions', positions.grad, _positions.grad, **tol)
    check_close('v_axes_u', axes_u.grad, _axes_u.grad, **tol)
    check_close('v_axes_v', axes_v.grad, _axes_v.grad, **tol)
    check_close('v_colors', colors.grad, _colors.grad, **tol)
    if dim_ch > 0:
        check_close('v_ch_coeffs', ch_coeffs.grad, _ch_coeffs.grad, **tol)
    check_close('v_opacities', opacities.grad, _opacities.grad, **tol)
    check_close('v_anisotropies', anisotropies.grad, _anisotropies.grad, **tol)
    check_close('v_depth_grad', depth_grads.grad, _depth_grads.grad, **tol)
    check_close('v_depth_normal_ref', depth_ref_im.grad, _depth_normal_ref.grad, **tol)
    check_close('v_background', background.grad, _background.grad, **tol)

    assert (positions.absgrad > 0).any()
    assert (positions.absgrad >= abs(positions.grad)[:,:2]).all()
    assert (ch_coeffs.absgrad > 0).any()
    # assert (ch_coeffs.absgrad >= abs(ch_coeffs.grad)).all()


if __name__ == "__main__":
    # torch.autograd.set_detect_anomaly(True)
    rasterize_simple.RETURN_IDX = True
    rasterize.RETURN_IDX = True
    test_rasterize()
