import pytest
import traceback
import torch
from torch.func import vjp  # type: ignore

from spirulae_splat.splat import _torch_impl
from spirulae_splat.splat import (
    project_gaussians,
    rasterize_gaussians_indices,
    rasterize_gaussians_sorted,
    rasterize_gaussians_simple_sorted,
)
import spirulae_splat.splat.cuda as _C
from spirulae_splat.splat._camera import _Camera

torch.manual_seed(41)

device = torch.device("cuda:0")


def check_close(name, a, b, atol=5e-6, rtol=1e-4):
    print(name, [*a.shape], [*b.shape], a.dtype, b.dtype)
    try:
        torch.testing.assert_close(a, b, atol=atol, rtol=rtol)
    except AssertionError:
        traceback.print_exc()
        diff = (torch.abs(a - b) / torch.abs(b).mean()).detach()
        print(f"{diff.max()=} {diff.mean()=}", end='\n\n')


@pytest.mark.skipif(not torch.cuda.is_available(), reason="No CUDA device")
def test_rasterize_sorted():

    num_points, H, W = 40, 20, 30
    # num_points, H, W = 20, 10, 15
    cx, cy = 0.45*W, 0.55*H
    fx, fy = 0.7*W, 0.8*W
    intrins = (fx, fy, cx, cy)
    cam = _Camera(H, W, "OPENCV", intrins)
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

    means3d = 0.8*torch.randn((num_points, 3), device=device)
    scales = 0.8*torch.exp(torch.randn((num_points, 2), device=device))
    quats = torch.randn((num_points, 4), device=device)
    quats /= torch.linalg.norm(quats, dim=-1, keepdim=True)
    # background = torch.rand(3, device=device, requires_grad=True)
    background = torch.zeros(3, device=device, requires_grad=True)
    _background = background.detach().clone().requires_grad_(True)

    colors = torch.randn((num_points, 3), device=device, requires_grad=True)
    opacities = (0.995 * torch.rand((num_points, 1), device=device)).requires_grad_(True)
    _colors = colors.detach().clone().requires_grad_(True)
    _opacities = opacities.detach().clone().requires_grad_(True)
    depth_ref_im = torch.randn((H, W, 1), device=device, requires_grad=True)
    _depth_normal_ref = depth_ref_im.detach().clone().requires_grad_(True)

    ch_degree_r, ch_degree_r_to_use = 3, 2
    ch_degree_phi, ch_degree_phi_to_use = 2, 1
    dim_ch = ch_degree_r * (2*ch_degree_phi+1)
    ch_coeffs = torch.randn((num_points, dim_ch, 3), device=device).requires_grad_(True)
    _ch_coeffs = ch_coeffs.detach().clone().requires_grad_(True)

    params = project_gaussians(
        means3d, scales, quats,
        viewmat, cam, clip_thresh,
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
        bounds,
        num_tiles_hit,
    ) = decode_params(params)
    (
        _positions,
        _axes_u,
        _axes_v,
        _bounds,
        _num_tiles_hit,
    ) = decode_params(params)

    num_intersects, sorted_indices = rasterize_gaussians_indices(
        positions, axes_u, axes_v,
        opacities, bounds, num_tiles_hit,
        cam
    )

    depth_reg_pairwise_factor = 0.7

    output = rasterize_gaussians_sorted(
        positions, axes_u, axes_v, colors,
        ch_degree_r, ch_degree_r_to_use,
        ch_degree_phi, ch_degree_phi_to_use,
        ch_coeffs,
        opacities,
        depth_ref_im, depth_reg_pairwise_factor,
        num_intersects, sorted_indices,
        cam
    )

    rgb_im, alpha_im = rasterize_gaussians_simple_sorted(
        positions, axes_u, axes_v, colors, opacities,
        num_intersects, sorted_indices,
        cam, background.detach().clone()
    )

    print("test consistency")
    # check_close('rgb', output[0], rgb_im)
    check_close('alpha', output[1], alpha_im)
    print()

    _output = _torch_impl.rasterize_gaussians_sorted(
        _positions, _axes_u, _axes_v, _colors,
        ch_degree_r, ch_degree_r_to_use,
        ch_degree_phi, ch_degree_phi_to_use,
        _ch_coeffs,
        _opacities,
        _depth_normal_ref, depth_reg_pairwise_factor,
        num_intersects, sorted_indices,
        intrins, H, W,
    )

    print("test forward")
    check_close('out_img', output[0], _output[0])
    check_close('out_alpha', output[1], _output[1])
    check_close('out_depth', output[2], _output[2])
    check_close('out_normal', output[3], _output[3])
    check_close('out_depth', output[4], _output[4])
    print()

    def fun(output):
        img, alpha, depth, normal, reg_depth = output
        img_r = torch.sin(img).norm(dim=2, keepdim=True) + 1
        alpha_r = torch.exp(torch.flip(alpha, [0, 1]))
        depth_r = torch.flip(torch.sin(0.5*depth), [0, 1])
        normal_r = torch.flip(torch.sigmoid(normal).norm(dim=2, keepdim=True), [0, 1])
        reg_r = torch.exp(-0.1*(reg_depth+1))
        return ((img_r+1) * (alpha_r+1) + (depth_r+1) + (normal_r+1) * (reg_r+1)).mean()
    fun(output).backward()
    fun(_output).backward()

    print("test backward")
    check_close('v_positions', positions.grad, _positions.grad)
    check_close('v_axes_u', axes_u.grad, _axes_u.grad)
    check_close('v_axes_v', axes_v.grad, _axes_v.grad)
    check_close('v_colors', colors.grad, _colors.grad)
    if dim_ch > 0:
        check_close('v_ch_coeffs', ch_coeffs.grad, _ch_coeffs.grad)
    check_close('v_opacities', opacities.grad, _opacities.grad)
    check_close('v_depth_normal_ref', depth_ref_im.grad, _depth_normal_ref.grad)
    # check_close('v_background', background.grad, _background.grad)

    assert (positions.absgrad > 0).any()
    assert (positions.absgrad >= abs(positions.grad)[:,:2]).all()
    # assert (ch_coeffs.absgrad > 0).any()
    # assert (ch_coeffs.absgrad >= abs(ch_coeffs.grad)).all()


if __name__ == "__main__":
    test_rasterize_sorted()
