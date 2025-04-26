import pytest
import traceback
import torch
from torch.func import vjp  # type: ignore

from spirulae_splat.splat import _torch_impl
from spirulae_splat.splat.project_gaussians import project_gaussians
from spirulae_splat.splat import (
    project_gaussians,
    rasterize_gaussians_indices,
    rasterize_gaussians_sorted,
    rasterize_gaussians_simplified_sorted,
    rasterize_gaussians_depth_sorted,
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
def test_rasterize_simplified_sorted():

    num_points, H, W = 40, 20, 30
    # num_points, H, W = 20, 10, 15
    # num_points, H, W = 1, 2, 3
    cx, cy = 0.45*W, 0.55*H
    fx, fy = 0.7*W, 0.8*W
    cam = _Camera(H, W, "OPENCV", (fx, fy, cx, cy))
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

    colors = torch.randn((num_points, 3), device=device, requires_grad=True)
    opacities = (0.995 * torch.rand((num_points, 1), device=device)).requires_grad_(True)
    _colors = colors.detach().clone().requires_grad_(True)
    _opacities = opacities.detach().clone().requires_grad_(True)

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

    output = rasterize_gaussians_simplified_sorted(
        positions, axes_u, axes_v, colors, opacities,
        num_intersects, sorted_indices, cam
    )

    output_r = rasterize_gaussians_sorted(
        positions, axes_u, axes_v, colors,
        0, 0, 0, 0, torch.zeros(len(positions), 0, 3).to(positions),
        opacities,
        torch.zeros((H, W, 1)).float().to(device), 1.0,
        num_intersects, sorted_indices,
        cam.intrins, H, W, cam.BLOCK_WIDTH,
    )

    output_d = rasterize_gaussians_depth_sorted(
        positions, axes_u, axes_v, opacities,
        sorted_indices,
        cam,
        "mean"
    )

    print("test consistency")
    check_close('rgb', output[0], output_r[0])
    check_close('alpha', output[1], output_r[1])
    check_close('depth', output[2], output_r[2])
    check_close('normal', output[3], output_r[3])
    check_close('reg_depth', output[4], output_r[4])
    check_close('depth1', output[2][..., 0:1]/(output[1]+1e-12), output_d, rtol=1e-3)
    print()

    _output = _torch_impl.rasterize_gaussians_simplified_sorted(
        _positions, _axes_u, _axes_v, _colors, _opacities,
        num_intersects, sorted_indices,
        cam.intrins, H, W,
    )

    print("test forward")
    check_close('out_img', output[0], _output[0])
    check_close('out_alpha', output[1], _output[1])
    check_close('out_depth', output[2], _output[2])
    check_close('out_normal', output[3], _output[3])
    check_close('out_depth_reg', output[4], _output[4])
    print()

    def fun(output):
        img, alpha, depth, normal, depth_reg = output
        img_r = torch.sin(img).norm(dim=2, keepdim=True) + 1
        normal_r = (torch.exp(normal)-0.75).norm(dim=2, keepdim=True)
        depth_r = torch.flip(torch.sigmoid(depth.norm(dim=2, keepdim=True)), [0, 1])
        reg_r = torch.exp(-0.1*(depth_reg+1))
        alpha_r = torch.exp(torch.flip(alpha, [0, 1]))
        return ((img_r+1) * (alpha_r+1) + (depth_r+1) * (reg_r+1) + normal_r).mean()
    fun(output).backward()
    fun(_output).backward()

    print("test backward")
    check_close('v_positions', positions.grad, _positions.grad)
    # print(axes_u.grad, _axes_u.grad)
    # print(axes_v.grad, _axes_v.grad)
    check_close('v_axes_u', axes_u.grad, _axes_u.grad)
    check_close('v_axes_v', axes_v.grad, _axes_v.grad)
    check_close('v_colors', colors.grad, _colors.grad)
    check_close('v_opacities', opacities.grad, _opacities.grad)

    assert (positions.absgrad > 0).any()
    assert (positions.absgrad >= abs(positions.grad)[:,:2]).all()


if __name__ == "__main__":
    # torch.autograd.set_detect_anomaly(True)
    test_rasterize_simplified_sorted()
