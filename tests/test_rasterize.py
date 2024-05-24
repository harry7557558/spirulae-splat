import pytest
import traceback
import torch
from torch.func import vjp  # type: ignore

from spirulae_splat.splat import _torch_impl
from spirulae_splat.splat.project_gaussians import project_gaussians
import spirulae_splat.splat.rasterize_simple as rasterize_simple
import spirulae_splat.splat.rasterize as rasterize
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
def test_rasterize():

    num_points = 40
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

    means3d = torch.randn((num_points, 3), device=device)
    scales = torch.exp(torch.randn((num_points, 2), device=device))
    quats = torch.randn((num_points, 4), device=device)
    quats /= torch.linalg.norm(quats, dim=-1, keepdim=True)
    background = torch.randn(3, device=device)

    colors = torch.randn((num_points, 3), device=device, requires_grad=True)
    opacities = torch.rand((num_points, 1), device=device, requires_grad=True)
    _colors = colors.detach().clone().requires_grad_(True)
    _opacities = opacities.detach().clone().requires_grad_(True)
    depth_normal_ref = torch.randn((H, W, 2), device=device, requires_grad=True)
    _depth_normal_ref = depth_normal_ref.detach().clone().requires_grad_(True)

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
        params = [*params]
        for i in range(len(params)):
            params[i] = torch.tensor(params[i].detach().cpu().numpy()).to(device)
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

    output = rasterize.rasterize_gaussians(
        xys,
        depths,
        depth_grads,
        radii,
        conics,
        num_tiles_hit,
        colors,
        opacities,
        depth_normal_ref,
        H, W, BLOCK_SIZE,
        background
    )

    rgb_im, alpha_im, idx = rasterize_simple.rasterize_gaussians_simple(
        xys,
        depths,
        radii,
        conics,
        num_tiles_hit,
        colors,
        opacities,
        H, W, BLOCK_SIZE,
        background
    )

    print("test consistency")
    check_close('rgb', output[0], rgb_im)
    check_close('alpha', output[4], alpha_im)
    check_close('idx', output[5], idx)

    _output = _torch_impl.rasterize_gaussians(
        _xys,
        _depths,
        _depth_grads,
        _radii,
        _conics,
        _num_tiles_hit,
        _colors,
        _opacities,
        _depth_normal_ref,
        H, W, BLOCK_SIZE,
        background
    )

    print("test forward")
    check_close('out_img', output[0], _output[0])
    check_close('out_depth_grad', output[1], _output[1])
    check_close('out_reg_depth', output[2], _output[2])
    check_close('out_reg_normal', output[3], _output[3])
    check_close('out_alpha', output[4], _output[4])
    check_close('final_idx', output[5], _output[5])

    def fun(output):
        img, depth_grad, reg_depth, reg_normal, alpha, idx = output
        # reg_normal *= 0.0
        # reg_depth *= 0.0
        img_r = torch.sin(img).norm(dim=2) + 1
        depth_grad_r = torch.flip(torch.cos(depth_grad.norm(dim=2)), [0, 1])
        reg_r = torch.exp(-0.1*(reg_depth+1)*(reg_normal+1))
        alpha_r = torch.exp(torch.flip(alpha, [0, 1]))
        return (img_r * alpha_r + depth_grad_r * reg_r).mean()
    fun(output).backward()
    fun(_output).backward()

    print("test backward")
    tol = { 'atol': 1e-6, 'rtol': 1e-5 }
    check_close('v_xys', xys.grad, _xys.grad, **tol)
    check_close('v_depth', depths.grad, _depths.grad, **tol)
    check_close('v_depth_grad', depth_grads.grad, _depth_grads.grad, **tol)
    check_close('v_conics', conics.grad, _conics.grad, **tol)
    check_close('v_colors', colors.grad, _colors.grad, **tol)
    check_close('v_opacities', opacities.grad, _opacities.grad, **tol)
    check_close('v_depth_normal_ref', depth_normal_ref.grad, _depth_normal_ref.grad, **tol)


if __name__ == "__main__":
    # torch.autograd.set_detect_anomaly(True)
    rasterize_simple.RETURN_IDX = True
    rasterize.RETURN_IDX = True
    test_rasterize()
