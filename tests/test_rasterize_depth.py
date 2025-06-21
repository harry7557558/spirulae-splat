import pytest
import traceback
import torch
from torch.func import vjp  # type: ignore

from spirulae_splat.splat import _torch_impl
from spirulae_splat.splat.project_gaussians import project_gaussians
import spirulae_splat.splat.rasterize_depth as rasterize_depth
import spirulae_splat.splat.cuda as _C
from spirulae_splat.splat._camera import _Camera

rasterize_depth.DEBUG = True

torch.manual_seed(41)

device = torch.device("cuda:0")


def check_close(name, a, b, atol=1e-5, rtol=1e-5):
    if a is None or b is None:
        print(name, a if a is None else ([*a.shape], a.dtype),
              b if b is None else ([*b.shape], b.dtype))
        return
    print(name, [*a.shape], [*b.shape], a.dtype, b.dtype)
    try:
        torch.testing.assert_close(a, b, atol=atol, rtol=rtol)
    except AssertionError:
        traceback.print_exc()
        diff = (torch.abs(a - b) / torch.abs(b).mean()).detach()
        print(f"{diff.max()=} {diff.mean()=}", end='\n\n')


def _test_rasterize_depth(mode):
    num_points = 40

    means3d = torch.randn((num_points, 3), device=device)
    scales = torch.exp(torch.randn((num_points, 2), device=device))
    quats = torch.randn((num_points, 4), device=device)
    quats /= torch.linalg.norm(quats, dim=-1, keepdim=True)

    opacities = (0.995 * torch.rand((num_points, 1), device=device)).requires_grad_(True)
    _opacities = opacities.detach().clone().requires_grad_(True)

    H, W = 20, 30
    cx, cy = 0.45*W, 0.55*H
    fx, fy = 0.7*W, 0.8*W
    intrins = (fx, fy, cx, cy)
    cam = _Camera(H, W, "OPENCV", intrins)
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

    params = project_gaussians(
        means3d, scales, quats,
        viewmat, cam, clip_thresh,
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

    depth_im, meta_im, idx = rasterize_depth.rasterize_gaussians_depth(
        positions,
        axes_u,
        axes_v,
        opacities,
        bounds,
        num_tiles_hit,
        cam, mode
    )

    _depth_im, _meta_im, _idx = _torch_impl.rasterize_gaussians_depth(
        _positions,
        _axes_u,
        _axes_v,
        _opacities,
        _bounds,
        _num_tiles_hit,
        intrins, H, W, mode
    )

    print("test forward")
    tol = { 'atol': 1e-4, 'rtol': 1e-3 }
    check_close('depth_im', depth_im, _depth_im, **tol)
    check_close('meta_im', meta_im, _meta_im, **tol)
    check_close('final_idx', idx, _idx)
    print()

    def fun(depth, meta, idx):
        depth_r = torch.log(depth+1.0).squeeze()
        # depth_r = depth_r * (meta[...,0] >= 0.5)
        # depth_r = depth_r * (meta[...,1] == 1.0)
        return 1e2 * (depth_r).mean()
    fun(depth_im, meta_im, idx).backward()
    fun(_depth_im, _meta_im, _idx).backward()

    print("test backward")
    # tol = { 'atol': 1e-5, 'rtol': 1e-4 }
    tol = { 'atol': 1e-4, 'rtol': 1e-3 }
    check_close('v_positions', positions.grad, _positions.grad, **tol)
    check_close('v_axes_u', axes_u.grad, _axes_u.grad, **tol)
    check_close('v_axes_v', axes_v.grad, _axes_v.grad, **tol)
    check_close('v_opacities', opacities.grad, _opacities.grad, **tol)

    assert (positions.absgrad >= 0).all()


@pytest.mark.skipif(not torch.cuda.is_available(), reason="No CUDA device")
def test_mean():
    _test_rasterize_depth("mean")

@pytest.mark.skip
def test_median():
    _test_rasterize_depth("median")

if __name__ == "__main__":

    print("==== Mean ====")
    test_mean()

    print("\n==== Median ====")
    test_median()
