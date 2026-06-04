import pytest
import torch
from torch.func import vjp  # type: ignore

from spirulae_splat.splat._torch_impl import quat_to_rotmat

from spirulae_splat.modules.core import Renderer

# commit 48df61993fed8e4742d0758b6c8c2b4e599c4124
from gsplat.rendering import rasterization as gsplat_rasterization

from spirulae_splat.splat.cuda import (
    ray_depth_to_linear_depth,
)

from utils import check_close, timeit

from typing import Literal

device = torch.device("cuda:0")


B, W, H = 4, 1440, 1080
N, SH_DEGREE, SH_DEGREE_TO_USE = 200000, 3, 3
PACKED = True
IS_FISHEYE = False
IS_ANTIALIASED = False
WITH_UT = False

def rasterize_ssplat(means, quats, scales, opacities, features_dc, features_sh, viewmats, Ks):
    renderer = Renderer(
        primitive="3dgut" if WITH_UT else ["3dgs", "mip"][IS_ANTIALIASED],
        splats_world=(means, quats, scales, opacities.unsqueeze(-1), features_dc, features_sh),
        cur_num_splats=len(means),
        packed=PACKED,
    )
    camera_model = ["pinhole", "fisheye"][IS_FISHEYE]
    # quats = torch.nn.functional.normalize(quats, dim=-1)  # affects gradient
    intrins = torch.stack([Ks[..., 0, 0], Ks[..., 1, 1], Ks[..., 0, 2], Ks[..., 1, 2]], dim=-1)  # fx, fy, cx, cy
    renderer.set_params(
        viewmats=viewmats,  # [C, 4, 4]
        intrins=intrins,  # [C, 4]
        width=W,
        height=H,
        sh_degree_to_use=SH_DEGREE_TO_USE,
        camera_model=camera_model,
    )
    renderer.forward()
    rgbd = renderer.render_colors
    Ts = renderer.render_Ts
    rgbd = [*rgbd[:2]]
    if WITH_UT:
        try:
            rgbd[1] = ray_depth_to_linear_depth(rgbd[1], camera_model, intrins)  # TODO: f(E[X]) != E[f(X)]
        except AttributeError:
            pass  # ray_depth_to_linear_depth_forward not compiled, depth comparison will differ
    return (*rgbd, 1.0 - Ts), renderer

def rasterize_gsplat(means, quats, scales, opacities, features_dc, features_sh, viewmats, Ks):
    rgbd, alpha, meta = gsplat_rasterization(
        means=means,
        quats=quats,
        scales=torch.exp(scales),
        opacities=torch.sigmoid(opacities).squeeze(-1),
        colors=torch.concatenate([features_dc.unsqueeze(1), features_sh], dim=1),  # TODO: slow
        viewmats=viewmats,  # [C, 4, 4]
        Ks=Ks,  # [C, 3, 3]
        width=W,
        height=H,
        sh_degree=SH_DEGREE_TO_USE,
        packed=PACKED and not WITH_UT,
        sparse_grad=False,
        rasterize_mode=["classic", "antialiased"][IS_ANTIALIASED],
        distributed=False,
        camera_model=["pinhole", "fisheye"][IS_FISHEYE],
        with_ut=WITH_UT,
        with_eval3d=WITH_UT,
        render_mode="RGB+ED",
    )
    rgb, depth = rgbd[..., :3], rgbd[..., 3:]
    # depth = torch.where(alpha == 0.0, depth.amax(), depth)
    return rgb, depth, alpha


def get_inputs(seed=44):
    torch.manual_seed(seed)

    means = torch.randn((N, 3), device=device)
    quats = torch.randn((N, 4), device=device)
    scales = torch.randn((N, 3), device=device)*0.7 - 5.0
    # means[0] -= 1e4
    opacities = torch.randn((N,), device=device)
    features_dc = torch.rand((N, 3), device=device)
    features_sh = 0.2 * torch.randn((N, (SH_DEGREE+1)**2-1, 3), device=device)

    quats = torch.nn.functional.normalize(quats)

    viewmats = torch.eye(4)[None].repeat(B, 1, 1)
    viewmats[:, :3, 3] = 1.0 * torch.randn((B, 3))
    viewmats[:, 2, 3] = torch.abs(viewmats[:, 2, 3])
    rotation = torch.randn(B, 4)
    rotation[:, 1:] *= 0.1
    rotation /= torch.norm(rotation, dim=-1, keepdim=True)
    viewmats[:, :3, :3] = quat_to_rotmat(rotation)
    viewmats = viewmats.contiguous().to(device)

    cx, cy = 0.5*W, 0.5*H
    fx, fy = 0.4*W, 0.4*W  # 0.4
    Ks = torch.tensor([[[fx, 0, cx], [0, fy, cy], [0, 0, 1]]]).repeat(B, 1, 1).to(device)
    Ks *= torch.exp(0.2*torch.randn_like(Ks))

    inputs = (means, quats, scales, opacities, features_dc, features_sh, viewmats, Ks)
    inputs = [torch.nn.Parameter(x.contiguous().clone()) for x in inputs]
    if WITH_UT:
        inputs = inputs[:-2] + [viewmats, Ks]
    return inputs



@pytest.mark.skipif(not torch.cuda.is_available(), reason="No CUDA device")
def test_rasterization():

    inputs = get_inputs()
    _inputs = get_inputs()

    outputs, renderer = rasterize_ssplat(*inputs)
    _outputs = rasterize_gsplat(*_inputs)

    print("test forward")
    tol = { 'atol': 1e-3, 'rtol': 1e-3 }
    check_close('rgb', outputs[0], _outputs[0], **tol)
    check_close('depth', outputs[1], _outputs[1], **tol)
    check_close('alpha', outputs[2], _outputs[2], **tol)
    print()

    if False:
        import matplotlib.pyplot as plt
        import numpy as np
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
        rgb = outputs[0].detach().cpu().numpy()
        _rgb = _outputs[0].detach().cpu().numpy()
        depth = outputs[1].detach().cpu().numpy()
        _depth = _outputs[1].detach().cpu().numpy()
        alpha = outputs[2].detach().cpu().numpy()
        _alpha = _outputs[2].detach().cpu().numpy()
        ax1.imshow(rgb[0])
        # ax2.imshow(alpha[0])
        # ax2.imshow(depth[0])
        ax2.imshow(np.log10(np.abs(depth[0]-_depth[0])))
        # ax2.imshow(rgb[1])
        # ax3.imshow(rgb[2])
        # ax4.imshow(rgb[3])
        ax3.imshow(_rgb[0])
        # ax4.imshow(_alpha[0])
        ax4.imshow(_depth[0])
        # ax4.imshow(_rgb[1])
        plt.show()
        # plt.savefig("/mnt/d/plot.png")
        exit(0)

    weights = [torch.randn_like(x.detach()) for x in _outputs]
    # weights[1] *= 0.0  # depth
    def fun(outputs):
        return sum([(w*o).sum() for w, o in zip(weights, outputs)])
    fun(_outputs).backward()

    renderer.zero_grad()
    renderer.backward((*weights[:-1], None), -weights[-1])
    grads = renderer.v_splats_world

    # print(inputs[0].grad)
    # print(_inputs[0].grad)

    print("test backward")
    tol = { 'atol': 5e-2, 'rtol': 5e-2 }
    check_close('means', grads[0], _inputs[0].grad, **tol)
    check_close('quats', grads[1], _inputs[1].grad, **tol)
    check_close('scales', grads[2], _inputs[2].grad, **tol)
    check_close('opacs', grads[3].squeeze(-1), _inputs[3].grad, **tol)
    check_close('features_dc', grads[4], _inputs[4].grad, **tol)
    check_close('features_sh', grads[5].reshape(N, -1, 3), _inputs[5].grad, **tol)
    check_close('viewmats', inputs[6].grad, _inputs[6].grad, **tol)
    print()

@pytest.mark.skipif(not torch.cuda.is_available(), reason="No CUDA device")
def profile_rasterization():

    inputs = get_inputs()
    _inputs = get_inputs()

    print("profile forward")
    timeit(lambda: rasterize_ssplat(*inputs), "ssplat forward")
    timeit(lambda: rasterize_gsplat(*_inputs), "gsplat forward")
    print()

    print("profile backward")
    outputs, renderer = rasterize_ssplat(*inputs)
    _outputs = rasterize_gsplat(*_inputs)

    renderer.zero_grad()

    weights = [torch.randn_like(x.detach()) for x in _outputs]
    def fun(outputs):
        return sum([(w*o).sum() for w, o in zip(weights, outputs)])
    loss = fun(outputs)
    _loss = fun(_outputs)

    timeit(lambda: renderer.backward((*weights[:-1], None), -weights[-1]), "ssplat backward", repeat=20)
    timeit(lambda: _loss.backward(retain_graph=True), "gsplat backward")
    print()



if __name__ == "__main__":
    import itertools

    N = 1000

    # packed x fisheye x (not_aa+no_ut, aa+no_ut, not_aa+with_ut)
    modes = [(False, False), (True, False), (False, True)]
    # modes = [(False, True)]
    combos = list(itertools.product([True, False], [False, True], modes))
    # modes = [(False, False), (True, False)]
    # combos = list(itertools.product([True, False], [False], modes))

    for packed_val, fisheye_val, (aa_val, ut_val) in combos:
        PACKED = packed_val
        IS_FISHEYE = fisheye_val
        IS_ANTIALIASED = aa_val
        WITH_UT = ut_val
        label = f"packed={PACKED} fisheye={IS_FISHEYE} antialiased={IS_ANTIALIASED} with_ut={WITH_UT}"
        print(f"=== {label} ===")
        try:
            test_rasterization()
        except Exception as e:
            import traceback
            traceback.print_exc()
        print()

    # N = 200000
    # profile_rasterization()
