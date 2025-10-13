import pytest
import torch
from torch.func import vjp  # type: ignore

from spirulae_splat.splat._torch_impl import quat_to_rotmat

from spirulae_splat.splat.rendering import rasterization as ssplat_rasterization
from gsplat.rendering import rasterization as gsplat_rasterization

from utils import check_close, timeit

device = torch.device("cuda:0")


B, W, H = 4, 1440, 1080
N, SH_DEGREE = 200000, 3
PACKED = False
IS_FISHEYE = False

def rasterize_ssplat(means, quats, scales, opacities, features_dc, features_sh, viewmats, Ks):
    rgbd, alpha, meta = ssplat_rasterization(
        means=means,
        quats=quats,
        scales=torch.exp(scales),
        opacities=opacities.squeeze(-1),
        colors_dc=features_dc,
        colors_sh=features_sh,
        viewmats=viewmats,  # [C, 4, 4]
        Ks=Ks,  # [C, 3, 3]
        width=W,
        height=H,
        sh_degree=SH_DEGREE,
        packed=PACKED,
        use_bvh=False,
        tile_size=16,
        absgrad=False,
        sparse_grad=False,
        rasterize_mode="classic",
        distributed=False,
        camera_model=["pinhole", "fisheye"][IS_FISHEYE],
        with_ut=IS_FISHEYE,
        with_eval3d=IS_FISHEYE,
        render_mode="RGB+D",
    )
    return rgbd[..., :3], rgbd[..., 3:], alpha

def rasterize_gsplat(means, quats, scales, opacities, features_dc, features_sh, viewmats, Ks):
    rgbd, alpha, meta = gsplat_rasterization(
        means=means,
        quats=quats,
        scales=torch.exp(scales),
        opacities=opacities.squeeze(-1),
        colors=torch.concatenate([features_dc.unsqueeze(1), features_sh], dim=1),  # TODO: slow
        viewmats=viewmats,  # [C, 4, 4]
        Ks=Ks,  # [C, 3, 3]
        width=W,
        height=H,
        sh_degree=SH_DEGREE,
        packed=PACKED,
        tile_size=16,
        absgrad=False,
        sparse_grad=False,
        rasterize_mode="classic",
        distributed=False,
        camera_model=["pinhole", "fisheye"][IS_FISHEYE],
        with_ut=IS_FISHEYE,
        with_eval3d=IS_FISHEYE,
        render_mode="RGB+D",
    )
    return rgbd[..., :3], rgbd[..., 3:], alpha


def get_inputs():
    torch.manual_seed(42)

    means = torch.randn((N, 3)).to(device)
    quats = torch.randn((N, 4)).to(device)
    scales = torch.randn((N, 3)).to(device)*0.7 - 5.0
    # means[0] -= 1e4
    opacities = torch.rand((N,)).to(device)
    features_dc = torch.rand((N, 3)).to(device)
    features_sh = 0.2 * torch.randn((N, (SH_DEGREE+1)**2-1, 3)).to(device)

    viewmats = torch.eye(4)[None].repeat(B, 1, 1)
    viewmats[:, :3, 3] = 1.0 * torch.randn((B, 3))
    viewmats[:, 2, 3] = torch.abs(viewmats[:, 2, 3])
    rotation = torch.randn(B, 4)
    rotation[:, 1:] *= 0.1
    rotation /= torch.norm(rotation, dim=-1, keepdim=True)
    viewmats[:, :3, :3] = quat_to_rotmat(rotation)
    viewmats = viewmats.contiguous().to(device)

    cx, cy = 0.5*W, 0.5*H
    fx, fy = 0.4*W, 0.4*W
    Ks = torch.tensor([[[fx, 0, cx], [0, fy, cy], [0, 0, 1]]]).repeat(B, 1, 1).to(device)
    Ks *= torch.exp(0.2*torch.randn_like(Ks))

    inputs = (means, quats, scales, opacities, features_dc, features_sh, viewmats, Ks)
    inputs = [torch.nn.Parameter(x.contiguous().clone()) for x in inputs]
    return inputs



@pytest.mark.skipif(not torch.cuda.is_available(), reason="No CUDA device")
def test_rasterization():

    inputs = get_inputs()
    _inputs = get_inputs()

    outputs = rasterize_ssplat(*inputs)
    _outputs = rasterize_gsplat(*_inputs)

    print("test forward")
    check_close('rgb', outputs[0], _outputs[0])
    check_close('depth', outputs[1], _outputs[1])
    check_close('alpha', outputs[2], _outputs[2])
    print()

    if False:
        import matplotlib.pyplot as plt
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
        rgb = outputs[0].detach().cpu().numpy()
        _rgb = _outputs[0].detach().cpu().numpy()
        ax1.imshow(rgb[0])
        # ax2.imshow(rgb[1])
        # ax3.imshow(rgb[2])
        # ax4.imshow(rgb[3])
        ax3.imshow(_rgb[0])
        # ax4.imshow(_rgb[1])
        plt.show()

    weights = [torch.randn_like(x.detach()) for x in _outputs]
    def fun(outputs):
        return sum([(w*o).sum() for w, o in zip(weights, outputs)])
    fun(outputs).backward()
    fun(_outputs).backward()

    # print(inputs[0].grad)
    # print(_inputs[0].grad)

    print("test backward")
    tol = { 'atol': 1e-3, 'rtol': 1e-4 }
    check_close('means', inputs[0].grad, _inputs[0].grad, **tol)
    check_close('quats', inputs[1].grad, _inputs[1].grad, **tol)
    check_close('scales', inputs[2].grad, _inputs[2].grad, **tol)
    check_close('opacs', inputs[3].grad, _inputs[3].grad, **tol)
    check_close('features_dc', inputs[4].grad, _inputs[4].grad, **tol)
    check_close('features_sh', inputs[5].grad, _inputs[5].grad, **tol)
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
    outputs = rasterize_ssplat(*inputs)
    _outputs = rasterize_gsplat(*_inputs)

    weights = [torch.randn_like(x.detach()) for x in _outputs]
    def fun(outputs):
        return sum([(w*o).sum() for w, o in zip(weights, outputs)])
    loss = fun(outputs)
    _loss = fun(_outputs)

    timeit(lambda: loss.backward(retain_graph=True), "ssplat backward", repeat=20)
    timeit(lambda: _loss.backward(retain_graph=True), "gsplat backward")
    print()



if __name__ == "__main__":

    N = 1000
    test_rasterization()
    print()

    N = 200000
    profile_rasterization()
