import pytest
import torch
from torch.func import vjp  # type: ignore

from spirulae_splat.splat._torch_impl import quat_to_rotmat

from spirulae_splat.splat.rendering import rasterization as ssplat_rasterization
from gsplat.rendering import rasterization as gsplat_rasterization

from spirulae_splat.splat.cuda import (
    ray_depth_to_linear_depth,
    svhash_create_initial_volume,
    svhash_get_voxels,
    svhash_split_voxels,
)

from utils import check_close, timeit

from typing import Literal

device = torch.device("cuda:0")

# # volume = svhash_create_initial_volume((-3, -3, -2), (3, 3, 2), (30, 30, 20))
# volume = svhash_create_initial_volume((-0.6, -0.6, -0.6), (0.6, 0.6, 0.6), (6, 6, 6))
# # volume = svhash_create_initial_volume((-0.6, -0.6, -0.6), (0.6, 0.6, 0.6), (1, 1, 1))
# # print(volume)
# voxels, voxel_indices = svhash_get_voxels(volume)
# print(voxels.shape, voxel_indices.shape)

# # densities_0 = 1.0*torch.exp(2.0*torch.randn(len(volume[0])).cuda())
# densities_0 = 1.0*torch.exp(2.0*torch.linspace(-1, 1, len(volume[0])).cuda())
# features_dc_0 = 0.5+1.0*voxels[:, :3]
# features_sh_0 = torch.zeros((len(features_dc_0), 0, 3)).cuda()

# torch.random.manual_seed(42)
# mask = torch.rand(len(voxels)).to(voxels.device) < 0.2
# print(mask)
# # print(volume)
# volume, cell_idx, vert_idx, vert_weight = svhash_split_voxels(volume, mask)
# # print(cell_idx)
# # print(vert_idx)
# # print(vert_weight)
# # print(vert_weight.sum(-1))
# # print(densities_0)
# densities_0 = torch.cat((
#     densities_0,
#     (densities_0[vert_idx.flatten()].reshape(*vert_idx.shape, *densities_0.shape[1:]) * vert_weight).sum(1)
#     # densities_0[vert_idx[:,0]]
# ), 0)
# # print(densities_0)
# # print(features_dc_0.shape)
# features_dc_0 = torch.cat((features_dc_0, features_dc_0[cell_idx]), 0)
# features_sh_0 = torch.cat((features_sh_0, features_sh_0[cell_idx]), 0)
# # # print(volume)
# voxels, voxel_indices = svhash_get_voxels(volume)
# # print(voxel_indices)
# # print(voxels.shape, voxel_indices.shape)

# # print(voxels)
# # exit(0)


B, W, H = 4, 1440, 1080
N, SH_DEGREE = 200000, 3
PACKED = False
IS_FISHEYE = False
IS_ANTIALIASED = False
WITH_UT = True

def rasterize_ssplat(means, quats, scales, opacities, features_dc, features_sh, viewmats, Ks):
    camera_model = ["pinhole", "fisheye"][IS_FISHEYE]
    # quats = torch.nn.functional.normalize(quats, dim=-1)
    intrins = torch.stack([Ks[..., 0, 0], Ks[..., 1, 1], Ks[..., 0, 2], Ks[..., 1, 2]], dim=-1)  # fx, fy, cx, cy
    rgbd, alpha, meta = ssplat_rasterization(
        primitive="3dgut" if WITH_UT else ["3dgs", "mip"][IS_ANTIALIASED],
        splat_params=(means, quats, scales, opacities.unsqueeze(-1), features_dc, features_sh),
        # primitive="opaque_triangle",
        # splat_params=(means, quats, scales+1.8, opacities.unsqueeze(-1).repeat(1, 2), features_dc, features_sh, features_dc.unsqueeze(-2).repeat(1, 2, 1)),
        # primitive="voxel",
        # splat_params=(torch.cat((means, 5.0*torch.exp(scales.mean(-1, True))), dim=-1), 2.0*torch.exp(opacities).unsqueeze(-1).repeat(1, 8), features_dc, features_sh),
        # splat_params=(voxels, 10.0*torch.exp(opacities)[voxel_indices], features_dc[:len(voxels)], features_sh[:len(voxels)]),
        # splat_params=(voxels, densities_0[voxel_indices], features_dc_0, features_sh_0),
        viewmats=viewmats,  # [C, 4, 4]
        intrins=intrins,  # [C, 4]
        width=W,
        height=H,
        packed=PACKED,
        use_bvh=False,
        sparse_grad=False,
        distributed=False,
        camera_model=camera_model,
        render_mode="RGB+D",
        # render_mode="RGB+D+N",
    )
    rgbd = [*rgbd[:2]]
    rgbd[1] = ray_depth_to_linear_depth(rgbd[1], camera_model, intrins)  # TODO: f(E[X]) != E[f(X)]
    return *rgbd, alpha

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
        sh_degree=SH_DEGREE,
        packed=PACKED,
        absgrad=False,
        sparse_grad=False,
        rasterize_mode=["classic", "antialiased"][IS_ANTIALIASED],
        distributed=False,
        camera_model=["pinhole", "fisheye"][IS_FISHEYE],
        with_ut=WITH_UT,
        with_eval3d=WITH_UT,
        render_mode="RGB+ED",
    )
    return rgbd[..., :3], rgbd[..., 3:], alpha


def get_inputs():
    torch.manual_seed(44)

    means = torch.randn((N, 3)).to(device)
    quats = torch.randn((N, 4)).to(device)
    scales = torch.randn((N, 3)).to(device)*0.7 - 5.0
    # means[0] -= 1e4
    opacities = torch.randn((N,)).to(device)
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

    outputs = rasterize_ssplat(*inputs)
    _outputs = rasterize_gsplat(*_inputs)

    print("test forward")
    tol = { 'atol': 1e-3, 'rtol': 1e-3 }
    check_close('rgb', outputs[0], _outputs[0], **tol)
    check_close('depth', outputs[1], _outputs[1], **tol)
    check_close('alpha', outputs[2], _outputs[2], **tol)
    print()

    if False:
        import matplotlib.pyplot as plt
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
        rgb = outputs[0].detach().cpu().numpy()
        _rgb = _outputs[0].detach().cpu().numpy()
        alpha = outputs[2].detach().cpu().numpy()
        _alpha = _outputs[2].detach().cpu().numpy()
        ax1.imshow(rgb[0])
        ax2.imshow(alpha[0])
        # ax2.imshow(rgb[1])
        # ax3.imshow(rgb[2])
        # ax4.imshow(rgb[3])
        ax3.imshow(_rgb[0])
        ax4.imshow(_alpha[0])
        # ax4.imshow(_rgb[1])
        # plt.show()
        plt.savefig("/mnt/d/plot.png")
        exit(0)

    weights = [torch.randn_like(x.detach()) for x in _outputs]
    weights[1] *= 0.0  # depth
    def fun(outputs):
        return sum([(w*o).sum() for w, o in zip(weights, outputs)])
    fun(outputs).backward()
    fun(_outputs).backward()

    # print(inputs[0].grad)
    # print(_inputs[0].grad)

    print("test backward")
    tol = { 'atol': 1e-2, 'rtol': 1e-2 }
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
