import torch
from torch import Tensor
import math
import numpy as np

from spirulae_splat.splat.cuda import _make_lazy_cuda_func

from functools import cache

@cache
def get_cubemap_faces_equirectangular():
    return torch.Tensor([
        [[1,0,0],[0,1,0],[0,0,1]],
        [[0,0,-1],[0,1,0],[1,0,0]],
        [[-1,0,0],[0,0,1],[0,1,0]],
        [[0,-1,0],[0,0,1],[-1,0,0]],
        [[1,0,0],[0,0,1],[0,-1,0]],
        [[1,0,0],[0,-1,0],[0,0,-1]],
    ]).float().cuda()

@torch.no_grad()
def warp_equirectangular_to_pinhole(
    camera_to_worlds: Tensor,  # [B, 3, 4],
    images: Tensor,  # [B, H, W, C]
):
    axes = get_cubemap_faces_equirectangular()
    b, h, w, c = images.shape
    out_shape = int(math.ceil(math.sqrt(h * w / len(axes))))
    images = _make_lazy_cuda_func("warp_image_equirectangular_to_pinhole")(
        images, axes, out_shape, out_shape
    )
    # print(torch.linalg.det(axes))
    # if b == 1:
    #     import matplotlib.pyplot as plt
    #     plt.imshow(images[0].reshape(-1, out_shape, c).detach().cpu().numpy())
    #     plt.show()

    # TODO: fused kernels
    camera_to_worlds = camera_to_worlds.clone()
    camera_to_worlds[:, :3, :3] = camera_to_worlds[:, :3, :3] * torch.tensor([[[1.0, -1.0, -1.0]]]).to(camera_to_worlds)
    camera_to_worlds = torch.cat((
        torch.einsum('cij,bkj->bcki', axes, camera_to_worlds[:, :3, :3]),
        camera_to_worlds[:, None, :3, 3:].repeat(1, len(axes), 1, 1)
    ), dim=-1)
    camera_to_worlds[:, :, :3, :3] = camera_to_worlds[:, :, :3, :3] * torch.tensor([[[[1.0, -1.0, -1.0]]]]).to(camera_to_worlds)
    intrins = torch.from_numpy(np.array([
        [out_shape/2, out_shape/2, out_shape/2, out_shape/2]
    ]*(b*len(axes)), dtype=np.float32))
    distortion_params = torch.zeros(b*len(axes), 10).float()
    return (
        camera_to_worlds.view(-1, *camera_to_worlds.shape[2:]),
        intrins,
        distortion_params,
        images.view(-1, *images.shape[2:]),
    )


@cache
def get_cubemap_faces_fisheye5():
    r2, r3, r6 = 2**0.5, 3**0.5, 6**0.5
    a0 = [1/r2, 1/r6, 1/r3]
    a1 = [-1/r2, 1/r6, 1/r3]
    a2 = [0, -2/r6, 1/r3]
    return torch.Tensor([
        # [a0, a1, a2],
        # [a1, a2, a0],
        # [a2, a0, a1],
        [[1,0,0],[0,1,0],[0,0,1]],
        [[0,1,0],[0,0,1],[1,0,0]],
        [[-1,0,0],[0,0,1],[0,1,0]],
        [[0,-1,0],[0,0,1],[-1,0,0]],
        [[1,0,0],[0,0,1],[0,-1,0]],
    ]).float().cuda()

@torch.no_grad()
def warp_wide_to_pinhole(
    camera_model: str,
    intrins: Tensor,  # [B, 4]
    dist_coeffs: Tensor,  # [B, 4]
    camera_to_worlds: Tensor,  # [B, 3, 4]
    images: Tensor,  # [B, H, W, C]
):
    if intrins.ndim == 1:
        intrins = intrins[None]
    if dist_coeffs.ndim == 1:
        dist_coeffs = dist_coeffs[None]
    if camera_to_worlds.ndim == 2:
        camera_to_worlds = camera_to_worlds[None]

    axes = get_cubemap_faces_fisheye5()
    b, h, w, c = images.shape
    out_shape = int(math.ceil(math.sqrt(h * w / len(axes))))
    images = _make_lazy_cuda_func("warp_image_wide_to_pinhole")(
        camera_model, intrins, dist_coeffs, images, axes, out_shape, out_shape
    )
    # print(torch.linalg.det(axes))
    # if b == 1:
    #     import matplotlib.pyplot as plt
    #     plt.imshow(images[0].reshape(-1, out_shape, c).detach().cpu().numpy())
    #     plt.show()

    # TODO: fused kernels
    camera_to_worlds = camera_to_worlds.clone()
    camera_to_worlds[:, :3, :3] = camera_to_worlds[:, :3, :3] * torch.tensor([[[1.0, -1.0, -1.0]]]).to(camera_to_worlds)
    camera_to_worlds = torch.cat((
        torch.einsum('cij,bkj->bcki', axes, camera_to_worlds[:, :3, :3]),
        camera_to_worlds[:, None, :3, 3:].repeat(1, len(axes), 1, 1)
    ), dim=-1)
    camera_to_worlds[:, :, :3, :3] = camera_to_worlds[:, :, :3, :3] * torch.tensor([[[[1.0, -1.0, -1.0]]]]).to(camera_to_worlds)
    intrins = torch.from_numpy(np.array([
        [out_shape/2, out_shape/2, out_shape/2, out_shape/2]
    ]*(b*len(axes)), dtype=np.float32))
    distortion_params = torch.zeros(b*len(axes), 10).float()
    return (
        camera_to_worlds.view(-1, *camera_to_worlds.shape[2:]),
        intrins,
        distortion_params,
        images.view(-1, *images.shape[2:]),
    )
