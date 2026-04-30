import torch
import numpy as np
from spirulae_splat.modules.camera import Cameras
from spirulae_splat.splat.rendering import rasterization
from spirulae_splat.viewer_legacy.utils import triangle_verts_to_quat_scale_mean

from spirulae_splat.splat.cuda import _make_lazy_cuda_func


def knn_dist(x: torch.Tensor, k: int = 4):
    if len(x) == 1:
        return 1.0
    if len(x) <= k:
        dists = torch.cdist(x, x).flatten()
        return torch.median(dists[dists > 0.0])
    x_np = x.detach().cpu().numpy()

    from scipy.spatial import cKDTree
    tree = cKDTree(x_np, balanced_tree=False)

    distances, indices = tree.query(x_np, k=k+1, workers=-1)

    # return distances[:, 1:].mean()
    # return np.median(distances[:, 1:], axis=-1).mean()
    return np.median(np.amax(distances[:, 1:], axis=-1))


@torch.no_grad()
def annotate_train_cameras(
    rgb: torch.Tensor,
    depths: torch.Tensor,
    alpha: torch.Tensor,
    view_camera: Cameras,
    cameras: Cameras,
    thumbnails: torch.Tensor,
    **kwargs
):
    if rgb.ndim == 4:
        rgb = rgb.squeeze(0)
    if depths.ndim == 4:
        depths = depths.squeeze(0)
    if alpha.ndim == 4:
        alpha = alpha.squeeze(0)

    key = '_annotation_size'
    if not hasattr(cameras, key):
        T = cameras.camera_to_worlds[:, :3, 3:4]
        size = 0.2 * knn_dist(T.squeeze(-1))
        setattr(cameras, key, size)
    else:
        size = getattr(cameras, key)
    # size = size * (kwargs.get("relative_scale", 1.0) or 1.0)
    # relative_scale = size / 0.05
    # print(relative_scale)

    R = view_camera.camera_to_worlds[:, :3, :3]  # 3 x 3
    T = view_camera.camera_to_worlds[:, :3, 3:4]  # 3 x 1
    # T = T * relative_scale
    R = R * torch.tensor([[[1.0, -1.0, -1.0]]])
    R_inv = R.transpose(-1, -2)
    T_inv = -torch.bmm(R_inv, T)
    view_viewmats = torch.eye(4, dtype=view_camera.camera_to_worlds.dtype)[None].repeat(len(view_camera), 1, 1)
    view_viewmats[:, :3, :3] = R_inv
    view_viewmats[:, :3, 3:4] = T_inv

    cameras = cameras.to(rgb.device)
    R = cameras.camera_to_worlds[:, :3, :3]  # 3 x 3
    T = cameras.camera_to_worlds[:, :3, 3:4]  # 3 x 1
    # T = T * relative_scale
    R = R * torch.tensor([[[1.0, -1.0, -1.0]]]).cuda()
    camera_to_worlds = torch.concat((R, T), dim=-1)

    key = "_is_fisheyes"
    if not hasattr(cameras, key):
        is_fisheyes = torch.Tensor([(ctype == "FISHEYE") for ctype in cameras.camera_type]).bool().cuda()
        setattr(cameras, key, is_fisheyes)
    else:
        is_fisheyes = getattr(cameras, key)

    from time import perf_counter
    torch.cuda.synchronize()
    time0 = perf_counter()
    # depths = depths * relative_scale
    rgb = _make_lazy_cuda_func("blit_train_cameras")(
        rgb, depths, alpha,
        view_camera.camera_type[0] == "FISHEYE",
        view_camera.intrins.cuda(),
        view_viewmats.cuda(),
        view_camera.distortion_params.cuda(),
        cameras.intrins,
        cameras.width,
        cameras.height,
        is_fisheyes,
        cameras.distortion_params,
        camera_to_worlds,
        thumbnails.cuda(),
        size,
        kwargs.get("show_training_cameras", False),
    )
    torch.cuda.synchronize()
    time1 = perf_counter()
    # print(1e3*(time1-time0), 'ms')
    return rgb
