import torch
import numpy as np
from spirulae_splat.modules.camera import Cameras, CameraType, _COLMAP_CAMERA_MODEL_TO_TYPE

from spirulae_splat.splat.cuda import _C


def knn_dist(x: torch.Tensor, k: int = 4):
    if len(x) == 1:
        return 1.0
    if len(x) <= k:
        dists = torch.cdist(x, x).flatten()
        return torch.median(dists[dists > 0.0])
    x_np = x.detach().cpu().numpy()

    # TODO: fix cameras at same location, e.g. pinhole converted from equirectangular

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

    if cameras.camera_type[0] == CameraType.EQUIRECTANGULAR.value:
        from spirulae_splat.modules.resample import warp_equirectangular_to_pinhole
        camera_to_worlds, intrins, distortion_params, thumbnails = \
             warp_equirectangular_to_pinhole(cameras.camera_to_worlds.cuda(), thumbnails.cuda())
        cameras = Cameras(
            intrins, distortion_params,
            thumbnails.shape[-3], thumbnails.shape[-2],
            camera_to_worlds, "PINHOLE"
        )
    elif kwargs.get("warp_to_pinhole", False):
        from spirulae_splat.modules.resample import warp_wide_to_pinhole
        intrins = cameras.intrins.clone()
        intrins[:, 0] *= thumbnails.shape[-2] / cameras.width
        intrins[:, 1] *= thumbnails.shape[-3] / cameras.height
        intrins[:, 2] *= thumbnails.shape[-2] / cameras.width
        intrins[:, 3] *= thumbnails.shape[-3] / cameras.height
        camera_to_worlds, intrins, distortion_params, thumbnails = \
             warp_wide_to_pinhole(
                cameras.camera_type[0],  # TODO: mixed fisheye/pinhole
                intrins.cuda(), cameras.distortion_params.cuda(), cameras.camera_to_worlds.cuda(),
                thumbnails.cuda()
            )
        cameras = Cameras(
            intrins, distortion_params,
            thumbnails.shape[-3], thumbnails.shape[-2],
            camera_to_worlds, "PINHOLE"
        )

    R = view_camera.camera_to_worlds[:, :3, :3]  # 3 x 3
    T = view_camera.camera_to_worlds[:, :3, 3:4]  # 3 x 1
    # T = T * relative_scale
    R = R * torch.tensor([[[1.0, -1.0, -1.0]]]).to(R)
    R_inv = R.transpose(-1, -2)
    T_inv = -torch.bmm(R_inv, T)
    view_viewmats = torch.eye(4, dtype=view_camera.camera_to_worlds.dtype)[None].repeat(len(view_camera), 1, 1)
    view_viewmats[:, :3, :3] = R_inv
    view_viewmats[:, :3, 3:4] = T_inv

    cameras = cameras.to("cuda")
    R = cameras.camera_to_worlds[:, :3, :3]  # 3 x 3
    T = cameras.camera_to_worlds[:, :3, 3:4]  # 3 x 1
    # T = T * relative_scale
    R = R * torch.tensor([[[1.0, -1.0, -1.0]]]).cuda()
    camera_to_worlds = torch.concat((R, T), dim=-1)

    # must match Common.cuh and projection_utils.slang
    camera_model_mapper = {
        "PINHOLE": 0,
        "FISHEYE": 1,
        "EQUISOLID": 2,
    }

    key = "_camera_models"
    if not hasattr(cameras, key):
        camera_models = torch.Tensor([camera_model_mapper[ctype.upper()] for ctype in cameras.camera_type]).int().cuda()
        setattr(cameras, key, camera_models)
    else:
        camera_models = getattr(cameras, key)

    # Kernel expects CUDA inputs. Coerce any CPU tensor to CUDA, then keep a
    # reference alive in _keepalive so its data_ptr stays valid through the C++
    # call (otherwise Python may free a temporary between _tv() and the kernel,
    # leading to non-deterministic crashes / garbled outputs).
    _keepalive = []
    def _tv(t):
        if t is None:
            return (0, 4, [0])
        if not t.is_cuda:
            t = t.cuda()
        t = t.contiguous()
        _keepalive.append(t)
        return (t.data_ptr(), t.element_size(), list(t.shape))

    h, w = rgb.shape[0], rgb.shape[1]
    # Output goes to CUDA (kernel writes); caller can move to CPU as needed.
    out_device = rgb.device if rgb.is_cuda else torch.device("cuda")
    out_rgb = torch.empty(h, w, 3, dtype=torch.uint8, device=out_device)

    _C.blit_train_cameras(
        _tv(rgb), _tv(depths), _tv(alpha),
        camera_model_mapper[view_camera.camera_type[0]],
        _tv(view_camera.intrins.cuda()),
        _tv(view_viewmats.cuda()),
        _tv(view_camera.distortion_params.cuda()),
        _tv(cameras.intrins),
        _tv(cameras.width),
        _tv(cameras.height),
        _tv(camera_models),
        _tv(cameras.distortion_params),
        _tv(camera_to_worlds),
        _tv(thumbnails.cuda()),
        size,
        kwargs.get("show_training_cameras", False),
        _tv(out_rgb),
    )
    return out_rgb
