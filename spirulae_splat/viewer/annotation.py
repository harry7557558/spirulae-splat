import torch
from spirulae_splat.modules.camera import Cameras
from spirulae_splat.splat.rendering import rasterization
from spirulae_splat.viewer_legacy.utils import triangle_verts_to_quat_scale_mean

from spirulae_splat.splat.cuda import _make_lazy_cuda_func

@torch.no_grad()
def annotate_train_cameras(
    rgb: torch.Tensor,
    depths: torch.Tensor,
    alpha: torch.Tensor,
    view_camera: Cameras,
    cameras: Cameras
):
    if rgb.ndim == 4:
        rgb = rgb.squeeze(0)
    if depths.ndim == 4:
        depths = depths.squeeze(0)
    if alpha.ndim == 4:
        alpha = alpha.squeeze(0)

    R = view_camera.camera_to_worlds[:, :3, :3]  # 3 x 3
    T = view_camera.camera_to_worlds[:, :3, 3:4]  # 3 x 1
    # if self.config.relative_scale is not None:
    #     T = T * self.config.relative_scale  # TODO
    R = R * torch.tensor([[[1.0, -1.0, -1.0]]])
    R_inv = R.transpose(-1, -2)
    T_inv = -torch.bmm(R_inv, T)
    view_viewmats = torch.eye(4, dtype=view_camera.camera_to_worlds.dtype)[None].repeat(len(view_camera), 1, 1)
    view_viewmats[:, :3, :3] = R_inv
    view_viewmats[:, :3, 3:4] = T_inv

    cameras = cameras.to(rgb.device)
    R = cameras.camera_to_worlds[:, :3, :3]  # 3 x 3
    T = cameras.camera_to_worlds[:, :3, 3:4]  # 3 x 1
    # if self.config.relative_scale is not None:
    #     T = T * self.config.relative_scale  # TODO
    R = R * torch.tensor([[[1.0, -1.0, -1.0]]]).cuda()
    camera_to_worlds = torch.concat((R, T), dim=-1)

    from time import perf_counter
    torch.cuda.synchronize()
    time0 = perf_counter()
    _make_lazy_cuda_func("blit_train_cameras")(
        rgb, depths, alpha,
        view_camera.camera_type[0] == "FISHEYE",
        view_camera.intrins.cuda(),
        view_viewmats.cuda(),
        view_camera.distortion_params.cuda(),
        cameras.intrins,
        cameras.width,
        cameras.height,
        camera_to_worlds,
    )
    torch.cuda.synchronize()
    time1 = perf_counter()
    print(1e3*(time1-time0), 'ms')
    return rgb

    if rgb.ndim == 3:
        rgb = rgb[None]
    if depths.ndim == 3:
        depths = depths[None]
    if alpha.ndim == 3:
        alpha = alpha[None]

    cameras = cameras.to(rgb.device)

    R = cameras.camera_to_worlds[:, :3, :3]  # 3 x 3
    T = cameras.camera_to_worlds[:, :3, 3:4]  # 3 x 1
    # if self.config.relative_scale is not None:
    #     T = T * self.config.relative_scale
    R = R * torch.tensor([[[1.0, -1.0, -1.0]]]).cuda()

    n = len(cameras)
    o = 0.1
    ot = 0.05
    cam_verts = 0.05 * torch.Tensor([
        # [[0, 0, 0], [-1, -1, 1], [1, -1, 1]],
        # [[0, 0, 0], [1, -1, 1], [1, 1, 1]],
        # [[0, 0, 0], [1, 1, 1], [-1, 1, 1]],
        # [[0, 0, 0], [-1, 1, 1], [-1, -1, 1]],

        # [[0, 0, 0], [1-o, -1, 1], [1, -1+o, 1]],
        # [[0, 0, 0], [1, 1-o, 1], [1-o, 1, 1]],
        # [[0, 0, 0], [-1+o, 1, 1], [-1, 1-o, 1]],
        # [[0, 0, 0], [-1, -1+o, 1], [-1+o, -1, 1]],

        [[-ot, 0, 0], [1-o, -1, 1], [1, -1+o, 1]],
        [[-ot, 0, 0], [0, o, 0], [1, -1+o, 1]],
        [[0, -ot, 0], [1, 1-o, 1], [1-o, 1, 1]],
        [[0, -ot, 0], [-o, 0, 0], [1-o, 1, 1]],
        [[ot, 0, 0], [-1+o, 1, 1], [-1, 1-o, 1]],
        [[ot, 0, 0], [0, -o, 0], [-1, 1-o, 1]],
        [[0, ot, 0], [-1, -1+o, 1], [-1+o, -1, 1]],
        [[0, ot, 0], [o, 0, 0], [-1+o, -1, 1]],

        [[1-o, -1, 1], [1, -1+o, 1], [1-o, 1, 1]],
        [[1, -1+o, 1], [1-o, 1, 1], [1, 1-o, 1]],
        [[1, 1-o, 1], [1-o, 1, 1], [-1, 1-o, 1]],
        [[1-o, 1, 1], [-1, 1-o, 1], [-1+o, 1, 1]],
        [[-1+o, 1, 1], [-1, 1-o, 1], [-1+o, -1, 1]],
        [[-1, 1-o, 1], [-1+o, -1, 1], [-1, -1+o, 1]],
        [[-1, -1+o, 1], [-1+o, -1, 1], [1, -1+o, 1]],
        [[-1+o, -1, 1], [1, -1+o, 1], [1-o, -1, 1]],

    ]).float().cuda()
    sx = (cameras.width / cameras.intrins[:, 0]).reshape(-1, 1)
    sy = (cameras.height / cameras.intrins[:, 1]).reshape(-1, 1)
    # sz = (sx * sy) ** (-1/3)
    sz = (sx * sy) ** -0.2
    R[:, :, 0] *= sx * sz
    R[:, :, 1] *= sy * sz
    R[:, :, 2] *= sz
    means = torch.einsum('aij,bkj->abki', R, cam_verts) + T[:, :, None, :].transpose(-1, -3)
    means = means.reshape(-1, 3, 3)
    quats, scales, means = triangle_verts_to_quat_scale_mean(means)
    n = len(means)
    opacities = torch.ones(n, 2).to(means)
    # opacities = torch.concatenate((torch.ones(n, 1).to(means), 0.5*torch.ones(n, 1).to(means)), dim=-1)
    features_dc = -1.3 * torch.ones(n, 3).to(means)
    features_ch = torch.zeros(n, 2, 3).to(means)
    features_sh = torch.zeros(n, 15, 3).to(means)
    splat_params = (means, quats, scales, opacities, features_dc, features_sh, features_ch)

    R = view_camera.camera_to_worlds[:, :3, :3]  # 3 x 3
    T = view_camera.camera_to_worlds[:, :3, 3:4]  # 3 x 1
    # if self.config.relative_scale is not None:
    #     T = T * self.config.relative_scale
    R = R * torch.tensor([[[1.0, -1.0, -1.0]]])
    R_inv = R.transpose(-1, -2)
    T_inv = -torch.bmm(R_inv, T)
    viewmats = torch.eye(4, dtype=view_camera.camera_to_worlds.dtype)[None].repeat(len(view_camera), 1, 1)
    viewmats[:, :3, :3] = R_inv
    viewmats[:, :3, 3:4] = T_inv

    sc = 4
    while max(rgb.shape[-2], rgb.shape[-3]) * sc > 1920 and sc > 1:
        sc -= 1
    rgbd, Ts, meta = rasterization(
        # primitive="3dgut",
        primitive="opaque_triangle",
        splat_params=[x.contiguous() for x in splat_params],
        viewmats=viewmats.cuda(),
        intrins=view_camera.intrins.cuda()*sc,
        width=rgb.shape[-2]*sc,
        height=rgb.shape[-3]*sc,
        camera_model=view_camera.camera_type[0],
    )
    rgb1, depths1 = rgbd[:2]
    rgb1 = torch.nn.functional.avg_pool2d(rgb1.permute(0, 3, 1, 2), sc).permute(0, 2, 3, 1)
    depths1 = torch.nn.functional.avg_pool2d(depths1.permute(0, 3, 1, 2), sc).permute(0, 2, 3, 1)
    Ts = torch.nn.functional.avg_pool2d(Ts.permute(0, 3, 1, 2), sc).permute(0, 2, 3, 1)
    rgb1 = torch.lerp(rgb1, rgb, Ts)

    return torch.where((depths < depths1) & (alpha > 1.0-Ts), rgb, rgb1)
