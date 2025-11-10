import numpy as np
import torch
import torch.nn.functional as F
import os
import json

from spirulae_splat.splat.rendering import rasterization
from spirulae_splat.splat.background_sh import render_background_sh
from spirulae_splat.splat.sh import spherical_harmonics
from spirulae_splat.splat._torch_impl import quat_mult
from spirulae_splat.splat.cuda import ray_depth_to_linear_depth

from spirulae_splat.viewer.utils import (
    quat_to_rotmat,
    rotmat_to_quat,
    quat_scale_to_triangle_verts,
    triangle_verts_to_quat_scale_mean,
)

from spirulae_splat.viewer.camera import Camera

from spirulae_splat.perf_timer import PerfTimer
timer = PerfTimer("render")

from typing import Literal


class SplatModel:
    def __init__(self, file_path: str):
        self.bgr = True
        self.flip_yz = False
        self.return_torch = False
        
        self.primitive = "3dgs"  # type: Literal["3dgs", "mip", "opaque_triangle"]

        self.dataparser_transform = np.eye(4)
        self.dataparser_scale = 1.0

        if os.path.isdir(file_path):
            for fn in os.listdir(file_path):
                if fn.endswith('.ckpt') or fn == 'config.yml':
                    file_path = os.path.join(file_path, fn)
                    break
        if file_path.endswith('.ckpt'):
            self.load_ckpt(file_path)
        elif file_path.endswith('config.yml'):
            self.load_config(file_path)
        else:
            raise ValueError("Must be .ckpt or config.yml")

    def load_ckpt(self, file_path, _from_load_config=False):
        if not _from_load_config:
            print("WARNING: ckpt file does not contain information about scene type. Assuming vanilla 3DGS. "
                "Use config.yml file instead of .ckpt file for more information.")
        checkpoint = torch.load(file_path, 'cpu', weights_only=False)
        pipeline = checkpoint['pipeline']

        self.gauss_params = {}
        for key in pipeline:
            if key.startswith("_model.gauss_params."):
                value = pipeline[key].to(torch.float32).cuda()
                key = key.split('.')[-1]
                self.gauss_params[key] = value
        self.background_color = pipeline['_model.background_color'].to(torch.float32).cuda()
        self.background_sh = pipeline['_model.background_sh'].to(torch.float32).cuda() \
            if '_model.background_sh' in pipeline else torch.zeros((0, 3))

        self.num_sh = self.features_sh.shape[1]
        self.sh_degree = {
            1: 0, 4: 1, 9: 2, 16: 3, 25: 4
        }[self.num_sh+1]
        self.background_sh_degree = {
            1: 0, 4: 1, 9: 2, 16: 3, 25: 4
        }[len(self.background_sh)+1]

    def load_config(self, file_path: str):
        save_dir = file_path[:file_path.rfind(os.path.sep)]
        ckpt_dir = os.path.join(save_dir, 'nerfstudio_models')
        for f in os.listdir(ckpt_dir):
            if f.endswith('.ckpt'):
                f = os.path.join(ckpt_dir, f)
                self.load_ckpt(f, _from_load_config=True)
                break

        content = open(file_path).read()
        if "primitive: mip" in content:
            self.primitive = "mip"
        elif "primitive: opaque_triangle" in content:
            self.primitive = "opaque_triangle"
        print("Primitive:", self.primitive)

        # load dataparser transforms
        dtr_path = os.path.join(save_dir, 'dataparser_transforms.json')
        with open(dtr_path, 'r') as fp:
            dtr = json.load(fp)
        self.dataparser_transform = np.concatenate((dtr['transform'], [[0, 0, 0, 1]]))
        self.dataparser_scale = dtr['scale']

    @torch.no_grad
    def change_frame(self, rot, tr, sc):
        from scipy.spatial.transform import Rotation
        from spirulae_splat.viewer.utils import rotate_sh_coeffs

        if isinstance(rot, np.ndarray):
            rot = torch.from_numpy(rot).to(self.means)
        if isinstance(tr, np.ndarray):
            tr = torch.from_numpy(tr).to(self.means)
        if isinstance(sc, np.ndarray):
            sc = torch.from_numpy(sc).to(self.means)

        self.gauss_params["means"] = (self.means * sc @ rot.T + tr)

        self.gauss_params["scales"] = (self.scales + np.log(sc))

        dq = rotmat_to_quat(rot)
        self.gauss_params["quats"] = quat_mult(dq, self.quats)

        self.gauss_params["features_sh"] = rotate_sh_coeffs(self.features_sh, rot, "gsplat")
        if self.background_sh_degree > 0:
            self.background_sh = rotate_sh_coeffs(self.background_sh[None], rot, "nerfstudio")[0]

    @torch.no_grad()
    def convert_to_input_frame(self, match: Literal['ply', 'json', None]=None):
        """Convert to the same coordinate frame as in input dataset"""

        if match == "ply":
            # TODO: load this dynamically from transforms.json
            # It's not always this for some datasets
            applied = np.array([
                [1, 0, 0, 0],
                [0, 0, 1, 0],
                [0, -1, 0, 0],
                [0, 0, 0, 1]
            ])  # to match sparse_pc.ply
        elif match == "json":
            applied = np.diag([1, -1, -1, 1])  # to match transforms.json
        else:
            applied = np.eye(4)

        transform = self.dataparser_transform
        transform = applied @ np.linalg.inv(transform)
        transform = torch.from_numpy(transform).to(self.means)
        rot = transform[:3, :3]
        tr = transform[:3, 3]
        sc = 1.0 / self.dataparser_scale

        self.change_frame(rot, tr, sc)

    def num_splats(self):
        return len(self.gauss_params["means"])

    @property
    def colors(self):
        if self.sh_degree > 0:
            return 0.28209479177387814*self.features_dc+0.5
        else:
            return self.features_dc
    @property
    def num_points(self):
        return self.means.shape[0]
    @property
    def means(self):
        return self.gauss_params["means"]
    @property
    def scales(self):
        return self.gauss_params["scales"]
    @property
    def quats(self):
        return self.gauss_params["quats"]
    @property
    def features_dc(self):
        return self.gauss_params["features_dc"]
    @property
    def features_sh(self):
        return self.gauss_params["features_sh"]
    @property
    def features_ch(self):
        return self.gauss_params["features_ch"]
    @property
    def opacities(self):
        return self.gauss_params["opacities"]

    def _get_background_image(self, camera: Camera, c2w: np.ndarray):
        """TODO: support distortion"""
        sh_degree = self.background_sh_degree
        if not (sh_degree > 0):
            return self.background_color.repeat(camera.h, camera.w, 1)
        sh_coeffs = torch.cat((self.background_color.unsqueeze(0), self.background_sh), dim=0)  # [(deg+1)^2, 3]

        if not self.flip_yz:
            c2w = c2w[..., :3, :3] @ np.diag([1, -1, -1])
        c2w = torch.from_numpy(c2w).to(sh_coeffs)

        fx, fy, cx, cy = camera.intrins
        Ks = torch.tensor([[fx, 0, cx], [0, fy, cy], [0, 0, 1]]).float().to(c2w)
        return render_background_sh(
            camera.w, camera.h,
            ['pinhole', 'fisheye'][camera.model == "OPENCV_FISHEYE"],
            Ks, c2w, sh_degree, sh_coeffs
        )[0]

    @torch.inference_mode()
    def _render(self, camera: Camera, c2w: np.ndarray, return_depth=False):

        timer.start()

        c2w = c2w.astype(np.float32)
        R = c2w[:3, :3]  # 3 x 3
        T = c2w[:3, 3:4]  # 3 x 1
        if self.flip_yz:
            R_edit = np.diag([1, -1, -1])
            R = R @ R_edit
        R_inv = R.T
        T_inv = -R_inv @ T
        viewmat = np.eye(4, dtype=np.float32)
        viewmat[:3, :3] = R_inv
        viewmat[:3, 3:4] = T_inv
        viewmat = torch.from_numpy(viewmat).to(self.means.device)
        ssplat_camera = camera._to_ssplat_camera()
        timer.mark("pre")

        if self.scales.shape[-1] == 3:
            fx, fy, cx, cy = ssplat_camera.intrins
            w, h = int(ssplat_camera.w), int(ssplat_camera.h)
            Ks = torch.tensor([[fx, 0, cx], [0, fy, cy], [0, 0, 1]]).float().to(viewmat)

            kwargs = {}
            dist_coeffs = torch.tensor([*ssplat_camera.dist_coeffs]).float().to(viewmat)
            is_fisheye = ssplat_camera.model == "OPENCV_FISHEYE"
            is_distorted = ssplat_camera.is_distorted() or any([x != 0 for x in dist_coeffs])
            if is_distorted:
                # TODO
                raise NotImplementedError()
                if is_fisheye:
                    kwargs['radial_coeffs'] = dist_coeffs[None]
                else:
                    kwargs["radial_coeffs"] = torch.concat((dist_coeffs[:2], 0.0*dist_coeffs[:2]))[None]
                    kwargs["tangential_coeffs"] = dist_coeffs[2:][None]

            rgbd, alpha, meta = rasterization(
                self.primitive,
                (self.means, F.normalize(self.quats, dim=-1), self.scales, self.opacities.squeeze(-1),
                    self.features_dc, self.features_sh)
                    if self.primitive in ['3dgs', 'mip'] else
                (self.means, F.normalize(self.quats, dim=-1), self.scales,
                    torch.concat([
                        torch.ones_like(self.opacities),
                        torch.ones_like(self.opacities)
                    ], dim=-1),
                    self.features_dc, self.features_sh, self.features_ch
                ),
                viewmats=viewmat[None].contiguous(),  # [C, 4, 4]
                Ks=Ks[None].contiguous(),  # [C, 3, 3]
                width=w,
                height=h,
                packed=False,
                use_bvh=False,
                absgrad=False,
                sparse_grad=False,
                distributed=False,
                camera_model=["pinhole", "fisheye"][is_fisheye],
                with_ut=True,
                with_eval3d=True,
                render_mode="RGB+ED",
                # render_mode="RGB+ED+N",
                **kwargs,
            )
            alpha = alpha[0]

            rgb = rgbd[0]
            if return_depth:
                depth = rgbd[1]
                depth = ray_depth_to_linear_depth(
                    depth,
                    ["pinhole", "fisheye"][is_fisheye],
                    Ks[None].contiguous(),
                    **kwargs
                )
                depth = torch.where(
                    alpha > 0.0, depth,
                    1.5*torch.amax(depth).detach()
                ).contiguous()

        else:
            raise NotImplementedError("2DGS is deprecated")

        # print(torch.mean(rgb[...,0]).item()*255.0, torch.mean(rgb[...,1]).item()*255.0, torch.mean(rgb[...,2]).item())
        # print(torch.amax(rgb, dim=(0,1))*255)

        # blend with background
        background = self.background_color.reshape((1, 1, 3))
        if self.background_sh_degree > 0 or True:
            background = self._get_background_image(camera, c2w)
            rgb = rgb + (1.0 - alpha) * background

        rgb = torch.clip(rgb, 0.0, 1.0)
        timer.mark("background")

        if return_depth:
            return rgb[0], depth[0]
        return rgb[0]

    @torch.no_grad()
    def render(self, camera, c2w, return_depth=False):
        c2w = camera.c2w_to_c2o(c2w)

        if return_depth:
            rgb, depth = self._render(camera, c2w, True)

            if self.return_torch:
                return rgb, depth

            im = rgb.cpu().numpy()
            im = (im*255).astype(np.uint8)
            depth = depth.cpu().numpy()
            return im, depth

        rgb = self._render(camera, c2w)
        if self.return_torch:
            return rgb

        im = rgb.cpu().numpy()
        im = (im*255).astype(np.uint8)
        return im
