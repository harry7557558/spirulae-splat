import numpy as np
import torch
import torch.nn.functional as F
import os

from spirulae_splat.splat import (
    project_gaussians,
    rasterize_gaussians_simple,
    rasterize_gaussians_depth,
    rasterize_gaussians_indices,
    rasterize_gaussians_simple_sorted,
    rasterize_gaussians_depth_sorted,
)
from spirulae_splat.splat.background_sh import render_background_sh
from spirulae_splat.splat.sh import spherical_harmonics

from spirulae_splat.viewer.camera import Camera

from spirulae_splat.perf_timer import PerfTimer
timer = PerfTimer("render")


class SplatModel:
    def __init__(self, file_path: str):
        if file_path.endswith('.ckpt'):
            self.load_ckpt(file_path)
        elif file_path.endswith('config.yml'):
            self.load_config(file_path)
        self.bgr = True
        self.flip_yz = False
        self.return_torch = False

    def load_ckpt(self, file_path):
        checkpoint = torch.load(file_path)
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
        self.num_ch = self.features_ch.shape[1]
        assert self.num_ch == 0, "Not implemented"
        self.sh_degree = {
            1: 0, 4: 1, 9: 2, 16: 3, 25: 4
        }[self.num_sh+1]
        self.ch_degree_r, self.ch_degree_phi = {
            0: (0, 0), 1: (1, 0), 2: (2, 0),
            3: (1, 1), 5: (1, 2), 6: (2, 1),
            7: (1, 3), 9: (3, 1), 10: (2, 2),
            14: (2, 3), 15: (3, 2), 21: (3, 3)
        }[self.num_ch]
        self.background_sh_degree = {
            1: 0, 4: 1, 9: 2, 16: 3, 25: 4
        }[len(self.background_sh)+1]

    def load_config(self, file_path: str):
        save_dir = file_path[:file_path.rfind(os.path.sep)]
        ckpt_dir = os.path.join(save_dir, 'nerfstudio_models')
        for f in os.listdir(ckpt_dir):
            if f.endswith('.ckpt'):
                f = os.path.join(ckpt_dir, f)
                self.load_ckpt(f)
                break

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
        viewmat = torch.from_numpy(c2w[:3, :3] * np.array([1,-1,-1], dtype=np.float32)).cuda()
        return render_background_sh(camera._to_ssplat_camera(), viewmat, sh_degree+1, sh_coeffs)

    @torch.inference_mode()
    def _render(self, camera: Camera, c2w: np.ndarray, return_depth=False, sort_per_pixel=True):

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
        ssplat_camera = camera._to_ssplat_camera()
        timer.mark("pre")

        quats_norm = F.normalize(self.quats, dim=-1)
        if self.sh_degree > 0:
            viewdirs = self.means.detach() - torch.from_numpy(c2w[:3, 3]).cuda()  # (N, 3)
            n = self.sh_degree
            rgbs = spherical_harmonics(n, viewdirs, self.features_dc, self.features_sh)
            rgbs = torch.relu(rgbs+0.5)  # type: ignore
        else:
            rgbs = self.features_dc
        if self.bgr:
            rgbs = rgbs[:, [2, 1, 0]]
        opacities = torch.sigmoid(self.opacities)
        timer.mark("map")

        (
            positions, axes_u, axes_v,
            bounds, num_tiles_hit
        ) = project_gaussians(  # type: ignore
            self.means,
            torch.exp(self.scales),
            quats_norm,
            torch.from_numpy(viewmat[:3, :]).cuda(),
            ssplat_camera,
        )  # type: ignore
        timer.mark("project")

        if sort_per_pixel:
            num_intersects, sorted_indices = rasterize_gaussians_indices(
                positions, axes_u, axes_v, opacities,
                bounds, num_tiles_hit,
                ssplat_camera
            )
            timer.mark("sort")

            rgb, alpha = rasterize_gaussians_simple_sorted(
                positions, axes_u, axes_v, rgbs, opacities,
                num_intersects, sorted_indices,
                ssplat_camera
            )
            timer.mark("rasterize")

        else:
            rgb, alpha = rasterize_gaussians_simple(
                positions,
                axes_u, axes_v,
                rgbs, opacities,
                bounds, num_tiles_hit, ssplat_camera,
            )
            timer.mark("rasterize")

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
            if sort_per_pixel:
                depth = rasterize_gaussians_depth_sorted(
                    positions, axes_u, axes_v, opacities,
                    sorted_indices,
                    ssplat_camera,
                    "median"
                )
            else:
                depth = rasterize_gaussians_depth(
                    positions, axes_u, axes_v, opacities,
                    bounds, num_tiles_hit,
                    ssplat_camera,
                    "median"
                )
            depth = torch.where(depth == 0.0, depth.max(), depth)
            timer.end("depth")
            return rgb, depth

        timer.end("depth")
        return rgb

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
