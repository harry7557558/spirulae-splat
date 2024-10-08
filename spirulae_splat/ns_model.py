# ruff: noqa: E741
# Copyright 2022 the Regents of the University of California, Nerfstudio Team and contributors. All rights reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
NeRF implementation that combines many recent advancements.
"""

import math
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Type, Union, Optional
from typing_extensions import Literal

import numpy as np
from scipy.spatial.transform import Rotation
import torch
from torch.nn import Parameter
from pytorch_msssim import SSIM

from spirulae_splat.splat._torch_impl import quat_to_rotmat
from spirulae_splat.splat.project_gaussians import project_gaussians
from spirulae_splat.splat.rasterize import rasterize_gaussians
from spirulae_splat.splat.rasterize_depth import rasterize_gaussians_depth
from spirulae_splat.splat.rasterize_simple import rasterize_gaussians_simple
from spirulae_splat.splat.sh import num_sh_bases, spherical_harmonics
from spirulae_splat.splat.relocation import compute_relocation, compute_relocation_split

from nerfstudio.cameras.camera_optimizers import (CameraOptimizer,
                                                  CameraOptimizerConfig)
from nerfstudio.cameras.cameras import Cameras
from nerfstudio.data.scene_box import OrientedBox
from nerfstudio.engine.callbacks import (TrainingCallback,
                                         TrainingCallbackAttributes,
                                         TrainingCallbackLocation)
from nerfstudio.engine.optimizers import Optimizers
# need following import for background color override
from nerfstudio.model_components import renderers
from nerfstudio.models.base_model import Model, ModelConfig
from nerfstudio.utils.colors import get_color
from nerfstudio.utils.rich_utils import CONSOLE
from nerfstudio.utils.math import components_from_spherical_harmonics


def random_quat_tensor(N):
    """
    Defines a random quaternion tensor of shape (N, 4)
    """
    u = torch.rand(N)
    v = torch.rand(N)
    w = torch.rand(N)
    return torch.stack(
        [
            torch.sqrt(1 - u) * torch.sin(2 * math.pi * v),
            torch.sqrt(1 - u) * torch.cos(2 * math.pi * v),
            torch.sqrt(u) * torch.sin(2 * math.pi * w),
            torch.sqrt(u) * torch.cos(2 * math.pi * w),
        ],
        dim=-1,
    )


def RGB2SH(rgb):
    """
    Converts from RGB values [0,1] to the 0th spherical harmonic coefficient
    """
    C0 = 0.28209479177387814
    return (rgb - 0.5) / C0


def SH2RGB(sh):
    """
    Converts from the 0th spherical harmonic coefficient to RGB values [0,1]
    """
    C0 = 0.28209479177387814
    return sh * C0 + 0.5


def resize_image(image: torch.Tensor, d: int):
    """
    Downscale images using the same 'area' method in opencv

    :param image shape [H, W, C]
    :param d downscale factor (must be 2, 4, 8, etc.)

    return downscaled image in shape [H//d, W//d, C]
    """
    import torch.nn.functional as tf

    weight = (1.0 / (d * d)) * torch.ones((1, 1, d, d), dtype=torch.float32, device=image.device)
    return tf.conv2d(image.permute(2, 0, 1)[:, None, ...], weight, stride=d).squeeze(1).permute(1, 2, 0)


class SaturateKeepGradient(torch.autograd.Function):
    @staticmethod
    def forward(ctx, x):
        return torch.clip(x, 0.0, 1.0)
    @staticmethod
    def backward(ctx, v_x):
        return v_x

def saturate_keep_gradient(x):
    return SaturateKeepGradient.apply(x)


@dataclass
class SpirulaeModelConfig(ModelConfig):

    _target: Type = field(default_factory=lambda: SpirulaeModel)

    warmup_length: int = 200
    """period of steps where refinement is turned off"""
    refine_every: int = 100
    """period of steps where gaussians are culled and densified"""
    resolution_schedule: int = 3000
    """training starts at 1/d resolution, every n steps this is doubled"""
    background_color: Literal["random", "black", "white"] = "white"
    """Whether to randomize the background color."""
    num_downscales: int = 2
    """at the beginning, resolution is 1/2^d, where d is this number"""
    use_mcmc: bool = False
    """use Markov-Chain Monte Carlo for gaussian control"""
    random_init: bool = False
    """whether to initialize the positions uniformly randomly (not SFM points)"""
    num_random: int = 20000
    """Number of gaussians to initialize if random init is used"""
    random_scale: float = 1.0
    """Position standard deviation to initialize random gaussians"""
    ssim_lambda: float = 0.2
    """weight of ssim loss"""
    output_depth_during_training: bool = False
    """If True, output depth during training. Otherwise, only output depth during evaluation."""
    use_camera_optimizer: bool = False
    """Whether to use camera optimizer"""
    camera_optimizer: CameraOptimizerConfig = field(default_factory=lambda: CameraOptimizerConfig(mode="SO3xR3"))
    """Config of the camera optimizer to use"""

    # classial control
    cull_alpha_thresh: float = 0.1
    """threshold of opacity for culling gaussians. One can set it to a lower value (e.g. 0.005) for higher quality."""
    cull_scale_thresh: float = 0.5
    """threshold of scale for culling huge gaussians"""
    cull_anisotropy_thresh: float = 10.0
    """threshold of quotient of scale for culling long thin gaussians"""
    cull_grad_thresh: float = 0.0001  # 0.0004 | 0.0001 | 0
    """threshold for culling gaussians with low visibility"""
    continue_cull_post_densification: bool = True
    """If True, continue to cull gaussians post refinement"""
    reset_alpha_every: int = 30
    """Every this many refinement steps, reset the alpha"""
    densify_xy_grad_thresh: float = 0.0006  # 0.001 | 0.0006
    """threshold of positional gradient norm for densifying gaussians"""
    densify_ch_grad_thresh: float = 3e-6  # 5e-6 | 3e-6 | 1e-6
    """threshold of cylindrical harmonics gradient norm for densifying gaussians"""
    densify_size_thresh: float = 0.01
    """below this size, gaussians are *duplicated*, otherwise split"""
    n_split_samples: int = 2
    """number of samples to split gaussians into"""
    stop_split_at: int = 15000
    """stop splitting at this step"""
    cull_screen_size: float = 0.15
    """if a gaussian is more than this percent of screen space, cull it"""
    split_screen_size: float = 0.05
    """if a gaussian is more than this percent of screen space, split it"""
    stop_screen_size_at: int = 4000
    """stop culling/splitting at this step WRT screen size of gaussians"""

    # MCMC control
    refine_start_iter: int = 500
    """start MCMC refinement at this number of steps"""
    refine_stop_iter: int = 25000
    """end MCMC refinement at this number of steps"""
    cap_max: int = 1000000
    """maximum number of splats for MCMC"""
    mcmc_split_splats: bool = True
    """whether to split splats into two splats with different positions"""
    noise_lr: float = 5e3
    """MCMC sampling noise learning rate"""
    min_opacity: float = 0.005
    """minimum opacity for MCMC relocation"""

    # representation
    use_anisotropy: bool = False
    """use anisotropy for splats"""
    sh_degree: int = 3
    """maximum degree of spherical harmonics to use"""
    sh_degree_interval: int = 1000
    """every n intervals turn on another sh degree"""
    ch_degree_r: int = 0  # 3 | 0
    """maximum radial degree of cylindrical harmonics to use"""
    ch_degree_phi: int = 0  # 3 | 0
    """maximum angular degree of cylindrical harmonics to use"""
    ch_degree_interval: int = 1000
    """every n intervals turn on another ch degree"""
    max_opacity: float = 0.995
    """maximum opacity of a gaussian, prevent numerical instability during backward"""
    train_background_color: bool = True
    """make background color trainable"""
    background_sh_degree: int = 3
    """enable background model"""

    # regularization
    use_scale_regularization: bool = True
    """If enabled, a scale regularization introduced in PhysGauss (https://xpandora.github.io/PhysGaussian/) is used for reducing huge spikey gaussians."""
    max_gauss_ratio: float = 5.0
    """threshold of ratio of gaussian max to min scale before applying regularization
    loss from the PhysGaussian paper
    """
    depth_reg_weight: float = 0.02
    """Weight for depth regularizer"""
    depth_reg_warmup: int = 0
    """warmup steps for depth regularizer"""
    normal_reg_weight: float = 0.04
    """Weight for normal regularizer"""
    normal_reg_per_splat_factor: float = 0.4
    """Factor of per-splat vs overall normal regularization, 0 to 1"""
    normal_reg_warmup: int = 0
    """warmup steps for normal regularizer"""
    reg_warmup_length: int = 2000
    """Warmup for depth and normal regularizers.
       IF THE INITIALIZATION IS RANDOM, only apply regularizers after this many steps."""
    mcmc_opacity_reg: float = 0.001
    """Opacity regularization from MCMC
       Lower usually gives more accurate geometry, original MCMC uses 0.01"""
    mcmc_scale_reg: float = 0.002
    """Scale regularization from MCMC, original MCMC uses 0.01"""


class SpirulaeModel(Model):
    """Template Model."""

    config: SpirulaeModelConfig

    def __init__(
        self,
        *args,
        seed_points: Optional[Tuple[torch.Tensor, torch.Tensor]] = None,
        **kwargs,
    ):
        self.seed_points = seed_points
        super().__init__(*args, **kwargs)

    def populate_modules(self):
        if self.seed_points is not None and not self.config.random_init:
            means = torch.nn.Parameter(self.seed_points[0])
            self.random_init = False
        else:
            means = torch.nn.Parameter(torch.randn((self.config.num_random, 3)) * self.config.random_scale)
            self.random_init = True
        self.xys_grad_norm = None
        self.ch_grad_norm = None
        self.max_2Dsize = None
        distances, indices = self.k_nearest_sklearn(means.data, 6)
        distances = torch.from_numpy(distances)
        # avg_dist = distances.mean(dim=-1, keepdim=True)
        # scales = torch.nn.Parameter(torch.log(avg_dist.repeat(1, 2)))
        points = means.data[indices] - means.data[:, None, :]
        points = points.cpu().numpy()
        U, S, Vt = np.linalg.svd(points)
        num_points = means.shape[0]
        for i in range(num_points):
            sorted_indices = np.argsort(-S[i])
            S[i] = S[i][sorted_indices]
            Vt[i] = Vt[i][sorted_indices]
        scales = np.log(1.5*S[:,:2]+1e-8)
        scales = 0.5*(scales+np.flip(scales,axis=1))
        scales = torch.nn.Parameter(torch.from_numpy(scales))
        # quats = torch.nn.Parameter(random_quat_tensor(num_points))
        quats = torch.nn.Parameter(torch.from_numpy(np.array(
            [Rotation.from_matrix(R.T).as_quat() for R in Vt],
            dtype=np.float32)))
        # colors
        dim_sh = num_sh_bases(self.config.sh_degree)
        dim_ch = self.config.ch_degree_r * (2*self.config.ch_degree_phi+1)

        if (
            self.seed_points is not None
            and not self.config.random_init
            # We can have colors without points.
            and self.seed_points[1].shape[0] > 0
        ):
            shs = torch.zeros((self.seed_points[1].shape[0], dim_sh, 3)).float().cuda()
            seed_color = self.seed_points[1] / 255
            if self.config.sh_degree > 0:
                shs[:, 0, :3] = RGB2SH(seed_color)
                shs[:, 1:, 3:] = 0.0
            else:
                shs[:, 0, :3] = seed_color
            if dim_ch > 0:
                shs *= 2
            chs = torch.zeros((self.seed_points[1].shape[0], dim_ch, 3)).float().cuda()
            features_dc = torch.nn.Parameter(shs[:, 0, :])
            features_sh = torch.nn.Parameter(shs[:, 1:, :])
            features_ch = torch.nn.Parameter(chs)
        else:
            features_dc = torch.nn.Parameter(torch.rand(num_points, 3))
            features_sh = torch.nn.Parameter(torch.zeros((num_points, dim_sh-1, 3)))
            features_ch = torch.nn.Parameter(torch.zeros((num_points, dim_ch, 3)))

        opacities = torch.nn.Parameter(torch.logit(0.1 * torch.ones(num_points, 1)))
        anisotropies = torch.nn.Parameter(torch.zeros((num_points, 2)))

        gauss_params = {
            "means": means,
            "scales": scales,
            "quats": quats,
            "features_dc": features_dc,
            "features_sh": features_sh,
            "features_ch": features_ch,
            "opacities": opacities,
            "anisotropies": anisotropies
        }
        self.gauss_params = torch.nn.ParameterDict(gauss_params)

        self.camera_optimizer: CameraOptimizer = self.config.camera_optimizer.setup(
            num_cameras=self.num_train_data, device="cpu"
        )

        # metrics
        from torchmetrics.image import PeakSignalNoiseRatio
        from torchmetrics.image.lpip import \
            LearnedPerceptualImagePatchSimilarity

        self.psnr = PeakSignalNoiseRatio(data_range=1.0)
        self.ssim = SSIM(data_range=1.0, size_average=True, channel=3)
        self.lpips = LearnedPerceptualImagePatchSimilarity(normalize=True)
        self.step = 0

        self.crop_box: Optional[OrientedBox] = None
        if self.config.background_color == "random":
            self.background_color = torch.tensor(
                [0.1490, 0.1647, 0.2157]
            )  # This color is the same as the default background color in Viser. This would only affect the background color when rendering.
        else:
            self.background_color = get_color(self.config.background_color)
        self.background_color = torch.nn.Parameter(self.background_color)
        if self.config.train_background_color:
            dim_sh = num_sh_bases(self.config.background_sh_degree)
            self.background_sh = torch.nn.Parameter(torch.zeros((dim_sh-1, 3)))
        else:
            self.background_sh = None

    @property
    def colors(self):
        if self.config.sh_degree > 0:
            return SH2RGB(self.features_dc)
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
    def scales_3d(self):
        scales = self.gauss_params["scales"]
        scales_thickness = torch.amin(scales, axis=1, keepdim=True) + math.log(0.001)
        scales_thickness += -torch.inf
        return torch.concat((scales, scales_thickness), dim=1)

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

    @property
    def anisotropies(self):
        return self.gauss_params["anisotropies"]

    def load_state_dict(self, dict, **kwargs):  # type: ignore
        # resize the parameters to match the new number of points
        self.step = 30000
        if "means" in dict:
            # For backwards compatibility, we remap the names of parameters from
            # means->gauss_params.means since old checkpoints have that format
            for p in ["means", "scales", "quats",
                      "features_dc", "features_sh", "features_ch",
                      "opacities", "anisotropies"]:
                dict[f"gauss_params.{p}"] = dict[p]
        newp = dict["gauss_params.means"].shape[0]
        for name, param in self.gauss_params.items():
            old_shape = param.shape
            new_shape = (newp,) + old_shape[1:]
            self.gauss_params[name] = torch.nn.Parameter(torch.zeros(new_shape, device=self.device))
        super().load_state_dict(dict, **kwargs)

    def k_nearest_sklearn(self, x: torch.Tensor, k: int):
        """
            Find k-nearest neighbors using sklearn's NearestNeighbors.
        x: The data tensor of shape [num_samples, num_features]
        k: The number of neighbors to retrieve
        """
        # Convert tensor to numpy array
        x_np = x.cpu().numpy()

        # Build the nearest neighbors model
        from sklearn.neighbors import NearestNeighbors

        nn_model = NearestNeighbors(n_neighbors=k + 1, algorithm="auto", metric="euclidean").fit(x_np)

        # Find the k-nearest neighbors
        distances, indices = nn_model.kneighbors(x_np)

        # Exclude the point itself from the result and return
        return distances[:, 1:].astype(np.float32), indices[:, 1:].astype(np.float32)

    def remove_from_optim(self, optimizer, deleted_mask, new_params):
        """removes the deleted_mask from the optimizer provided"""
        assert len(new_params) == 1
        assert isinstance(optimizer, torch.optim.Adam), "Only works with Adam"

        param = optimizer.param_groups[0]["params"][0]
        param_state = optimizer.state[param]
        del optimizer.state[param]

        # Modify the state directly without deleting and reassigning.
        if "exp_avg" in param_state:
            param_state["exp_avg"] = param_state["exp_avg"][~deleted_mask]
            param_state["exp_avg_sq"] = param_state["exp_avg_sq"][~deleted_mask]

        # Update the parameter in the optimizer's param group.
        del optimizer.param_groups[0]["params"][0]
        del optimizer.param_groups[0]["params"]
        optimizer.param_groups[0]["params"] = new_params
        optimizer.state[new_params[0]] = param_state

    def remove_from_all_optim(self, optimizers, deleted_mask):
        param_groups = self.get_gaussian_param_groups()
        for group, param in param_groups.items():
            if group not in self.gauss_params:
                continue
            self.remove_from_optim(optimizers.optimizers[group], deleted_mask, param)
        torch.cuda.empty_cache()

    def dup_in_optim(self, optimizer, dup_mask, new_params, n=2):
        """adds the parameters to the optimizer"""
        assert isinstance(optimizer, torch.optim.Adam), "Only works with Adam"
        param = optimizer.param_groups[0]["params"][0]
        param_state = optimizer.state[param]
        if "exp_avg" in param_state:
            repeat_dims = (n,) + tuple(1 for _ in range(param_state["exp_avg"].dim() - 1))
            param_state["exp_avg"] = torch.cat(
                [
                    param_state["exp_avg"],
                    torch.zeros_like(param_state["exp_avg"][dup_mask.squeeze()]).repeat(*repeat_dims),
                ],
                dim=0,
            )
            param_state["exp_avg_sq"] = torch.cat(
                [
                    param_state["exp_avg_sq"],
                    torch.zeros_like(param_state["exp_avg_sq"][dup_mask.squeeze()]).repeat(*repeat_dims),
                ],
                dim=0,
            )
        del optimizer.state[param]
        optimizer.state[new_params[0]] = param_state
        optimizer.param_groups[0]["params"] = new_params
        del param

    def dup_in_all_optim(self, optimizers, dup_mask, n):
        param_groups = self.get_gaussian_param_groups()
        for group, param in param_groups.items():
            if group not in self.gauss_params:
                continue
            self.dup_in_optim(optimizers.optimizers[group], dup_mask, param, n)

    def reset_optim(self, optimizer, sampled_idxs):
        """removes the deleted_mask from the optimizer provided"""
        assert isinstance(optimizer, torch.optim.Adam), "Only works with Adam"

        param = optimizer.param_groups[0]["params"][0]
        param_state = optimizer.state[param]
        del optimizer.state[param]

        if "exp_avg" in param_state:
            param_state["exp_avg"][sampled_idxs] = 0.0
            param_state["exp_avg_sq"][sampled_idxs] = 0.0

    def reset_all_optim(self, optimizers, sampled_idxs):
        param_groups = self.get_gaussian_param_groups()
        for group, param in param_groups.items():
            if group not in self.gauss_params:
                continue
            self.reset_optim(optimizers.optimizers[group], sampled_idxs)
        torch.cuda.empty_cache()

    def set_crop(self, crop_box: Optional[OrientedBox]):
        self.crop_box = crop_box

    def set_background(self, background_color: torch.Tensor):
        assert background_color.shape == (3,)
        # self.background_color = background_color

    @torch.no_grad()
    def after_train(self, step: int):
        assert step == self.step
        if self.max_2Dsize is None:
            self.max_2Dsize = torch.zeros_like(self.radii, dtype=torch.float32)
        # to save some training time, we no longer need to update those stats post refinement
        if self.step >= self.config.stop_split_at:
            return
        if not hasattr(self.positions, 'absgrad'):
            return
        if not hasattr(self.features_ch, 'absgrad'):
            return
        # keep track of a moving average of grad norms
        visible_mask = (self.radii > 0).flatten()
        assert self.positions.absgrad is not None  # type: ignore
        grads = self.positions.absgrad.detach().norm(dim=-1)  # type: ignore
        ch_grads = self.features_ch.absgrad.detach().flatten()
        # print(f"grad norm min {grads.min().item()} max {grads.max().item()} mean {grads.mean().item()} size {grads.shape}")
        if self.xys_grad_norm is None or self.ch_grad_norm is None:
            self.xys_grad_norm = grads
            self.ch_grad_norm = ch_grads
            self.vis_counts = torch.ones_like(self.xys_grad_norm)
        else:
            assert self.vis_counts is not None
            self.vis_counts[visible_mask] = self.vis_counts[visible_mask] + 1
            self.xys_grad_norm[visible_mask] = grads[visible_mask] + self.xys_grad_norm[visible_mask]
            self.ch_grad_norm[visible_mask] = ch_grads[visible_mask] + self.ch_grad_norm[visible_mask]

        # update the max screen size, as a ratio of number of pixels
        # newradii = self.radii.detach()[self.depth_sort_i][visible_mask]
        newradii = self.radii.detach()[visible_mask]
        self.max_2Dsize[visible_mask] = torch.maximum(
            self.max_2Dsize[visible_mask],
            newradii / float(max(self.last_size[0], self.last_size[1])),
        )

    @torch.no_grad()
    def refinement_after(self, optimizers: Optimizers, step):
        assert step == self.step
        if self.step <= self.config.warmup_length:
            return
        # Offset all the opacity reset logic by refine_every so that we don't
        # save checkpoints right when the opacity is reset (saves every 2k)
        # then cull
        # only split/cull if we've seen every image since opacity reset
        reset_interval = self.config.reset_alpha_every * self.config.refine_every
        do_densification = (
            self.step < self.config.stop_split_at
            and self.step % reset_interval > self.num_train_data + self.config.refine_every
        )
        # assert self.xys_grad_norm is not None and self.vis_counts is not None and self.max_2Dsize is not None
        if self.xys_grad_norm is None or self.ch_grad_norm is None \
            or self.vis_counts is None or self.max_2Dsize is None:
            do_densification = False
            self.avg_xy_grad_norm = None
        else:
            self.avg_xy_grad_norm = (self.xys_grad_norm / self.vis_counts) * 0.5# * max(self.last_size[0], self.last_size[1])
            self.avg_ch_grad_norm = self.ch_grad_norm / self.vis_counts
        # with torch.no_grad():
        #     print("avg_grad_norm", torch.mean(self.avg_grad_norm).item(), torch.median(self.avg_grad_norm).item())
        #     assert False
        if do_densification:
            # then we densify
            high_grads = (self.avg_xy_grad_norm > self.config.densify_xy_grad_thresh).squeeze()
            high_grads |= (self.avg_ch_grad_norm > self.config.densify_ch_grad_thresh).squeeze()
            splits = (self.scales.exp().max(dim=-1).values > self.config.densify_size_thresh).squeeze()
            if self.step < self.config.stop_screen_size_at:
                splits |= (self.max_2Dsize > self.config.split_screen_size).squeeze()
            splits &= high_grads
            nsamps = self.config.n_split_samples
            split_params = self.split_splats(splits, nsamps)

            dups = (self.scales.exp().max(dim=-1).values <= self.config.densify_size_thresh).squeeze()
            dups &= high_grads
            dup_params = self.dup_splats(dups)
            for name, param in self.gauss_params.items():
                self.gauss_params[name] = torch.nn.Parameter(
                    torch.cat([param.detach(), split_params[name], dup_params[name]], dim=0)
                )

            # append zeros to the max_2Dsize tensor
            self.max_2Dsize = torch.cat(
                [
                    self.max_2Dsize,
                    torch.zeros_like(split_params["scales"][:, 0]),
                    torch.zeros_like(dup_params["scales"][:, 0]),
                ],
                dim=0,
            )

            split_idcs = torch.where(splits)[0]
            self.dup_in_all_optim(optimizers, split_idcs, nsamps)

            dup_idcs = torch.where(dups)[0]
            self.dup_in_all_optim(optimizers, dup_idcs, 1)

            # After a guassian is split into two new gaussians, the original one should also be pruned.
            splits_mask = torch.cat(
                (
                    splits,
                    torch.zeros(
                        nsamps * splits.sum() + dups.sum(),
                        device=self.device,
                        dtype=torch.bool,
                    ),
                )
            )

            deleted_mask = self.cull_splats(splits_mask)
        elif self.step >= self.config.stop_split_at and self.config.continue_cull_post_densification:
            deleted_mask = self.cull_splats()
        else:
            # if we donot allow culling post refinement, no more gaussians will be pruned.
            deleted_mask = None

        if deleted_mask is not None:
            self.remove_from_all_optim(optimizers, deleted_mask)

        if self.step < self.config.stop_split_at and self.step % reset_interval == self.config.refine_every:
            # Reset value is set to be twice of the cull_alpha_thresh
            reset_value = self.config.cull_alpha_thresh * 2.0
            self.opacities.data = torch.clamp(
                self.opacities.data,
                max=torch.logit(torch.tensor(reset_value, device=self.device)).item(),
            )
            # reset the exp of optimizer
            optim = optimizers.optimizers["opacities"]
            param = optim.param_groups[0]["params"][0]
            param_state = optim.state[param]
            param_state["exp_avg"] = torch.zeros_like(param_state["exp_avg"])
            param_state["exp_avg_sq"] = torch.zeros_like(param_state["exp_avg_sq"])

        self.xys_grad_norm = None
        self.vis_counts = None
        self.avg_xy_grad_norm = None
        self.avg_ch_grad_norm = None
        self.max_2Dsize = None

    @torch.no_grad()
    def cull_splats(self, extra_cull_mask: Optional[torch.Tensor] = None):
        """
        This function deletes gaussians with under a certain opacity threshold
        extra_cull_mask: a mask indicates extra gaussians to cull besides existing culling criterion
        """
        n_bef = self.num_points
        # cull transparent ones
        culls = (torch.sigmoid(self.opacities) < self.config.cull_alpha_thresh).squeeze()
        below_alpha_count = torch.sum(culls).item()
        # cull ones with low visibility
        low_grad_count = 0
        if self.avg_xy_grad_norm is not None:
            low_grads = torch.zeros_like(culls)
            low_grads[:len(self.avg_xy_grad_norm)] = (self.avg_xy_grad_norm < self.config.cull_grad_thresh).squeeze()
            culls |= low_grads
            low_grad_count = torch.sum(low_grads).item()
        # others
        toobigs_count = 0
        toothins_count = 0
        if extra_cull_mask is not None:
            culls = culls | extra_cull_mask
        if self.step > self.config.refine_every * self.config.reset_alpha_every:
            # cull huge ones
            toobigs = (torch.exp(self.scales).max(dim=-1).values > self.config.cull_scale_thresh).squeeze()
            if self.step < self.config.stop_screen_size_at:
                # cull big screen space
                assert self.max_2Dsize is not None
                toobigs = toobigs | (self.max_2Dsize > self.config.cull_screen_size).squeeze()
            culls = culls | toobigs
            toobigs_count = torch.sum(toobigs).item()
            # cull long/thin ones
            toothins = (
                torch.exp(self.scales).max(dim=-1).values /
                torch.exp(self.scales).min(dim=-1).values > \
                    self.config.cull_anisotropy_thresh).squeeze()
            culls = culls | toothins
            toothins_count = torch.sum(toothins).item()
        for name, param in self.gauss_params.items():
            self.gauss_params[name] = torch.nn.Parameter(param[~culls])

        CONSOLE.log(
            f"Culled {n_bef - self.num_points} gaussians "
            f"({below_alpha_count} below alpha thresh, {low_grad_count} low visibility, {toobigs_count} too bigs, {toothins_count} too thins, {self.num_points} remaining)"
        )

        return culls

    @torch.no_grad()
    def split_splats(self, split_mask, samps):
        """
        This function splits gaussians that are too large
        """
        n_splits = split_mask.sum().item()
        CONSOLE.log(f"Splitting {split_mask.sum().item()/self.num_points} gaussians: {n_splits}/{self.num_points}")
        centered_samples = torch.randn((samps * n_splits, 3), device=self.device)  # Nx3 of axis-aligned scales
        scaled_samples = (
            torch.exp(self.scales_3d[split_mask].repeat(samps, 1)) * \
                # torch.exp(0.3*centered_samples)
                (centered_samples)
        )  # how these scales are rotated
        quats = self.quats[split_mask] / self.quats[split_mask].norm(dim=-1, keepdim=True)  # normalize them first
        rots = quat_to_rotmat(quats.repeat(samps, 1))  # how these scales are rotated
        rotated_samples = torch.bmm(rots, scaled_samples[..., None]).squeeze()
        new_means = 0.4*rotated_samples + self.means[split_mask].repeat(samps, 1)
        # step 2, sample new colors
        new_features_dc = self.features_dc[split_mask].repeat(samps, 1)
        new_features_sh = self.features_sh[split_mask].repeat(samps, 1, 1)
        new_features_ch = 0.0 * self.features_ch[split_mask].repeat(samps, 1, 1)
        # step 3, sample new opacities
        new_opacities = self.opacities[split_mask].repeat(samps, 1)
        # new_opacities = 1.0-torch.sqrt(1.0-torch.clip(new_opacities,0.,1.))
        # self.opacities[split_mask] = 1.0-torch.sqrt(1.0-torch.clip(self.opacities[split_mask],0.,1.))
        new_anisotropies = self.anisotropies[split_mask].repeat(samps, 1)
        # step 4, sample new scales
        size_fac = 1.6
        new_scales = torch.log(torch.exp(self.scales[split_mask]) / size_fac).repeat(samps, 1)
        self.scales[split_mask] = torch.log(torch.exp(self.scales[split_mask]) / size_fac)
        # step 5, sample new quats
        new_quats = self.quats[split_mask].repeat(samps, 1)
        out = {
            "means": new_means,
            "features_dc": new_features_dc,
            "features_sh": new_features_sh,
            "features_ch": new_features_ch,
            "opacities": new_opacities,
            "anisotropies": new_anisotropies,
            "scales": new_scales,
            "quats": new_quats,
        }
        for name, param in self.gauss_params.items():
            if name not in out:
                out[name] = param[split_mask].repeat(samps, 1)
        return out

    @torch.no_grad()
    def dup_splats(self, dup_mask):
        """
        This function duplicates gaussians that are too small
        """
        n_dups = dup_mask.sum().item()
        CONSOLE.log(f"Duplicating {dup_mask.sum().item()/self.num_points} gaussians: {n_dups}/{self.num_points}")
        new_dups = {}
        for name, param in self.gauss_params.items():
            new_dups[name] = param[dup_mask]
        new_opacities = new_dups['opacities']
        # new_opacities = 1.0 - torch.sqrt(1.0-torch.clip(new_opacities,0.,1.))
        new_dups['opacities'] = new_opacities
        self.opacities[dup_mask] = new_opacities
        return new_dups

    @torch.no_grad()
    def mcmc_after_train(self, optimizers: Optimizers, step: int):
        optimizer = optimizers.optimizers['means']
        for param_group in optimizer.param_groups:
            lr = param_group['lr']
            break
        self.mcmc_add_noise_to_splats(lr)

    @torch.no_grad()
    def mcmc_refinement_after(self, optimizers: Optimizers, step: int):
        assert step == self.step
        if self.step <= self.config.refine_start_iter or \
            self.step > self.config.refine_stop_iter:
            return
        num_relocated = self.mcmc_relocate_splats(optimizers)
        CONSOLE.log(f"Step {step}: Relocated {num_relocated} splats.")
        num_new = self.mcmc_add_new_splats(optimizers)
        CONSOLE.log(f"Step {step}: Added {num_new} splats.")

    @torch.no_grad()
    def mcmc_relocate_splats(self, optimizers: Optimizers) -> int:
        dead_mask = torch.sigmoid(self.opacities.squeeze()) <= self.config.min_opacity
        dead_indices = dead_mask.nonzero(as_tuple=True)[0]
        alive_indices = (~dead_mask).nonzero(as_tuple=True)[0]
        num_gs = len(dead_indices)
        if num_gs <= 0:
            return num_gs
        
        eps = torch.finfo(torch.float32).eps
        probs = torch.sigmoid(self.opacities)[alive_indices].squeeze()
        probs = probs / (probs.sum() + eps)
        sampled_idxs = torch.multinomial(probs, num_gs, replacement=True)
        sampled_idxs = alive_indices[sampled_idxs]

        if self.config.mcmc_split_splats:
            new_position_offsets, new_opacities, new_scales = compute_relocation_split(
                self.positions[sampled_idxs],
                self.quats[sampled_idxs],
                torch.sigmoid(self.opacities)[sampled_idxs],
                torch.exp(self.scales)[sampled_idxs],
            )
            new_opacities = torch.clamp(
                new_opacities, max=0.9, min=self.config.min_opacity)
            self.opacities[sampled_idxs] = torch.logit(new_opacities)
            self.scales[sampled_idxs] = torch.log(new_scales)
            for name, param in self.gauss_params.items():
                self.gauss_params[name][dead_indices] = param[sampled_idxs]
            self.means[dead_indices] += new_position_offsets
            self.means[sampled_idxs] -= new_position_offsets

        else:
            new_opacities, new_scales = compute_relocation(
                torch.sigmoid(self.opacities)[sampled_idxs],
                torch.exp(self.scales)[sampled_idxs],
                torch.bincount(sampled_idxs)[sampled_idxs] + 1,
            )
            new_opacities = torch.clamp(
                new_opacities, max=0.9, min=self.config.min_opacity)
            self.opacities[sampled_idxs] = torch.logit(new_opacities)
            self.scales[sampled_idxs] = torch.log(new_scales)
            for name, param in self.gauss_params.items():
                self.gauss_params[name][dead_indices] = param[sampled_idxs]

        self.reset_all_optim(optimizers, dead_indices)
        self.reset_all_optim(optimizers, sampled_idxs)

        return num_gs

    @torch.no_grad()
    def mcmc_add_new_splats(self, optimizers: Optimizers) -> int:
        current_num_points = len(self.opacities)
        target_num = min(self.config.cap_max, int(1.05 * current_num_points)+1)
        num_gs = max(0, target_num - current_num_points)
        if num_gs <= 0:
            return num_gs

        # Sample for new splats
        eps = torch.finfo(torch.float32).eps
        probs = torch.sigmoid(self.opacities).squeeze()
        probs = probs / (probs.sum() + eps)
        sampled_idxs = torch.multinomial(probs, num_gs, replacement=True)

        if self.config.mcmc_split_splats:
            new_position_offsets, new_opacities, new_scales = compute_relocation_split(
                self.positions[sampled_idxs],
                self.quats[sampled_idxs],
                torch.sigmoid(self.opacities)[sampled_idxs],
                torch.exp(self.scales)[sampled_idxs],
            )
            new_opacities = torch.clamp(
                new_opacities, max=0.9, min=self.config.min_opacity)
            self.opacities[sampled_idxs] = torch.logit(new_opacities)
            self.scales[sampled_idxs] = torch.log(new_scales)
            for name, param in self.gauss_params.items():
                self.gauss_params[name] = torch.concatenate((param, param[sampled_idxs]))
            self.means[sampled_idxs] += new_position_offsets
            self.means[current_num_points:] -= new_position_offsets

        else:
            new_opacities, new_scales = compute_relocation(
                opacities=torch.sigmoid(self.opacities)[sampled_idxs],
                scales=torch.exp(self.scales)[sampled_idxs],
                ratios=torch.bincount(sampled_idxs)[sampled_idxs] + 1,
            )
            new_opacities = torch.clamp(
                new_opacities, max=0.9, min=self.config.min_opacity)
            self.opacities[sampled_idxs] = torch.logit(new_opacities)
            self.scales[sampled_idxs] = torch.log(new_scales)
            for name, param in self.gauss_params.items():
                self.gauss_params[name] = torch.concatenate((param, param[sampled_idxs]))

        self.dup_in_all_optim(optimizers, sampled_idxs, 2)
        return num_gs

    @torch.no_grad()
    def mcmc_add_noise_to_splats(self, last_lr):
        opacities = torch.sigmoid(self.opacities)
        scales = torch.exp(self.scales_3d)
        scales[:, 2] = 0.5 * torch.fmin(scales[:,0], scales[:,1])

        R = quat_to_rotmat(self.quats)  # (..., 3, 3)
        M = R * scales[..., None, :]  # (..., 3, 3)
        covars = torch.bmm(M, M.transpose(-1, -2))  # (..., 3, 3)

        def op_sigmoid(x, k=100, x0=1.0-self.config.min_opacity):
            return 1 / (1 + torch.exp(-k * (x - x0)))

        noise = (
            torch.randn_like(self.means)
            * (op_sigmoid(1 - opacities))
            * self.config.noise_lr
            * last_lr
        )
        noise = torch.bmm(covars, noise.unsqueeze(-1)).squeeze(-1)
        self.gauss_params['means'] = self.gauss_params['means'] + noise

    def get_training_callbacks(
        self, training_callback_attributes: TrainingCallbackAttributes
    ) -> List[TrainingCallback]:
        cbs = []
        cbs.append(TrainingCallback([TrainingCallbackLocation.BEFORE_TRAIN_ITERATION], self.step_cb))
        # return cbs
        if self.config.use_mcmc:
            cbs.append(
                TrainingCallback(
                    [TrainingCallbackLocation.AFTER_TRAIN_ITERATION],
                    self.mcmc_after_train,
                    args=[training_callback_attributes.optimizers],
                )
            )
            cbs.append(
                TrainingCallback(
                    [TrainingCallbackLocation.AFTER_TRAIN_ITERATION],
                    self.mcmc_refinement_after,
                    update_every_num_iters=self.config.refine_every,
                    args=[training_callback_attributes.optimizers],
                )
            )
        else:
            # The order of these matters
            cbs.append(
                TrainingCallback(
                    [TrainingCallbackLocation.AFTER_TRAIN_ITERATION],
                    self.after_train,
                )
            )
            cbs.append(
                TrainingCallback(
                    [TrainingCallbackLocation.AFTER_TRAIN_ITERATION],
                    self.refinement_after,
                    update_every_num_iters=self.config.refine_every,
                    args=[training_callback_attributes.optimizers],
                )
            )
        return cbs

    def step_cb(self, step):
        self.step = step

    def get_gaussian_param_groups(self) -> Dict[str, List[Parameter]]:
        # Here we explicitly use the means, scales as parameters so that the user can override this function and
        # specify more if they want to add more optimizable params to gaussians.
        return {
            name: [self.gauss_params[name]]
            for name in [
                "means", "scales", "quats",
                "features_dc", "features_sh", "features_ch",
                "opacities"
            ] + ["anisotropies"] * self.config.use_anisotropy
        }

    def get_param_groups(self) -> Dict[str, List[Parameter]]:
        """Obtain the parameter groups for the optimizers

        Returns:
            Mapping of different parameter groups
        """
        gps = self.get_gaussian_param_groups()
        if self.config.train_background_color:
            gps["background_color"] = [self.background_color]
            if self.background_sh is not None:
                gps["background_sh"] = [self.background_sh]
        if self.config.use_camera_optimizer:
            self.camera_optimizer.get_param_groups(param_groups=gps)
        return gps

    def _get_downscale_factor(self):
        if self.training:
            return 2 ** max(
                (self.config.num_downscales - self.step // self.config.resolution_schedule),
                0,
            )
        else:
            return 1

    def _downscale_if_required(self, image):
        d = self._get_downscale_factor()
        if d > 1:
            return resize_image(image, d)
        return image

    def get_background_image(self, camera: Cameras, mask: Optional[torch.Tensor]=None):
        if not isinstance(camera, Cameras):
            print("Called get_background_image with not a camera")
            return {}
        assert camera.shape[0] == 1, "Only one camera at a time"
        camera_downscale = self._get_downscale_factor()
        camera.rescale_output_resolution(1 / camera_downscale)

        W, H = int(camera.width.item()), int(camera.height.item())

        background = self.background_color.repeat(H, W, 1)
        sh_degree = self.config.background_sh_degree
        if not (sh_degree > 0):
            return background

        # camera.generate_rays seems to be slow, write from scratch
        if False:
            ray_bundle = camera.generate_rays(
                camera_indices=0, keep_shape=False, disable_distortion=True)
            if mask is not None:
                mask_indices = torch.where(mask.flatten())
                ray_bundle = ray_bundle[mask_indices]
            if not len(ray_bundle) > 0:
                return background
            directions = ray_bundle.directions
        else:
            device = background.device
            y, x = torch.meshgrid(torch.arange(H), torch.arange(W), indexing="ij")
            x = x.flatten().float().to(device) + 0.5
            y = y.flatten().float().to(device) + 0.5
            if mask is not None:
                mask_indices = torch.where(mask.flatten())
                x, y = x[mask_indices], y[mask_indices]
            if not len(x) > 0:
                return background
            fx, fy = camera.fx[0].item(), camera.fy[0].item()
            cx, cy = camera.cx[0].item(), camera.cy[0].item()
            coord = torch.stack([(x - cx) / fx, -(y - cy) / fy, -torch.ones_like(x)], -1)
            rotation = camera.camera_to_worlds[0][:3, :3]
            directions = torch.matmul(coord, rotation.T)
            norm = torch.maximum(torch.linalg.norm(directions, dim=-1, keepdims=True), torch.tensor([1e-6]).to(device))
            directions = directions / norm

        sh_components = components_from_spherical_harmonics(sh_degree+1, directions)
        sh_coeffs = torch.cat((self.background_color.unsqueeze(0), self.background_sh), dim=0)
        bg_flat = torch.matmul(sh_components, sh_coeffs)
        bg_flat = torch.clamp(bg_flat+0.5, min=0.0)

        if mask is not None:
            background = background.view(-1, 3)
            background[mask_indices] = bg_flat
        else:
            background = bg_flat

        return background.view(H, W, 3)

    @staticmethod
    def get_empty_outputs(width: int, height: int, background: torch.Tensor) -> Dict[str, Union[torch.Tensor, List]]:
        rgb = background.repeat(height, width, 1)
        depth = background.new_ones(*rgb.shape[:2], 1) * 10
        alpha = background.new_zeros(*rgb.shape[:2], 1)
        depth_grad_1_vis = background.repeat(height, width, 1)
        depth_grad_2_vis = background.repeat(height, width, 1)
        reg_depth = background.new_zeros(*rgb.shape[:2], 1)
        reg_normal = background.new_zeros(*rgb.shape[:2], 1)
        depth_grad_1 = background.repeat(height, width, 1)
        depth_grad_2 = background.repeat(height, width, 1)
        return {
            "rgb": rgb,
            "depth": depth,
            "depth_grad_1_vis": depth_grad_1_vis,
            "depth_grad_2_vis": depth_grad_2_vis,
            "reg_depth": reg_depth,
            "reg_normal": reg_normal,
            "alpha": alpha,
            "depth_grad_1": depth_grad_1,
            "depth_grad_2": depth_grad_2,
            "background": background
        }

    def get_outputs(self, camera: Cameras) -> Dict[str, Union[torch.Tensor, List]]:
        """Takes in a Ray Bundle and returns a dictionary of outputs.

        Args:
            ray_bundle: Input bundle of rays. This raybundle should have all the
            needed information to compute the outputs.

        Returns:
            Outputs of model. (ie. rendered colors)
        """
        if not isinstance(camera, Cameras):
            print("Called get_outputs with not a camera")
            return {}
        assert camera.shape[0] == 1, "Only one camera at a time"

        # get the background color
        if self.training and self.config.use_camera_optimizer:
            optimized_camera_to_world = self.camera_optimizer.apply_to_camera(camera)[0, ...]
        else:
            optimized_camera_to_world = camera.camera_to_worlds[0, ...]

        if self.crop_box is not None and not self.training:
            crop_ids = self.crop_box.within(self.means).squeeze()
            if crop_ids.sum() == 0:
                return self.get_empty_outputs(
                    int(camera.width.item()), int(camera.height.item()),
                    self.background_color.detach())
        else:
            crop_ids = None
        camera_downscale = self._get_downscale_factor()
        camera.rescale_output_resolution(1 / camera_downscale)
        # shift the camera to center of scene looking at center
        R = optimized_camera_to_world[:3, :3]  # 3 x 3
        T = optimized_camera_to_world[:3, 3:4]  # 3 x 1

        # flip the z and y axes to align with gsplat conventions
        R_edit = torch.diag(torch.tensor([1, -1, -1], device=self.device, dtype=R.dtype))
        R = R @ R_edit
        # analytic matrix inverse to get world2camera matrix
        R_inv = R.T
        T_inv = -R_inv @ T
        viewmat = torch.eye(4, device=R.device, dtype=R.dtype)
        viewmat[:3, :3] = R_inv
        viewmat[:3, 3:4] = T_inv
        # calculate the FOV of the camera given fx and fy, width and height
        cx = camera.cx.item()
        cy = camera.cy.item()
        W, H = int(camera.width.item()), int(camera.height.item())
        self.last_size = (H, W)

        if crop_ids is not None:
            opacities_crop = self.opacities[crop_ids]
            anisotropies_crop = self.anisotropies[crop_ids]
            means_crop = self.means[crop_ids]
            features_dc_crop = self.features_dc[crop_ids]
            features_sh_crop = self.features_sh[crop_ids]
            features_ch_crop = self.features_ch[crop_ids]
            scales_crop = self.scales[crop_ids]
            quats_crop = self.quats[crop_ids]
        else:
            opacities_crop = self.opacities
            anisotropies_crop = self.anisotropies
            means_crop = self.means
            features_dc_crop = self.features_dc
            features_sh_crop = self.features_sh
            features_ch_crop = self.features_ch
            scales_crop = self.scales
            quats_crop = self.quats
        colors_crop = torch.cat((features_dc_crop[:, None, :], features_sh_crop), dim=1)

        quats_norms = torch.norm(quats_crop, p=2, dim=1)
        quats_norms = quats_norms.unsqueeze(1)
        quats_norms = torch.where(quats_norms == 0,
                                  torch.tensor([1.0],device=quats_crop.device),
                                  quats_norms)
        quats_crop = quats_crop / quats_norms

        if self.config.sh_degree > 0:
            viewdirs = means_crop.detach() - optimized_camera_to_world.detach()[:3, 3]  # (N, 3)
            n = min(self.step // self.config.sh_degree_interval, self.config.sh_degree)
            rgbs = spherical_harmonics(n, viewdirs, colors_crop)  # input unnormalized viewdirs
            rgbs = torch.clamp(rgbs + 0.5, min=0.0)  # type: ignore
        else:
            rgbs = colors_crop[:, 0, :]

        # print(self.config.sh_degree, rgbs.shape)

        BLOCK_WIDTH = 16  # this controls the tile size of rasterization, 16 is a good default
        intrins = (camera.fx.item(), camera.fy.item(), cx, cy)
        (
            positions, axes_u, axes_v,
            depth_grads,
            bounds, num_tiles_hit
        ) = project_gaussians(  # type: ignore
            means_crop,
            torch.exp(scales_crop),
            quats_crop,
            viewmat.squeeze()[:3, :],
            intrins,
            H, W,
            BLOCK_WIDTH,
        )  # type: ignore

        # rescale the camera back to original dimensions before returning
        camera.rescale_output_resolution(camera_downscale)

        opacities = self.config.max_opacity * torch.sigmoid(opacities_crop)

        depth_im_ref = rasterize_gaussians_depth(
            positions, axes_u, axes_v,
            opacities, anisotropies_crop,
            bounds, num_tiles_hit,
            intrins, H, W, BLOCK_WIDTH
        )
        depth_im_ref = torch.where(
            depth_im_ref > 0.0, depth_im_ref,
            torch.amax(depth_im_ref).detach()
        ).contiguous()

        # finite difference gradient
        depth_normal_grad_ref = torch.zeros_like(depth_im_ref).repeat(1,1,2)
        depth_normal_grad_ref[1:-1,1:-1] = torch.concatenate((
            (depth_im_ref[1:-1,2:]-depth_im_ref[1:-1,:-2])/2,
            (depth_im_ref[2:,1:-1]-depth_im_ref[:-2,1:-1])/2
        ), axis=2)
        depth_normal_grad_ref = depth_normal_grad_ref.contiguous()
        depth_normal_grad_ref_norm = torch.norm(depth_normal_grad_ref, dim=2, keepdim=True)
        depth_normal_grad_ref_normalized = torch.where(
            depth_normal_grad_ref_norm > 0,
            depth_normal_grad_ref / (depth_normal_grad_ref_norm+1e-4),
            0.0).contiguous()
        depth_ref_im = torch.concatenate(
            (depth_normal_grad_ref_normalized, depth_im_ref), dim=-1)

        # main rasterization
        background_color = self.background_color
        if self.config.background_sh_degree > 0:
            background_color = torch.zeros_like(background_color)
        ch_degree = self.step // self.config.ch_degree_interval
        (
            rgb, depth_im,
            reg_depth, reg_normal,
            alpha
        ) = rasterize_gaussians(  # type: ignore
            positions, axes_u, axes_v,
            rgbs,
            self.config.ch_degree_r, min(ch_degree, self.config.ch_degree_r),
            self.config.ch_degree_phi, min(ch_degree, self.config.ch_degree_phi),
            features_ch_crop,
            opacities, anisotropies_crop,
            depth_grads, depth_ref_im,
            bounds, num_tiles_hit,
            intrins, H, W, BLOCK_WIDTH,
            background_color
        )  # type: ignore
        alpha = alpha[..., None]
        rgb = torch.clamp(rgb, max=1.0)  # type: ignore
        depth_grad_im, depth_im = depth_im[...,0:2], depth_im[...,2:3]
        depth_grad_im = torch.where(alpha > 0, depth_grad_im / alpha, 0.0).contiguous()
        depth_im = torch.where(alpha > 0, depth_im / alpha, depth_im.detach().max()).contiguous()

        depth_grad_im_norm = torch.norm(depth_grad_im, dim=2, keepdim=True)
        depth_grad_im_normalized = torch.where(
            depth_grad_im_norm > 0,
            depth_grad_im / (depth_grad_im_norm+1e-4),
            0.0).contiguous()

        depth_grad_im_vis = None
        depth_grad_ref_vis = None
        anisotropy_vis = None
        if self.config.output_depth_during_training or not self.training:
            def depth_grad_vis(im):
                im_norm = torch.norm(im, dim=2, keepdim=True)
                hue_im = torch.atan2(im[...,1:2], im[...,0:1])
                return torch.clip(
                    torch.tanh(im_norm * 0.5*torch.numel(depth_im)**0.5) * \
                        (torch.concat((
                        torch.cos(hue_im),
                        torch.cos(hue_im-2.0*np.pi/3.0),
                        torch.cos(hue_im+2.0*np.pi/3.0)
                    ), axis=2)*0.5+0.5) * 2.0, 0.0, 1.0)
            with torch.no_grad():
                depth_grad_im_vis = depth_grad_vis(depth_grad_im)
                depth_grad_ref_vis = depth_grad_vis(0.1*depth_normal_grad_ref)

                pad_zero = torch.zeros_like(opacities)
                anisotropy_norm = torch.norm(anisotropies_crop, dim=1, keepdim=True)
                meta, _ = rasterize_gaussians_simple(
                    positions, axes_u, axes_v,
                    torch.concatenate((anisotropy_norm, pad_zero, pad_zero), dim=1),
                    opacities, anisotropies_crop,
                    bounds, num_tiles_hit,
                    intrins, H, W, BLOCK_WIDTH
                )
                anisotropy_vis = meta[..., 0:1]

                # print(torch.amin(depth_grad_norm_im).item(), torch.mean(depth_grad_norm_im).item(), torch.amax(depth_grad_norm_im).item())
                # assert (reg_depth > -1e-6).all()
                # print(torch.mean(reg_depth).item(), torch.mean(reg_normal).item())
                # pass
                # depth_grad_reg = self.depth_grads.abs().mean()
                # print(depth_grad_reg.item())

                # print(self.features_ch.mean().item(), self.features_ch.std().item())
                # print(self.anisotropies.mean().item(), self.anisotropies.std().item())

        self.intrins = intrins
        self.positions = positions
        self.radii = num_tiles_hit**0.5 / 2 * BLOCK_WIDTH

        # blend with background
        if self.config.background_sh_degree > 0:
            mask = (alpha < 0.98)
            background = self.get_background_image(camera, mask)
            rgb = rgb + (1.0 - alpha) * background

        # clamp pixel to between 0 and 1
        # rgb = torch.clip(rgb, 0.0, 1.0)
        rgb = saturate_keep_gradient(rgb)

        # do this for L2 loss
        # reg_depth = torch.sqrt(torch.relu(reg_depth+0.01))

        return {
            "rgb": rgb,
            "depth": depth_im_ref,
            # "depth_diff": torch.abs(depth_im-depth_ref),
            # "alpha_diff": torch.abs(alpha-alpha_ref),
            "depth_grad_1_vis": depth_grad_im_vis,
            "depth_grad_2_vis": depth_grad_ref_vis,
            "reg_depth": reg_depth.unsqueeze(2),
            "reg_normal": reg_normal.unsqueeze(2),
            "anisotropy_vis": anisotropy_vis,
            "alpha": alpha,
            "depth_grad_1": depth_grad_im_normalized,
            "depth_grad_2": depth_normal_grad_ref_normalized,
            "background": self.background_color
        }  # type: ignore

    def get_gt_img(self, image: torch.Tensor):
        """Compute groundtruth image with iteration dependent downscale factor for evaluation purpose

        Args:
            image: tensor.Tensor in type uint8 or float32
        """
        if image.dtype == torch.uint8:
            image = image.float() / 255.0
        gt_img = self._downscale_if_required(image)
        return gt_img.to(self.device)

    def composite_with_background(self, image, background) -> torch.Tensor:
        """Composite the ground truth image with a background color when it has an alpha channel.

        Args:
            image: the image to composite
            background: the background color
        """
        if image.shape[2] == 4:
            alpha = image[..., -1].unsqueeze(-1).repeat((1, 1, 3))
            return alpha * image[..., :3] + (1 - alpha) * background
        else:
            return image

    def get_metrics_dict(self, outputs, batch) -> Dict[str, torch.Tensor]:
        """Compute and returns metrics.

        Args:
            outputs: the output to compute loss dict to
            batch: ground truth batch corresponding to outputs
        """
        gt_rgb = self.composite_with_background(self.get_gt_img(batch["image"]), outputs["background"])
        metrics_dict = {}
        predicted_rgb = outputs["rgb"]
        metrics_dict["psnr"] = self.psnr(predicted_rgb, gt_rgb)

        metrics_dict["gaussian_count"] = self.num_points

        if self.config.use_camera_optimizer:
            self.camera_optimizer.get_metrics_dict(metrics_dict)
        return metrics_dict

    def get_loss_dict(self, outputs, batch, metrics_dict=None) -> Dict[str, torch.Tensor]:
        """Computes and returns the losses dict.

        Args:
            outputs: the output to compute loss dict to
            batch: ground truth batch corresponding to outputs
            metrics_dict: dictionary of metrics, some of which we can use for loss
        """
        gt_img_rgba = self.get_gt_img(batch["image"])
        gt_img = self.composite_with_background(gt_img_rgba, outputs["background"])
        pred_img = outputs["rgb"]

        alpha_loss = 0.0
        if gt_img_rgba.shape[2] == 4:
            alpha = gt_img_rgba[..., -1].unsqueeze(-1)
            alpha_loss = alpha_loss + torch.relu(outputs['alpha']-alpha).mean()
            print(alpha_loss.item())

        # Set masked part of both ground-truth and rendered image to black.
        # This is a little bit sketchy for the SSIM loss.
        if "mask" in batch:
            # batch["mask"] : [H, W, 1]
            mask = self._downscale_if_required(batch["mask"])
            mask = mask.to(self.device)
            assert mask.shape[:2] == gt_img.shape[:2] == pred_img.shape[:2]
            gt_img = gt_img * mask
            pred_img = pred_img * mask
            # compute alpha loss
            alpha_loss = alpha_loss + torch.relu(outputs['alpha']-mask).mean()

        Ll1 = torch.abs(gt_img - pred_img).mean()
        simloss = 1 - self.ssim(gt_img.permute(2, 0, 1)[None, ...], pred_img.permute(2, 0, 1)[None, ...])
        if self.config.use_scale_regularization and self.step % 10 == 0:
            scale_exp = torch.exp(self.scales)
            scale_reg = (
                torch.maximum(
                    scale_exp.amax(dim=-1) / scale_exp.amin(dim=-1),
                    torch.tensor(self.config.max_gauss_ratio),
                )
                - self.config.max_gauss_ratio
            )
            scale_reg = 0.1 * scale_reg.mean()
        else:
            scale_reg = torch.tensor(0.0).to(self.device)

        # depth regularizer
        alpha = outputs['alpha']
        reg_depth = outputs["reg_depth"]
        reg_normal = outputs["reg_normal"]
        weight_depth_reg = self.config.depth_reg_weight * \
            min(self.step / max(self.config.depth_reg_warmup, 1), 1)
        weight_normal_reg = self.config.normal_reg_weight * \
            min(self.step / max(self.config.normal_reg_warmup, 1), 1)
        if self.step < self.config.reg_warmup_length and self.random_init:
            weight_depth_reg, weight_normal_reg = 0.0, 0.0
        depth_reg = weight_depth_reg * reg_depth.sum() / alpha.sum()
        normal_reg_per_splat = weight_normal_reg * reg_normal.sum() / alpha.sum()

        # normal regularizer
        depth_grad_1 = outputs['depth_grad_1']
        depth_grad_2 = outputs['depth_grad_2']
        dot = (depth_grad_1 * depth_grad_2).sum(dim=2,keepdim=True)
        normal_reg_s = 1.0-dot
        normal_reg_overall = weight_normal_reg * (normal_reg_s*alpha).sum() / alpha.sum()
        normal_reg = normal_reg_overall + (normal_reg_per_splat-normal_reg_overall) * \
            self.config.normal_reg_per_splat_factor
        
        # MCMC regularizers
        use_mcmc = self.config.use_mcmc
        mcmc_opacity_reg = use_mcmc * self.config.mcmc_opacity_reg * \
            torch.abs(torch.sigmoid(self.opacities)).mean()
        mcmc_scale_reg = use_mcmc * self.config.mcmc_scale_reg * \
            torch.abs(torch.exp(self.scales)).mean()

        # regularizations for parameters
        quat_norm_reg = 0.1 * (torch.log(self.quats.norm(dim=-1)+0.001)**2).mean()

        loss_dict = {
            "main_loss": (1 - self.config.ssim_lambda) * Ll1 + self.config.ssim_lambda * simloss,
            "alpha_loss": alpha_loss,
            "scale_reg": scale_reg,
            "depth_reg": depth_reg,
            "normal_reg": normal_reg,
            "quat_reg": quat_norm_reg,
            'mcmc_opacity_reg': mcmc_opacity_reg,
            'mcmc_scale_reg': mcmc_scale_reg,
        }

        if self.training and self.config.use_camera_optimizer:
            # Add loss from camera optimizer
            self.camera_optimizer.get_loss_dict(loss_dict)

        return loss_dict

    @torch.no_grad()
    def get_outputs_for_camera(self, camera: Cameras, obb_box: Optional[OrientedBox] = None) -> Dict[str, torch.Tensor]:
        """Takes in a camera, generates the raybundle, and computes the output of the model.
        Overridden for a camera-based gaussian model.

        Args:
            camera: generates raybundle
        """
        assert camera is not None, "must provide camera to gaussian model"
        self.set_crop(obb_box)
        outs = self.get_outputs(camera.to(self.device))
        return outs  # type: ignore

    def get_image_metrics_and_images(
        self, outputs: Dict[str, torch.Tensor], batch: Dict[str, torch.Tensor]
    ) -> Tuple[Dict[str, float], Dict[str, torch.Tensor]]:
        """Writes the test image outputs.

        Args:
            image_idx: Index of the image.
            step: Current step.
            batch: Batch of data.
            outputs: Outputs of the model.

        Returns:
            A dictionary of metrics.
        """
        gt_rgb = self.composite_with_background(self.get_gt_img(batch["image"]), outputs["background"])
        predicted_rgb = outputs["rgb"]

        combined_rgb = torch.cat([gt_rgb, predicted_rgb], dim=1)

        # Switch images from [H, W, C] to [1, C, H, W] for metrics computations
        gt_rgb = torch.moveaxis(gt_rgb, -1, 0)[None, ...]
        predicted_rgb = torch.moveaxis(predicted_rgb, -1, 0)[None, ...]

        psnr = self.psnr(gt_rgb, predicted_rgb)
        ssim = self.ssim(gt_rgb, predicted_rgb)
        lpips = self.lpips(gt_rgb, predicted_rgb)

        # all of these metrics will be logged as scalars
        metrics_dict = {"psnr": float(psnr.item()), "ssim": float(ssim)}  # type: ignore
        metrics_dict["lpips"] = float(lpips)
        metrics_dict["gaussian_count"] = float(self.num_points)

        images_dict = {"img": combined_rgb}

        return metrics_dict, images_dict
