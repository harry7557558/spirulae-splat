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
2DGS implementation that combines many recent advancements.
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
import torch.nn.functional as F
from pytorch_msssim import SSIM
from fused_ssim import fused_ssim

from spirulae_splat.splat._torch_impl import depth_map, depth_inv_map
from spirulae_splat.splat import (
    project_gaussians,
    rasterize_gaussians_simple,
    rasterize_gaussians_depth,
    rasterize_gaussians,
    rasterize_gaussians_simplified,
    rasterize_gaussians_indices,
    rasterize_gaussians_simple_sorted,
    rasterize_gaussians_depth_sorted,
    rasterize_gaussians_sorted,
    rasterize_gaussians_simplified_sorted,
    depth_to_points,
    depth_to_normal,
    BLOCK_WIDTH
)
from spirulae_splat.splat.background_sh import render_background_sh
from spirulae_splat.splat.sh import num_sh_bases, spherical_harmonics
from spirulae_splat.strategy import DefaultStrategy, MCMCStrategy
from spirulae_splat.splat._camera import _Camera

from nerfstudio.cameras.camera_optimizers import (CameraOptimizer,
                                                  CameraOptimizerConfig)
from nerfstudio.cameras.cameras import Cameras, CameraType
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

from nerfstudio.model_components.lib_bilagrid import BilateralGrid, color_correct, slice, total_variation_loss

from spirulae_splat.perf_timer import PerfTimer
timerr = PerfTimer("get_outputs", ema_tau=100)
timerl = PerfTimer("get_loss", ema_tau=100)


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
    weight = (1.0 / (d * d)) * torch.ones((1, 1, d, d), dtype=torch.float32, device=image.device)
    return F.conv2d(image.float().permute(2, 0, 1)[:, None, ...], weight, stride=d).squeeze(1).permute(1, 2, 0).to(image)


class SaturateKeepGradient(torch.autograd.Function):
    @staticmethod
    def forward(ctx, x, xmin, xmax):
        return torch.clip(x, xmin, xmax)
    @staticmethod
    def backward(ctx, v_x):
        return v_x, None, None

def saturate_keep_gradient(x, xmin=None, xmax=None):
    return SaturateKeepGradient.apply(x, xmin, xmax)


class MaskGradient(torch.autograd.Function):
    @staticmethod
    def forward(ctx, x, mask):
        ctx.mask = mask
        return x
    @staticmethod
    def backward(ctx, v_x):
        return v_x * ctx.mask, None


@dataclass
class SpirulaeModelConfig(ModelConfig):

    _target: Type = field(default_factory=lambda: SpirulaeModel)

    warmup_length: int = 500
    """period of steps where refinement is turned off"""
    stop_refine_at: int = 27000
    """period of steps where refinement is turned off"""
    refine_every: int = 100
    """period of steps where gaussians are culled and densified"""
    resolution_schedule: int = 3000
    """training starts at 1/d resolution, every n steps this is doubled"""
    background_color: Literal["random", "black", "white", "gray"] = "gray"
    """Whether to randomize the background color."""
    num_downscales: int = 2
    """at the beginning, resolution is 1/2^d, where d is this number"""
    use_mcmc: bool = False
    """use Markov-Chain Monte Carlo for gaussian control
        Disable per-pixel sorting if you use this"""
    random_init: bool = False
    """whether to initialize the positions uniformly randomly (not SFM points)"""
    num_random: int = 20000
    """Number of gaussians to initialize if random init is used"""
    random_scale: float = 1.0
    """Position standard deviation to initialize random gaussians"""
    ssim_lambda: float = 0.4  # 0.2
    """weight of ssim loss"""
    use_camera_optimizer: bool = False
    """Whether to use camera optimizer"""
    camera_optimizer: CameraOptimizerConfig = field(default_factory=lambda: CameraOptimizerConfig(mode="SO3xR3"))
    """Config of the camera optimizer to use"""
    kernel_radius: float = 1.0
    """Radius of the splatting kernel, 3.0 for Gaussian and 1.0 for polynomial"""
    loss_scale: float = 0.003
    """Scaling for loss values to normalize gradient"""
    target_absgrad: float = 5.5e-6
    """Will auto tune loss_scale to fit absgrad to this number, higher gives more splats"""
    relative_scale: Optional[float] = None
    """Manually set scale when a scene is poorly scaled by nerfstudio
        (e.g. Zip-NeRF dataset, very large-scale scenes across multiple street blocks)"""

    # classial control
    cull_alpha_thresh: float = 0.1
    """threshold of opacity for culling gaussians. One can set it to a lower value (e.g. 0.005) for higher quality."""
    cull_scale_thresh: float = 1.0
    """threshold of scale for culling huge gaussians"""
    cull_grad_thresh: float = 0.0  # 3e-4 | 1e-4 | 1e-5 | 0.0
    """threshold for culling gaussians with low visibility"""
    continue_cull_post_densification: bool = True
    """If True, continue to cull gaussians post refinement"""
    reset_alpha_every: int = 30
    """Every this many refinement steps, reset the alpha"""
    densify_xy_grad_thresh: float = 0.003  # 0.003 | 0.001
    """threshold of positional gradient norm for densifying gaussians"""
    # densify_ch_grad_thresh: float = 3e-6  # 5e-6 | 3e-6 | 1e-6
    # """threshold of cylindrical harmonics gradient norm for densifying gaussians"""
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
    mcmc_warmup_length: int = 500
    """start MCMC refinement at this number of steps"""
    mcmc_stop_refine_at: int = 25000
    """end MCMC refinement at this number of steps"""
    mcmc_cap_max: int = 100000
    """maximum number of splats for MCMC, dataset-specific tuning required"""
    mcmc_split_splats: bool = True
    """whether to split splats into two splats with different positions"""
    mcmc_noise_lr: float = 5e5
    """MCMC sampling noise learning rate"""
    mcmc_min_opacity: float = 0.005
    """minimum opacity for MCMC relocation"""

    # representation
    use_per_pixel_sorting: bool = False
    """enable per-pixel sorting, recommend disable this when using MCMC"""
    per_pixel_sorting_warmup: int = 2000
    """use per pixel sorting only after this number of steps"""
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
    adaptive_exposure_mode: Optional[Literal[
        "linear", "log_linear", "channel", "log_channel", "affine", "log_affine"
        ]] = "log_channel"
    """Adaptive exposure mode
        linear: gt ~ k * pred, with scalar k
        log_linear: log(gt) ~ k * log(pred) + b, with scalar k and b
        channel: linear but per channel
        log_channel: log_linear but per channel
        affine: gt ~ A * pred, with 3x3 matrix A
        log_affine: log(gt) ~ A * log(pred) + b, with 3x3 matrix A and vec3 b
       Additional coefficients (k,A,b) are optimized online using linear least squares
       May degrade performance at boundary when segmentation mask supervision is used, at least when tested on SAM2.1
    """
    adaptive_exposure_start_iter: int = 1000
    """Start adaptive exposure at this number of steps"""
    use_bilateral_grid: bool = False
    """If True, use bilateral grid to handle the ISP changes in the image space. This technique was introduced in the paper 'Bilateral Guided Radiance Field Processing' (https://bilarfpro.github.io/).
       Makes training much slower - TODO: fused bilagrid in CUDA"""
    grid_shape: Tuple[int, int, int] = (16, 16, 8)
    """Shape of the bilateral grid (X, Y, W)"""

    use_3dgs: bool = False
    """Use 3DGS instead of 2DGS, with limited support for existing features
        Tested with:
         - MCMC strategy
         - Depth supervision with alpha_supervision_weight (0.02,0.01) instead of (0.002,0)
         - Use bilateral grid, disable adaptive exposure
        """

    # regularization
    use_scale_regularization: bool = True
    """If enabled, a scale regularization introduced in PhysGauss (https://xpandora.github.io/PhysGaussian/) is used for reducing huge spikey gaussians."""
    max_gauss_ratio: float = 10.0
    """threshold of ratio of gaussian max to min scale before applying regularization
    loss from the PhysGaussian paper
    """
    depth_mode: Literal["mean", "median"] = "mean"
    """Depth rendering mode, use mean for stable training and median for high meshing resolution"""
    depth_reg_weight: float = 0.02
    """Weight for depth regularizer"""
    depth_reg_pairwise_factor: float = 1.0  # 0.7, 1.0
    """Factor of pairwise vs depth fitting depth regularization, 0 to 1"""
    depth_reg_warmup: int = 12000
    """warmup steps for depth regularizer"""
    normal_reg_weight: float = 0.04
    """Weight for normal regularizer"""
    normal_reg_warmup: int = 12000
    """warmup steps for normal regularizer"""
    reg_warmup_length: int = 4000
    """Warmup for depth and normal regularizers.
       only apply regularizers after this many steps."""
    alpha_loss_weight: int = 0.01
    """Weight for alpha, if mask is provided."""
    mcmc_opacity_reg: float = 0.01  # 0.01 in original paper
    """Opacity regularization from MCMC
       Lower usually gives more accurate geometry"""
    mcmc_scale_reg: float = 0.01  # 0.01 in original paper
    """Scale regularization from MCMC"""
    exposure_reg_image: float = 0.1
    """Between 0 and 1; For exposure regularization, include this fraction of L1 loss between GT image and image before exposure adjustment"""
    exposure_reg_param: float = 0.002
    """Make sure image look right in auto exposure mode,
       Use an L2 cost to match exposure parameter to the parameter at no exposure adjustment;
       (may be summed/averaged across a batch)"""

    # supervision using a foundation depth model
    # enable these by setting `depth_model` in data manager config
    depth_supervision_start_iter: int = 1000
    """Start using foundation model depth at this number of steps"""
    depth_distortion_depth_degree: int = -1  # 3
    """Hyperparameter for depth distortion model, controls depth embedding, see code for details
        Larger gives more parameters in depth distortion model, -1 to disable"""
    depth_distortion_uv_degree: int = -1  # 1
    """Hyperparameter for depth distortion model, controls image space embedding, see code for details
        Larger gives more parameters in depth distortion model, -1 to disable"""
    depth_supervision_weight: float = 0.5
    """Weight for depth supervision by comparing rendered depth with depth predicted by a foundation model"""
    normal_supervision_weight: float = 0.1
    """Weight for normal supervision by comparing gradient of rendered depth with gradient of depth predicted by a foundation model"""
    alpha_supervision_weight: float = 0.002
    """Weight for alpha supervision by rendered alpha with alpha predicted by a foundation model
        Useful for removing floaters from sky for outdoor scenes"""
    alpha_supervision_weight_under: float = 0.0
    """Similar to alpha_supervision_weight, but applies when renderer opacity is lower than reference opacity"""


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

        if self.config.use_3dgs:
            global fully_fused_projection
            global rasterize_to_pixels
            global isect_tiles
            global isect_offset_encode
            from spirulae_splat.splat.gsplat_3dgs_wrapper import (
                fully_fused_projection,
                rasterize_to_pixels,
                isect_tiles,
                isect_offset_encode,
            )

    def populate_modules(self):
        if self.seed_points is not None and not self.config.random_init:
            means = self.seed_points[0]
            if self.config.relative_scale is not None:
                means *= self.config.relative_scale
            means = torch.nn.Parameter(means)
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
        points = means.data[indices] - means.data[:, None, :]
        points = points.cpu().numpy()
        U, S, Vt = np.linalg.svd(points)
        num_points = means.shape[0]
        for i in range(num_points):
            sorted_indices = np.argsort(-S[i])
            S[i] = S[i][sorted_indices]
            Vt[i] = Vt[i][sorted_indices]
        scales = S
        if not self.config.use_3dgs:
            scales = S[:, :2]
        scales = np.log(1.5*scales+1e-8)
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

        gauss_params = {
            "means": means,
            "scales": scales,
            "quats": quats,
            "features_dc": features_dc,
            "features_sh": features_sh,
            "features_ch": features_ch,
            "opacities": opacities,
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

        if self.config.use_bilateral_grid:
            self.bil_grids = BilateralGrid(
                num=self.num_train_data,
                grid_X=self.config.grid_shape[0],
                grid_Y=self.config.grid_shape[1],
                grid_W=self.config.grid_shape[2],
            )

        self.crop_box: Optional[OrientedBox] = None
        if self.config.background_color == "random":
            self.background_color = torch.tensor(
                [0.1490, 0.1647, 0.2157]
            )  # This color is the same as the default background color in Viser. This would only affect the background color when rendering.
        elif self.config.background_color == "gray":
            self.background_color = torch.tensor([0.5, 0.5, 0.5])
        else:
            self.background_color = get_color(self.config.background_color)
        self.background_color = torch.nn.Parameter(self.background_color)
        if self.config.train_background_color:
            dim_sh = num_sh_bases(self.config.background_sh_degree)
            self.background_sh = torch.nn.Parameter(torch.zeros((dim_sh-1, 3)))
        else:
            self.background_sh = None

        self._train_batch_size = 1
        self._set_strategy()

    def _set_strategy(self):
        # Strategy for GS densification

        # MCMC mode
        if self.config.use_mcmc:
            self.strategy = MCMCStrategy(
                cap_max=self.config.mcmc_cap_max,
                noise_lr=self.config.mcmc_noise_lr,
                refine_start_iter=self.config.mcmc_warmup_length,
                refine_stop_iter=self.config.mcmc_stop_refine_at,
                refine_every=self.config.refine_every,
                min_opacity=self.config.mcmc_min_opacity,
                is_3dgs=self.config.use_3dgs,
            )
            self.strategy_state = self.strategy.initialize_state()
            return

        # classical mode
        reset_every = self.config.reset_alpha_every * self.config.refine_every
        pause_refine_after_reset = min(
            self.num_train_data//self._train_batch_size+self.config.refine_every,
            int(0.8*reset_every)
        )
        self.strategy = DefaultStrategy(
            prune_opa=self.config.cull_alpha_thresh,
            grow_grad2d=self.config.densify_xy_grad_thresh,
            prune_grad2d=self.config.cull_grad_thresh,
            grow_scale3d=self.config.densify_size_thresh,
            grow_scale2d=self.config.split_screen_size,
            prune_scale3d=self.config.cull_scale_thresh,
            prune_scale2d=self.config.cull_screen_size,
            refine_scale2d_stop_iter=self.config.stop_screen_size_at,
            refine_start_iter=self.config.warmup_length,
            refine_stop_iter=self.config.stop_refine_at,
            split_stop_iter=self.config.stop_split_at,
            reset_every=reset_every,
            refine_every=self.config.refine_every,
            pause_refine_after_reset=pause_refine_after_reset,
            absgrad=True,
            revised_opacity=False,
            verbose=True,
        )
        self.strategy_state = self.strategy.initialize_state(scene_scale=1.0)

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

    def load_state_dict(self, dict, **kwargs):  # type: ignore
        # resize the parameters to match the new number of points
        self.step = 30000
        if "means" in dict:
            # For backwards compatibility, we remap the names of parameters from
            # means->gauss_params.means since old checkpoints have that format
            for p in ["means", "scales", "quats",
                      "features_dc", "features_sh", "features_ch",
                      "opacities"]:
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

    def set_crop(self, crop_box: Optional[OrientedBox]):
        self.crop_box = crop_box

    def set_background(self, background_color: torch.Tensor):
        assert background_color.shape == (3,)
        # self.background_color = background_color

    def step_cb(self, optimizers: Optimizers, step):
        self.step = step
        self.optimizers = optimizers.optimizers

    def step_post_backward(self, step):
        assert step == self.step
        if self.config.use_mcmc:
            self.strategy.step_post_backward(
                params=self.gauss_params,
                optimizers=self.optimizers,
                state=self.strategy_state,
                step=self.step,
                info=self.info,
                lr=self.optimizers['means'].param_groups[0]['lr']
            )
        else:
            self.strategy.step_post_backward(
                params=self.gauss_params,
                optimizers=self.optimizers,
                state=self.strategy_state,
                step=self.step,
                info=self.info,
                packed=False,
            )
        if False:
            # adaptive densification threshold
            # TODO: do this again after resolution change
            if not hasattr(self, '_loss_scale_n'):
                self._loss_scale_n = 0
                self._loss_scale_item = 0.0
            step_unit = self.config.reset_alpha_every*self.config.refine_every
            if step_unit/3 <= step < step_unit:
                absgrad = self.info['means2d'].absgrad.mean().item()
                loss_scale = self.config.target_absgrad / (absgrad / self.config.loss_scale)
                self._loss_scale_n += 1
                self._loss_scale_item = self._loss_scale_item + (loss_scale-self._loss_scale_item)/self._loss_scale_n
                # self._loss_scale_item = self._loss_scale_item + (loss_scale-self._loss_scale_item)/(0.01*step_unit)
                if self._loss_scale_n > 100:
                    self.config.loss_scale = self._loss_scale_item
                # print(self._loss_scale_n, self._loss_scale_item)
        return
        # for splatfacto: around 1.5 and around 1e-6 @ 1k iters
        print(self.info['radii'].float().mean().item(),
              self.info['means2d'].absgrad.mean().item())

    def get_training_callbacks(
        self, training_callback_attributes: TrainingCallbackAttributes
    ) -> List[TrainingCallback]:
        cbs = []
        cbs.append(
            TrainingCallback(
                [TrainingCallbackLocation.BEFORE_TRAIN_ITERATION],
                self.step_cb,
                args=[training_callback_attributes.optimizers],
            )
        )
        cbs.append(
            TrainingCallback(
                [TrainingCallbackLocation.AFTER_TRAIN_ITERATION],
                self.step_post_backward,
            )
        )
        return cbs

    def get_gaussian_param_groups(self) -> Dict[str, List[Parameter]]:
        # Here we explicitly use the means, scales as parameters so that the user can override this function and
        # specify more if they want to add more optimizable params to gaussians.
        return {
            name: [self.gauss_params[name]]
            for name in [
                "means", "scales", "quats",
                "features_dc", "features_sh", "features_ch",
                "opacities"
            ]
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
        if self.config.use_bilateral_grid:
            gps["bilateral_grid"] = list(self.bil_grids.parameters())
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

    def _apply_bilateral_grid(self, rgb: torch.Tensor, cam_idx: int, H: int, W: int) -> torch.Tensor:
        # make xy grid
        grid_y, grid_x = torch.meshgrid(
            torch.linspace(0, 1.0, H, device=self.device),
            torch.linspace(0, 1.0, W, device=self.device),
            indexing="ij",
        )
        grid_xy = torch.stack([grid_x, grid_y], dim=-1).unsqueeze(0)

        rgb = torch.clip(rgb, 0.0, 1.0)

        out = slice(
            bil_grids=self.bil_grids,
            rgb=rgb,
            xy=grid_xy,
            grid_idx=torch.tensor(cam_idx, device=self.device, dtype=torch.long),
        )
        return out["rgb"]

    def get_background_image(self, camera: Cameras, ssplat_camera: _Camera):
        if not isinstance(camera, Cameras):
            print("Called get_background_image with not a camera")
            return {}
        assert camera.shape[0] == 1, "Only one camera at a time"
        camera_downscale = self._get_downscale_factor()
        camera.rescale_output_resolution(1 / camera_downscale)

        W, H = int(camera.width.item()), int(camera.height.item())

        sh_degree = self.config.background_sh_degree
        if not self.config.train_background_color or not (sh_degree > 0):
            return self.background_color.repeat(H, W, 1)

        rotation = camera.camera_to_worlds[0][:3, :3]
        sh_coeffs = torch.cat((self.background_color.unsqueeze(0), self.background_sh), dim=0)  # [(deg+1)^2, 3]
        return render_background_sh(ssplat_camera, rotation, sh_degree, sh_coeffs)

    @staticmethod
    def get_empty_outputs(width: int, height: int, background: torch.Tensor) -> Dict[str, Union[torch.Tensor, List]]:
        return {}

    def get_outputs(self, camera: Cameras) -> Dict[str, Union[torch.Tensor, List]]:
        """Takes in a camera and returns a dictionary of outputs."""
        if isinstance(camera, list):
            # TODO: support multi-GPU
            return [self.get_outputs(c) for c in camera]
        if not isinstance(camera, Cameras):
            print("Called get_outputs with not a camera")
            return {}
        assert camera.shape[0] == 1, "Only one camera at a time"

        timerr.start()

        if self.training and self.config.use_camera_optimizer:
            optimized_camera_to_world = self.camera_optimizer.apply_to_camera(camera)[0, ...]
        else:
            optimized_camera_to_world = camera.camera_to_worlds[0, ...]
        optimized_camera_to_world = optimized_camera_to_world.cpu()
        timerr.mark(".")  # 70us

        camera_downscale = self._get_downscale_factor()
        camera.rescale_output_resolution(1 / camera_downscale)
        W, H = int(camera.width.item()), int(camera.height.item())
        intrins = (camera.fx.item(), camera.fy.item(), camera.cx.item(), camera.cy.item())
        camera.rescale_output_resolution(camera_downscale)
        if camera.camera_type.item() == CameraType.FISHEYE.value:
            dist_coeffs = tuple(camera.distortion_params.cpu().numpy().flatten().tolist())
            ssplat_camera = _Camera(H, W, "OPENCV_FISHEYE", intrins, dist_coeffs)
        else:
            ssplat_camera = _Camera(H, W, "", intrins)
        timerr.mark(".")  # 400us

        R = optimized_camera_to_world[:3, :3]  # 3 x 3
        T = optimized_camera_to_world[:3, 3:4]  # 3 x 1
        if self.config.relative_scale is not None:
            T = T * self.config.relative_scale
        R = R * torch.tensor([[1.0, -1.0, -1.0]])
        R_inv = R.T
        T_inv = -R_inv @ T
        viewmat = torch.eye(4, dtype=optimized_camera_to_world.dtype)
        viewmat[:3, :3] = R_inv
        viewmat[:3, 3:4] = T_inv
        self.last_size = (H, W)
        timerr.mark("pre")  # 80us

        quats = F.normalize(self.quats, dim=-1)
        opacities = self.config.max_opacity * torch.sigmoid(self.opacities)
        timerr.mark("map")  # 300us-450us

        if self.config.sh_degree > 0:
            viewdirs = self.means.detach() - optimized_camera_to_world[:3, 3].cuda()  # (N, 3)
            n = min(self.step // self.config.sh_degree_interval, self.config.sh_degree)
            rgbs = spherical_harmonics(n, viewdirs, self.features_dc, self.features_sh)  # input unnormalized viewdirs
            rgbs = rgbs + 0.5
            # comment this - discourage below zero, make it compatible with existing viewers
            # rgbs = torch.relu(rgbs)
        else:
            rgbs = self.features_dc
        timerr.mark("sh")  # 200us-350us

        if not self.config.use_3dgs:
            # 2DGS rendering

            (
                positions, axes_u, axes_v,
                bounds, num_tiles_hit
            ) = project_gaussians(  # type: ignore
                self.means,
                torch.exp(self.scales),
                quats,
                viewmat[:3, :].cuda(),
                ssplat_camera,
            )  # type: ignore
            timerr.mark("project")  # 200us-350us

            use_per_pixel_sorting = (
                self.config.use_per_pixel_sorting and
                self.step >= self.config.per_pixel_sorting_warmup
            )
            if use_per_pixel_sorting:
                raster_indices = rasterize_gaussians_indices(
                    positions, axes_u, axes_v, opacities,
                    bounds, num_tiles_hit,
                    ssplat_camera
                )
                rasterize_depth = rasterize_gaussians_depth_sorted
                rasterize = rasterize_gaussians_sorted
                rasterize_simplified = rasterize_gaussians_simplified_sorted
            else:
                raster_indices = (bounds, num_tiles_hit)
                rasterize_depth = rasterize_gaussians_depth
                rasterize = rasterize_gaussians
                rasterize_simplified = rasterize_gaussians_simplified

            # slower but more capable two-pass rendering
            if False or (
                self.config.depth_reg_pairwise_factor < 1.0 or \
                self.config.ch_degree_r * (2*self.config.ch_degree_phi+1) > 0 or \
                self.config.depth_mode != "mean"
            ):
                assert not ssplat_camera.is_distorted(), \
                    "Fisheye camera is not supported for CH or median depth"

                depth_im_ref = rasterize_depth(
                    positions, axes_u, axes_v, opacities,
                    *((raster_indices[1],) if use_per_pixel_sorting else raster_indices),
                    ssplat_camera,
                    self.config.depth_mode
                )
                depth_im_ref = torch.where(
                    depth_im_ref > 0.0, depth_im_ref,
                    torch.amax(depth_im_ref).detach()
                ).contiguous()
                timerr.mark("depth")  # ?us

                # main rasterization
                ch_degree = self.step // self.config.ch_degree_interval
                (rgb, alpha, depth_im, normal_im, reg_depth) \
                = rasterize(  # type: ignore
                    positions, axes_u, axes_v, rgbs,
                    self.config.ch_degree_r, min(ch_degree, self.config.ch_degree_r),
                    self.config.ch_degree_phi, min(ch_degree, self.config.ch_degree_phi),
                    self.features_ch,
                    opacities,
                    depth_im_ref,
                    # background_color,
                    self.config.depth_reg_pairwise_factor,
                    *raster_indices,
                    ssplat_camera
                )  # type: ignore
                timerr.mark("render")  # ?us

            # fast one-pass rendering
            else:
                (rgb, alpha, depth_im, normal_im, reg_depth) \
                = rasterize_simplified(
                    positions, axes_u, axes_v, rgbs, opacities,
                    *raster_indices, ssplat_camera
                )
                depth_im_ref = torch.where(
                    alpha > 0.0, depth_im[..., :1] / alpha,
                    torch.amax(depth_im[..., :1]).detach()
                ).contiguous()
                timerr.mark("render")  # 750us-1200us
        
            radii = self.config.kernel_radius/3.0 * \
                      num_tiles_hit.unsqueeze(0)**0.5 / 2 * BLOCK_WIDTH
            means2d = positions

        else:
            # 3DGS rendering

            fx, fy, cx, cy = ssplat_camera.intrins
            Ks = torch.tensor([[fx, 0, cx], [0, fy, cy], [0, 0, 1]]).float()
            (
                camera_ids, gaussian_ids, radii,
                means2d, depths, conics, compensations,
            ) = fully_fused_projection(
                self.means, None, quats, torch.exp(self.scales),
                viewmat[None].cuda(), Ks[None].cuda(), W, H,
                eps2d=0.3, packed=True,
            )

            tile_size = 16
            tile_width = math.ceil(W / float(tile_size))
            tile_height = math.ceil(H / float(tile_size))
            tiles_per_gauss, isect_ids, flatten_ids = isect_tiles(
                means2d, radii, depths,
                tile_size, tile_width, tile_height,
                packed=True, n_cameras=1,
                camera_ids=camera_ids, gaussian_ids=gaussian_ids,
            )
            isect_offsets = isect_offset_encode(isect_ids, 1, tile_width, tile_height)

            rgbd = torch.cat((
                rgbs[gaussian_ids],
                depth_map(depths)[..., None]
            ), dim=-1)
            opacs = opacities[gaussian_ids].squeeze(-1)

            rgbd, alpha = rasterize_to_pixels(
                means2d, conics, rgbd, opacs,
                W, H, tile_size,
                isect_offsets, flatten_ids,
                backgrounds=None, packed=True, absgrad=True,
            )
            rgbd = rgbd[0]
            alpha = alpha[0]

            rgb = rgbd[..., :3]
            depth_im_ref = torch.where(
                alpha > 0.0, rgbd[..., 3:] / alpha,
                torch.amax(rgbd[..., 3:]).detach()
            ).contiguous()

            radii = radii.unsqueeze(0)

        if self.training:
            self.info = {
                "width": W,
                "height": H,
                "n_cameras": camera.shape[0],
                "radii": radii,
                "means2d": means2d,
            }
        timerr.mark("post")  # 100us-200us

        # blend with background
        background = self.background_color.reshape((1, 1, 3))
        if self.config.background_sh_degree > 0 or True:
            background = self.get_background_image(camera, ssplat_camera)
            rgb = rgb + (1.0 - alpha) * background
        timerr.mark("bkg")  # 300us-450us

        # apply bilateral grid
        if self.config.use_bilateral_grid and self.training:
            if camera.metadata is not None and "cam_idx" in camera.metadata:
                rgb = rgb.unsqueeze(0)
                rgb = self._apply_bilateral_grid(rgb, camera.metadata["cam_idx"], H, W)
                rgb = rgb.squeeze(0)

        # saturate pixel value
        else:
            if not self.training:
                rgb = torch.clip(rgb, 0.0, 1.0)
            else:
                rgb = torch.relu(rgb)
                # rgb = saturate_keep_gradient(rgb, 0.0, None)
        timerr.mark("bilagrid")  # 300us-450us

        # normal regularization
        depth_im_ref = depth_inv_map(depth_im_ref)
        depth_normal, alpha_diffused = depth_to_normal(depth_im_ref, ssplat_camera, None, True, alpha)
        if not self.config.use_3dgs:
            normal_im = F.normalize(normal_im, dim=-1)
            reg_normal = 1.0 - (depth_normal * normal_im).sum(-1, True)
        timerr.end("normal_reg")  # 600us-900us
        # -> ?us-?us median depth, 3000us-4500us one pass

        # pack outputs
        outputs = {
            "rgb": rgb,
            "depth": depth_im_ref,
            "depth_normal": 0.5+0.5*depth_normal,
        }
        if not self.config.use_3dgs:
            outputs["render_normal"] = 0.5+0.5*normal_im
            outputs["reg_depth"] = torch.sqrt(torch.relu(reg_depth*alpha)+1e-8)
            outputs["reg_normal"] = torch.sqrt(torch.relu(reg_normal*alpha)+1e-8) * alpha_diffused
        outputs["alpha"] = alpha
        outputs["background"] = background

        if not self.training and not self.config.use_3dgs and use_per_pixel_sorting:
            intersects = raster_indices[0].float() / raster_indices[1].shape[-1]
            outputs["num_intersects"] = intersects.reshape((H, W, 1)).repeat(1, 1, 3)

        if self.training:
            outputs["ssplat_camera"] = ssplat_camera

        # convert linear depth to ray depth, for correct gl_z_buf_depth in Viser
        if not self.training:
            undist_map = ssplat_camera.get_undist_map(always=True)
            distances = torch.sqrt((undist_map*undist_map).sum(-1, True) + 1.0)
            outputs["depth"] = outputs["depth"] * distances

        return outputs

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
        return {}  # TODO
        gt_rgb = self.composite_with_background(self.get_gt_img(batch["image"]), outputs["background"])
        metrics_dict = {}
        predicted_rgb = outputs["rgb"]
        metrics_dict["psnr"] = self.psnr(predicted_rgb, gt_rgb)

        metrics_dict["gaussian_count"] = self.num_points

        if self.config.use_camera_optimizer:
            self.camera_optimizer.get_metrics_dict(metrics_dict)
        return metrics_dict

    def depth_supervision(self, ref_depth, pred_depth, pred_alpha, _uv_cache={}):
        # details: https://github.com/harry7557558/Graphics/blob/master/mapping/relative_depth_matching/depth_fitting_01.ipynb

        # resize depth
        h, w, _ = pred_depth.shape
        if ref_depth.shape != pred_depth.shape:
            ref_depth = F.interpolate(
                ref_depth.reshape((1, 1, *ref_depth.shape[:2])),
                size=(h, w), mode='bilinear'
            ).reshape(pred_depth.shape)
        alpha = (ref_depth > 0.0).squeeze(-1)

        # generate UV embeddings
        if not self.config.depth_distortion_uv_degree >= 0:
            uv = None

        elif (h, w) not in _uv_cache:
            u = (torch.arange(w, dtype=torch.float32)+0.5)/w * 2.0 - 1.0
            v = (torch.arange(h, dtype=torch.float32)+0.5)/h * 2.0 - 1.0
            u = u[None, :].float().cuda().repeat((h, 1))
            v = v[:, None].float().cuda().repeat((1, w))
            u, v = u.flatten(), v.flatten()

            uv = []
            for i in range(self.config.depth_distortion_uv_degree+1):
                for j in range(self.config.depth_distortion_uv_degree+1):
                    uv.append(torch.cos(math.pi/2*i*u)*torch.cos(math.pi/2*j*v))
                    if j != 0:
                        uv.append(torch.cos(math.pi/2*i*u)*torch.sin(math.pi/2*j*v))
                    if i != 0:
                        uv.append(torch.sin(math.pi/2*i*u)*torch.cos(math.pi/2*j*v))
                    if i != 0 and j != 0:
                        uv.append(torch.sin(math.pi/2*i*u)*torch.sin(math.pi/2*j*v))

            uv = torch.stack(uv)
            _uv_cache[(h, w)] = uv

        else:
            uv = _uv_cache[(h, w)]

        # depth distortion fitting
        if self.config.depth_supervision_weight > 0.0 or self.config.normal_supervision_weight > 0.0:
            # normalize depths - inputs are metric linear depth
            x, y = ref_depth.flatten(), pred_depth.flatten()
            # x, y = pred_depth.flatten(), ref_depth.flatten()

            m = alpha.flatten().float()
            m_sum = alpha.sum().item()

            if (
                self.config.depth_distortion_depth_degree >= 0 and
                self.config.depth_distortion_uv_degree >= 0
            ):
                # normalize by mean and log
                x = x / ((x*m).sum() / m_sum)
                y = y / ((y*m).sum() / m_sum)
                x, y = depth_map(x), depth_map(y)

                # generate depth embeddings
                A = []
                for k in range(self.config.depth_distortion_depth_degree+1):
                    ed = x**k / math.factorial(k) * torch.exp(-x)
                    A.append(ed * uv)
                A = torch.concatenate(A)

                # linear least squares
                mat = (alpha.reshape(1,-1) * A) @ A.T
                vec = A @ (alpha.flatten() * y)
                c = torch.linalg.solve(mat, vec)
                corr_depth = torch.relu(c @ A)
                depth_diff = (corr_depth - y).reshape((h, w))

            else:
                # normalize log by mean and std
                x, y = torch.log(torch.relu(x) + 1e-6), torch.log(torch.relu(y) + 1e-6)
                x = x - (x*m).sum() / m_sum
                y = y - (y*m).sum() / m_sum
                x = x / torch.sqrt(((x*x)*m).sum() / m_sum)
                y = y / torch.sqrt(((y*y)*m).sum() / m_sum)

                # distortion disabled, just use it
                depth_diff = (x - y).reshape((h, w))
        # TODO: sometimes it's overfitting mean depth, need per-splat depth supervision

        # depth loss
        depth_loss = 0.0
        if self.config.depth_supervision_weight > 0.0:
            # depth_loss = ((depth_diff**2)*alpha).mean()
            depth_loss = (torch.abs(depth_diff)*alpha).mean()

        # normal loss
        normal_loss = 0.0
        if self.config.normal_supervision_weight > 0.0:
            normal_diff_x = (depth_diff[1:, :] - depth_diff[:-1, :]) * (alpha[1:, :] & alpha[:-1, :])
            normal_diff_y = (depth_diff[:, 1:] - depth_diff[:, :-1]) * (alpha[:, 1:] & alpha[:, :-1])
            # normal_loss = (normal_diff_x**2).mean() + (normal_diff_y**2).mean()
            normal_loss = normal_diff_x.abs().mean() + normal_diff_y.abs().mean()

        # alpha loss
        pred_alpha, alpha = pred_alpha.float().flatten(), alpha.float().flatten()
        alpha_loss = self.get_alpha_loss(pred_alpha, alpha)
        alpha_loss_over = 0.0
        if self.config.alpha_supervision_weight_under > 0.0:
            alpha_loss_over = self.get_alpha_loss(1.0-pred_alpha, 1.0-alpha)

        return (
            self.config.depth_supervision_weight * depth_loss,
            self.config.normal_supervision_weight * normal_loss,
            self.config.alpha_supervision_weight * alpha_loss +
                self.config.alpha_supervision_weight_under * alpha_loss_over
        )

    def exposure_correction(
            self,
            pred_img: torch.Tensor, gt_img: torch.Tensor,
            alpha_mask: Optional[torch.Tensor]=None
        ):
        exposure_param_reg = 0.0
        eeps = 0.5/255.0

        total_alpha = max(alpha_mask.sum().item(), 1.0) if alpha_mask is not None else 1.0
        def immean(img: torch.Tensor, dim: Optional[int]=None, keepdim=False):
            if alpha_mask is None:
                return img.mean(dim, keepdim)
            mask = alpha_mask.reshape(*img.shape[:-1], 1)
            return (img*mask).sum(dim, keepdim) / total_alpha

        # gt ~ k * pred
        if self.config.adaptive_exposure_mode == "linear":
            k = (immean(gt_img)+eeps) / (immean(pred_img)+eeps)
            pred_img_e = k * (pred_img+eeps) - eeps
            # print(k.item())
            exposure_param_reg = (k-1.0)**2

        # log(gt) ~ k * log(pred) + b
        elif self.config.adaptive_exposure_mode == "log_linear":
            x = torch.log(pred_img+eeps)
            y = torch.log(gt_img+eeps)
            x_mean, y_mean = immean(x), immean(y)
            x2_mean, xy_mean = immean(x**2), immean(x*y)
            m = (xy_mean - x_mean * y_mean) / (x2_mean - x_mean**2)
            b = y_mean - m * x_mean
            pred_img_e = torch.exp(m * x + b) - eeps
            # print(m.item(), b.item())
            exposure_param_reg = (m-1.0)**2 + b**2

        # gt ~ k * pred, per channel
        elif self.config.adaptive_exposure_mode == "channel":
            x = pred_img.reshape(-1, 3)
            y = gt_img.reshape(-1, 3)
            k = (immean(y, 0, True)+eeps) / (immean(x, 0, True)+eeps)
            pred_img_e = k * (pred_img+eeps) - eeps
            # print(k.detach())
            exposure_param_reg = torch.mean((k-1.0)**2)

        # log(gt) ~ k * log(pred) + b, per channel
        elif self.config.adaptive_exposure_mode == "log_channel":
            x = torch.log(pred_img + eeps).reshape(-1, 3)
            y = torch.log(gt_img + eeps).reshape(-1, 3)
            x_mean, y_mean = immean(x, 0, True), immean(y, 0, True)
            x2_mean, xy_mean = immean(x**2, 0, True), immean(x*y, 0, True)
            m = (xy_mean - x_mean * y_mean) / (x2_mean - x_mean**2)
            b = y_mean - m * x_mean
            pred_img_e = torch.exp(m * x + b).reshape(pred_img.shape) - eeps
            # print(m.detach(), b.detach())
            exposure_param_reg = torch.mean((m-1.0)**2 + b**2)

        # gt ~ A * pred
        elif self.config.adaptive_exposure_mode == "affine":
            x = (pred_img + eeps).reshape(-1, 3)
            y = (gt_img + eeps).reshape(-1, 3)
            if alpha_mask is not None:
                x = x * alpha_mask.reshape(-1, 1)
                y = y * alpha_mask.reshape(-1, 1)
            xTy, xTx = x.T @ y, x.T @ x
            xTy, xTx = xTy.cpu(), xTx.cpu()
            try:
                xTx = xTx + (1e-6*torch.trace(xTx).item()) * torch.eye(3)
                A = torch.linalg.inv(xTx) @ xTy
            except torch._C._LinAlgError:
                print("Warning: torch._C._LinAlgError")
                print(xTx)
                pred_img_e = pred_img.detach()
            else:
                A = A.to(x.device)
                pred_img_e = (x @ A).reshape(pred_img.shape) - eeps
                exposure_param_reg = torch.mean((A-torch.eye(3).to(A))**2)

        # log(gt) ~ A * log(pred) + b
        elif self.config.adaptive_exposure_mode == "log_affine":
            x = torch.log(pred_img + eeps).reshape(-1, 3)
            y = torch.log(gt_img + eeps).reshape(-1, 3)
            if alpha_mask is not None:
                x = x * alpha_mask.reshape(-1, 1)
                y = y * alpha_mask.reshape(-1, 1)
            mean_x = immean(x, 0, True)
            mean_y = immean(y, 0, True)
            centered_x = x - mean_x
            centered_y = y - mean_y
            xTy = centered_x.T @ centered_y
            xTx = centered_x.T @ centered_x
            xTy, xTx = xTy.cpu(), xTx.cpu()
            try:
                xTx = xTx + (1e-6*torch.trace(xTx).item()) * torch.eye(3)
                A = torch.linalg.inv(xTx) @ xTy
            except torch._C._LinAlgError:
                print("Warning: torch._C._LinAlgError")
                print(xTx)
                pred_img_e = pred_img.detach()
            else:
                A = A.to(x.device)
                B = mean_y - mean_x @ A
                pred_img_e = torch.exp(x @ A + B).reshape(pred_img.shape) - eeps
                exposure_param_reg = 0.5*(torch.mean((A-torch.eye(3).to(A))**2) + torch.mean(B**2))

        else:
            raise ValueError("Invalid adaptive_exposure_mode:", self.config.adaptive_exposure_mode)
        
        timerl.mark("exposure")  # 200us-500us

        if alpha_mask is not None:
            pred_img_e = pred_img + (pred_img_e-pred_img) * alpha_mask
        return pred_img_e, exposure_param_reg

    def get_alpha_loss(self, x, y):
        """Compute asymmetric loss for alpha, penalize only when first argument is greater than second argument

        Args:
            x: alpha output by the model
            y: reference alpha with same shape as x
        """
        # return torch.relu(x-y).mean()
        # return F.binary_cross_entropy(x, y, reduction="mean")
        return F.binary_cross_entropy(torch.fmax(x, y), y, reduction="mean")

    def get_loss_dict(self, outputs, batch, metrics_dict=None) -> Dict[str, torch.Tensor]:
        """Computes and returns the losses dict.

        Args:
            outputs: the output to compute loss dict to
            batch: ground truth batch corresponding to outputs
            metrics_dict: dictionary of metrics, some of which we can use for loss
        """
        if isinstance(outputs, list) and isinstance(batch, list):
            assert len(outputs) == len(batch)
            if self._train_batch_size != len(outputs):
                self._train_batch_size = len(outputs)
                self._set_strategy()
            losses = []
            for outputs_i, batch_i in zip(outputs, batch):
                losses.append(self.get_loss_dict(outputs_i, batch_i))
            loss = losses[0]
            for loss_i in losses[1:]:
                for key, value in loss_i.items():
                    loss[key] = loss[key] + value
            for key in loss:
                loss[key] = loss[key] / self._train_batch_size
            return loss

        # mask out of bound (e.g. fisheye circle)
        camera_mask = None
        ssplat_camera = outputs["ssplat_camera"]  # type: _Camera
        if ssplat_camera.is_distorted():
            undist_map = ssplat_camera.get_undist_map()
            camera_mask = torch.isfinite(undist_map.sum(-1, True))
            if not camera_mask.all():
                for key in ['rgb', 'depth', 'alpha', 'background']:
                    outputs[key] = MaskGradient.apply(outputs[key], camera_mask)

        timerl.start()
        gt_img_rgba = self.get_gt_img(batch["image"])
        gt_img = self.composite_with_background(gt_img_rgba, outputs["background"])
        pred_img = outputs["rgb"]

        # alpha channel for bounded objects - apply a cost on rendered alpha
        alpha_loss = 0.0
        if gt_img_rgba.shape[2] == 4:
            alpha = gt_img_rgba[..., -1].unsqueeze(-1)
            alpha_loss = alpha_loss + self.get_alpha_loss(outputs['alpha'], alpha)

        # separate mask for dynamic objects, text, etc.
        # simply don't consider it when evaluating loss, unless theres's no alpha channel, where a cost is applied
        mask = None
        if "mask" in batch:
            # batch["mask"] : [H, W, 1]
            mask = self._downscale_if_required(batch["mask"])
            mask = mask.float().to(self.device)
            assert mask.shape[:2] == gt_img.shape[:2] == pred_img.shape[:2]
            # can be little bit sketchy for the SSIM loss
            gt_img = torch.lerp(outputs["background"], gt_img, mask)
            # pred_img = torch.lerp(outputs["background"], pred_img, mask)
            if isinstance(alpha_loss, float) and alpha_loss == 0.0:
                alpha_loss = alpha_loss + self.get_alpha_loss(outputs['alpha'], mask)
        timerl.mark("alpha")  # ~100us

        # depth supervision
        depth_supervision_loss, normal_supervision_loss, alpha_supervision_loss = 0.0, 0.0, 0.0
        if "depth" in batch and self.step > self.config.depth_supervision_start_iter and \
            (self.config.depth_supervision_weight > 0.0 or
             self.config.normal_supervision_weight > 0.0 or
             self.config.alpha_supervision_weight > 0.0 or
             self.config.alpha_supervision_weight_under > 0.0):
            depth_supervision_loss, normal_supervision_loss, alpha_supervision_loss \
                = self.depth_supervision(batch["depth"].cuda(), outputs["depth"], outputs["alpha"])
            timerl.mark("depth")  # ?us

        # correct exposure
        pred_img_e = pred_img
        exposure_param_reg = 0.0
        if self.config.adaptive_exposure_mode is not None and \
            self.step > self.config.adaptive_exposure_start_iter:
            # mask, note that alpha_mask is defined in "mask out of bound" step
            alpha_mask = camera_mask
            if mask is not None and alpha_mask is not None:
                alpha_mask = alpha_mask & mask
            # call function
            pred_img_e, exposure_param_reg = self.exposure_correction(pred_img, gt_img, alpha_mask)

        # image loss
        if False:
            pred_img_e = saturate_keep_gradient(pred_img_e, 0.0, 1.0)
            pred_img = saturate_keep_gradient(pred_img, 0.0, 1.0)
        else:
            pred_img_e = torch.clip(pred_img_e, 0.0, 1.0)
            pred_img = torch.clip(pred_img, 0.0, 1.0)
        Ll1_e = torch.abs(gt_img - pred_img_e).mean()
        Ll1 = torch.abs(gt_img - pred_img).mean()
        simloss = torch.zeros_like(Ll1_e)
        if self.config.ssim_lambda > 0.0:
            gt_img_bchw = gt_img.permute(2, 0, 1).unsqueeze(0)
            pred_img_bchw = pred_img_e.permute(2, 0, 1).unsqueeze(0).contiguous()
            # simloss = 1 - self.ssim(pred_img_bchw, gt_img_bchw)
            # simloss = 1 - fused_ssim(pred_img_bchw, gt_img_bchw, padding="valid")
            simloss = 1 - fused_ssim(pred_img_bchw, gt_img_bchw, padding="same")
        timerl.mark("image")  # 750us-3000us, 100us-300us without ssim

        # scale regularization
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
        timerl.mark("scale")  # <100us

        # depth and normal regularizers
        depth_reg, normal_reg = 0.0, 0.0
        if not self.config.use_3dgs:
            alpha = outputs['alpha']
            reg_depth = outputs["reg_depth"]
            reg_normal = outputs["reg_normal"]
            weight_depth_reg = self.config.depth_reg_weight * \
                min(self.step / max(self.config.depth_reg_warmup, 1), 1)
            weight_normal_reg = self.config.normal_reg_weight * \
                min(self.step / max(self.config.normal_reg_warmup, 1), 1)
            if self.step < self.config.reg_warmup_length:
                weight_depth_reg, weight_normal_reg = 0.0, 0.0
            depth_reg = weight_depth_reg * reg_depth.mean()
            normal_reg = weight_normal_reg * reg_normal[1:-1, 1:-1].mean()
            timerl.mark("reg")  # ~100us

        # MCMC regularizers
        use_mcmc = self.config.use_mcmc
        mcmc_opacity_reg = use_mcmc * self.config.mcmc_opacity_reg * \
            torch.abs(torch.sigmoid(self.opacities)).mean()
        mcmc_scale_reg = use_mcmc * self.config.mcmc_scale_reg * \
            torch.abs(torch.exp(self.scales)).mean()

        # regularizations for parameters
        quat_norm = self.quats.norm(dim=-1)
        quat_norm_reg = 0.1 * (quat_norm-1.0-torch.log(quat_norm)).mean()
        timerl.mark("mcmc")  # ~100us

        loss_dict = {
            "main_loss": torch.lerp(torch.lerp(Ll1_e, simloss, self.config.ssim_lambda), Ll1, self.config.exposure_reg_image),
            "depth_ref_loss": depth_supervision_loss,
            "normal_ref_loss": normal_supervision_loss,
            "alpha_ref_loss": alpha_supervision_loss,
            "alpha_loss": self.config.alpha_loss_weight * alpha_loss,
            "scale_reg": scale_reg,
            "depth_reg": depth_reg,
            "normal_reg": normal_reg,
            "quat_reg": quat_norm_reg,
            'mcmc_opacity_reg': mcmc_opacity_reg,
            'mcmc_scale_reg': mcmc_scale_reg,
            "exposure_param_reg": self.config.exposure_reg_param * exposure_param_reg,
        }
        for key, value in loss_dict.items():
            loss_dict[key] = self.config.loss_scale * value
        timerl.mark("scale")  # <100us

        if self.config.use_bilateral_grid:
            loss_dict["tv_loss"] = 10 * total_variation_loss(self.bil_grids.grids)

        if self.training and self.config.use_camera_optimizer:
            # Add loss from camera optimizer
            self.camera_optimizer.get_loss_dict(loss_dict)
        timerl.end("camera")  # <100us -> 1600us-4200us, 900us-1400us without ssim

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
