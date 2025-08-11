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
    depth_to_normal,
    BLOCK_WIDTH,
)
from spirulae_splat.splat.utils import resize_image
from spirulae_splat.splat.background_sh import render_background_sh
from spirulae_splat.splat.sh import num_sh_bases, spherical_harmonics
from spirulae_splat.strategy import DefaultStrategy, MCMCStrategy
from spirulae_splat.splat._camera import _Camera

from spirulae_splat.modules.training_losses import SplatTrainingLosses

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



class SaturateKeepGradient(torch.autograd.Function):
    @staticmethod
    def forward(ctx, x, xmin, xmax):
        return torch.clip(x, xmin, xmax)
    @staticmethod
    def backward(ctx, v_x):
        return v_x, None, None

def saturate_keep_gradient(x, xmin=None, xmax=None):
    return SaturateKeepGradient.apply(x, xmin, xmax)



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
    ssim_lambda: float = 0.4
    """weight of ssim loss; 0.2 for optimal PSNR, higher for better visual quality"""
    ssim_warmup: int = 0
    """warmup of ssim loss"""
    use_camera_optimizer: bool = False
    """Whether to use camera optimizer"""
    camera_optimizer: CameraOptimizerConfig = field(default_factory=lambda: CameraOptimizerConfig(mode="SO3xR3"))
    """Config of the camera optimizer to use"""
    kernel_radius: float = 1.0
    """Radius of the splatting kernel, 3.0 for Gaussian and 1.0 for polynomial"""
    relative_scale: Optional[float] = None
    """Manually set scale when a scene is poorly scaled by nerfstudio
        (e.g. Zip-NeRF dataset, very large-scale scenes across multiple street blocks)"""
    compute_depth_normal: bool = True
    """Compute normal from depth. Required for 2DGS and supervision. Disabling this can reduce VRAM usage and speed up training."""

    # classial control
    cull_alpha_thresh: float = 0.1
    """threshold of opacity for culling gaussians. One can set it to a lower value (e.g. 0.005) for higher quality."""
    cull_scale_thresh: float = 1.0
    """threshold of world scale for culling huge gaussians"""
    split_scale_thresh: float = float('inf')
    """threshold of world scale for splitting huge gaussians"""
    cull_grad_thresh: float = 0.0  # 3e-4 | 1e-4 | 1e-5 | 0.0
    """threshold for culling gaussians with low visibility"""
    continue_cull_post_densification: bool = True
    """If True, continue to cull gaussians post refinement"""
    reset_alpha_every: int = 30
    """Every this many refinement steps, reset the alpha"""
    densify_xy_grad_thresh: float = 0.0005
    """threshold of positional gradient norm for densifying gaussians"""
    densify_size_thresh: float = 0.01
    """below this size, gaussians are *duplicated*, otherwise split"""
    stop_split_at: int = 15000
    """stop splitting at this step"""
    cull_screen_size: float = 0.15
    """if a gaussian is more than this fraction of screen space, cull it"""
    split_screen_size: float = 0.05
    """if a gaussian is more than this fraction of screen space, split it"""
    stop_screen_size_at: int = 4000
    """stop culling/splitting at this step WRT screen size of gaussians"""

    # MCMC control
    mcmc_warmup_length: int = 500
    """start MCMC refinement at this number of steps"""
    mcmc_cap_max: int = 100000
    """maximum number of splats for MCMC, dataset-specific tuning required"""
    mcmc_noise_lr: float = 5e5
    """MCMC sampling noise learning rate"""
    mcmc_min_opacity: float = 0.005
    """minimum opacity for MCMC relocation"""
    mcmc_prob_grad_weight: float = 0.0
    """weight of position gradient used in sampling Gaussians to relocate/add to, uses only opacity if 0 and only gradient of 1"""
    relocate_screen_size: float = float('inf')
    """if a gaussian is more than this fraction of screen space, relocate it
        Useful for fisheye with 3DGUT, may drop PSNR for conventional cameras
        For likely better quality, use mcmc_max_screen_size instead"""
    mcmc_max_screen_size: float = 0.15
    """if a gaussian is more than this fraction of screen space, clip scale and increase opacity
        Intended to be an MCMC-friendly alternative of relocate_screen_size"""
    mcmc_max_world_size: float = float('inf')
    """if a gaussian is more than this of world space, clip scale
        Useful if you see huge floaters at a distance in large indoor space"""

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
    adaptive_exposure_warmup: int = 1000
    """Start adaptive exposure at this number of steps"""
    use_bilateral_grid: bool = False
    """If True, use bilateral grid to handle the ISP changes in the image space. This technique was introduced in the paper 'Bilateral Guided Radiance Field Processing' (https://bilarfpro.github.io/).
       Makes training much slower - TODO: fused bilagrid in CUDA"""
    bilagrid_shape: Tuple[int, int, int] = (16, 16, 8)
    """Shape of the bilateral grid (X, Y, W)"""

    use_3dgs: bool = False
    """Use 3DGS instead of 2DGS, with limited support for existing features
        Tested with MCMC strategy, 3.0 kernel radius, 1.0 loss scale
    """

    # regularization
    scale_regularization_weight: float = 0.0
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
    """warmup steps for depth regularizer, regularization weight ramps up"""
    normal_reg_weight: float = 0.04
    """Weight for normal regularizer"""
    normal_reg_warmup: int = 12000
    """warmup steps for normal regularizer, regularization weight ramps up"""
    alpha_reg_weight: float = 0.025
    """Weight for alpha regularizer (encourage alpha to go to either 0 or 1)"""
    alpha_reg_warmup: int = 12000
    """warmup steps for alpha regularizer, regularization weight ramps up"""
    reg_warmup_length: int = 4000
    """Warmup steps for depth, normal, and alpha regularizers.
       only apply regularizers after this many steps."""
    alpha_loss_weight: float = 0.01
    """Weight for alpha, if mask is provided.
       Set this to 0.0 to use masks to ignore distractors (e.g. people and cars, area outside fisheye circle)
       Set this to a positive value to remove background (e.g. sky, background around centered object)
       See scripts/SAM2-GUI for what I use to generate masks"""
    mcmc_opacity_reg: float = 0.01  # 0.01 in original paper
    """Opacity regularization from MCMC
       Lower usually gives more accurate geometry"""
    mcmc_scale_reg: float = 0.01  # 0.01 in original paper
    """Scale regularization from MCMC"""
    erank_reg: float = 0.01
    """erank regularization weight, for 3DGS only -
        see *Effective Rank Analysis and Regularization for Enhanced 3D Gaussian Splatting, Hyung et al.*"""
    erank_reg_s3: float = 0.05
    """erank regularization weight for smallest dimension, for 3DGS only"""
    erank_reg_warmup: int = 1000
    """only apply erank regularization after this many steps"""
    exposure_reg_image: float = 0.1
    """Between 0 and 1; For exposure regularization, include this fraction of L1 loss between GT image and image before exposure adjustment"""
    exposure_reg_param: float = 0.002
    """Make sure image look right in auto exposure mode,
       Use an L2 cost to match exposure parameter to the parameter at no exposure adjustment;
       (may be summed/averaged across a batch)"""
    randomize_background: bool = False
    """Randomize background color during training to discourage transparency. This will override background model."""

    # supervision using a foundation depth model
    # enable these by setting `depth_model` in data manager config
    supervision_warmup: int = 0
    """Start using foundation model depth at this number of steps"""
    depth_distortion_depth_degree: int = -1  # 3
    """Hyperparameter for depth distortion model, controls depth embedding, see code for details
        Larger gives more parameters in depth distortion model, -1 to disable"""
    depth_distortion_uv_degree: int = -1  # 1
    """Hyperparameter for depth distortion model, controls image space embedding, see code for details
        Larger gives more parameters in depth distortion model, -1 to disable"""
    depth_supervision_weight: float = 0.01
    """Weight for depth supervision by comparing rendered depth with depth predicted by a foundation model"""
    normal_supervision_weight: float = 0.01
    """Weight for normal supervision by comparing normal from rendered depth with normal from depth predicted by a foundation model"""
    alpha_supervision_weight: float = 0.005
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
            global gsplat
            import gsplat

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
        Vt[:,:,2] *= np.linalg.det(Vt)[:,None]
        num_points = means.shape[0]
        for i in range(num_points):
            if self.config.use_3dgs or True:
                S[i] = np.prod(S[i])**(1/len(S[i])) * np.ones(S[i].shape)
                continue
            sorted_indices = np.argsort(-S[i])
            S[i] = S[i][sorted_indices]
            Vt[i] = Vt[i][sorted_indices]
        scales = S
        if not self.config.use_3dgs:
            scales = S[:, :2]
        scales = np.log(1.5*scales/self.config.kernel_radius+1e-8)
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
            shs = torch.zeros((self.seed_points[1].shape[0], dim_sh, 3)).float()
            seed_color = self.seed_points[1] / 255
            if self.config.sh_degree > 0:
                shs[:, 0, :3] = RGB2SH(seed_color)
                shs[:, 1:, 3:] = 0.0
            else:
                shs[:, 0, :3] = seed_color
            if dim_ch > 0:
                shs *= 2
            chs = torch.zeros((self.seed_points[1].shape[0], dim_ch, 3)).float()
            features_dc = torch.nn.Parameter(shs[:, 0, :].contiguous())
            features_sh = torch.nn.Parameter(shs[:, 1:, :].contiguous())
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

        self.training_losses = SplatTrainingLosses(self.config, self.num_train_data)

        if not self.config.compute_depth_normal and not self.config.use_3dgs:
            raise ValueError("--pipeline.model.compute_depth_normal must be enabled for 2DGS")

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
            current_num = len(self.means)
            final_num = self.config.mcmc_cap_max
            warmup_steps = (self.config.resolution_schedule * (self.config.num_downscales+1)
                            - self.config.mcmc_warmup_length) / self.config.refine_every
            grow_factor = max(final_num / current_num, 1.0) ** (1/warmup_steps)
            # grow_factor = min(grow_factor, 1.05)
            # grow_factor = max(grow_factor, 1.05)
            grow_factor = 1.05

            self.strategy = MCMCStrategy(
                cap_max=self.config.mcmc_cap_max,
                noise_lr=self.config.mcmc_noise_lr,
                refine_start_iter=self.config.mcmc_warmup_length,
                refine_stop_iter=self.config.stop_refine_at,
                refine_every=self.config.refine_every,
                grow_factor=grow_factor,
                min_opacity=self.config.mcmc_min_opacity,
                prob_grad_weight=self.config.mcmc_prob_grad_weight,
                is_3dgs=self.config.use_3dgs,
                relocate_scale2d=self.config.relocate_screen_size,
                max_scale2d=self.config.mcmc_max_screen_size,
                max_scale3d=self.config.mcmc_max_world_size,
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
            split_scale3d=self.config.split_scale_thresh,
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
            gps["bilateral_grid"] = list(self.training_losses.bil_grids.parameters())
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

    def get_background_image(self, camera: Cameras, ssplat_camera: _Camera):
        if not isinstance(camera, Cameras):
            print("Called get_background_image with not a camera")
            return {}
        assert camera.shape[0] == 1, "Only one camera at a time"
        camera_downscale = self._get_downscale_factor()
        camera.rescale_output_resolution(1 / camera_downscale)

        W, H = int(camera.width.item()), int(camera.height.item())

        if self.config.randomize_background:
            # return torch.rand_like(self.background_color).repeat(H, W, 1)
            return torch.rand((H, W, 3), device=self.background_color.device)

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
        device = self.means.device

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
            ssplat_camera = _Camera(H, W, "OPENCV_FISHEYE", intrins, dist_coeffs, device=device)
        else:
            ssplat_camera = _Camera(H, W, "", intrins, device=device)
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
            viewdirs = self.means.detach() - T.squeeze(-1).to(device)  # (N, 3)
            n = min(self.step // self.config.sh_degree_interval, self.config.sh_degree)
            rgbs = spherical_harmonics(n, viewdirs, self.features_dc, self.features_sh)  # input unnormalized viewdirs
            rgbs = rgbs + 0.5
            # comment this - discourage below zero, make it compatible with existing viewers
            # rgbs = torch.relu(rgbs)
        else:
            rgbs = self.features_dc
        timerr.mark("sh")  # 200us-350us

        max_depth_scale = 2.0 if self.training else 1.0

        if not self.config.use_3dgs:
            # 2DGS rendering

            (
                positions, axes_u, axes_v,
                bounds, num_tiles_hit
            ) = project_gaussians(  # type: ignore
                self.means,
                torch.exp(self.scales),
                quats,
                viewmat[:3, :].to(device),
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
                    max_depth_scale*torch.amax(depth_im_ref).detach()
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
                    max_depth_scale*torch.amax(depth_im[..., :1]).detach()
                ).contiguous()
                timerr.mark("render")  # 750us-1200us
        
            radii = torch.stack([
                bounds[:, 2] - bounds[:, 0],
                bounds[:, 3] - bounds[:, 1],
            ]).T.unsqueeze(0)
            means2d = positions
            depths = positions[..., 2]

        else:
            # 3DGS rendering

            fx, fy, cx, cy = ssplat_camera.intrins
            Ks = torch.tensor([[fx, 0, cx], [0, fy, cy], [0, 0, 1]]).float()
            
            kwargs = {}
            is_fisheye = ssplat_camera.model == "OPENCV_FISHEYE"
            if is_fisheye:
                if not self.config.use_mcmc:
                    raise ValueError("3DGS training with fisheye camera is currently only supported for MCMC.")
                dist_coeffs = torch.tensor(ssplat_camera.dist_coeffs).float().to(device)
                if len(dist_coeffs) == 4:
                    kwargs['radial_coeffs'] = dist_coeffs[None]
                elif len(dist_coeffs) == 6:
                    kwargs["radial_coeffs"] = dist_coeffs[:4][None]
                    kwargs["tangential_coeffs"] = dist_coeffs[4:][None]
                    # TODO: make sure GSplat 3DGUT actually supports this??
                else:
                    raise ValueError("Only support fisheye with 4 or 6 distortion coefficients")

            rgbd, alpha, meta = gsplat.rasterization(
                means=self.means,
                quats=quats,
                scales=torch.exp(self.scales),
                opacities=opacities.squeeze(-1),
                colors=rgbs,
                viewmats=viewmat[None].to(device),  # [C, 4, 4]
                Ks=Ks[None].to(device),  # [C, 3, 3]
                width=W,
                height=H,
                packed=False,
                absgrad=(not self.config.use_mcmc),
                sparse_grad=False,
                rasterize_mode="classic",
                distributed=False,
                camera_model=["pinhole", "fisheye"][is_fisheye],
                with_ut=is_fisheye,
                with_eval3d=is_fisheye,
                render_mode="RGB+D",
                **kwargs,
            )
            rgbd = rgbd[0]
            alpha = alpha[0]

            rgb = rgbd[..., :3]
            depth_im_ref = torch.where(
                alpha > 0.0, rgbd[..., 3:] / alpha,
                max_depth_scale*torch.amax(rgbd[..., 3:]).detach()
            ).contiguous()

            means2d = meta["means2d"]
            if self.config.use_mcmc:
                means2d = self.means
            radii = meta["radii"]
            depths = meta["depths"]

        if self.training:
            self.info = {
                "width": W,
                "height": H,
                "n_cameras": camera.shape[0],
                "radii": radii,
                "means2d": means2d,
                "depths": depths,
            }
        timerr.mark("post")  # 100us-200us

        # blend with background
        background = self.get_background_image(camera, ssplat_camera)
        rgb = rgb + (1.0 - alpha) * background
        timerr.mark("bkg")  # 300us-450us

        # apply bilateral grid
        if self.config.use_bilateral_grid and self.training:
            if camera.metadata is not None and "cam_idx" in camera.metadata:
                rgb = rgb.unsqueeze(0)
                rgb = self.training_losses.apply_bilateral_grid(rgb, camera.metadata["cam_idx"], H, W)
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
        if not self.config.use_3dgs:
            depth_im_ref = depth_inv_map(depth_im_ref)
        if self.config.compute_depth_normal or not self.training:
            depth_normal, alpha_diffused = depth_to_normal(depth_im_ref, ssplat_camera, None, True, alpha)
        if not self.config.use_3dgs:
            normal_im = F.normalize(normal_im, dim=-1)
            reg_normal = 1.0 - (depth_normal * normal_im).sum(-1, True)
            reg_normal = torch.nan_to_num(reg_normal, 0.0, 0.0, 0.0)
        timerr.end("normal_reg")  # 600us-900us
        # -> ?us-?us median depth, 3000us-4500us one pass

        # pack outputs
        outputs = {
            "rgb": rgb,
            "depth": depth_im_ref,
        }
        if self.config.compute_depth_normal or not self.training:
            outputs["depth_normal"] = depth_normal
        if not self.config.use_3dgs:
            outputs["render_normal"] = 0.5+0.5*normal_im
            outputs["reg_depth"] = torch.sqrt(torch.relu(reg_depth*alpha)+1e-8)
            outputs["reg_normal"] = torch.sqrt(torch.relu(reg_normal*alpha_diffused)+1e-8)
        outputs["alpha"] = alpha
        outputs["background"] = background

        if not self.training and not self.config.use_3dgs and use_per_pixel_sorting:
            intersects = raster_indices[0].float() / raster_indices[1].shape[-1]
            outputs["num_intersects"] = intersects.reshape((H, W, 1)).repeat(1, 1, 3)
        if not self.training:
            outputs["alpha"] = outputs["alpha"].reshape((H, W, 1)).repeat(1, 1, 3)

        if self.training:
            outputs["ssplat_camera"] = ssplat_camera

        # convert linear depth to ray depth, for correct gl_z_buf_depth in Viser
        if not self.training:
            undist_map = ssplat_camera.get_undist_map(always=True)
            distances = torch.sqrt((undist_map*undist_map).sum(-1, True) + 1.0)
            outputs["depth"] = outputs["depth"] * distances
            outputs["depth"] = torch.clip(outputs["depth"], max=torch.quantile(outputs["depth"], 0.99))
            if self.config.relative_scale is not None:
                outputs["depth"] /= self.config.relative_scale
            if "depth_normal" in outputs:
                outputs["depth_normal"] = 0.5+0.5*outputs["depth_normal"]

        return outputs

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

    def get_loss_dict(self, outputs, batch, metrics_dict=None) -> Dict[str, torch.Tensor]:
        """Computes and returns the losses dict.

        Args:
            outputs: the output to compute loss dict to
            batch: ground truth batch corresponding to outputs
            metrics_dict: dictionary of metrics, some of which we can use for loss
        """

        # Per-image losses
        if isinstance(outputs, list) and isinstance(batch, list):
            assert len(outputs) == len(batch)
            if self._train_batch_size != len(outputs):
                self._train_batch_size = len(outputs)
                self._set_strategy()
            losses = []
            for outputs_i, batch_i in zip(outputs, batch):
                losses.append(self.training_losses(self.step, batch_i, outputs_i))
            loss_dict = losses[0]
            for loss_i in losses[1:]:
                for key, value in loss_i.items():
                    loss_dict[key] = loss_dict[key] + value
            for key in loss_dict:
                loss_dict[key] = loss_dict[key] / self._train_batch_size

        else:
            loss_dict = self.training_losses(self.step, batch, outputs)

        # Per-splat losses
        self.training_losses.get_static_losses(
            self.step,
            self.quats, self.scales, self.opacities,
            loss_dict
        )

        # Camera optimizer loss
        if self.training and self.config.use_camera_optimizer:
            self.camera_optimizer.get_loss_dict(loss_dict)

        self.print_loss_dict(loss_dict)
        return loss_dict

    def print_loss_dict(self, losses: Dict[str, torch.Tensor], _max_vals={}):

        # get VRAM usage (only supports single GPU at this time)
        free, total = torch.cuda.mem_get_info(self.means.device)
        used = (total - free) / 1024**3
        used_percentage = (1 - free/total)*100
        mem_stats = f"{used:.2f}\N{ZERO WIDTH SPACE}GB {used_percentage:.0f}%"

        def fmt(key: str, s: float, decimals=None) -> str:
            if s == 0.0:
                return '~'

            l = float(losses[key]) / s
            if not math.isfinite(l):  # not finite
                return str(l)

            if key not in _max_vals or self.step % 1000 == 0:
                _max_vals[key] = abs(l)
            _max_vals[key] = max(_max_vals[key], abs(l))
            if _max_vals[key] == 0.0:
                return '~'

            if decimals is None:
                decimals = int(max(-math.log10(0.001*_max_vals[key]), 0))
            s = f"{{:.{decimals}f}}".format(l)
            if s.startswith('0.'):
                s = s[1:]
            return s

        mcmc_reg = (self.config.use_mcmc and self.step < self.config.stop_refine_at)
        chunks = [
            f"[N] {len(self.opacities)} {mem_stats}",
            f"[C] {fmt('main_loss', 1.0)} "
            f"{fmt('alpha_loss', self.config.alpha_loss_weight)} "
            f"{fmt('psnr', 1.0, 2)} "
            f"{fmt('ssim', 1.0, 3)}",
            f"[S] {fmt('depth_ref_loss', self.config.depth_supervision_weight)} "
            f"{fmt('normal_ref_loss', self.config.normal_supervision_weight)} "
            f"{fmt('alpha_ref_loss', self.config.alpha_supervision_weight)}",
            f"[G] {fmt('depth_reg', self.training_losses.get_2dgs_reg_weights()[0])} "
            f"{fmt('normal_reg', self.training_losses.get_2dgs_reg_weights()[1])} "
            f"{fmt('alpha_reg', self.training_losses.get_alpha_reg_weight())}",
            f"[M] {fmt('mcmc_opacity_reg', self.config.mcmc_opacity_reg * mcmc_reg)} "
            f"{fmt('mcmc_scale_reg', self.config.mcmc_scale_reg * mcmc_reg)}",
            f"[R] {fmt('erank_reg', max(self.config.erank_reg_s3, self.config.erank_reg))} "
            f"{fmt('scale_reg', self.config.scale_regularization_weight)}",
            f"[E] {fmt('tv_loss', 10.0)} "
            f"{fmt('exposure_param_reg', self.config.exposure_reg_param)}",
        ]
        chunks = [c for c in chunks if any(char.isdigit() for char in c)]
        CONSOLE.print(' '.join(chunks).replace('\n', '') + "    ", end="\r")

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
        return {}, {}   # TODO
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
