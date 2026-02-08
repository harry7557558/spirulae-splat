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
3DGS implementation that combines many recent advancements.
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

from spirulae_splat.splat._torch_impl import quat_to_rotmat
from spirulae_splat.splat import (
    rasterization,
    depth_to_normal,
    BLOCK_WIDTH,
)
from spirulae_splat.splat.utils import resize_image, _TORCH_COMPILE_ARGS
from spirulae_splat.splat.background_sh import render_background_sh
from spirulae_splat.splat.sh import num_sh_bases, spherical_harmonics
from spirulae_splat.strategy import DefaultStrategy, MCMCStrategy, OpaqueStrategy, SVRasterStrategy
from spirulae_splat.splat._camera import _Camera

from spirulae_splat.modules.training_losses import SplatTrainingLosses
from spirulae_splat.splat.cuda._wrapper_per_pixel import blend_background, linear_rgb_to_srgb, get_color_transform_matrix
from spirulae_splat.splat.cuda import (
    svhash_create_initial_volume,
    svhash_get_voxels
)

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

    primitive: Literal["3dgs", "mip", "3dgut", "3dgut_sv", "opaque_triangle", "voxel"] = "3dgut"
    """Splat primitive to use"""
    fit: Literal["rgb", "depth", "normal"] = "rgb"
    """Fit RGB image by default, fit depth/normal for geometry reconstruction
        Will apply the same as RGB for depth and normal (e.g. bilagrid, SSIM, unless you disable them)"""
    fit_normal_depth_factor: float = 0.2
    """Probability to use depth normal instead of rendered normal in normal fitting mode"""
    fit_normal_depth_warmup_steps: int = 3000
    """Warmup for probability to use depth normal instead of rendered normal in normal fitting mode"""

    num_iterations: int = 30000
    """number of training iterations, should be consistent with --max_num_iterations"""
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
    num_downscales: int = 0
    """at the beginning, resolution is 1/2^d, where d is this number"""
    use_mcmc: bool = True
    """use Markov-Chain Monte Carlo for gaussian control
        Disable per-pixel sorting if you use this"""
    random_init: bool = False
    """whether to initialize the positions uniformly randomly (not SFM points)"""
    num_random: int = 200000
    """Number of gaussians to initialize if random init is used"""
    random_scale: float = 1.0
    """Position standard deviation to initialize random gaussians"""
    ssim_lambda: float = 0.4
    """Weight of ssim loss; 0.2 for optimal PSNR, higher for better visual quality"""
    lpips_lambda: float = 0.0
    """Weight of lpips loss for better perceptual quality; Note that this can make training much slower"""
    use_camera_optimizer: bool = False
    """Whether to use camera optimizer
        Note: this only works well in patch batching mode"""
    camera_optimizer: CameraOptimizerConfig = field(default_factory=lambda: CameraOptimizerConfig(mode="SO3xR3"))
    """Config of the camera optimizer to use"""
    kernel_radius: float = 3.0
    """Radius of the splatting kernel, 3.0 for Gaussian and 0.5 for triangle"""
    relative_scale: Optional[float] = None
    """Manually set scale when a scene is poorly scaled by nerfstudio
        (e.g. Zip-NeRF dataset, very large-scale scenes across multiple street blocks)"""
    compute_depth_normal: bool = True
    """Compute normal from depth. Required for 2DGS and supervision. Disabling this can reduce VRAM usage and speed up training."""
    packed: bool = False
    """Pack projection outputs, reduce VRAM usage at large batch size but can be slightly slower"""
    use_bvh: bool = False
    """Use BVH for splat-patch intersection test, may be faster when batching large number of small patches"""
    supersampling: int = 1
    """Antialiasing by rendering at higher resolution and downsampling to a lower resolution, as per triangle splatting +"""

    # classial control
    cull_alpha_thresh: float = 0.005
    """threshold of opacity for culling gaussians. One can set it to a lower value (e.g. 0.005) for higher quality."""
    cull_scale_thresh: float = 0.15
    """threshold of world scale for culling huge gaussians"""
    split_scale_thresh: float = 0.05
    """threshold of world scale for splitting huge gaussians"""
    cull_grad_thresh: float = 3e-4  # 3e-4 | 1e-4 | 1e-5 | 0.0
    """threshold for culling gaussians with low visibility"""
    continue_cull_post_densification: bool = True
    """If True, continue to cull gaussians post refinement"""
    reset_alpha_every: int = 30
    """Every this many refinement steps, reset the alpha"""
    densify_xy_grad_thresh: float = 0.005
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
    mcmc_cap_max: int = 1_000_000
    """maximum number of splats for MCMC, dataset-specific tuning required"""
    mcmc_noise_lr: float = 5e5
    """MCMC sampling noise learning rate"""
    mcmc_min_opacity: float = 0.005
    """minimum opacity for MCMC relocation"""
    mcmc_growth_factor: float = 1.05
    """multiply number of splats by this number at every refinement"""
    mcmc_prob_grad_weight: float = 0.0
    """weight of position gradient used in sampling Gaussians to relocate/add to, uses only opacity if 0 and only gradient of 1"""
    relocate_screen_size: float = float('inf')
    """if a gaussian is more than this fraction of screen space, relocate it
        Useful for fisheye with 3DGUT, may drop PSNR for conventional cameras
        For likely better quality, use mcmc_max_screen_size instead"""
    mcmc_max_screen_size: float = 0.15
    """if a gaussian is more than this fraction of screen space, clip scale and increase opacity
        Intended to be an MCMC-friendly alternative of relocate_screen_size"""
    mcmc_max_screen_size_clip_hardness: float = 1.1
    """clip hardness for Gaussians with large screen space size, between 1 and infinity, larger is harder"""
    mcmc_max_world_size: float = float('inf')
    """if a gaussian is more than this of world space, clip scale
        Useful if you see huge floaters at a distance in large indoor space"""

    # representation
    sh_degree: int = 3
    """maximum degree of spherical harmonics to use"""
    sh_degree_interval: int = 1000
    """every n intervals turn on another sh degree; TODO: deprecated?"""
    num_sv: int = 4  # 8
    """number of spherical voronoi to use"""
    train_background_color: bool = True
    """make background color trainable"""
    background_sh_degree: int = 4
    """enable background model"""
    adaptive_exposure_mode: Optional[Literal[
        "linear", "log_linear", "channel", "log_channel", "affine", "log_affine"
        ]] = None
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
    use_bilateral_grid: bool = True
    """If True, use bilateral grid to handle the ISP changes in the image space.
        This technique was introduced in the paper 'Bilateral Guided Radiance Field Processing' (https://bilarfpro.github.io/)."""
    bilagrid_shape: Tuple[int, int, int] = (16, 16, 8)
    """Shape of the bilateral grid (X, Y, W)"""
    use_bilateral_grid_for_geometry: bool = True
    """If True, use bilateral grid for depth and normal (e.g. AI generated biased ones)"""
    bilagrid_shape_geometry: Tuple[int, int, int] = (8, 8, 4)
    """Shape of the bilateral grid for depth and normal (X, Y, W)"""
    use_ppisp: bool = False
    """If True, use the per-pixel ISP model (PPISP) to handle per-pixel color distortions."""
    ppisp_reg_exposure_mean: float = 1.0
    """Encourage exposure mean ~ 0 to resolve SH <-> exposure ambiguity in PPISP."""
    ppisp_reg_vig_center: float = 0.02
    """Encourage vignetting optical center near image center in PPISP."""
    ppisp_reg_vig_non_pos: float = 0.01
    """Penalize positive vignetting alpha coefficients in PPISP (should be <= 0)."""
    ppisp_reg_vig_channel_var: float = 0.1
    """Encourage similar vignetting across RGB channels in PPISP."""
    ppisp_reg_color_mean: float = 1.0
    """Encourage color correction mean ~ 0 across frames in PPISP."""
    ppisp_reg_crf_channel_var: float = 0.1
    """Encourage similar CRF parameters across RGB channels in PPISP."""
    use_linear_color_space: bool = False
    """Whether to assume training images are in linear color space."""
    image_color_space: Literal[None, "ACES2065-1", "ACEScg", "Rec.2020", "AdobeRGB", "DCI-P3"] = None
    """Color space of input images. If set, this only affects displayed images during training, not the final output.
        Note that tonemap is not applied for displayed images."""

    use_3dgs: bool = True
    """Must be True, kept for backward compatibility"""

    # regularization
    scale_regularization_weight: float = 0.0
    """If enabled, a scale regularization introduced in PhysGauss (https://xpandora.github.io/PhysGaussian/) is used for reducing huge spikey gaussians."""
    max_gauss_ratio: float = 10.0
    """threshold of ratio of gaussian max to min scale before applying regularization
    loss from the PhysGaussian paper
    """
    depth_mode: Literal["mean", "median"] = "mean"
    """Depth rendering mode, use mean for stable training and median for high meshing resolution"""
    depth_reg_weight: float = 0.0
    """Weight for depth distortion regularizer"""
    normal_distortion_reg_weight: float = 0.0
    """Weight for normal distortion regularizer"""
    rgb_distortion_reg_weight: float = 0.0
    """Weight for rgb distortion regularizer"""
    distortion_reg_warmup: int = 6000
    """warmup steps for depth regularizer, regularization weight ramps up"""
    normal_reg_weight: float = 0.04
    """Weight for normal regularizer"""
    normal_reg_warmup: int = 6000
    """warmup steps for normal regularizer, regularization weight ramps up"""
    alpha_reg_weight: float = 0.025
    """Weight for alpha regularizer (encourage alpha to go to either 0 or 1)
        Recommend using with --pipeline.model.cull_screen_size for better results"""
    alpha_reg_warmup: int = 12000
    """warmup steps for alpha regularizer, regularization weight ramps up"""
    reg_warmup_length: int = 3000
    """Warmup steps for depth, normal, and alpha regularizers.
       only apply regularizers after this many steps."""
    apply_loss_for_mask: bool = False
    """Set this to False to use masks to ignore distractors (e.g. people and cars, area outside fisheye circle, over exposure)
       Set this to True to remove background (e.g. sky, background outside centered object)"""
    enable_sky_masking: bool = True
    """If enabled, sky from depth map will be used for masking. Alpha loss will always be applied for sky."""
    alpha_loss_weight: float = 0.01
    """Loss weight for alpha, applies when rendered alpha is above reference alpha"""
    alpha_loss_weight_under: float = 0.005
    """Loss weight for alpha, applies when rendered alpha is below reference alpha"""
    mcmc_opacity_reg: float = 0.01  # 0.01 in original paper
    """Opacity regularization from MCMC
       Lower usually gives more accurate geometry"""
    mcmc_scale_reg: float = 0.01  # 0.01 in original paper
    """Scale regularization from MCMC"""
    erank_reg: float = 0.01
    """erank regularization weight, for 3DGS only -
        see *Effective Rank Analysis and Regularization for Enhanced 3D Gaussian Splatting, Hyung et al.*"""
    erank_reg_s3: float = 0.01
    """erank regularization weight for smallest dimension, for 3DGS only"""
    erank_reg_warmup: int = 1000
    """only apply erank regularization after this many steps"""
    exposure_reg_image: float = 0.1
    """Between 0 and 1; For exposure regularization, include this fraction of L1 loss between GT image and image before exposure adjustment"""
    exposure_reg_param: float = 0.002
    """Make sure image look right in auto exposure mode,
       Use an L2 cost to match exposure parameter to the parameter at no exposure adjustment;
       (may be summed/averaged across a batch)"""
    randomize_background: Literal[True, False, "opaque-only", "non-sky-only"] = False
    """Use random noise for background color during training to discourage transparency."""

    # supervision using a foundation depth model
    # enable these by setting `depth_model` in data manager config
    supervision_warmup: int = 1000
    """Start using foundation model depth at this number of steps"""
    depth_distortion_depth_degree: int = -1  # 3
    """Hyperparameter for depth distortion model, controls depth embedding, see code for details
        Larger gives more parameters in depth distortion model, -1 to disable"""
    depth_distortion_uv_degree: int = -1  # 1
    """Hyperparameter for depth distortion model, controls image space embedding, see code for details
        Larger gives more parameters in depth distortion model, -1 to disable"""
    depth_supervision_weight: float = 0.0
    """Weight for depth supervision by comparing rendered depth with depth predicted by a foundation model
        Warn that this can reduce quality if AI generated depth is heavily biased"""
    normal_supervision_weight: float = 0.01
    """Weight for normal supervision by comparing normal from rendered depth with normal from depth predicted by a foundation model"""

    overfit_score_aggregation_mode: Literal['max', 'min', 'mean'] = 'max'
    """Mode to aggregate multiple overfitting objectives.
        Use max for more aggressive early stopping, min for more conservative early stopping, and mean for something in between."""
    validation_loss_average_window: int = 500
    """Window to calculate moving average validation loss for early stop"""
    early_stop_patience: int = 1000
    """Stop training if overfitting score remains positive for this number of iterations"""
    early_stop_warmup: int = 12000
    """Warmup steps for early stop, will not early stop before this number of steps
        Recommend setting this number no less than regularization warmups"""

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
            if len(self.seed_points[0]) > self.config.mcmc_cap_max:
                indices = torch.randperm(len(self.seed_points[0]))[:self.config.mcmc_cap_max]
                self.seed_points = [t[indices] for t in self.seed_points]
            means = self.seed_points[0]
            if self.config.relative_scale is not None:
                means *= self.config.relative_scale
            self.random_init = False
        else:
            means = torch.randn((self.config.num_random, 3)) * self.config.random_scale
            self.random_init = True
        self.xys_grad_norm = None
        self.ch_grad_norm = None
        self.max_2Dsize = None

        scale_init = 0.1 if self.config.use_mcmc else 1.0  # per original papers
        opacity_init = 0.5 if self.config.use_mcmc else 0.1  # per original papers
        if self.config.use_mcmc and self.config.mcmc_max_screen_size < 1.0:
            scale_init, opacity_init = 0.5, 0.1
        if self.config.train_background_color or self.config.randomize_background:
            # scale_init, opacity_init = 1.0, 0.1
            scale_init, opacity_init = 0.5, 0.1

        if self.config.primitive in ["voxel"]:
            means_mean, means_std = torch.mean(means, 0), torch.std(means, 0)
            means_extend = 4.0 * means_std
            pos_min = (means_mean - means_extend).detach().cpu().numpy()
            pos_max = (means_mean + means_extend).detach().cpu().numpy()
            num_voxels = len(means)
            unit_size = 1.0 * np.cbrt(np.prod(pos_max - pos_min) / num_voxels)
            pos_diff = np.ceil((pos_max - pos_min) / unit_size).astype(np.int32)
            self.svhash = svhash_create_initial_volume(pos_min, pos_max, pos_diff)
            self.voxels, self.voxel_indices = svhash_get_voxels(self.svhash)
            num_verts = int(len(self.svhash[0]))
            num_points = int(len(self.voxels))
            density_init = 1.0 * opacity_init / np.cbrt(np.prod(means_std.detach().cpu().numpy()))
            densities = density_init * torch.ones(num_verts)
            # densities = torch.where(densities > 1.1, densities, 1.1 * (torch.log(densities) + 0.904689820196))

        elif self.config.primitive in ["3dgs", "mip", "3dgut"] or True:
            distances, indices = self.k_nearest_sklearn(means.data, 4)
            num_points = means.shape[0]
            scales = torch.sqrt(torch.from_numpy(distances**2).mean(-1, keepdim=True)).repeat(1, 3).float()
            scales = torch.log(scale_init * scales / (self.config.kernel_radius/3.0) + 1e-8)
            quats = F.normalize(torch.randn((num_points, 4)), dim=-1)
            opacities = torch.logit(opacity_init * torch.ones(num_points, 1))

        else:
            distances, indices = self.k_nearest_sklearn(means.data, 6)
            distances = torch.from_numpy(distances)
            # avg_dist = distances.mean(dim=-1, keepdim=True)
            points = means.data[indices] - means.data[:, None, :]
            U, S, Vt = np.linalg.svd(points)
            Vt[:,:,2] *= np.linalg.det(Vt)[:,None]
            num_points = means.shape[0]
            S = np.prod(S, axis=-1, keepdims=True)**(1/3) * np.ones(S.shape)
            scales = S
            scales = np.log(1.5*scales/self.config.kernel_radius+1e-8)
            scales = torch.from_numpy(scales.astype(np.float32))
            quats = torch.from_numpy(np.array(
                Rotation.from_matrix(np.transpose(Vt, axes=(0, 2, 1))).as_quat(),
                dtype=np.float32))
            opacities = torch.logit(0.1 * torch.ones(num_points, 1))

        # if self.config.primitive in ["opaque_triangle"]:
        if self.config.primitive in []:
            M = quat_to_rotmat(quats) * torch.exp(scales[..., None, :])
            means = means.unsqueeze(-2) + torch.bmm(M, torch.tensor([[
                [1, 0, 0], [-0.5, 0.75**0.5, 0], [-0.5, -0.75**0.5, 0]
            ]]).to(M).repeat(len(M), 1, 1))

        gauss_params = {}
        if self.config.primitive in ["3dgs", "mip", "3dgut", "3dgut_sv", "opaque_triangle"]:
            gauss_params['means'] = torch.nn.Parameter(means)
            gauss_params['quats'] = torch.nn.Parameter(quats)
            gauss_params['scales'] = torch.nn.Parameter(scales)
            gauss_params['opacities'] = torch.nn.Parameter(opacities)
        elif self.config.primitive in ["voxel"]:
            gauss_params['densities'] = torch.nn.Parameter(densities)

        # colors
        dim_sh = num_sh_bases(self.config.sh_degree)

        if (
            self.seed_points is not None
            and not self.config.random_init
            # We can have colors without points.
            and self.seed_points[1].shape[0] > 0
        ):
            seed_color = self.seed_points[1] / 255
            if len(seed_color) != num_points:
                seed_color = seed_color.mean(0, True).repeat(num_points, 1)
        else:
            seed_color = torch.rand(num_points, 3)

        if self.config.primitive in ["3dgs", "mip", "3dgut", "opaque_triangle", "voxel"]:
            shs = torch.zeros((num_points, dim_sh, 3)).float()
            shs[:, 0, :3] = (seed_color - 0.5) / 0.28209479177387814
            if self.config.sh_degree > 0:
                shs[:, 1:, 3:] = 0.0
            gauss_params['features_dc'] = torch.nn.Parameter(shs[:, 0, :].contiguous())
            gauss_params['features_sh'] = torch.nn.Parameter(shs[:, 1:, :].contiguous())

        elif self.config.primitive in ["3dgut_sv"]:
            sv_colors = seed_color.unsqueeze(1).repeat(1, self.config.num_sv, 1)
            sv_sites = torch.zeros_like(sv_colors)
            gauss_params['sv_colors'] = torch.nn.Parameter(sv_colors.contiguous())
            gauss_params['sv_sites'] = torch.nn.Parameter(sv_sites.contiguous())

        if self.config.primitive == "opaque_triangle":
            gauss_params['features_ch'] = torch.nn.Parameter(torch.zeros((num_points, 2, 3)))

        self.gauss_params = torch.nn.ParameterDict(gauss_params)

        if self.config.use_camera_optimizer:
            self.camera_optimizer: CameraOptimizer = self.config.camera_optimizer.setup(
                num_cameras=self.num_train_data, device="cpu"
            )

        self.training_losses = SplatTrainingLosses(self.config, self.num_train_data)

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
        if self.config.use_linear_color_space:
            self.background_color = self.background_color ** 2.2
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
        if self.config.primitive in ['voxel']:
            self.strategy = SVRasterStrategy(
            )
            self.strategy_state = self.strategy.initialize_state()
            return

        # opaque triangle mode
        if self.config.primitive == "opaque_triangle":
            current_num = len(self.means)
            final_num = self.config.mcmc_cap_max

            min_warmup_steps = self.config.refine_every * self.config.reset_alpha_every
            num_steps_until_full = math.log(max(final_num, current_num) / current_num) / math.log(self.config.mcmc_growth_factor) * \
                self.config.refine_every + (self.config.warmup_length + self.config.refine_every)
            warmup_steps_0 = min(2.0*max(min_warmup_steps, num_steps_until_full), self.config.num_iterations/3)

            self.strategy = OpaqueStrategy(
                cap_max=self.config.mcmc_cap_max,
                noise_lr=self.config.mcmc_noise_lr,
                refine_start_iter=self.config.mcmc_warmup_length,
                warmup_steps=warmup_steps_0,
                refine_stop_iter=self.config.num_iterations,
                refine_every=self.config.refine_every,
                grow_factor=self.config.mcmc_growth_factor,
                min_opacity=self.config.mcmc_min_opacity,
                relocate_scale2d=self.config.relocate_screen_size,
                max_scale2d=self.config.mcmc_max_screen_size,
                max_scale2d_clip_hardness=self.config.mcmc_max_screen_size_clip_hardness,
                max_scale3d=self.config.mcmc_max_world_size,
                geometry_optimizer_stop_iter=self.config.stop_refine_at
            )
            self.strategy_state = self.strategy.initialize_state()
            return

        # MCMC mode
        if self.config.use_mcmc:
            current_num = len(self.means)
            final_num = self.config.mcmc_cap_max
            warpup_steps = math.log(max(final_num, current_num) / current_num) / math.log(self.config.mcmc_growth_factor) * \
                self.config.refine_every + (self.config.warmup_length + self.config.refine_every)
            self.mcmc_num_steps_until_full = warpup_steps

            self.strategy = MCMCStrategy(
                cap_max=self.config.mcmc_cap_max,
                noise_lr=self.config.mcmc_noise_lr,
                refine_start_iter=self.config.mcmc_warmup_length,
                refine_stop_iter=self.config.stop_refine_at,
                refine_every=self.config.refine_every,
                grow_factor=self.config.mcmc_growth_factor,
                min_opacity=self.config.mcmc_min_opacity,
                prob_grad_weight=self.config.mcmc_prob_grad_weight,
                is_3dgs=self.config.use_3dgs,
                relocate_scale2d=self.config.relocate_screen_size,
                max_scale2d=self.config.mcmc_max_screen_size,
                max_scale2d_clip_hardness=self.config.mcmc_max_screen_size_clip_hardness,
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
            kernel_radius=self.config.kernel_radius,
            absgrad=True,
            revised_opacity=False,
            verbose=True,
        )
        self.strategy_state = self.strategy.initialize_state(scene_scale=1.0)

    @property
    def colors(self):
        if self.config.primitive == "3dgut_sv":
            return self.sv_colors.mean(-2) + 0.5
        C0 = 0.28209479177387814
        return self.features_dc * C0 + 0.5

    @property
    def num_points(self):
        if self.means is not None:
            return self.means.shape[0]
        return self.features_dc.shape[0]

    @property
    def means(self):
        return self.gauss_params.get("means", None)

    @property
    def scales(self):
        return self.gauss_params.get("scales", None)
    
    @property
    def quats(self):
        return self.gauss_params.get("quats", None)

    @property
    def features_dc(self):
        return self.gauss_params.get("features_dc", None)

    @property
    def features_sh(self):
        return self.gauss_params.get("features_sh", None)

    @property
    def features_ch(self):
        return self.gauss_params.get("features_ch", None)

    @property
    def sv_sites(self):
        return self.gauss_params.get("sv_sites", None)

    @property
    def sv_colors(self):
        return self.gauss_params.get("sv_colors", None)

    @property
    def opacities(self):
        return self.gauss_params.get("opacities", None)

    @property
    def densities(self):
        return self.gauss_params.get("densities", None)

    def load_state_dict(self, dict, **kwargs):  # type: ignore
        # resize the parameters to match the new number of points
        self.step = self.config.num_iterations
        for p in ["means", "scales", "quats",
                    "features_dc", "features_sh", "features_ch",
                    "sv_sites", "sv_colors",
                    "opacities", "densities"]:
            if p in dict:
                dict[f"gauss_params.{p}"] = dict[p]
        newp = dict["gauss_params.features_dc"].shape[0]
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
        if self.config.primitive == 'voxel':
            self.svhash, self.voxels, self.voxel_indices = \
            self.strategy.step_post_backward(
                self.svhash, self.voxels, self.voxel_indices,
                params=self.gauss_params,
                optimizers=self.optimizers,
                state=self.strategy_state,
                step=self.step,
                info=self.info,
                packed=(self.config.packed or self.config.use_bvh),
            )
        elif self.config.use_mcmc:
            self.strategy.step_post_backward(
                params=self.gauss_params,
                optimizers=self.optimizers,
                state=self.strategy_state,
                step=self.step,
                info=self.info,
                lr=self.optimizers['means'].param_groups[0]['lr'],
                packed=(self.config.packed or self.config.use_bvh),
            )
        else:
            self.strategy.step_post_backward(
                params=self.gauss_params,
                optimizers=self.optimizers,
                state=self.strategy_state,
                step=self.step,
                info=self.info,
                packed=(self.config.packed or self.config.use_bvh),
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
                "sv_sites", "sv_colors",
                "opacities", "densities"
            ] if name in self.gauss_params
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
        if self.config.use_bilateral_grid_for_geometry:
            gps["bilateral_grid_geometry"] = \
                list(self.training_losses.bil_grids_depth.parameters()) + \
                list(self.training_losses.bil_grids_normal.parameters())
        if self.config.use_ppisp:
            gps["ppisp"] = [self.training_losses.ppisp_params]

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

    def get_background_image(self, camera: Cameras, c2w: torch.Tensor, intrins: torch.Tensor):
        if not isinstance(camera, Cameras):
            print("Called get_background_image with not a camera")
            return {}

        W, H = int(camera.width[0].item()), int(camera.height[0].item())

        if self.config.randomize_background == True:
            # return torch.rand_like(self.background_color).repeat(H, W, 1)
            return torch.rand((len(camera), H, W, 3), device=self.background_color.device)

        sh_degree = self.config.background_sh_degree
        if not self.config.train_background_color or not (sh_degree > 0):
            return self.background_color[None].repeat(len(camera), H, W, 1)

        sh_coeffs = torch.cat((self.background_color.unsqueeze(0), self.background_sh), dim=0)  # [(deg+1)^2, 3]
        return render_background_sh(
            camera.width[0].item(), camera.height[0].item(),
            ['pinhole', 'fisheye'][camera.camera_type[0].item() == CameraType.FISHEYE.value],
            intrins, c2w[..., :3, :3], sh_degree, sh_coeffs
        )

    @staticmethod
    def get_empty_outputs(width: int, height: int, background: torch.Tensor) -> Dict[str, Union[torch.Tensor, List]]:
        return {}

    def get_outputs(self, camera: Cameras, val: bool=False) -> Dict[str, Union[torch.Tensor, List]]:
        """Takes in a camera and returns a dictionary of outputs."""

        if isinstance(camera, Tuple) and len(camera) == 2:
            train_outputs = self.get_outputs(camera[0])
            # with torch.no_grad():
            val_outputs = self.get_outputs(camera[1], val=True)
            return train_outputs, val_outputs

        if not isinstance(camera, Cameras):
            print("Called get_outputs with not a camera")
            return {}

        device = self.means.device if self.means is not None else self.features_dc.device

        if self.training and self.config.use_camera_optimizer:
            camera.metadata['cam_idx'] = camera.metadata['cam_idx'].flatten()
            optimized_camera_to_world = self.camera_optimizer.apply_to_camera(camera)
        else:
            optimized_camera_to_world = camera.camera_to_worlds
        optimized_camera_to_world = optimized_camera_to_world.cpu()

        camera_downscale = self._get_downscale_factor()
        if camera_downscale > 1:
            camera.rescale_output_resolution(1 / camera_downscale)

        # TODO: separate different sizes/intrins
        W, H = int(camera.width[0].item()), int(camera.height[0].item())

        R = optimized_camera_to_world[:, :3, :3]  # 3 x 3
        T = optimized_camera_to_world[:, :3, 3:4]  # 3 x 1
        if self.config.relative_scale is not None:
            T = T * self.config.relative_scale
        R = R * torch.tensor([[[1.0, -1.0, -1.0]]])
        R_inv = R.transpose(-1, -2)
        T_inv = -torch.bmm(R_inv, T)
        viewmats = torch.eye(4, dtype=optimized_camera_to_world.dtype)[None].repeat(len(camera), 1, 1)
        viewmats[:, :3, :3] = R_inv
        viewmats[:, :3, 3:4] = T_inv

        max_depth_scale = 2.0 if self.training else 1.0

        assert self.config.use_3dgs, "2DGS is deprecated"

        # Call GSplat for 3DGS rendering
        intrins = torch.concatenate((camera.fx, camera.fy, camera.cx, camera.cy), dim=1)
        
        kwargs = {'actual_width': W, 'actual_height': H}
        if 'actual_width' in camera.metadata:
            kwargs['actual_width'] = int(camera.metadata['actual_width'] + 0.5)
        if 'actual_height' in camera.metadata:
            kwargs['actual_height'] = int(camera.metadata['actual_height'] + 0.5)
        if 'patch_offsets' in camera.metadata:
            kwargs['patch_offsets'] = camera.metadata['patch_offsets']
        is_fisheye = (camera.camera_type[0].item() == CameraType.FISHEYE.value)
        if not self.training:
            pass
            # is_fisheye = True
            # kwargs['dist_coeffs'] = torch.tensor([[0.0259, 0.0082, 0.0002, -0.0013, -0.0012, -0.0008, 0.0000, 0.0000, -0.0006, -0.0001]]).float().cuda()

            # TODO: investigate why this uses a ton of VRAM
            # kwargs['dist_coeffs'] = torch.tensor([[-0.29, 0.07, 0, 0, 1e-5, -1e-3, 0, 0, 0, 0]]).float().cuda()

        if camera.distortion_params is not None and camera.distortion_params.any():
            kwargs['dist_coeffs'] = camera.distortion_params

        if 'visibility_masks' in camera.metadata:
            visibility_masks = camera.metadata['visibility_masks']
            visibility_masks = self._downscale_if_required(visibility_masks)
            kwargs['masks'] = visibility_masks.bool()

        optimized_camera_to_world = optimized_camera_to_world.to(device)
        viewmats = viewmats.to(device)
        intrins = intrins.to(device)

        TILE_SIZE = 64
        gh, gw = (H+TILE_SIZE-1) // TILE_SIZE, (W+TILE_SIZE-1) // TILE_SIZE
        def split_into_tiles(viewmat, intrins):
            dh, dw = torch.meshgrid(torch.arange(gh)*TILE_SIZE, torch.arange(gw)*TILE_SIZE)
            intrins = intrins.clone().repeat(gh*gw, 1)
            viewmat = viewmat.clone().repeat(gh*gw, 1, 1)
            intrins[:, 2] -= dw.flatten().to(intrins)
            intrins[:, 3] -= dh.flatten().to(intrins)
            return viewmat, intrins

        def merge_tiles(im):
            im = im.reshape(gh, gw, TILE_SIZE, TILE_SIZE, -1).permute(0, 2, 1, 3, 4).reshape(1, gh*TILE_SIZE, gw*TILE_SIZE, -1)
            im = im[:, :camera.height[0].item(), :camera.width[0].item(), :]
            return im

        # if not self.training:
        #     viewmats_0, intrins_0 = viewmats, intrins
        #     viewmats, intrins = split_into_tiles(viewmats, intrins)
        #     if 'dist_coeffs' in kwargs:
        #         kwargs['dist_coeffs'] = kwargs['dist_coeffs'].repeat(gh*gw, 1, 1)
        #     W = TILE_SIZE
        #     H = TILE_SIZE

        if self.config.primitive in ['3dgs', 'mip', '3dgut']:
            splat_params = (
                self.means, F.normalize(self.quats, dim=-1), self.scales,
                self.opacities.squeeze(-1),
                self.features_dc, self.features_sh
            )
        if self.config.primitive in ['3dgut_sv']:
            splat_params = (
                self.means, F.normalize(self.quats, dim=-1), self.scales,
                self.opacities.squeeze(-1),
                self.sv_sites, self.sv_colors
            )
        elif self.config.primitive in ['opaque_triangle']:
            # hardness = min(max(self.step / self.config.stop_refine_at, 0.1), 1.0)
            opacity_floor = self.strategy.get_opacity_floor(self.step)
            hardness = self.strategy.get_hardness(self.step)
            splat_params = (
                self.means, F.normalize(self.quats, dim=-1), self.scales,
                #  hardness * torch.ones_like(self.opacities.squeeze(-1))
                #  hardness + (1.0-hardness) * torch.sigmoid(self.opacities.squeeze(-1))
                torch.concat([
                    self.strategy.map_opacities(self.step, self.opacities),
                    hardness * torch.ones_like(self.opacities)
                ], dim=-1),
                self.features_dc, self.features_sh, self.features_ch
            )
        elif self.config.primitive in ['voxel']:
            # print(self.svhash)
            # torch.cuda.synchronize()
            # print(torch.amin(voxel_indices), torch.amax(voxel_indices), self.densities.shape)
            # print(voxels)
            # print(voxel_indices)
            # exit(0)
            splat_params = (
                self.voxels, self.densities[self.voxel_indices],
                self.features_dc, self.features_sh
            )
            # print([x.shape for x in splat_params])
        if val:
            splat_params = tuple([(p.detach() if isinstance(p, torch.Tensor) else p) for p in splat_params])

        use_bvh = self.config.use_bvh and self.training and not val
        rgbd, alpha, meta = rasterization(
            self.config.primitive,
            splat_params,
            # "voxel",
            # (torch.cat((self.means, 20.0*torch.exp(self.scales.mean(-1, True))), dim=-1), torch.exp(self.opacities).repeat(1, 8), self.features_dc, self.features_sh),
            # # (torch.cat((self.means, 0.025*torch.ones_like(self.scales.mean(-1, True))), dim=-1), torch.exp(self.opacities).repeat(1, 8), self.features_dc, self.features_sh),

            # (self.means, hardness * torch.ones_like(self.opacities.squeeze(-1))),
            viewmats=viewmats,  # [C, 4, 4]
            intrins=intrins * self.config.supersampling,  # [C, 4]
            width=W * self.config.supersampling,
            height=H * self.config.supersampling,
            packed=(self.config.packed or use_bvh),
            use_bvh=(use_bvh),
            # packed=True,
            # use_bvh=True,
            relative_scale=self.config.relative_scale,
            sparse_grad=False,
            distributed=False,
            camera_model=["pinhole", "fisheye"][is_fisheye],
            render_mode="RGB+D" if self.config.primitive in ['3dgs', 'mip', '3dgut', '3dgut_sv'] else "RGB+D+N",
            output_distortion=any([c != 0.0 for c in self.training_losses.get_2dgs_reg_weights()[0]]),
            **kwargs,
        )
        if self.config.supersampling != 1:
            rgbd = [resize_image(im, self.config.supersampling) for im in rgbd]
            alpha = resize_image(alpha, self.config.supersampling)

        # if not self.training:
        #     W = gw*TILE_SIZE
        #     H = gh*TILE_SIZE
        #     # print(rgbd.shape, alpha.shape)
        #     rgbd = [merge_tiles(comp) for comp in rgbd]
        #     alpha = merge_tiles(alpha)
        #     for key in ['rgb_distortion', 'depth_distortion', 'normal_distortion']:
        #         if key in meta:
        #             meta[key] = merge_tiles(meta[key])
        #     W, H = camera.width[0].item(), camera.height[0].item()
        #     viewmats, intrins = viewmats_0, intrins_0
        #     if 'dist_coeffs' in kwargs:
        #         kwargs['dist_coeffs'] = kwargs['dist_coeffs'][:1]

        if self.config.compute_depth_normal or not self.training:
            depth_im_ref = torch.where(
                alpha > 0.0, rgbd[1],
                max_depth_scale*torch.amax(rgbd[1]).detach()
            ).contiguous()
        else:
            depth_im_ref = None

        # normals
        render_normal = None
        if len(rgbd) > 2:
            render_normal = torch.where(alpha > 0.0, F.normalize(rgbd[2], dim=-1), rgbd[2])

        depth_normal = None
        if self.config.fit == "depth_normal" or not self.training:
            depth_normal = depth_to_normal(
                depth_im_ref, ["pinhole", "fisheye"][is_fisheye], intrins, **kwargs
            )

        if self.config.fit == "rgb":
            rgb = rgbd[0]
        elif self.config.fit == "depth":
            rgb = depth_im_ref
            rgb = rgb / rgb.mean()
            rgb = (rgb / (1.0 + rgb)).repeat(1, 1, 1, 3)
        elif self.config.fit == "normal":
            assert depth_normal is not None, "Depth normal map not available"
            factor = self.config.fit_normal_depth_factor * min(
                -1 + 2 * self.step / self.config.fit_normal_depth_warmup_steps, 1.0)
            if render_normal is None:
                factor = 1.0
            if np.random.random() < factor:# and self.training:
                rgb = 0.5+0.5*depth_normal
            else:
                rgb = 0.5+0.5*render_normal

        means2d = meta["means2d"]
        if self.config.use_mcmc:
            means2d = self.means
        radii = meta["radii"]
        depths = meta["depths"]

        if self.training:
            self.info = {
                "width": kwargs["actual_width"],
                "height": kwargs["actual_height"],
                "n_cameras": camera.shape[0],
                "n_train": self.num_train_data,
                "radii": radii,
                "means2d": means2d,
                "depths": depths,
            }
            if 'patch_offsets' in camera.metadata:
                self.info['patch_offsets'] = camera.metadata['patch_offsets']
            if 'actual_images_per_batch' in camera.metadata:
                self.info['n_cameras'] = camera.metadata['actual_images_per_batch']
            for key in ['gaussian_ids', 'camera_ids', 'max_blending', 'bvh_time']:
                if key in meta:
                    self.info[key] = meta[key]

        # blend with background
        if self.config.fit == "rgb":
            background = self.get_background_image(camera, optimized_camera_to_world, intrins)
            # rgb = torch.clip(rgb + (1.0 - alpha) * background, 0.0, 1.0)
            rgb = blend_background(rgb, alpha, background)
        else:
            background = torch.zeros_like(rgb)

        # visualize PPISP for debugging
        if not self.training and False:
            ppisp_param = self.training_losses.ppisp_params[0:1]
            # ppisp_param = ppisp_param.detach() + math.sin(self.step / 10.0)
            from spirulae_splat.modules.training_losses import apply_ppisp
            rgb = apply_ppisp(rgb, ppisp_param, intrins)

        # pack outputs
        outputs = {
            "rgb": rgb,
        }
        if depth_im_ref is not None:
            outputs["depth"] = depth_im_ref
        if render_normal is not None:
            outputs["normal"] = render_normal
        if depth_normal is not None:
            outputs["depth_normal"] = depth_normal
        outputs["alpha"] = alpha
        outputs["background"] = background

        for key in ['rgb_distortion', 'depth_distortion', 'normal_distortion']:
            if key in meta:
                value = meta[key]
                if value is not None:
                    if self.config.supersampling != 1:
                        value = resize_image(value, self.config.supersampling)
                    if not self.training:
                        value = torch.sqrt(value + (1/255)**2) - (1/255)
                    outputs[key] = value

        if not self.training:
            outputs["alpha"] = outputs["alpha"].reshape((H, W, 1)).repeat(1, 1, 3)
            outputs["background"] = torch.clip(outputs["background"], min=0.0, max=1.0)

        # convert linear depth to ray depth, for correct gl_z_buf_depth in Viser
        if not self.training:
            # undist_map = ssplat_camera.get_undist_map(always=True)
            # distances = torch.sqrt((undist_map*undist_map).sum(-1, True) + 1.0)
            # outputs["depth"] = outputs["depth"] * distances
            # outputs["depth"] = torch.clip(outputs["depth"], max=torch.quantile(outputs["depth"], 0.99))
            if self.config.relative_scale is not None:
                outputs["depth"] /= self.config.relative_scale
            if "normal" in outputs:
                outputs["normal"] = 0.5+0.5*outputs["normal"]
            if "depth_normal" in outputs:
                outputs["depth_normal"] = 0.5+0.5*outputs["depth_normal"]
            if self.config.use_linear_color_space or self.config.image_color_space != None:
                outputs["rgb_raw"] = outputs["rgb"]
                for key in ['rgb', 'background']:
                    if key not in outputs:
                        continue
                    if self.config.image_color_space != None:
                        color_transform = get_color_transform_matrix(self.config.image_color_space)
                        outputs[key] = torch.matmul(outputs[key], color_transform.T).clip(0, 1)
                    if self.config.use_linear_color_space:
                        outputs[key] = linear_rgb_to_srgb(outputs[key])
            for key in outputs:
                outputs[key] = outputs[key].squeeze(0)

        if self.training:
            kwargs["intrins"] = intrins
            kwargs["camera_model"] = ["pinhole", "fisheye"][is_fisheye]
            outputs["camera"] = camera
            outputs["camera_intrins"] = kwargs
        # if not self.training and True:
        #     outputs['ray'] = merge_tiles(meta['intersection_count'].float().reshape(-1, 1, 1, 1).repeat(1, TILE_SIZE, TILE_SIZE, 1))

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

        # Per-splat losses (CUDA async)
        static_losses = {}
        self.training_losses.get_static_losses(
            self.step,
            self.quats, self.scales, self.opacities,
            static_losses
        )

        # Camera optimizer loss (depends on nerfstudio, usually CPU)
        if self.training and self.config.use_camera_optimizer:
            self.camera_optimizer.get_loss_dict(static_losses)

        # Per-image losses (CUDA with .item() that involves GPU-CPU sync)
        self.info['num_train_data'] = self.num_train_data

        # Outputs and batch are tuples if eval data is provided
        # self.training_losses will call CUDA functions to compute training loss
        total_val_loss = None
        if isinstance(outputs, tuple) and isinstance(batch, tuple):
            loss_dict = self.training_losses(self.step, batch[0], outputs[0], self.info)
            # for key, value in outputs[1].items():
            #     if isinstance(value, torch.Tensor):
            #         outputs[1][key] = value.detach()
            val_loss_dict = self.training_losses(self.step, batch[1], outputs[1], self.info, True)
            total_val_loss = torch.stack([
                x for x in val_loss_dict.values() if isinstance(x, torch.Tensor)
            ]).sum()
        else:
            loss_dict = self.training_losses(self.step, batch, outputs, self.info)

        # Total train and validation losses
        with torch.no_grad():
            total_train_loss = torch.stack([
                x for x in loss_dict.values() if isinstance(x, torch.Tensor)
            ]).sum()
            # prevent double backward
            if isinstance(total_train_loss, torch.Tensor):
                total_train_loss = total_train_loss.item()

        if total_val_loss is not None:
            # will have gradient to e.g. bilagrid, camera poses, but not Gaussian params
            loss_dict['val'] = total_val_loss
            # prevent double backward
            if isinstance(total_val_loss, torch.Tensor):
                total_val_loss = total_val_loss.item()

        # Add static losses
        loss_dict.update(static_losses)

        # Store train/val loss running stats
        if not hasattr(self, 'train_loss_history'):
            self.train_loss_history = []
        self.train_loss_history.append(total_train_loss)
        if not hasattr(self, 'val_loss_history'):
            self.val_loss_history = []
        if not hasattr(self, 'val_lpips_history'):
            self.val_lpips_history = []
        if total_val_loss is not None:
            while len(self.val_loss_history) < len(self.train_loss_history):
                self.val_loss_history.append(total_val_loss)
            while len(self.val_lpips_history) < len(self.train_loss_history):
                self.val_lpips_history.append(val_loss_dict['lpips_val'])
            # TODO: possbily linear interpolation?
        sw_width = self.config.validation_loss_average_window
        if len(self.val_loss_history) >= sw_width:
            # TODO: O(1) moving average update
            self.total_train_loss = sum(self.train_loss_history[-sw_width:]) / sw_width
            self.total_val_loss = sum(self.val_loss_history[-sw_width:]) / sw_width
            self.total_val_lpips = sum(self.val_lpips_history[-sw_width:]) / sw_width

        # Compute overfitting score
        if not hasattr(self, 'overfit_count'):
            self.overfit_count = 0
        if len(self.val_loss_history) >= 2*sw_width:
            # note that training loss can increase at e.g. regularization weight scheduling
            prev_total_train_loss = sum(self.train_loss_history[-2*sw_width:-sw_width]) / sw_width
            prev_total_val_loss = sum(self.val_loss_history[-2*sw_width:-sw_width]) / sw_width
            prev_total_val_lpips = sum(self.val_lpips_history[-2*sw_width:-sw_width]) / sw_width
            train_loss_increase = self.total_train_loss - prev_total_train_loss
            val_loss_increase = self.total_val_loss - prev_total_val_loss
            val_lpips_increase = self.total_val_lpips - prev_total_val_lpips
            overfit_scores = []
            overfit_scores.append(val_loss_increase - max(train_loss_increase, 0))
            overfit_scores.append(val_lpips_increase)
            overfit_score = max(overfit_scores) if self.config.overfit_score_aggregation_mode == 'max' else \
                            min(overfit_scores) if self.config.overfit_score_aggregation_mode == 'min' else \
                            sum(overfit_scores) / len(overfit_scores)
            # overfit_score = min(val_loss_increase - max(train_loss_increase, 0), val_lpips_increase)
            if not hasattr(self, 'overfit_score_history'):
                self.overfit_score_history = []
            self.overfit_score_history.append(overfit_score)
            # Early stop of overfits
            if len(self.overfit_score_history) >= sw_width:
                self.overfit_score = sum(self.overfit_score_history[-sw_width:]) / sw_width
                if self.overfit_score > 0.0:
                    self.overfit_count += 1
                    if self.overfit_count > self.config.early_stop_patience and \
                            self.step > self.config.early_stop_warmup + self.config.early_stop_patience//2:
                        print("\n\n\n\n\n\nOverfitting detected. Quality is becoming worse with more training. Quitting.\n")
                        exit(0)
                else:
                    self.overfit_count = 0

        self.print_loss_dict(loss_dict)

        # Reduce before return -
        # nerfstudio will later use functools.reduce(torch.add, ...), which is slow
        loss_dict = {'total': torch.stack([
            x for x in loss_dict.values() if isinstance(x, torch.Tensor)
        ]).sum()}
        return loss_dict

    def print_loss_dict(self, losses: Dict[str, torch.Tensor], _storage={}):

        # print can take a few ms depending on system, do in a separate thread
        if 'print_queue' not in _storage and False:
            import threading
            import sys
            import time

            print_queue = []

            def printer():
                while True:
                    if len(print_queue) == 0:
                        continue
                    msg = print_queue[0]
                    print_queue.clear()
                    sys.stdout.write(msg)
                    sys.stdout.flush()
                    time.sleep(0.1)

            threading.Thread(target=printer, daemon=True).start()

            _storage['print_queue'] = print_queue

        # caching
        skip = 10
        if self.step % skip != 0 and 'chunks' in _storage:
            chunks = _storage['chunks']
            # CONSOLE.print(chunks, end='')  # slow
            print(chunks, end='')
            # _storage['print_queue'].append(chunks)
            return

        # text formatting (instead of using rich.CONSOLE that's slow)
        def bracket(s: str):
            return f"\033[1m[\033[m{s}\033[1m]\033[m"

        def orange(s: str):
            return f"\033[33m{s}\033[m"

        def boldcyan(s: str):
            return f"\033[1;36m{s}\033[m"

        def redbkg(s: str, threshold: float = 1.0):
            return f"\033[1;41m{s.replace('\033[m', '')}\033[m" if threshold >= 1.0 else s

        # get VRAM usage (only supports single GPU at this time)
        device = self.means.device if self.means is not None else self.features_dc.device
        free, total = torch.cuda.mem_get_info(device)
        used = (total - free) / 1024**3
        used_percentage = (1 - free/total)*100
        mem_stats = boldcyan(f"{used:.2f}") + f"\N{ZERO WIDTH SPACE}GB " + boldcyan(f"{used_percentage:.0f}") + "%"

        if hasattr(self, 'total_val_loss'):
            losses['val_total'] = self.total_val_loss
        if hasattr(self, 'total_train_loss'):
            losses['train_total'] = self.total_train_loss
        if hasattr(self, 'total_val_lpips'):
            losses['lpips_val'] = self.total_val_lpips
        if hasattr(self, 'overfit_score'):
            losses['overfit_score'] = self.overfit_score

        def fmt(key: str, s: float, decimals=None) -> str:
            if s == 0.0 or key not in losses:
                return '~'

            l = losses[key]
            if isinstance(l, torch.Tensor):
                l = l.detach().item()
            l = l / s
            if not math.isfinite(l):  # not finite
                return str(l)

            if key not in _storage or self.step % 1000 == 0:
                _storage[key] = abs(l)
            _storage[key] = max(_storage[key], abs(l))
            if _storage[key] == 0.0:
                return '~'

            if decimals is None:  # 3 sig figs
                decimals = int(max(-math.log10(0.001*abs(l)), 0)) if l != 0.0 else 0
            s = f"{{:.{decimals}f}}".format(l)
            if s.startswith('0.'):
                s = s[1:]
            return boldcyan(s)

        reg_mcmc = (self.config.use_mcmc and self.step < self.config.stop_refine_at)
        reg_2dgs = self.training_losses.get_2dgs_reg_weights()
        opacity_floor = (
            f"{bracket('OpacFloor')} {self.strategy.get_opacity_floor(self.step):.3f} "
            f"{bracket('Hardness')} {self.strategy.get_hardness(self.step):.3f}"
        ).replace('0.', '.') \
            if self.config.primitive == "opaque_triangle" else ""
        chunks = [
            f"{bracket('N')} {boldcyan(self.num_points)}",
            f"{redbkg(bracket('Mem'), used_percentage/90)} {mem_stats}",
            f"{bracket('Train')} {orange('loss')}={fmt('image_loss', 1.0)} "
            f"{orange('psnr')}={fmt('psnr', 1.0, 2)} "
            f"{orange('ssim')}={fmt('ssim', 1.0, 3)}" + \
                f" lpips={fmt('lpips', 1.0, 3)}" * ('lpips' in losses),
            "                \n",
            f"{bracket('RefLoss')} {orange('depth')}={fmt('depth_ref_loss', self.config.depth_supervision_weight, 3)} "
            f"{orange('normal')}={fmt('normal_ref_loss', self.config.normal_supervision_weight, 3)} "
            f"{orange('alpha')}={fmt('alpha_ref_loss', 0.5*(self.config.alpha_loss_weight+self.config.alpha_loss_weight_under), 4)}",
            f"{bracket('DistLoss')} {orange('depth')}={fmt('depth_dist_reg', reg_2dgs[0][0], 3)} "
            f"{orange('normal')}={fmt('normal_dist_reg', reg_2dgs[0][1], 3)} "
            f"{orange('rgb')}={fmt('rgb_dist_reg', reg_2dgs[0][2], 3)}",
            "                \n",
            f"{bracket('ImReg')} {orange('normal')}={fmt('normal_reg', reg_2dgs[1], 3)} "
            f"{orange('alpha')}={fmt('alpha_reg', self.training_losses.get_alpha_reg_weight(), 3)}",
            f"{bracket('SplatReg')} {orange('opac')}={fmt('mcmc_opacity_reg', self.config.mcmc_opacity_reg * reg_mcmc, 3)} "
            f"{orange('scale')}={fmt('mcmc_scale_reg', self.config.mcmc_scale_reg * reg_mcmc, 4)} "
            f"{orange('erank')}={fmt('erank_reg', max(self.config.erank_reg_s3, self.config.erank_reg, 3))} "
            f"{orange('aniso')}={fmt('scale_reg', self.config.scale_regularization_weight, 3)}",
            "                \n",
            f"{bracket('BilagridTVLoss')} {orange('rgb')}={fmt('tv_loss', 10.0)} "
            f"{orange('depth')}={fmt('tv_loss_depth', 10.0)} "
            f"{orange('normal')}={fmt('tv_loss_normal', 10.0)}",
            f"{bracket('PPISPReg')} {orange('eμ')}={fmt('ppisp_reg_exposure_mean', self.config.ppisp_reg_exposure_mean)}"
            "                \n"
            f"{orange('vc')}={fmt('ppisp_reg_vig_center', self.config.ppisp_reg_vig_center)} "
            f"{orange('v+')}={fmt('ppisp_reg_vig_non_pos', self.config.ppisp_reg_vig_non_pos)} "
            f"{orange('vσ')}={fmt('ppisp_reg_vig_channel_var', self.config.ppisp_reg_vig_channel_var)} "
            f"{orange('cμ')}={fmt('ppisp_reg_color_mean', self.config.ppisp_reg_color_mean)}",
            f"{orange('rσ')}={fmt('ppisp_reg_crf_channel_var', self.config.ppisp_reg_crf_channel_var)}",
            "                \n",
            f"{bracket('Validation')} {orange('train')}={fmt('train_total', 1.0)} "
            f"{orange('val')}={fmt('val_total', 1.0)} "
            f"{orange('lpips')}={fmt('lpips_val', 1.0)} "
            f"{orange('overfit')}={fmt('overfit_score', 1.0)} "
            f"{redbkg(f'{self.overfit_count}/{self.config.early_stop_patience}',
                self.overfit_count/(0.8*self.config.early_stop_patience))}"
        ]
        additional_chunks = []
        if len(opacity_floor) > 0:
            additional_chunks.append(opacity_floor)
        if 'bvh_time' in self.info:
            losses['bvh_time'] = self.info['bvh_time']
            additional_chunks.append(f"{bracket('BvhTime')} {fmt('bvh_time', 1.0, 1)} ms")
        if len(additional_chunks) > 0:
            chunks.append("                \n")
            chunks.extend(additional_chunks)
        chunks.append("                \n")
        chunks = ' '.join(chunks).replace('\n ', '\n')
        chunks += "\033[F"*(chunks.count('\n'))
        print(chunks, end='')
        # _storage['print_queue'].append(chunks)
        _storage['chunks'] = chunks

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
        # return {}, {}   # TODO
        # gt_rgb = self.composite_with_background(self.get_gt_img(batch["image"]), outputs["background"])
        gt_rgb = batch["image"]  # TODO
        if gt_rgb.dtype == torch.uint8:
            gt_rgb = gt_rgb.float() / 255.0
        predicted_rgb = outputs["rgb"]

        combined_rgb = torch.cat([gt_rgb, predicted_rgb], dim=1)

        # Switch images from [H, W, C] to [1, C, H, W] for metrics computations
        gt_rgb = torch.moveaxis(gt_rgb, -1, 0)[None, ...]
        predicted_rgb = torch.moveaxis(predicted_rgb, -1, 0)[None, ...]

        # metrics
        from torchmetrics.image import PeakSignalNoiseRatio
        from torchmetrics.image.lpip import \
            LearnedPerceptualImagePatchSimilarity

        self.psnr = PeakSignalNoiseRatio(data_range=1.0).to(self.device)
        self.ssim = SSIM(data_range=1.0, size_average=True, channel=3).to(self.device)
        self.lpips = LearnedPerceptualImagePatchSimilarity(normalize=True).to(self.device)

        psnr = self.psnr(gt_rgb, predicted_rgb)
        ssim = self.ssim(gt_rgb, predicted_rgb)
        lpips = self.lpips(gt_rgb, predicted_rgb)

        # all of these metrics will be logged as scalars
        metrics_dict = {"psnr": float(psnr.item()), "ssim": float(ssim)}  # type: ignore
        metrics_dict["lpips"] = float(lpips)
        metrics_dict["gaussian_count"] = float(self.num_points)

        images_dict = {"img": combined_rgb}

        return metrics_dict, images_dict
