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

import spirulae_splat
from spirulae_splat.splat._torch_impl import quat_to_rotmat
from spirulae_splat.modules.core import Renderer
from spirulae_splat.splat.utils import resize_image
from spirulae_splat.splat.background_sh import render_background_sh
from spirulae_splat.splat.sh import num_sh_bases, spherical_harmonics
from spirulae_splat.strategy import MCMCStrategy, OpaqueStrategy, SVRasterStrategy

from spirulae_splat.modules.training_losses import SplatTrainingLosses
from spirulae_splat.modules.optimizer import get_scheduled_lr
from spirulae_splat.splat.cuda._wrapper_per_pixel import (
    blend_background,
    blend_background_noise,
    rgb_to_srgb,
    get_color_transform_matrix,
    _make_lazy_cuda_func
)
from spirulae_splat.splat.cuda._wrapper_projection import scatter_max
from spirulae_splat.splat.cuda import (
    svhash_create_initial_volume,
    svhash_get_voxels
)

from spirulae_splat.modules.camera import Cameras, CameraType


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
class SpirulaeSplatModelConfig:

    # Representation
    primitive: Literal["3dgs", "mip", "3dgut", "3dgut_sv", "opaque_triangle", "voxel"] = "3dgut"
    """Splat primitive to use"""
    sh_degree: int = 3
    """Maximum degree of spherical harmonics to use."""
    num_sv: int = 4  # 8
    """Number of spherical voronoi to use."""
    background_color: Literal["black", "white", "gray"] = "black"
    """Whether to randomize the background color."""
    train_background_color: bool = False
    """Whether to train background color."""
    background_sh_degree: int = 4
    """SH degree for background color."""
    randomize_background: Literal[True, False, "opaque-only", "non-sky-only"] = False
    """Use random noise for background color during training to discourage transparency."""
    randomize_background_warmup: int = 2000
    """Number of steps to warmup background randomization. This applies when randomize_background is True"""
    randomize_background_pre_warmup: float = 0.25
    """Weight of randomization at start of training (0 to 1). Higher value reduce the chance of washing away splat opacities near the beginning of training."""
    kernel_radius: float = 3.0
    """Radius of the splatting kernel, 3.0 for Gaussian and 0.5 for triangle"""

    # Training loss
    relative_scale: Optional[float] = None
    """Manually set scale when a scene is poorly scaled, i.e. increase this for large datasets.
        If not set, will use a scale agnostic optimizer. To prevent this, set it to 1.0."""
    l2_lambda: float = 0.0
    """Weight of L2 loss, default 0.0"""
    ssim_lambda: float = 0.2
    """Weight of ssim loss; 0.2 for academic baseline, higher for potentially more high-frequency details, lower for less blurry background in outdoor scenes"""
    lpips_lambda: float = 0.0
    """Weight of lpips loss for better perceptual quality; Note that this can make training much slower"""
    num_loss_scales: int = 0
    """Number of scales for image loss. For multi-scale loss, image is downscaled by 2 this number of times, and losses are averaged across scales. Improves convergence for high-resolution images."""
    use_camera_optimizer: bool = False
    """Whether to use camera optimizer
        Note: this only works well in patch batching mode"""
    # camera_optimizer: CameraOptimizerConfig = field(default_factory=lambda: CameraOptimizerConfig(mode="SO3xR3"))
    # """Config of the camera optimizer to use"""  # TODO
    packed: bool = True
    """Pack projection outputs, reduce VRAM usage at large batch size but can be slightly slower"""
    use_bvh: bool = False
    """Use BVH for splat-patch intersection test, may be faster when batching large number of small patches"""
    compute_hessian_diagonal: Literal[None, "position", "all"] = None
    """What parameter sets to compute an approximation of Hessian diagonal as well as a Jacobian-residual product in backward pass. Required for second-order optimizer."""
    optimizer_offload: Literal[None, "sh", "all"] = None
    """Whether to offload optimizer momentum to CPU to save VRAM. This is only supported for Adam optimizer."""
    supersampling: int = 1
    """Antialiasing by rendering at higher resolution and downsampling to a lower resolution, as per triangle splatting +"""
    resolution_schedule: int = 3000
    """training starts at 1/d resolution, every n steps this is doubled"""
    num_downscales: int = 0
    """At the beginning, resolution is 1/2^d, where d is this number"""

    # Densification
    use_mcmc: bool = True
    """Must be True for 3DGS methods."""
    preallocate_splat_tensors: bool = True
    """Whether to pre-allocate Gaussian attribute tensors to cap_max to avoid OOM during densification"""
    cap_max: int = 1_000_000
    """maximum number of splats, dataset-specific tuning required"""
    refine_every: int = 100
    """Densify every this number of steps"""
    refine_start_iter: int = 500
    """Start densification at this number of steps"""
    refine_stop_num_iter: int = 5000
    """Stop densification at this number of steps before maximum number of training iterations"""
    noise_lr: float = 5e5
    """Scalar for MCMC-style noise injection"""
    min_opacity: float = 0.005
    """Minimum Gaussian opacity before relocation"""
    growth_factor: float = 1.05
    """Multiply number of splats by this number at each densification step"""
    relocate_heuristic_weight: float = 1.0
    """Weight of gradient used in sampling Gaussians to relocate/add to.
        If 0.0, use only opacity; If 1.0, use other heuristics (see strategy/mcmc.py for details)."""
    use_edge_aware_score: bool = True
    """Whether to use edge aware score to guide densification.
        If True, it computes edge aware score following https://arxiv.org/abs/2603.08661
        Note that this is only active when relocate_heuristic_weight is nonzero"""
    use_loss_map: bool = True
    """Whether to use loss map to guide densification.
        Note that this is only active when relocate_heuristic_weight is nonzero."""
    use_long_axis_split: bool = True
    """whether to use long-axis split described in https://arxiv.org/abs/2508.12313 for relocation and sample add.
        When combined with relocate_heuristic_weight=1.0, this can give less blurry background details for unbounded outdoor scenes."""
    relocate_screen_size: float = float('inf')
    """if a gaussian is more than this fraction of screen space, relocate it
        Useful for fisheye with 3DGUT, may drop PSNR for conventional cameras
        For likely better quality, use max_screen_size instead"""
    max_screen_size: float = 0.3
    """if a gaussian is more than this fraction of screen space, clip scale and increase opacity
        Intended to be an MCMC-friendly alternative of relocate_screen_size"""
    max_screen_size_clip_hardness: float = 1.5
    """clip hardness for Gaussians with large screen space size, between 1 and infinity, larger is harder"""
    max_world_size: float = float('inf')
    """if a gaussian is more than this of world space, clip scale
        Useful if you see huge floaters at a distance in large indoor space"""
    reset_alpha_every: int = 30
    """Every this many refinement steps, reset the alpha. Only applies for opaque triangle splatting."""

    # Exposure/WB correction
    use_bilateral_grid: bool = True
    """If True, use bilateral grid to handle the ISP changes in the image space.
        This technique was introduced in the paper 'Bilateral Guided Radiance Field Processing' (https://bilarfpro.github.io/)."""
    bilagrid_shape: Tuple[int, int, int] = (16, 16, 8)
    """Shape of the bilateral grid, typically `16 16 8`, or `8 8 4` for scenes with low-texture surfaces."""
    bilagrid_type: Literal["affine", "ppisp", "loglinear"] = "loglinear"
    """What the bilateral grid predicts.
        affine: 4x3 matrix per original bilateral grid.
        ppisp: PPISP exposure and color parameters, generally gives less color shift but can be less numerically stable.
        loglinear: 3x3 linear transformation matrix with log-encoded diagonals, balances color shift and numerical stability.
    """
    use_bilateral_grid_for_geometry: bool = True
    """If True, use bilateral grid for depth and normal (e.g. AI generated biased ones)"""
    bilagrid_shape_geometry: Tuple[int, int, int] = (8, 8, 4)
    """Shape of the bilateral grid for depth and normal (X, Y, W)"""
    bilagrid_tv_loss_weight: float = 10.0
    """Total variation loss weight for bilateral grid used for radiance"""
    bilagrid_mean_reg_weight: float = 10.0
    """Regularization to discourage bilateral grid color shift"""
    optimize_bilagrid_frequencies: bool = False
    """Whether to optimize bilagrid parameters in frequency domain instead of time domain"""
    bilagrid_tv_loss_weight_geometry: float = 10.0
    """Total variation loss weight for bilateral grid used for geometry"""
    use_ppisp: bool = False
    """If True, use the PPISP model (https://research.nvidia.com/labs/sil/projects/ppisp/) to handle per-pixel color distortions."""
    ppisp_param_type: Literal["original", "rqs"] = "rqs"
    """Parameterization for PPISP. "original" implements the original paper,
        "rqs" uses a parameterization that is more friendly to optimization and can produce better results in darker areas."""
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

    # Linear and wide-gamut
    image_color_is_linear: bool = False
    """Whether to assume training images are in linear color space."""
    image_color_gamut: Literal[None, "ACES2065-1", "ACEScg", "Rec.2020", "AdobeRGB", "DCI-P3"] = None
    """Color gamut of input images. If None, Rec.709 will be used. Note that tonemap is not applied."""
    splat_color_is_linear: Literal[True, False, None] = None
    """Whether to train splats in linear color space. If None, will use same as images."""
    splat_color_gamut: Literal["Rec.709", "ACES2065-1", "ACEScg", "Rec.2020", "AdobeRGB", "DCI-P3", None] = None
    """Color gamut of trained splats. If None, will use same as images. Note that tonemap is not applied."""
    convert_initial_point_cloud_color: Literal[True, False, None] = None
    """If True, this will assume color in initial point cloud is sRGB, and convert if images are in a linear or wide-gamut color space."""

    # Regularization
    suppress_initial_scales: bool = False
    """Whether to suppress scales during initialization to discourage large floaters in vacant areas"""
    scale_regularization_weight: float = 0.0
    """If enabled, a scale regularization introduced in PhysGauss (https://xpandora.github.io/PhysGaussian/) is used for reducing huge spikey gaussians."""
    max_gauss_ratio: float = 10.0
    """Threshold of ratio of gaussian max to min scale before applying regularization loss from the PhysGaussian paper"""
    depth_distortion_reg: float = 0.0
    """Weight for depth distortion regularizer"""
    normal_distortion_reg: float = 0.0
    """Weight for normal distortion regularizer"""
    rgb_distortion_reg: float = 0.0
    """Weight for rgb distortion regularizer"""
    distortion_reg_warmup: int = 6000
    """warmup steps for depth regularizer, regularization weight ramps up"""
    normal_reg_weight: float = 0.04
    """Weight for normal regularizer"""
    normal_reg_warmup: int = 6000
    """warmup steps for normal regularizer, regularization weight ramps up"""
    alpha_reg_weight: float = 0.0
    """Weight for alpha regularizer (encourage alpha to go to either 0 or 1)"""
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
    mcmc_opacity_reg: float = 0.001  # 0.01 in original paper
    """Encourage low opacity to aid densification, per MCMC."""
    mcmc_scale_reg: float = 0.01  # 0.01 in original paper
    """Encourage low scale, per MCMC."""
    erank_reg: float = 0.0
    """erank regularization weight, for 3DGS only -
        see *Effective Rank Analysis and Regularization for Enhanced 3D Gaussian Splatting, Hyung et al.*"""
    erank_reg_s3: float = 0.0
    """erank regularization weight for smallest dimension, for 3DGS only"""
    quat_norm_reg: float = 0.01
    """Weight to regularize quaternion norm to identity"""

    # supervision using a foundation depth model
    # enable these by setting `depth_model` in data manager config
    supervision_warmup: int = 0
    """Start using foundation model depth at this number of steps"""
    depth_supervision_weight: float = 0.0
    """Weight for depth supervision by comparing rendered depth with depth predicted by a foundation model
        Warn that this can reduce quality if AI generated depth is heavily biased"""
    normal_supervision_weight: float = 0.01
    """Weight for normal supervision by comparing normal from rendered depth with normal from depth predicted by a foundation model"""

    # Validation
    overfit_score_aggregation_mode: Literal['max', 'min', 'mean'] = 'min'
    """Mode to aggregate multiple overfitting objectives.
        Use max for more aggressive early stopping, min for more conservative early stopping, and mean for something in between."""
    validation_loss_average_window: int = 500
    """Window to calculate moving average validation loss for early stop"""
    early_stop_patience: int = 1000
    """Stop training if overfitting score remains positive for this number of iterations"""
    early_stop_warmup: int = 12000
    """Warmup steps for early stop, will not early stop before this number of steps
        Recommend setting this number no less than regularization warmups"""


class SpirulaeSplatModel(torch.nn.Module):
    """Template Model."""

    config: SpirulaeSplatModelConfig

    def __init__(
        self,
        trainer_config: 'spirulae_splat.modules.trainer.TrainerConfig',
        seed_points: Optional[Tuple[torch.Tensor, torch.Tensor]] = None,
        cameras: Optional[Cameras] = None,
    ):
        super().__init__()

        self.trainer_config = trainer_config
        self.config = trainer_config.model  # type: SpirulaeSplatModelConfig

        self.seed_points = (seed_points['points3D_xyz'], seed_points['points3D_rgb'])
        self.cameras = cameras

        self.num_train_data = len(self.cameras)

        if CameraType.EQUIRECTANGULAR.value in cameras.camera_type:
            assert len(set(cameras.camera_type)) == 1, "Mixed equirectangular and pinhole/fisheye is not supported"
            self.num_train_data *= 6  # TODO

        self.info = {}

        self.populate_modules()

        if self.config.primitive in ['3dgs', 'mip', '3dgut']:
            splat_params = (
                self.means, self.quats, self.scales,
                self.opacities,
                self.features_dc, self.features_sh
            )
        if self.config.primitive in ['3dgut_sv']:
            splat_params = (
                self.means, self.quats, self.scales,
                self.opacities,
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

        self.renderer = Renderer(
            self.config.primitive,
            splat_params,
            self.seed_points[0].shape[0]  # TODO: voxel
        )

    def populate_modules(self):
        if self.seed_points is not None:
            if len(self.seed_points[0]) > self.config.cap_max:
                indices = torch.randperm(len(self.seed_points[0]))[:self.config.cap_max]
                self.seed_points = [t[indices] for t in self.seed_points]
            means = self.seed_points[0]
            if self.config.relative_scale is not None:
                means *= self.config.relative_scale
        else:
            raise ValueError("No seed point found")
        self.xys_grad_norm = None
        self.ch_grad_norm = None
        self.max_2Dsize = None

        scale_init, opacity_init = 0.1, 0.5  # per MCMC paper
        if self.config.use_mcmc and self.config.max_screen_size < 1.0:
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
            from time import perf_counter
            time0 = perf_counter()
            distances, indices = self.k_nearest_neighbor(means.data, 4)
            time1 = perf_counter()
            # print("k_nearest_neighbor time:", time1-time0)
            num_points = means.shape[0]
            scales = torch.sqrt(torch.from_numpy(distances**2).mean(-1, keepdim=True)).repeat(1, 3).float()
            scales = torch.log(scale_init * scales / (self.config.kernel_radius/3.0) + 1e-8)
            if self.config.suppress_initial_scales:
                scales = self.suppress_initial_scales(means, scales)
            quats = F.normalize(torch.randn((num_points, 4)), dim=-1)
            opacities = torch.logit(opacity_init * torch.ones(num_points, 1))

        else:
            distances, indices = self.k_nearest_neighbor(means.data, 6)
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
            # We can have colors without points.
            and self.seed_points[1].shape[0] > 0
        ):
            seed_color = self.seed_points[1] / 255
            if len(seed_color) != num_points:
                seed_color = seed_color.mean(0, True).repeat(num_points, 1)
        else:
            seed_color = torch.rand(num_points, 3)

        if self.config.splat_color_gamut is None:
            self.config.splat_color_gamut = self.config.image_color_gamut
        elif self.config.convert_initial_point_cloud_color is None:
            self.config.convert_initial_point_cloud_color = True
        if self.config.splat_color_gamut == "Rec.709":
            self.config.splat_color_gamut = None
        if self.config.splat_color_is_linear is None:
            self.config.splat_color_is_linear = self.config.image_color_is_linear
        elif self.config.convert_initial_point_cloud_color is None:
            self.config.convert_initial_point_cloud_color = True
        if self.config.convert_initial_point_cloud_color is None:
            self.config.convert_initial_point_cloud_color = False
        if self.config.convert_initial_point_cloud_color:
            if self.config.splat_color_is_linear or self.config.splat_color_gamut is not None:
                seed_color = torch.where(seed_color < 0.055, seed_color / 12.92, torch.relu((seed_color + 0.055) / 1.055) ** 2.4)
            if self.config.splat_color_gamut is not None:
                color_transform = get_color_transform_matrix(self.config.image_color_gamut, device=seed_color.device)
                seed_color = seed_color @ torch.linalg.inv(color_transform).T
                if not self.config.splat_color_is_linear:
                    seed_color = torch.where(seed_color < 0.0031308, 12.92 * seed_color, 1.055 * torch.relu(seed_color) ** (1/2.4) - 0.055)

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

        new_num_points = num_points
        if self.config.use_mcmc and self.config.preallocate_splat_tensors:
            new_num_points = max(num_points, self.config.cap_max)
            for key, value in gauss_params.items():
                if not isinstance(value, torch.Tensor) or value.shape[0] != num_points:
                    continue
                if (value == 0).all():
                    value = torch.zeros(new_num_points, *value.shape[1:], device=value.device, dtype=value.dtype)
                else:
                    value = torch.concat((value, torch.zeros(new_num_points-num_points, *value.shape[1:], device=value.device, dtype=value.dtype)))
                gauss_params[key] = torch.nn.Parameter(value)

        self.gauss_params = torch.nn.ParameterDict(gauss_params)

        optim_info = {}
        if self.config.use_mcmc and self.config.preallocate_splat_tensors:
            optim_info['num_splats'] = num_points
        for key, value in self.gauss_params.items():
            if isinstance(value, torch.Tensor) and value.shape[0] == new_num_points:
                value.optim_info = {**optim_info}

        if self.config.use_camera_optimizer:
            self.camera_optimizer = self.config.camera_optimizer.setup(
                num_cameras=self.num_train_data, device="cpu"
            )

        self.training_losses = SplatTrainingLosses(self.config, self.num_train_data)

        self.step = 0

        if self.config.background_color == "gray":
            self.background_color = torch.tensor([0.5, 0.5, 0.5])
        elif self.config.background_color == "black":
            self.background_color = torch.tensor([0.0, 0.0, 0.0])
        elif self.config.background_color == "white":
            self.background_color = torch.tensor([1.0, 1.0, 1.0])
        if self.config.splat_color_is_linear:
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
            final_num = self.config.cap_max

            min_warmup_steps = self.config.refine_every * self.config.reset_alpha_every
            num_steps_until_full = math.log(max(final_num, current_num) / current_num) / math.log(self.config.growth_factor) * \
                self.config.refine_every + (self.config.refine_start_iter + self.config.refine_every)
            warmup_steps_0 = min(2.0*max(min_warmup_steps, num_steps_until_full), self.trainer_config.num_iterations/3)

            self.strategy = OpaqueStrategy(
                cap_max=self.config.cap_max,
                noise_lr=self.config.noise_lr,
                refine_start_iter=self.config.refine_start_iter,
                warmup_steps=warmup_steps_0,
                refine_stop_iter=self.trainer_config.num_iterations,
                refine_every=self.config.refine_every,
                grow_factor=self.config.growth_factor,
                min_opacity=self.config.min_opacity,
                relocate_scale2d=self.config.relocate_screen_size,
                max_scale2d=self.config.max_screen_size,
                max_scale2d_clip_hardness=self.config.max_screen_size_clip_hardness,
                max_scale3d=self.config.max_world_size,
                geometry_optimizer_stop_iter=self.trainer_config.num_iterations-self.config.refine_stop_num_iter
            )
            self.strategy_state = self.strategy.initialize_state()
            return

        # MCMC mode
        if self.config.use_mcmc:
            current_num = len(self.means)
            final_num = self.config.cap_max
            warpup_steps = math.log(max(final_num, current_num) / current_num) / math.log(self.config.growth_factor) * \
                self.config.refine_every + (self.config.refine_start_iter + self.config.refine_every)
            self.mcmc_num_steps_until_full = warpup_steps

            self.strategy = MCMCStrategy(
                cap_max=self.config.cap_max,
                noise_lr=self.config.noise_lr,
                refine_start_iter=self.config.refine_start_iter,
                refine_stop_iter=self.trainer_config.num_iterations-self.config.refine_stop_num_iter,
                refine_every=self.config.refine_every,
                grow_factor=self.config.growth_factor,
                min_opacity=self.config.min_opacity,
                prob_grad_weight=self.config.relocate_heuristic_weight,
                use_long_axis_split=self.config.use_long_axis_split,
                is_3dgs=True,
                relocate_scale2d=self.config.relocate_screen_size,
                max_scale2d=self.config.max_screen_size,
                max_scale2d_clip_hardness=self.config.max_screen_size_clip_hardness,
                max_scale3d=self.config.max_world_size,
            )
            self.strategy_state = self.strategy.initialize_state()
            return

    def _get_gauss_param(self, key: str):
        tensor = self.gauss_params.get(key, None)
        return tensor
        if hasattr(tensor, "optim_info") and 'num_splats' in tensor.optim_info:
            if len(tensor) <= tensor.optim_info['num_splats']:
                return tensor
            result = tensor[:tensor.optim_info['num_splats']]
            if hasattr(tensor, 'optim_info'):
                result.optim_info = tensor.optim_info
            return result
        return tensor

    @property
    def colors(self):
        if self.config.primitive == "3dgut_sv":
            return self.sv_colors.mean(-2) + 0.5
        C0 = 0.28209479177387814
        return self.features_dc * C0 + 0.5

    @property
    def num_points(self):
        param = self.means if self.means is not None else self.features_dc
        if hasattr(param, 'optim_info') and 'num_splats' in param.optim_info:
            return param.optim_info['num_splats']
        return param.shape[0]

    @property
    def means(self):
        return self._get_gauss_param("means")

    @property
    def scales(self):
        return self._get_gauss_param("scales")
    
    @property
    def quats(self):
        return self._get_gauss_param("quats")

    @property
    def features_dc(self):
        return self._get_gauss_param("features_dc")

    @property
    def features_sh(self):
        return self._get_gauss_param("features_sh")

    @property
    def features_ch(self):
        return self._get_gauss_param("features_ch")

    @property
    def sv_sites(self):
        return self._get_gauss_param("sv_sites")

    @property
    def sv_colors(self):
        return self._get_gauss_param("sv_colors")

    @property
    def opacities(self):
        param = self._get_gauss_param("opacities")
        # # TODO: properly do this in projection
        # if hasattr(param, 'optim_info') and 'num_splats' in param.optim_info:
        #     num_splats = param.optim_info['num_splats']
        #     if num_splats < param.shape[0]:
        #         param.data[num_splats:] = -10.0  # logit(1/255) ~ -5.54
        return param

    @property
    def densities(self):
        # TODO: num_splats
        return self._get_gauss_param("densities")

    def load_state_dict(self, dict, **kwargs):  # type: ignore
        # resize the parameters to match the new number of points
        self.step = self.trainer_config.num_iterations
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

    def k_nearest_neighbor(self, x: torch.Tensor, k: int):
        """
        Find k-nearest neighbors using scipy.spatial.cKDTree.
        Uses all available CPU cores.
        """
        x_np = x.detach().cpu().numpy()

        # Build the tree
        # balanced_tree=False can speed up build time for large datasets
        from scipy.spatial import cKDTree
        tree = cKDTree(x_np, balanced_tree=False)

        # Query the tree
        # workers=-1 uses all available CPU cores
        distances, indices = tree.query(x_np, k=k+1, workers=-1)

        # Exclude the point itself (first column)
        return distances[:, 1:].astype(np.float32), indices[:, 1:].astype(np.int64)

    def suppress_initial_scales(self, means: torch.Tensor, scales: torch.Tensor):
        assert CameraType.EQUIRECTANGULAR.value not in self.cameras, "TODO: Not Implemented"

        R = self.cameras.camera_to_worlds[:, :3, :3]  # 3 x 3
        T = self.cameras.camera_to_worlds[:, :3, 3:4]  # 3 x 1
        if self.config.relative_scale is not None:
            T = T * self.config.relative_scale
        R = R * torch.tensor([[[1.0, -1.0, -1.0]]])
        R_inv = R.transpose(-1, -2)
        T_inv = -torch.bmm(R_inv, T)
        viewmats = torch.eye(4, dtype=self.cameras.camera_to_worlds.dtype)[None].repeat(len(self.cameras), 1, 1)
        viewmats[:, :3, :3] = R_inv
        viewmats[:, :3, 3:4] = T_inv

        log_scales = -_make_lazy_cuda_func("cov_scale_init")(
            *[x.cuda() if x is not None else x for x in (
                means,
                self.cameras.camera_type == CameraType.FISHEYE.value,
                torch.concatenate((self.cameras.width, self.cameras.height), dim=-1).int(),
                torch.concatenate((self.cameras.fx, self.cameras.fy, self.cameras.cx, self.cameras.cy), dim=1),
                viewmats,
                self.cameras.distortion_params
            )]
        )

        original_mean = scales.mean().item()

        log_scales += (original_mean - log_scales.mean().item())
        log_scales = log_scales.repeat(1, 3)
        # return log_scales

        scales = torch.fmin(log_scales, scales.cuda())
        scales += (original_mean - scales.mean().item())
        return scales

    # def step_cb(self, optimizers: Dict, step: int):
    #     self.step = step
    #     self.optimizers = optimizers

    def step_cb(self, step: int):
        self.step = step

        self.info = {}
        self.renderer.zero_grad()

    def step_post_backward(self):
        return  # TODO
        if self.config.primitive == 'voxel':
            raise NotImplementedError()
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
                lr=get_scheduled_lr(self.optimizers['means'].param_groups[0], self.step+1),
                packed=(self.config.packed or self.config.use_bvh),
            )
        else:
            raise NotImplementedError()
            self.strategy.step_post_backward(
                params=self.gauss_params,
                optimizers=self.optimizers,
                state=self.strategy_state,
                step=self.step,
                info=self.info,
                packed=(self.config.packed or self.config.use_bvh),
            )

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
            gps["bilateral_grid_depth"] = list(self.training_losses.bil_grids_depth.parameters())
            gps["bilateral_grid_normal"] = list(self.training_losses.bil_grids_normal.parameters())
        if self.config.use_ppisp:
            gps["ppisp"] = [self.training_losses.ppisp_params]
        gps['_dummy'] = [self.training_losses._dummy]

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

    def blend_background(
            self, camera: Cameras, c2w: torch.Tensor, intrins: torch.Tensor,
            W: Optional[int] = None, H: Optional[int] = None, is_fisheye: Optional[bool] = None,
            rgb: Optional[torch.Tensor] = None, transmittance: Optional[torch.Tensor] = None
        ):
        if not isinstance(camera, Cameras):
            print("Called blend_background with not a camera")
            return {}

        if W is None or H is None:
            W, H = int(camera.width[0].item()), int(camera.height[0].item())
        if is_fisheye is None:
            is_fisheye = (camera.camera_type[0] == "FISHEYE")
        sh_degree = self.config.background_sh_degree

        if self.config.randomize_background == True:
            # return torch.rand_like(self.background_color).repeat(H, W, 1)
            if rgb is None or transmittance is None:
                return None
            randomize_weight = min(self.step / max(self.config.randomize_background_warmup, 1), 1)
            randomize_weight = 1.0 - (1.0 - self.config.randomize_background_pre_warmup) * (1.0 - randomize_weight)
            return blend_background_noise(self.config.splat_color_is_linear, rgb, transmittance, randomize_weight)

        elif not self.config.train_background_color and self.config.background_color == "black":
            if rgb is None or transmittance is None:
                return None
            return rgb.clip_(0.0, 1.0)

        elif not self.config.train_background_color or not (sh_degree > 0):
            background = self.background_color[None].repeat(len(camera), H, W, 1)

        else:
            sh_coeffs = torch.cat((self.background_color.unsqueeze(0), self.background_sh), dim=0)  # [(deg+1)^2, 3]
            background = render_background_sh(
                W, H,
                "FISHEYE" if is_fisheye else "PINHOLE",
                intrins, c2w[..., :3, :3], sh_degree, sh_coeffs
            )

        if rgb is None or transmittance is None:
            return background
        return blend_background(rgb, transmittance, background)

    @staticmethod
    def get_empty_outputs(width: int, height: int, background: torch.Tensor) -> Dict[str, Union[torch.Tensor, List]]:
        return {}

    def forward(self):
        raise NotImplementedError()

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
            camera.rescale(1 / camera_downscale)

        # TODO: separate different sizes/intrins
        W, H = int(camera.width[0].item()), int(camera.height[0].item())
        is_fisheye = (camera.camera_type[0] == "FISHEYE")

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

        intrins = camera.intrins
        
        kwargs = {'actual_width': W, 'actual_height': H}
        if 'actual_width' in camera.metadata:
            kwargs['actual_width'] = int(camera.metadata['actual_width'] + 0.5)
        if 'actual_height' in camera.metadata:
            kwargs['actual_height'] = int(camera.metadata['actual_height'] + 0.5)
        if 'patch_offsets' in camera.metadata:
            kwargs['patch_offsets'] = camera.metadata['patch_offsets']
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
            if visibility_masks is None:
                visibility_masks = torch.zeros((len(camera), H, W, 1), dtype=torch.bool, device=device)
            visibility_masks = self._downscale_if_required(visibility_masks)
            kwargs['max_blending_masks'] = visibility_masks.bool()

        if 'accum_weight_map' in camera.metadata:
            accum_weight_map = camera.metadata['accum_weight_map']
            if accum_weight_map is not None:
                accum_weight_map = self._downscale_if_required(accum_weight_map)
                kwargs['accum_weight_map'] = accum_weight_map

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
            im = im[:, :H, :W, :]
            return im

        # if not self.training:
        #     viewmats_0, intrins_0 = viewmats, intrins
        #     viewmats, intrins = split_into_tiles(viewmats, intrins)
        #     if 'dist_coeffs' in kwargs:
        #         kwargs['dist_coeffs'] = kwargs['dist_coeffs'].repeat(gh*gw, 1, 1)
        #     W = TILE_SIZE
        #     H = TILE_SIZE

        # TODO
        # if val:
        #     splat_params = tuple([(p.detach() if isinstance(p, torch.Tensor) else p) for p in splat_params])

        # rendering
        use_bvh = self.config.use_bvh and self.training and not val
        self.renderer.set_params(
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
            camera_model=["pinhole", "fisheye"][is_fisheye],
            output_distortion=any([c != 0.0 for c in self.training_losses.get_2dgs_reg_weights()[0]]),
            compute_hessian_diagonal=self.config.compute_hessian_diagonal,
            **kwargs,
        )
        self.renderer.forward()
        rgbd = self.renderer.render_colors
        Ts = self.renderer.render_Ts
        meta = self.renderer.meta
        # torch.cuda.empty_cache()
        if self.config.supersampling != 1:
            rgbd = [resize_image(im, self.config.supersampling) for im in rgbd]
            Ts = resize_image(Ts, self.config.supersampling)

        # if not self.training:
        #     W = gw*TILE_SIZE
        #     H = gh*TILE_SIZE
        #     # print(rgbd.shape, Ts.shape)
        #     rgbd = [merge_tiles(comp) for comp in rgbd]
        #     Ts = merge_tiles(Ts)
        #     for key in ['rgb_distortion', 'depth_distortion', 'normal_distortion']:
        #         if key in meta:
        #             meta[key] = merge_tiles(meta[key])
        #     viewmats, intrins = viewmats_0, intrins_0
        #     if 'dist_coeffs' in kwargs:
        #         kwargs['dist_coeffs'] = kwargs['dist_coeffs'][:1]

        depth_im_ref = rgbd[1]

        # normals
        render_normal = None
        if len(rgbd) > 2:
            render_normal = torch.where(Ts < 1.0, F.normalize(rgbd[2], dim=-1), rgbd[2])

        depth_normal = None

        rgb = rgbd[0]

        radii = meta["radii"]
        depths = meta["depths"]

        if self.training:
            radii = scatter_max(self.opacities.flatten(), radii, meta["gaussian_ids"])
            if 'radii' in self.info:
                self.info['radii'] = torch.fmax(self.info['radii'], radii)
            else:
                self.info['radii'] = radii

            self.info.update({
                "width": kwargs["actual_width"],
                "height": kwargs["actual_height"],
                "n_cameras": len(camera),
                "n_train": self.num_train_data,
                "means2d": self.means,
                "depths": depths,
                "backward_info": self.renderer.backward_info,
            })
            if 'patch_offsets' in camera.metadata:
                self.info['patch_offsets'] = camera.metadata['patch_offsets']
            if 'actual_images_per_batch' in camera.metadata:
                self.info['n_cameras'] = camera.metadata['actual_images_per_batch']
            for key in ['gaussian_ids', 'camera_ids', 'max_blending', 'bvh_time']:
                if key in meta:
                    self.info[key] = meta[key]

        # blend with background
        rgb = self.blend_background(camera, optimized_camera_to_world, intrins, W, H, is_fisheye, rgb, Ts)

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
        if self.training:
            outputs["backward_info"] = self.renderer.backward_info
        if depth_im_ref is not None:
            outputs["depth"] = depth_im_ref
        if render_normal is not None:
            outputs["normal"] = render_normal
        if depth_normal is not None:
            outputs["depth_normal"] = depth_normal

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
            outputs["alpha"] = (1.0 - Ts).reshape((H, W, 1)).repeat(1, 1, 3)
            background = self.blend_background(camera, optimized_camera_to_world, intrins, W, H, is_fisheye)
            if background is not None:
                outputs["background"] = torch.clip(background, min=0.0, max=1.0)
        else:
            outputs["transmittance"] = Ts

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
            if self.config.splat_color_is_linear or self.config.splat_color_gamut != None:
                outputs["rgb_raw"] = outputs["rgb"]
                for key in ['rgb', 'background']:
                    if key not in outputs:
                        continue
                    if self.config.splat_color_is_linear or self.config.splat_color_gamut != None:
                        color_matrix = get_color_transform_matrix(self.config.splat_color_gamut)
                        outputs[key] = rgb_to_srgb(outputs[key], self.config.splat_color_is_linear, color_matrix).clip(0, 1)
            for key in outputs:
                outputs[key] = outputs[key].squeeze(0)

        if self.training:
            kwargs["intrins"] = intrins
            kwargs["camera_model"] = ["pinhole", "fisheye"][is_fisheye]
            outputs["camera"] = camera
            outputs["camera_intrins"] = kwargs
        # if not self.training and True:
        #     outputs['ray'] = merge_tiles(meta['intersection_count'].float().reshape(-1, 1, 1, 1).repeat(1, TILE_SIZE, TILE_SIZE, 1))

        # setup override required for optimizers
        # TODO: more reliable way than setattr tensor?
        if self.training:
            optim_info = {
                'radii': self.info['radii'],
            }
            for key, value in self.gauss_params.items():
                if isinstance(value, torch.Tensor) and key in self.gauss_params:
                    optim_info[key] = value
            for key, value in self.gauss_params.items():
                if isinstance(value, torch.Tensor) and key in self.gauss_params:
                    if not hasattr(value, 'optim_info'):  # can happen during eval
                        value.optim_info = {}
                    value.optim_info.update(optim_info)

            if self.config.relative_scale is None and "means" in self.gauss_params:
                self.means.optim_info['optimizer_override'] = "fused_adam_scale_agnostic_mean"
            if "quats" in self.gauss_params:
                self.quats.optim_info['optimizer_override'] = "fused_adam_riemannian_quat"
            if self.config.splat_color_is_linear:
                if "features_dc" in self.gauss_params and 'opacities' in self.gauss_params:
                    # self.features_dc.optim_info['optimizer_override'] = "fused_adam_linear_rgb_optim"
                    self.features_dc.optim_info['optimizer_override'] = "fused_adamtr_linear_rgb_optim"
                if "features_sh" in self.gauss_params and "features_dc" in self.gauss_params and 'opacities' in self.gauss_params:
                    self.features_sh.optim_info['optimizer_override'] = "fused_adamtr_linear_rgb_sh_optim"
            elif self.config.splat_color_gamut is not None:
                if "features_dc" in self.gauss_params and 'opacities' in self.gauss_params:
                    self.features_dc.optim_info['optimizer_override'] = "fused_adamtr_rgb_optim"
                if "features_sh" in self.gauss_params and "features_dc" in self.gauss_params and 'opacities' in self.gauss_params:
                    self.features_sh.optim_info['optimizer_override'] = "fused_adamtr_rgb_sh_optim"
            if self.config.optimizer_offload == "sh" and "features_sh" in self.gauss_params:
                if hasattr(self.features_sh, "optimizer_override"):
                    raise ValueError("Optimizer offloading is not supported for linear color space")
                self.features_sh.optim_info['optimizer_offload'] = True
            elif self.config.optimizer_offload == "all":
                for key, value in self.gauss_params.items():
                    if isinstance(value, torch.Tensor) and not (hasattr(value, 'optim_info') and "optimizer_override" in value.optim_info):
                        value.optim_info['optimizer_offload'] = True
                if self.config.use_bilateral_grid:
                    self.training_losses.bil_grids.grids.optim_info = {'optimizer_offload': True}
                if self.config.use_bilateral_grid_for_geometry:
                    self.training_losses.bil_grids_depth.grids.optim_info = {'optimizer_offload': True}
                    self.training_losses.bil_grids_normal.grids.optim_info = {'optimizer_offload': True}

        # return outputs

        # Debug densification
        if not self.training and self.step > 1:
        # if self.step > 1:

            param_to_vis = self.renderer.densify_accum_buffer[:, 0]
            param_to_vis = param_to_vis / param_to_vis.mean()
            # param_to_vis = torch.log10(param_to_vis + 1e-30)

            indices = _make_lazy_cuda_func("weighted_sample_without_replacement")(
                -1, self.renderer.densify_accum_buffer, None, len(param_to_vis) // 20, self.step
            )
            param_to_vis = torch.zeros_like(param_to_vis)
            param_to_vis[indices] = 1

            param_to_vis = (param_to_vis - 0.5) / 0.28
            param_to_vis = param_to_vis.unsqueeze(-1).repeat(1, 3)
            param_to_vis = torch.concatenate((param_to_vis, torch.zeros_like(self.means[len(param_to_vis):])), 0)

            old_splats_world = self.renderer.splats_world
            self.renderer.splats_world = (
                self.means, self.quats, self.scales,
                self.opacities,
                param_to_vis, torch.zeros_like(self.features_sh)
            )
            self.renderer.set_params(
                viewmats=viewmats,  # [C, 4, 4]
                intrins=intrins * self.config.supersampling,  # [C, 4]
                width=W * self.config.supersampling,
                height=H * self.config.supersampling,
                packed=(self.config.packed or use_bvh),
                use_bvh=(use_bvh),
                relative_scale=self.config.relative_scale,
                camera_model=["pinhole", "fisheye"][is_fisheye],
                output_distortion=any([c != 0.0 for c in self.training_losses.get_2dgs_reg_weights()[0]]),
                compute_hessian_diagonal=self.config.compute_hessian_diagonal,
                **kwargs,
            )
            self.renderer.forward()
            rgbd = self.renderer.render_colors
            Ts = self.renderer.render_Ts
            meta = self.renderer.meta
            self.renderer.splats_world = old_splats_world
            outputs['refinement_score'] = rgbd[0][0, :, :, :].mean(dim=-1, keepdim=True)

        return outputs

    def backward(self, outputs, loss_grads):
        if not self.training or loss_grads is None:
            return
        if 'backward_info' not in outputs:
            return

        v_render_rgb = loss_grads[0]
        v_render_depth = loss_grads[2] if len(loss_grads) > 2 else None
        v_render_normal = loss_grads[4] if len(loss_grads) > 4 else None
        v_render_Ts = loss_grads[7] if len(loss_grads) > 7 else None
        v_rgb_distortion = loss_grads[8] if len(loss_grads) > 8 else None
        v_depth_distortion = loss_grads[9] if len(loss_grads) > 9 else None
        v_normal_distortion = loss_grads[10] if len(loss_grads) > 10 else None

        if v_render_Ts is None:
            v_render_Ts = torch.zeros_like(outputs['transmittance'])

        if self.config.primitive in ['opaque_triangle']:
            render_grads = (v_render_rgb, v_render_depth, v_render_normal)
        else:
            render_grads = (v_render_rgb, v_render_depth)

        self.renderer.backward(
            render_grads,
            v_render_Ts,
            v_rgb_distortion,
            v_depth_distortion,
            v_normal_distortion,
        )

    def optim_step(self):
        self.renderer.optim_step(self.step, self.config, self.trainer_config.optimizer)
        if self.step < self.trainer_config.num_iterations-self.config.refine_stop_num_iter:
            self.renderer.densify_step(self.step, self.config, None)
        # self.step_post_backward()

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

    def get_loss_dict(self, outputs, batch, batch_size: int, no_static_losses=False) -> Dict[str, torch.Tensor]:
        """Computes and returns the losses dict.

        Args:
            outputs: the output to compute loss dict to
            batch: ground truth batch corresponding to outputs
            metrics_dict: dictionary of metrics, some of which we can use for loss
        """

        static_losses = {}
        if not no_static_losses:

            # Per-splat losses (CUDA async)
            self.training_losses.get_static_losses(
                self.step,
                self.quats, self.scales, self.opacities,
                static_losses,
                self.info['backward_info']
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
            loss_dict, loss_grad = self.training_losses(self.step, batch[0], outputs[0], self.info)
            # for key, value in outputs[1].items():
            #     if isinstance(value, torch.Tensor):
            #         outputs[1][key] = value.detach()
            val_loss_dict, val_loss_grad = self.training_losses(self.step, batch[1], outputs[1], self.info, val=True)
            total_val_loss = torch.stack([
                x for x in val_loss_dict.values() if isinstance(x, torch.Tensor)
            ]).sum()
        else:
            loss_dict, loss_grad = self.training_losses(self.step, batch, outputs, self.info)

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

        # Apply scalar in batching mode
        if batch_size != 1:
            for key, value in loss_dict.items():
                if isinstance(value, torch.Tensor):
                    loss_dict[key] = value * (1 / batch_size)

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

        if not no_static_losses:
            # TODO: print averaged loss dict in split_batch mode
            self.print_loss_dict(loss_dict)

        return loss_dict, loss_grad

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
        mem_stats = boldcyan(f"{used:.2f}") + f"\N{ZERO WIDTH SPACE}GiB " + boldcyan(f"{used_percentage:.0f}") + "%"

        if hasattr(self, 'total_val_loss'):
            losses['val_total'] = self.total_val_loss
        if hasattr(self, 'total_train_loss'):
            losses['train_total'] = self.total_train_loss
        if hasattr(self, 'total_val_lpips'):
            losses['lpips_val'] = self.total_val_lpips
        if hasattr(self, 'overfit_score'):
            losses['overfit_score'] = self.overfit_score

        def fmt(key: str, s: float, decimals=None, sigfigs=3) -> str:
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

            if (abs(l) < 0.1**(sigfigs+1) or abs(l) >= 1.0-(0.1**sigfigs)) and decimals is None:
                s = f"{{:.{sigfigs}g}}".format(l)
            else:
                if decimals is None:  # 3 sig figs
                    decimals = int(max(-math.log10((0.1**sigfigs)*abs(l)), 0)) if l != 0.0 else 0
                s = f"{{:.{decimals}f}}".format(l)
            if s.startswith('0.'):
                s = s[1:]
            if s.startswith('-0.'):
                s = '-' + s[2:]
            return boldcyan(s)

        reg_mcmc = (self.config.use_mcmc and self.step < self.trainer_config.num_iterations-self.config.refine_stop_num_iter)
        reg_2dgs = self.training_losses.get_2dgs_reg_weights()
        opacity_floor = (
            f"{bracket('OpacFloor')} {self.strategy.get_opacity_floor(self.step):.3f} "
            f"{bracket('Hardness')} {self.strategy.get_hardness(self.step):.3f}"
        ).replace('0.', '.') \
            if self.config.primitive == "opaque_triangle" else ""
        chunks = [
            f"{boldcyan(self.step)} "
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
            f"{bracket('SplatReg')} {orange('o')}={fmt('mcmc_opacity_reg', self.config.mcmc_opacity_reg * reg_mcmc, 3)} "
            f"{orange('s')}={fmt('mcmc_scale_reg', self.config.mcmc_scale_reg * reg_mcmc, 4)} "
            f"{orange('q')}={fmt('quat_norm_reg', self.config.quat_norm_reg, sigfigs=1)} "
            f"{orange('erank')}={fmt('erank_reg', max(self.config.erank_reg_s3, self.config.erank_reg))} "
            f"{orange('aniso')}={fmt('scale_reg', self.config.scale_regularization_weight)}",
            "                \n",
            f"{bracket('Bilagrid')} {orange('reg')}={fmt('bilagrid_mean_reg', self.config.bilagrid_mean_reg_weight)} "
            f"{orange('tv.rgb')}={fmt('tv_loss', self.config.bilagrid_tv_loss_weight)} "
            f"{orange('depth')}={fmt('tv_loss_depth', self.config.bilagrid_tv_loss_weight_geometry)} "
            f"{orange('normal')}={fmt('tv_loss_normal', self.config.bilagrid_tv_loss_weight_geometry)}"
            "                \n",
            f"{bracket('PPISP')} {orange('eμ')}={fmt('ppisp_reg_exposure_mean', self.config.ppisp_reg_exposure_mean)} "
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

    def set_gradient_accumulation_steps(self, gradient_accumulation_step: int, _trainer=[]):
        if len(trainer) == 0:
            # TODO: some less hacky way?
            import inspect
            trainer = inspect.stack()[4].frame.f_locals['self']
            _trainer.append(trainer)
        else:
            trainer = _trainer[0]
        trainer.gradient_accumulation_steps  # type: dict[str, int]

    @torch.no_grad()
    def get_outputs_for_camera(self, camera: Cameras) -> Dict[str, torch.Tensor]:
        """Takes in a camera, generates the raybundle, and computes the output of the model.
        Overridden for a camera-based gaussian model.

        Args:
            camera: generates raybundle
        """
        assert camera is not None, "must provide camera to gaussian model"
        outs = self.get_outputs(camera.to(self.device))
        return outs  # type: ignore

    @torch.inference_mode()
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
        # TODO: linear and wide-gamut color spaces
        gt_rgb = gt_rgb[..., :3]  # TODO: RGBA
        gt_rgb = gt_rgb.squeeze(0)

        predicted_rgb = outputs["rgb"]
        predicted_rgb = torch.clip(predicted_rgb, 0.0, 1.0)

        combined_rgb = torch.cat([gt_rgb, predicted_rgb], dim=1)

        # Switch images from [H, W, C] to [1, C, H, W] for metrics computations
        gt_rgb = torch.moveaxis(gt_rgb, -1, 0)[None, ...]
        predicted_rgb = torch.moveaxis(predicted_rgb, -1, 0)[None, ...]

        # metrics
        from torchmetrics.image import PeakSignalNoiseRatio
        from torchmetrics.image.lpip import LearnedPerceptualImagePatchSimilarity
        from pytorch_msssim import SSIM
        from torchmetrics.image import StructuralSimilarityIndexMeasure  # used by gsplat, not numerically identical to above

        if not hasattr(self, 'psnr'):
            self.psnr = PeakSignalNoiseRatio(data_range=1.0).cuda()
            self.ssim = SSIM(data_range=1.0, size_average=True, channel=3).cuda()
            self.ssim_torchmetrics = StructuralSimilarityIndexMeasure(data_range=1.0).cuda()
            self.lpips_alex = LearnedPerceptualImagePatchSimilarity(net_type="alex", normalize=True).cuda()
            self.lpips_vgg = LearnedPerceptualImagePatchSimilarity(net_type="vgg", normalize=False).cuda()

        l1 = torch.abs(gt_rgb - predicted_rgb).mean()
        psnr = self.psnr(gt_rgb, predicted_rgb)
        ssim = self.ssim(gt_rgb, predicted_rgb)
        ssim_torchmetrics = self.ssim_torchmetrics(gt_rgb, predicted_rgb)
        lpips_alex = self.lpips_alex(gt_rgb, predicted_rgb)
        lpips_vgg = self.lpips_vgg(gt_rgb, predicted_rgb)

        metrics_dict = {
            "l1": float(l1),
            "psnr": float(psnr),
            "ssim_pytorch_msssim": float(ssim),
            "ssim_torchmetrics": float(ssim_torchmetrics),
            "lpips_alex": float(lpips_alex),
            "lpips_vgg": float(lpips_vgg),
            "gaussian_count": float(self.num_points),
        }

        images_dict = {"img": combined_rgb}

        return metrics_dict, images_dict
