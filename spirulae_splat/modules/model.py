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
from spirulae_splat.modules.core import Renderer

from spirulae_splat.modules.verbose import TrainingVerbose
from spirulae_splat.modules._profile import PROFILE_TRAIN_STEP
from spirulae_splat.splat.cuda._wrapper_per_pixel import (
    blend_background,
    blend_background_noise,
    rgb_to_srgb,
    get_color_transform_matrix,
    _make_lazy_cuda_func
)
from spirulae_splat.splat.cuda import (
    _C,
)

from spirulae_splat.modules.camera import Cameras, CameraType
from spirulae_splat.splat import depth_to_normal




@dataclass
class SpirulaeSplatModelConfig:

    # Representation
    primitive: Literal["3dgs", "mip", "3dgut"] = "3dgut"
    """Splat primitive to use"""
    sh_degree: int = 3
    """Maximum degree of spherical harmonics to use."""
    sh_degree_warmup_every: int = 1000
    """Increase SH degree every this number of iterations"""
    background_mode: Literal["black", "noise", "sh"] = "black"
    """Background mode, black per convention, noise to discourage transparency, sh for skybox."""
    background_noise_warmup: int = 2000
    """Number of steps to warmup background noise. This applies when background_mode is noise"""
    background_noise_pre_warmup: float = 0.25
    """Weight of background noise at start of training (0 to 1). Higher value reduce the chance of washing away splat opacities near the beginning of training."""
    background_sh_degree: int = 4
    """SH degree for background color, only used when background_mode is sh."""

    # Training loss
    relative_scale: Optional[float] = None
    """Manually set scale when a scene is poorly scaled, i.e. increase this for large datasets.
        If not set, will use a scale agnostic optimizer. To prevent this, set it to 1.0."""
    l1_weight: float = 1.0
    """Weight of L1 loss, default 1.0"""
    l2_weight: float = 0.0
    """Weight of L2 loss, default 0.0"""
    ssim_lambda: float = 0.2
    """Weight of ssim loss; 0.2 for academic baseline, higher for potentially more high-frequency details, lower for less blurry background in outdoor scenes"""
    # YUV (BT.601) per-pixel supervision. These are raw, relative weights —
    # they sit on the same "element-wise RGB" budget as l2_lambda and the
    # implicit RGB-L1 (whose raw weight is 1.0); the loss helper normalizes
    # the full {RGB-L1, RGB-L2, Y-L1, Y-L2, U-L2, V-L2} group to sum to
    # (1 - ssim_lambda), preserving the user's RGB-vs-YUV split.
    l1_weight_y: float = 0.0
    """Weight of per-pixel BT.601 luma (Y) L1 loss."""
    l2_weight_y: float = 0.0
    """Weight of per-pixel BT.601 luma (Y) L2 loss."""
    l2_weight_u: float = 0.0
    """Weight of per-pixel BT.601 chroma U L2 loss."""
    l2_weight_v: float = 0.0
    """Weight of per-pixel BT.601 chroma V L2 loss."""
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
    use_fused_proj_bwd_optim: bool = True
    """Whether to use fused projection backward and optimizer.
        More memory efficient for large number of Gaussians, with slight performance hit."""
    sh_quantization_level: int = 1
    """SH quantization level: a single int that selects one of two
        (param bits, optim bits) configurations.
            0 = off          : 32-bit param, fp32 optim
            1 = light        : 16-bit param, 8-bit packed optim (2 B / cell)
        Collapsing the prior independent param+optim bit controls into a
        single level minimizes the FPBO kernel instantiations."""
    compute_hessian_diagonal: Literal[None, "position", "all"] = None
    """What parameter sets to compute an approximation of Hessian diagonal as well as a Jacobian-residual product in backward pass. Required for second-order optimizer."""
    optimizer_offload: Literal[None, "sh", "all"] = None
    """Whether to offload optimizer momentum to CPU to save VRAM. This is only supported for Adam optimizer."""
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
    noise_lr: float = 80.0  # 5e5 * 1.6e-4
    """MCMC-like noise injection magnitude at start of training"""
    noise_lr_final: float = 0.8  # 5e5 * 1.6e-6
    """MCMC-like noise injection magnitude at end of training"""
    min_opacity: float = 0.005
    """Minimum Gaussian opacity before relocation"""
    growth_factor: float = 1.05
    """Multiply number of splats by this number at each densification step"""
    use_revised_densification: bool = True
    """Whether to use revised densification instead of original MCMC."""
    densify_score_mode: Literal["mean", "max", "median", "geom"] = "mean"
    """How to accumulate per-splat scores across iterations for densification.
        "mean":   running mean of |w|.
        "max":    running max of |w|.
        "median": running median of |w| (approximation).
        "geom":   running geometric mean of |w|."""
    use_edge_aware_score: bool = False
    """Whether to use edge aware score to guide densification.
        If True, it computes edge aware score following https://arxiv.org/abs/2603.08661
        Note that this is only active when use_revised_densification"""
    use_loss_map: bool = True
    """Whether to use loss map to guide densification.
        Note that this is only active when use_revised_densification."""
    use_structure_only_loss_map: bool = True
    """When True, the densification loss map is filled only from the SSIM
        structure term (no per-pixel L1/L2 / auxiliary supervisory terms and
        no SSIM luminance/contrast). Intended to bias densification toward
        pattern/edge mismatches instead of brightness/contrast errors.
        Affects loss_map only; training gradients and scalar losses unchanged."""
    use_long_axis_split: bool = True
    """whether to use long-axis split described in https://arxiv.org/abs/2508.12313 for relocation and sample add.
        When combined with use_revised_densification, this can give less blurry background details for unbounded outdoor scenes."""
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
    bilagrid_optim_bits: int = 8
    """Bit depth for the RGB bilagrid optimizer state quantization. 32 = off
        (fp32 g1/g2 for Adam, fp32 accum for AdaGrad). 4 or 8 = packed
        QuantizedAdamState / QuantizedTensorLog at BITS."""
    bilagrid_geometry_optim_bits: int = 8
    """Bit depth for the depth+normal bilagrid optimizer state quantization.
        Same semantics as bilagrid_optim_bits."""
    use_adagrad_bilagrid_optim: bool = False
    """Use AdaGrad (lr_decay=0, weight_decay=0, initial_accumulator_value=0,
       eps=1e-15) instead of Adam for all bilateral-grid parameters (RGB +
       depth + normal). When True, the bilagrid LR fields read from
       OptimizerConfig switch to ``bilagrid_adagrad_*_lr``. If
       ``bilagrid_*_optim_bits`` is 4 or 8 the per-cell AdaGrad accumulator
       is stored block-wise log-quantized at that bit depth
       (QuantizedTensorLog<BITS, 256>)."""
    bilagrid_tv_loss_weight: float = 10.0
    """Total variation loss weight for bilateral grid used for radiance"""
    bilagrid_shift_reg_weight: float = 0.0
    """RGB bilagrid color-shift regularizer. Penalizes the
        dataset-wide mean of sign(post-bilagrid - pre-bilagrid) per channel:
        R = w * ||EMA[mean_p sign(post - pre)]||^2. The gradient is injected
        on the POST-bilagrid v_render_rgb buffer and flows through the
        bilagrid vjp into the grid params (and as a small leak into the
        splats). 0 disables. Typical values: 0.01--1.0; tune w.r.t. main loss
        scale."""
    bilagrid_shift_reg_ema_period: int = 750
    """EMA decay length for the bilagrid color-shift regularizer, in BATCHES.
        beta = max(0, 1 - 1/period). Should be roughly the number of batches
        per epoch so the EMA estimates the dataset-wide mean. Ignored when
        bilagrid_shift_reg_weight = 0."""
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
    scale_init: Optional[float] = None
    """Initial scale. If not set, auto decide"""
    opacity_init: Optional[float] = None
    """Initial opacity. If not set, auto decide"""
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
    alpha_loss_weight_under: float = 0.0
    """Loss weight for alpha, applies when rendered alpha is below reference alpha"""
    opacity_reg: float = 0.01  # 0.01 for MCMC
    """Encourage low opacity to aid densification, per MCMC."""
    scale_reg: float = 0.01  # 0.01 for MCMC
    """Encourage low scale, per MCMC."""
    opacity_decay: float = 0.0  # 0.004 in MRNF
    """Decay opacity to aid densification, per MRNF."""
    scale_decay: float = 0.0  # 0.002 in MRNF
    """Decay scale to aid densification, per MRNF."""
    erank_reg: float = 0.0
    """erank regularization weight, for 3DGS only -
        see *Effective Rank Analysis and Regularization for Enhanced 3D Gaussian Splatting, Hyung et al.*"""
    erank_reg_s3: float = 0.0
    """erank regularization weight for smallest dimension, for 3DGS only"""
    quat_norm_reg: float = 0.01
    """Weight to regularize quaternion norm to identity"""
    sh_reg: float = 0.001
    """Regularize SH magnitude to find a balance between bilagrid/PPISP and improve generalizability."""
    overexposure_reg: float = 0.0
    """Image-space L2 penalty on rendered RGB outside [0, 1]:
       L = w * mean(max(-x, x-1, 0)^2) over all pixels and channels of the
       raw rendered image (pre-bilagrid / pre-PPISP / pre-color-space).
       When non-zero a dedicated CUDA kernel adds dL/dx directly into the
       rendered-RGB gradient; the scalar loss value is never materialized."""

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

        # Pre-init guess at the post-split camera count for the bilagrid table.
        # When use_cpp_data_manager is on, Trainer._setup_cpp_data_manager
        # overrides this with the exact per-camera-K sum, so mixed datasets
        # work even when this guess is approximate.
        if trainer_config.datamanager.use_cpp_data_manager:
            # Worst-case upper bound: every camera split 6x. The trainer
            # downscales to the real per-camera sum right after model init.
            self.num_train_data *= 6
        elif CameraType.EQUIRECTANGULAR.value in cameras.camera_type:
            assert len(set(cameras.camera_type)) == 1, "Mixed equirectangular and pinhole/fisheye is not supported"
            self.num_train_data *= 6  # TODO
        elif trainer_config.datamanager.warp_to_pinhole:
            self.num_train_data *= 5  # TODO

        self.info = {}

        self.g1_bilagrid = None
        self.g2_bilagrid = None
        self.g1_ppisp = None
        self.g2_ppisp = None

        self.populate_modules()

        if self.config.primitive in ['3dgs', 'mip', '3dgut']:
            splat_params = (
                self.means, self.quats, self.scales,
                self.opacities,
                self.features_dc, self.features_sh
            )

        self.core = Renderer(
            self.config.primitive,
            [x.cpu().contiguous() for x in splat_params],
            self.seed_points[0].shape[0],
            packed=(self.config.packed or self.config.use_bvh),
            use_bvh=self.config.use_bvh,
            use_fused_proj_bwd_optim=self.config.use_fused_proj_bwd_optim,
            sh_quantization_level=self.config.sh_quantization_level,
        )

        # Engine-side background blending. populate_modules already created the
        # background_color / background_sh Parameters; mode is chosen from the
        # config and matched in the C++ engine state. Without one of these, the
        # engine_train_step path skips blending entirely and the model never
        # learns to fill transparent regions (legacy `blend_background` in
        # model.py runs only for eval/viewer rendering, not training).
        #
        # `_engine_bg_mode` is consulted by self.blend_background() to fetch
        # the rendered skybox / mean-color image from the C++ side (otherwise
        # the dropdown would render a stale Python image from the un-trained
        # self.background_color / background_sh Parameters).
        self._engine_bg_mode = "none"
        if self.config.background_mode == "noise":
            self.core.engine_init_background_noise(
                bool(self.config.splat_color_is_linear))
            self._engine_bg_mode = "noise"
        elif self.config.background_mode == "sh":
            self.core.engine_init_background_sh(
                int(self.config.background_sh_degree),
                splat_color_is_linear=bool(self.config.splat_color_is_linear)
            )
            self._engine_bg_mode = "sh"

        # Engine-side linear / wide-gamut color space. The splat side (pred
        # RGB) runs every forward and inverts through the vjp on the loss
        # gradient before raster bwd. The image side (GT RGB) runs once per
        # upload inside set_training_data. Either side disabled -> identity.
        splat_cs_on  = bool(self.config.splat_color_is_linear) or (self.config.splat_color_gamut is not None)
        image_cs_on  = bool(self.config.image_color_is_linear) or (self.config.image_color_gamut is not None)
        splat_matrix = get_color_transform_matrix(self.config.splat_color_gamut) if splat_cs_on else None
        image_matrix = get_color_transform_matrix(self.config.image_color_gamut) if image_cs_on else None
        self.core.engine_init_color_space(
            splat_enabled=splat_cs_on,
            splat_is_linear=bool(self.config.splat_color_is_linear),
            splat_color_matrix=splat_matrix,
            image_enabled=image_cs_on,
            image_is_linear=bool(self.config.image_color_is_linear),
            image_color_matrix=image_matrix,
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

        scale_init, opacity_init = self.config.scale_init, self.config.opacity_init
        if scale_init is None:
            scale_init = 0.5
        if opacity_init is None:
            opacity_init = 0.1

        if self.config.primitive in ["3dgs", "mip", "3dgut"]:
            from time import perf_counter
            time0 = perf_counter()
            distances, indices = self.k_nearest_neighbor(means.data, 4)
            time1 = perf_counter()
            # print("k_nearest_neighbor time:", time1-time0)
            num_points = means.shape[0]
            scales = torch.sqrt(torch.from_numpy(distances**2).mean(-1, keepdim=True)).repeat(1, 3).float()
            scales = torch.log(scale_init * scales + 1e-8)
            if self.config.suppress_initial_scales:
                scales = self.suppress_initial_scales(means, scales)
            quats = F.normalize(torch.randn((num_points, 4)), dim=-1)
            opacities = torch.logit(opacity_init * torch.ones(num_points, 1))

        gauss_params = {}
        if self.config.primitive in ["3dgs", "mip", "3dgut"]:
            gauss_params['means'] = torch.nn.Parameter(means)
            gauss_params['quats'] = torch.nn.Parameter(quats)
            gauss_params['scales'] = torch.nn.Parameter(scales)
            gauss_params['opacities'] = torch.nn.Parameter(opacities)

        # colors
        dim_sh = (self.config.sh_degree + 1) ** 2

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

        if self.config.primitive in ["3dgs", "mip", "3dgut"]:
            shs = torch.zeros((num_points, dim_sh, 3)).float()
            shs[:, 0, :3] = (seed_color - 0.5) / 0.28209479177387814
            if self.config.sh_degree > 0:
                shs[:, 1:, 3:] = 0.0
            gauss_params['features_dc'] = torch.nn.Parameter(shs[:, 0, :].contiguous())
            gauss_params['features_sh'] = torch.nn.Parameter(shs[:, 1:, :].contiguous())

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

        # self.gauss_params = torch.nn.ParameterDict(gauss_params)
        self.gauss_params = gauss_params

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

        self._bilagrid_rgb_init = False
        self._bilagrid_depth_init = False
        self._bilagrid_normal_init = False
        # PPISP: also lazy. Allocated on first iteration where config + LR is
        # positive (RGB always has a supervision path via rendering loss).
        self._ppisp_init = False
        self.training_verboser = TrainingVerbose()

        self.step = 0

        self._train_batch_size = 1

        torch.cuda.empty_cache()

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
    def opacities(self):
        param = self._get_gauss_param("opacities")
        return param

    def load_state_dict(self, dict, **kwargs):  # type: ignore
        # resize the parameters to match the new number of points
        self.step = self.trainer_config.num_iterations
        for p in ["means", "scales", "quats",
                    "features_dc", "features_sh",
                    "opacities"]:
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
        self.core.zero_grad()

    def get_2dgs_reg_weights(self):
        """Returns (depth_reg, normal_dist_reg, rgb_dist_reg) and normal_reg,
        each scaled by distortion_reg_warmup progress. Lifted from training_losses.py."""
        factor = min(self.step / max(self.config.distortion_reg_warmup, 1), 1)
        weight_depth_reg = self.config.depth_distortion_reg * factor
        weight_normal_dist_reg = self.config.normal_distortion_reg * factor
        weight_rgb_dist_reg = self.config.rgb_distortion_reg * factor
        weight_normal_reg = self.config.normal_reg_weight * factor
        return (weight_depth_reg, weight_normal_dist_reg, weight_rgb_dist_reg), weight_normal_reg

    def _maybe_init_bilagrid(self, batch):
        """Lazily allocate each bilagrid type the first time all of its
        enablement conditions are met. Called per training iteration.

        A bilagrid is allocated only when:
          - the corresponding config flag is True;
          - the base learning rate (config) is positive — a 0 LR means the
            grid would never update, so allocation is pure waste;
          - for depth / normal: the relevant supervision weight is positive AND
            the current batch actually contains that modality (gt depth /
            normal). RGB always has a gradient path through the rendering loss,
            so it doesn't gate on a supervision weight.

        Once a type is allocated it stays for the rest of training; LR / TV
        scheduling is handled per-step in engine_train_step.
        """
        cfg = self.config
        optim_cfg = self.trainer_config.optimizer

        # RGB: config flag + positive base LR. Always has supervision via the
        # rendering loss, and RGB is always present in training batches.
        if (not self._bilagrid_rgb_init
                and cfg.use_bilateral_grid
                and optim_cfg.bilagrid_lr > 0.0):
            X, Y, W_g = cfg.bilagrid_shape  # (grid_X, grid_Y, grid_W) → (W, H, L)
            self.core.engine_init_bilagrid(
                self.num_train_data,
                rgb_type=cfg.bilagrid_type,
                rgb_LHW=(W_g, Y, X),
                rgb_optim_bits=cfg.bilagrid_optim_bits,
                use_adagrad=cfg.use_adagrad_bilagrid_optim,
            )
            self._bilagrid_rgb_init = True

        # Depth: config flag + LR + supervision weight + dataset has depth.
        if (not self._bilagrid_depth_init
                and cfg.use_bilateral_grid_for_geometry
                and optim_cfg.bilagrid_depth_lr > 0.0
                and cfg.depth_supervision_weight > 0.0
                and batch.get('depth', None) is not None):
            X, Y, W_g = cfg.bilagrid_shape_geometry
            self.core.engine_init_bilagrid(
                self.num_train_data,
                depth_LHW=(W_g, Y, X),
                geometry_optim_bits=cfg.bilagrid_geometry_optim_bits,
                use_adagrad=cfg.use_adagrad_bilagrid_optim,
            )
            self._bilagrid_depth_init = True

        # Normal: config flag + LR + supervision weight + dataset has normals.
        if (not self._bilagrid_normal_init
                and cfg.use_bilateral_grid_for_geometry
                and optim_cfg.bilagrid_normal_lr > 0.0
                and cfg.normal_supervision_weight > 0.0
                and batch.get('normal', None) is not None):
            X, Y, W_g = cfg.bilagrid_shape_geometry
            self.core.engine_init_bilagrid(
                self.num_train_data,
                normal_LHW=(W_g, Y, X),
                geometry_optim_bits=cfg.bilagrid_geometry_optim_bits,
                use_adagrad=cfg.use_adagrad_bilagrid_optim,
            )
            self._bilagrid_normal_init = True

        # PPISP (RGB only): config flag + positive base LR. Always has supervision
        # via the rendering loss; no dataset-side gating needed.
        if (not self._ppisp_init
                and cfg.use_ppisp
                and optim_cfg.ppisp_lr > 0.0):
            self.core.engine_init_ppisp(
                self.num_train_data,
                param_type=cfg.ppisp_param_type
            )
            self._ppisp_init = True

    def step_post_backward(self):
        return  # TODO

    def get_gaussian_param_groups(self) -> Dict[str, List[Parameter]]:
        # Here we explicitly use the means, scales as parameters so that the user can override this function and
        # specify more if they want to add more optimizable params to gaussians.
        return {
            name: [self.gauss_params[name]]
            for name in [
                "means", "scales", "quats",
                "features_dc", "features_sh",
                "opacities"
            ] if name in self.gauss_params
        }

    def get_param_groups(self) -> Dict[str, List[Parameter]]:
        """Obtain the parameter groups for the optimizers

        Returns:
            Mapping of different parameter groups
        """
        raise NotImplementedError()

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

    def blend_background(
            self, camera: Cameras, c2w: torch.Tensor, intrins: torch.Tensor,
            W: Optional[int] = None, H: Optional[int] = None, camera_model: str = None,
            rgb: Optional[torch.Tensor] = None, transmittance: Optional[torch.Tensor] = None
        ):
        if not isinstance(camera, Cameras):
            print("Called blend_background with not a camera")
            return {}

        if W is None or H is None:
            W, H = int(camera.width[0].item()), int(camera.height[0].item())
        if camera_model is None:
            camera_model = camera.camera_type[0]
        sh_degree = self.config.background_sh_degree

        if self.config.background_mode not in ['noise', 'sh']:
            return None

        if getattr(self, "_engine_bg_mode", "none") != "none":
            if rgb is not None and transmittance is not None:
                return rgb
            B = int(transmittance.shape[0]) if transmittance is not None else len(camera)
            bg_buf = torch.empty(B, H, W, 3, dtype=torch.float32)
            if self.core.engine_copy_background_to_host(bg_buf):
                return bg_buf
            return None

        raise NotImplementedError()
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
                camera_model.upper(),
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

    def get_outputs(self, camera: Cameras, val: bool=False, want_keys=None) -> Dict[str, Union[torch.Tensor, List]]:
        """Takes in a camera and returns a dictionary of outputs.

        `want_keys` is an optional iterable of keys the caller is interested in.
        Currently used to skip the debug SH/refinement_score renders and the
        depth_normal kernel in the eval/viewer path when they weren't requested.
        Passing None (default) keeps the legacy "render everything" behavior.
        """
        _want = set(want_keys) if want_keys is not None else None

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
        camera_model = camera.camera_type[0].upper()

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
        # if 'patch_offsets' in camera.metadata:
        #     kwargs['patch_offsets'] = camera.metadata['patch_offsets']
        if not self.training:
            pass
            # camera_model = "FISHEYE"
            # kwargs['dist_coeffs'] = torch.tensor([[0.0259, 0.0082, 0.0002, -0.0013, -0.0012, -0.0008, 0.0000, 0.0000, -0.0006, -0.0001]]).float().cuda()

            # TODO: investigate why this uses a ton of VRAM
            # kwargs['dist_coeffs'] = torch.tensor([[-0.29, 0.07, 0, 0, 1e-5, -1e-3, 0, 0, 0, 0]]).float().cuda()

        if camera.distortion_params is not None and camera.distortion_params.any():
            kwargs['dist_coeffs'] = camera.distortion_params

        if 'accum_weight_map' in camera.metadata:
            accum_weight_map = camera.metadata['accum_weight_map']
            if accum_weight_map is not None:
                kwargs['accum_weight_map'] = accum_weight_map

        # Engine path: keep on CPU, C++ does H→D copy
        optimized_camera_to_world = optimized_camera_to_world.cpu()
        viewmats = viewmats.cpu().contiguous()
        intrins = intrins.cpu().contiguous()

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
        self.core.set_params(
            viewmats=viewmats,  # [C, 4, 4]
            intrins=intrins,  # [C, 4]
            width=W,
            height=H,
            sh_degree_to_use=self.step // max(self.config.sh_degree_warmup_every, 1),
            # packed=True,
            # use_bvh=True,
            relative_scale=self.config.relative_scale,
            camera_model=camera_model,
            output_distortion=any([c != 0.0 for c in self.get_2dgs_reg_weights()[0]]),
            compute_hessian_diagonal=self.config.compute_hessian_diagonal,
            **kwargs,
        )
        self.core.forward()
        rgbdn = self.core.render_colors
        Ts = self.core.render_Ts
        meta = self.core.meta
        # torch.cuda.empty_cache()

        # if not self.training:
        #     W = gw*TILE_SIZE
        #     H = gh*TILE_SIZE
        #     # print(rgbdn.shape, Ts.shape)
        #     rgbdn = [merge_tiles(comp) for comp in rgbdn]
        #     Ts = merge_tiles(Ts)
        #     for key in ['rgb_distortion', 'depth_distortion', 'normal_distortion']:
        #         if key in meta:
        #             meta[key] = merge_tiles(meta[key])
        #     viewmats, intrins = viewmats_0, intrins_0
        #     if 'dist_coeffs' in kwargs:
        #         kwargs['dist_coeffs'] = kwargs['dist_coeffs'][:1]

        if len(rgbdn) == 2:
            rgb, depth_im_ref = rgbdn
            render_normal = None
        else:
            rgb, depth_im_ref, render_normal = rgbdn

        # normals
        if render_normal is not None:  # TODO: fused kernel
            render_normal = torch.where(Ts < 1.0, F.normalize(render_normal, dim=-1), render_normal)

        depth_normal = None
        if (depth_im_ref is not None and not self.training
                and (_want is None or 'depth_normal' in _want)):
            # depth_im_ref is CPU (from engine_copy_render_to_host). Move to CUDA for the kernel.
            depth_cuda = depth_im_ref.cuda().contiguous() if not depth_im_ref.is_cuda else depth_im_ref.contiguous()
            intrins_cuda = intrins.cuda().contiguous() if not intrins.is_cuda else intrins.contiguous()
            if camera.distortion_params is not None:
                dist_coeffs_cuda = camera.distortion_params.cuda().contiguous()
            else:
                dist_coeffs_cuda = torch.zeros(len(camera), 10, dtype=torch.float32, device='cuda')
            depth_normal_cuda = torch.empty(
                *depth_cuda.shape[:-1], 3, dtype=torch.float32, device='cuda')
            _C.depth_to_normal_forward(
                camera.camera_type[0].upper(),
                self.core._tv(intrins_cuda),
                self.core._tv(dist_coeffs_cuda),
                self.config.primitive not in ['3dgs', 'mip'],  # is_ray_depth
                self.core._tv(depth_cuda),
                self.core._tv(depth_normal_cuda),
            )
            depth_normal = depth_normal_cuda.cpu()

        # radii = meta["radii"]
        # depths = meta["depths"]

        if self.training:
            # self.info['radii'] = radii

            self.info.update({
                "width": kwargs["actual_width"],
                "height": kwargs["actual_height"],
                "n_cameras": len(camera),
                "n_train": self.num_train_data,
                "backward_info": self.core.backward_info,
            })
            if 'patch_offsets' in camera.metadata:
                self.info['patch_offsets'] = camera.metadata['patch_offsets']
            if 'actual_images_per_batch' in camera.metadata:
                self.info['n_cameras'] = camera.metadata['actual_images_per_batch']
            for key in ['gaussian_ids', 'camera_ids', 'bvh_time']:
                if key in meta:
                    self.info[key] = meta[key]

        # pack outputs
        outputs = {
            "rgb": rgb,
        }
        if self.training:
            outputs["backward_info"] = self.core.backward_info
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
                    if not self.training:
                        value = torch.sqrt(value + (1/255)**2) - (1/255)
                    outputs[key] = value

        if not self.training:
            outputs["alpha"] = (1.0 - Ts).reshape((H, W, 1)).repeat(1, 1, 3)
            background = self.blend_background(camera, optimized_camera_to_world, intrins, W, H, camera_model)
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
            # Color-space conversion (linear / wide-gamut -> sRGB) for the
            # rendered RGB + skybox is done entirely on the C++ side:
            # `outputs["rgb"]` comes from engine_copy_render_to_host which
            # writes the sRGB post-conversion render; the pre-conversion
            # render is in self.core.render_rgb_raw. `outputs["background"]`
            # comes from engine_copy_background_to_host which also applies
            # the conversion before D->H. We just clip & wire here.
            if self.config.splat_color_is_linear or self.config.splat_color_gamut is not None:
                outputs["rgb_raw"] = self.core.render_rgb_raw
                outputs["rgb"] = outputs["rgb"].clip(0, 1)
                if "background" in outputs:
                    outputs["background"] = outputs["background"].clip(0, 1)
            for key in outputs:
                outputs[key] = outputs[key].squeeze(0)

        if self.training:
            kwargs["intrins"] = intrins
            kwargs["camera_model"] = camera_model
            outputs["camera"] = camera
            outputs["camera_intrins"] = kwargs
        # if not self.training and True:
        #     outputs['ray'] = merge_tiles(meta['intersection_count'].float().reshape(-1, 1, 1, 1).repeat(1, TILE_SIZE, TILE_SIZE, 1))

        # setup override required for optimizers
        # TODO: more reliable way than setattr tensor?
        if self.training and False:
            optim_info = {
                # 'radii': self.info['radii'],
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
        # if not hasattr(self.core, 'densify_accum_buffer'):
        #     return outputs

        # Debug SH: render with DC zeroed out to visualize SH-only contribution.
        # Skipped on the viewer hot path when the caller didn't request the
        # `sh` buffer (want_keys gating).
        if not self.training and (_want is None or 'sh' in _want):
            zero_dc = torch.zeros_like(self.features_dc)
            sh_rgb = self.core.engine_debug_forward(
                override_features_dc=zero_dc,
                override_sh_degree=self.step // max(self.config.sh_degree_warmup_every, 1)
            )
            outputs['sh'] = sh_rgb[0]  # first camera

        # Debug densification: visualize per-splat accum weight as color
        if not self.training and self.step > 1 and (_want is None or 'refinement_score' in _want):
            accum_weight = self.core.engine_get_accum_weight()
            param_to_vis = accum_weight
            mean_val = param_to_vis.mean()
            if mean_val > 0:
                param_to_vis = param_to_vis / mean_val
            param_to_vis = (param_to_vis - 0.5) / 0.28
            param_to_vis = param_to_vis.unsqueeze(-1).repeat(1, 3)
            # Pad to max_num_splats
            if len(param_to_vis) < len(self.means):
                param_to_vis = torch.cat((param_to_vis, torch.zeros(
                    len(self.means) - len(param_to_vis), 3, device=self.means.device)), 0)
            vis_rgb = self.core.engine_debug_forward(
                override_features_dc=param_to_vis.contiguous(),
                override_sh_degree=0
            )
            outputs['refinement_score'] = vis_rgb[0].mean(dim=-1, keepdim=True)

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

        render_grads = (v_render_rgb, v_render_depth, v_render_normal)

        self.core.backward(
            render_grads,
            v_render_Ts,
            v_rgb_distortion,
            v_depth_distortion,
            v_normal_distortion,
        )

    def optim_step(self):
        max_steps = self.trainer_config.num_iterations
        optim_config = self.trainer_config.optimizer

        # Engine path: optimizer + densification via core -> Engine.cpp
        # Splats stay on device — no D→H copy between iterations.
        if self.core.primitive in ['3dgs', 'mip', '3dgut'] and not self.core.use_bvh:
            self.core.engine_optim_step(self.step, max_steps, self.config, optim_config)
            self.core.engine_densify_step(self.step, max_steps, self.config)
            return
        raise NotImplementedError()

        # Legacy Python path
        self.core.optim_step(self.config, optim_config, self.step, max_steps)
        self.training_losses.optim_step(optim_config, self.step, max_steps)
        self.core.densify_step(self.step, max_steps, self.config, None)
        # self.step_post_backward()

    def verbose(self):

        if self.step % self.training_verboser.cache_skip == 0:

            # Per splat losses
            # hyperparams = (
            #     np.sign(self.config.opacity_reg) ** 0,
            #     np.sign(self.config.scale_reg) ** 0,
            #     np.sign(self.config.max_gauss_ratio),
            #     np.sign(self.config.scale_regularization_weight),
            #     np.sign(self.config.erank_reg),
            #     np.sign(self.config.erank_reg_s3),
            #     np.sign(self.config.quat_norm_reg)
            # )
            # losses = _make_lazy_cuda_func("compute_per_splat_losses_forward")(
            #     self.scales, self.opacities, self.quats,
            #     *hyperparams
            # )
            # (opac, scale, aniso, erank, quatnorm) = losses.cpu().numpy().tolist()
            # self.training_verboser.add_metric('opac', opac, last_only=True)
            # self.training_verboser.add_metric('scale', scale, last_only=True)
            # if aniso > 0:
            #     self.training_verboser.add_metric('aniso', aniso, last_only=True)
            # if erank > 0:
            #     self.training_verboser.add_metric('erank', erank, last_only=True)

            # Static losses removed; bilagrid TV losses now come from the C++
            # engine_train_step loss_dict each iteration (see verbose metrics).

            # VRAM breakdown (C++ pool only — PyTorch-side tensors are tracked
            # separately by torch.cuda.memory_allocated). Each pool key gets
            # mapped to one of five buckets by its prefix:
            #   splat       : world params + grads + optim states (per-Gaussian)
            #   splat x img : projection outputs/gradients, tile intersection
            #   image       : image-space tensors (renders, GT copies, grads)
            #   appearance  : appearance + geometry grids and ppisp merged
            #   other       : everything else (camera tables, tiny scratches)
            from spirulae_splat.splat.cuda import _C
            buckets = {'splat': 0.0, 'splat x img': 0.0, 'image': 0.0,
                       'appearance': 0.0, 'viewer': 0.0, 'other': 0.0}
            GiB = 1024 ** 3
            for key, _used, cap in _C.engine_get_pool_breakdown():
                gib = cap / GiB
                if key.startswith('world.') \
                        or key.startswith('eng.v_') \
                        or key.startswith('eng.g1_') \
                        or key.startswith('eng.g2_') \
                        or key.startswith('eng.sh_quant_') \
                        or key in ('eng.radii', 'eng.accum_buffer',
                                   'eng.bias_correction_steps',
                                   'eng.quant_bounds_sh'):
                    buckets['splat'] += gib
                elif key.startswith('proj.') \
                        or key.startswith('isect.') \
                        or key.startswith('isect_post.') \
                        or key.startswith('fused_proj_bwd.') \
                        or key.startswith('raster_bwd.'):
                    buckets['splat x img'] += gib
                elif key.startswith('render.') or key.startswith('renders.') \
                        or key in ('eng.v_rgb', 'eng.v_depth', 'eng.v_Ts',
                                   'eng.v_depth_normal', 'eng.v_ref_depth',
                                   'eng.v_ref_normal',
                                   'eng.depth_normal', 'eng.loss_map',
                                   'eng.bg_sky.v_Ts_scratch',
                                   'gt.rgb', 'gt.normal', 'gt.staging_u8'):
                    buckets['image'] += gib
                elif key.startswith('eng.bg.') or key.startswith('eng.ppisp.'):
                    buckets['appearance'] += gib
                elif key.startswith('vis.') or key.startswith('viewer.'):
                    buckets['viewer'] += gib
                else:
                    buckets['other'] += gib
            # DeviceScratch (workspace for cub::DeviceRadixSort etc.) is part
            # of the splat x image pipeline (used during tile intersection
            # sort + projection backward sort).
            buckets['splat x img'] += _C.engine_get_scratch_bytes() / GiB

            self.training_verboser.add_metric('splat_vram', buckets['splat'], last_only=True)
            self.training_verboser.add_metric('image_vram', buckets['image'], last_only=True)
            self.training_verboser.add_metric('splat_x_image_vram', buckets['splat x img'], last_only=True)
            self.training_verboser.add_metric('appearance_vram', buckets['appearance'], last_only=True)
            self.training_verboser.add_metric('viewer_vram', buckets['viewer'], last_only=True)
            self.training_verboser.add_metric('other_vram', buckets['other'], last_only=True)

        self.training_verboser.verbose(self.step, self.trainer_config.num_iterations)

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

    @staticmethod
    def _tv(tensor: torch.Tensor):
        if tensor is None:
            return (0, 4, [0])
        assert tensor.is_contiguous(), f"Tensor must be contiguous, got strides {tensor.stride()}"
        return (tensor.data_ptr(), tensor.element_size(), list(tensor.shape))

    def _build_loss_weights(self, step):
        """Build the per-pixel loss_weights array passed to the engine.

        The "element-wise RGB" group {RGB-L1, RGB-L2, Y-L1, Y-L2, U-L2, V-L2}
        is normalized so its total equals ``1 - ssim_lambda``, and the split
        between the RGB-only sub-group and the YUV sub-group preserves the
        user's raw config ratio. RGB-internal L1/L2 split (l2_lambda) and the
        YUV-internal ratios are also preserved.

        With all four YUV weights at 0 (default), this reduces exactly to the
        old (1 - l2_lambda) * (1 - ssim_lambda), l2_lambda * (1 - ssim_lambda)
        split.

        Order must match per_pixel_losses.slang::LossWeightIndex.
        """
        cfg = self.config
        dist_factor = min(step / max(cfg.distortion_reg_warmup, 1), 1.0)
        reg_active = float(step >= cfg.reg_warmup_length)
        sup_active = float(step > cfg.supervision_warmup)
        alpha_reg_factor = cfg.alpha_reg_weight * min(
            step / max(cfg.alpha_reg_warmup, 1), 1.0)

        # Raw element-wise RGB-group weights (relative). RGB-L1 fills the slot
        # left by L2 inside the RGB sub-budget of 1.0; YUV weights are
        # configured directly by the user.
        w_rgb_l1_raw = max(0.0, cfg.l1_weight)
        w_rgb_l2_raw = max(0.0, cfg.l2_weight)
        w_y_l1_raw   = max(0.0, cfg.l1_weight_y)
        w_y_l2_raw   = max(0.0, cfg.l2_weight_y)
        w_u_l2_raw   = max(0.0, cfg.l2_weight_u)
        w_v_l2_raw   = max(0.0, cfg.l2_weight_v)
        total_raw = (w_rgb_l1_raw + w_rgb_l2_raw
                     + w_y_l1_raw + w_y_l2_raw + w_u_l2_raw + w_v_l2_raw)
        # When all weights are 0 (degenerate), keep zeros to avoid NaN.
        scale = ((1.0 - cfg.ssim_lambda) / total_raw) if total_raw > 0.0 else 0.0

        return [
            w_rgb_l1_raw * scale,                              # RgbSupL1
            w_rgb_l2_raw * scale,                              # RgbSupL2
            w_y_l1_raw   * scale,                              # YSupL1
            w_y_l2_raw   * scale,                              # YSupL2
            w_u_l2_raw   * scale,                              # USupL2
            w_v_l2_raw   * scale,                              # VSupL2
            sup_active * cfg.depth_supervision_weight,         # DepthSup
            sup_active * cfg.normal_supervision_weight,        # NormalSup
            cfg.apply_loss_for_mask * cfg.alpha_loss_weight,   # AlphaSup
            cfg.apply_loss_for_mask * cfg.alpha_loss_weight_under,  # AlphaSupUnder
            reg_active * cfg.normal_reg_weight * dist_factor,  # NormalReg
            reg_active * alpha_reg_factor,                     # AlphaReg
            reg_active * cfg.rgb_distortion_reg * dist_factor,  # RgbDistReg
            reg_active * cfg.depth_distortion_reg * dist_factor,  # DepthDistReg
            reg_active * cfg.normal_distortion_reg * dist_factor,  # NormalDistReg
        ]

    def _engine_get_loss_grad(self, outputs, batch, batch_size):
        """Engine path: compute loss + backward entirely in C++."""

        # --- Prepare ground truth (CPU, C++ does H→D) ---
        gt_rgb = batch["image"].cpu()
        if gt_rgb.dtype == torch.uint8:
            gt_rgb = gt_rgb.float() / 255.0
        elif gt_rgb.dtype == torch.uint16:
            gt_rgb = gt_rgb.float() / 65535.0
        elif gt_rgb.dtype == torch.float16:
            gt_rgb = gt_rgb.float()
        gt_rgb = gt_rgb[..., :3].contiguous()

        gt_depth = batch.get('depth', None)
        if gt_depth is not None:
            gt_depth = gt_depth.cpu().float()
            if len(gt_depth.shape) == 3:
                gt_depth = gt_depth.unsqueeze(-1)
            gt_depth = gt_depth.contiguous()

        gt_normal = batch.get('normal', None)
        if gt_normal is not None:
            if gt_normal.dtype == torch.uint8:
                gt_normal = gt_normal.float() / (255/2) - 1.0
            gt_normal = gt_normal.cpu().float().contiguous()

        # External mask -> gt_alpha. The slang kernel derives RGB / depth /
        # normal / alpha-sup masks internally from this + the gt sentinels.
        gt_alpha = batch.get('mask', None)
        if gt_alpha is not None:
            if gt_alpha.dtype != torch.bool:
                gt_alpha = gt_alpha.cpu().float() > 0.5
            gt_alpha = gt_alpha.contiguous()

        self.core.engine_set_training_data(gt_rgb, gt_depth, gt_normal, gt_alpha)

        # --- Loss weights ---
        step = self.step
        cfg = self.config

        loss_weights = self._build_loss_weights(step)
        w_ssim = cfg.ssim_lambda
        num_loss_scales = cfg.num_loss_scales + 1
        compute_loss_map = cfg.use_loss_map or (cfg.compute_hessian_diagonal is not None)
        structure_only_loss_map = cfg.use_structure_only_loss_map

        # --- Compute loss + backward via core (gradients managed by C++ pool) ---
        loss_dict = self.core.engine_compute_loss_backward(
            step, loss_weights, w_ssim, num_loss_scales, compute_loss_map,
            structure_only_loss_map,
            overexposure_reg_weight=cfg.overexposure_reg,
        )

        for key, value in loss_dict.items():
            self.training_verboser.add_metric(key, value)
        self.training_verboser.add_metric("num_splats", self.core.cur_num_splats, last_only=True)
        self.training_verboser.add_metric("max_num_splats", self.core.max_num_splats, last_only=True)
        self.training_verboser.add_metric("num_sh", min(self.core.sh_degree_to_use, self.config.sh_degree), last_only=True)
        self.training_verboser.add_metric("max_num_sh", self.config.sh_degree, last_only=True)

        self._engine_backward_done = True
        return None

    def engine_train_step(self, camera, batch):
        """Fused training step. Tensors passed as-is (CPU or CUDA) — C++ side detects
        device location and avoids redundant copies. Minimizes H↔D transfers."""
        # ---- profiling probes (gated by PROFILE_TRAIN_STEP) ----
        if PROFILE_TRAIN_STEP:
            _prof = getattr(self, "_ets_prof", None)
            if _prof is None:
                self._ets_prof = _prof = {
                    "n": 0, "warmup": 10,
                    "cam_prep": 0, "gt_prep": 0, "weights_prep": 0,
                    "presync": 0, "ccall": 0, "ccall_after_sync": 0, "post": 0,
                }
            from time import perf_counter_ns as _t
            _t0 = _t()
        # --- Camera params: keep on whatever device they came from ---
        if self.config.use_camera_optimizer:
            camera.metadata['cam_idx'] = camera.metadata['cam_idx'].flatten()
            optimized_camera_to_world = self.camera_optimizer.apply_to_camera(camera)
        else:
            optimized_camera_to_world = camera.camera_to_worlds

        camera_downscale = self._get_downscale_factor()
        if camera_downscale > 1:
            camera.rescale(1 / camera_downscale)

        W, H = int(camera.width[0].item()), int(camera.height[0].item())
        camera_model = camera.camera_type[0].upper()

        R = optimized_camera_to_world[:, :3, :3]
        T = optimized_camera_to_world[:, :3, 3:4]
        if self.config.relative_scale is not None:
            T = T * self.config.relative_scale
        R = R * torch.tensor([[[1.0, -1.0, -1.0]]]).to(R)
        R_inv = R.transpose(-1, -2)
        T_inv = -torch.bmm(R_inv, T)
        viewmats = torch.eye(4, dtype=optimized_camera_to_world.dtype,
                              device=optimized_camera_to_world.device)[None].repeat(len(camera), 1, 1)
        viewmats[:, :3, :3] = R_inv
        viewmats[:, :3, 3:4] = T_inv
        viewmats = viewmats.contiguous()
        intrins = camera.intrins.contiguous()

        dist_coeffs = None
        if camera.distortion_params is not None and camera.distortion_params.any():
            dist_coeffs = camera.distortion_params.contiguous()
        if dist_coeffs is None:
            dist_coeffs = torch.zeros(len(camera), 10, dtype=torch.float32, device=viewmats.device)

        if PROFILE_TRAIN_STEP: _t1 = _t()
        # --- GT data: hand raw bytes to the C++ engine ---
        # For uint8 / uint16 inputs, the C++ side does the float conversion on
        # GPU (3-4x smaller H->D payload). The 4 per-pixel mask buffers
        # (gt_rgb_mask, gt_depth_mask, gt_normal_mask, gt_alpha_mask) are gone:
        # they're derived inside per_pixel_losses.slang from the gt_depth = 0
        # sentinel, the gt_normal sum > -2.366 sentinel, and gt_alpha presence.
        gt_rgb = batch["image"]
        if gt_rgb.dtype == torch.float16:
            gt_rgb = gt_rgb.float()
        gt_rgb = gt_rgb[..., :3].contiguous()

        gt_depth = batch.get('depth', None)
        if gt_depth is not None:
            # Keep original dtype (float32 / uint16). C++ casts on GPU.
            if gt_depth.dtype == torch.float16:
                gt_depth = gt_depth.float()
            if len(gt_depth.shape) == 3:
                gt_depth = gt_depth.unsqueeze(-1)
            gt_depth = gt_depth.contiguous()

        gt_normal = batch.get('normal', None)
        if gt_normal is not None:
            # Keep original dtype (float32 / uint8). C++ does the 2x/255 - 1
            # scaling on GPU when uint8.
            if gt_normal.dtype == torch.float16:
                gt_normal = gt_normal.float()
            gt_normal = gt_normal.contiguous()

        # External mask -> gt_alpha (bool / uint8). Drives RGB mask AND alpha
        # supervision target in the slang kernel.
        gt_alpha = batch.get('mask', None)
        if gt_alpha is not None:
            # Coerce to bool so the C++ side hands a 1 byte/pixel buffer to the
            # kernel; this is small enough to keep as-is rather than convert.
            if gt_alpha.dtype != torch.bool:
                gt_alpha = gt_alpha > 0.5
            gt_alpha = gt_alpha.contiguous()

        if PROFILE_TRAIN_STEP: _t2 = _t()
        # --- Loss weights ---
        step = self.step
        cfg = self.config
        loss_weights = self._build_loss_weights(step)
        w_ssim = cfg.ssim_lambda
        num_loss_scales = cfg.num_loss_scales + 1
        compute_loss_map = cfg.use_loss_map or (cfg.compute_hessian_diagonal is not None)
        structure_only_loss_map = cfg.use_structure_only_loss_map

        sh_degree_to_use = step // max(cfg.sh_degree_warmup_every, 1)
        max_steps = self.trainer_config.num_iterations

        # --- Bilagrid: lazy per-type init on first iteration that needs it ---
        self._maybe_init_bilagrid(batch)
        optim_cfg = self.trainer_config.optimizer
        max_steps_lr = optim_cfg.max_steps if optim_cfg.max_steps is not None else max_steps
        # Pass the full per-image cam_idx tensor [C_batch] so the C++ kernels
        # can handle batches with mixed camera indices via indirect grid/param
        # lookup (no per-image gather, no per-image kernel launches).
        needs_cam_indices = (
            self._bilagrid_rgb_init or self._bilagrid_depth_init or
            self._bilagrid_normal_init or self._ppisp_init)
        if (needs_cam_indices and camera.metadata is not None
                and 'cam_idx' in camera.metadata):
            bilagrid_cam_indices = camera.metadata['cam_idx'].flatten().to(torch.int32).contiguous()
        else:
            bilagrid_cam_indices = None
        # When AdaGrad is enabled for bilagrids, swap to the dedicated AdaGrad
        # LR schedule names; the scheduler logic itself is unchanged.
        _bg_lr_keys = (
            ('bilagrid_adagrad',        'bilagrid_adagrad_depth',  'bilagrid_adagrad_normal')
            if cfg.use_adagrad_bilagrid_optim else
            ('bilagrid',                'bilagrid_depth',          'bilagrid_normal'))
        if self._bilagrid_rgb_init:
            bilagrid_lr_rgb = optim_cfg.get_scheduled_lr(_bg_lr_keys[0], step, max_steps_lr)
            bilagrid_tv_weight_rgb = cfg.bilagrid_tv_loss_weight
            bilagrid_shift_reg_weight_rgb = cfg.bilagrid_shift_reg_weight
            _period = max(int(cfg.bilagrid_shift_reg_ema_period), 1)
            bilagrid_shift_reg_beta_rgb = max(0.0, 1.0 - 1.0 / _period)
        else:
            bilagrid_lr_rgb = 0.0
            bilagrid_tv_weight_rgb = 0.0
            bilagrid_shift_reg_weight_rgb = 0.0
            bilagrid_shift_reg_beta_rgb = 0.0
        if self._bilagrid_depth_init:
            bilagrid_lr_depth = optim_cfg.get_scheduled_lr(_bg_lr_keys[1], step, max_steps_lr)
            bilagrid_tv_weight_depth = cfg.bilagrid_tv_loss_weight_geometry
        else:
            bilagrid_lr_depth = 0.0
            bilagrid_tv_weight_depth = 0.0
        if self._bilagrid_normal_init:
            bilagrid_lr_normal = optim_cfg.get_scheduled_lr(_bg_lr_keys[2], step, max_steps_lr)
            bilagrid_tv_weight_normal = cfg.bilagrid_tv_loss_weight_geometry
        else:
            bilagrid_lr_normal = 0.0
            bilagrid_tv_weight_normal = 0.0

        # Background blending per-step args.
        if cfg.background_mode == "noise":
            rw = min(step / max(cfg.background_noise_warmup, 1), 1.0)
            bg_randomize_weight = 1.0 - (1.0 - cfg.background_noise_pre_warmup) * (1.0 - rw)
            bg_lr_dc = 0.0
            bg_lr_sh = 0.0
        elif cfg.background_mode == "sh":
            bg_randomize_weight = 0.0
            bg_lr_dc = optim_cfg.get_scheduled_lr('background_dc', step, max_steps_lr)
            bg_lr_sh = optim_cfg.get_scheduled_lr('background_sh',    step, max_steps_lr)
        else:
            bg_randomize_weight = 0.0
            bg_lr_dc = 0.0
            bg_lr_sh = 0.0
        bg_seed = int(step) & 0x7FFFFFFF

        # PPISP: zero LR and reg weights when not yet initialized — C++ side
        # treats this as a no-op.
        if self._ppisp_init:
            ppisp_lr = optim_cfg.get_scheduled_lr('ppisp', step, max_steps_lr)
            ppisp_reg_exposure_mean   = cfg.ppisp_reg_exposure_mean
            ppisp_reg_vig_center      = cfg.ppisp_reg_vig_center
            ppisp_reg_vig_non_pos     = cfg.ppisp_reg_vig_non_pos
            ppisp_reg_vig_channel_var = cfg.ppisp_reg_vig_channel_var
            ppisp_reg_color_mean      = cfg.ppisp_reg_color_mean
            ppisp_reg_crf_channel_var = cfg.ppisp_reg_crf_channel_var
        else:
            ppisp_lr = 0.0
            ppisp_reg_exposure_mean = 0.0
            ppisp_reg_vig_center = 0.0
            ppisp_reg_vig_non_pos = 0.0
            ppisp_reg_vig_channel_var = 0.0
            ppisp_reg_color_mean = 0.0
            ppisp_reg_crf_channel_var = 0.0

        if PROFILE_TRAIN_STEP:
            _t3 = _t()
            # Sync before the C++ call so we can attribute its time cleanly:
            # presync = time waited for prior queued GPU work to finish;
            # ccall   = pure C++ engine_train_step (kernel launches; after
            #           option 3 the per-step D->H is async so this should be
            #           sub-ms, with gpu_drain growing instead).
            torch.cuda.synchronize()
            _t3b = _t()
        # --- Fused C++ training step ---
        loss_dict = self.core.engine_train_step(
            step, max_steps,
            sh_degree_to_use,
            W, H, camera_model,
            viewmats, intrins, dist_coeffs,
            gt_rgb, gt_depth, gt_normal, gt_alpha,
            loss_weights, w_ssim, num_loss_scales, compute_loss_map,
            structure_only_loss_map,
            self.config, self.trainer_config.optimizer,
            bilagrid_cam_indices=bilagrid_cam_indices,
            bilagrid_lr_rgb=bilagrid_lr_rgb,
            bilagrid_lr_depth=bilagrid_lr_depth,
            bilagrid_lr_normal=bilagrid_lr_normal,
            bilagrid_tv_weight_rgb=bilagrid_tv_weight_rgb,
            bilagrid_tv_weight_depth=bilagrid_tv_weight_depth,
            bilagrid_tv_weight_normal=bilagrid_tv_weight_normal,
            bilagrid_shift_reg_weight_rgb=bilagrid_shift_reg_weight_rgb,
            bilagrid_shift_reg_beta_rgb=bilagrid_shift_reg_beta_rgb,
            ppisp_lr=ppisp_lr,
            ppisp_reg_exposure_mean=ppisp_reg_exposure_mean,
            ppisp_reg_vig_center=ppisp_reg_vig_center,
            ppisp_reg_vig_non_pos=ppisp_reg_vig_non_pos,
            ppisp_reg_vig_channel_var=ppisp_reg_vig_channel_var,
            ppisp_reg_color_mean=ppisp_reg_color_mean,
            ppisp_reg_crf_channel_var=ppisp_reg_crf_channel_var,
            bg_lr_dc=bg_lr_dc,
            bg_lr_sh=bg_lr_sh,
            bg_randomize_weight=bg_randomize_weight,
            bg_seed=bg_seed,
            overexposure_reg_weight=cfg.overexposure_reg,
        )

        if PROFILE_TRAIN_STEP: _t4 = _t()
        # --- Verbose metrics ---
        for key, value in loss_dict.items():
            self.training_verboser.add_metric(key, value)
        self.training_verboser.add_metric("num_splats", self.core.cur_num_splats, last_only=True)
        self.training_verboser.add_metric("max_num_splats", self.core.max_num_splats, last_only=True)
        self.training_verboser.add_metric("num_sh", min(sh_degree_to_use, self.config.sh_degree), last_only=True)
        self.training_verboser.add_metric("max_num_sh", self.config.sh_degree, last_only=True)

        if PROFILE_TRAIN_STEP:
            _t5 = _t()
            if _prof["warmup"] > 0:
                _prof["warmup"] -= 1
            else:
                _prof["n"] += 1
                _prof["cam_prep"]        += (_t1 - _t0)
                _prof["gt_prep"]         += (_t2 - _t1)
                _prof["weights_prep"]    += (_t3 - _t2)
                _prof["presync"]         += (_t3b - _t3)
                _prof["ccall"]           += (_t4 - _t3b)
                _prof["post"]            += (_t5 - _t4)
            if _prof["n"] >= 100:
                n = _prof["n"]
                print(
                    f"[PROF_ETS n={n}] "
                    f"cam={_prof['cam_prep']/n/1e6:.3f}ms "
                    f"gt={_prof['gt_prep']/n/1e6:.3f}ms "
                    f"weights={_prof['weights_prep']/n/1e6:.3f}ms "
                    f"presync={_prof['presync']/n/1e6:.3f}ms "
                    f"ccall={_prof['ccall']/n/1e6:.3f}ms "
                    f"post={_prof['post']/n/1e6:.3f}ms",
                    flush=True,
                )
                for k in ("cam_prep", "gt_prep", "weights_prep", "presync", "ccall", "post"):
                    _prof[k] = 0
                _prof["n"] = 0

    def engine_train_step_managed(self, step: int):
        """DataManager-driven fused training step.

        The engine pulls the next batch (camera params + RGB/mask/depth/normal)
        from the C++ DataManager installed by Trainer; this method only computes
        the per-step Python config (loss weights, LR schedule, bilagrid / PPISP /
        background args) and dispatches to ``core.engine_train_step_managed``.

        Required precondition: ``Trainer._setup_cpp_data_manager()`` has been
        called, installing widths/heights/viewmats/intrins/dist_coeffs on the
        engine side at training-frame scale. This entrypoint
        therefore does NOT touch ``camera`` / ``batch`` — the on-engine state
        already reflects all per-camera baking.
        """
        cfg = self.config

        # Loss weights -- identical schedule to engine_train_step.
        loss_weights = self._build_loss_weights(step)
        w_ssim = cfg.ssim_lambda
        num_loss_scales = cfg.num_loss_scales + 1
        compute_loss_map = cfg.use_loss_map or (cfg.compute_hessian_diagonal is not None)
        structure_only_loss_map = cfg.use_structure_only_loss_map

        sh_degree_to_use = step // max(cfg.sh_degree_warmup_every, 1)
        max_steps = self.trainer_config.num_iterations
        optim_cfg = self.trainer_config.optimizer
        max_steps_lr = optim_cfg.max_steps if optim_cfg.max_steps is not None else max_steps

        # Bilagrid lazy init — feed sentinel-bearing dict so _maybe_init_bilagrid
        # picks up depth/normal availability without an actual batch tensor.
        synthetic_batch = {}
        if getattr(self, '_dm_has_depths', False):
            synthetic_batch['depth'] = True
        if getattr(self, '_dm_has_normals', False):
            synthetic_batch['normal'] = True
        self._maybe_init_bilagrid(synthetic_batch)

        _bg_lr_keys = (
            ('bilagrid_adagrad',        'bilagrid_adagrad_depth',  'bilagrid_adagrad_normal')
            if cfg.use_adagrad_bilagrid_optim else
            ('bilagrid',                'bilagrid_depth',          'bilagrid_normal'))
        if self._bilagrid_rgb_init:
            bilagrid_lr_rgb = optim_cfg.get_scheduled_lr(_bg_lr_keys[0], step, max_steps_lr)
            bilagrid_tv_weight_rgb = cfg.bilagrid_tv_loss_weight
            bilagrid_shift_reg_weight_rgb = cfg.bilagrid_shift_reg_weight
            _period = max(int(cfg.bilagrid_shift_reg_ema_period), 1)
            bilagrid_shift_reg_beta_rgb = max(0.0, 1.0 - 1.0 / _period)
        else:
            bilagrid_lr_rgb = 0.0
            bilagrid_tv_weight_rgb = 0.0
            bilagrid_shift_reg_weight_rgb = 0.0
            bilagrid_shift_reg_beta_rgb = 0.0
        if self._bilagrid_depth_init:
            bilagrid_lr_depth = optim_cfg.get_scheduled_lr(_bg_lr_keys[1], step, max_steps_lr)
            bilagrid_tv_weight_depth = cfg.bilagrid_tv_loss_weight_geometry
        else:
            bilagrid_lr_depth = 0.0
            bilagrid_tv_weight_depth = 0.0
        if self._bilagrid_normal_init:
            bilagrid_lr_normal = optim_cfg.get_scheduled_lr(_bg_lr_keys[2], step, max_steps_lr)
            bilagrid_tv_weight_normal = cfg.bilagrid_tv_loss_weight_geometry
        else:
            bilagrid_lr_normal = 0.0
            bilagrid_tv_weight_normal = 0.0

        # Background per-step args.
        if cfg.background_mode == "noise":
            rw = min(step / max(cfg.background_noise_warmup, 1), 1.0)
            bg_randomize_weight = 1.0 - (1.0 - cfg.background_noise_pre_warmup) * (1.0 - rw)
            bg_lr_dc = 0.0
            bg_lr_sh = 0.0
        elif cfg.background_mode == "sh":
            bg_randomize_weight = 0.0
            bg_lr_dc = optim_cfg.get_scheduled_lr('background_dc', step, max_steps_lr)
            bg_lr_sh = optim_cfg.get_scheduled_lr('background_sh',    step, max_steps_lr)
        else:
            bg_randomize_weight = 0.0
            bg_lr_dc = 0.0
            bg_lr_sh = 0.0
        bg_seed = int(step) & 0x7FFFFFFF

        # PPISP per-step args.
        if self._ppisp_init:
            ppisp_lr = optim_cfg.get_scheduled_lr('ppisp', step, max_steps_lr)
            ppisp_reg_exposure_mean   = cfg.ppisp_reg_exposure_mean
            ppisp_reg_vig_center      = cfg.ppisp_reg_vig_center
            ppisp_reg_vig_non_pos     = cfg.ppisp_reg_vig_non_pos
            ppisp_reg_vig_channel_var = cfg.ppisp_reg_vig_channel_var
            ppisp_reg_color_mean      = cfg.ppisp_reg_color_mean
            ppisp_reg_crf_channel_var = cfg.ppisp_reg_crf_channel_var
        else:
            ppisp_lr = 0.0
            ppisp_reg_exposure_mean = 0.0
            ppisp_reg_vig_center = 0.0
            ppisp_reg_vig_non_pos = 0.0
            ppisp_reg_vig_channel_var = 0.0
            ppisp_reg_color_mean = 0.0
            ppisp_reg_crf_channel_var = 0.0

        loss_dict = self.core.engine_train_step_managed(
            step, max_steps,
            sh_degree_to_use,
            loss_weights, w_ssim, num_loss_scales,
            compute_loss_map, structure_only_loss_map,
            self.config, self.trainer_config.optimizer,
            bilagrid_lr_rgb=bilagrid_lr_rgb,
            bilagrid_lr_depth=bilagrid_lr_depth,
            bilagrid_lr_normal=bilagrid_lr_normal,
            bilagrid_tv_weight_rgb=bilagrid_tv_weight_rgb,
            bilagrid_tv_weight_depth=bilagrid_tv_weight_depth,
            bilagrid_tv_weight_normal=bilagrid_tv_weight_normal,
            bilagrid_shift_reg_weight_rgb=bilagrid_shift_reg_weight_rgb,
            bilagrid_shift_reg_beta_rgb=bilagrid_shift_reg_beta_rgb,
            ppisp_lr=ppisp_lr,
            ppisp_reg_exposure_mean=ppisp_reg_exposure_mean,
            ppisp_reg_vig_center=ppisp_reg_vig_center,
            ppisp_reg_vig_non_pos=ppisp_reg_vig_non_pos,
            ppisp_reg_vig_channel_var=ppisp_reg_vig_channel_var,
            ppisp_reg_color_mean=ppisp_reg_color_mean,
            ppisp_reg_crf_channel_var=ppisp_reg_crf_channel_var,
            bg_lr_dc=bg_lr_dc,
            bg_lr_sh=bg_lr_sh,
            bg_randomize_weight=bg_randomize_weight,
            bg_seed=bg_seed,
            overexposure_reg_weight=cfg.overexposure_reg,
        )

        for key, value in loss_dict.items():
            self.training_verboser.add_metric(key, value)
        self.training_verboser.add_metric("num_splats",     self.core.cur_num_splats, last_only=True)
        self.training_verboser.add_metric("max_num_splats", self.core.max_num_splats, last_only=True)
        self.training_verboser.add_metric("num_sh",         min(sh_degree_to_use, self.config.sh_degree), last_only=True)
        self.training_verboser.add_metric("max_num_sh",     self.config.sh_degree, last_only=True)

    def get_loss_grad(self, outputs, batch, batch_size: int) -> Dict[str, torch.Tensor]:
        """Computes and returns the losses dict."""

        # Engine path: set training data, compute loss + backward in C++
        if self.core.primitive in ['3dgs', 'mip', '3dgut'] and not self.core.use_bvh:
            return self._engine_get_loss_grad(outputs, batch, batch_size)

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
        # gt_rgb = self.composite_with_background(self.get_gt_img(batch["image"]), outputs["background"])
        gt_rgb = batch["image"].cuda()  # TODO
        if gt_rgb.dtype == torch.uint8:
            gt_rgb = gt_rgb.float() / 255.0
        # TODO: linear and wide-gamut color spaces
        gt_rgb = gt_rgb[..., :3]  # TODO: RGBA
        gt_rgb = gt_rgb.squeeze(0)

        predicted_rgb = outputs["rgb"].cuda()
        predicted_rgb = torch.clip(predicted_rgb, 0.0, 1.0)

        # combined_rgb = torch.cat([gt_rgb, predicted_rgb], dim=1)
        combined_rgb = predicted_rgb

        from fused_bilagrid import color_correct
        corrected_rgb = color_correct(predicted_rgb, gt_rgb)

        # Switch images from [H, W, C] to [1, C, H, W] for metrics computations
        gt_rgb = torch.moveaxis(gt_rgb, -1, 0)[None, ...]
        predicted_rgb = torch.moveaxis(predicted_rgb, -1, 0)[None, ...]
        corrected_rgb = torch.moveaxis(corrected_rgb, -1, 0)[None, ...]

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

        cc_l1 = torch.abs(gt_rgb - corrected_rgb).mean()
        cc_psnr = self.psnr(gt_rgb, corrected_rgb)
        cc_ssim = self.ssim(gt_rgb, corrected_rgb)
        cc_ssim_torchmetrics = self.ssim_torchmetrics(gt_rgb, corrected_rgb)
        cc_lpips_alex = self.lpips_alex(gt_rgb, corrected_rgb)
        cc_lpips_vgg = self.lpips_vgg(gt_rgb, corrected_rgb)

        metrics_dict = {
            "l1": float(l1),
            "psnr": float(psnr),
            # "ssim_pytorch_msssim": float(ssim),
            # "ssim_torchmetrics": float(ssim_torchmetrics),
            "ssim": float(ssim_torchmetrics),
            "lpips_alex": float(lpips_alex),
            "lpips_vgg": float(lpips_vgg),
            "cc_l1": float(cc_l1),
            "cc_psnr": float(cc_psnr),
            # "cc_ssim_pytorch_msssim": float(cc_ssim),
            # "cc_ssim_torchmetrics": float(cc_ssim_torchmetrics),
            "cc_ssim": float(cc_ssim_torchmetrics),
            "cc_lpips_alex": float(cc_lpips_alex),
            "cc_lpips_vgg": float(cc_lpips_vgg),
            "gaussian_count": float(self.core.cur_num_splats),
        }

        images_dict = {"img": combined_rgb}

        return metrics_dict, images_dict
