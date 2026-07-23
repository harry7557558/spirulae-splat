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
    primitive: Literal["3dgs", "mip", "3dgut"] = "3dgs"
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
    loss_scale_min_pixels: int = 1920
    """If positive, overrides num_loss_scales per image based on resolution, in units of pixels.
        num_loss_scales is chosen so the smallest image dimension is halved down toward (but not below) this many pixels.
        e.g. with 2000: min dim 1999 -> num_loss_scales=0, 2000 -> 1, 4000 -> 2, 8000 -> 3, etc.
        Adapts per training step, so datasets with mixed image resolutions get the right count per image automatically."""
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
    split_batch: bool = True
    """Split the camera batch into one-camera sub-batches inside the C++ train step.
        Per-splat grads accumulate via atomicAdd across sub-batches; a single
        optim+densify pass at the end consumes the accumulator with grad_scale = 1/B.
        Drops peak VRAM for the immediate projection / rasterization buffers by
        roughly 1/B. Per-image grad magnitude vs regularization weight stays
        batch-size invariant. Not compatible with use_fused_proj_bwd_optim or with
        the warped train-step path."""
    quantization_level: int = 1
    """SH quantization level: a single int that selects one of two
        (param bits, optim bits) configurations.
            0 = off          : 32-bit param, fp32 optim
            1 = light        : 16-bit param, 8-bit packed optim (2 B / cell)
        Collapsing the prior independent param+optim bit controls into a
        single level minimizes the FPBO kernel instantiations."""
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
    min_init_fraction: float = 0.0
    """minimum fraction of splats out of cap_max at initialization"""
    refine_every: int = 100
    """Densify every this number of steps"""
    refine_start_iter: int = 500
    """Start densification at this number of steps"""
    refine_stop_num_iter: int = 5000
    """Stop densification at this number of steps before maximum number of training iterations"""
    refine_stop_iter: int = 25000
    """Densification runs until max(this, num_iterations - refine_stop_num_iter).
    Without this floor, runs shorter than refine_stop_num_iter would never densify at all
    (num_iterations - refine_stop_num_iter goes negative), which confuses users."""
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
    densify_score_blend_world_grad: float = 0.0
    """Blend weight `w` in [0, 1] between the image-space loss score and the
        world-space gradient score for densification. The per-step score is
        (image-space accum_weight)^(1-w) * (||dL/dmean_world|| * max post-exp
        world scale)^w. The world-grad term favors world-space-large splats
        (e.g. distant background in unbounded outdoor scenes) that the
        image-space score under-weights; the geometric blend is invariant to
        each score's global scale so no cross-normalization is needed.
        0 (default) = image-space score only, identical cost and behavior to
        before. 1 = world-grad score only; the per-pixel densification loss
        map (densify_loss_map_mode) is skipped entirely. In between, both
        scores are computed (one extra float per splat of VRAM)."""
    densify_loss_map_mode: Literal[
        "none",
        "loss_full",
        "ssim_full",
        "ssim_cs",
        "ssim_structure",
        "edge_aware",
        "robust_edge_aware",
    ] = "ssim_structure"
    """What gets accumulated into the per-pixel densification loss map. The
        loss map is read by raster bwd to weight the per-splat accum_weight.
        Only active when use_revised_densification. Modes:
        "none":              no loss map (uniform alpha*T accumulation).
        "loss_full":         per-pixel L1/L2 + auxiliary supervisory terms +
                             full SSIM (luminance*contrast*structure).
        "ssim_full":         full SSIM only.
        "ssim_cs":           contrast*structure SSIM (no luminance).
        "ssim_structure":    structure-only SSIM, biases toward pattern/edge
                             mismatches and ignores brightness/contrast errors.
        "edge_aware":        canny edge magnitude of GT rgb (Plenoxels-style,
                             https://arxiv.org/abs/2603.08661). Biases
                             densification toward GT edges directly, regardless
                             of how well the splats already reconstruct them.
        "robust_edge_aware": RobustNeRF-style Tukey biweight on the BT.601
                             luma of |render - GT|, capped at the per-image
                             `densify_robust_edge_aware_quantile`, then canny.
                             Near-zero where the render already matches GT,
                             zeroed past the quantile cutoff so distractor
                             pixels (people/cars/operator) don't pull splats
                             toward them, and luminance-shift tolerant since
                             a global DC residual has no spatial gradient.
        For num_loss_scales > 0 the map is computed per scale and the per-scale
        results are averaged (matches the multi-scale loss accumulation).
        Affects loss_map only; training gradients and scalar losses unchanged."""
    densify_robust_edge_aware_quantile: float = 0.9
    """Per-image quantile of the luma residual used as the Tukey biweight
        cutoff in `robust_edge_aware` mode. Pixels whose residual exceeds this
        quantile get zero densification weight. Lower values are more
        aggressive outlier rejection (good for distractor-heavy datasets);
        higher are more permissive (good for clean datasets where real edges
        may produce large residuals). Ignored unless
        `densify_loss_map_mode == "robust_edge_aware"`."""
    use_long_axis_split: bool = True
    """whether to use long-axis split described in https://arxiv.org/abs/2508.12313 for relocation and sample add.
        When combined with use_revised_densification, this can give less blurry background details for unbounded outdoor scenes."""
    long_axis_split_opacity_k: Tuple[float, float, int] = (0.6, 0.6, 4500)
    """opacity split factor `k` for long-axis split, as (initial, final, warmup_steps).
        Each split child keeps opacity `logit^-1(k / (1 + exp(-logit_opacity) - k))`; `k` is linearly
        scheduled from `initial` to `final` over the first `warmup_steps` training steps, then held at `final`.
        Larger `k` preserves more opacity per child (denser, sharper); smaller `k` fades children faster."""
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
    bilagrid_type: Literal["affine", "ppisp", "loglinear"] = "ppisp"
    """What the bilateral grid predicts.
        affine: 4x3 matrix per original bilateral grid.
        ppisp: PPISP exposure and color parameters, generally gives less color shift but can be less numerically stable.
        loglinear: 3x3 linear transformation matrix with log-encoded diagonals, balances color shift and numerical stability.
    """
    use_bilateral_grid_for_geometry: bool = True
    """If True, use bilateral grid for depth and normal (e.g. AI generated biased ones)"""
    bilagrid_shape_geometry: Tuple[int, int, int] = (8, 8, 4)
    """Shape of the bilateral grid for depth and normal (X, Y, W)"""
    use_adagrad_bilagrid_optim: bool = True
    """Use AdaGrad (lr_decay=0, weight_decay=0, initial_accumulator_value=0,
       eps=1e-15) instead of Adam for all bilateral-grid parameters (RGB +
       depth + normal). When True, the bilagrid LR fields read from
       OptimizerConfig switch to ``bilagrid_adagrad_*_lr``. Bilagrid bit
       depths are coupled to `quantization_level`: level 0 = fp32; level 1 =
       16-bit value + 8x2-bit optimizer state across all three bilagrid types."""
    bilagrid_tv_loss_weight: float = 10.0
    """Total variation loss weight for bilateral grid used for radiance"""
    color_shift_reg_weight: float = 0.0
    """Color-shift regularizer for the combined bilagrid + PPISP color
        transform. Penalizes the dataset-wide mean of sign(post - pre) per
        channel, where `pre` is the splat-side rendered RGB (input to whichever
        of bilagrid / PPISP runs first) and `post` is the final image fed to
        the photometric loss: R = w * ||EMA[mean_p sign(post - pre)]||^2. The
        gradient is injected on the POST-transforms v_render_rgb buffer and
        flows through each transform's vjp into its parameters (and as a small
        leak into the splats). Active when at least one of bilagrid_rgb / PPISP
        is enabled; 0 disables. Typical values: 0.01--1.0."""
    color_shift_reg_ema_period: int = 750
    """EMA decay length for the color-shift regularizer, in BATCHES.
        beta = max(0, 1 - 1/period). Should be roughly the number of batches
        per epoch so the EMA estimates the dataset-wide mean. Ignored when
        color_shift_reg_weight = 0."""
    bilagrid_tv_loss_weight_geometry: float = 10.0
    """Total variation loss weight for bilateral grid used for geometry"""
    use_ppisp: bool = True
    """If True, use the PPISP model (https://research.nvidia.com/labs/sil/projects/ppisp/) to handle per-pixel color distortions."""
    ppisp_param_type: Literal["original", "rqs", "no_crf"] = "no_crf"
    """Parameterization for PPISP. "original" implements the original paper,
        "rqs" uses a parameterization that is more friendly to optimization and can produce better results in darker areas,
        "no_crf" omits the tone-curve (CRF) stage entirely and just clamps colors to [0,1] after exposure/vignetting/color correction."""
    use_adagrad_ppisp_optim: bool = True
    """Use unscheduled AdaGrad (lr_decay=0, weight_decay=0,
       initial_accumulator_value=0, eps=1e-15) instead of Adam for the PPISP
       parameter table. When True, the PPISP LR reads from
       ``OptimizerConfig.ppisp_adagrad_lr`` (constant) instead of the scheduled
       ``ppisp_lr``. No quantization path either way."""
    apply_ppisp_before_bilagrid: bool = True
    """When True, the PPISP forward runs BEFORE the RGB bilagrid (and PPISP
       backward runs AFTER bilagrid backward), i.e. render -> PPISP -> bilagrid
       -> loss. Otherwise: render -> bilagrid -> PPISP -> loss. Only meaningful
       when both ``use_bilateral_grid`` and ``use_ppisp`` are enabled."""
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
    reg_warmup_length: int = 0
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
    input_depth_is_ray_depth: bool = False
    """Whether the input/supervision depth maps store ray depth (Euclidean distance along the camera
        ray) rather than linear (z) depth. The rasterizer renders ray depth, so when this is False
        (the common case, e.g. most foundation-model depths) the GT depth is converted from linear to
        ray depth in place on the GPU before the depth bilateral grid / loss. Set True for depth maps
        already in ray depth, e.g. >180deg fisheye captures where linear depth is ill-defined."""
    normal_supervision_weight: float = 0.01
    """Weight for normal supervision by comparing normal from rendered depth with normal from depth predicted by a foundation model"""

    # Median-depth losses. The median depth is the depth where accumulated
    # transmittance crosses 1/2 (sharper than the expected/mean depth). These
    # are emitted by the rasterizer only when at least one weight below is > 0.
    mean_median_depth_weight: float = 0.0
    """L1 between the mean (expected) depth and the median depth, where both are nonzero."""
    median_depth_normal_reg_weight: float = 0.0
    """normal_loss between the normal from the median depth and the normal from the mean depth."""
    median_normal_supervision_weight: float = 0.0
    """normal_loss between the normal from the median depth and the reference (foundation-model) normal."""
    median_render_normal_reg_weight: float = 0.0
    """normal_loss between the normal from the median depth and the rendered normal (placeholder until render_normal exists)."""
    median_warmup: int = 6000
    """Linear warmup length (steps) shared by all four median-depth loss weights: each ramps 0->full over this many steps."""

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


# Map densify_loss_map_mode literal -> int passed to the C++ engine. Must
# stay in sync with DensifyLossMapMode in PerPixelLoss.cuh.
_DENSIFY_LOSS_MAP_MODE_TO_INT = {
    "none":              0,
    "loss_full":         1,
    "ssim_full":         2,
    "ssim_cs":           3,
    "ssim_structure":    4,
    "edge_aware":        5,
    "robust_edge_aware": 6,
}


class SpirulaeSplatModel(torch.nn.Module):
    """Template Model."""

    config: SpirulaeSplatModelConfig

    def __init__(
        self,
        trainer_config: 'spirulae_splat.modules.trainer.TrainerConfig',
        seed_points: Optional[Tuple[torch.Tensor, torch.Tensor]] = None,
        cameras: Optional[Cameras] = None,
        attach_engine: bool = False,
        num_train_data: Optional[int] = None,
    ):
        """attach_engine=True binds to a world the engine already holds.

        The native TrainerSession (src/app/TrainerCore.cpp) does the seeding,
        the background / color-space init and the DataManager setup, so in
        attach mode this class contributes only what stays in Python: host
        parameter tensors at the right shapes (the destination for
        engine_copy_splats_to_host, and the layout resume_adapt targets), the
        resolved color config, and the eval / viewer render paths. Repeating
        any of the init here would overwrite the session's work.
        """
        super().__init__()

        self.trainer_config = trainer_config
        self.config = trainer_config.model  # type: SpirulaeSplatModelConfig
        self._attach_engine = bool(attach_engine)

        self.seed_points = (seed_points['points3D_xyz'], seed_points['points3D_rgb'])
        self.cameras = cameras

        if num_train_data is not None:
            # Exact POST-split camera count, from the session's camera bake.
            # Sizes the bilagrid / PPISP tables.
            self.num_train_data = int(num_train_data)
        else:
            self.num_train_data = len(self.cameras)
            # Worst-case upper bound: every camera split 6x. Callers that know
            # the real post-split count pass num_train_data instead.
            self.num_train_data *= 6

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
            split_batch=self.config.split_batch,
            quantization_level=self.config.quantization_level,
            upload=not self._attach_engine,
        )

        if self._attach_engine:
            # setup_engine() already ran seeding, background and color space.
            # Pull the device world into the host parameters so gauss_params,
            # checkpointing and resume see the real splats.
            self._engine_bg_mode = self.config.background_mode \
                if self.config.background_mode in ("noise", "sh") else "none"
            self.core.cur_num_splats = int(_C.engine_get_cur_num_splats())
            self.core.engine_sync_splats_to_host()
            return

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

    def _resolve_color_config(self):
        """Resolve the Optional color-space fields to concrete values.

        Port-mate of TrainerCore's resolve_color(): both read the same raw
        config and must land on the same answer, so this runs in attach mode
        too -- the engine side is already resolved, but the Python eval /
        viewer paths read these fields.
        """
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

    def _seed_count(self):
        """(num_points, capacity) without doing the seeding.

        Mirrors the resolution in populate_modules below and in TrainerCore's
        seed_splats(); attach mode needs the shapes, not the values.
        """
        min_init = max(int(min(self.config.min_init_fraction, 1.0) * self.config.cap_max), 1)
        min_init = min(min_init, self.config.cap_max)
        n_src = len(self.seed_points[0])
        num_points = min(n_src, self.config.cap_max) if n_src >= min_init else min_init
        capacity = num_points
        if self.config.use_mcmc and self.config.preallocate_splat_tensors:
            capacity = max(num_points, self.config.cap_max)
        return num_points, capacity

    def _populate_modules_attached(self):
        """Allocate host parameters at the engine's layout, no seeding."""
        self._resolve_color_config()
        num_points, capacity = self._seed_count()
        engine_n = int(_C.engine_get_cur_num_splats())
        if engine_n != num_points:
            raise RuntimeError(
                f"attach_engine: the engine holds {engine_n} splats but this "
                f"config resolves to {num_points}. The session and the model "
                f"disagree about seeding -- they must be built from the same "
                f"config and the same point cloud.")

        dim_sh = (self.config.sh_degree + 1) ** 2
        z = lambda *shape: torch.nn.Parameter(torch.zeros(*shape))
        self.gauss_params = {
            'means':       z(capacity, 3),
            'quats':       z(capacity, 4),
            'scales':      z(capacity, 3),
            'opacities':   z(capacity, 1),
            'features_dc': z(capacity, 3),
            'features_sh': z(capacity, dim_sh - 1, 3),
        }
        optim_info = {'num_splats': num_points} \
            if (self.config.use_mcmc and self.config.preallocate_splat_tensors) else {}
        for value in self.gauss_params.values():
            value.optim_info = {**optim_info}

        self.xys_grad_norm = None
        self.ch_grad_norm = None
        self.max_2Dsize = None
        # Set by the session's setup_engine(); the trainer copies its RunState
        # over right after construction.
        self._bilagrid_rgb_init = False
        self._bilagrid_depth_init = False
        self._bilagrid_normal_init = False
        self._ppisp_init = False
        self.training_verboser = TrainingVerbose()
        self.step = 0
        self._train_batch_size = 1

    def populate_modules(self):
        if self._attach_engine:
            return self._populate_modules_attached()
        if self.seed_points is not None:
            min_init = max(int(min(self.config.min_init_fraction, 1.0) * self.config.cap_max), 1)
            min_init = min(min_init, self.config.cap_max)
            if len(self.seed_points[0]) < min_init:
                distances, indices_1 = self.k_nearest_neighbor(self.seed_points[0], 2)
                indices = torch.arange(min_init) % len(self.seed_points[0])
                self.seed_points = [
                    torch.lerp(t[indices], t[indices_1[indices, 0]], 0.5 * torch.randn_like(t[indices]))
                    if i == 0 else t[indices]
                    for (i, t) in enumerate(self.seed_points)
                ]
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
            if (seed_color == seed_color[0][0]).all():
                seed_color = torch.rand_like(seed_color)
        else:
            seed_color = torch.rand(num_points, 3)

        self._resolve_color_config()
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
        if CameraType.EQUIRECTANGULAR.value in self.cameras:
            # cov_scale_init projects each splat through a pinhole/fisheye model
            # to bound its initial scale by its pixel footprint; there is no
            # equirectangular path yet. Skip the suppression rather than crash so
            # direct equirectangular training still works (this is only an init
            # heuristic). TODO: equirectangular cov_scale_init.
            print("suppress_initial_scales: skipped (equirectangular not supported)")
            return scales

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
                use_adagrad=cfg.use_adagrad_bilagrid_optim,
            )
            self._bilagrid_normal_init = True

        # PPISP (RGB only): config flag + positive base LR. Always has supervision
        # via the rendering loss; no dataset-side gating needed.
        _ppisp_lr_init = (optim_cfg.ppisp_adagrad_lr
                          if cfg.use_adagrad_ppisp_optim else optim_cfg.ppisp_lr)
        if (not self._ppisp_init
                and cfg.use_ppisp
                and _ppisp_lr_init > 0.0):
            self.core.engine_init_ppisp(
                self.num_train_data,
                param_type=cfg.ppisp_param_type,
                use_adagrad=cfg.use_adagrad_ppisp_optim,
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
            # Enable the eval/viewer distortion render when any distortion
            # regularizer is configured, or on-demand when the viewer requests a
            # distortion buffer. Which channels are actually produced is decided
            # in C++ from the primitive's rendered channels (normal is honored
            # only once a normal-rendering primitive exists).
            output_distortion=any(c != 0.0 for c in self.get_2dgs_reg_weights()[0])
                or (_want is not None and any('distortion' in (k or '') for k in _want)),
            output_median=any([
                self.config.mean_median_depth_weight > 0.0,
                self.config.median_depth_normal_reg_weight > 0.0,
                self.config.median_normal_supervision_weight > 0.0,
                self.config.median_render_normal_reg_weight > 0.0,
            ]),
            **kwargs,
        )
        self.core.forward()
        rgbdn = self.core.render_colors
        Ts = self.core.render_Ts
        median_im = getattr(self.core, 'render_median', None)
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
                True,  # is_ray_depth
                self.core._tv(depth_cuda),
                self.core._tv(depth_normal_cuda),
            )
            depth_normal = depth_normal_cuda.cpu()

        # Median depth -> normal, same path as depth_normal above.
        median_normal = None
        if (median_im is not None and not self.training
                and (_want is None or 'normal_median' in _want)):
            median_cuda = median_im.cuda().contiguous() if not median_im.is_cuda else median_im.contiguous()
            intrins_cuda = intrins.cuda().contiguous() if not intrins.is_cuda else intrins.contiguous()
            if camera.distortion_params is not None:
                dist_coeffs_cuda = camera.distortion_params.cuda().contiguous()
            else:
                dist_coeffs_cuda = torch.zeros(len(camera), 10, dtype=torch.float32, device='cuda')
            median_normal_cuda = torch.empty(
                *median_cuda.shape[:-1], 3, dtype=torch.float32, device='cuda')
            _C.depth_to_normal_forward(
                camera.camera_type[0].upper(),
                self.core._tv(intrins_cuda),
                self.core._tv(dist_coeffs_cuda),
                True,  # is_ray_depth
                self.core._tv(median_cuda),
                self.core._tv(median_normal_cuda),
            )
            median_normal = median_normal_cuda.cpu()

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
        if median_im is not None:
            outputs["depth_median"] = median_im
        if median_normal is not None:
            outputs["normal_median"] = median_normal

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
                if "depth_median" in outputs:
                    outputs["depth_median"] = outputs["depth_median"] / self.config.relative_scale
            if "normal" in outputs:
                outputs["normal"] = 0.5+0.5*outputs["normal"]
            if "depth_normal" in outputs:
                outputs["depth_normal"] = 0.5+0.5*outputs["depth_normal"]
            if "normal_median" in outputs:
                outputs["normal_median"] = 0.5+0.5*outputs["normal_median"]
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

    # ---- The training path that used to live here is gone ----------------
    # engine_train_step / engine_train_step_managed / _build_loss_weights /
    # optim_step (and the unreachable _engine_get_loss_grad / get_loss_grad
    # that hung off them) were a second implementation of what
    # src/app/TrainerCore.cpp does -- its own header comment described itself
    # as a line-by-line port of them. `Trainer` now drives the C++ one through
    # `_C.TrainerSession` (docs/restructure-proposal.md §4.3), and
    # tests/python/step_config_golden.json freezes the values the two were
    # proven to agree on. Per-step training logic belongs in
    # TrainerCore::build_step_config; do not re-add it here.

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
            # separately by torch.cuda.memory_allocated). Each pool buffer's
            # bucket is now the authoritative category tag carried by its
            # pool-slot metadata (PoolSlots.h), returned per-buffer by the
            # engine — no more guessing from the key prefix here:
            #   splat       : world params + grads + optim states (per-Gaussian)
            #   splat x img : projection outputs/gradients, tile intersection
            #   image       : image-space tensors (renders, GT copies, grads)
            #   appearance  : bilagrid / background-SH / PPISP / color transforms
            #   viewer      : interactive viewer + visualizer caches/scratch
            #   other       : everything else (camera tables, tiny scratches)
            from spirulae_splat.splat.cuda import _C
            buckets = {'splat': 0.0, 'splat x img': 0.0, 'image': 0.0,
                       'appearance': 0.0, 'viewer': 0.0, 'other': 0.0}
            GiB = 1024 ** 3
            for key, category, _used, cap in _C.engine_get_pool_breakdown_categorized():
                buckets[category] = buckets.get(category, 0.0) + cap / GiB
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
