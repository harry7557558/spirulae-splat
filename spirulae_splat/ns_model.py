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
from torch.nn.functional import normalize
from pytorch_msssim import SSIM

from spirulae_splat.splat._torch_impl import quat_to_rotmat
from spirulae_splat.splat import (
    project_gaussians,
    rasterize_gaussians_simple,
    rasterize_gaussians_depth,
    rasterize_gaussians,
    rasterize_gaussians_simplified,
    depth_to_points,
    depth_to_normal,
)
from spirulae_splat.splat.background_sh import render_background_sh
from spirulae_splat.splat.sh import num_sh_bases, spherical_harmonics
from spirulae_splat.strategy import DefaultStrategy

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
    import torch.nn.functional as tf

    weight = (1.0 / (d * d)) * torch.ones((1, 1, d, d), dtype=torch.float32, device=image.device)
    return tf.conv2d(image.float().permute(2, 0, 1)[:, None, ...], weight, stride=d).squeeze(1).permute(1, 2, 0).to(image)


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
    ssim_lambda: float = 0.4  # 0.2
    """weight of ssim loss"""
    output_depth_during_training: bool = False
    """If True, output depth during training. Otherwise, only output depth during evaluation."""
    use_camera_optimizer: bool = True
    """Whether to use camera optimizer"""
    camera_optimizer: CameraOptimizerConfig = field(default_factory=lambda: CameraOptimizerConfig(mode="SO3xR3"))
    """Config of the camera optimizer to use"""
    kernel_radius: float = 1.0
    """Radius of the splatting kernel, 3.0 for Gaussian and 1.0 for polynomial"""
    loss_scale: float = 0.003
    """Scaling for loss values to normalize gradient"""
    target_absgrad: float = 5.5e-6
    """Will auto tune loss_scale to fit absgrad to this number, higher gives more splats"""

    # classial control
    cull_alpha_thresh: float = 0.1
    """threshold of opacity for culling gaussians. One can set it to a lower value (e.g. 0.005) for higher quality."""
    cull_scale_thresh: float = 0.5
    """threshold of scale for culling huge gaussians"""
    cull_anisotropy_thresh: float = np.inf
    """threshold of quotient of scale for culling long thin gaussians"""
    cull_grad_thresh: float = 0.0001  # 0.0003 | 0.0001
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
    """
    adaptive_exposure_start_iter: int = 1000
    """Start adaptive exposure at this number of steps"""

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
    mcmc_opacity_reg: float = 0.001
    """Opacity regularization from MCMC
       Lower usually gives more accurate geometry, original MCMC uses 0.01"""
    mcmc_scale_reg: float = 0.002
    """Scale regularization from MCMC, original MCMC uses 0.01"""
    exposure_reg_splat: float = 0.0
    """Make sure image look right in auto exposure mode, match mean splat color with mean image color"""
    exposure_reg_image: float = 0.1
    """Between 0 and 1; For exposure regularization, include this fraction of L1 loss between GT image and image before exposure adjustment"""
    exposure_reg_param: float = 0.002
    """Make sure image look right in auto exposure mode,
       Use an L2 cost to match exposure parameter to the parameter at no exposure adjustment;
       (may be summed/averaged across a batch)"""


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
        if not self.config.use_anisotropy:
            anisotropies.requires_grad_(False)

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

        # Strategy for GS densification
        reset_every = self.config.reset_alpha_every * self.config.refine_every
        pause_refine_after_reset = min(self.num_train_data, reset_every//2) + self.config.refine_every
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
            refine_stop_iter=self.config.stop_split_at,
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
        if not self.config.use_anisotropy:
            self.gauss_params['anisotropies'].requires_grad_(False)
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

    def get_background_image(self, camera: Cameras):
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

        fx, fy = camera.fx[0].item(), camera.fy[0].item()
        cx, cy = camera.cx[0].item(), camera.cy[0].item()
        rotation = camera.camera_to_worlds[0][:3, :3]
        sh_coeffs = torch.cat((self.background_color.unsqueeze(0), self.background_sh), dim=0)  # [(deg+1)^2, 3]
        return render_background_sh(W, H, (fx, fy, cx, cy), rotation, sh_degree+1, sh_coeffs, 16)

    @staticmethod
    def get_empty_outputs(width: int, height: int, background: torch.Tensor) -> Dict[str, Union[torch.Tensor, List]]:
        return {}

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

        timerr.start()

        if self.training and self.config.use_camera_optimizer:
            optimized_camera_to_world = self.camera_optimizer.apply_to_camera(camera)[0, ...]
        else:
            optimized_camera_to_world = camera.camera_to_worlds[0, ...]
        timerr.mark(".")  # 500us

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
        W, H = int(camera.width.item()), int(camera.height.item())
        intrins = (camera.fx.item(), camera.fy.item(), camera.cx.item(), camera.cy.item())
        camera.rescale_output_resolution(camera_downscale)
        timerr.mark(".")  # 450us-800us

        R = optimized_camera_to_world[:3, :3]  # 3 x 3
        T = optimized_camera_to_world[:3, 3:4]  # 3 x 1
        R_edit = torch.diag(torch.tensor([1, -1, -1], device=self.device, dtype=R.dtype))
        R = R @ R_edit
        R_inv = R.T
        T_inv = -R_inv @ T
        viewmat = torch.eye(4, device=R.device, dtype=R.dtype)
        viewmat[:3, :3] = R_inv
        viewmat[:3, 3:4] = T_inv
        self.last_size = (H, W)
        timerr.mark("pre")  # 300us-500us

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
        timerr.mark("crop")  # 80us-200us

        quats_crop = normalize(quats_crop, dim=-1)
        timerr.mark("norm")  # 100us-250us

        if self.config.sh_degree > 0:
            viewdirs = means_crop.detach() - optimized_camera_to_world.detach()[:3, 3]  # (N, 3)
            n = min(self.step // self.config.sh_degree_interval, self.config.sh_degree)
            rgbs = spherical_harmonics(n, viewdirs, colors_crop)  # input unnormalized viewdirs
            rgbs = torch.clamp(rgbs + 0.5, min=0.0)  # type: ignore
        else:
            rgbs = colors_crop[:, 0, :]
        opacities = self.config.max_opacity * torch.sigmoid(opacities_crop)
        timerr.mark("color")  # 200us-300us

        BLOCK_WIDTH = 16
        (
            positions, axes_u, axes_v,
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
        timerr.mark("project")  # 150us-250us

        # slower but more capable two-pass rendering
        if False or (
            self.config.depth_reg_pairwise_factor < 1.0 or \
            self.config.ch_degree_r * (2*self.config.ch_degree_phi+1) > 0 or \
            self.config.depth_mode != "mean"
        ):

            depth_im_ref = rasterize_gaussians_depth(
                positions, axes_u, axes_v,
                opacities, anisotropies_crop,
                bounds, num_tiles_hit,
                intrins, H, W, BLOCK_WIDTH,
                self.config.depth_mode
            )
            depth_im_ref = torch.where(
                depth_im_ref > 0.0, depth_im_ref,
                torch.amax(depth_im_ref).detach()
            ).contiguous()
            timerr.mark("depth")  # 700us-1200us

            # main rasterization
            ch_degree = self.step // self.config.ch_degree_interval
            (rgb, alpha, depth_im, normal_im, reg_depth) \
             = rasterize_gaussians(  # type: ignore
                positions, axes_u, axes_v,
                rgbs,
                self.config.ch_degree_r, min(ch_degree, self.config.ch_degree_r),
                self.config.ch_degree_phi, min(ch_degree, self.config.ch_degree_phi),
                features_ch_crop,
                opacities, anisotropies_crop,
                depth_im_ref,
                # background_color,
                self.config.depth_reg_pairwise_factor,
                bounds, num_tiles_hit,
                intrins, H, W, BLOCK_WIDTH,
            )  # type: ignore
            timerr.mark("render")  # 650us-1400us

        # fast one-pass rendering
        else:
            (rgb, alpha, depth_im, normal_im, reg_depth) \
             = rasterize_gaussians_simplified(
                positions, axes_u, axes_v,
                rgbs,
                opacities, anisotropies_crop,
                bounds, num_tiles_hit,
                intrins, H, W, BLOCK_WIDTH,
            )
            depth_im_ref = torch.where(
                alpha > 0.0, depth_im[..., :1] / alpha,
                torch.amax(depth_im[..., :1]).detach()
            ).contiguous()
            timerr.mark("render")  # 800us-1500us

        anisotropy_vis = None
        if self.config.output_depth_during_training or not self.training:
        # if False:
            with torch.no_grad():

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

        else:
            self.info = {
                "width": W,
                "height": H,
                "n_cameras": camera.shape[0],
                "radii": self.config.kernel_radius/3.0 * \
                      num_tiles_hit.unsqueeze(0)**0.5 / 2 * BLOCK_WIDTH,
                "means2d": positions,
            }
        timerr.mark("post")  # 150us-300us

        # blend with background
        background = self.background_color.reshape((1, 1, 3))
        if self.config.background_sh_degree > 0 or True:
            background = self.get_background_image(camera)
            rgb = rgb + (1.0 - alpha) * background

        # saturate pixel value
        if not self.training:
            rgb = torch.clip(rgb, 0.0, 1.0)
        else:
            rgb = torch.clip(rgb, 0.0, None)
            # rgb = saturate_keep_gradient(rgb, 0.0, None)

        timerr.mark("bkg")  # 250us-450us

        # normal regularization
        normal_im = normalize(normal_im, dim=-1)
        depth_normal, alpha_diffused = depth_to_normal(depth_im_ref, None, intrins, True, alpha)
        reg_normal = 1.0 - (depth_normal * normal_im).sum(-1, True)
        timerr.end("normal_reg")  # 500us-900us
        # -> 4000us-7000us median depth, 3500us-6000us one pass

        return {
            "rgb": rgb,
            "depth": depth_im_ref,
            "depth_normal": 0.5+0.5*depth_normal,
            "render_normal": 0.5+0.5*normal_im,
            "reg_depth": torch.sqrt(torch.relu(reg_depth*alpha)+1e-8) / depth_im_ref.clip(1e-3),
            "reg_normal": torch.sqrt(torch.relu(reg_normal*alpha)+1e-8) * alpha_diffused,
            # "reg_depth": reg_depth / depth_im_ref.clip(1e-3)**2,
            # "reg_normal": reg_normal * alpha_diffused,
            "anisotropy_vis": anisotropy_vis,
            "alpha": alpha,
            "background": background,
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
        timerl.start()
        gt_img_rgba = self.get_gt_img(batch["image"])
        gt_img = self.composite_with_background(gt_img_rgba, outputs["background"])
        pred_img = outputs["rgb"]

        # alpha channel for bounded objects - apply a cost on rendered alpha
        alpha_loss = 0.0
        if gt_img_rgba.shape[2] == 4:
            alpha = gt_img_rgba[..., -1].unsqueeze(-1)
            alpha_loss = alpha_loss + torch.relu(outputs['alpha']-alpha).mean()
            print(alpha_loss.item())

        # separate mask for dynamic objects, text, etc. - simply don't consider it when evaluating loss
        if "mask" in batch:
            # batch["mask"] : [H, W, 1]
            mask = self._downscale_if_required(batch["mask"])
            mask = mask.float().to(self.device)
            assert mask.shape[:2] == gt_img.shape[:2] == pred_img.shape[:2]
            # can be little bit sketchy for the SSIM loss
            gt_img = torch.lerp(outputs["background"], gt_img, mask)
            pred_img = torch.lerp(outputs["background"], pred_img, mask)
        timerl.mark("alpha")  # ~100us

        # handle exposure
        exposure_reg = 0.0
        exposure_param_reg = 0.0
        pred_img_e = pred_img
        if self.config.adaptive_exposure_mode is not None and \
            self.step > self.config.adaptive_exposure_start_iter:
            eeps = 0.5/255.0

            # gt ~ k * pred
            if self.config.adaptive_exposure_mode == "linear":
                k = (gt_img.mean()+eeps) / (pred_img.mean()+eeps)
                pred_img_e = k * (pred_img+eeps) - eeps
                # print(k.item())
                exposure_param_reg = (k-1.0)**2

            # log(gt) ~ k * log(pred) + b
            elif self.config.adaptive_exposure_mode == "log_linear":
                x = torch.log(pred_img+eeps)
                y = torch.log(gt_img+eeps)
                x_mean, y_mean = x.mean(), y.mean()
                x2_mean, xy_mean = (x**2).mean(), (x*y).mean()
                m = (xy_mean - x_mean * y_mean) / (x2_mean - x_mean**2)
                b = y_mean - m * x_mean
                pred_img_e = torch.exp(m * x + b) - eeps
                # print(m.item(), b.item())
                exposure_param_reg = (m-1.0)**2 + b**2

            # gt ~ k * pred, per channel
            if self.config.adaptive_exposure_mode == "channel":
                x = pred_img.reshape(-1, 3)
                y = gt_img.reshape(-1, 3)
                k = (y.mean(0, True)+eeps) / (x.mean(0, True)+eeps)
                pred_img_e = k * (pred_img+eeps) - eeps
                # print(k.detach())
                exposure_param_reg = torch.mean((k-1.0)**2)

            # log(gt) ~ k * log(pred) + b, per channel
            if self.config.adaptive_exposure_mode == "log_channel":
                x = torch.log(pred_img + eeps).reshape(-1, 3)
                y = torch.log(gt_img + eeps).reshape(-1, 3)
                x_mean, y_mean = x.mean(0, True), y.mean(0, True)
                x2_mean, xy_mean = (x**2).mean(0, True), (x*y).mean(0, True)
                m = (xy_mean - x_mean * y_mean) / (x2_mean - x_mean**2)
                b = y_mean - m * x_mean
                pred_img_e = torch.exp(m * x + b).reshape(pred_img_e.shape) - eeps
                # print(m.detach(), b.detach())
                exposure_param_reg = torch.mean((m-1.0)**2 + b**2)

            # gt ~ A * pred
            elif self.config.adaptive_exposure_mode == "affine":
                x = (pred_img + eeps).reshape(-1, 3)
                y = (gt_img + eeps).reshape(-1, 3)
                xTy, xTx = x.T @ y, x.T @ x
                try:
                    A = torch.linalg.inv(xTx) @ xTy
                except torch._C._LinAlgError:
                    pass
                else:
                    pred_img_e = (x @ A).reshape(pred_img.shape) - eeps
                    exposure_param_reg = torch.mean((A-torch.eye(3).to(A))**2)

            # log(gt) ~ A * log(pred) + b
            elif self.config.adaptive_exposure_mode == "log_affine":
                x = torch.log(pred_img + eeps).reshape(-1, 3)
                y = torch.log(gt_img + eeps).reshape(-1, 3)
                mean_x = torch.mean(x, dim=0, keepdim=True)
                mean_y = torch.mean(y, dim=0, keepdim=True)
                centered_x = x - mean_x
                centered_y = y - mean_y
                xTy = centered_x.T @ centered_y
                xTx = centered_x.T @ centered_x
                try:
                    A = torch.linalg.inv(xTx) @ xTy
                except torch._C._LinAlgError:
                    pass
                else:
                    B = mean_y - mean_x @ A
                    pred_img_e = torch.exp(x @ A + B).reshape(pred_img.shape) - eeps
                    exposure_param_reg = 0.5*(torch.mean((A-torch.eye(3).to(A))**2) + torch.mean(B**2))
            
            timerl.mark("exposure")  # 200us-500us

            # keep good value for natural image
            if self.config.exposure_reg_splat > 0.0:
                if not hasattr(self, '_exposure_info') or self._exposure_info_n < 4096:
                    # get mean and covariance of image colors
                    gt_color_flat = gt_img[..., :3].reshape((-1, 3))
                    gt_mean = gt_color_flat.mean(0).detach()
                    gt_mean2 = torch.cov(gt_color_flat.T).detach()
                    if not hasattr(self, '_exposure_info'):
                        self._exposure_info_n = 1
                        self._exposure_info_gt_mean = gt_mean
                        self._exposure_info_gt_mean2 = gt_mean2
                    else:
                        self._exposure_info_n += 1
                        self._exposure_info_gt_mean = torch.lerp(
                            self._exposure_info_gt_mean, gt_mean,
                            1.0 / self._exposure_info_n
                        )
                        self._exposure_info_gt_mean2 = torch.lerp(
                            self._exposure_info_gt_mean2, gt_mean2,
                            1.0 / self._exposure_info_n
                        )
                splat_colors = self.features_dc
                if self.config.sh_degree > 0:
                    splat_colors = SH2RGB(splat_colors)
                splat_colors = torch.nan_to_num(
                    splat_colors, nan=0.0, posinf=0.0, neginf=0.0)
                splat_weights = torch.exp(self.scales.mean(-1, True)) * torch.sigmoid(self.opacities)
                splat_color_mean = (splat_colors*splat_weights).sum(0) / splat_weights.sum(0)
                splat_color_mean2 = torch.cov(splat_colors.T, aweights=splat_weights.flatten())
                diff1 = torch.abs(splat_color_mean - self._exposure_info_gt_mean).mean()
                diff2 = torch.abs(splat_color_mean2 - self._exposure_info_gt_mean2).mean()
                exposure_reg = self.config.exposure_reg_splat * (diff1+diff2)
                # print(splat_color_mean.detach(), diff1.item(), diff2.item())
                timerl.mark("exposure_reg")

        # image loss
        if False:
            pred_img_e = saturate_keep_gradient(pred_img_e, 0.0, 1.0)
            pred_img = saturate_keep_gradient(pred_img, 0.0, 1.0)
        else:
            pred_img_e = torch.clip(pred_img_e, 0.0, 1.0)
            pred_img = torch.clip(pred_img, 0.0, 1.0)
        Ll1_e = torch.abs(gt_img - pred_img_e).mean()
        Ll1 = torch.abs(gt_img - pred_img).mean()
        simloss = 0.0
        if self.config.ssim_lambda > 0.0:
            simloss = 1 - self.ssim(gt_img.permute(2, 0, 1)[None, ...], pred_img_e.permute(2, 0, 1)[None, ...])
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
        quat_norm_reg = 0.1 * (torch.log(self.quats.norm(dim=-1)+0.001)**2).mean()
        timerl.mark("mcmc")  # ~100us

        loss_dict = {
            "main_loss": torch.lerp(torch.lerp(Ll1_e, simloss, self.config.ssim_lambda), Ll1, self.config.exposure_reg_image),
            "alpha_loss": alpha_loss,
            "scale_reg": scale_reg,
            "depth_reg": depth_reg,
            "normal_reg": normal_reg,
            "quat_reg": quat_norm_reg,
            'mcmc_opacity_reg': mcmc_opacity_reg,
            'mcmc_scale_reg': mcmc_scale_reg,
            "exposure_reg": exposure_reg,
            "exposure_param_reg": self.config.exposure_reg_param * exposure_param_reg,
        }
        for key, value in loss_dict.items():
            loss_dict[key] = self.config.loss_scale * value
        timerl.mark("scale")  # <100us

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
