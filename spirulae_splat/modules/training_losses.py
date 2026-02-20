import torch
import torch.nn.functional as F
import math
import random

from spirulae_splat.modules.supervision import SupervisionLosses
from spirulae_splat.modules.exposure_correction import ExposureCorrection

from spirulae_splat.splat._camera import _Camera
from spirulae_splat.splat.utils import resize_image, _TORCH_COMPILE_ARGS

from spirulae_splat.splat.cuda import (
    _C,
    _make_lazy_cuda_func,
    ray_depth_to_linear_depth
)

from spirulae_splat.splat.cuda._wrapper_per_pixel import (
    depth_to_normal,
    linear_rgb_to_srgb,
    get_color_transform_matrix,
    apply_ppisp
)

from spirulae_splat.splat.cuda._wrapper_projection import add_gradient_component

from fused_bilagrid import BilateralGrid, fused_bilagrid_sample, total_variation_loss

from typing import List, Optional, Literal

import spirulae_splat.modules.enhancer

# from torchmetrics.image.lpip import LearnedPerceptualImagePatchSimilarity
from spirulae_splat.modules.lpips import LearnedPerceptualImagePatchSimilarity

from nerfstudio.cameras.cameras import Cameras, CameraType


class _BackwardMetadataInjector(torch.autograd.Function):
    @staticmethod
    def forward(ctx, x):
        ctx.meta = {}
        def set_meta(key, meta):
            nonlocal ctx
            ctx.meta[key] = meta
        return x, set_meta
    @staticmethod
    def backward(ctx, v_x, *kwargs):
        for key, value in ctx.meta.items():
            setattr(v_x, key, value)
        del ctx.meta
        return v_x


class _MaskGradient(torch.autograd.Function):
    @staticmethod
    def forward(ctx, x, mask):
        ctx.mask = mask
        return x
    @staticmethod
    def backward(ctx, v_x):
        return v_x * ctx.mask, None


class FusedSSIM(torch.autograd.Function):
    @staticmethod
    def forward(ctx, img1, img2, train=True, return_ssim_map=False):
        ssim, ssim_map, dm_dmu1, dm_dsigma1_sq, dm_dsigma12 = \
            _make_lazy_cuda_func("fused_ssim_forward")(img1, img2, train, return_ssim_map)

        ctx.save_for_backward(img1.detach(), img2, dm_dmu1, dm_dsigma1_sq, dm_dsigma12)

        return ssim, ssim_map

    @staticmethod
    def backward(ctx, dL_dssim, dL_dmap):
        img1, img2, dm_dmu1, dm_dsigma1_sq, dm_dsigma12 = ctx.saved_tensors
        grad = _make_lazy_cuda_func("fused_ssim_backward")(img1, img2, dL_dssim, dm_dmu1, dm_dsigma1_sq, dm_dsigma12)
        return grad, None, None, None


class Dct3D(torch.autograd.Function):
    @staticmethod
    def forward(ctx, data: torch.Tensor, norm: Literal[None, "ortho"]="ortho", type: int=1):
        if norm != "ortho" or type != 1:
            raise NotImplementedError()
        return _make_lazy_cuda_func("dct3d_type1_ortho")(data)

    @staticmethod
    def backward(ctx, v_data: torch.Tensor):
        return _make_lazy_cuda_func("dct3d_type1_ortho")(v_data)


class _ComputePerSplatLosses(torch.autograd.Function):

    @staticmethod
    def forward(
        ctx,
        scales, opacities, quats,
        mcmc_opacity_reg: float,
        mcmc_scale_reg: float,
        max_gauss_ratio: float,
        scale_regularization_weight: float,
        erank_reg: float,
        erank_reg_s3: float,
        quat_norm_reg_weight: float,
        compute_hessian_diagonal: Literal[None, "position", "all"] = None
    ):

        hyperparams = (
            mcmc_opacity_reg,
            mcmc_scale_reg,
            max_gauss_ratio,
            scale_regularization_weight,
            erank_reg,
            erank_reg_s3,
            quat_norm_reg_weight
        )

        losses = _make_lazy_cuda_func("compute_per_splat_losses_forward")(
            scales, opacities, quats,
            *hyperparams
        )

        ctx.hyperparams = hyperparams
        ctx.save_for_backward(scales, opacities, quats)
        ctx.compute_hessian_diagonal = (compute_hessian_diagonal == "all")

        return losses

    @staticmethod
    def backward(ctx, v_losses):

        hyperparams = ctx.hyperparams
        scales, opacities, quats = ctx.saved_tensors

        if ctx.compute_hessian_diagonal:
            v_inputs, vr_inputs, h_inputs = \
            _make_lazy_cuda_func("compute_per_splat_losses_backward_with_hessian_diagonal")(
                scales, opacities, quats,
                v_losses,
                *hyperparams
            )
            for v, vr, h in zip(v_inputs, vr_inputs, h_inputs):
                add_gradient_component(v, 'gradr', vr)
                add_gradient_component(v, 'hess', h)
                # print(v.shape, torch.isfinite(h).float().mean().item(), torch.nan_to_num(v, 0, 0, 0).mean().item(), torch.nan_to_num(vr, 0, 0, 0).mean().item(), torch.nan_to_num(h, 0, 0, 0).mean().item())
        else:
            v_inputs = _make_lazy_cuda_func("compute_per_splat_losses_backward")(
                scales, opacities, quats,
                v_losses,
                *hyperparams
            )
        return (*v_inputs, *([None]*len(hyperparams)), None)


class _ComputePerPixelLosses(torch.autograd.Function):

    @staticmethod
    def forward(
        ctx,
        render_rgb: Optional[torch.Tensor],
        ref_rgb: Optional[torch.Tensor],
        render_depth: Optional[torch.Tensor],
        ref_depth: Optional[torch.Tensor],
        render_normal: Optional[torch.Tensor],
        depth_normal: Optional[torch.Tensor],
        ref_normal: Optional[torch.Tensor],
        render_alpha: Optional[torch.Tensor],
        rgb_dist: Optional[torch.Tensor],
        depth_dist: Optional[torch.Tensor],
        normal_dist: Optional[torch.Tensor],
        ref_alpha: Optional[torch.Tensor],
        mask: Optional[torch.Tensor],
        depth_mask: Optional[torch.Tensor],
        normal_mask: Optional[torch.Tensor],
        alpha_mask: Optional[torch.Tensor],
        weights: List[float],
        num_train_images: int = -1,
        camera_indices: Optional[torch.Tensor] = None,
        return_loss_map: bool = False
    ):
        if not isinstance(camera_indices, torch.Tensor):
            num_train_images = -1
            camera_indices = None
        
        tensors = (
            render_rgb,
            ref_rgb,
            render_depth,
            ref_depth,
            render_normal,
            depth_normal,
            ref_normal,
            render_alpha,
            rgb_dist,
            depth_dist,
            normal_dist,
            ref_alpha,
            mask,
            depth_mask,
            normal_mask,
            alpha_mask
        )

        losses, raw_losses, loss_map = _C.compute_per_pixel_losses_forward(
            *tensors, weights,
            num_train_images, camera_indices,
            return_loss_map
        )
        # print(losses)
        # print(raw_losses[0].detach().cpu().numpy().tolist())
        # print(raw_losses[1].detach().cpu().numpy().tolist())

        ctx.weights = weights
        ctx.num_train_images = num_train_images
        ctx.save_for_backward(*tensors, raw_losses, camera_indices)

        # print('loss_map:', loss_map.mean().item(), loss_map.median().item())

        return losses, loss_map

    @staticmethod
    def backward(ctx, v_losses, v_loss_map):
        # warn that loss map is not differentiable
        grads = _C.compute_per_pixel_losses_backward(
            *ctx.saved_tensors[:-1],
            ctx.weights,
            v_losses,
            ctx.num_train_images,
            ctx.saved_tensors[-1],
        )
        return *grads, *([None]*(len(ctx.needs_input_grad)-len(grads)))


class _ComputePPISPRegularization(torch.autograd.Function):
    @staticmethod
    def forward(
        ctx,
        param_type: str,
        ppisp_params: torch.Tensor,  # [B, PPISP_NUM_PARAMS]
        loss_weights: List[float],
    ):
        losses, raw_losses = _C.compute_ppsip_regularization_forward(
            ppisp_params,
            loss_weights,
            param_type
        )

        ctx.loss_weights = loss_weights
        ctx.save_for_backward(ppisp_params, raw_losses)
        ctx.param_type = param_type

        return losses

    @staticmethod
    def backward(ctx, v_losses):
        ppisp_params, raw_losses = ctx.saved_tensors

        v_ppisp_params = _C.compute_ppsip_regularization_backward(
            ppisp_params,
            ctx.loss_weights,
            raw_losses,
            v_losses.contiguous(),
            ctx.param_type
        )

        return None, v_ppisp_params, None


DEFAULT_BILAGRID_PARAMS = torch.Tensor([
    1, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, 1, 0
]).float().cuda()

DEFAULT_PPISP_PARAMS = [
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.013658988289535046, 0.013658988289535046, 0.37816452980041504, 0.0, 0.013658988289535046, 0.013658988289535046, 0.37816452980041504, 0.0, 0.013658988289535046, 0.013658988289535046, 0.37816452980041504, 0.0]
]

DEFAULT_PPISP_PARAMS_RQS = [
    [0.0] * 39
]

class SplatTrainingLosses(torch.nn.Module):

    def __init__(self, config, num_training_data):
        super().__init__()

        self.step = 0
        self.config = config
        self.num_train_data = num_training_data

        self.supervision_losses = SupervisionLosses(config)
        self.exposure_correction = ExposureCorrection(config)

        if self.config.use_bilateral_grid:
            self.bil_grids = BilateralGrid(
                num=self.num_train_data,
                grid_X=self.config.bilagrid_shape[0],
                grid_Y=self.config.bilagrid_shape[1],
                grid_W=self.config.bilagrid_shape[2],
            )
            if self.config.optimize_bilagrid_frequencies:
                self.bil_grids.grids.data = \
                    Dct3D.apply(self.bil_grids.grids.data.cuda()).to(self.bil_grids.grids.data.device)
        if self.config.use_bilateral_grid_for_geometry:
            self.bil_grids_depth = BilateralGrid(
                num=self.num_train_data,
                grid_X=self.config.bilagrid_shape_geometry[0],
                grid_Y=self.config.bilagrid_shape_geometry[1],
                grid_W=self.config.bilagrid_shape_geometry[2],
            )
            self.bil_grids_normal = BilateralGrid(
                num=self.num_train_data,
                grid_X=self.config.bilagrid_shape_geometry[0],
                grid_Y=self.config.bilagrid_shape_geometry[1],
                grid_W=self.config.bilagrid_shape_geometry[2],
            )
        if self.config.use_ppisp:
            self.ppisp_params = torch.Tensor(
                DEFAULT_PPISP_PARAMS_RQS if self.config.ppisp_param_type == "rqs" else DEFAULT_PPISP_PARAMS
            ).repeat(self.num_train_data, 1)
            self.ppisp_params = torch.nn.Parameter(self.ppisp_params)

        self.lpips_dtype = torch.bfloat16 if torch.cuda.is_bf16_supported() else torch.float32
        # self.lpips_dtype = torch.float32
        if self.config.lpips_lambda > 0.0:
            self.lpips = LearnedPerceptualImagePatchSimilarity(
                net_type="vgg", normalize=True
                # net_type="alex", normalize=True
            ).to(self.lpips_dtype)

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

    def get_gt_img(self, image: torch.Tensor):
        """Compute groundtruth image with iteration dependent downscale factor for evaluation purpose

        Args:
            image: tensor.Tensor in type uint8, uint16, or float32
        """
        if image.dtype == torch.uint8:
            image = image.float() / 255.0
        elif image.dtype == torch.uint16:
            image = image.float() / 65535.0
        gt_img = self._downscale_if_required(image)
        return gt_img

    def apply_bilateral_grid(self, bilagrid, rgb: torch.Tensor, cam_idx: int, **kwargs) -> torch.Tensor:
        """rgb must be clamped to 0-1"""
        try:
            grid_idx = torch.tensor(cam_idx, device=rgb.device, dtype=torch.long).flatten()
            grids = bilagrid.grids[grid_idx]
            if self.config.optimize_bilagrid_frequencies:
                grids = Dct3D.apply(grids)
            out = fused_bilagrid_sample(
                grids, coords=None, rgb=rgb.unsqueeze(1),
                actual_width=kwargs.get('actual_width', None),
                actual_height=kwargs.get('actual_height', None),
                patch_offsets=kwargs.get('patch_offsets', None),
            ).squeeze(1)
        except TypeError:
            raise RuntimeError(
                "\033[93m"
                "You are likely using an outdated version of fused_bilagrid. "
                "Please update to latest fused_bilagrid from https://github.com/harry7557558/fused-bilagrid"
                "\033[0m"
            )
        return out

    @staticmethod
    def get_visibility_masks(batch, device=torch.device("cuda")):
        masks = None

        if "mask" in batch:
            batch_mask = batch['mask'].to(device)
            if len(batch_mask.shape) == 3:
                batch_mask = batch_mask.unsqueeze(-1)
            masks = batch_mask

        if 'depth' in batch:
            gt_depth = batch['depth'].to(device)
            if len(gt_depth.shape) == 3:
                gt_depth = gt_depth.unsqueeze(-1)

            none_sky_mask = gt_depth < torch.amax(gt_depth, dim=(1,2,3), keepdims=True).detach()
            masks = masks & none_sky_mask if masks is not None else none_sky_mask

        return masks

    def get_2dgs_reg_weights(self):
        factor = min(self.step / max(self.config.distortion_reg_warmup, 1), 1)
        weight_depth_reg = self.config.depth_reg_weight * factor
        weight_normal_dist_reg = self.config.normal_distortion_reg_weight * factor
        weight_rgb_dist_reg = self.config.rgb_distortion_reg_weight * factor
        weight_normal_reg = self.config.normal_reg_weight * factor
        return (weight_depth_reg, weight_normal_dist_reg, weight_rgb_dist_reg), weight_normal_reg

    def get_alpha_reg_weight(self):
        return self.config.alpha_reg_weight * \
            min(self.step / max(self.config.alpha_reg_warmup, 1), 1)

    def forward(self, step: int, batch, outputs, meta={}, val=False):
        self.step = step

        device = outputs['rgb'].device
        camera = outputs["camera"]
        intrins = torch.concatenate((camera.fx, camera.fy, camera.cx, camera.cy), dim=1)

        # If reference image is empty, AI generate from rendered image
        if "image" not in batch:
            image_shape = outputs["rgb"].shape
            batch["image"] = spirulae_splat.modules.enhancer.infer(outputs["rgb"].detach()).detach()
            # batch["image"] = outputs["rgb"].detach()
            B, H, W, C = batch["image"].shape
            for key, value in outputs.items():
                if isinstance(value, torch.Tensor) and value.shape[:3] == image_shape[:3]:
                    outputs[key] = value[:, :H, :W, :]
            camera.height = H * torch.ones_like(camera.height)
            camera.width = W * torch.ones_like(camera.width)
            if random.random() < 0.1:
                from PIL import Image
                Image.fromarray((255*torch.clip(batch["image"].detach()[0,:,:,:3],0,1)).byte().cpu().numpy()).save("/mnt/d/temp.png")

        # Load rendering outputs
        pred_rgb = outputs["rgb"]
        pred_depth = outputs["depth"] if 'depth' in outputs else None
        pred_normal = outputs["normal"] if 'normal' in outputs else None
        pred_depth_normal = outputs["depth_normal"] if 'depth_normal' in outputs else None
        pred_alpha = outputs["alpha"] if 'alpha' in outputs else None

        gt_rgb, gt_depth, gt_normal, gt_alpha = None, None, None, None  # for loss
        gt_rgb_mask, gt_depth_mask, gt_normal_mask, gt_alpha_mask = None, None, None, None  # for masking
        none_sky_mask = None

        # load alpha
        if "mask" in batch:
            batch_mask = self._downscale_if_required(batch['mask'].to(device).float()) > 0.5
            gt_rgb_mask = batch_mask
            if self.config.apply_loss_for_mask:
                gt_alpha = batch_mask

        # load depth
        if 'depth' in batch:
            gt_depth = batch['depth'].to(device)
            if len(gt_depth.shape) == 3:
                gt_depth = gt_depth.unsqueeze(-1)
            gt_depth = self._downscale_if_required(gt_depth)
            gt_depth_mask = (gt_depth != 0.0)

            # sky mask
            none_sky_mask = gt_depth < torch.amax(gt_depth, dim=(1,2,3), keepdims=True).detach()
            gt_alpha_mask = gt_depth_mask

            # apply bilagrid
            if self.config.use_bilateral_grid_for_geometry and \
                    (camera.metadata is not None and "cam_idx" in camera.metadata):
                # TODO: fused kernel
                # TODO: might not be the best way to use RGB bilagrid
                B, H, W, C = gt_depth.shape
                gt_depth = gt_depth * (
                    gt_depth_mask.float().sum(dim=(1, 2, 3), keepdims=True) /
                    (gt_depth * gt_depth_mask.float()).sum(dim=(1, 2, 3), keepdims=True)  # TODO: fix zero division
                )
                gt_depth = gt_depth / (gt_depth + 1.0)
                gt_depth = self.apply_bilateral_grid(
                    self.bil_grids_depth,
                    gt_depth.repeat(1, 1, 1, 3), camera.metadata["cam_idx"],
                    **meta
                )[..., :1]
                gt_depth = gt_depth / (1.0 - gt_depth).clip(max=0.999)

        # load normal
        if 'normal' in batch:
            gt_normal = batch['normal']
            if gt_normal.dtype == torch.uint8:
                gt_normal = gt_normal.float() / (255/2) - 1.0
            gt_normal = self._downscale_if_required(gt_normal.to(device))
            gt_normal_mask = (gt_normal.sum(-1, True) > -2.366)  # background is (-1, -1, -1)

            # apply bilagrid
            if self.config.use_bilateral_grid_for_geometry and \
                    (camera.metadata is not None and "cam_idx" in camera.metadata):
                # TODO: fused kernel
                # TODO: might not be the best way to use RGB bilagrid
                B, H, W, C = gt_normal.shape
                gt_normal = self.apply_bilateral_grid(
                    self.bil_grids_normal,
                    0.5+0.5*gt_normal, camera.metadata["cam_idx"],
                    **meta
                ) * 2.0 - 1.0

        # mask sky
        if none_sky_mask is not None:
            # apply to depth mask (already there if sky mask is there)
            gt_depth_mask = gt_depth_mask & none_sky_mask
            # apply to normal desk
            if gt_normal_mask is None:
                gt_normal_mask = none_sky_mask
            else:
                gt_normal_mask = gt_normal_mask & none_sky_mask
            if not self.config.enable_sky_masking:
                none_sky_mask = None
        if none_sky_mask is not None:
            # apply loss to discourage opacity
            if gt_alpha is not None:
                gt_alpha = gt_alpha & none_sky_mask
            else:
                gt_alpha = none_sky_mask
            # apply to rgb mask
            if not (self.config.apply_loss_for_mask or self.config.train_background_color):
                if gt_rgb_mask is not None:
                    gt_rgb_mask = gt_rgb_mask & none_sky_mask
                else:
                    gt_rgb_mask = none_sky_mask

        # load RGB
        if self.config.fit == "rgb":
            gt_img_rgba = self.get_gt_img(batch["image"].to(device))
        elif self.config.fit == "depth":
            gt_img_rgba = gt_depth
        elif self.config.fit in ["normal"]:
            gt_img_rgba = 0.5+0.5*F.normalize(gt_normal, dim=-1)
        gt_rgb = gt_img_rgba[..., :3]

        # update alpha if image is RGBA
        if gt_img_rgba.shape[-1] == 4 and self.config.alpha_loss_weight > 0.0:
            alpha = gt_img_rgba[..., -1].unsqueeze(-1)
            gt_rgb_mask = gt_rgb_mask & alpha if gt_rgb_mask is not None else alpha
            if self.config.apply_loss_for_mask:
                gt_alpha = gt_alpha & alpha if gt_alpha is not None else alpha

        # replace parts of background with random noise to discourage transparency
        background = outputs["background"]
        if self.config.randomize_background == "opaque-only":
            background_mask = gt_rgb_mask
            if none_sky_mask is not None:
                background_mask = none_sky_mask if background_mask is None else \
                    background_mask & none_sky_mask
            if background_mask is not None:
                background = torch.where(background_mask, torch.rand_like(background), background)
        if self.config.randomize_background == "non-sky-only":
            if none_sky_mask is not None:
                background = torch.where(none_sky_mask, torch.rand_like(background), background)

        # do this to make SSIM/LPIPS happier + encourage sky transparency
        if gt_rgb_mask is not None and False:
            gt_rgb = torch.where(gt_rgb_mask, gt_rgb, background)
        if gt_rgb_mask is not None or none_sky_mask is not None:
            background_mask = gt_rgb_mask
            if none_sky_mask is not None and ('max_blending' in meta or random.random() < 0.1):
                background_mask = none_sky_mask if background_mask is None else \
                    background_mask & none_sky_mask
            if background_mask is not None:
                pred_rgb = torch.where(background_mask, pred_rgb, background)

        # convert linear RGB to sRGB if needed
        if self.config.image_color_space != None:
            color_transform = get_color_transform_matrix(self.config.image_color_space)
            pred_rgb = torch.matmul(pred_rgb, color_transform.T)
            gt_rgb = torch.matmul(gt_rgb, color_transform.T)
            # don't clip; exposure correction and loss should ideally handle out-of-gamut colors
        if self.config.use_linear_color_space:
            pred_rgb = linear_rgb_to_srgb(pred_rgb)
            gt_rgb = linear_rgb_to_srgb(gt_rgb)

        # apply exposure correction
        if self.config.use_bilateral_grid and self.config.fit == "rgb" and \
                camera.metadata is not None and "cam_idx" in camera.metadata:
            pred_rgb = self.apply_bilateral_grid(
                self.bil_grids,
                pred_rgb, camera.metadata["cam_idx"],
                **meta
            )
        if self.config.use_ppisp and self.config.fit == "rgb" and \
                camera.metadata is not None and "cam_idx" in camera.metadata:
            indices = torch.tensor(camera.metadata["cam_idx"]).flatten().to(device)
            ppisp_param = self.ppisp_params[indices, :]
            pred_rgb = apply_ppisp(
                self.config.ppisp_param_type,
                pred_rgb, ppisp_param,
                intrins=intrins,
                actual_image_width=camera.metadata.get('actual_width', None),
                actual_image_height=camera.metadata.get('actual_height', None),
            )

        # correct exposure (deprecated)
        if self.config.adaptive_exposure_mode is not None and \
            self.step > self.config.adaptive_exposure_warmup:
            raise NotImplementedError("Adaptive exposure is deprecated. Use bilateral grid instead.")

        # ssim loss
        ssim, ssim_map = FusedSSIM.apply(
            pred_rgb.contiguous(),
            gt_rgb.contiguous(),
            True,  # TODO: detect torch.no_grad()
            self.config.compute_hessian_diagonal is not None
        )

        # call fused kernel to compute loss
    
        (weight_depth_dist_reg, weight_normal_dist_reg, weight_rgb_dist_reg), weight_normal_reg = \
            self.get_2dgs_reg_weights()

        if pred_depth_normal is None and pred_depth is not None and (pred_normal is not None or gt_normal is not None):
            is_fisheye = (camera.camera_type[0].item() == CameraType.FISHEYE.value)
            pred_depth_normal = depth_to_normal(
                pred_depth, ["pinhole", "fisheye"][is_fisheye],
                intrins, camera.distortion_params
            )

        losses, loss_map = _ComputePerPixelLosses.apply(
            pred_rgb,
            gt_rgb,
            pred_depth,
            gt_depth,
            pred_normal,
            pred_depth_normal,
            gt_normal,
            pred_alpha,
            outputs['rgb_distortion'] if 'rgb_distortion' in outputs else None,
            outputs['depth_distortion'] if 'depth_distortion' in outputs else None,
            outputs['normal_distortion'] if 'normal_distortion' in outputs else None,
            gt_alpha,
            gt_rgb_mask,
            gt_depth_mask,
            gt_normal_mask,
            gt_alpha_mask,
            [
                # RGB supervision
                (1.0 - self.config.ssim_lambda) * (1.0 - self.config.lpips_lambda),
                # depth supervison
                float(self.step > self.config.supervision_warmup) *
                    self.config.depth_supervision_weight,
                # normal supervision
                float(self.step > self.config.supervision_warmup) *
                    self.config.normal_supervision_weight,
                # alpha supervision (over and under)
                self.config.alpha_loss_weight,
                0.0 if gt_alpha is None else self.config.alpha_loss_weight_under,
                # normal regularization
                float(self.step > self.config.reg_warmup_length) *
                    weight_normal_reg,
                # alpha regularization
                float(self.step >= self.config.reg_warmup_length) *
                    self.get_alpha_reg_weight(),
                # distortion regularizations (RGB, depth, normal)
                float(self.step >= self.config.reg_warmup_length) * weight_rgb_dist_reg,
                float(self.step >= self.config.reg_warmup_length) * weight_depth_dist_reg,
                float(self.step >= self.config.reg_warmup_length) * weight_normal_dist_reg,
            ],
            meta.get("num_train_data", -1),
            camera.metadata.get('cam_idx', None),
            self.config.compute_hessian_diagonal is not None
        )
        (
            rgb_l1, rgb_psnr,
            depth_supervision_loss, normal_supervision_loss, alpha_supervision_loss,
            normal_reg, alpha_reg,
            rgb_dist_reg, depth_dist_reg, normal_dist_reg
        ) = losses
        if loss_map is not None and ssim_map is not None:
            w_ssim = 1.0 - (1.0 - self.config.ssim_lambda) * (1.0 - self.config.lpips_lambda)
            loss_map += w_ssim * ssim_map
            if 'backward_metadata_injector' in outputs:
                outputs['backward_metadata_injector']("loss_map", loss_map)  # to be able to get it in backward
        # if loss_map is not None and self.step % 100 == 0:
        #     import matplotlib.pyplot as plt
        #     plt.imshow(loss_map[0].detach().cpu().numpy())
        #     plt.show()

        image_loss = (1.0 - self.config.ssim_lambda) * rgb_l1 + self.config.ssim_lambda * (1.0 - ssim)

        # LPIPS for training
        if self.config.lpips_lambda > 0.0:
            lpips = self.lpips(
                pred_rgb.permute(0, 3, 1, 2).to(self.lpips_dtype).clip(0, 1),
                gt_rgb.permute(0, 3, 1, 2).to(self.lpips_dtype).clip(0, 1)
            ).float()
            image_loss = torch.lerp(image_loss, lpips, self.config.lpips_lambda)

        # LPIPS for validation
        if val:
            if not hasattr(self, 'lpips_val'):
                self.lpips_val = LearnedPerceptualImagePatchSimilarity(
                    net_type="alex", normalize=True
                ).to(dtype=self.lpips_dtype, device=pred_rgb.device)
                # .to(memory_format=torch.channels_last)
            with torch.no_grad():
                # .contiguous(memory_format=torch.channels_last)
                lpips_val = self.lpips_val(
                    pred_rgb.permute(0, 3, 1, 2).to(self.lpips_dtype).clip(0, 1),
                    gt_rgb.permute(0, 3, 1, 2).to(self.lpips_dtype).clip(0, 1)
                ).float()

        # metrics, readable from console during training
        with torch.no_grad():
            # list_cap_max = self.num_train_data
            list_cap_max = self.config.refine_every
            if not hasattr(self, '_running_metrics'):
                self._running_metrics = { 'psnr': [], 'ssim': [], 'lpips': [] }
            psnr_list = self._running_metrics['psnr']
            ssim_list = self._running_metrics['ssim']
            psnr = rgb_psnr.item()
            ssim = ssim.item()
            psnr_list.append(psnr)
            ssim_list.append(ssim)
            if len(psnr_list) > list_cap_max:
                del psnr_list[0]
                del ssim_list[0]
            psnr = sum(psnr_list) / len(psnr_list)
            ssim = sum(ssim_list) / len(ssim_list)

            if self.config.lpips_lambda > 0.0:
                lpips_list = self._running_metrics['lpips']
                lpips = lpips.item()
                lpips_list.append(lpips)
                if len(lpips_list) > list_cap_max:
                    del lpips_list[0]
                lpips = sum(lpips_list) / len(lpips_list)

        loss_dict = {
            # [C] RGB and alpha
            "image_loss": image_loss,
            "psnr": float(psnr),
            "ssim": float(ssim),
            # [S] supervision
            "depth_ref_loss": depth_supervision_loss,
            "normal_ref_loss": normal_supervision_loss,
            "alpha_ref_loss": alpha_supervision_loss,
            # [G] 2DGS
            "normal_reg": normal_reg,
            "alpha_reg": alpha_reg,
            "depth_dist_reg": depth_dist_reg,
            "normal_dist_reg": normal_dist_reg,
            "rgb_dist_reg": rgb_dist_reg,
        }
        if self.config.lpips_lambda > 0.0:
            loss_dict['lpips'] = float(lpips)
        if val:
            loss_dict['lpips_val'] = float(lpips_val)

        return loss_dict

    # @torch.compile(**_TORCH_COMPILE_ARGS)
    def _erank_reg(self, gauss_scales):
        scales = torch.exp(2.0*gauss_scales)
        if self.config.use_3dgs:
            x, y, z = torch.unbind(scales, -1)
            s = x + y + z
            s1 = torch.fmax(torch.fmax(x, y), z) / s
            s3 = torch.fmin(torch.fmin(x, y), z) / s
            s2 = 1 - s1 - s3
            r = torch.exp(-s1*torch.log(s1) - s2*torch.log(s2) - s3*torch.log(s3))
            reg = torch.relu(-torch.log(r - 0.99999))
            return self.config.erank_reg * reg.mean() + \
                self.config.erank_reg_s3 * s3.mean()
        else:
            x, y = torch.unbind(scales, -1)
            s = x + y
            s1 = torch.fmax(x, y) / s
            s2 = torch.fmin(x, y) / s
            r = torch.exp(-s1*torch.log(s1) - s2*torch.log(s2))
            reg = torch.relu(-torch.log(r - 0.99999))
            return self.config.erank_reg * reg.mean()

    # @torch.compile(**_TORCH_COMPILE_ARGS)
    def _scale_reg(self, gauss_scales):
        scale_exp = torch.exp(gauss_scales)
        scale_reg = (
            torch.maximum(
                scale_exp.amax(dim=-1) / scale_exp.amin(dim=-1),
                torch.tensor(self.config.max_gauss_ratio),
            )
            - self.config.max_gauss_ratio
        )
        return self.config.scale_regularization_weight * scale_reg.mean()

    # @torch.compile(**_TORCH_COMPILE_ARGS)
    def _mcmc_opac_reg(self, gauss_opacities):
        mcmc_opacity_reg = torch.sigmoid(gauss_opacities).mean()
        return self.config.mcmc_opacity_reg * mcmc_opacity_reg

    # @torch.compile(**_TORCH_COMPILE_ARGS)
    def _mcmc_scale_reg(self, gauss_scales):
        mcmc_scale_reg = torch.exp(gauss_scales).mean()
        # mcmc_scale_reg = gauss_scales.mean()
        # mcmc_scale_reg = torch.where(gauss_scales < 0, torch.exp(gauss_scales), gauss_scales+1).mean()
        return self.config.mcmc_scale_reg * mcmc_scale_reg

    # @torch.compile(**_TORCH_COMPILE_ARGS)
    def _quat_norm_reg(self, gauss_quats):
        # quat_norm = gauss_quats.norm(dim=-1)
        quat_norm = torch.sqrt((gauss_quats**2).sum(dim=-1))
        quat_norm_reg = 0.01 * (quat_norm-1.0-torch.log(quat_norm)).mean()
        return quat_norm_reg

    def get_static_losses(self, step: int, gauss_quats, gauss_scales, gauss_opacities, loss_dict, _use_torch_impl=False):
        """Separately process losses that are not dependent on images"""
        self.step = step

        # bilagrid total variation loss
        if self.config.use_bilateral_grid:
            bilagrid = self.bil_grids.grids
            if self.config.optimize_bilagrid_frequencies:
                bilagrid = Dct3D.apply(bilagrid)
        if self.config.use_bilateral_grid:
            loss_dict["tv_loss"] = self.config.bilagrid_tv_loss_weight * total_variation_loss(bilagrid)
        if self.config.use_bilateral_grid_for_geometry:
            if self.config.depth_supervision_weight > 0.0:  # do this because bilagrid backward is expensive, especially in patched mode
                loss_dict["tv_loss_depth"] = self.config.bilagrid_tv_loss_weight_geometry * total_variation_loss(self.bil_grids_depth.grids)
            loss_dict["tv_loss_normal"] = self.config.bilagrid_tv_loss_weight_geometry * total_variation_loss(self.bil_grids_normal.grids)

        # bilagrid regularization loss
        if self.config.use_bilateral_grid:
            bilagrid_mean = torch.mean(bilagrid, dim=(0, 2, 3, 4))
            bilagrid_mean_reg = F.mse_loss(bilagrid_mean, DEFAULT_BILAGRID_PARAMS)
            loss_dict['bilagrid_mean_reg'] = self.config.bilagrid_mean_reg_weight * bilagrid_mean_reg

        # PPISP regularization loss
        if self.config.use_ppisp:
            ppisp_loss_weights = [
                self.config.ppisp_reg_exposure_mean,
                self.config.ppisp_reg_vig_center,
                self.config.ppisp_reg_vig_non_pos,
                self.config.ppisp_reg_vig_channel_var,
                self.config.ppisp_reg_color_mean,
                self.config.ppisp_reg_crf_channel_var,
            ]
            ppisp_reg_loss = _ComputePPISPRegularization.apply(
                self.config.ppisp_param_type,
                self.ppisp_params,
                ppisp_loss_weights
            )
            loss_dict['ppisp_reg_exposure_mean'] = ppisp_reg_loss[0]
            loss_dict['ppisp_reg_vig_center'] = ppisp_reg_loss[1]
            loss_dict['ppisp_reg_vig_non_pos'] = ppisp_reg_loss[2]
            loss_dict['ppisp_reg_vig_channel_var'] = ppisp_reg_loss[3]
            loss_dict['ppisp_reg_color_mean'] = ppisp_reg_loss[4]
            loss_dict['ppisp_reg_crf_channel_var'] = ppisp_reg_loss[5]

        if self.config.primitive == "voxel":
            return loss_dict

        if not _use_torch_impl:
            losses = _ComputePerSplatLosses.apply(
                gauss_scales, gauss_opacities, gauss_quats,
                self.config.mcmc_opacity_reg * float(self.config.use_mcmc),
                self.config.mcmc_scale_reg * float(self.config.use_mcmc),
                self.config.max_gauss_ratio,
                self.config.scale_regularization_weight,
                self.config.erank_reg * float(self.step >= self.config.erank_reg_warmup),
                self.config.erank_reg_s3 * float(self.step >= self.config.erank_reg_warmup),
                0.01,
                self.config.compute_hessian_diagonal
            )
            loss_dict['mcmc_opacity_reg'] = losses[0]
            loss_dict['mcmc_scale_reg'] = losses[1]
            loss_dict['scale_reg'] = losses[2]
            loss_dict['erank_reg'] = losses[3]
            loss_dict['quat_norm_reg'] = losses[4]
            return loss_dict

        # MCMC regularizers
        mcmc_opacity_reg, mcmc_scale_reg = 0.0, 0.0
        if self.config.use_mcmc:
            if self.config.mcmc_opacity_reg > 0.0:
                mcmc_opacity_reg = self._mcmc_opac_reg(gauss_opacities)
            if self.config.mcmc_scale_reg > 0.0:
                mcmc_scale_reg = self._mcmc_scale_reg(gauss_scales)
        loss_dict['mcmc_opacity_reg'] = mcmc_opacity_reg
        loss_dict['mcmc_scale_reg'] = mcmc_scale_reg

        # scale regularization
        scale_reg = 0.0
        if self.config.scale_regularization_weight > 0.0:
            scale_reg = self._scale_reg(gauss_scales)
        loss_dict['scale_reg'] = scale_reg

        # erank regularizers
        erank_reg = 0.0
        if self.step >= self.config.erank_reg_warmup and (
            self.config.erank_reg > 0.0 or self.config.erank_reg_s3 > 0.0
        ):
            erank_reg = self._erank_reg(gauss_scales)
        loss_dict['erank_reg'] = erank_reg

        # regularizations for parameters
        loss_dict['quat_norm_reg'] = self._quat_norm_reg(gauss_quats)

        return loss_dict
