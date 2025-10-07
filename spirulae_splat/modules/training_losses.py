import torch
import math

from spirulae_splat.modules.supervision import SupervisionLosses
from spirulae_splat.modules.exposure_correction import ExposureCorrection

from spirulae_splat.splat._camera import _Camera
from spirulae_splat.splat.utils import resize_image, _TORCH_COMPILE_ARGS

from spirulae_splat.splat.cuda import _C

from fused_ssim import fused_ssim

from fused_bilagrid import BilateralGrid, slice, total_variation_loss


class _MaskGradient(torch.autograd.Function):
    @staticmethod
    def forward(ctx, x, mask):
        ctx.mask = mask
        return x
    @staticmethod
    def backward(ctx, v_x):
        return v_x * ctx.mask, None


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
        quat_norm_reg_weight: float
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

        losses = _C.compute_per_splat_losses_forward(
            scales, opacities, quats,
            *hyperparams
        )

        ctx.hyperparams = hyperparams
        ctx.save_for_backward(scales, opacities, quats)

        return losses

    @staticmethod
    def backward(ctx, v_losses):

        hyperparams = ctx.hyperparams
        scales, opacities, quats = ctx.saved_tensors

        v_inputs = _C.compute_per_splat_losses_backward(
            scales, opacities, quats,
            v_losses,
            *hyperparams
        )
        return (*v_inputs, *([None]*len(hyperparams)))



class SplatTrainingLosses(torch.nn.Module):

    def __init__(self, config, num_training_data):
        super().__init__()

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
            image: tensor.Tensor in type uint8 or float32
        """
        if image.dtype == torch.uint8:
            image = image.float() / 255.0
        gt_img = self._downscale_if_required(image)
        return gt_img

    def apply_bilateral_grid(self, rgb: torch.Tensor, cam_idx: int, H: int, W: int) -> torch.Tensor:
        """rgb must be clamped to 0-1"""
        out = slice(
            bil_grids=self.bil_grids, rgb=rgb, xy=None,
            grid_idx=torch.tensor(cam_idx, device=rgb.device, dtype=torch.long)[:,None],
        )
        return out["rgb"]

    # @torch.compile(**_TORCH_COMPILE_ARGS)
    def composite_with_background(self, image, background) -> torch.Tensor:
        """Composite the ground truth image with a background color when it has an alpha channel.

        Args:
            image: the image to composite
            background: the background color
        """
        if image.shape[-1] == 4:
            alpha = image[..., -1].unsqueeze(-1).repeat((1, 1, 1, 3))
            return alpha * image[..., :3] + (1 - alpha) * background
        else:
            return image

    def get_2dgs_reg_weights(self):
        weight_depth_reg = self.config.depth_reg_weight * \
            min(self.step / max(self.config.depth_reg_warmup, 1), 1)
        weight_normal_reg = self.config.normal_reg_weight * \
            min(self.step / max(self.config.normal_reg_warmup, 1), 1)
        return weight_depth_reg, weight_normal_reg

    def get_alpha_reg_weight(self):
        return self.config.alpha_reg_weight * \
            min(self.step / max(self.config.alpha_reg_warmup, 1), 1)

    # @torch.compile(**_TORCH_COMPILE_ARGS)
    def _image_loss_0(self, gt_img, pred_img, pred_img_e):
        pred_img_e = torch.clip(pred_img_e, 0.0, 1.0)
        pred_img = torch.clip(pred_img, 0.0, 1.0)
        Ll1_e = torch.abs(gt_img - pred_img_e).mean()
        Ll1 = torch.abs(gt_img - pred_img).mean()
        Ll2_e = ((gt_img-pred_img_e)**2).mean()

        gt_img_bchw = gt_img.permute(0, 3, 1, 2).contiguous()
        pred_img_bchw = pred_img_e.permute(0, 3, 1, 2).contiguous()
        return gt_img_bchw, pred_img_bchw, Ll1_e, Ll1, Ll2_e

    # @torch.compile(dynamic=False)
    def _image_loss_1(self, Ll1_e, Ll1, ssim, ssim_lambda, exposure_reg_image):
        return torch.lerp(torch.lerp(Ll1_e, 1-ssim, ssim_lambda), Ll1, exposure_reg_image)

    def image_loss(self, gt_img, pred_img, pred_img_e, exposure_reg_image):
        gt_img_bchw, pred_img_bchw, Ll1_e, Ll1, Ll2_e = self._image_loss_0(gt_img, pred_img, pred_img_e)
        # ssim = fused_ssim(pred_img_bchw, gt_img_bchw, padding="valid")
        ssim = fused_ssim(pred_img_bchw, gt_img_bchw, padding="same")

        ssim_lambda = self.config.ssim_lambda * min(self.step/max(self.config.ssim_warmup,1), 1)
        return self._image_loss_1(Ll1_e, Ll1, ssim, ssim_lambda, exposure_reg_image), Ll2_e, ssim

    # @torch.compile(**_TORCH_COMPILE_ARGS)
    def alpha_reg(self, alpha):
        weight_alpha_reg = self.get_alpha_reg_weight()
        if self.config.randomize_background:
            reg_alpha = 1.0 - alpha**2  # push to 1
        else:
            # reg_alpha = torch.log(4.0*torch.clip(alpha*(1.0-alpha), min=1e-2))
            reg_alpha = 4.0*alpha*(1.0-alpha)
        return weight_alpha_reg * reg_alpha.mean()

    def forward(self, step: int, batch, outputs):
        self.step = step

        # mask out of bound (e.g. fisheye circle)
        camera_mask = None
        # TODO
        # ssplat_camera = outputs["ssplat_camera"]  # type: _Camera
        # if ssplat_camera.is_distorted():
        #     undist_map = ssplat_camera.get_undist_map()
        #     camera_mask = torch.isfinite(undist_map.sum(-1, True))
        #     if not camera_mask.all():
        #         for key in ['rgb', 'depth', 'alpha', 'background']:
        #             if key in outputs:
        #                 outputs[key] = _MaskGradient.apply(outputs[key], camera_mask)
        device = outputs['rgb'].device

        gt_img_rgba = self.get_gt_img(batch["image"].to(device))
        gt_img = self.composite_with_background(gt_img_rgba, outputs["background"])
        pred_img = outputs["rgb"]

        # alpha channel for bounded objects - apply a cost on rendered alpha
        alpha_loss = 0.0
        if gt_img_rgba.shape[-1] == 4:
            alpha = gt_img_rgba[..., -1].unsqueeze(-1)
            alpha_loss = alpha_loss + SupervisionLosses.get_alpha_loss(outputs['alpha'], alpha)

        # separate mask for dynamic objects, text, etc.
        # simply don't consider it when evaluating loss, unless theres's no alpha channel, where a cost is applied
        mask = None
        if "mask" in batch:
            # batch["mask"] : [H, W, 1]
            mask = self._downscale_if_required(batch["mask"])
            mask = mask.float().to(gt_img.device)
            assert mask.shape[:-1] == gt_img.shape[:-1] == pred_img.shape[:-1]
            # can be little bit sketchy for the SSIM loss
            gt_img = torch.lerp(outputs["background"], gt_img, mask)
            # pred_img = torch.lerp(outputs["background"], pred_img, mask)
            if isinstance(alpha_loss, float) and alpha_loss == 0.0:
                alpha_loss = alpha_loss + SupervisionLosses.get_alpha_loss(outputs['alpha'], mask)

        alpha_loss = self.config.alpha_loss_weight * alpha_loss

        # depth supervision
        depth_supervision_loss, normal_supervision_loss, alpha_supervision_loss = 0.0, 0.0, 0.0
        if "depth" in batch and self.step > self.config.supervision_warmup and \
            (self.config.depth_supervision_weight > 0.0 or
             self.config.normal_supervision_weight > 0.0 or
             self.config.alpha_supervision_weight > 0.0 or
             self.config.alpha_supervision_weight_under > 0.0):
            if batch["depth"].ndim == 2:
                batch["depth"] = batch["depth"].unsqueeze(-1)
            depth = self._downscale_if_required(batch["depth"].to(device))
            depth_supervision_loss, normal_supervision_loss, alpha_supervision_loss \
                = self.supervision_losses(
                    depth, ssplat_camera, outputs["depth"],
                    outputs["depth_normal"] if self.config.compute_depth_normal else None,
                    outputs["alpha"]
                )

        # correct exposure
        pred_img_e = pred_img
        exposure_param_reg = 0.0
        exposure_reg_image = 0.0
        if self.config.adaptive_exposure_mode is not None and \
            self.step > self.config.adaptive_exposure_warmup:
            exposure_reg_image = self.config.exposure_reg_image
            # mask, note that alpha_mask is defined in "mask out of bound" step
            alpha_mask = camera_mask
            if mask is not None and alpha_mask is not None:
                alpha_mask = alpha_mask & mask
            # call function
            pred_img_e, exposure_param_reg = self.exposure_correction(pred_img, gt_img, alpha_mask)

        # image loss
        image_loss, mse, ssim = self.image_loss(gt_img, pred_img, pred_img_e, exposure_reg_image)

        # depth and normal regularizers
        depth_reg, normal_reg = 0.0, 0.0
        if not self.config.use_3dgs and self.step >= self.config.reg_warmup_length:
            reg_depth = outputs["reg_depth"]
            reg_normal = outputs["reg_normal"]
            weight_depth_reg, weight_normal_reg = self.get_2dgs_reg_weights()
            depth_reg = weight_depth_reg * reg_depth.mean()
            normal_reg = weight_normal_reg * reg_normal[1:-1, 1:-1].mean()

        # alpha regularizer
        alpha_reg = 0.0
        if self.step >= self.config.reg_warmup_length:
            alpha = outputs['alpha']
            alpha_reg = self.alpha_reg(alpha)

        # metrics, readable from console during training
        with torch.no_grad():
            if not hasattr(self, '_running_metrics'):
                self._running_metrics = { 'psnr': [], 'ssim': [] }
            psnr_list = self._running_metrics['psnr']
            ssim_list = self._running_metrics['ssim']
            psnr = -10.0 * math.log10(mse.item())
            ssim = ssim.item()
            psnr_list.append(psnr)
            ssim_list.append(ssim)
            if len(psnr_list) > self.num_train_data:
                del psnr_list[0]
                del ssim_list[0]
            psnr = sum(psnr_list) / len(psnr_list)
            ssim = sum(ssim_list) / len(ssim_list)

        loss_dict = {
            # [C] RGB and alpha
            "image_loss": image_loss,
            "alpha_loss": alpha_loss,
            "psnr": float(psnr),
            "ssim": float(ssim),
            # [S] supervision
            "depth_ref_loss": depth_supervision_loss,
            "normal_ref_loss": normal_supervision_loss,
            "alpha_ref_loss": alpha_supervision_loss,
            # [G] 2DGS
            "depth_reg": depth_reg,
            "normal_reg": normal_reg,
            "alpha_reg": alpha_reg,
            # [E] exposure
            "tv_loss": 0.0,  # see get_per_splat_losses()
            "exposure_param_reg": exposure_param_reg,
        }

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
        bilagrid_tv_loss = 0.0
        if self.config.use_bilateral_grid:
            bilagrid_tv_loss = 10 * total_variation_loss(self.bil_grids.grids)
        loss_dict["tv_loss"] = bilagrid_tv_loss

        if not _use_torch_impl:
            losses = _ComputePerSplatLosses.apply(
                gauss_scales, gauss_opacities, gauss_quats,
                self.config.mcmc_opacity_reg * float(self.config.use_mcmc),
                self.config.mcmc_scale_reg * float(self.config.use_mcmc),
                self.config.max_gauss_ratio,
                self.config.scale_regularization_weight,
                self.config.erank_reg * float(self.step >= self.config.erank_reg_warmup),
                self.config.erank_reg_s3 * float(self.step >= self.config.erank_reg_warmup),
                0.01
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
