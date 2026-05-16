import torch
import torch.nn.functional as F
import math
import random

# from spirulae_splat.modules.supervision import SupervisionLosses
# from spirulae_splat.modules.exposure_correction import ExposureCorrection

from spirulae_splat.splat._camera import _Camera
from spirulae_splat.splat.utils import resize_image, _TORCH_COMPILE_ARGS

from spirulae_splat.splat.cuda import (
    _C,
    _make_lazy_cuda_func,
    ray_depth_to_linear_depth,
    depth_normal_loss,
)

from spirulae_splat.splat.cuda._wrapper_per_pixel import (
    depth_to_normal,
    depth_normal_loss,
    rgb_to_srgb,
    get_color_transform_matrix,
    apply_ppisp
)
from spirulae_splat.modules.edge_detector import detect_edge

from spirulae_splat.splat.cuda._wrapper_projection import add_gradient_component

from spirulae_splat.modules.optimizer import OptimizerConfig

try:
    from fused_bilagrid import (
        BilateralGrid,
        BilateralGridPPISP,
        BilateralGridLoglinear,
        BilateralGridDepth,
        BilateralGridNormal,
        fused_bilagrid_sample,
        fused_bilagrid_ppisp_sample,
        fused_bilagrid_loglinear_sample,
        fused_bilagrid_depth_sample,
        fused_bilagrid_normal_sample,
        total_variation_loss,
        channel_mean
    )
    from fused_bilagrid import _C as fused_bilagrid_C
except:
    raise RuntimeError(
        "\033[93m"
        "You are likely using an incompatible version of fused_bilagrid. "
        "Please install latest fused_bilagrid from https://github.com/harry7557558/fused-bilagrid, branch `dev`."
        "\033[0m"
    )

from typing import List, Optional, Literal, Any, Union

import spirulae_splat.modules.enhancer

# from torchmetrics.image.lpip import LearnedPerceptualImagePatchSimilarity
from spirulae_splat.modules.lpips import LearnedPerceptualImagePatchSimilarity


class _MaskGradient(torch.autograd.Function):
    @staticmethod
    def forward(ctx, x, mask):
        ctx.mask = mask
        return x
    @staticmethod
    def backward(ctx, v_x):
        return v_x * ctx.mask, None


class _Dummy:
    def __init__(self):
        pass

class _SoftDetach(torch.autograd.Function):
    """Detach from PyTorch autograd pipeline to allow managing gradient manually"""
    @staticmethod
    def forward(ctx, dummy: torch.Tensor, tensors: _Dummy):
        ctx.set_materialize_grads(False)
        for tensor in tensors.tensors:
            if isinstance(tensor, torch.Tensor):
                assert tensor.is_leaf
            else:
                assert tensor is None
        ctx.save_for_backward(dummy)
        return tensors
    @staticmethod
    def backward(ctx, v_tensors):
        print(v_tensors)
        return ctx.saved_tensors[0], None


class _MemoryEfficientBilagridFetch(torch.autograd.Function):
    @staticmethod
    def forward(ctx, dummy: torch.Tensor, bilagrid: List, idx=None):
        ctx.set_materialize_grads(False)
        ctx.bilagrid = bilagrid[0]
        assert ctx.bilagrid.grids.is_leaf
        ctx.save_for_backward(dummy, idx)
        if idx is None:
            return ctx.bilagrid.grids
        return ctx.bilagrid.grids[idx]
    @staticmethod
    def backward(ctx, v_sample):
        """Accumulate gradient in place"""
        bilagrid = ctx.bilagrid
        dummy, idx = ctx.saved_tensors
        if v_sample is None:
            return dummy, None, None
        if idx is None:
            if bilagrid.grids.grad is None:
                bilagrid.grids.grad = v_sample
            else:
                bilagrid.grids.grad.add_(v_sample)
        else:
            if bilagrid.grids.grad is None:
                bilagrid.grids.grad = torch.zeros_like(bilagrid.grids.data)
            bilagrid.grids.grad.index_add_(0, idx, v_sample)
        return dummy, None, None


class _BilagridFusedTotalVariationLoss(torch.autograd.Function):
    @staticmethod
    def forward(ctx, dummy: torch.Tensor, bilagrid: List):
        ctx.set_materialize_grads(False)
        ctx.bilagrid = bilagrid[0]
        assert ctx.bilagrid.grids.is_leaf
        ctx.save_for_backward(dummy)
        return fused_bilagrid_C.tv_loss_forward(ctx.bilagrid.grids)

    @staticmethod
    def backward(ctx, v_output):
        """Accumulate gradient in place"""
        bilagrid = ctx.bilagrid
        (dummy,) = ctx.saved_tensors
        if v_output is None:
            return dummy, None, None
        if bilagrid.grids.grad is None:
            bilagrid.grids.grad = fused_bilagrid_C.tv_loss_backward(bilagrid.grids, v_output.contiguous())
        else:
            fused_bilagrid_C.tv_loss_backward_inplace(bilagrid.grids, v_output.contiguous(), bilagrid.grids.grad)
        return dummy, None, None


class _BilagridFusedRegularization(torch.autograd.Function):
    @staticmethod
    def forward(ctx, dummy: torch.Tensor, bilagrid: List):
        ctx.set_materialize_grads(False)
        ctx.bilagrid = bilagrid
        assert ctx.bilagrid.grids.is_leaf
        ctx.save_for_backward(dummy)
        return (
            fused_bilagrid_C.tv_loss_forward(ctx.bilagrid.grids),
            fused_bilagrid_C.channel_mean_forward(ctx.bilagrid.grids)
        )

    @staticmethod
    def backward(ctx, v_tv_loss, v_channel_mean):
        """Accumulate gradient in place"""
        bilagrid = ctx.bilagrid
        (dummy,) = ctx.saved_tensors
        if v_tv_loss is not None:
            if bilagrid.grids.grad is None:
                bilagrid.grids.grad = fused_bilagrid_C.tv_loss_backward(bilagrid.grids, v_tv_loss.contiguous())
            else:
                fused_bilagrid_C.tv_loss_backward_inplace(bilagrid.grids, v_tv_loss.contiguous(), bilagrid.grids.grad)
        if v_channel_mean is not None:
            if bilagrid.grids.grad is None:
                bilagrid.grids.grad = fused_bilagrid_C.channel_mean_backward(bilagrid.grids, v_channel_mean.contiguous())
            else:
                fused_bilagrid_C.channel_mean_backward_inplace(bilagrid.grids, v_channel_mean.contiguous(), bilagrid.grids.grad)
        return dummy, None, None


class FusedSSIM(torch.autograd.Function):
    @staticmethod
    def forward(ctx, img1, img2, train=True, return_ssim_map=False, is_l1=False):
        ssim, ssim_map, dm_dmu1, dm_dsigma1_sq, dm_dsigma12 = \
            _make_lazy_cuda_func("fused_ssim_forward")(img1, img2, train, return_ssim_map, is_l1)
        if not train:
            dm_dmu1, dm_dsigma1_sq, dm_dsigma12 = None, None, None
            if is_l1:
                raise NotImplementedError()

        ctx.save_for_backward(img1.detach(), img2, dm_dmu1, dm_dsigma1_sq, dm_dsigma12)

        return ssim, ssim_map

    @staticmethod
    def backward(ctx, dL_dssim, dL_dmap):
        # TODO: support gradient to ref_rgb
        img1, img2, dm_dmu1, dm_dsigma1_sq, dm_dsigma12 = ctx.saved_tensors
        grad = _make_lazy_cuda_func("fused_ssim_backward")(img1, img2, dL_dssim, dm_dmu1, dm_dsigma1_sq, dm_dsigma12)
        return grad, None, None, None, None


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
        opacity_reg: float,
        scale_reg: float,
        max_gauss_ratio: float,
        scale_regularization_weight: float,
        erank_reg: float,
        erank_reg_s3: float,
        quat_norm_reg_weight: float,
        compute_hessian_diagonal: Literal[None, "position", "all"] = None,
        backward_info: Optional[dict] = None,
    ):

        hyperparams = (
            opacity_reg,
            scale_reg,
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
        ctx.backward_info = backward_info

        return losses

    @staticmethod
    def backward(ctx, v_losses):

        hyperparams = ctx.hyperparams
        scales, opacities, quats = ctx.saved_tensors

        if ctx.compute_hessian_diagonal:
            assert ctx.backward_info is not None
            v_inputs, vr_inputs, h_inputs = \
            _make_lazy_cuda_func("compute_per_splat_losses_backward_with_hessian_diagonal")(
                scales, opacities, quats,
                v_losses,
                *hyperparams
            )
            for key, v, vr, h in zip('scales opacities quats'.split(), v_inputs, vr_inputs, h_inputs):
                add_gradient_component(ctx.backward_info, key+'.gradr', vr)
                add_gradient_component(ctx.backward_info, key+'.hess', h)
                # print(v.shape, torch.isfinite(h).float().mean().item(), torch.nan_to_num(v, 0, 0, 0).mean().item(), torch.nan_to_num(vr, 0, 0, 0).mean().item(), torch.nan_to_num(h, 0, 0, 0).mean().item())
        else:
            v_inputs = _make_lazy_cuda_func("compute_per_splat_losses_backward")(
                scales, opacities, quats,
                v_losses,
                *hyperparams
            )
        return (*v_inputs, *([None]*len(hyperparams)), None, None)


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

class SplatTrainingLosses:

    def __init__(self, config: 'spirulae_splat.modules.model.SpirulaeModelConfig', num_training_data):

        self.step = 0
        self.config = config
        self.num_train_data = num_training_data

        # self.supervision_losses = SupervisionLosses(config)
        # self.exposure_correction = ExposureCorrection(config)

        if self.config.use_bilateral_grid:
            self.bilagrid = {
                'affine': BilateralGrid,
                'ppisp': BilateralGridPPISP,
                'loglinear': BilateralGridLoglinear
            }[self.config.bilagrid_type](
                num=self.num_train_data,
                grid_X=self.config.bilagrid_shape[0],
                grid_Y=self.config.bilagrid_shape[1],
                grid_W=self.config.bilagrid_shape[2],
            ).cuda()
            if self.config.optimize_bilagrid_frequencies:
                raise NotImplementedError()
                self.bilagrid.grids.data = \
                    Dct3D.apply(self.bilagrid.grids.data.cuda()).to(self.bilagrid.grids.data.device)
        if self.config.use_bilateral_grid_for_geometry:
            # TODO: some way to avoid introducing VRAM overhead when geometry is not provided
            self.bilagrid_depth = BilateralGridDepth(
                num=self.num_train_data,
                grid_X=self.config.bilagrid_shape_geometry[0],
                grid_Y=self.config.bilagrid_shape_geometry[1],
                grid_W=self.config.bilagrid_shape_geometry[2],
            ).cuda()
            self.bilagrid_normal = BilateralGridNormal(
                num=self.num_train_data,
                grid_X=self.config.bilagrid_shape_geometry[0],
                grid_Y=self.config.bilagrid_shape_geometry[1],
                grid_W=self.config.bilagrid_shape_geometry[2],
            ).cuda()
        if self.config.use_ppisp:
            self.ppisp_params = torch.Tensor(
                DEFAULT_PPISP_PARAMS_RQS if self.config.ppisp_param_type == "rqs" else DEFAULT_PPISP_PARAMS
            ).repeat(self.num_train_data, 1).cuda()
            self.ppisp_params = torch.nn.Parameter(self.ppisp_params)

        self._dummy = torch.nn.Parameter(torch.zeros(1, dtype=torch.float32))

        self.lpips_dtype = torch.bfloat16 if torch.cuda.is_bf16_supported() else torch.float32
        # self.lpips_dtype = torch.float32
        if self.config.lpips_lambda > 0.0:
            self.lpips = LearnedPerceptualImagePatchSimilarity(
                net_type="vgg", normalize=True
                # net_type="alex", normalize=True
            ).to(self.lpips_dtype).cuda()

    def soft_detach(self, tensors):
        dummy = _Dummy()
        dummy.tensors = tensors
        dummy = _SoftDetach.apply(self._dummy, dummy)
        return dummy.tensors

    def _get_downscale_factor(self):
        return 1   ## TODO
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

    @torch.no_grad()
    def get_gt_img(self, image: torch.Tensor):
        """Compute groundtruth image with iteration dependent downscale factor for evaluation purpose

        Args:
            image: tensor.Tensor in type uint8, uint16, float16, or float32
        """
        if image.dtype == torch.uint8:
            # image = image.float() / 255.0
            image = _make_lazy_cuda_func("uint8_image_to_float")(image)
        elif image.dtype == torch.uint16:
            # image = image.float() / 65535.0
            image = _make_lazy_cuda_func("uint16_image_to_float")(image)
        elif image.dtype == torch.float16:
            image = image.float()
        gt_img = self._downscale_if_required(image)
        return gt_img

    def apply_bilateral_grid(
            self,
            bilagrid: Union[BilateralGrid, BilateralGridDepth, BilateralGridNormal, BilateralGridLoglinear, BilateralGridPPISP],
            rgb: torch.Tensor,
            cam_idx: int,
            bilagrid_type: Optional[str] = None,
            **kwargs
        ) -> torch.Tensor:
        """rgb must be clamped to 0-1"""
        try:
            grid_idx = cam_idx
            if isinstance(grid_idx, int):
                assert False
                grid_idx = torch.tensor(grid_idx, device=rgb.device, dtype=torch.long).flatten()
            grids = bilagrid.grids[grid_idx]
            if self.config.optimize_bilagrid_frequencies:
                raise NotImplementedError()
                grids = Dct3D.apply(grids)
            out = {
                'affine': fused_bilagrid_sample,
                'ppisp': fused_bilagrid_ppisp_sample,
                'loglinear': fused_bilagrid_loglinear_sample,
                'depth': fused_bilagrid_depth_sample,
                'normal': fused_bilagrid_normal_sample,
            }[self.config.bilagrid_type if bilagrid_type is None else bilagrid_type](
                grids, None, rgb.unsqueeze(1),
                actual_width=kwargs.get('width', None),
                actual_height=kwargs.get('height', None),
                patch_offsets=kwargs.get('patch_offsets', None),
            ).squeeze(1)
        except TypeError:
            raise RuntimeError(
                "\033[93m"
                "You are likely using an incompatible version of fused_bilagrid. "
                "Please install latest fused_bilagrid from https://github.com/harry7557558/fused-bilagrid, branch `dev`."
                "\033[0m"
            )
        # if (self.step+1) % 100 == 0 and grids.shape[1] == 9:
        #     import matplotlib.pyplot as plt
        #     print(grids.shape)
        #     print(torch.std(grids))
        #     plt.imshow(out[0].detach().cpu().numpy())
        #     plt.show()
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
        weight_depth_reg = self.config.depth_distortion_reg * factor
        weight_normal_dist_reg = self.config.normal_distortion_reg * factor
        weight_rgb_dist_reg = self.config.rgb_distortion_reg * factor
        weight_normal_reg = self.config.normal_reg_weight * factor
        return (weight_depth_reg, weight_normal_dist_reg, weight_rgb_dist_reg), weight_normal_reg

    def get_alpha_reg_weight(self):
        return self.config.alpha_reg_weight * \
            min(self.step / max(self.config.alpha_reg_warmup, 1), 1)

    def __call__(self, step: int, batch, outputs, meta={}, val=False, return_grad=True):
        self.step = step

        device = outputs['rgb'].device
        camera = outputs["camera"]
        intrins = camera.intrins

        # If reference image is empty, AI generate from rendered image
        if "image" not in batch:
            raise NotImplementedError()
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
        pred_transmittance = outputs["transmittance"] if 'transmittance' in outputs else None

        gt_rgb, gt_depth, gt_normal, gt_alpha = None, None, None, None  # for loss
        gt_rgb_mask, gt_depth_mask, gt_normal_mask, gt_alpha_mask = None, None, None, None  # for masking
        none_sky_mask = None

        with torch.no_grad():
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

            # load normal
            if 'normal' in batch:
                gt_normal = batch['normal']
                if gt_normal.dtype == torch.uint8:
                    gt_normal = gt_normal.float() / (255/2) - 1.0
                gt_normal = self._downscale_if_required(gt_normal.to(device))
                gt_normal_mask = (gt_normal.sum(-1, True) > -2.366)  # background is (-1, -1, -1)

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
            gt_img_rgba = self.get_gt_img(batch["image"].to(device))
            gt_rgb = gt_img_rgba if gt_img_rgba.shape[-1] == 3 else gt_img_rgba[..., :3]

            # update alpha if image is RGBA
            if gt_img_rgba.shape[-1] == 4 and self.config.alpha_loss_weight > 0.0:
                alpha = (gt_img_rgba[..., -1].unsqueeze(-1) > 0.5)
                gt_img_rgba = gt_img_rgba[..., :3]
                gt_rgb_mask = gt_rgb_mask & alpha if gt_rgb_mask is not None else alpha
                if self.config.apply_loss_for_mask:
                    gt_alpha = gt_alpha & alpha if gt_alpha is not None else alpha

            # convert to sRGB if needed
            # don't clip; exposure correction and loss should ideally handle out-of-gamut colors
            if self.config.image_color_is_linear or self.config.image_color_gamut != None:
                color_matrix = get_color_transform_matrix(self.config.image_color_gamut)
                gt_rgb = rgb_to_srgb(gt_rgb, self.config.image_color_is_linear, color_matrix)

        # apply bilagrid for geometry
        if self.config.use_bilateral_grid_for_geometry and \
                (camera.metadata is not None and "cam_idx" in camera.metadata):
            if gt_depth is not None:
                gt_depth_pre_bilagrid = gt_depth
                gt_depth = self.apply_bilateral_grid(
                    self.bilagrid_depth,
                    gt_depth, camera.metadata["cam_idx"],
                    bilagrid_type="depth",
                    **meta
                )
            if gt_normal is not None:
                gt_normal_pre_bilagrid = gt_normal
                gt_normal = self.apply_bilateral_grid(
                    self.bilagrid_normal,
                    gt_normal, camera.metadata["cam_idx"],
                    bilagrid_type="normal",
                    **meta
                )

        # replace parts of background with random noise to discourage transparency
        background = outputs["background"] if 'background' in outputs else gt_rgb
        if self.config.randomize_background == "opaque-only":
            raise NotImplementedError("TODO")
            background_mask = gt_rgb_mask
            if none_sky_mask is not None:
                background_mask = none_sky_mask if background_mask is None else \
                    background_mask & none_sky_mask
            if background_mask is not None:
                background = torch.where(background_mask, torch.rand_like(background), background)
        if self.config.randomize_background == "non-sky-only":
            raise NotImplementedError("TODO")
            if none_sky_mask is not None:
                background = torch.where(none_sky_mask, torch.rand_like(background), background)

        # do this to make SSIM/LPIPS happier + encourage sky transparency
        if gt_rgb_mask is not None or none_sky_mask is not None:
            # raise NotImplementedError("TODO")
            # background_mask = gt_rgb_mask
            # if none_sky_mask is not None and ('max_blending' in meta or random.random() < 0.1):
            #     background_mask = none_sky_mask if background_mask is None else \
            #         background_mask & none_sky_mask
            # if background_mask is not None:
            #     pred_rgb = torch.where(background_mask, pred_rgb, background)
            pass

        # convert to sRGB if needed
        # don't clip; exposure correction and loss should ideally handle out-of-gamut colors
        if self.config.splat_color_is_linear or self.config.splat_color_gamut != None:
            raise NotImplementedError("TODO")
            color_matrix = get_color_transform_matrix(self.config.splat_color_gamut)
            pred_rgb = rgb_to_srgb(pred_rgb, self.config.splat_color_is_linear, color_matrix)

        # edge detector to guide densification
        # TODO: multi resolution
        edge_map = None
        if self.config.use_edge_aware_score and self.config.relocate_heuristic_weight != 0.0:
            edge_map = detect_edge(gt_rgb, gt_rgb_mask)  # no_grad

        # Apply bilateral grid
        pred_rgb_pre_bilagrid = None
        if self.config.use_bilateral_grid and \
                camera.metadata is not None and "cam_idx" in camera.metadata:
            pred_rgb_pre_bilagrid = pred_rgb
            pred_rgb = self.apply_bilateral_grid(
                self.bilagrid,
                pred_rgb, camera.metadata["cam_idx"],
                bilagrid_type=self.config.bilagrid_type,
                **meta
            )
            # if step >= 0 and step % 100 == 0:
            #     import matplotlib.pyplot as plt
            #     fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
            #     ax1.imshow(pred_rgb_pre_bilagrid[0].detach().cpu().numpy())
            #     ax2.imshow(pred_rgb[0].detach().cpu().numpy())
            #     ax3.imshow(gt_rgb[0].detach().cpu().numpy())
            #     plt.show()

        # Apply PPISP
        pred_rgb_pre_ppisp = None
        if self.config.use_ppisp and \
                camera.metadata is not None and "cam_idx" in camera.metadata:
            pred_rgb_pre_ppisp = pred_rgb
            indices = torch.tensor(camera.metadata["cam_idx"]).flatten().to(device)
            pred_rgb = apply_ppisp(
                self.config.ppisp_param_type,
                pred_rgb, self.ppisp_params[indices, :],
                intrins=intrins,
                actual_image_width=camera.metadata.get('actual_width', None),
                actual_image_height=camera.metadata.get('actual_height', None),
            )

        # Depth to normal
        need_depth_to_normal_grad = False
        if pred_depth_normal is None and pred_depth is not None and (pred_normal is not None or gt_normal is not None):
            is_ray_depth = (self.config.primitive not in ['3dgs', 'mip'])
            pred_depth_normal = depth_to_normal(
                pred_depth, camera.camera_type[0].upper(),
                intrins, camera.distortion_params,
                is_ray_depth
            )
            need_depth_to_normal_grad = True

        # handle multi-resolution loss

        # TODO: fix depth supervision
        # pred_depth is linear depth for 3dgs and mip primitives and ray depth otherwise
        # gt_depth is linear depth for perspective cameras and ray depth for fisheye cameras

        all_images = [[
            pred_rgb,
            gt_rgb,
            pred_depth,
            gt_depth,
            pred_normal,
            pred_depth_normal,
            gt_normal,
            pred_transmittance,
            outputs['rgb_distortion'] if 'rgb_distortion' in outputs else None,
            outputs['depth_distortion'] if 'depth_distortion' in outputs else None,
            outputs['normal_distortion'] if 'normal_distortion' in outputs else None,
            gt_alpha,
            gt_rgb_mask,
            gt_depth_mask,
            gt_normal_mask,
            gt_alpha_mask,
        ]]

        (weight_depth_dist_reg, weight_normal_dist_reg, weight_rgb_dist_reg), weight_normal_reg = \
            self.get_2dgs_reg_weights()

        loss_weights = [
            # RGB supervision
            (1.0 - self.config.l2_lambda) * (1.0 - self.config.ssim_lambda) * (1.0 - self.config.lpips_lambda),
            self.config.l2_lambda * (1.0 - self.config.ssim_lambda) * (1.0 - self.config.lpips_lambda),
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
        ]
        w_ssim = self.config.ssim_lambda * (1.0 - self.config.lpips_lambda)

        needs_input_grad = [
            True,  # pred_rgb
            False,  # gt_rgb
            True,  # pred_depth
            self.config.use_bilateral_grid_for_geometry,  # gt_depth
            True,  # pred_normal,
            True,  # pred_depth_normal,
            self.config.use_bilateral_grid_for_geometry,  # gt_normal,
            True,  # pred_transmittance,
            True, True, True,  # distortion, auto False is None
        ]

        compute_loss_map = self.config.use_loss_map or (self.config.compute_hessian_diagonal is not None)

        if not hasattr(self, '_v_losses'):
            NUM_LOSSES = 10
            v_losses = torch.ones(NUM_LOSSES, device="cuda", dtype=torch.float32)
            if v_losses.numel() > 1:
                v_losses[1] = 0.0  # psnr
            self._v_losses = v_losses
        else:
            v_losses = self._v_losses

        losses, loss_map, grads, (psnr, ssim) = _make_lazy_cuda_func("compute_multi_scale_per_pixel_losses")(
            self.config.num_loss_scales + 1,
            # 2,
            *all_images[0],
            loss_weights,
            w_ssim,
            v_losses,
            needs_input_grad,
            meta.get("num_train_data", -1),
            camera.metadata.get('cam_idx', None),
            compute_loss_map
        )
        grads = [*grads]

        # Depth to normal backward
        if need_depth_to_normal_grad:
            grads[2] += _make_lazy_cuda_func("depth_to_normal_backward")(
                camera.camera_type[0].upper(), intrins, camera.distortion_params,
                is_ray_depth, pred_depth, grads[5]
            )

        # PPISP backward
        if pred_rgb_pre_ppisp is not None:
            indices = torch.tensor(camera.metadata["cam_idx"]).flatten().to(device)  # TODO: use pre-computed
            grads[0], v_ppisp = _make_lazy_cuda_func("ppisp_backward")(
                pred_rgb_pre_ppisp, self.ppisp_params[indices, :],
                intrins,
                # camera.metadata.get('actual_width', None),
                # camera.metadata.get('actual_height', None),
                pred_rgb[0].shape[2], pred_rgb[0].shape[1],  # TODO
                grads[0], self.config.ppisp_param_type
            )
            if not hasattr(self, 'v_ppisp') or self.v_ppisp is None:
                ppisp_loss_weights = [
                    self.config.ppisp_reg_exposure_mean,
                    self.config.ppisp_reg_vig_center,
                    self.config.ppisp_reg_vig_non_pos,
                    self.config.ppisp_reg_vig_channel_var,
                    self.config.ppisp_reg_color_mean,
                    self.config.ppisp_reg_crf_channel_var,
                ]
                ppisp_regs, ppisp_raw_regs = _C.compute_ppsip_regularization_forward(
                    self.ppisp_params,
                    ppisp_loss_weights,
                    self.config.ppisp_param_type
                )
                self.v_ppisp = _C.compute_ppsip_regularization_backward(
                    self.ppisp_params,
                    ppisp_loss_weights,
                    ppisp_raw_regs,
                    torch.ones_like(ppisp_regs),
                    self.config.ppisp_param_type
                )
            self.v_ppisp[camera.metadata["cam_idx"]] += v_ppisp

        # Bilagrid backward
        if pred_rgb_pre_bilagrid is not None:
            v_bilagrid, grads[0] = {
                'affine': fused_bilagrid_C.bilagrid_uniform_sample_backward,
                'ppisp': fused_bilagrid_C.bilagrid_ppisp_uniform_sample_backward,
                'loglinear': fused_bilagrid_C.bilagrid_loglinear_uniform_sample_backward,
            }[self.config.bilagrid_type](
                self.bilagrid.grids[camera.metadata["cam_idx"]],  # TODO: use pre-computed
                pred_rgb_pre_bilagrid.unsqueeze(1),
                grads[0].unsqueeze(1),
                1, 8, 8, 5  # TODO
            )
            pred_rgb_pre_bilagrid = None
            grads[0] = grads[0].squeeze(1)
            if not hasattr(self, 'v_bilagrid') or self.v_bilagrid is None:
                # TODO: fuse this into optimizer
                self.v_bilagrid = fused_bilagrid_C.tv_loss_backward(
                    self.bilagrid.grids, torch.full((1,), self.config.bilagrid_tv_loss_weight)
                )
                # TODO
                # fused_bilagrid_C.channel_mean_backward_inplace(
                #     self.bilagrid.grids, v_channel_mean.contiguous(), self.v_bilagrid)
            self.v_bilagrid[camera.metadata["cam_idx"]] += v_bilagrid

        # Bilagrid geometry backward
        if self.config.use_bilateral_grid_for_geometry and \
                (camera.metadata is not None and "cam_idx" in camera.metadata):
            # TODO: we don't need gradient to gt_geometry
            if gt_normal is not None:
                v_bilagrid, v_geom = fused_bilagrid_C.bilagrid_normal_uniform_sample_backward(
                    self.bilagrid_normal.grids[camera.metadata["cam_idx"]],  # TODO: use pre-computed
                    gt_normal_pre_bilagrid.unsqueeze(1),
                    grads[5].unsqueeze(1),
                    1, 8, 8, 5  # TODO
                )
                gt_normal_pre_bilagrid = None
                if not hasattr(self, 'v_bilagrid_normal') or self.v_bilagrid_normal is None:
                    # TODO: fuse this into optimizer
                    self.v_bilagrid_normal = fused_bilagrid_C.tv_loss_backward(
                        self.bilagrid_normal.grids, torch.full((1,), self.config.bilagrid_tv_loss_weight_geometry)
                    )
                self.v_bilagrid_normal[camera.metadata["cam_idx"]] += v_bilagrid
            if gt_depth is not None:
                v_bilagrid, v_geom = fused_bilagrid_C.bilagrid_depth_uniform_sample_backward(
                    self.bilagrid_depth.grids[camera.metadata["cam_idx"]],  # TODO: use pre-computed
                    gt_depth_pre_bilagrid.unsqueeze(1),
                    fused_bilagrid_C.compute_depth_scalars(gt_depth_pre_bilagrid, False),  # TODO: use pre-computed
                    grads[2].unsqueeze(1),
                    1, 8, 8, 5  # TODO
                )
                gt_depth_pre_bilagrid = None
                if not hasattr(self, 'v_bilagrid_depth') or self.v_bilagrid_depth is None:
                    # TODO: fuse this into optimizer
                    self.v_bilagrid_depth = fused_bilagrid_C.tv_loss_backward(
                        self.bilagrid_depth.grids, torch.full((1,), self.config.bilagrid_tv_loss_weight_geometry)
                    )
                self.v_bilagrid_depth[camera.metadata["cam_idx"]] += v_bilagrid
            v_geom = None

        (
            rgb_loss, rgb_psnr,
            depth_supervision_loss, normal_supervision_loss, alpha_supervision_loss,
            normal_reg, alpha_reg,
            rgb_dist_reg, depth_dist_reg, normal_dist_reg
        ) = losses
        if loss_map is not None:
            if 'backward_info' in outputs:
                outputs['backward_info']['loss_map'] = loss_map  # to be able to get it in backward
        if edge_map is not None:
            if 'backward_info' in outputs:
                outputs['backward_info']['accum_weight_map'] = edge_map
        elif self.config.relocate_heuristic_weight != 0.0 and loss_map is not None:
            if 'backward_info' in outputs:
                outputs['backward_info']['accum_weight_map'] = loss_map
        # if self.step % 100 == 0:
        #     import matplotlib.pyplot as plt
        #     fig, ((ax0, ax1), (ax2, ax3)) = plt.subplots(2, 2)
        #     ax0.imshow(gt_rgb[0].detach().cpu().numpy())
        #     if gt_rgb_mask is not None:
        #         ax1.imshow(gt_rgb_mask[0].detach().cpu().numpy())
        #     if loss_map is not None:
        #         ax2.imshow(loss_map[0].detach().cpu().numpy())
        #     if edge_map is not None:
        #         ax3.imshow(edge_map[0].detach().cpu().numpy())
        #     plt.tight_layout()
        #     plt.show()

        # note that rgb_loss is already multipled by (1.0 - self.config.ssim_lambda) * (1.0 - self.config.lpips_lambda)
        image_loss = rgb_loss.item() + w_ssim * (1.0 - ssim)

        # LPIPS for training
        if self.config.lpips_lambda > 0.0:
            raise NotImplementedError()
            if self.config.compute_hessian_diagonal is not None:
                raise NotImplementedError()
            if self.config.num_loss_scales != 0:
                raise NotImplementedError()
            lpips = self.lpips(
                pred_rgb.permute(0, 3, 1, 2).to(self.lpips_dtype).clip(0, 1),
                gt_rgb.permute(0, 3, 1, 2).to(self.lpips_dtype).clip(0, 1)
            ).float()
            image_loss = image_loss + self.config.lpips_lambda * lpips

        # LPIPS for validation
        if val:
            raise NotImplementedError()
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

        # metrics, readable from terminal during training
        with torch.no_grad():
            # list_cap_max = self.num_train_data
            list_cap_max = self.config.refine_every
            if not hasattr(self, '_running_metrics'):
                self._running_metrics = { 'psnr': [], 'ssim': [], 'lpips': [] }
            psnr_list = self._running_metrics['psnr']
            ssim_list = self._running_metrics['ssim']
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
            raise NotImplementedError()
            loss_dict['lpips_val'] = float(lpips_val)

        return loss_dict, grads

    def get_static_losses(self, step: int, gauss_quats, gauss_scales, gauss_opacities, loss_dict, backward_info: Optional[dict] = None):
        """Separately process losses that are not dependent on images"""
        self.step = step

        # bilagrid regularization loss
        if self.config.use_bilateral_grid:
            if self.config.optimize_bilagrid_frequencies:
                raise NotImplementedError()
                bilagrid = Dct3D.apply(bilagrid)
                tv_loss_fun = total_variation_loss
            tv_loss, bilagrid_mean = _BilagridFusedRegularization.apply(self._dummy, self.bilagrid)

            # total variation loss
            loss_dict["tv_loss"] = self.config.bilagrid_tv_loss_weight * tv_loss

            # channel mean
            if self.config.bilagrid_type == "affine":
                bilagrid_mean_reg = F.mse_loss(bilagrid_mean, DEFAULT_BILAGRID_PARAMS)
            elif self.config.bilagrid_type in ["ppisp", "loglinear"]:
                bilagrid_mean_reg = F.mse_loss(bilagrid_mean, torch.zeros_like(bilagrid_mean))
            loss_dict['bilagrid_mean_reg'] = self.config.bilagrid_mean_reg_weight * bilagrid_mean_reg

        # bilagrid regularization loss for geometry
        if self.config.use_bilateral_grid_for_geometry:
            if self.config.depth_supervision_weight > 0.0:  # do this because bilagrid backward is expensive, especially in patched mode
                loss_dict["tv_loss_depth"] = self.config.bilagrid_tv_loss_weight_geometry * total_variation_loss(self.bilagrid_depth.grids)
            loss_dict["tv_loss_normal"] = self.config.bilagrid_tv_loss_weight_geometry * total_variation_loss(self.bilagrid_normal.grids)

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

        losses = _ComputePerSplatLosses.apply(
            gauss_scales, gauss_opacities, gauss_quats,
            self.config.opacity_reg * float(self.config.use_mcmc),
            self.config.scale_reg * float(self.config.use_mcmc),
            self.config.max_gauss_ratio,
            self.config.scale_regularization_weight,
            self.config.erank_reg,
            self.config.erank_reg_s3,
            self.config.quat_norm_reg,
            self.config.compute_hessian_diagonal,
            backward_info,
        )
        loss_dict['opacity_reg'] = losses[0]
        loss_dict['scale_reg'] = losses[1]
        loss_dict['scale_reg'] = losses[2]
        loss_dict['erank_reg'] = losses[3]
        loss_dict['quat_norm_reg'] = losses[4]
        return loss_dict

    def optim_step(self, optim_config: OptimizerConfig, step: int, max_steps: int):

        if hasattr(self, 'v_bilagrid') and self.v_bilagrid is not None:
            assert self.bilagrid.grids.numel() == self.v_bilagrid.numel()
            if not hasattr(self, 'g1_bilagrid') is None:
                self.g1_bilagrid = torch.zeros_like(self.bilagrid.grids)
                self.g2_bilagrid = torch.zeros_like(self.bilagrid.grids)
            lr = optim_config.get_scheduled_lr("bilagrid", self.step, max_steps)
            _make_lazy_cuda_func("fused_adam_with_steps")(
                self.bilagrid.grids, self.v_bilagrid, self.g1_bilagrid, self.g2_bilagrid,
                lr, step+1, 0.0, 0.0
            )
            self.v_bilagrid = None

        if hasattr(self, 'v_ppisp') and self.v_ppisp is not None:
            if not hasattr(self, 'g1_ppisp') is None:
                self.g1_ppisp = torch.zeros_like(self.ppisp_params)
                self.g2_ppisp = torch.zeros_like(self.ppisp_params)
            lr = optim_config.get_scheduled_lr("ppisp", self.step, max_steps)
            _make_lazy_cuda_func("fused_adam_with_steps")(
                self.ppisp_params, self.v_ppisp, self.g1_ppisp, self.g2_ppisp,
                lr, step+1, 0.0, 0.0
            )
            self.v_ppisp = None

        if hasattr(self, 'v_bilagrid_depth') and self.v_bilagrid_depth is not None:
            assert self.bilagrid_depth.grids.numel() == self.v_bilagrid_depth.numel()
            if not hasattr(self, 'g1_bilagrid_depth') is None:
                self.g1_bilagrid_depth = torch.zeros_like(self.bilagrid_depth.grids)
                self.g2_bilagrid_depth = torch.zeros_like(self.bilagrid_depth.grids)
            lr = optim_config.get_scheduled_lr("bilagrid_depth", self.step, max_steps)
            _make_lazy_cuda_func("fused_adam_with_steps")(
                self.bilagrid_depth.grids, self.v_bilagrid_depth, self.g1_bilagrid_depth, self.g2_bilagrid_depth,
                lr, step+1, 0.0, 0.0
            )
            self.v_bilagrid_depth = None

        if hasattr(self, 'v_bilagrid_normal') and self.v_bilagrid_normal is not None:
            assert self.bilagrid_normal.grids.numel() == self.v_bilagrid_normal.numel()
            if not hasattr(self, 'g1_bilagrid_depth') is None:
                self.g1_bilagrid_normal = torch.zeros_like(self.bilagrid_normal.grids)
                self.g2_bilagrid_normal = torch.zeros_like(self.bilagrid_normal.grids)
            lr = optim_config.get_scheduled_lr("bilagrid_normal", self.step, max_steps)
            _make_lazy_cuda_func("fused_adam_with_steps")(
                self.bilagrid_normal.grids, self.v_bilagrid_normal, self.g1_bilagrid_normal, self.g2_bilagrid_normal,
                lr, step+1, 0.0, 0.0
            )
            self.v_bilagrid_normal = None
