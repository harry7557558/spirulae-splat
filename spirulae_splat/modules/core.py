import torch
import torch.nn.functional as F
from torch import Tensor
from typing import Dict, Optional, Tuple, Literal
import math
import numpy as np


import spirulae_splat

from spirulae_splat.splat.cuda._wrapper import (
    intersect_splat_tile,
    fully_fused_projection,
    fully_fused_projection_hetero,
    rasterize_to_pixels,
    spherical_harmonics,
)

from spirulae_splat.splat.cuda import (
    _C,
    _make_lazy_cuda_func,
)

from spirulae_splat.modules.optimizer import OptimizerConfig


class Renderer:

    def __init__(
        self,
        primitive: Literal["3dgs", "mip", "3dgut", "3dgut_sv", "opaque_triangle", "voxel"],
        splats_world: tuple[Tensor],  # means, quats, scales, opacities
        cur_num_splats: int,
        packed: bool = True,
        use_bvh: bool = False,
        use_fused_proj_bwd_optim: bool = False,
        quantize_sh_optim: bool = True,
    ):
        for tensor in splats_world:
            assert tensor.is_contiguous(), "Tensor must be contiguous"

        features_dc, features_sh, sv_sites, sv_colors = None, None, None, None
        if primitive in ["3dgs", "mip", "3dgut", "3dgut_sv", "opaque_triangle"]:
            if primitive in ["3dgs", "mip", "3dgut"]:
                assert len(splats_world) == 6, "3DGS requires 6 params (means, quats, scales, opacities, color, sh)"
                means, quats, scales, opacities, features_dc, features_sh = splats_world
            elif primitive in ["3dgut_sv"]:
                assert len(splats_world) == 6, "3DGS SV requires 6 params (means, quats, scales, opacities, sv_dir, sv_color)"
                means, quats, scales, opacities, sv_sites, sv_colors = splats_world
            else:
                assert len(splats_world) == 7, "Opaque triangle requires 4 params (means, quats, scales, opacities, color, ch, sh)"
                means, quats, scales, opacities, features_dc, features_sh, features_ch = splats_world
            assert len(means.shape) == 2, "means must have at 2 dimensions"
            N = means.shape[-2]
            device = means.device
            assert means.shape == (N, 3), means.shape
            assert quats.shape == (N, 4), quats.shape
            assert scales.shape == (N, 3), scales.shape
            if primitive in ["3dgs", "mip", "3dgut", "3dgut_sv"]:
                assert opacities.shape == (N, 1), opacities.shape
            else:
                assert opacities.shape == (N, 2), opacities.shape
                assert features_ch.shape == (N, 2, 3), features_ch.shape
        elif primitive in ["voxel"]:
            assert len(splats_world) == 4, "Voxel requires 4 params (pos_sizes, densities, features_dc, features_sh)"
            pos_sizes, densities, features_dc, features_sh = splats_world
            N = pos_sizes.shape[-2]
            device = pos_sizes.device
            assert pos_sizes.shape == (N, 4), pos_sizes.shape
            assert densities.shape == (N, 8), densities.shape
        else:
            raise ValueError(f"Invalid primitive ({primitive})")
        if features_dc is not None and features_sh is not None:
            assert features_dc.shape == (N, 3), features_dc.shape
            assert features_sh.shape == (N, 0, 3) or \
                features_sh.shape == (N, 3, 3) or \
                features_sh.shape == (N, 8, 3) or \
                features_sh.shape == (N, 15, 3) or \
                features_sh.shape == (N, 24, 3), features_sh.shape
        if sv_sites is not None and sv_colors is not None:
            num_sv = sv_sites.shape[-2]
            assert sv_sites.shape == (N, num_sv, 3)
            assert sv_colors.shape == (N, num_sv, 3)

        self.device = device
        self.cur_num_splats = cur_num_splats
        self.max_num_splats = N

        self.primitive = primitive
        self.splats_world = splats_world


        self.packed = packed
        self.use_bvh = use_bvh
        self.use_fused_proj_bwd_optim = use_fused_proj_bwd_optim
        self.quantize_sh_optim = quantize_sh_optim

        # Set by the trainer when train_frame != "normalized". Equal to
        # 1 / dataparser_scale_factor — how much bigger the training frame is
        # than the would-be normalized frame. Used to rescale means_lr,
        # scale_reg, max_world_size, and noise_lr in the engine_* calls.
        self.train_frame_scale = 1.0

        # Engine path: upload splats to device pool once at init (idempotent on C++ side).
        if self.primitive in ['3dgs', 'mip', '3dgut'] and not self.use_bvh:
            _C.set_data_3dgs(
                self.cur_num_splats,
                *[self._tv(t) for t in self.splats_world]
            )

    def set_params(
        self,
        viewmats: Tensor,  # [..., C, 4, 4]
        intrins: Tensor,  # [..., C, 4]
        width: int,
        height: int,
        sh_degree_to_use: int = 4,
        output_distortion: bool = False,
        compute_hessian_diagonal: Literal[None, "position", "all"] = None,
        relative_scale: Optional[float] = None,
        backgrounds: Optional[Tensor] = None,
        accum_weight_map: Optional[Tensor] = None,
        max_blending_masks: Optional[Tensor] = None,
        camera_model: Literal["pinhole", "ortho", "fisheye", "equisolid"] = "pinhole",
        dist_coeffs: Optional[Tensor] = None,  # [..., C, 10]
        actual_width: int = None,
        actual_height: int = None
    ):
        self.viewmats = viewmats.cpu().contiguous()
        self.intrins = intrins.cpu().contiguous()
        self.width = width
        self.height = height
        self.sh_degree_to_use = sh_degree_to_use
        self.output_distortion = output_distortion
        self.compute_hessian_diagonal = compute_hessian_diagonal
        self.relative_scale = relative_scale
        self.backgrounds = backgrounds
        self.accum_weight_map = accum_weight_map
        self.max_blending_masks = max_blending_masks
        self.camera_model = camera_model
        self.dist_coeffs = dist_coeffs.cpu() if dist_coeffs is not None else None
        self.actual_width = actual_width
        self.actual_height = actual_height
        if self.dist_coeffs is not None:
            self.dist_coeffs = self.dist_coeffs.contiguous()

    @staticmethod
    def _tv(tensor: torch.Tensor):
        """Convert torch.Tensor to (data_ptr, element_size, shape) for C++ Engine calls.
        Works with both CPU and CUDA tensors."""
        if tensor is None:
            return (0, 4, [0])
        assert tensor.is_contiguous(), f"Tensor must be contiguous, got strides {tensor.stride()}"
        return (tensor.data_ptr(), tensor.element_size(), list(tensor.shape))

    def engine_forward(self):
        """Runs projection + tile intersection + rasterization via C++ Engine."""
        C = self.viewmats.shape[0]
        H, W = self.height, self.width

        dist_coeffs = self.dist_coeffs
        if dist_coeffs is None:
            dist_coeffs = torch.zeros(C, 10, dtype=torch.float32)

        _C.set_data_3dgs(
            self.cur_num_splats,
            *[self._tv(t) for t in self.splats_world]
        )
        _C.set_camera_params(
            self.width, self.height,
            self.camera_model.upper(),
            self._tv(self.viewmats),
            self._tv(self.intrins),
            self._tv(dist_coeffs)
        )
        _C.engine_forward_3dgs(
            self.primitive,
            self.sh_degree_to_use,
            self.packed,
        )

        rgb = torch.empty(C, H, W, 3, dtype=torch.float32)
        depth = torch.empty(C, H, W, 1, dtype=torch.float32)
        Ts = torch.empty(C, H, W, 1, dtype=torch.float32)
        _C.engine_copy_render_to_host(self._tv(rgb), self._tv(depth), self._tv(Ts))

        self.render_colors = (rgb, depth)
        self.render_Ts = Ts

    def engine_debug_forward(self, override_features_dc=None, override_sh_degree=-1):
        """Re-render with custom features_dc and/or sh_degree for debugging.
        Returns CPU RGB tensor [C, H, W, 3]."""
        C = self.viewmats.shape[0]
        H, W = self.height, self.width
        out_rgb = torch.zeros(C, H, W, 3, dtype=torch.float32)
        dc_tv = self._tv(override_features_dc) if override_features_dc is not None else (0, 4, [0])
        _C.engine_debug_forward(dc_tv, override_sh_degree, self._tv(out_rgb))
        return out_rgb

    def engine_get_accum_weight(self):
        """Copy per-splat accum_buffer from C++ pool (D→H), compute weight as col0/col1."""
        N = _C.engine_get_cur_num_splats()
        buf = torch.zeros(N, 2, dtype=torch.float32)
        _C.engine_copy_accum_buffer(self._tv(buf))
        weight = torch.where(buf[:, 1] != 0, buf[:, 0] / buf[:, 1], torch.zeros_like(buf[:, 0]))
        return weight

    def engine_set_training_data(self, gt_rgb, gt_depth=None, gt_normal=None,
                                  gt_alpha=None):
        # 4-mask buffers (rgb/depth/normal/alpha) gone — derived in the slang
        # kernel from gt_alpha / gt_depth=0 / gt_normal sentinel.
        _C.set_training_data(
            self._tv(gt_rgb), self._tv(gt_depth), self._tv(gt_normal), self._tv(gt_alpha)
        )

    def engine_compute_loss_backward(self, step, loss_weights, w_ssim, num_loss_scales, compute_loss_map):
        """Compute loss + rasterization backward + projection backward in C++.
        Gradients are managed by C++ pool."""
        loss_dict = _C.engine_compute_loss_backward(
            step, loss_weights, w_ssim, num_loss_scales, compute_loss_map
        )
        return loss_dict

    def _build_optim_config(self, step, max_steps, model_config, optim_config):
        """Build a C++ OptimConfig with all per-step LRs + regularization
        weights already resolved (incl. train_frame_scale / alpha)."""
        if optim_config.max_steps is not None:
            max_steps = optim_config.max_steps
        alpha = self.train_frame_scale
        means_lr = optim_config.get_scheduled_lr("means", step, max_steps)
        if not optim_config.use_scale_agnostic_mean:
            means_lr *= alpha

        c = _C.OptimConfig()
        c.lr_means       = means_lr
        c.lr_quats       = optim_config.get_scheduled_lr("quats",       step, max_steps)
        c.lr_scales      = optim_config.get_scheduled_lr("scales",      step, max_steps)
        c.lr_opacities   = optim_config.get_scheduled_lr("opacities",   step, max_steps)
        c.lr_features_dc = optim_config.get_scheduled_lr("features_dc", step, max_steps)
        c.lr_features_sh = optim_config.get_scheduled_lr("features_sh", step, max_steps)
        c.max_gauss_ratio             = model_config.max_gauss_ratio
        c.scale_regularization_weight = model_config.scale_regularization_weight
        c.mcmc_opacity_reg_weight     = model_config.opacity_reg
        c.mcmc_scale_reg_weight       = model_config.scale_reg / alpha
        c.erank_reg_weight            = model_config.erank_reg
        c.erank_reg_weight_s3         = model_config.erank_reg_s3
        c.quat_norm_reg_weight        = model_config.quat_norm_reg
        c.sh_reg_weight               = model_config.sh_reg
        c.use_scale_agnostic_mean        = optim_config.use_scale_agnostic_mean
        c.quantize_sh                    = self.quantize_sh_optim
        c.use_per_splat_bias_correction  = optim_config.use_per_splat_bias_correction
        return c

    def _build_densify_config(self, model_config):
        alpha = self.train_frame_scale
        noise_lr_scalar = 1.0 if model_config.relocate_heuristic_weight >= 1.0 else alpha
        c = _C.DensifyConfig()
        c.refine_start_iter             = model_config.refine_start_iter
        c.refine_stop_num_iter          = model_config.refine_stop_num_iter
        c.refine_every                  = model_config.refine_every
        c.growth_factor                 = model_config.growth_factor
        c.min_opacity                   = model_config.min_opacity
        c.max_screen_size               = model_config.max_screen_size
        c.max_screen_size_clip_hardness = model_config.max_screen_size_clip_hardness
        c.max_world_size                = model_config.max_world_size * alpha
        c.noise_lr                      = model_config.noise_lr * noise_lr_scalar
        c.noise_lr_final                = model_config.noise_lr_final * noise_lr_scalar
        c.relocate_heuristic_weight     = model_config.relocate_heuristic_weight
        return c

    def engine_optim_step(self, step, max_steps, model_config, optim_config):
        """Run optimizer step via C++ Engine. All tensors managed by C++ pool."""
        _C.engine_optim_step(
            step,
            self._build_optim_config(step, max_steps, model_config, optim_config),
        )

    def engine_init_bilagrid(self, n_grids, rgb_type=None, rgb_LHW=None,
                              depth_LHW=None, normal_LHW=None,
                              quantize_rgb_optim=False,
                              quantize_geometry_optim=False):
        """One-time allocation of C++ bilagrid storage. Pass `None` for any type
        to leave it disabled. Must be called after set_data_3dgs (i.e., after
        engine_forward has been invoked at least once).

        quantize_*_optim: store the Adam g1/g2 moments as uint8 + per-block
        float4 bounds (same scheme as SH quantization). Separate per-axis
        because the RGB grid is much larger than the geometry grids in
        photogrammetry-scale datasets."""
        if rgb_type is not None and rgb_LHW is not None:
            L, H, W = rgb_LHW
            _C.engine_init_bilagrid_rgb(n_grids, rgb_type, L, H, W,
                                        bool(quantize_rgb_optim))
        if depth_LHW is not None:
            L, H, W = depth_LHW
            _C.engine_init_bilagrid_depth(n_grids, L, H, W,
                                          bool(quantize_geometry_optim))
        if normal_LHW is not None:
            L, H, W = normal_LHW
            _C.engine_init_bilagrid_normal(n_grids, L, H, W,
                                           bool(quantize_geometry_optim))

    def engine_bilagrid_forward(self, cam_indices):
        """Apply bilagrid forward for the current batch. cam_indices: [C_batch]
        int32 tensor of per-image camera-table indices (or None for identity).
        In-place on rendered RGB and on the GT depth/normal buffers."""
        _C.engine_bilagrid_forward(self._tv(cam_indices))

    @staticmethod
    def _build_bilagrid_step_config(lr_rgb, lr_depth, lr_normal,
                                     tv_weight_rgb, tv_weight_depth, tv_weight_normal):
        c = _C.BilagridStepConfig()
        c.lr_rgb           = float(lr_rgb)
        c.lr_depth         = float(lr_depth)
        c.lr_normal        = float(lr_normal)
        c.tv_weight_rgb    = float(tv_weight_rgb)
        c.tv_weight_depth  = float(tv_weight_depth)
        c.tv_weight_normal = float(tv_weight_normal)
        return c

    def engine_bilagrid_optim_step(self, step, lr_rgb, lr_depth, lr_normal,
                                    tv_weight_rgb, tv_weight_depth, tv_weight_normal):
        """Adam step + TV-loss regularization for each enabled bilagrid type."""
        _C.engine_bilagrid_optim_step(step, self._build_bilagrid_step_config(
            lr_rgb, lr_depth, lr_normal,
            tv_weight_rgb, tv_weight_depth, tv_weight_normal))

    def engine_init_ppisp(self, n_grids, param_type="original"):
        """One-time allocation of C++ PPISP per-camera parameter table. Seeded
        with the type's default values. Must be called after set_data_3dgs."""
        _C.engine_init_ppisp(int(n_grids), str(param_type))

    def engine_ppisp_forward(self, cam_indices):
        """Apply PPISP forward in place on the current rendered RGB.
        cam_indices: [C_batch] int32 tensor (or None for identity)."""
        _C.engine_ppisp_forward(self._tv(cam_indices))

    @staticmethod
    def _build_ppisp_step_config(lr,
                                 reg_exposure_mean, reg_vig_center,
                                 reg_vig_non_pos, reg_vig_channel_var,
                                 reg_color_mean, reg_crf_channel_var):
        c = _C.PpispStepConfig()
        c.lr = float(lr)
        # Order must match enum PPISPRegLossIndex
        # (ExposureMean, VignettingCenter, VignettingNonPositivity,
        #  VignettingChannelVariance, ColorMean, CRFChannelVariance).
        c.reg_weights = [
            float(reg_exposure_mean), float(reg_vig_center),
            float(reg_vig_non_pos),   float(reg_vig_channel_var),
            float(reg_color_mean),    float(reg_crf_channel_var),
        ]
        return c

    def engine_ppisp_optim_step(self, step, lr,
                                 reg_exposure_mean, reg_vig_center,
                                 reg_vig_non_pos, reg_vig_channel_var,
                                 reg_color_mean, reg_crf_channel_var):
        """Adam step over PPISP params + folds in the 6 regularization losses."""
        _C.engine_ppisp_optim_step(int(step), self._build_ppisp_step_config(
            lr, reg_exposure_mean, reg_vig_center,
            reg_vig_non_pos, reg_vig_channel_var,
            reg_color_mean, reg_crf_channel_var))

    def engine_save_checkpoint(self, output_dir, step, full_dump=False):
        """Write checkpoint files into ``output_dir`` (created if missing).
        Always writes ``splat.ply`` (cur_num_splats only, NaN/low-opacity-filtered),
        ``bilagrid_*.npy`` / ``ppisp.npy`` (when enabled) and ``meta.txt``.
        When ``full_dump=True`` additionally writes a ``full/`` subfolder with
        every world parameter (max_num_splats) plus Adam optimizer state."""
        import os
        os.makedirs(str(output_dir), exist_ok=True)
        _C.engine_save_checkpoint(str(output_dir), bool(full_dump), int(step))

    def engine_train_step(self, step, max_steps,
                          # Forward
                          sh_degree_to_use,
                          width, height, camera_model,
                          viewmats, intrins, dist_coeffs,
                          # GT data. Per-pixel masks are derived in the slang
                          # kernel; gt_alpha drives the RGB mask + alpha-sup
                          # target. apply_loss_for_mask folded into weights
                          # on the model.py side.
                          gt_rgb, gt_depth, gt_normal, gt_alpha,
                          # Loss config
                          loss_weights, w_ssim, num_loss_scales, compute_loss_map,
                          # Configs
                          model_config, optim_config,
                          # Bilagrid (pass cam_indices + lrs + tv weights;
                          # ignored if engine_init_bilagrid_* was never called).
                          # cam_indices: [C_batch] int32 tensor, or None.
                          bilagrid_cam_indices=None,
                          bilagrid_lr_rgb=0.0, bilagrid_lr_depth=0.0, bilagrid_lr_normal=0.0,
                          bilagrid_tv_weight_rgb=0.0, bilagrid_tv_weight_depth=0.0,
                          bilagrid_tv_weight_normal=0.0,
                          # PPISP (ignored if engine_init_ppisp was never called).
                          # Reuses bilagrid_cam_idx for the camera selector.
                          ppisp_lr=0.0,
                          ppisp_reg_exposure_mean=0.0, ppisp_reg_vig_center=0.0,
                          ppisp_reg_vig_non_pos=0.0, ppisp_reg_vig_channel_var=0.0,
                          ppisp_reg_color_mean=0.0, ppisp_reg_crf_channel_var=0.0):
        """Single fused training step (set_camera + set_gt + fwd + loss/bwd + optim + densify).
        All input tensors are CPU; returns loss_dict for verbose."""
        max_steps_lr = optim_config.max_steps if optim_config.max_steps is not None else max_steps

        cfg = _C.EngineStepConfig()
        cfg.loss.weights          = loss_weights
        cfg.loss.w_ssim           = float(w_ssim)
        cfg.loss.num_loss_scales  = int(num_loss_scales)
        cfg.loss.compute_loss_map = bool(compute_loss_map)
        cfg.optim    = self._build_optim_config(step, max_steps_lr, model_config, optim_config)
        cfg.densify  = self._build_densify_config(model_config)
        cfg.bilagrid = self._build_bilagrid_step_config(
            bilagrid_lr_rgb, bilagrid_lr_depth, bilagrid_lr_normal,
            bilagrid_tv_weight_rgb, bilagrid_tv_weight_depth, bilagrid_tv_weight_normal)
        cfg.ppisp    = self._build_ppisp_step_config(
            ppisp_lr, ppisp_reg_exposure_mean, ppisp_reg_vig_center,
            ppisp_reg_vig_non_pos, ppisp_reg_vig_channel_var,
            ppisp_reg_color_mean, ppisp_reg_crf_channel_var)

        result = _C.engine_train_step(
            step, max_steps,
            self.primitive, sh_degree_to_use, self.packed,
            width, height, camera_model.upper(),
            self._tv(viewmats), self._tv(intrins), self._tv(dist_coeffs),
            self._tv(gt_rgb), self._tv(gt_depth), self._tv(gt_normal), self._tv(gt_alpha),
            self._tv(bilagrid_cam_indices),
            cfg,
        )
        # Update Python-side cur_num_splats (densification happened in C++)
        num_added = int(result.pop("num_added", 0))
        result.pop("cur_num_splats", None)
        result.pop("max_num_splats", None)
        self.cur_num_splats += num_added
        self.sh_degree_to_use = sh_degree_to_use
        return result

    def engine_densify_step(self, step, max_steps, model_config):
        """Run densification step via C++ Engine. All tensors managed by C++ pool."""
        num_added = _C.engine_densify_step(step, max_steps,
                                           self._build_densify_config(model_config))
        self.cur_num_splats += num_added

    def engine_sync_splats_to_host(self):
        """Copy device splat data back to CPU PyTorch parameters after optim+densify."""
        _C.engine_copy_splats_to_host(*[self._tv(t) for t in self.splats_world])

    def zero_grad(self):
        return
        raise NotImplementedError("Use engine instead. Code below for reference.")

        # Engine path: gradients, optimizer state, radii managed by C++ pool
        if self.primitive in ['3dgs', 'mip', '3dgut'] and not self.use_bvh:
            return

        if hasattr(self, 'v_splats_world'):
            for tensor in self.v_splats_world:
                if tensor is not None:
                    tensor.zero_()
        else:
            if self.use_fused_proj_bwd_optim:
                if self.primitive in ['3dgs', 'mip']:
                    self.v_splats_world = [None] * len(self.splats_world)
                elif self.primitive in ['3dgut', '3dgut_sv']:
                    self.v_splats_world = [
                        torch.zeros_like(self.splats_world[0]),  # means
                        torch.zeros_like(self.splats_world[1]),  # quats
                        torch.zeros_like(self.splats_world[2]),  # scales
                        None,  # opacities
                        None,  # features_dc
                        None,  # features_sh
                    ]
                else:
                    raise NotImplementedError()
            else:
                self.v_splats_world = [
                    torch.zeros_like(x) for x in self.splats_world
                ]

        def alloc_optim_state():
            if self.quantize_sh_optim:
                return [torch.zeros_like(x) for x in self.splats_world[:-1]] + \
                    [torch.zeros(self.splats_world[-1].shape, device=self.splats_world[-1].device, dtype=torch.uint8)]
            return [torch.zeros_like(x) for x in self.splats_world]

        if not hasattr(self, 'g1_splats_world'):
            self.g1_splats_world = alloc_optim_state()
        if not hasattr(self, 'g2_splats_world'):
            self.g2_splats_world = alloc_optim_state()

        if self.quantize_sh_optim:
            if not hasattr(self, 'quant_bounds_sh'):
                BLOCK_SIZE = 256
                n = (self.splats_world[-1].numel() + BLOCK_SIZE-1) // BLOCK_SIZE
                self.quant_bounds_sh = torch.zeros((n, 4), dtype=torch.float32, device=self.splats_world[-1].device)

        if hasattr(self, 'radii'):
            # _make_lazy_cuda_func("set_zero")(self.radii)
            self.radii = self.radii.zero_()

    def projection_forward(self):
        raise NotImplementedError("Use engine instead. Code below for reference.")
        if self.primitive not in ["3dgs", "mip", "3dgut", "3dgut_sv"]:
            raise NotImplementedError()

        proj_returns = _make_lazy_cuda_func(
            f"projection_{self.primitive}_packed_forward" if self.packed else
            f"projection_{self.primitive}_forward"
        )(
            self.cur_num_splats,
            self.sh_degree_to_use,
            self.splats_world,
            self.viewmats, self.intrins, self.width, self.height,
            self.camera_model.upper(), self.dist_coeffs
        )
        if self.packed:
            camera_ids, gaussian_ids, aabb, sorting_depths, radii, splats_proj = proj_returns
        else:
            aabb, sorting_depths, radii, splats_proj = proj_returns
            camera_ids, gaussian_ids = None, None

        self.aabb = aabb  # xyxy
        self.depths = sorting_depths
        self.splats_proj = splats_proj
        self.camera_ids, self.gaussian_ids = camera_ids, gaussian_ids

        if not hasattr(self, 'radii'):
            self.radii = radii
        self.radii = torch.fmax(self.radii, radii)

    def rasterize_forward(self):
        raise NotImplementedError("Use engine instead. Code below for reference.")
        if self.primitive in ['3dgs', 'mip']:
            (
                self.render_colors, self.render_Ts, self.render_last_ids,
                 render2_outputs, distortion_outputs  # TODO
            ) = _make_lazy_cuda_func(
                f"rasterization_{self.primitive}_forward"
            )(
                self.splats_world, self.splats_proj, self.gaussian_ids,
                self.width, self.height, self.isect_offsets, self.flatten_ids,
                self.output_distortion
            )
        else:
            (
                self.render_colors, self.render_Ts, self.render_last_ids,
                render2_outputs, distortion_outputs  # TODO
            ) = _make_lazy_cuda_func(
                f"rasterization_{self.primitive}_forward"
            )(
                self.splats_world, self.splats_proj, self.gaussian_ids,
                self.viewmats, self.intrins, self.camera_model.upper(), self.dist_coeffs,
                self.width, self.height, self.isect_offsets, self.flatten_ids,
                self.output_distortion
            )

    def rasterize_backward(self):
        raise NotImplementedError("Use engine instead. Code below for reference.")
        if self.compute_hessian_diagonal:
            raise NotImplementedError()
            v_splats, vr_splats, h_splats, accum_weight = _make_lazy_cuda_func(
                f"rasterization_{['3dgs', 'mip'][ctx.antialiased]}_backward_with_hessian_diagonal"
            )(
                self.max_num_splats,
                (means2d, depths, conics, opacities, colors),
                ctx.backward_info.get('gaussian_ids', None),
                backgrounds, masks,
                width, height, isect_offsets, flatten_ids, render_Ts, last_ids,
                None, None,
                ctx.backward_info['loss_map'] if ctx.compute_hessian_diagonal else None,
                ctx.backward_info.get('accum_weight_map', None),
                (v_render_rgbs.contiguous(), v_render_depths.contiguous()),
                v_render_Ts.contiguous(),
                None
            )
            v_means2d, v_depths, v_conics, v_opacities, v_colors = v_splats
            for key, v, vr, h in zip('means2d depths conics proj_opacities colors'.split(), v_splats, vr_splats, h_splats):
                add_gradient_component(ctx.backward_info, key+'.gradr', vr)
                add_gradient_component(ctx.backward_info, key+'.hess', h)
        else:
            if self.primitive in ['3dgs', 'mip']:
                self.v_splats_world, self.v_splats_proj, accum_weight = _make_lazy_cuda_func(
                    f"rasterization_{self.primitive}_backward"
                )(
                    self.max_num_splats,
                    self.splats_world, self.splats_proj, self.gaussian_ids,
                    self.width, self.height, self.isect_offsets, self.flatten_ids, self.render_Ts, self.render_last_ids,
                    self.backward_info.get('accum_weight_map', None),
                    self.v_render_colors,
                    self.v_render_Ts,
                    self.v_splats_world if hasattr(self, 'v_splats_world') else None,
                    None,
                )
            else:
                cuda_return = _make_lazy_cuda_func(
                    f"rasterization_{self.primitive}_backward"
                )(
                    self.max_num_splats,
                    self.splats_world, self.splats_proj, self.gaussian_ids,
                    self.viewmats, self.intrins, self.camera_model.upper(), self.dist_coeffs,
                    self.width, self.height, self.isect_offsets, self.flatten_ids, self.render_Ts, self.render_last_ids,
                    # (render_rgbs, render_depths) if ctx.output_distortion else None,
                    # (render2_rgbs, render2_depths) if ctx.output_distortion else None,
                    # ctx.backward_info['loss_map'] if ctx.compute_hessian_diagonal else None,
                    None, None, None,  # TODO
                    self.backward_info.get('accum_weight_map', None),
                    self.v_render_colors,
                    self.v_render_Ts,
                    # (v_distortion_rgbs.contiguous(), v_distortion_depths.contiguous()) if ctx.output_distortion else None,
                    # ctx.needs_input_grad[13]
                    None,
                    self.v_splats_world if hasattr(self, 'v_splats_world') else None,
                    None,
                    False,  # TODO
                )
                if self.compute_hessian_diagonal:
                    v_splats, v_viewmats, vr_splats, h_splats, accum_weight = cuda_return
                else:
                    self.v_splats_world, self.v_splats_proj, v_viewmats, accum_weight = cuda_return
        if accum_weight is not None:
            self.backward_info['accum_weight'] = accum_weight

        # TODO
        # self.v_backgrounds = None
        # if self.needs_input_grad[5]:
        #     self.v_backgrounds = (torch.cat(self.v_render_colors, dim=-1) * \
        #                      self.v_render_Ts.float()).sum(dim=(-3, -2))

    def projection_backward(self):
        raise NotImplementedError("Use engine instead. Code below for reference.")
        if self.use_fused_proj_bwd_optim:
            return

        if self.compute_hessian_diagonal == "all":
            raise NotImplementedError()
            # TODO: in place gradient accumulation
            assert ctx.backward_info is not None
            vr_proj = (*[
                torch.zeros_like(x)
                for x in (means, quats, scales, opacities, features_dc)
            ], None)
            h_proj = (*[
                torch.zeros_like(x)
                for x in (means, quats, scales, opacities, features_dc)
            ], None)
        elif self.compute_hessian_diagonal == "position":
            raise NotImplementedError()
            vr_proj = torch.zeros_like(means)
            h_proj = torch.zeros_like(means)
        # v_viewmats = torch.zeros_like(viewmats) if ctx.needs_input_grad[7] else None
        self.v_viewmats = None  # TODO

        if self.compute_hessian_diagonal is not None:
            raise NotImplementedError()
            if ctx.primitive not in ["3dgs", "mip", "3dgut"]:
                raise NotImplementedError()
            raster_return_keys = 'means2d depths conics proj_opacities colors'.split() \
                if ctx.primitive in ['3dgs', 'mip'] else 'depths proj_scales proj_opacities colors'.split()
            assert len(v_proj_returns) == len(raster_return_keys), ([x.shape for x in v_proj_returns], raster_return_keys)
            proj_return_keys = 'means quats scales opacities features_dc'.split() \
                if ctx.compute_hessian_diagonal == "all" else ['means']
            _make_lazy_cuda_func(
                f"projection_{ctx.primitive}_backward_with_hessian_diagonal" if ctx.compute_hessian_diagonal == "all"
                else f"projection_{ctx.primitive}_backward_with_position_hessian_diagonal"
            )(
                (means, quats, scales, opacities, features_dc, features_sh),
                viewmats, intrins, ctx.width, ctx.height, ctx.camera_model_type, ctx.dist_coeffs,
                camera_ids, gaussian_ids, aabb,
                [x.contiguous() for x in v_proj_returns],
                [ctx.backward_info[key+'.gradr'] for key in raster_return_keys],
                [ctx.backward_info[key+'.hess'] for key in raster_return_keys],
                v_proj, v_viewmats, vr_proj, h_proj,
            )
            for key in raster_return_keys:
                if key not in proj_return_keys:
                    del ctx.backward_info[key+'.gradr']
                    del ctx.backward_info[key+'.hess']
            for key in proj_return_keys:
                if 'proj_'+key+'.gradr' not in ctx.backward_info and 'proj_'+key+'.hess' not in ctx.backward_info:
                    continue
                assert 'proj_'+key+'.gradr' in ctx.backward_info
                assert 'proj_'+key+'.hess' in ctx.backward_info
                if gaussian_ids is not None:  # packed
                    assert ctx.primitive in ["3dgut"]
                    # ref_tensor = {"means": means, "quats": quats}[key]
                    # ctx.backward_info[key+'.gradr'] = scatter_add(ref_tensor, ctx.backward_info['proj_'+key+'.gradr'], gaussian_ids)
                    # ctx.backward_info[key+'.hess'] = scatter_add(ref_tensor, ctx.backward_info['proj_'+key+'.hess'], gaussian_ids)
                    ctx.backward_info[key+'.gradr'] = ctx.backward_info['proj_'+key+'.gradr']
                    ctx.backward_info[key+'.hess'] = ctx.backward_info['proj_'+key+'.hess']
                    del ctx.backward_info['proj_'+key+'.gradr']
                    del ctx.backward_info['proj_'+key+'.hess']
            if ctx.compute_hessian_diagonal == "position":
                vr_proj, h_proj = [vr_proj], [h_proj]
            for key, x, vr, h in zip(proj_return_keys, ctx.saved_tensors, vr_proj, h_proj):
                add_gradient_component(ctx.backward_info, key+'.gradr', vr)
                add_gradient_component(ctx.backward_info, key+'.hess', h)
                assert hasattr(x, 'optim_info')
                x.optim_info['gradr'] = ctx.backward_info[key+'.gradr']
                x.optim_info['hess'] = ctx.backward_info[key+'.hess']
                del ctx.backward_info[key+'.gradr']
                del ctx.backward_info[key+'.hess']
        else:
            _make_lazy_cuda_func(
                f"projection_{self.primitive}_backward"
            )(
                self.cur_num_splats,
                self.sh_degree_to_use,
                self.splats_world,
                self.viewmats, self.intrins, self.width, self.height, self.camera_model.upper(), self.dist_coeffs,
                self.camera_ids, self.gaussian_ids, self.aabb,
                self.v_splats_proj,
                self.v_splats_world, self.v_viewmats,
            )

    def forward(self):

        self.meta = {}

        C = self.viewmats.shape[-3]
        assert self.viewmats.shape == (C, 4, 4), self.viewmats.shape
        assert self.intrins.shape == (C, 4), self.intrins.shape

        self.backward_info = {}

        # Engine path: unified C++ forward for 3dgs/mip/3dgut
        if self.primitive in ['3dgs', 'mip', '3dgut'] and not self.use_bvh:
            self.engine_forward()
            self.meta.update({
                "camera_ids": None, "gaussian_ids": None, "depths": None,
                "width": self.width, "height": self.height, "n_cameras": C,
            })
            return

        raise NotImplementedError("Use engine instead. Code below for reference.")
        if self.use_bvh:
            raise NotImplementedError()
            # raise NotImplementedError()
            assert self.packed, "BVH must be packed"
            assert B == 1, "Not support batching"
            from time import perf_counter
            torch.cuda.synchronize()
            time0 = perf_counter()
            intersection_count_map, intersection_splat_id = intersect_splat_tile(
                self.primitive,
                self.splats_world,
                self.width,
                self.height,
                self.viewmats.contiguous(),
                self.intrins.contiguous(),
                self.camera_model,
                self.dist_coeffs,
                1.0 if self.relative_scale is None else self.relative_scale
            )
            torch.cuda.synchronize()
            time1 = perf_counter()
            bvh_time = 1e3*(time1-time0)
            proj_results = fully_fused_projection_hetero(
                self.primitive,
                self.splats_world,
                self.viewmats,
                self.intrins,
                self.actual_width,
                self.actual_height,
                self.width,
                self.height,
                intersection_count_map,
                intersection_splat_id,
                camera_model=self.camera_model,
                dist_coeffs=self.dist_coeffs,
                backward_info=self.backward_info,
            )

        else:
            self.projection_forward()

        self.meta.update(
            {
                # global batch and camera ids
                "camera_ids": self.camera_ids,
                # local gaussian_ids
                "gaussian_ids": self.gaussian_ids,
                # "radii": radii,
                "depths": self.depths,
                # "normals": normals,
            }
        )
        if self.use_bvh:
            raise NotImplementedError()
            self.meta['bvh_time'] = bvh_time
        # if heterogeneous:
        #     meta.update({
        #         "intersection_count": intersection_count_map[1:]-intersection_count_map[:-1]
        #     })

        # Identify intersecting tiles
        TILE_SIZE = 16
        tile_width = math.ceil(self.width / float(TILE_SIZE))
        tile_height = math.ceil(self.height / float(TILE_SIZE))
        intersect_tile_splat_params = None
        if self.primitive in ['3dgs', 'mip']:
            intersect_tile_splat_params = (self.splats_proj[0], self.splats_proj[2], self.splats_proj[3])
        if True:  # Note that GUT is already an approximation
            if self.primitive in ['3dgut', '3dgut_sv']:
                intersect_tile_splat_params = (None, self.splats_proj[0], self.splats_proj[1])
        isect_ids, flatten_ids, isect_offsets = _make_lazy_cuda_func(f"intersect_tile")(
            self.aabb,
            self.depths,
            intersect_tile_splat_params,
            C, self.intrins, self.width, self.height,
            self.camera_ids if self.packed else None,
        )
        self.isect_offsets = isect_offsets.reshape((C, tile_height, tile_width))
        self.flatten_ids = flatten_ids

        self.meta.update(
            {
                "radii": self.radii,
                "tile_width": tile_width,
                "tile_height": tile_height,
                # "tiles_per_gauss": tiles_per_gauss,
                # "isect_ids": isect_ids,
                # "flatten_ids": flatten_ids,
                # "isect_offsets": isect_offsets,
                "width": self.width,
                "height": self.height,
                "tile_size": TILE_SIZE,
                "n_cameras": C,
            }
        )

        self.rasterize_forward()

        # render_colors[1] is depth, map to expected depth
        if len(self.render_colors) > 1:
            self.raw_depth = self.render_colors[1]
            mapped_depth = _make_lazy_cuda_func("rendered_depth_to_expected_depth_forward")(self.raw_depth, self.render_Ts)
            self.render_colors = (
                self.render_colors[0],
                mapped_depth,
                *self.render_colors[2:]
            )

        self.render_colors = self.render_colors
        self.render_Ts = self.render_Ts

    def fused_proj_bwd_optim_step(
        self,
        model_config: 'spirulae_splat.modules.model.SpirulaeSplatModelConfig',
        optim_config: OptimizerConfig,
        step: int,
        max_steps: int
    ):
        raise NotImplementedError("Use engine instead. Code below for reference.")
        if self.primitive not in ["3dgs", "mip", "3dgut", "3dgut_sv"]:
            raise NotImplementedError()

        if optim_config.max_steps is not None:
            max_steps = optim_config.max_steps

        if optim_config.use_per_splat_bias_correction:
            if not hasattr(self, 'optim_bias_correction_step'):
                self.optim_bias_correction_step = torch.ones(
                    self.max_num_splats, dtype=torch.int32, device=self.radii.device)
            else:
                self.optim_bias_correction_step += 1
            bias_correction_step = self.optim_bias_correction_step
        else:
            bias_correction_step = step + 1

        g_optim = [
            [*self.g1_splats_world], [*self.g2_splats_world],
            None,  # sh_packed (AoS joint (u, sqrt_g2) bytes when quantize)
            None,  # sh_quant_bounds
        ]
        if self.quantize_sh_optim:
            # Legacy non-engine path doesn't materialize sh_packed; engine
            # training path uses the C++ QuantizedAdamState class directly.
            g_optim[0][-1] = None
            g_optim[1][-1] = None
            g_optim[2] = self.sh_packed if hasattr(self, "sh_packed") else None
            g_optim[3] = self.quant_bounds_sh

        _make_lazy_cuda_func(f"fused_projection_bwd_optimizer_{self.primitive}")(
            self.cur_num_splats,
            self.splats_world,
            self.viewmats, self.intrins, self.width, self.height, self.camera_model.upper(), self.dist_coeffs,
            self.camera_ids if self.packed else None,
            self.gaussian_ids if self.packed else None,
            self.aabb,
            self.v_splats_world, None, None,
            self.v_splats_proj, None, None,
            *g_optim,
            self.radii,
            optim_config.get_scheduled_lr("means", step, max_steps),
            optim_config.get_scheduled_lr("quats", step, max_steps),
            optim_config.get_scheduled_lr("scales", step, max_steps),
            optim_config.get_scheduled_lr("opacities", step, max_steps),
            optim_config.get_scheduled_lr("features_dc", step, max_steps),
            optim_config.get_scheduled_lr("features_sh", step, max_steps),
            model_config.max_gauss_ratio,
            model_config.scale_regularization_weight,
            model_config.opacity_reg,
            model_config.scale_reg,
            # 0.0, 0.0,  # TODO
            model_config.erank_reg,
            model_config.erank_reg_s3,
            model_config.quat_norm_reg,
            model_config.sh_reg,
            optim_config.use_scale_agnostic_mean,
            bias_correction_step
        )

    def optim_step(
        self,
        model_config: 'spirulae_splat.modules.model.SpirulaeSplatModelConfig',
        optim_config: OptimizerConfig,
        step: int,
        max_steps: int
    ):
        raise NotImplementedError("Use engine instead. Code below for reference.")
        if self.use_fused_proj_bwd_optim:
            self.fused_proj_bwd_optim_step(model_config, optim_config, step, max_steps)
            return

        if self.primitive not in ["3dgs", "mip", "3dgut", "3dgut_sv"]:
            raise NotImplementedError()

        if optim_config.max_steps is not None:
            max_steps = optim_config.max_steps

        if optim_config.use_per_splat_bias_correction:
            if not hasattr(self, 'optim_bias_correction_step'):
                self.optim_bias_correction_step = torch.ones(
                    self.max_num_splats, dtype=torch.int32, device=self.radii.device)
            else:
                self.optim_bias_correction_step += 1
            bias_correction_step = self.optim_bias_correction_step
        else:
            bias_correction_step = step + 1

        # geometry, includes regularization and MCMC add noise
        _make_lazy_cuda_func("fused_optim_3dgs_geometry")(
            self.cur_num_splats,
            self.splats_world[0],
            self.v_splats_world[0],
            self.g1_splats_world[0],
            self.g2_splats_world[0],
            self.splats_world[1],
            self.v_splats_world[1],
            self.g1_splats_world[1],
            self.g2_splats_world[1],
            self.splats_world[2],
            self.v_splats_world[2],
            self.g1_splats_world[2],
            self.g2_splats_world[2],
            self.splats_world[3],
            self.v_splats_world[3],
            self.g1_splats_world[3],
            self.g2_splats_world[3],
            self.radii,
            optim_config.get_scheduled_lr("means", step, max_steps),
            optim_config.get_scheduled_lr("quats", step, max_steps),
            optim_config.get_scheduled_lr("scales", step, max_steps),
            optim_config.get_scheduled_lr("opacities", step, max_steps),
            model_config.max_gauss_ratio,
            model_config.scale_regularization_weight,
            model_config.opacity_reg,
            model_config.scale_reg,
            # 0.0, 0.0,  # TODO
            model_config.erank_reg,
            model_config.erank_reg_s3,
            model_config.quat_norm_reg,
            optim_config.use_scale_agnostic_mean,
            bias_correction_step
        )

        _make_lazy_cuda_func("fused_adam_with_steps")(
            self.cur_num_splats,
            self.splats_world[4],
            self.v_splats_world[4],
            self.g1_splats_world[4],
            self.g2_splats_world[4],
            optim_config.get_scheduled_lr("features_dc", step, max_steps),
            bias_correction_step,
            model_config.sh_reg, 0.5 / 0.28209479177387814
        )

        if self.quantize_sh_optim:
            # TODO: probably better to use structure of array for SH here
            _make_lazy_cuda_func("fused_adam_with_steps_8bit")(
                self.cur_num_splats,
                self.splats_world[5],
                self.v_splats_world[5],
                self.g1_splats_world[5],
                self.g2_splats_world[5],
                self.quant_bounds_sh,
                optim_config.get_scheduled_lr("features_sh", step, max_steps),
                bias_correction_step,
                model_config.sh_reg, 0.0
            )
            # print(self.quant_bounds_sh)
        else:
            _make_lazy_cuda_func("fused_adam_with_steps")(
                self.cur_num_splats,
                self.splats_world[5],
                self.v_splats_world[5],
                self.g1_splats_world[5],
                self.g2_splats_world[5],
                optim_config.get_scheduled_lr("features_sh", step, max_steps),
                bias_correction_step,
                model_config.sh_reg, 0.0
            )

        return
        g1_sh = self.g1_splats_world[5][:self.cur_num_splats] * self.cur_num_splats / (1 - 0.9 ** (step+1))
        g2_sh = self.g2_splats_world[5][:self.cur_num_splats] * self.cur_num_splats**2 / (1 - 0.999 ** (step+1))
        g2_sh **= 0.5
        print(g1_sh.std().item(), g2_sh.mean().item(), g2_sh.std().item())

    def densify_step(
        self,
        step: int,
        max_steps: int,
        model_config: 'spirulae_splat.modules.model.SpirulaeSplatModelConfig',
        optim_config: OptimizerConfig
    ):
        raise NotImplementedError("Use engine instead. Code below for reference.")
        densify_ongoing = (step < max_steps - model_config.refine_stop_num_iter)
        densify_step = densify_ongoing and (step > model_config.refine_start_iter and step % model_config.refine_every == 0)
        use_revised_densification = (model_config.relocate_heuristic_weight >= 1.0)

        # Clip large splats
        progress = (step+0.5) / max_steps
        if np.isfinite(model_config.max_screen_size) or np.isfinite(model_config.max_world_size):
            _make_lazy_cuda_func("densify_clip_scale")(
                self.cur_num_splats,
                self.radii,
                self.splats_world[2],  # scales
                # self.splats_world[3],  # opacs
                None,
                model_config.max_screen_size,
                model_config.max_screen_size_clip_hardness,# ** (progress**2),
                model_config.max_world_size,
            )

        # Update Densification Score

        if densify_ongoing and use_revised_densification:

            if not hasattr(self, 'densify_accum_buffer'):
                self.densify_accum_buffer = torch.zeros(
                    self.max_num_splats, 2, dtype=torch.float32, device=self.radii.device)
            if 'accum_weight' in self.backward_info:
                # is_max_mode = True  # better for clean datasets with uneven coverage
                is_max_mode = False  # better for noisy/high resolution datasets
                densify_score = self.backward_info['accum_weight']
                # densify_scale = self.splats_world[2]
                densify_scale = None
                # densify_opac = None
                densify_opac = self.splats_world[3]
            else:
                is_max_mode = False
                densify_score = self.v_splats_world[3] # v_opacs
                densify_scale = None
                densify_opac = self.splats_world[3]
                if densify_score is None:  # fused backward and optimizer mode
                    densify_score = densify_opac

            _make_lazy_cuda_func("densify_update_weight")(
                self.cur_num_splats,
                self.radii,
                densify_scale,
                densify_opac,
                densify_score,
                self.densify_accum_buffer,
                is_max_mode
            )

        # Apply Densification

        if densify_step:
            torch.cuda.empty_cache()

        if densify_step and use_revised_densification:

            # relocation
            _make_lazy_cuda_func("relocate_splats_with_long_axis_split")(
                self.cur_num_splats,
                model_config.min_opacity,
                *self.splats_world,
                *self.g1_splats_world,
                *self.g2_splats_world,
                self.densify_accum_buffer,
                getattr(self, 'optim_bias_correction_step', None),
                2*step+0
            )

            # add more splats
            n_target = min(self.max_num_splats, int(model_config.growth_factor * self.cur_num_splats))
            num_add = max(0, n_target - self.cur_num_splats)
            _make_lazy_cuda_func("add_splats_with_long_axis_split")(
                self.cur_num_splats,
                num_add,
                *self.splats_world,
                *self.g1_splats_world,
                *self.g2_splats_world,
                self.densify_accum_buffer,
                getattr(self, 'optim_bias_correction_step', None),
                2*step+1
            )
            self.cur_num_splats += num_add

        elif densify_step:

            # mcmc relocation
            _make_lazy_cuda_func("relocate_splats_mcmc")(
                self.cur_num_splats,
                model_config.min_opacity,
                *self.splats_world,
                *self.g1_splats_world,
                *self.g2_splats_world,
                getattr(self, 'optim_bias_correction_step', None),
                2*step+0
            )

            # mcmc sample add
            n_target = min(self.max_num_splats, int(model_config.growth_factor * self.cur_num_splats))
            num_add = max(0, n_target - self.cur_num_splats)
            _make_lazy_cuda_func("add_splats_mcmc")(
                self.cur_num_splats,
                num_add,
                model_config.min_opacity,
                *self.splats_world,
                *self.g1_splats_world,
                *self.g2_splats_world,
                getattr(self, 'optim_bias_correction_step', None),
                2*step+1
            )
            self.cur_num_splats += num_add

        if densify_step:
            if hasattr(self, 'densify_accum_buffer'):
                _make_lazy_cuda_func("set_zero")(self.densify_accum_buffer)

            torch.cuda.empty_cache()

        # Add MCMC noise

        if model_config.opacity_decay != 0.0 or model_config.scale_decay != 0.0:
            raise NotImplementedError()

        if model_config.noise_lr > 0.0 and model_config.noise_lr_final > 0.0:
            noise_scalar = model_config.noise_lr * (model_config.noise_lr_final / model_config.noise_lr) ** progress

            if model_config.relocate_heuristic_weight >= 1.0:
                _make_lazy_cuda_func("revised_add_noise")(
                    self.cur_num_splats,
                    noise_scalar,
                    self.radii,
                    self.splats_world[0],  # means
                    self.splats_world[2],  # quats
                    self.splats_world[1],  # scales
                    self.splats_world[3],  # opacs
                )
            else:
                _make_lazy_cuda_func("mcmc_add_noise")(
                    self.cur_num_splats,
                    noise_scalar,
                    self.splats_world[0],  # means
                    self.splats_world[2],  # quats
                    self.splats_world[1],  # scales
                    self.splats_world[3],  # opacs
                )
