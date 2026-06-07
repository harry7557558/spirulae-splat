import torch
import torch.nn.functional as F
from torch import Tensor
from typing import Dict, Optional, Tuple, Literal
import math
import numpy as np


import spirulae_splat

from spirulae_splat.splat.cuda import (
    _C,
    _make_lazy_cuda_func,
)

from spirulae_splat.modules.optimizer import OptimizerConfig


def engine_reset():
    """Wipe the C++ EngineState singleton and free all engine-pool / scratch
    VRAM. Must be called when swapping training datasets in the same Python
    process (e.g. ss_benchmark looping over scenes); without it the new run
    inherits the previous run's world splats, camera table, bilagrid / PPISP /
    background config, optimizer moments, and color-space matrices."""
    _C.engine_reset()


class Renderer:

    def __init__(
        self,
        primitive: Literal["3dgs", "mip", "3dgut"],
        splats_world: tuple[Tensor],  # means, quats, scales, opacities
        cur_num_splats: int,
        packed: bool = True,
        use_bvh: bool = False,
        use_fused_proj_bwd_optim: bool = False,
        sh_quantization_level: int = 0,
    ):
        for tensor in splats_world:
            assert tensor.is_contiguous(), "Tensor must be contiguous"

        features_dc, features_sh = None, None
        if primitive in ["3dgs", "mip", "3dgut"]:
            assert len(splats_world) == 6, "3DGS requires 6 params (means, quats, scales, opacities, color, sh)"
            means, quats, scales, opacities, features_dc, features_sh = splats_world
            assert len(means.shape) == 2, "means must have at 2 dimensions"
            N = means.shape[-2]
            device = means.device
            assert means.shape == (N, 3), means.shape
            assert quats.shape == (N, 4), quats.shape
            assert scales.shape == (N, 3), scales.shape
            assert opacities.shape == (N, 1), opacities.shape
        else:
            raise ValueError(f"Invalid primitive ({primitive})")
        if features_dc is not None and features_sh is not None:
            assert features_dc.shape == (N, 3), features_dc.shape
            assert features_sh.shape == (N, 0, 3) or \
                features_sh.shape == (N, 3, 3) or \
                features_sh.shape == (N, 8, 3) or \
                features_sh.shape == (N, 15, 3) or \
                features_sh.shape == (N, 24, 3), features_sh.shape

        self.device = device
        self.cur_num_splats = cur_num_splats
        self.max_num_splats = N

        self.primitive = primitive
        self.splats_world = splats_world


        self.packed = packed
        self.use_bvh = use_bvh
        self.use_fused_proj_bwd_optim = use_fused_proj_bwd_optim
        # Single SH quantization level. 0=off, 1=light (16-bit value + 8-bit
        # optim), 2=heavy (8-bit value + 8-bit optim). Collapses the FPBO
        # dispatcher's bit-pair axes into a single template parameter so the
        # generated kernel-instantiation .cu files distribute work evenly and
        # there are 3x fewer wrappers compared to enumerating each bit-pair.
        if sh_quantization_level not in (0, 1, 2):
            raise ValueError(
                f"sh_quantization_level must be 0 (off), 1 (light), or 2 (heavy); "
                f"got {sh_quantization_level}")
        self.sh_quantization_level = sh_quantization_level
        # Derived per-component bits. Used by non-FPBO paths
        # (fused_adam_step_quantized, EngineState) and by EngineCheckpoint
        # metadata. The FPBO dispatcher uses the level directly.
        _level_to_optim_bits = {0: 32, 1: 8, 2: 8}
        _level_to_value_bits = {0: 32, 1: 16, 2: 8}
        self.sh_optim_bits = _level_to_optim_bits[sh_quantization_level]
        self.sh_value_bits = _level_to_value_bits[sh_quantization_level]

        # Set by engine_init_color_space; consulted by the forward path so it
        # can skip allocating a separate rgb_raw host buffer when no color
        # space is active. Default is False (no color space).
        self._has_engine_color_space = False

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
        # rgb_raw: pre-color-space-conversion render (linear / wide-gamut).
        # The C++ side writes into this buffer only when an engine color
        # space is configured; otherwise leave rgb_raw aliased to rgb so we
        # avoid a redundant D->H of the identical buffer.
        rgb_raw = torch.empty(C, H, W, 3, dtype=torch.float32) \
            if self._has_engine_color_space else rgb
        _C.engine_copy_render_to_host(
            self._tv(rgb), self._tv(depth), self._tv(Ts), self._tv(rgb_raw)
        )

        self.render_colors = (rgb, depth)
        self.render_Ts = Ts
        self.render_rgb_raw = rgb_raw

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
        return buf[:, 0]

    def engine_set_training_data(self, gt_rgb, gt_depth=None, gt_normal=None,
                                  gt_alpha=None):
        # 4-mask buffers (rgb/depth/normal/alpha) gone — derived in the slang
        # kernel from gt_alpha / gt_depth=0 / gt_normal sentinel.
        _C.set_training_data(
            self._tv(gt_rgb), self._tv(gt_depth), self._tv(gt_normal), self._tv(gt_alpha)
        )

    def engine_compute_loss_backward(self, step, loss_weights, w_ssim,
                                     num_loss_scales, compute_loss_map,
                                     structure_only_loss_map,
                                     overexposure_reg_weight=0.0):
        """Compute loss + rasterization backward + projection backward in C++.
        Gradients are managed by C++ pool."""
        loss_dict = _C.engine_compute_loss_backward(
            step, loss_weights, w_ssim, num_loss_scales, compute_loss_map,
            bool(structure_only_loss_map),
            float(overexposure_reg_weight),
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
        c.sh_quantization_level          = self.sh_quantization_level
        c.sh_optim_bits                  = self.sh_optim_bits
        c.sh_value_bits                  = self.sh_value_bits
        c.use_per_splat_bias_correction  = optim_config.use_per_splat_bias_correction
        c.use_fused_proj_bwd_optim       = self.use_fused_proj_bwd_optim
        c.color_is_linear = model_config.splat_color_is_linear
        c.use_color_trust_region = model_config.splat_color_is_linear
        c.eps_tr = 1e-6 * 0.01 ** (step / max_steps)  # TODO: make this configurable
        return c

    def _build_densify_config(self, model_config):
        alpha = self.train_frame_scale
        noise_lr_scalar = 1.0 if model_config.use_revised_densification else alpha
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
        c.use_revised_densification     = model_config.use_revised_densification
        c.score_mode                    = {
            "mean":   0,
            "max":    1,
            "median": 2,
            "geom":   3,
        }[str(getattr(model_config, "densify_score_mode", "mean"))]
        return c

    def engine_optim_step(self, step, max_steps, model_config, optim_config):
        """Run optimizer step via C++ Engine. All tensors managed by C++ pool."""
        _C.engine_optim_step(
            step,
            self._build_optim_config(step, max_steps, model_config, optim_config),
        )

    def engine_init_bilagrid(self, n_grids, rgb_type=None, rgb_LHW=None,
                              depth_LHW=None, normal_LHW=None,
                              rgb_optim_bits=32,
                              geometry_optim_bits=32,
                              use_adagrad=False):
        """One-time allocation of C++ bilagrid storage. Pass `None` for any type
        to leave it disabled. Must be called after set_data_3dgs (i.e., after
        engine_forward has been invoked at least once).

        rgb_optim_bits / geometry_optim_bits: bit depth for optimizer-state
        quantization. 32 = no quantization (fp32 g1/g2 for Adam, fp32 accum for
        AdaGrad). 4 or 8 = packed QuantizedAdamState<BITS, 256> (Adam) or
        QuantizedTensorLog<BITS, 256> (AdaGrad). The RGB grid is much larger
        than the geometry grids in photogrammetry-scale datasets, so the depth
        is configurable per-axis.

        use_adagrad: select AdaGrad over Adam for this bilagrid. Mirrored to
        all three sub-types (RGB / depth / normal) in this call — bilagrids
        share the optimizer choice on the Python side."""
        for label, b in (("rgb", rgb_optim_bits), ("geometry", geometry_optim_bits)):
            if b not in (32, 8, 4):
                raise ValueError(
                    f"engine_init_bilagrid: {label}_optim_bits must be 32 (off), "
                    f"4, or 8; got {b}")
        if rgb_type is not None and rgb_LHW is not None:
            L, H, W = rgb_LHW
            _C.engine_init_bilagrid_rgb(n_grids, rgb_type, L, H, W,
                                        int(rgb_optim_bits),
                                        bool(use_adagrad))
        if depth_LHW is not None:
            L, H, W = depth_LHW
            _C.engine_init_bilagrid_depth(n_grids, L, H, W,
                                          int(geometry_optim_bits),
                                          bool(use_adagrad))
        if normal_LHW is not None:
            L, H, W = normal_LHW
            _C.engine_init_bilagrid_normal(n_grids, L, H, W,
                                           int(geometry_optim_bits),
                                           bool(use_adagrad))

    def engine_bilagrid_forward(self, cam_indices):
        """Apply bilagrid forward for the current batch. cam_indices: [C_batch]
        int32 tensor of per-image camera-table indices (or None for identity).
        In-place on rendered RGB and on the GT depth/normal buffers."""
        _C.engine_bilagrid_forward(self._tv(cam_indices))

    @staticmethod
    def _build_bilagrid_step_config(lr_rgb, lr_depth, lr_normal,
                                     tv_weight_rgb, tv_weight_depth, tv_weight_normal,
                                     shift_reg_weight_rgb=0.0,
                                     shift_reg_beta_rgb=0.0):
        c = _C.BilagridStepConfig()
        c.lr_rgb           = float(lr_rgb)
        c.lr_depth         = float(lr_depth)
        c.lr_normal        = float(lr_normal)
        c.tv_weight_rgb    = float(tv_weight_rgb)
        c.tv_weight_depth  = float(tv_weight_depth)
        c.tv_weight_normal = float(tv_weight_normal)
        c.shift_reg_weight_rgb = float(shift_reg_weight_rgb)
        c.shift_reg_beta_rgb   = float(shift_reg_beta_rgb)
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

    def engine_init_background_noise(self, splat_color_is_linear):
        _C.engine_init_background_noise(bool(splat_color_is_linear))

    def engine_init_background_sh(self, sh_degree, splat_color_is_linear):
        _C.engine_init_background_sh(int(sh_degree), bool(splat_color_is_linear))

    def engine_init_color_space(
        self,
        splat_enabled, splat_is_linear, splat_color_matrix,
        image_enabled, image_is_linear, image_color_matrix,
    ):
        # Flatten to row-major float lists. Empty lists when the side is
        # disabled — the C++ side ignores the matrix in that case.
        def _flat(m):
            if m is None: return []
            t = m.detach().contiguous().cpu().float()
            return t.flatten().tolist()
        _C.engine_init_color_space(
            bool(splat_enabled), bool(splat_is_linear), _flat(splat_color_matrix),
            bool(image_enabled), bool(image_is_linear), _flat(image_color_matrix),
        )
        # Cached so the forward path can decide whether to allocate a
        # separate rgb_raw host buffer (avoids a redundant D->H copy when
        # the engine has no color space).
        self._has_engine_color_space = bool(splat_enabled)

    @staticmethod
    def _build_background_step_config(lr_dc, lr_sh, randomize_weight, seed):
        c = _C.BackgroundStepConfig()
        c.lr_dc            = float(lr_dc)
        c.lr_sh            = float(lr_sh)
        c.randomize_weight = float(randomize_weight)
        c.seed             = int(seed) & 0xFFFFFFFF
        return c

    def engine_background_optim_step(self, step, lr_dc, lr_sh):
        """Adam step over background SH coefficients. No-op for noise mode or
        when train_color was False at init."""
        cfg = self._build_background_step_config(lr_dc, lr_sh, 0.0, 0)
        _C.engine_background_optim_step(int(step), cfg)

    def engine_copy_background_to_host(self, out_image):
        """Fill a host (..., H, W, 3) float buffer with the engine's current
        background image (skybox for SH mode, uniform mean color for noise
        mode). Returns True if the engine has an active background mode."""
        return bool(_C.engine_copy_background_to_host(self._tv(out_image)))

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
                          structure_only_loss_map,
                          # Configs
                          model_config, optim_config,
                          # Bilagrid (pass cam_indices + lrs + tv weights;
                          # ignored if engine_init_bilagrid_* was never called).
                          # cam_indices: [C_batch] int32 tensor, or None.
                          bilagrid_cam_indices=None,
                          bilagrid_lr_rgb=0.0, bilagrid_lr_depth=0.0, bilagrid_lr_normal=0.0,
                          bilagrid_tv_weight_rgb=0.0, bilagrid_tv_weight_depth=0.0,
                          bilagrid_tv_weight_normal=0.0,
                          # RGB bilagrid color-shift regularizer (design 1).
                          # 0 disables. beta is the EMA decay.
                          bilagrid_shift_reg_weight_rgb=0.0,
                          bilagrid_shift_reg_beta_rgb=0.0,
                          # PPISP (ignored if engine_init_ppisp was never called).
                          # Reuses bilagrid_cam_idx for the camera selector.
                          ppisp_lr=0.0,
                          ppisp_reg_exposure_mean=0.0, ppisp_reg_vig_center=0.0,
                          ppisp_reg_vig_non_pos=0.0, ppisp_reg_vig_channel_var=0.0,
                          ppisp_reg_color_mean=0.0, ppisp_reg_crf_channel_var=0.0,
                          # Background (ignored if engine_init_background_* was
                          # never called).
                          bg_lr_dc=0.0, bg_lr_sh=0.0,
                          bg_randomize_weight=0.0, bg_seed=0,
                          # Image-space overexposure regularization (model.py
                          # config field). Zero -> kernel not launched.
                          overexposure_reg_weight=0.0):
        """Single fused training step (set_camera + set_gt + fwd + loss/bwd + optim + densify).
        All input tensors are CPU; returns loss_dict for verbose."""
        max_steps_lr = optim_config.max_steps if optim_config.max_steps is not None else max_steps

        cfg = _C.EngineStepConfig()
        cfg.loss.weights          = loss_weights
        cfg.loss.w_ssim           = float(w_ssim)
        cfg.loss.num_loss_scales  = int(num_loss_scales)
        cfg.loss.compute_loss_map = bool(compute_loss_map)
        cfg.loss.structure_only_loss_map = bool(structure_only_loss_map)
        cfg.loss.overexposure_reg_weight = float(overexposure_reg_weight)
        cfg.optim    = self._build_optim_config(step, max_steps_lr, model_config, optim_config)
        cfg.densify  = self._build_densify_config(model_config)
        cfg.bilagrid = self._build_bilagrid_step_config(
            bilagrid_lr_rgb, bilagrid_lr_depth, bilagrid_lr_normal,
            bilagrid_tv_weight_rgb, bilagrid_tv_weight_depth, bilagrid_tv_weight_normal,
            shift_reg_weight_rgb=bilagrid_shift_reg_weight_rgb,
            shift_reg_beta_rgb=bilagrid_shift_reg_beta_rgb)
        cfg.ppisp    = self._build_ppisp_step_config(
            ppisp_lr, ppisp_reg_exposure_mean, ppisp_reg_vig_center,
            ppisp_reg_vig_non_pos, ppisp_reg_vig_channel_var,
            ppisp_reg_color_mean, ppisp_reg_crf_channel_var)
        cfg.background = self._build_background_step_config(
            bg_lr_dc, bg_lr_sh, bg_randomize_weight, bg_seed)

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

    def engine_train_step_managed(self, step, max_steps,
                                  sh_degree_to_use,
                                  loss_weights, w_ssim, num_loss_scales,
                                  compute_loss_map, structure_only_loss_map,
                                  model_config, optim_config,
                                  bilagrid_lr_rgb=0.0, bilagrid_lr_depth=0.0,
                                  bilagrid_lr_normal=0.0,
                                  bilagrid_tv_weight_rgb=0.0,
                                  bilagrid_tv_weight_depth=0.0,
                                  bilagrid_tv_weight_normal=0.0,
                                  bilagrid_shift_reg_weight_rgb=0.0,
                                  bilagrid_shift_reg_beta_rgb=0.0,
                                  ppisp_lr=0.0,
                                  ppisp_reg_exposure_mean=0.0, ppisp_reg_vig_center=0.0,
                                  ppisp_reg_vig_non_pos=0.0, ppisp_reg_vig_channel_var=0.0,
                                  ppisp_reg_color_mean=0.0, ppisp_reg_crf_channel_var=0.0,
                                  bg_lr_dc=0.0, bg_lr_sh=0.0,
                                  bg_randomize_weight=0.0, bg_seed=0,
                                  overexposure_reg_weight=0.0):
        """DataManager-driven fused training step.

        Same EngineStepConfig as ``engine_train_step`` but no per-batch tensors —
        the C++ engine pulls the next batch from its installed DataManager
        (see ``_C.engine_setup_data_manager``) and dispatches the fused
        forward + loss/bwd + optim + densify pipeline. Returns the same
        loss_dict shape."""
        max_steps_lr = optim_config.max_steps if optim_config.max_steps is not None else max_steps

        cfg = _C.EngineStepConfig()
        cfg.loss.weights          = loss_weights
        cfg.loss.w_ssim           = float(w_ssim)
        cfg.loss.num_loss_scales  = int(num_loss_scales)
        cfg.loss.compute_loss_map = bool(compute_loss_map)
        cfg.loss.structure_only_loss_map = bool(structure_only_loss_map)
        cfg.loss.overexposure_reg_weight = float(overexposure_reg_weight)
        cfg.optim    = self._build_optim_config(step, max_steps_lr, model_config, optim_config)
        cfg.densify  = self._build_densify_config(model_config)
        cfg.bilagrid = self._build_bilagrid_step_config(
            bilagrid_lr_rgb, bilagrid_lr_depth, bilagrid_lr_normal,
            bilagrid_tv_weight_rgb, bilagrid_tv_weight_depth, bilagrid_tv_weight_normal,
            shift_reg_weight_rgb=bilagrid_shift_reg_weight_rgb,
            shift_reg_beta_rgb=bilagrid_shift_reg_beta_rgb)
        cfg.ppisp    = self._build_ppisp_step_config(
            ppisp_lr, ppisp_reg_exposure_mean, ppisp_reg_vig_center,
            ppisp_reg_vig_non_pos, ppisp_reg_vig_channel_var,
            ppisp_reg_color_mean, ppisp_reg_crf_channel_var)
        cfg.background = self._build_background_step_config(
            bg_lr_dc, bg_lr_sh, bg_randomize_weight, bg_seed)

        result = _C.engine_train_step_managed(
            step, max_steps,
            self.primitive, sh_degree_to_use, self.packed,
            cfg,
        )
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

    def projection_forward(self):
        raise NotImplementedError("Use engine instead. Code below for reference.")
        if self.primitive not in ["3dgs", "mip", "3dgut"]:
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
            if self.primitive in ['3dgut']:
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
