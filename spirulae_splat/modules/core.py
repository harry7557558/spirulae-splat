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
        split_batch: bool = False,
        quantization_level: int = 0,
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
        # When True, engine_train_step splits the camera batch into one-camera
        # sub-batches, atomic-adds per-splat grads across them, and runs a
        # single optim+densify pass at the end (with grad_scale = 1/B inside
        # the splat Adam kernels). Reduces peak VRAM by ~1/B for the
        # immediate projection / rasterization buffers. Not compatible with
        # use_fused_proj_bwd_optim or the warped train-step path.
        self.split_batch = split_batch
        # Single SH quantization level. 0=off, 1=light (16-bit value + 8-bit
        # optim). Collapses the FPBO dispatcher's bit-pair axes into a single
        # template parameter so the generated kernel-instantiation .cu files
        # distribute work evenly.
        if quantization_level not in (0, 1):
            raise ValueError(
                f"quantization_level must be 0 (off) or 1 (light); "
                f"got {quantization_level}")
        self.quantization_level = quantization_level
        # Derived per-component bits. Used by non-FPBO paths
        # (fused_adam_step_quantized, EngineState) and by EngineCheckpoint
        # metadata. The FPBO dispatcher uses the level directly.
        _level_to_optim_bits = {0: 32, 1: 8}
        _level_to_value_bits = {0: 32, 1: 16}
        # Non-SH Adam state (means, quats, scales, opacities, features_dc).
        # Level 1 quantizes to 16-bit (joint u/log_s = 4 B/cell); FPBO-only.
        _level_to_non_sh_optim_bits = {0: 32, 1: 16}
        self.sh_optim_bits     = _level_to_optim_bits[quantization_level]
        self.sh_value_bits     = _level_to_value_bits[quantization_level]
        self.non_sh_optim_bits = _level_to_non_sh_optim_bits[quantization_level]

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
        output_median: bool = False,
        compute_hessian_diagonal: Literal[None, "position", "all"] = None,
        relative_scale: Optional[float] = None,
        accum_weight_map: Optional[Tensor] = None,
        max_blending_masks: Optional[Tensor] = None,
        camera_model: Literal["pinhole", "ortho", "fisheye", "equisolid", "equirectangular"] = "pinhole",
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
        self.output_median = output_median
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
        output_median = getattr(self, "output_median", False)
        # DistortionType (mirror Primitive.cuh): 0=None, 1=D, 2=RGB_D, 3=DN, 4=RGB_DN.
        # For eval/viewer request the full set the primitive renders (RGB_D for
        # all current primitives) so any distortion channel is viewable.
        output_distortion = getattr(self, "output_distortion", False) \
            and self.primitive in ("3dgut", "3dgs", "mip")
        dist_type = 2 if output_distortion else 0  # RGB_D when on, else None
        _C.engine_forward_3dgs(
            self.primitive,
            self.sh_degree_to_use,
            self.packed,
            output_median,
            dist_type,
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
        median = torch.empty(C, H, W, 1, dtype=torch.float32) \
            if output_median else None
        _C.engine_copy_render_to_host(
            self._tv(rgb), self._tv(depth), self._tv(Ts), self._tv(rgb_raw),
            self._tv(median) if median is not None else (0, 0, [])
        )

        self.render_colors = (rgb, depth)
        self.render_Ts = Ts
        self.render_rgb_raw = rgb_raw
        self.render_median = median

        # Distortion image D = W*S - C^2 (per channel). RGB_D primitives emit
        # rgb + depth only. Stored as (rgb_dist [C,H,W,3], depth_dist [C,H,W,1]).
        if output_distortion:
            rgb_dist = torch.empty(C, H, W, 3, dtype=torch.float32)
            depth_dist = torch.empty(C, H, W, 1, dtype=torch.float32)
            _C.engine_copy_distortion_to_host(
                self._tv(rgb_dist), self._tv(depth_dist))
            self.render_distortion = (rgb_dist, depth_dist)
        else:
            self.render_distortion = None

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
                                     loss_map_mode,
                                     robust_edge_aware_quantile,
                                     overexposure_reg_weight=0.0,
                                     color_shift_reg_weight=0.0,
                                     color_shift_reg_beta=0.0,
                                     loss_scale_min_pixels=0):
        """Compute loss + rasterization backward + projection backward in C++.
        Gradients are managed by C++ pool."""
        loss_dict = _C.engine_compute_loss_backward(
            step, loss_weights, w_ssim, num_loss_scales,
            int(loss_scale_min_pixels), compute_loss_map,
            int(loss_map_mode),
            float(robust_edge_aware_quantile),
            float(overexposure_reg_weight),
            float(color_shift_reg_weight),
            float(color_shift_reg_beta),
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
        c.quantization_level          = self.quantization_level
        c.sh_optim_bits                  = self.sh_optim_bits
        c.sh_value_bits                  = self.sh_value_bits
        c.non_sh_optim_bits              = self.non_sh_optim_bits
        c.use_per_splat_bias_correction  = optim_config.use_per_splat_bias_correction
        c.use_fused_proj_bwd_optim       = self.use_fused_proj_bwd_optim
        c.write_densify_world_grad_score = (
            float(getattr(model_config, "densify_score_blend_world_grad", 0.0)) > 0.0
            and model_config.use_revised_densification)
        c.split_batch     = self.split_batch
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
        c.score_blend_world_grad        = float(getattr(
            model_config, "densify_score_blend_world_grad", 0.0))
        k_init, k_final, k_warmup = model_config.long_axis_split_opacity_k
        c.las_split_opacity_k_init      = float(k_init)
        c.las_split_opacity_k_final     = float(k_final)
        c.las_split_opacity_k_warmup    = int(k_warmup)
        return c

    def engine_optim_step(self, step, max_steps, model_config, optim_config):
        """Run optimizer step via C++ Engine. All tensors managed by C++ pool."""
        _C.engine_optim_step(
            step,
            self._build_optim_config(step, max_steps, model_config, optim_config),
        )

    def engine_init_bilagrid(self, n_grids, rgb_type=None, rgb_LHW=None,
                              depth_LHW=None, normal_LHW=None,
                              use_adagrad=False):
        """One-time allocation of C++ bilagrid storage. Pass `None` for any type
        to leave it disabled. Must be called after set_data_3dgs (i.e., after
        engine_forward has been invoked at least once).

        Bilagrid bit depths are coupled to `self.quantization_level`:
            level 0: 32-bit grids + fp32 optim state.
            level 1: 16-bit linearly-quantized grids (QuantizedTensor<16, 256>
                     packed bytes + per-256-cell float2 bounds) AND 8x2-bit
                     packed optim state (QuantizedAdamState<8, 256> / Adam,
                     QuantizedTensorLog<8, 256> / AdaGrad).
        Both axes apply uniformly across RGB / depth / normal -- the bilagrid
        path is small enough that per-type tuning doesn't pay off.

        use_adagrad: select AdaGrad over Adam for this bilagrid. Mirrored to
        all three sub-types (RGB / depth / normal) in this call -- bilagrids
        share the optimizer choice on the Python side."""
        if self.quantization_level == 0:
            optim_bits = 32
            value_bits = 32
        else:
            optim_bits = 8
            value_bits = 16
        if rgb_type is not None and rgb_LHW is not None:
            L, H, W = rgb_LHW
            _C.engine_init_bilagrid_rgb(n_grids, rgb_type, L, H, W,
                                        int(optim_bits),
                                        int(value_bits),
                                        bool(use_adagrad))
        if depth_LHW is not None:
            L, H, W = depth_LHW
            _C.engine_init_bilagrid_depth(n_grids, L, H, W,
                                          int(optim_bits),
                                          int(value_bits),
                                          bool(use_adagrad))
        if normal_LHW is not None:
            L, H, W = normal_LHW
            _C.engine_init_bilagrid_normal(n_grids, L, H, W,
                                           int(optim_bits),
                                           int(value_bits),
                                           bool(use_adagrad))

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

    def engine_init_ppisp(self, n_grids, param_type="original",
                           use_adagrad=False):
        """One-time allocation of C++ PPISP per-camera parameter table. Seeded
        with the type's default values. Must be called after set_data_3dgs.

        use_adagrad: select unscheduled AdaGrad (lr_decay=0, weight_decay=0,
        initial_accumulator_value=0, eps=1e-15) over Adam. PPISP has no
        quantization path -- state stays fp32 either way."""
        _C.engine_init_ppisp(int(n_grids), str(param_type), bool(use_adagrad))

    def engine_ppisp_forward(self, cam_indices):
        """Apply PPISP forward in place on the current rendered RGB.
        cam_indices: [C_batch] int32 tensor (or None for identity)."""
        _C.engine_ppisp_forward(self._tv(cam_indices))

    @staticmethod
    def _build_ppisp_step_config(lr,
                                 reg_exposure_mean, reg_vig_center,
                                 reg_vig_non_pos, reg_vig_channel_var,
                                 reg_color_mean, reg_crf_channel_var,
                                 run_before_bilagrid=False):
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
        # Forward order: False -> bilagrid then PPISP (default); True -> PPISP
        # then bilagrid. Backward order is auto-inverted on the C++ side.
        c.run_before_bilagrid = bool(run_before_bilagrid)
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
        Always writes ``splat.ply`` (cur_num_splats only, NaN/low-opacity-filtered)
        and ``state.tar`` -- a POSIX tar of ``state.json`` (runtime/validation
        manifest) plus one flat typed ``.npy`` per saved pool buffer, named by its
        DevicePool slot. ``full_dump=False`` saves the Always slots (appearance /
        inference params); ``full_dump=True`` also saves the Resume slots (world
        raw params + all optimizer state) needed to resume training. The buffer
        set is selected by SaveClass metadata, so it stays in sync automatically."""
        import os
        os.makedirs(str(output_dir), exist_ok=True)
        _C.engine_save_checkpoint(str(output_dir), bool(full_dump), int(step))

    def engine_load_checkpoint(self, input_dir):
        """Restore engine state from ``input_dir``/state.tar (resume training).

        Preconditions (established by the caller): the engine skeleton is
        configured -- world seeded/allocated at max_num_splats (done by
        ``set_data_3dgs`` during Renderer construction), and any bilagrid / PPISP
        channels present in the checkpoint already initialized. This restores the
        world raw params, optimizer state, densify aux, and appearance state from
        the Resume-class .npy payloads and returns the saved training ``step``.
        The engine's cur_num_splats is set by the loader; the caller mirrors it
        onto ``self.cur_num_splats`` from the checkpoint's state.json. Requires a
        ``save_full_checkpoint=True`` dump."""
        return int(_C.engine_load_checkpoint(str(input_dir)))

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
                          loss_map_mode, robust_edge_aware_quantile,
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
                          ppisp_reg_color_mean=0.0, ppisp_reg_crf_channel_var=0.0,
                          # When True, PPISP forward runs BEFORE bilagrid (and
                          # PPISP bwd runs AFTER bilagrid bwd). Default False:
                          # bilagrid -> PPISP. Ignored when only one is enabled.
                          apply_ppisp_before_bilagrid=False,
                          # Background (ignored if engine_init_background_* was
                          # never called).
                          bg_lr_dc=0.0, bg_lr_sh=0.0,
                          bg_randomize_weight=0.0, bg_seed=0,
                          # Image-space overexposure regularization (model.py
                          # config field). Zero -> kernel not launched.
                          overexposure_reg_weight=0.0,
                          # Combined bilagrid + PPISP color-shift regularizer
                          # (design 1). Active when at least one of bilagrid_rgb
                          # / PPISP is enabled. 0 disables.
                          color_shift_reg_weight=0.0,
                          color_shift_reg_beta=0.0):
        """Single fused training step (set_camera + set_gt + fwd + loss/bwd + optim + densify).
        All input tensors are CPU; returns loss_dict for verbose."""
        max_steps_lr = optim_config.max_steps if optim_config.max_steps is not None else max_steps

        cfg = _C.EngineStepConfig()
        cfg.loss.weights          = loss_weights
        cfg.loss.w_ssim           = float(w_ssim)
        cfg.loss.num_loss_scales  = int(num_loss_scales)
        cfg.loss.loss_scale_min_pixels = int(getattr(model_config, "loss_scale_min_pixels", 0))
        cfg.loss.compute_loss_map = bool(compute_loss_map)
        cfg.loss.loss_map_mode = int(loss_map_mode)
        cfg.loss.robust_edge_aware_quantile = float(robust_edge_aware_quantile)
        cfg.loss.overexposure_reg_weight = float(overexposure_reg_weight)
        cfg.loss.color_shift_reg_weight = float(color_shift_reg_weight)
        cfg.loss.color_shift_reg_beta   = float(color_shift_reg_beta)
        # When False, GT depth is linear (z) depth and is converted to ray
        # depth in place on upload (set_training_data), to match the ray depth
        # the rasterizer renders.
        cfg.loss.input_depth_is_ray_depth = bool(model_config.input_depth_is_ray_depth)
        cfg.optim    = self._build_optim_config(step, max_steps_lr, model_config, optim_config)
        cfg.densify  = self._build_densify_config(model_config)
        cfg.bilagrid = self._build_bilagrid_step_config(
            bilagrid_lr_rgb, bilagrid_lr_depth, bilagrid_lr_normal,
            bilagrid_tv_weight_rgb, bilagrid_tv_weight_depth, bilagrid_tv_weight_normal)
        cfg.ppisp    = self._build_ppisp_step_config(
            ppisp_lr, ppisp_reg_exposure_mean, ppisp_reg_vig_center,
            ppisp_reg_vig_non_pos, ppisp_reg_vig_channel_var,
            ppisp_reg_color_mean, ppisp_reg_crf_channel_var,
            run_before_bilagrid=apply_ppisp_before_bilagrid)
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
                                  compute_loss_map, loss_map_mode,
                                  robust_edge_aware_quantile,
                                  model_config, optim_config,
                                  bilagrid_lr_rgb=0.0, bilagrid_lr_depth=0.0,
                                  bilagrid_lr_normal=0.0,
                                  bilagrid_tv_weight_rgb=0.0,
                                  bilagrid_tv_weight_depth=0.0,
                                  bilagrid_tv_weight_normal=0.0,
                                  ppisp_lr=0.0,
                                  ppisp_reg_exposure_mean=0.0, ppisp_reg_vig_center=0.0,
                                  ppisp_reg_vig_non_pos=0.0, ppisp_reg_vig_channel_var=0.0,
                                  ppisp_reg_color_mean=0.0, ppisp_reg_crf_channel_var=0.0,
                                  apply_ppisp_before_bilagrid=False,
                                  bg_lr_dc=0.0, bg_lr_sh=0.0,
                                  bg_randomize_weight=0.0, bg_seed=0,
                                  overexposure_reg_weight=0.0,
                                  color_shift_reg_weight=0.0,
                                  color_shift_reg_beta=0.0):
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
        cfg.loss.loss_scale_min_pixels = int(getattr(model_config, "loss_scale_min_pixels", 0))
        cfg.loss.compute_loss_map = bool(compute_loss_map)
        cfg.loss.loss_map_mode = int(loss_map_mode)
        cfg.loss.robust_edge_aware_quantile = float(robust_edge_aware_quantile)
        cfg.loss.overexposure_reg_weight = float(overexposure_reg_weight)
        cfg.loss.color_shift_reg_weight = float(color_shift_reg_weight)
        cfg.loss.color_shift_reg_beta   = float(color_shift_reg_beta)
        cfg.optim    = self._build_optim_config(step, max_steps_lr, model_config, optim_config)
        cfg.densify  = self._build_densify_config(model_config)
        cfg.bilagrid = self._build_bilagrid_step_config(
            bilagrid_lr_rgb, bilagrid_lr_depth, bilagrid_lr_normal,
            bilagrid_tv_weight_rgb, bilagrid_tv_weight_depth, bilagrid_tv_weight_normal)
        cfg.ppisp    = self._build_ppisp_step_config(
            ppisp_lr, ppisp_reg_exposure_mean, ppisp_reg_vig_center,
            ppisp_reg_vig_non_pos, ppisp_reg_vig_channel_var,
            ppisp_reg_color_mean, ppisp_reg_crf_channel_var,
            run_before_bilagrid=apply_ppisp_before_bilagrid)
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
        # Engine backward zeroes engine().grad.* at the start of every call
        # (_alloc_grad_buffers), so there is nothing to clear here.
        return

    def backward(self, v_render_colors, v_render_Ts):
        """vjp from supplied output cotangents (the backward seed), bypassing
        the loss. Mirrors the removed per-kernel backward:

            v_render_colors = (v_rgb [C,H,W,3], v_depth [C,H,W,1], ...)
                              extra / None entries are ignored.
            v_render_Ts     = cotangent w.r.t. transmittance [C,H,W,1]
                              (for alpha = 1 - T, pass -v_alpha).

        Per-splat gradients are pulled back to host into self.v_splats_world
        as (means, quats, scales, opacities[N,1], features_dc, features_sh[N,K,3]).

        NOTE: v_depth is applied to the engine's raw rasterized depth
        (sum of per-fragment weighted depths), which is what engine_forward
        returns. gsplat's "ED" is that divided by alpha, so depth gradients are
        only directly comparable once that normalization is accounted for (or
        when the depth cotangent is zero). The 3dgs / mip path also does not
        produce a viewmats gradient.
        """
        assert self.primitive in ['3dgs', 'mip', '3dgut'] and not self.use_bvh, \
            "backward() is only implemented for the engine path"
        C = self.viewmats.shape[0]
        H, W = self.height, self.width

        v_rgb = v_render_colors[0].contiguous()
        if len(v_render_colors) > 1 and v_render_colors[1] is not None:
            v_depth = v_render_colors[1].contiguous()
        else:
            v_depth = torch.zeros(C, H, W, 1, dtype=torch.float32)
        v_Ts = v_render_Ts.contiguous()

        _C.engine_backward_from_render_grad(
            self._tv(v_rgb), self._tv(v_depth), self._tv(v_Ts)
        )
        self._read_back_grads()

    def _read_back_grads(self):
        """Copy engine().grad.* (D->H) into self.v_splats_world."""
        N = self.max_num_splats
        K = self.splats_world[5].shape[1] if len(self.splats_world) >= 6 else 0
        v_means       = torch.empty(N, 3, dtype=torch.float32)
        v_quats       = torch.empty(N, 4, dtype=torch.float32)
        v_scales      = torch.empty(N, 3, dtype=torch.float32)
        v_opacities   = torch.empty(N, 1, dtype=torch.float32)
        v_features_dc = torch.empty(N, 3, dtype=torch.float32)
        v_features_sh = torch.empty(N, K, 3, dtype=torch.float32)
        _C.engine_copy_grads_to_host(
            self._tv(v_means), self._tv(v_quats), self._tv(v_scales),
            self._tv(v_opacities), self._tv(v_features_dc), self._tv(v_features_sh),
        )
        self.v_splats_world = (
            v_means, v_quats, v_scales, v_opacities, v_features_dc, v_features_sh
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
            if getattr(self, "render_distortion", None) is not None:
                rgb_dist, depth_dist = self.render_distortion
                self.meta["rgb_distortion"] = rgb_dist
                self.meta["depth_distortion"] = depth_dist
            return
