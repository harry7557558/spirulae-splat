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
        cur_num_splats: int
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
            assert len(means.shape) >= 2, "means must have at least 2 dimensions"
            batch_dims = means.shape[:-2]
            N = means.shape[-2]
            device = means.device
            assert means.shape == batch_dims + (N, 3), means.shape
            assert quats.shape == batch_dims + (N, 4), quats.shape
            assert scales.shape == batch_dims + (N, 3), scales.shape
            if primitive in ["3dgs", "mip", "3dgut", "3dgut_sv"]:
                assert opacities.shape == batch_dims + (N, 1), opacities.shape
            else:
                assert opacities.shape == batch_dims + (N, 2), opacities.shape
                assert features_ch.shape == batch_dims + (N, 2, 3), features_ch.shape
        elif primitive in ["voxel"]:
            assert len(splats_world) == 4, "Voxel requires 4 params (pos_sizes, densities, features_dc, features_sh)"
            pos_sizes, densities, features_dc, features_sh = splats_world
            batch_dims = pos_sizes.shape[:-2]
            N = pos_sizes.shape[-2]
            device = pos_sizes.device
            assert pos_sizes.shape == batch_dims + (N, 4), pos_sizes.shape
            assert densities.shape == batch_dims + (N, 8), densities.shape
        else:
            raise ValueError(f"Invalid primitive ({primitive})")
        if features_dc is not None and features_sh is not None:
            assert features_dc.shape == batch_dims + (N, 3), features_dc.shape
            assert features_sh.shape == batch_dims + (N, 0, 3) or \
                features_sh.shape == batch_dims + (N, 3, 3) or \
                features_sh.shape == batch_dims + (N, 8, 3) or \
                features_sh.shape == batch_dims + (N, 15, 3), features_sh.shape
        if sv_sites is not None and sv_colors is not None:
            num_sv = sv_sites.shape[-2]
            assert sv_sites.shape == batch_dims + (N, num_sv, 3)
            assert sv_colors.shape == batch_dims + (N, num_sv, 3)

        self.batch_dims = batch_dims
        self.device = device
        self.cur_num_splats = cur_num_splats
        self.max_num_splats = N

        self.primitive = primitive
        self.splats_world = splats_world

    def set_params(
        self,
        viewmats: Tensor,  # [..., C, 4, 4]
        intrins: Tensor,  # [..., C, 4]
        width: int,
        height: int,
        packed: bool = True,
        use_bvh: bool = False,
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
        self.viewmats = viewmats
        self.intrins = intrins
        self.width = width
        self.height = height
        self.packed = packed
        self.use_bvh = use_bvh
        self.output_distortion = output_distortion
        self.compute_hessian_diagonal = compute_hessian_diagonal
        self.relative_scale = relative_scale
        self.backgrounds = backgrounds
        self.accum_weight_map = accum_weight_map
        self.max_blending_masks = max_blending_masks
        self.camera_model = camera_model
        self.dist_coeffs = dist_coeffs
        self.actual_width = actual_width
        self.actual_height = actual_height
        if self.dist_coeffs is not None:
            self.dist_coeffs = self.dist_coeffs.contiguous()

    def zero_grad(self):
        if hasattr(self, 'v_splats_world'):
            for tensor in self.v_splats_world:
                _make_lazy_cuda_func("set_zero")(tensor)
        else:
            self.v_splats_world = [
                torch.zeros_like(x) for x in self.splats_world
            ]

        if not hasattr(self, 'g1_splats_world'):
            self.g1_splats_world = [
                torch.zeros_like(x) for x in self.splats_world
            ]
        if not hasattr(self, 'g2_splats_world'):
            self.g2_splats_world = [
                torch.zeros_like(x) for x in self.splats_world
            ]

        if hasattr(self, 'radii'):
            _make_lazy_cuda_func("set_zero")(self.radii)

    def projection_forward(self):
        if self.primitive not in ["3dgs", "mip", "3dgut", "3dgut_sv"]:
            raise NotImplementedError()

        proj_returns = _make_lazy_cuda_func(
            f"projection_{self.primitive}_packed_forward" if self.packed else
            f"projection_{self.primitive}_forward"
        )(
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
                v_splats_world, self.v_splats_proj, accum_weight = _make_lazy_cuda_func(
                    f"rasterization_{self.primitive}_backward"
                )(
                    self.max_num_splats,
                    self.splats_world, self.splats_proj, self.gaussian_ids,
                    self.width, self.height, self.isect_offsets, self.flatten_ids, self.render_Ts, self.render_last_ids,
                    self.backward_info.get('accum_weight_map', None),
                    self.v_render_colors,
                    self.v_render_Ts,
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
                    None, False  # TODO
                )
                if self.compute_hessian_diagonal:
                    v_splats, v_viewmats, vr_splats, h_splats, accum_weight = cuda_return
                else:
                    v_splats_world, self.v_splats_proj, v_viewmats, accum_weight = cuda_return
            # TODO: fuse
            if not hasattr(self, 'v_splats_world'):
                self.v_splats_world = v_splats_world
            else:
                for i in range(len(v_splats_world)):
                    self.v_splats_world[i] += v_splats_world[i].reshape(self.v_splats_world[i].shape)
        if accum_weight is not None:
            self.backward_info['accum_weight'] = accum_weight

        # TODO
        # self.v_backgrounds = None
        # if self.needs_input_grad[5]:
        #     self.v_backgrounds = (torch.cat(self.v_render_colors, dim=-1) * \
        #                      self.v_render_Ts.float()).sum(dim=(-3, -2))

    def projection_backward(self):

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
                self.splats_world,
                self.viewmats, self.intrins, self.width, self.height, self.camera_model.upper(), self.dist_coeffs,
                self.camera_ids, self.gaussian_ids, self.aabb,
                self.v_splats_proj,
                self.v_splats_world, self.v_viewmats,
            )

    def forward(self):

        if self.cur_num_splats < self.max_num_splats:
            self.splats_world[3].data[self.cur_num_splats:] = -10.0  # TODO

        self.meta = {}

        B = math.prod(self.batch_dims)
        C = self.viewmats.shape[-3]
        I = B * C
        assert self.viewmats.shape == self.batch_dims + (C, 4, 4), self.viewmats.shape
        assert self.intrins.shape == self.batch_dims + (C, 4), self.intrins.shape

        self.backward_info = {}
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
            I, self.intrins, self.width, self.height,
            self.camera_ids if self.packed else None,
        )
        self.isect_offsets = isect_offsets.reshape(self.batch_dims + (C, tile_height, tile_width))
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
                "n_batches": B,
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

    def backward(
        self,
        v_render_colors: Tuple[torch.Tensor, torch.Tensor, torch.Tensor],
        v_render_Ts: torch.Tensor,
        v_rgb_distortion: Optional[torch.Tensor] = None,
        v_depth_distortion: Optional[torch.Tensor] = None,
        v_normal_distortion: Optional[torch.Tensor] = None,
    ):
        assert len(v_render_colors) == 3, "v_render_colors must contain RGB, depth, and normal"
        for tensor in v_render_colors:
            assert tensor is None or tensor.is_contiguous()
        assert v_render_Ts.is_contiguous()

        if len(self.render_colors) > 1:

            v_raw_depth, v_t_from_depth = _make_lazy_cuda_func("rendered_depth_to_expected_depth_backward")(
                self.raw_depth, self.render_Ts, v_render_colors[1])

            v_render_Ts = v_render_Ts + v_t_from_depth
            v_render_colors = (
                v_render_colors[0],
                v_raw_depth,
                *v_render_colors[2:]
            )

        self.v_render_colors = v_render_colors
        self.v_render_Ts = v_render_Ts

        self.rasterize_backward()
        self.projection_backward()

        

    def optim_step(
        self,
        model_config: 'spirulae_splat.modules.model.SpirulaeSplatModelConfig',
        optim_config: OptimizerConfig,
        step: int,
        max_steps: int
    ):

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
            model_config.noise_lr,
            model_config.min_opacity,
            model_config.max_gauss_ratio,
            model_config.scale_regularization_weight,
            model_config.opacity_reg,
            model_config.scale_reg,
            # 0.0, 0.0,  # TODO
            model_config.erank_reg,
            model_config.erank_reg_s3,
            model_config.quat_norm_reg,
            0.0 if step % model_config.refine_every != 0 else
                (1.0 - (step+1) / max_steps) * model_config.opacity_decay,
            1.0 if step % model_config.refine_every != 0 else
                1.0 - (1.0 - (step+1) / max_steps) * model_config.scale_decay,
            optim_config.use_scale_agnostic_mean,
            bias_correction_step
        )

        _make_lazy_cuda_func("fused_adam_with_steps")(
            self.splats_world[4],
            self.v_splats_world[4],
            self.g1_splats_world[4],
            self.g2_splats_world[4],
            optim_config.get_scheduled_lr("features_dc", step, max_steps),
            bias_correction_step,
            model_config.sh_reg, 0.5 / 0.28209479177387814
        )

        _make_lazy_cuda_func("fused_adam_with_steps")(
            self.splats_world[5],
            self.v_splats_world[5],
            self.g1_splats_world[5],
            self.g2_splats_world[5],
            optim_config.get_scheduled_lr("features_sh", step, max_steps),
            bias_correction_step,
            model_config.sh_reg, 0.0
        )

    def densify_step(
        self,
        step: int, 
        max_steps: int,
        model_config: 'spirulae_splat.modules.model.SpirulaeSplatModelConfig',
        optim_config: OptimizerConfig
    ):
        # clip large splats
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

        if step >= max_steps - model_config.refine_stop_num_iter:
            return

        # update accumulation weight
        if model_config.relocate_heuristic_weight >= 1.0:

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

            _make_lazy_cuda_func("densify_update_weight")(
                self.cur_num_splats,
                self.radii,
                densify_scale,
                densify_opac,
                densify_score,
                self.densify_accum_buffer,
                is_max_mode
            )

        if step % model_config.refine_every != 0:
            return
        if step <= model_config.refine_start_iter:
            return

        if model_config.relocate_heuristic_weight >= 1.0:

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
        
        else:

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
        
        if hasattr(self, 'densify_accum_buffer'):
            _make_lazy_cuda_func("set_zero")(self.densify_accum_buffer)
