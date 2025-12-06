import math
import numpy as np
from dataclasses import dataclass
from typing import Any, Dict, Union, Literal, List, Tuple, Dict

import torch
from torch import Tensor

from .base import Strategy
from .ops import (
    _update_param_with_optimizer,
    inject_noise_to_position, relocate, sample_add, remove, split, _multinomial_sample,
    relocate_opaque_triangles, sample_add_opaque_triangles
)

from spirulae_splat.splat.cuda import (
    svhash_get_voxels,
    svhash_split_voxels,
)



@dataclass
class SVRasterStrategy(Strategy):
    """Strategy that loosely based on the paper:

    `Sparse Voxels Rasterization: Real-time High-fidelity Radiance Field Rendering <https://arxiv.org/abs/2412.04459>`_

    """

    prune_opa: float = 0.005
    grow_grad2d: float = 0.0002
    prune_grad2d: float = 0.0
    grow_scale3d: float = 0.01
    grow_scale2d: float = 0.05
    prune_scale3d: float = 0.1
    prune_scale2d: float = 0.15
    split_scale3d: float = float('inf')
    refine_scale2d_stop_iter: int = 0
    refine_start_iter: int = 500
    refine_stop_iter: int = 25_000
    split_stop_iter: int = 15_000
    reset_every: int = 3000
    refine_every: int = 100
    pause_refine_after_reset: int = 0
    kernel_radius: float = 3.0
    absgrad: bool = False
    revised_opacity: bool = False
    verbose: bool = False

    def initialize_state(self, scene_scale: float = 1.0) -> Dict[str, Any]:
        state = {"density_grads": None, "count": None, "scene_scale": scene_scale}
        return state

    def check_sanity(
        self,
        params: Union[Dict[str, torch.nn.Parameter], torch.nn.ParameterDict],
        optimizers: Dict[str, torch.optim.Optimizer],
    ):
        super().check_sanity(params, optimizers)
        # The following keys are required for this strategy.
        for key in ["densities"]:
            assert key in params, f"{key} is required in params but missing."

    def step_pre_backward(
        self,
        params: Union[Dict[str, torch.nn.Parameter], torch.nn.ParameterDict],
        optimizers: Dict[str, torch.optim.Optimizer],
        state: Dict[str, Any],
        step: int,
        info: Dict[str, Any],
    ):
        """Callback function to be executed before the `loss.backward()` call."""
        assert (
            self.key_for_gradient in info
        ), "The 2D means of the Gaussians is required but missing."
        info[self.key_for_gradient].retain_grad()

    def step_post_backward(
        self,
        svhash,
        voxels,
        voxel_indices,
        params: Union[Dict[str, torch.nn.Parameter], torch.nn.ParameterDict],
        optimizers: Dict[str, torch.optim.Optimizer],
        state: Dict[str, Any],
        step: int,
        info: Dict[str, Any],
        packed: bool = False,
    ):
        """Callback function to be executed after the `loss.backward()` call."""

        self._update_state(voxel_indices, params, state, info, packed=packed)

        if (
            step > self.refine_start_iter and step < self.refine_stop_iter
            # and step % self.refine_every == 0
            and step % self.refine_every == max(self.refine_every // 10, 1)
            and step % self.reset_every >= self.pause_refine_after_reset
        ):
            # grow voxels
            if step < self.split_stop_iter:
                svhash = self._grow_voxels(svhash, params, optimizers, state, step)
                voxels, voxel_indices = svhash_get_voxels(svhash)

            # prune voxels
            # TODO

            # reset running stats
            state["density_grads"].zero_()
            state["count"].zero_()
            torch.cuda.empty_cache()

        return svhash, voxels, voxel_indices

    def _update_state(
        self,
        voxel_indices,
        params: Union[Dict[str, torch.nn.Parameter], torch.nn.ParameterDict],
        state: Dict[str, Any],
        info: Dict[str, Any],
        packed: bool = False,
    ):
        for key in [
            "width",
            "height",
            "n_cameras",
            "radii",
        ] + ["gaussian_ids"] * packed:
            assert key in info, f"{key} is required but missing."

        device = params['densities'].device

        # initialize state on the first run
        n_voxels = len(params['features_dc'])

        if state["density_grads"] is None:
            state["density_grads"] = torch.zeros(n_voxels, device=device)
        if state["count"] is None:
            state["count"] = torch.zeros(n_voxels, device=device)

        # update the running state
        if packed:
            voxel_ids = info["gaussian_ids"]  # [nnz]
        else:
            sel = torch.amax(info["radii"], dim=-1) > 0.0  # [C, N]
            voxel_ids = torch.where(sel)[1]  # [nnz]

        # TODO: compute absgrad density per paper
        voxel_densities = torch.exp(2.0*params['densities'])[voxel_indices]
        voxel_density_grads = params['densities'].grad[voxel_indices]
        grads = (voxel_densities * torch.abs(voxel_density_grads)).mean(-1)

        # update the running state
        if not packed:
            grads = grads[voxel_ids]  # [nnz]
        state["density_grads"].index_add_(0, voxel_ids, grads)
        state["count"].index_add_(
            0, voxel_ids, torch.ones_like(voxel_ids, dtype=torch.float32)
        )

    @torch.no_grad()
    def _grow_voxels(
        self,
        svhash,
        params: Union[Dict[str, torch.nn.Parameter], torch.nn.ParameterDict],
        optimizers: Dict[str, torch.optim.Optimizer],
        state: Dict[str, Any],
        step: int,
    ):
        count = state["count"]
        grads = state["density_grads"] / count.clamp_min(1)
        device = grads.device

        grad_threshold = max(torch.quantile(grads, 1-0.5/8), self.grow_grad2d)
        for i in range(20):
            print(grad_threshold)

        is_split = (grads > grad_threshold)

        svhash, cell_idx, vert_idx, vert_weight = svhash_split_voxels(svhash, is_split)

        def param_fn(name: str, p: Tensor) -> Tensor:
            # per vertex attributes
            if name in ["densities"]:
                return torch.cat((
                    p,
                    (p[vert_idx.flatten()].reshape(len(vert_idx), 8, *p.shape[1:]) * vert_weight).sum(1)
                ), 0)
            # per cell attributes
            else:
                return torch.cat((p, p[cell_idx]), 0)

        def optimizer_fn(key: str, v: Tensor) -> Tensor:
            return param_fn(key, v)

        _update_param_with_optimizer(param_fn, optimizer_fn, params, optimizers)

        for k, v in state.items():
            if isinstance(v, torch.Tensor):
                state[k] = torch.cat((v, v[cell_idx]), 0)

        return svhash


    @torch.no_grad()
    def _prune_voxels(
        self,
        params: Union[Dict[str, torch.nn.Parameter], torch.nn.ParameterDict],
        optimizers: Dict[str, torch.optim.Optimizer],
        state: Dict[str, Any],
        step: int,
    ) -> int:
        # TODO
        return
