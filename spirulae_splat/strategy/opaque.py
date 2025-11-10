import math
import numpy as np
from dataclasses import dataclass
from typing import Any, Dict, Union, Literal, List, Tuple, Dict

import torch
from torch import Tensor

from .base import Strategy
from .ops import (
    inject_noise_to_position, relocate, sample_add, remove, split, _multinomial_sample,
    relocate_opaque_triangles, sample_add_opaque_triangles
)


@dataclass
class OpaqueStrategy(Strategy):
    """Strategy that follows the paper:

    `3D Gaussian Splatting as Markov Chain Monte Carlo <https://arxiv.org/abs/2404.09591>`_

    This strategy will:

    - Periodically teleport GSs with low opacity to a place that has high opacity.
    - Periodically introduce new GSs sampled based on the opacity distribution.
    - Periodically perturb the GSs locations.

    Args:
        cap_max (int): Maximum number of GSs. Default to 1_000_000.
        noise_lr (float): MCMC samping noise learning rate. Default to 5e5.
        refine_start_iter (int): Start refining GSs after this iteration. Default to 500.
        refine_stop_iter (int): Stop refining GSs after this iteration. Default to 25_000.
        refine_every (int): Refine GSs every this steps. Default to 100.
        min_opacity (float): GSs with opacity below this value will be pruned. Default to 0.005.
        verbose (bool): Whether to print verbose information. Default to False.

    Examples:

        >>> from gsplat import MCMCStrategy, rasterization
        >>> params: Dict[str, torch.nn.Parameter] | torch.nn.ParameterDict = ...
        >>> optimizers: Dict[str, torch.optim.Optimizer] = ...
        >>> strategy = MCMCStrategy()
        >>> strategy.check_sanity(params, optimizers)
        >>> strategy_state = strategy.initialize_state()
        >>> for step in range(1000):
        ...     render_image, render_alpha, info = rasterization(...)
        ...     loss = ...
        ...     loss.backward()
        ...     strategy.step_post_backward(params, optimizers, strategy_state, step, info, lr=1e-3)

    """

    cap_max: int = 1_000_000
    noise_lr: float = 5e5
    refine_start_iter: int = 500
    warmup_steps: int = 6_000
    refine_stop_iter: int = 30_000
    refine_every: int = 100
    geometry_optimizer_stop_iter: int = 29000
    grow_factor: float = 1.05
    min_opacity: float = 0.005
    final_opacity_floor: float = 0.9999
    final_hardness: float = 0.99
    hard_prune_opacity: float = 0.2
    gradual_prune_opacity: float = 0.45
    relocate_scale2d: float = float('inf')
    max_scale2d: float = float('inf')
    max_scale3d: float = float('inf')
    verbose: bool = False
    key_for_gradient: Literal["means2d", "gradient_2dgs"] = "means2d"

    def initialize_state(self) -> Dict[str, Any]:
        """Initialize and return the running state for this strategy."""
        return { "radii": None, "max_blending": None, "n_train_seen": 0 }

    def get_opacity_floor(self, step):
        a = (step - self.warmup_steps) / (self.refine_stop_iter - self.warmup_steps)
        return 1.0 - (1.0 - self.final_opacity_floor) ** max(min(a, 1.0), 0.0)

    def get_hardness(self, step):
        a = (step - self.warmup_steps) / (self.refine_stop_iter - self.warmup_steps)
        return 1.0 - (1.0 - self.final_hardness) ** max(min(a, 1.0), 0.0)

    def map_opacities(self, step, opacities):
        opacity_floor = self.get_opacity_floor(step)
        return opacity_floor + (1.0-opacity_floor) * torch.sigmoid(opacities)

    @torch.no_grad()
    def _update_state(
        self,
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
            "max_blending",
            self.key_for_gradient,
            "depths",
            "n_train",
        ] + ["gaussian_ids"] * packed:
            assert key in info, f"{key} is required but missing."

        grad_info = info[self.key_for_gradient]

        # initialize state on the first run
        n_gaussian = len(list(params.values())[0])

        if state["radii"] is None:
            assert "radii" in info, "radii is required but missing."
            state["radii"] = torch.zeros(n_gaussian, device=grad_info.device).float()

        if state["max_blending"] is None:
            assert "max_blending" in info, "max_blending is required but missing."
            state["max_blending"] = torch.zeros(n_gaussian, device=grad_info.device).float()

        state["n_train_seen"] += info["n_cameras"]

        # update the running state
        if packed:
            # grads is [nnz, 2]
            gs_ids = info["gaussian_ids"]  # [nnz]
            radii = info["radii"]  # [nnz]
        else:
            # grads is [C, N, 2]
            sel = info["radii"] > 0.0  # [C, N]
            gs_ids = torch.where(sel)[1]  # [nnz]
            radii = info["radii"][sel]  # [nnz]

        # Should be ideally using scatter max
        normalized_radii = radii / float(max(info["width"], info["height"]))
        state["radii"][gs_ids] = torch.maximum(
            state["radii"][gs_ids],
            normalized_radii
        )
        state["max_blending"] = torch.maximum(
            state["max_blending"],
            info["max_blending"]
        )

        # large splats in screen space
        # clip scale while increase opacity to encourage being relocated to
        if np.isfinite(self.max_scale2d):
            assert not packed
            # TODO: optionally, actually do anisotropic scale in 3d
            oversize_factor = torch.clip(normalized_radii / self.max_scale2d, min=1.0, max=1.1)
            oversize_factor = torch.log(oversize_factor).unsqueeze(-1)
            scales, opacities = params['scales'].data, params['opacities'].data
            scales[gs_ids] -= oversize_factor
            opacities[gs_ids] += scales.shape[-1] * oversize_factor
            opacities[gs_ids] = torch.clip(opacities[gs_ids], max=5.0)  # sigmoid(5.0)=0.993
            state["radii"][gs_ids] *= (normalized_radii <= self.max_scale2d).float()

        # large splats in world space
        # clip scale, without increasing opacity (which causes problems with background removal)
        if np.isfinite(self.max_scale3d):
            params['scales'].data.clip_(max=math.log(self.max_scale3d))

    def _get_probs(self, state, params, step):
        # opacs = self.map_opacities(step, params["opacities"].flatten())
        opacs = torch.sigmoid(params["opacities"].flatten())
        return opacs
        # scales = torch.exp(params["scales"][..., :2].mean(-1))
        scales = torch.nan_to_num(state["radii"], 0.0, 0.0, 0.0).clip(min=self.max_scale2d/2)
        return opacs * scales
        # return opacs * scales**0.5

    def check_sanity(
        self,
        params: Union[Dict[str, torch.nn.Parameter], torch.nn.ParameterDict],
        optimizers: Dict[str, torch.optim.Optimizer],
    ):
        """Sanity check for the parameters and optimizers.

        Check if:
            * `params` and `optimizers` have the same keys.
            * Each optimizer has exactly one param_group, corresponding to each parameter.
            * The following keys are present: {"means", "scales", "quats", "opacities"}.

        Raises:
            AssertionError: If any of the above conditions is not met.

        .. note::
            It is not required but highly recommended for the user to call this function
            after initializing the strategy to ensure the convention of the parameters
            and optimizers is as expected.
        """

        super().check_sanity(params, optimizers)
        # The following keys are required for this strategy.
        for key in ["means", "scales", "quats", "opacities"]:
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
        pass

    def step_post_backward(
        self,
        params: Union[Dict[str, torch.nn.Parameter], torch.nn.ParameterDict],
        optimizers: Dict[str, torch.optim.Optimizer],
        state: Dict[str, Any],
        step: int,
        info: Dict[str, Any],
        lr: float,
    ):
        """Callback function to be executed after the `loss.backward()` call.

        Args:
            lr (float): Learning rate for "means" attribute of the GS.
        """
        self._update_state(params, state, info)

        if step == self.warmup_steps:
            state["max_blending"] = torch.zeros_like(state["max_blending"])
        if step > self.warmup_steps and state["n_train_seen"] > info["n_train"]:
            state["n_train_seen"] = 0
            self._relocate_gs_max_blending(params, optimizers, state, step)
            state["max_blending"] = torch.zeros_like(state["max_blending"])

        if step > self.geometry_optimizer_stop_iter:
            for key in optimizers:
                if key in ['means', 'scales', 'quats']:
                    optimizers[key].param_groups[0]['lr'] = 0.0

        if (
            step >= self.refine_every
            and step % self.refine_every == 0
        ):
            # relocate splats
            n_relocated_gs = self._relocate_gs(params, optimizers, state, step)
            if self.verbose:
                print(f"Step {step}: Relocated {n_relocated_gs} GSs.")

        if (
            step < self.refine_stop_iter
            and step > self.refine_start_iter
            and step % self.refine_every == 0
        ):
            # add new splats
            n_new_gs = self._add_new_gs(params, optimizers, state, step)
            if self.verbose:
                print(
                    f"Step {step}: Added {n_new_gs} GSs. "
                    f"Now having {len(params['means'])} GSs."
                )

            torch.cuda.empty_cache()

        # add noise to GSs
        scalar = lr * self.noise_lr #* min(step/max(self.refine_start_iter,1), 1)
        sc = max(1 - self.get_opacity_floor(step) / self.gradual_prune_opacity, 0.0)
        inject_noise_to_position(
            params=params, optimizers=optimizers, state={}, scaler=scalar*sc,
            min_opacity=self.min_opacity*sc,
            # opacities=self.map_opacities(step, params["opacities"].flatten())
            opacities=self.map_opacities(0, params["opacities"].flatten())
        )

    @torch.no_grad()
    def _relocate_gs(
        self,
        params: Union[Dict[str, torch.nn.Parameter], torch.nn.ParameterDict],
        optimizers: Dict[str, torch.optim.Optimizer],
        state: Dict[str, Any],
        step: int
    ) -> int:
        opacity_floor = self.get_opacity_floor(step)
        opacities = self.map_opacities(step, params["opacities"].flatten())
        relocate_mask = ~(torch.isfinite(opacities))

        # relocate huge splats
        # TODO: logically it's better to split in this case?
        if np.isfinite(self.relocate_scale2d):
            relocate_mask |= (state["radii"] > self.relocate_scale2d)
        state["radii"] *= 0

        # relocate low opacity, as in original MCMC
        if step > self.refine_start_iter and step < self.refine_stop_iter:
            # min_opacity = opacity_floor + self.min_opacity * (1.0 - opacity_floor)
            min_opacity = self.min_opacity
            relocate_mask |= (opacities <= min_opacity)

        # hard prune
        if step // self.refine_every == self.warmup_steps // self.refine_every:
            relocate_mask |= (opacities < self.hard_prune_opacity)

        # get rid of floaters to prepare for opacity increase
        if opacity_floor > 0.0 and opacity_floor < self.gradual_prune_opacity:
            a = opacity_floor + self.min_opacity * max(1 - opacity_floor / self.gradual_prune_opacity, 0.0)
            relocate_mask |= (opacities < a)

        # relocate low opacity ones to high opacity ones
        if False:
            relocate_quantile = (self.grow_factor - 1) * (1 - opacity_floor)
            relocate_mask |= (opacities < torch.quantile(opacities, relocate_quantile))

        n_gs = relocate_mask.sum().item()
        if n_gs > 0:
            relocate_opaque_triangles(
                params=params,
                optimizers=optimizers,
                state=state,
                mask=relocate_mask,
                probs=self._get_probs(state, params, step),
            )
        return n_gs

    @torch.no_grad()
    def _add_new_gs(
        self,
        params: Union[Dict[str, torch.nn.Parameter], torch.nn.ParameterDict],
        optimizers: Dict[str, torch.optim.Optimizer],
        state: Dict[str, Any],
        step: int
    ) -> int:
        current_n_points = len(params["means"])
        n_target = min(self.cap_max, int(self.grow_factor * current_n_points))
        n_gs = max(0, n_target - current_n_points)
        if not (n_gs > 0):
            return 0
        
        sample_add_opaque_triangles(
            params=params,
            optimizers=optimizers,
            state=state,
            n=(n_gs+2)//3,
            probs=self._get_probs(state, params, step),
        )
        return n_gs

    @torch.no_grad()
    def _relocate_gs_max_blending(
        self,
        params: Union[Dict[str, torch.nn.Parameter], torch.nn.ParameterDict],
        optimizers: Dict[str, torch.optim.Optimizer],
        state: Dict[str, Any],
        step: int
    ) -> int:
        # prune_threshold = min(self.hard_prune_opacity, self.get_opacity_floor(step))
        prune_threshold = self.min_opacity
        relocate_mask = (state["max_blending"] < prune_threshold)

        n_gs = relocate_mask.sum().item()
        if n_gs > 0:
            relocate_opaque_triangles(
                params=params,
                optimizers=optimizers,
                state=state,
                mask=relocate_mask,
                probs=self._get_probs(state, params, step),
            )
        return n_gs
