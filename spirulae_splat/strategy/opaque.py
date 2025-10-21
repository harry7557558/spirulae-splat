import math
import numpy as np
from dataclasses import dataclass
from typing import Any, Dict, Union, Literal, List, Tuple, Dict

import torch
from torch import Tensor

from .base import Strategy
from .ops import inject_noise_to_position, relocate, sample_add, remove, split, _multinomial_sample


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
    warmup_steps_0: int = 6_000
    warmup_steps_1: int = 15_000
    refine_stop_iter: int = 30_000
    refine_every: int = 100
    grow_factor: float = 1.05
    min_opacity: float = 0.005
    prune_opacity: float = 0.2
    relocate_scale2d: float = float('inf')
    max_scale2d: float = float('inf')
    max_scale3d: float = float('inf')
    prob_grad_weight: float = 0.0
    kernel_radius: float = 1.0
    verbose: bool = False
    key_for_gradient: Literal["means2d", "gradient_2dgs"] = "means2d"

    def initialize_state(self) -> Dict[str, Any]:
        """Initialize and return the running state for this strategy."""
        if self.prob_grad_weight > 0.0:
            return { "grad2d": None, "count": None, "radii": None }
        return { "radii": None }

    def get_opacity_floor(self, step):
        a = (step - self.warmup_steps_0) / (self.refine_stop_iter - self.warmup_steps_0)
        return max(min(a, 1.0), 0.0)

    def get_hardness(self, step):
        a = (step - self.warmup_steps_0) / (self.refine_stop_iter - self.warmup_steps_0)
        return max(min(a, 1.0), 0.0)

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
            self.key_for_gradient,
            "depths",
        ] + ["gaussian_ids"] * packed:
            assert key in info, f"{key} is required but missing."

        grad_info = info[self.key_for_gradient]

        # initialize state on the first run
        n_gaussian = len(list(params.values())[0])

        if self.prob_grad_weight > 0.0:
            if state["grad2d"] is None:
                state["grad2d"] = torch.zeros(n_gaussian, device=grad_info.device)
            if state["count"] is None:
                state["count"] = torch.zeros(n_gaussian, device=grad_info.device)
        if state["radii"] is None:
            assert "radii" in info, "radii is required but missing."
            state["radii"] = torch.zeros(n_gaussian, device=grad_info.device).float()

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

        # large splats in screen space
        # clip scale while increase opacity to encourage being relocated to
        if np.isfinite(self.max_scale2d):
            assert not packed
            # TODO: optionally, actually do anisotropic scale in 3d
            oversize_factor = torch.clip(normalized_radii / self.max_scale2d, min=1.0)
            oversize_factor = torch.log(oversize_factor).unsqueeze(-1)
            scales, opacities = params['scales'].data, params['opacities'].data
            scales[gs_ids] -= oversize_factor
            opacities[gs_ids] += scales.shape[-1] * oversize_factor
            opacities[gs_ids] = torch.clip(opacities[gs_ids], max=5.0)  # sigmoid(5.0)=0.993

        # large splats in world space
        # clip scale, without increasing opacity (which causes problems with background removal)
        if np.isfinite(self.max_scale3d):
            params['scales'].data.clip_(max=math.log(self.max_scale3d))

        if not (self.prob_grad_weight > 0.0):
            return

        # normalize grads to [-1, 1] screen space
        if not hasattr(grad_info, 'grad') or grad_info.grad is None:
            print("Error: grad not found")
            return
        grads = grad_info.grad.clone()
        if grads.shape[-1] == 3:  # world space, the 3DGUT case
            sel = sel[..., :1] & (info["depths"] > 0.01).unsqueeze(-1)  # [C, N]
            gs_ids = torch.where(sel)[1]  # [nnz]
            grads = grads.norm(dim=-1, keepdim=True)
            grads = grads * info["depths"].unsqueeze(-1)  # to screen space
        else:
            grads[..., 0] *= info["width"] / 2.0 * info["n_cameras"]
            grads[..., 1] *= info["height"] / 2.0 * info["n_cameras"]
        if len(grads.shape) == 2:
            assert info["n_cameras"] == 1
            grads = grads.unsqueeze(0)

        # update the running state
        if not packed:
            grads = grads[sel]  # [nnz, 2]
        if len(grads.shape) == 2:
            grads = grads.norm(dim=-1)
        state["grad2d"].index_add_(0, gs_ids, torch.relu(grads))
        state["count"].index_add_(
            0, gs_ids, torch.ones_like(gs_ids, dtype=torch.float32)
        )

    def _get_probs(self, state, params, step):
        if not (self.prob_grad_weight > 0.0):
            return None
        opacs = self.map_opacities(step, params["opacities"])
        grads = state["grad2d"] / state["count"].clamp_min(1)
        grads = ((0.5+torch.sort(grads)[1]) / len(grads)).unsqueeze(-1)
        return opacs * (1.0-self.prob_grad_weight) + grads * self.prob_grad_weight

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

        if (
            step >= self.refine_every
            and step % self.refine_every == 0
        ):
            # prune splats
            opacity_floor = self.get_opacity_floor(step)
            if opacity_floor > 0.0 and opacity_floor < self.prune_opacity:
                self._prune_gs(params, optimizers, state, step)

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
            # TODO: if after warmup_steps_1, split instead of creating overlapping splats?
            n_new_gs = self._add_new_gs(params, optimizers, state, step)
            if self.verbose:
                print(
                    f"Step {step}: Added {n_new_gs} GSs. "
                    f"Now having {len(params['means'])} GSs."
                )

            if self.prob_grad_weight > 0.0:
                state["grad2d"] *= 0.0
                state["count"] *= 0.0
            torch.cuda.empty_cache()

        # add noise to GSs
        scalar = lr * self.noise_lr #* min(step/max(self.refine_start_iter,1), 1)
        inject_noise_to_position(
            params=params, optimizers=optimizers, state={}, scaler=scalar,
            min_opacity=self.min_opacity,
            opacities=self.map_opacities(step, params["opacities"].flatten())
        )

    @torch.no_grad()
    def _relocate_gs(
        self,
        params: Union[Dict[str, torch.nn.Parameter], torch.nn.ParameterDict],
        optimizers: Dict[str, torch.optim.Optimizer],
        state: Dict[str, Any],
        step: int
    ) -> int:
        opacities = self.map_opacities(step, params["opacities"].flatten())
        relocate_mask = ~(torch.isfinite(opacities))

        # relocate huge splats
        # TODO: logically it's better to split in this case?
        if step < self.refine_stop_iter and np.isfinite(self.relocate_scale2d):
            relocate_mask |= (state["radii"] > self.relocate_scale2d)
        state["radii"] *= 0

        # relocate low opacity, as in original MCMC
        if step > self.refine_start_iter and step < self.refine_stop_iter:
            relocate_mask |= (opacities <= self.min_opacity)

        n_gs = relocate_mask.sum().item()
        if n_gs > 0:
            relocate(
                params=params,
                optimizers=optimizers,
                state=state,
                mask=relocate_mask,
                is_3dgs=True,
                probs=self._get_probs(state, params, step),
                min_opacity=self.min_opacity,
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
        
        # TODO: split triangle into two
        if step > self.warmup_steps_1 and False:
            probs = self.map_opacities(step, params["opacities"].flatten())
            samples = _multinomial_sample(probs, n_gs, replacement=False)
            is_split = torch.zeros(len(probs), device=probs.device, dtype=torch.bool)
            is_split[samples] = True
            split(
                params=params,
                optimizers=optimizers,
                state=state,
                mask=is_split,
                std_scale=self.kernel_radius/3.0,
                revised_opacity=False,
            )
        else:
            sample_add(
                params=params,
                optimizers=optimizers,
                state=state,
                n=n_gs,
                probs=self._get_probs(state, params, step),
                min_opacity=self.min_opacity,
            )
        return n_gs

    @torch.no_grad()
    def _prune_gs(
        self,
        params: Union[Dict[str, torch.nn.Parameter], torch.nn.ParameterDict],
        optimizers: Dict[str, torch.optim.Optimizer],
        state: Dict[str, Any],
        step: int
    ) -> int:
        opacity_floor = self.get_opacity_floor(step)
        opacities = self.map_opacities(step, params["opacities"].flatten())

        if True:
            # TODO: make it work if opacity floor grows faster than prune opacity
            a = (step - self.warmup_steps_0) / (self.warmup_steps_1 - self.warmup_steps_0)
            is_prune = opacities < self.prune_opacity * max(min(a, 1.0), 0.0)
        else:
            a = min(opacity_floor + self.min_opacity, self.get_opacity_floor(self.warmup_steps_1))
            is_prune = opacities < a

        n_prune = is_prune.sum().item()
        if n_prune > 0:
            remove(params=params, optimizers=optimizers, state=state, mask=is_prune)

        return n_prune
