import math
from dataclasses import dataclass
from typing import Any, Dict, Union, Literal

import torch
from torch import Tensor

from .base import Strategy
from .ops import inject_noise_to_position, relocate, sample_add


@dataclass
class MCMCStrategy(Strategy):
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
    refine_stop_iter: int = 25_000
    refine_every: int = 100
    grow_factor: float = 1.05
    min_opacity: float = 0.005
    relocate_scale2d: float = float('inf')
    prob_grad_weight: float = 0.0
    is_3dgs: bool = False
    verbose: bool = False
    key_for_gradient: Literal["means2d", "gradient_2dgs"] = "means2d"

    def initialize_state(self) -> Dict[str, Any]:
        """Initialize and return the running state for this strategy."""
        return { "grad2d": None, "count": None, "radii": None }

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

        if state["grad2d"] is None:
            state["grad2d"] = torch.zeros(n_gaussian, device=grad_info.device)
        if state["count"] is None:
            state["count"] = torch.zeros(n_gaussian, device=grad_info.device)
        if state["radii"] is None:
            assert "radii" in info, "radii is required but missing."
            state["radii"] = torch.zeros(n_gaussian, device=grad_info.device)

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
        state["radii"][gs_ids] = torch.maximum(
            state["radii"][gs_ids],
            # normalize radii to [0, 1] screen space
            radii / float(max(info["width"], info["height"])),
        )

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

    def _get_probs(self, state, params):
        if not (self.prob_grad_weight > 0.0):
            return None
        opacs = torch.sigmoid(params["opacities"])
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

        # move to the correct device
        if (
            step < self.refine_stop_iter
            and step > self.refine_start_iter
            and step % self.refine_every == 0
        ):
            # teleport GSs
            n_relocated_gs = self._relocate_gs(params, optimizers, state)
            if self.verbose:
                print(f"Step {step}: Relocated {n_relocated_gs} GSs.")

            # add new GSs
            n_new_gs = self._add_new_gs(params, optimizers, state)
            if self.verbose:
                print(
                    f"Step {step}: Added {n_new_gs} GSs. "
                    f"Now having {len(params['means'])} GSs."
                )

            torch.cuda.empty_cache()

        # add noise to GSs
        scalar = lr * self.noise_lr #* min(step/max(self.refine_start_iter,1), 1)
        inject_noise_to_position(
            params=params, optimizers=optimizers, state={}, scaler=scalar,
            min_opacity = self.min_opacity
        )

    @torch.no_grad()
    def _relocate_gs(
        self,
        params: Union[Dict[str, torch.nn.Parameter], torch.nn.ParameterDict],
        optimizers: Dict[str, torch.optim.Optimizer],
        state: Dict[str, Any],
    ) -> int:
        opacities = torch.sigmoid(params["opacities"].flatten())
        relocate_mask = opacities <= self.min_opacity
        relocate_mask |= state["radii"] > self.relocate_scale2d
        n_gs = relocate_mask.sum().item()
        if n_gs > 0:
            relocate(
                params=params,
                optimizers=optimizers,
                state=state,
                mask=relocate_mask,
                is_3dgs=self.is_3dgs,
                probs=self._get_probs(state, params),
                min_opacity=self.min_opacity,
            )
        return n_gs

    @torch.no_grad()
    def _add_new_gs(
        self,
        params: Union[Dict[str, torch.nn.Parameter], torch.nn.ParameterDict],
        optimizers: Dict[str, torch.optim.Optimizer],
        state: Dict[str, Any],
    ) -> int:
        current_n_points = len(params["means"])
        n_target = min(self.cap_max, int(self.grow_factor * current_n_points))
        n_gs = max(0, n_target - current_n_points)
        if n_gs > 0:
            sample_add(
                params=params,
                optimizers=optimizers,
                state=state,
                n=n_gs,
                probs=self._get_probs(state, params),
                min_opacity=self.min_opacity,
            )
        return n_gs
