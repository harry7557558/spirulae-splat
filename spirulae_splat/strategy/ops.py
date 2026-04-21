import numpy as np
from typing import Callable, Dict, List, Tuple, Union, Optional, Literal

import torch
import torch.nn.functional as F
from torch import Tensor

from spirulae_splat.splat._torch_impl import quat_to_rotmat, quat_scale_to_covar_preci

from spirulae_splat.splat.cuda._wrapper_projection import _make_lazy_cuda_func

from spirulae_splat.viewer_legacy.utils import (
    quat_scale_to_triangle_verts,
    split_triangles,
    split_triangles_4,
    triangle_verts_to_quat_scale_mean
)


def get_param_attr(params, attr):
    param = params[attr]
    if hasattr(param, 'optim_info') and 'num_splats' in param.optim_info:
        param = param[:param.optim_info['num_splats']]
    return param


def get_param_grad(params, attr):
    param = params[attr]
    if hasattr(param, 'optim_info') and 'num_splats' in param.optim_info:
        return param.grad[:param.optim_info['num_splats']]
    return param.grad


def get_param_gradr(params, attr):
    param = params[attr]
    if not hasattr(param, 'optim_info') or 'gradr' not in param.optim_info:
        return None
    gradr = param.optim_info['gradr']
    if 'num_splats' in param.optim_info:
        return gradr[:param.optim_info['num_splats']]
    return gradr


import math

N_BINOMS = 51
BINOMS = torch.zeros((N_BINOMS, N_BINOMS))
for n in range(N_BINOMS):
    for k in range(n + 1):
        BINOMS[n, k] = math.comb(n, k)

@torch.no_grad()
def _multinomial_sample(weights: Tensor, n: int, replacement: bool = True) -> Tensor:
    """Sample from a distribution using torch.multinomial or numpy.random.choice.

    This function adaptively chooses between `torch.multinomial` and `numpy.random.choice`
    based on the number of elements in `weights`. If the number of elements exceeds
    the torch.multinomial limit (2^24), it falls back to using `numpy.random.choice`.

    Args:
        weights (Tensor): A 1D tensor of weights for each element.
        n (int): The number of samples to draw.
        replacement (bool): Whether to sample with replacement. Default is True.

    Returns:
        Tensor: A 1D tensor of sampled indices.
    """
    num_elements = weights.size(0)
    weights = torch.nan_to_num(weights, 0.0, 0.0, 0.0)

    if num_elements <= 2**24:
        # Use torch.multinomial for elements within the limit
        return torch.multinomial(weights, n, replacement=replacement)
    else:
        # Fallback to numpy.random.choice for larger element spaces
        weights = weights / weights.sum()
        weights_np = weights.detach().cpu().numpy()
        sampled_idxs_np = np.random.choice(
            num_elements, size=n, p=weights_np, replace=replacement
        )
        sampled_idxs = torch.from_numpy(sampled_idxs_np)

        # Return the sampled indices on the original device
        return sampled_idxs.to(weights.device)


@torch.no_grad()
def _update_param_with_optimizer(
    param_fn: Callable[[str, Tensor], Tensor],
    optimizer_fn: Callable[[str, Tensor], Tensor],
    params: Union[Dict[str, torch.nn.Parameter], torch.nn.ParameterDict],
    optimizers: Dict[str, torch.optim.Optimizer],
    names: Union[List[str], None] = None,
):
    """Update the parameters and the state in the optimizers with defined functions.

    Args:
        param_fn: A function that takes the name of the parameter and the parameter itself,
            and returns the new parameter.
        optimizer_fn: A function that takes the key of the optimizer state and the state value,
            and returns the new state value.
        params: A dictionary of parameters.
        optimizers: A dictionary of optimizers, each corresponding to a parameter.
        names: A list of key names to update. If None, update all. Default: None.
    """
    if names is None:
        # If names is not provided, update all parameters
        names = list(params.keys())

    for name in names:
        param = params[name]
        new_param = param_fn(name, param)
        params[name] = new_param
        if name not in optimizers:
            if name.startswith('_'):
                continue
            assert not param.requires_grad, (
                f"Optimizer for {name} is not found, but the parameter is trainable."
                f"Got requires_grad={param.requires_grad}"
            )
            continue
        optimizer = optimizers[name]
        for i in range(len(optimizer.param_groups)):
            param_state = optimizer.state[param]
            del optimizer.state[param]
            for key in param_state.keys():
                if key not in ["step", "step1", "step2"]:
                    v = param_state[key]
                    if hasattr(param, 'optim_info'):
                        v.optim_info = param.optim_info
                    param_state[key] = optimizer_fn(key, v)
            optimizer.param_groups[i]["params"] = [new_param]
            optimizer.state[new_param] = param_state


@torch.no_grad()
def duplicate(
    params: Union[Dict[str, torch.nn.Parameter], torch.nn.ParameterDict],
    optimizers: Dict[str, torch.optim.Optimizer],
    state: Dict[str, Tensor],
    mask: Tensor,
):
    """Inplace duplicate the Gaussian with the given mask.

    Args:
        params: A dictionary of parameters.
        optimizers: A dictionary of optimizers, each corresponding to a parameter.
        mask: A boolean mask to duplicate the Gaussians.
    """
    device = mask.device
    sel = torch.where(mask)[0]

    def param_fn(name: str, p: Tensor) -> Tensor:
        p_new = torch.nn.Parameter(torch.cat([p, p[sel]]))
        if hasattr(p, 'optim_info'):
            p_new.optim_info = p.optim_info
        return p_new

    def optimizer_fn(key: str, v: Tensor) -> Tensor:
        v_new = torch.cat([v, torch.zeros((len(sel), *v.shape[1:]), device=device)])
        if hasattr(v, 'optim_info'):
            v_new.optim_info = v.optim_info
        return v_new

    # update the parameters and the state in the optimizers
    _update_param_with_optimizer(param_fn, optimizer_fn, params, optimizers)
    # update the extra running state
    for k, v in state.items():
        if isinstance(v, torch.Tensor):
            state[k] = torch.cat((v, v[sel]))


@torch.no_grad()
def split(
    params: Union[Dict[str, torch.nn.Parameter], torch.nn.ParameterDict],
    optimizers: Dict[str, torch.optim.Optimizer],
    state: Dict[str, Tensor],
    mask: Tensor,
    std_scale: float = 1.0,
    revised_opacity: bool = False,
):
    """Inplace split the Gaussian with the given mask.

    Args:
        params: A dictionary of parameters.
        optimizers: A dictionary of optimizers, each corresponding to a parameter.
        mask: A boolean mask to split the Gaussians.
        revised_opacity: Whether to use revised opacity formulation
          from arXiv:2404.06109. Default: False.
    """
    device = mask.device
    sel = torch.where(mask)[0]
    rest = torch.where(~mask)[0]

    scales = torch.exp(get_param_attr(params, "scales")[sel])
    scales_3d = scales
    if scales.shape[-1] == 2:
        sz = 0.1*torch.fmin(scales[...,0:1], scales[...,1:2])
        scales_3d = torch.cat((scales, sz), dim=-1) / 3.0
    quats = F.normalize(get_param_attr(params, "quats")[sel], dim=-1)
    rotmats = quat_to_rotmat(quats)  # [N, 3, 3]
    samples = std_scale * torch.einsum(
        "nij,nj,bnj->bni",
        rotmats,
        scales_3d,
        torch.randn(2, len(scales_3d), 3, device=device),
    )  # [2, N, 3]

    def param_fn(name: str, p: Tensor) -> Tensor:
        repeats = [2] + [1] * (p.dim() - 1)
        if name == "means":
            p_split = (p[sel] + samples).reshape(-1, 3)  # [2N, 3]
        elif name == "scales":
            p_split = torch.log(scales / 1.6).repeat(2, 1)  # [2N, 3]
        elif name == "opacities" and revised_opacity:
            new_opacities = 1.0 - torch.sqrt(1.0 - torch.sigmoid(p[sel]))
            p_split = torch.logit(new_opacities).repeat(repeats)  # [2N]
        else:
            p_split = p[sel].repeat(repeats)
        p_new = torch.nn.Parameter(torch.cat([p[rest], p_split]))
        if hasattr(p, 'optim_info'):
            p_new.optim_info = p.optim_info
        return p_new

    def optimizer_fn(key: str, v: Tensor) -> Tensor:
        v_split = torch.zeros((2 * len(sel), *v.shape[1:]), device=device)
        v_new = torch.cat([v[rest], v_split])
        if hasattr(v, 'optim_info'):
            v_new.optim_info = v.optim_info
        return v_new

    # update the parameters and the state in the optimizers
    _update_param_with_optimizer(param_fn, optimizer_fn, params, optimizers)
    # update the extra running state
    for k, v in state.items():
        if isinstance(v, torch.Tensor):
            repeats = [2] + [1] * (v.dim() - 1)
            v_new = v[sel].repeat(repeats)
            state[k] = torch.cat((v[rest], v_new))


@torch.no_grad()
def remove(
    params: Union[Dict[str, torch.nn.Parameter], torch.nn.ParameterDict],
    optimizers: Dict[str, torch.optim.Optimizer],
    state: Dict[str, Tensor],
    mask: Tensor,
):
    """Inplace remove the Gaussian with the given mask.

    Args:
        params: A dictionary of parameters.
        optimizers: A dictionary of optimizers, each corresponding to a parameter.
        mask: A boolean mask to remove the Gaussians.
    """
    sel = torch.where(~mask)[0]

    def param_fn(name: str, p: Tensor) -> Tensor:
        p_new = torch.nn.Parameter(p[sel], requires_grad=p.requires_grad)
        if hasattr(p, 'optim_info'):
            p_new.optim_info = p.optim_info
        return p_new

    def optimizer_fn(key: str, v: Tensor) -> Tensor:
        v_new = v[sel]
        if hasattr(v, 'optim_info'):
            v_new.optim_info = v.optim_info
        return v_new

    # update the parameters and the state in the optimizers
    _update_param_with_optimizer(param_fn, optimizer_fn, params, optimizers)
    # update the extra running state
    for k, v in state.items():
        if isinstance(v, torch.Tensor):
            state[k] = v[sel]


@torch.no_grad()
def reset_opa(
    params: Union[Dict[str, torch.nn.Parameter], torch.nn.ParameterDict],
    optimizers: Dict[str, torch.optim.Optimizer],
    state: Dict[str, Tensor],
    value: float,
):
    """Inplace reset the opacities to the given post-sigmoid value.

    Args:
        params: A dictionary of parameters.
        optimizers: A dictionary of optimizers, each corresponding to a parameter.
        value: The value to reset the opacities
    """

    def param_fn(name: str, p: Tensor) -> Tensor:
        if name == "opacities":
            p_new = torch.clamp(p, max=torch.logit(torch.tensor(value)).item())
            p_new = torch.nn.Parameter(p_new, requires_grad=p.requires_grad)
            if hasattr(p, 'optim_info'):
                p_new.optim_info = p.optim_info
            return p_new
        else:
            raise ValueError(f"Unexpected parameter name: {name}")

    def optimizer_fn(key: str, v: Tensor) -> Tensor:
        v_new = torch.zeros_like(v)
        if hasattr(v, 'optim_info'):
            v_new.optim_info = v.optim_info
        return v_new

    # update the parameters and the state in the optimizers
    _update_param_with_optimizer(
        param_fn, optimizer_fn, params, optimizers, names=["opacities"]
    )


@torch.no_grad()
def relocate(
    params: Union[Dict[str, torch.nn.Parameter], torch.nn.ParameterDict],
    optimizers: Dict[str, torch.optim.Optimizer],
    state: Dict[str, Tensor],
    mask: Tensor,
    is_3dgs: bool,
    probs: Optional[Tensor]=None,
    min_opacity: float = 0.005,
):
    """Inplace relocate some dead Gaussians to the lives ones.

    Args:
        params: A dictionary of parameters.
        optimizers: A dictionary of optimizers, each corresponding to a parameter.
        mask: A boolean mask to indicates which Gaussians are dead.
    """
    opacities = torch.sigmoid(get_param_attr(params, "opacities"))

    dead_indices = mask.nonzero(as_tuple=True)[0]
    alive_indices = (~mask).nonzero(as_tuple=True)[0]
    n = len(dead_indices)
    if n == 0:
        return

    # Sample for new GSs
    eps = torch.finfo(torch.float32).eps
    if probs is None:
        probs = opacities
    probs = probs[alive_indices].flatten()  # ensure its shape is [N,]
    sampled_idxs = _multinomial_sample(probs, n, replacement=True)
    sampled_idxs = alive_indices[sampled_idxs]
    new_opacities, new_scales = _make_lazy_cuda_func("compute_relocation")(
        opacities[sampled_idxs],
        torch.exp(get_param_attr(params, "scales"))[sampled_idxs],
        torch.bincount(sampled_idxs).to(torch.int32)[sampled_idxs] + 1,
        BINOMS.to(opacities.device),
        len(BINOMS)
    )
    new_opacities = torch.clamp(new_opacities, max=1.0 - eps, min=min_opacity)

    def param_fn(name: str, p: Tensor) -> Tensor:
        if name == "opacities":
            p[sampled_idxs] = torch.logit(new_opacities)
        elif name == "scales":
            p[sampled_idxs] = torch.log(new_scales)
        p[dead_indices] = p[sampled_idxs]
        return p

    def optimizer_fn(key: str, v: Tensor) -> Tensor:
        v[sampled_idxs] = 0
        return v

    # update the parameters and the state in the optimizers
    _update_param_with_optimizer(param_fn, optimizer_fn, params, optimizers)
    # update the extra running state
    for k, v in state.items():
        if isinstance(v, torch.Tensor):
            v[sampled_idxs] = 0


@torch.no_grad()
def relocate_long_axis_split(
    params: Union[Dict[str, torch.nn.Parameter], torch.nn.ParameterDict],
    optimizers: Dict[str, torch.optim.Optimizer],
    state: Dict[str, Tensor],
    mask: Tensor,
    probs: Optional[Tensor]=None,
):
    """Inplace relocate some dead Gaussians to the lives ones.

    Args:
        params: A dictionary of parameters.
        optimizers: A dictionary of optimizers, each corresponding to a parameter.
        mask: A boolean mask to indicates which Gaussians are dead.
    """
    opacities = torch.sigmoid(get_param_attr(params, "opacities"))

    dead_indices = mask.nonzero(as_tuple=True)[0]
    alive_indices = (~mask).nonzero(as_tuple=True)[0]
    n = len(dead_indices)
    if n == 0:
        return
    if n >= len(alive_indices):
        return  # training likely has become unstable

    # Sample for new GSs
    if probs is None:
        probs = opacities
    probs = probs[alive_indices].flatten()  # ensure its shape is [N,]
    # sampled_idxs = _multinomial_sample(probs, n, replacement=False)
    sampled_idxs = torch.argsort(probs)[-n:]
    sampled_idxs = alive_indices[sampled_idxs]
    new_scales, new_opacities, mean_offsets = _make_lazy_cuda_func("long_axis_split")(
        "3dgs",
        get_param_attr(params, "scales")[sampled_idxs],
        get_param_attr(params, "opacities")[sampled_idxs],
        get_param_attr(params, "quats")[sampled_idxs]
    )

    def param_fn(name: str, p: Tensor) -> Tensor:
        if name == "means":
            p[dead_indices] = p[sampled_idxs] + mean_offsets
            p[sampled_idxs] -= mean_offsets
        elif name == "scales":
            p[sampled_idxs] = new_scales
            p[dead_indices] = new_scales
        elif name == "opacities":
            p[sampled_idxs] = new_opacities
            p[dead_indices] = new_opacities
        else:
            p[dead_indices] = p[sampled_idxs]
        return p

    def optimizer_fn(key: str, v: Tensor) -> Tensor:
        v[sampled_idxs] = 0
        return v

    # update the parameters and the state in the optimizers
    _update_param_with_optimizer(param_fn, optimizer_fn, params, optimizers)
    # update the extra running state
    for k, v in state.items():
        if isinstance(v, torch.Tensor):
            v[sampled_idxs] = 0


@torch.no_grad()
def relocate_opaque_triangles(
    params: Union[Dict[str, torch.nn.Parameter], torch.nn.ParameterDict],
    optimizers: Dict[str, torch.optim.Optimizer],
    state: Dict[str, Tensor],
    mask: Tensor,
    probs: Tensor
):
    """Inplace relocate some dead triangles to the lives ones.

    Args:
        params: A dictionary of parameters.
        optimizers: A dictionary of optimizers, each corresponding to a parameter.
        mask: A boolean mask to indicates which triangles are dead.
    """

    dead_indices = mask.nonzero(as_tuple=True)[0]
    alive_indices = (~mask).nonzero(as_tuple=True)[0]
    n = len(dead_indices)
    ns = (n + 2) // 3

    # Sample for new triangles
    probs = probs[alive_indices].flatten()  # ensure its shape is [N,]
    sampled_idxs = _multinomial_sample(probs, ns, replacement=True)
    sampled_idxs = alive_indices[sampled_idxs]

    triangles = quat_scale_to_triangle_verts(
        get_param_attr(params, "quats")[sampled_idxs],
        get_param_attr(params, "scales")[sampled_idxs],
        get_param_attr(params, "means")[sampled_idxs]
    )
    # TODO: better way to handle triangles sampled >2 times?
    # split_1, split_2 = split_triangles(triangles)  # TODO: doesn't look nice, possibly bug?
    split_2, split_1 = split_triangles_4(triangles)
    quats1, scales1, means1 = triangle_verts_to_quat_scale_mean(split_1[:n])
    quats2, scales2, means2 = triangle_verts_to_quat_scale_mean(split_2)

    def param_fn(name: str, p: Tensor) -> Tensor:
        if name == "quats":
            p[dead_indices] = quats1
            p[sampled_idxs] = quats2
        elif name == "scales":
            p[dead_indices] = scales1
            p[sampled_idxs] = scales2
        elif name == "means":
            p[dead_indices] = means1
            p[sampled_idxs] = means2
        else:
            p[dead_indices] = p[sampled_idxs.repeat(3)[:n]]
        return p

    def optimizer_fn(key: str, v: Tensor) -> Tensor:
        v[dead_indices] = 0
        v[sampled_idxs] = 0
        return v

    # update the parameters and the state in the optimizers
    _update_param_with_optimizer(param_fn, optimizer_fn, params, optimizers)
    # update the extra running state
    for k, v in state.items():
        if isinstance(v, torch.Tensor):
            v[dead_indices] = 0
            v[sampled_idxs] = 0


@torch.no_grad()
def sample_add(
    params: Union[Dict[str, torch.nn.Parameter], torch.nn.ParameterDict],
    optimizers: Dict[str, torch.optim.Optimizer],
    state: Dict[str, Tensor],
    n: int,
    probs: Optional[Tensor]=None,
    min_opacity: float = 0.005
):
    # TODO: get this right for opaque triangle splatting
    opacities = torch.sigmoid(get_param_attr(params, "opacities"))

    eps = torch.finfo(torch.float32).eps
    if probs is None:
        probs = opacities
    probs = probs.flatten()

    sampled_idxs = _multinomial_sample(probs, n, replacement=True)
    new_opacities, new_scales = _make_lazy_cuda_func("compute_relocation")(
        opacities[sampled_idxs],
        torch.exp(get_param_attr(params, "scales"))[sampled_idxs],
        torch.bincount(sampled_idxs).to(torch.int32)[sampled_idxs] + 1,
        BINOMS.to(opacities.device),
        len(BINOMS)
    )
    new_opacities = torch.clamp(new_opacities, max=1.0 - eps, min=min_opacity)

    def param_fn(name: str, p: Tensor) -> Tensor:
        if name == "opacities":
            p[sampled_idxs] = torch.logit(new_opacities)
        elif name == "scales":
            p[sampled_idxs] = torch.log(new_scales)
        if hasattr(p, 'optim_info') and 'num_splats' in p.optim_info:
            num_splats = p.optim_info['num_splats']
            assert num_splats + len(sampled_idxs) <= len(p)
            p[num_splats:num_splats+len(sampled_idxs)] = p[sampled_idxs]
            return p
        p_new = torch.nn.Parameter(torch.cat([p, p[sampled_idxs]]), requires_grad=p.requires_grad)
        if hasattr(p, 'optim_info'):
            p_new.optim_info = p.optim_info
        return p_new

    def optimizer_fn(key: str, v: Tensor) -> Tensor:
        if hasattr(v, 'optim_info') and 'num_splats' in v.optim_info:
            num_splats = v.optim_info['num_splats']
            assert num_splats <= len(v)
            v[num_splats:].zero_()
            return v
        v_new = torch.cat([v, torch.zeros((len(sampled_idxs), *v.shape[1:]), device=v.device)])
        if hasattr(v, 'optim_info'):
            v_new.optim_info = v.optim_info
        return v_new

    # update the parameters and the state in the optimizers
    _update_param_with_optimizer(param_fn, optimizer_fn, params, optimizers)
    # update the extra running state
    for k, v in state.items():
        v_new = torch.zeros((len(sampled_idxs), *v.shape[1:]), device=v.device)
        if isinstance(v, torch.Tensor):
            state[k] = torch.cat((v, v_new))


@torch.no_grad()
def sample_add_long_axis_split(
    params: Union[Dict[str, torch.nn.Parameter], torch.nn.ParameterDict],
    optimizers: Dict[str, torch.optim.Optimizer],
    state: Dict[str, Tensor],
    n: int,
    probs: Optional[Tensor]=None
):
    # TODO: get this right for opaque triangle splatting
    opacities = torch.sigmoid(get_param_attr(params, "opacities"))

    if probs is None:
        probs = opacities
    probs = probs.flatten()

    if False:
        split(
            params=params,
            optimizers=optimizers,
            state=state,
            mask=(probs > torch.sort(probs.flatten(), descending=True)[0][n]),
            std_scale=1.0,
            revised_opacity=False,
        )
        return

    # sampled_idxs = _multinomial_sample(probs, n, replacement=False)
    sampled_idxs = torch.argsort(probs)[-n:]
    new_scales, new_opacities, mean_offsets = _make_lazy_cuda_func("long_axis_split")(
        "3dgs",
        get_param_attr(params, "scales")[sampled_idxs],
        get_param_attr(params, "opacities")[sampled_idxs],
        get_param_attr(params, "quats")[sampled_idxs]
    )

    def param_fn(name: str, p: Tensor) -> Tensor:
        if name == "means":
            p_cat = p[sampled_idxs] + mean_offsets
            p[sampled_idxs] -= mean_offsets
        elif name == "scales":
            p[sampled_idxs] = new_scales
            p_cat = new_scales
        elif name == "opacities":
            p[sampled_idxs] = new_opacities
            p_cat = new_opacities
        else:
            p_cat = p[sampled_idxs]
        if hasattr(p, 'optim_info') and 'num_splats' in p.optim_info:
            num_splats = p.optim_info['num_splats']
            assert num_splats + len(p_cat) <= len(p)
            p[num_splats:num_splats+len(p_cat)] = p_cat
            return p
        p_new = torch.nn.Parameter(torch.cat([p, p_cat]), requires_grad=p.requires_grad)
        if hasattr(p, 'optim_info'):
            p_new.optim_info = p.optim_info
        return p_new

    def optimizer_fn(key: str, v: Tensor) -> Tensor:
        if hasattr(v, 'optim_info') and 'num_splats' in v.optim_info:
            num_splats = v.optim_info['num_splats']
            assert num_splats <= len(v)
            v[num_splats:].zero_()
            return v
        v_new = torch.cat([v, torch.zeros((len(sampled_idxs), *v.shape[1:]), device=v.device)])
        if hasattr(v, 'optim_info'):
            v_new.optim_info = v.optim_info
        return v_new

    # update the parameters and the state in the optimizers
    _update_param_with_optimizer(param_fn, optimizer_fn, params, optimizers)
    # update the extra running state
    for k, v in state.items():
        v_new = torch.zeros((len(sampled_idxs), *v.shape[1:]), device=v.device)
        if isinstance(v, torch.Tensor):
            state[k] = torch.cat((v, v_new))


@torch.no_grad()
def sample_add_opaque_triangles(
    params: Union[Dict[str, torch.nn.Parameter], torch.nn.ParameterDict],
    optimizers: Dict[str, torch.optim.Optimizer],
    state: Dict[str, Tensor],
    n: int,
    probs: Tensor
):
    eps = torch.finfo(torch.float32).eps
    probs = probs.flatten()
    sampled_idxs = _multinomial_sample(probs, n, replacement=True)

    triangles = quat_scale_to_triangle_verts(
        get_param_attr(params, "quats")[sampled_idxs],
        get_param_attr(params, "scales")[sampled_idxs],
        get_param_attr(params, "means")[sampled_idxs]
    )
    # TODO: better way to handle triangles sampled >2 times?
    # split_1, split_2 = split_triangles(triangles)  # TODO: doesn't look nice, possibly bug?
    split_1, split_2 = split_triangles_4(triangles)
    quats1, scales1, means1 = triangle_verts_to_quat_scale_mean(split_1)
    quats2, scales2, means2 = triangle_verts_to_quat_scale_mean(split_2)

    def param_fn(name: str, p: Tensor) -> Tensor:
        if name == "quats":
            p[sampled_idxs] = quats1
            p_cat = quats2
        elif name == "scales":
            p[sampled_idxs] = scales1
            p_cat = scales2
        elif name == "means":
            p[sampled_idxs] = means1
            p_cat = means2
        else:
            p_cat = p[sampled_idxs].repeat(3, *([1]*(len(p.shape)-1)))
        if hasattr(p, 'optim_info') and 'num_splats' in p.optim_info:
            num_splats = p.optim_info['num_splats']
            assert num_splats + len(p_cat) <= len(p)
            p[num_splats:num_splats+len(p_cat)] = p_cat
            return p
        p_new = torch.nn.Parameter(torch.cat([p, p_cat]), requires_grad=p.requires_grad)
        if hasattr(p, 'optim_info'):
            p_new.optim_info = p.optim_info
        return p_new

    def optimizer_fn(key: str, v: Tensor) -> Tensor:
        if hasattr(v, 'optim_info') and 'num_splats' in v.optim_info:
            num_splats = v.optim_info['num_splats']
            assert num_splats <= len(v)
            v[num_splats:].zero_()
            return v
        v[sampled_idxs] = 0
        v_new = torch.cat([v, torch.zeros((3*len(sampled_idxs), *v.shape[1:]), device=v.device)])
        if hasattr(v, 'optim_info'):
            v_new.optim_info = v.optim_info
        return v_new

    # update the parameters and the state in the optimizers
    _update_param_with_optimizer(param_fn, optimizer_fn, params, optimizers)
    # update the extra running state
    for k, v in state.items():
        if isinstance(v, torch.Tensor):
            v[sampled_idxs] = 0
            v_new = torch.zeros((3*len(sampled_idxs), *v.shape[1:]), device=v.device)
            state[k] = torch.cat((v, v_new))


@torch.no_grad()
def inject_noise_to_position(
    primitive: Literal["3dgs", "mip", "3dgut", "opaque_triangle"],
    params: Union[Dict[str, torch.nn.Parameter], torch.nn.ParameterDict],
    optimizers: Dict[str, torch.optim.Optimizer],
    state: Dict[str, Tensor],
    scaler: float,
    min_opacity: float,
    opacities: Optional[Tensor] = None
):
    if scaler <= 0.0 or min_opacity <= 0.0:
        return

    if opacities is None:
        opacities = torch.sigmoid(get_param_attr(params, "opacities").flatten())

    _make_lazy_cuda_func("mcmc_add_noise")(
        primitive,
        scaler,
        min_opacity,
        get_param_attr(params, "means"),
        get_param_attr(params, "scales"),
        get_param_attr(params, "quats"),
        opacities
    )
