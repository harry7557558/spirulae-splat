from dataclasses import dataclass, asdict
from typing import Any, Dict, List, Optional, Type, Literal

import torch
from torch.optim.optimizer import Optimizer

from spirulae_splat.splat.cuda import (
    _C,
    _make_lazy_cuda_func,
)


def get_scheduled_lr(param_group, step):
    lr = param_group['lr']
    lr_final = param_group['lr_final']
    lr_pre_warmup = param_group['lr_pre_warmup']
    warmup_steps = param_group['warmup_steps']
    max_steps = param_group['max_steps']
    scheduled_lr = lr
    if lr_final is not None:
        scheduled_lr = lr * (lr_final / lr) ** min(step / max_steps, 1.0)
    if warmup_steps is not None:
        scheduled_lr = min(scheduled_lr, lr_pre_warmup + (lr - lr_pre_warmup) * min(step / warmup_steps, 1.0))
    return scheduled_lr


class FusedAdam(Optimizer):
    """
    Fully fused CUDA implementation of Adam optimizer.
    """
    
    def __init__(
        self,
        params,
        lr: float = 1e-3,
        lr_final: Optional[float] = None,
        lr_pre_warmup: float = 0.0,
        betas: tuple = (0.9, 0.999),
        eps: float = 1e-8,
        tr: Optional[float] = 1e-6,
        tr_final: Optional[float] = 1e-8,
        warmup_steps: Optional[int] = None,
        max_steps: Optional[int] = 30000,
        **kwargs
    ):
        if lr < 0.0:
            raise ValueError(f"Invalid learning rate: {lr}")
        if eps < 0.0:
            raise ValueError(f"Invalid epsilon value: {eps}")
        if not 0.0 <= betas[0] < 1.0:
            raise ValueError(f"Invalid beta parameter at index 0: {betas[0]}")
        if not 0.0 <= betas[1] < 1.0:
            raise ValueError(f"Invalid beta parameter at index 1: {betas[1]}")

        # print(kwargs)
        if 'weight_decay' in kwargs and kwargs['weight_decay'] != 0.0:
            raise NotImplementedError("FusedAdam currently only supports weight_decay=0")
        
        defaults = dict(
            lr=lr, lr_final=lr_final, lr_pre_warmup=lr_pre_warmup, betas=betas, eps=eps,
            tr=tr, tr_final=tr_final, warmup_steps=warmup_steps, max_steps=max_steps
        )
        super(FusedAdam, self).__init__(params, defaults)

    @torch.no_grad()
    def step(self, closure=None):
        """Performs a single optimization step.
        
        Args:
            closure: A closure that reevaluates the model and returns the loss.
        """
        loss = None
        if closure is not None:
            with torch.enable_grad():
                loss = closure()
        
        for group in self.param_groups:
            beta1, beta2 = group['betas']
            lr = group['lr']
            lr_final = group['lr_final']
            lr_pre_warmup = group['lr_pre_warmup']
            eps = group['eps']

            tr = group['tr']
            tr_final = group['tr_final']
            warmup_steps = group['warmup_steps']
            max_steps = group['max_steps']
            
            # Collect all tensors for this group
            params_with_grad = []
            grads = []
            exp_avgs = []
            exp_avg_sqs = []
            steps = []
            
            for p in group['params']:
                if p.grad is None:
                    continue
                
                state = self.state[p]
                
                # State initialization
                if len(state) == 0:
                    state['step'] = 0
                    if hasattr(p, "optim_info") and p.optim_info.get("optimizer_offload", False):
                        state['exp_avg'] = torch.zeros_like(p, memory_format=torch.contiguous_format, device="cpu")
                        state['exp_avg_sq'] = torch.zeros_like(p, memory_format=torch.contiguous_format, device="cpu")
                    else:
                        state['exp_avg'] = torch.zeros_like(p, memory_format=torch.contiguous_format)
                        state['exp_avg_sq'] = torch.zeros_like(p, memory_format=torch.contiguous_format)
                p.exp_avg = state['exp_avg']
                p.exp_avg_sq = state['exp_avg_sq']
                
                state['step'] += 1
                scheduled_lr = get_scheduled_lr(group, state['step'])

                if not (hasattr(p, "optim_info") and "optimizer_override" in p.optim_info):
                    # Use default Adam
                    params_with_grad.append(
                        p[:p.optim_info['num_splats']] if hasattr(p, 'optim_info') and 'num_splats' in p.optim_info else p
                    )
                    grads.append(p.grad)
                    exp_avgs.append(state['exp_avg'])
                    exp_avg_sqs.append(state['exp_avg_sq'])
                    steps.append(state['step'])
                    continue

                additional_params = []
                eps_params = [eps]
                eps_tr = tr * (tr_final / tr) ** min(state['step'] / max_steps, 1.0)
                if p.optim_info['optimizer_override'] == "fused_adam_linear_rgb_optim":
                    pass
                elif p.optim_info['optimizer_override'] in ["fused_adamtr_linear_rgb_optim", "fused_adamtr_rgb_optim"]:
                    additional_params = [p.optim_info['opacities']]
                    eps_params.append(eps_tr)
                elif p.optim_info['optimizer_override'] in ["fused_adamtr_linear_rgb_sh_optim", "fused_adamtr_rgb_sh_optim"]:
                    additional_params = [p.optim_info['features_dc'], p.optim_info['opacities']]
                    eps_params.append(eps_tr)
                elif p.optim_info['optimizer_override'] in ["fused_adam_scale_agnostic_mean"]:
                    scales = p.optim_info['scales']
                    radii = p.optim_info['radii']
                    # if hasattr(p, 'optim_info') and 'num_splats' in p.optim_info:
                    #     scales = scales[:p.optim_info['num_splats']]
                    #     radii = radii[:p.optim_info['num_splats']]
                    # additional_params = [scales - scales.mean(), p.optim_info['quats'], p.optim_info['opacities'], radii / radii.mean()]
                    additional_params = [scales, p.optim_info['quats'], p.optim_info['opacities'], radii]
                    # print(torch.median((torch.exp(scales)/radii.unsqueeze(-1).repeat(1,3))[radii>0]))  # 0.2-0.3 for Mip-NeRF 360
                    eps_params.append(eps_tr)

                _make_lazy_cuda_func(p.optim_info['optimizer_override'])(
                    p[:p.optim_info['num_splats']] if hasattr(p, 'optim_info') and 'num_splats' in p.optim_info else p,
                    p.grad,
                    state['exp_avg'],
                    state['exp_avg_sq'],
                    *additional_params,
                    scheduled_lr,
                    beta1,
                    beta2,
                    *eps_params,
                    state['step']
                )

            if len(params_with_grad) == 0:
                continue
            
            # Use the step count from the first parameter (all should be the same)
            step = steps[0]
            scheduled_lr = get_scheduled_lr(group, step)
            
            # Launch fused CUDA kernel for all tensors
            _make_lazy_cuda_func("fused_adam_multi")(
                params_with_grad,
                grads,
                exp_avgs,
                exp_avg_sqs,
                scheduled_lr,
                beta1,
                beta2,
                eps,
                step
            )
        
        return loss


@dataclass
class FusedAdamOptimizerConfig:
    """Basic optimizer config with Adam"""

    _target: Type = FusedAdam

    lr: float = 1e-3
    lr_final: Optional[float] = None
    warmup_steps: Optional[int] = None
    lr_pre_warmup: float = 0.0

    betas: tuple = (0.9, 0.999)
    eps: float = 1e-8

    tr: Optional[float] = 1e-6
    tr_final: Optional[float] = 1e-8
    max_steps: Optional[int] = 30000



class Fused3DGS2Tr(Optimizer):
    """
    Fully fused CUDA implementation of 3DGS^2-TR optimizer for Gaussian means.
    https://arxiv.org/abs/2602.00395
    """
    
    def __init__(
        self,
        params,
        mode: Optional[Literal["mean", "scale", "opacity", "quat"]],
        lr: float = 1e-6,
        lr_final: Optional[float] = None,
        lr_pre_warmup: float = 0.0,
        betas: tuple = (0.9, 0.999),
        eps: float = 1e-8,
        warmup_steps: Optional[int] = None,
        max_steps: Optional[int] = 30000,
        # eps_tr: float = 1e-6,
        **kwargs
    ):
        if mode is None:
            raise ValueError("Mode must be set")
        if lr < 0.0:
            raise ValueError(f"Invalid learning rate: {lr}")
        if eps < 0.0:
            raise ValueError(f"Invalid epsilon value: {eps}")
        # if eps_tr < 0.0:
        #     raise ValueError(f"Invalid epsilon value: {eps_tr}")
        if not 0.0 <= betas[0] < 1.0:
            raise ValueError(f"Invalid beta parameter at index 0: {betas[0]}")
        if not 0.0 <= betas[1] < 1.0:
            raise ValueError(f"Invalid beta parameter at index 1: {betas[1]}")

        if 'weight_decay' in kwargs and kwargs['weight_decay'] != 0.0:
            raise NotImplementedError("Fused3DGS2Tr currently only supports weight_decay=0")
        
        defaults = dict(
            mode=mode,
            lr=lr, lr_final=lr_final, lr_pre_warmup=lr_pre_warmup, betas=betas, eps=eps,
            warmup_steps=warmup_steps, max_steps=max_steps
        )
        super(Fused3DGS2Tr, self).__init__(params, defaults)
    
    @torch.no_grad()
    def step(self, closure=None):
        """Performs a single optimization step.
        
        Args:
            closure: A closure that reevaluates the model and returns the loss.
        """
        loss = None
        if closure is not None:
            with torch.enable_grad():
                loss = closure()
        
        for group in self.param_groups:
            mode = group['mode']
            beta1, beta2 = group['betas']
            lr = group['lr']
            lr_final = group['lr_final']
            lr_pre_warmup = group['lr_pre_warmup']
            eps = group['eps']
            warmup_steps = group['warmup_steps']
            max_steps = group['max_steps']
            
            for p in group['params']:
                if p.grad is None:
                    continue
                
                # if mode == "scale":
                #     print(p.grad.mean().item(), p.gradr.mean().item(), p.hess.mean().item())
                
                state = self.state[p]
                
                # State initialization
                if len(state) == 0:
                    state['step1'] = 0
                    state['step2'] = 0
                    # Exponential moving average of gradient values
                    state['exp_avg'] = torch.zeros_like(p, memory_format=torch.contiguous_format)
                    # Exponential moving average of squared gradient values
                    state['exp_avg_sq'] = torch.zeros_like(p, memory_format=torch.contiguous_format)
                p.exp_avg = state['exp_avg']
                p.exp_avg_sq = state['exp_avg_sq']
                
                state['step1'] += 1
                state['step2'] += 1
                scheduled_lr = get_scheduled_lr(group, state['step1'])

                additional_params = []
                if mode == "mean":
                    additional_params = [p.optim_info['scales'], p.optim_info['quats'], p.optim_info['opacities']]
                elif mode == "scale":
                    additional_params = [p.optim_info['opacities']]
                elif mode == "color":
                    additional_params = [p.optim_info['opacities']]
                elif mode == "quat":
                    additional_params = [p.optim_info['scales'], p.optim_info['opacities']]

                _make_lazy_cuda_func(f"fused_3dgs2tr_{mode}_optim")(
                    p[:p.optim_info['num_splats']] if hasattr(p, 'optim_info') and 'num_splats' in p.optim_info else p,
                    # p.grad,
                    p.optim_info['gradr'],  # gradient residual product
                    p.optim_info['hess'],
                    *additional_params,
                    state['exp_avg'],
                    state['exp_avg_sq'],
                    # scheduled_lr,
                    1.0,
                    # min(state['step1'] / 3000.0, 1.0) ** 2,
                    beta1,
                    beta2,
                    eps,
                    # eps_tr,
                    scheduled_lr,
                    state['step1'],
                    state['step2']
                )
                # zero_grad
                p.optim_info['gradr'].zero_()
                p.optim_info['hess'].zero_()

        return loss


@dataclass
class FusedNewtonOptimizerConfig:
    """Basic optimizer config with Newton"""

    _target: Type = Fused3DGS2Tr

    lr: float = 1e-6
    lr_final: Optional[float] = None
    lr_pre_warmup: float = 0.0
    warmup_steps: Optional[int] = None
    max_steps: Optional[int] = 30000

    betas: tuple = (0.9, 0.999)
    eps: float = 1e-8

    mode: Optional[Literal["mean", "scale", "opacity", "quat"]] = None


def create_optimizers(model: torch.nn.Module, config: Dict):
    optimizers = {}
    for param_key, param in model.get_param_groups().items():
        for key, optim in config.items():
            if key != param_key:
                continue
            optim = optim._target(param, **asdict(optim))
            if param_key in optimizers:
                raise RuntimeError(f"Ambiguous optimizer names for {param_key}")
            optimizers[param_key] = optim
            break
        if param_key not in optimizers:
            raise RuntimeError(f"No optimizer found for {param_key}")
    return optimizers


@dataclass
class OptimizerConfig:
    max_steps: Optional[int] = None
    use_scale_agnostic_mean: bool = True
    use_per_splat_bias_correction: bool = False

    # MCMC
    means_lr: float = 1.6e-4
    means_lr_final: Optional[float] = 1.6e-6
    scales_lr: float = 0.005
    scales_lr_final: Optional[float] = None
    quats_lr: float = 0.0005
    opacities_lr: float = 0.05
    features_dc_lr: float = 0.0025
    features_sh_lr: float = 0.0025 / 20

    # MRNF
    # means_lr: float = 1.28e-4
    # means_lr_final: Optional[float] = 1.6e-7
    # scales_lr: float = 0.02
    # scales_lr_final: Optional[float] = 0.005
    # quats_lr: float = 0.0015
    # opacities_lr: float = 0.025
    # features_dc_lr: float = 0.005
    # features_sh_lr: float = 0.005 / 20

    features_ch_lr: float = 0.0025 / 5
    sv_sites_lr: float = 0.01
    sv_colors_lr: float = 0.0005
    densities_lr: float = 0.05
    densities_lr_final: Optional[float] = 0.0005
    background_color_lr: float = 0.0025
    background_sh_lr: float = 0.0025 / 5
    bilagrid_lr: float = 2e-3
    bilagrid_lr_final: Optional[float] = 1e-4
    bilagrid_lr_warmup: int = 1000
    bilagrid_depth_lr: float = 2e-3
    bilagrid_depth_lr_final: Optional[float] = 1e-4
    bilagrid_depth_lr_warmup: int = 2000
    bilagrid_normal_lr: float = 5e-4
    bilagrid_normal_lr_final: Optional[float] = 4e-5
    bilagrid_normal_lr_warmup: int = 2000
    ppisp_lr: float = 2e-3
    ppisp_lr_final: Optional[float] = 2e-5
    ppisp_lr_warmup: int = 500   # TODO: pre-warmup
    camera_opt_lr: float = 1e-4
    camera_opt_lr_final: Optional[float] = 5e-7
    camera_opt_lr_warmup: int = 1000
