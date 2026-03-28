from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Type, Literal

import torch
from torch.optim.optimizer import Optimizer

from nerfstudio.engine.optimizers import OptimizerConfig

from spirulae_splat.splat.cuda import (
    _C,
    _make_lazy_cuda_func,
)

class FusedAdam(Optimizer):
    """
    Fully fused CUDA implementation of Adam optimizer.
    """
    
    def __init__(
        self,
        params,
        lr: float = 1e-3,
        betas: tuple = (0.9, 0.999),
        eps: float = 1e-8,
        tr: Optional[float] = 1e-6,
        tr_final: Optional[float] = 1e-8,
        tr_max_steps: Optional[int] = 30000,
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
        
        defaults = dict(lr=lr, betas=betas, eps=eps, tr=tr, tr_final=tr_final, tr_max_steps=tr_max_steps)
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
            eps = group['eps']

            tr = group['tr']
            tr_final = group['tr_final']
            tr_max_steps = group['tr_max_steps']
            
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
                eps_tr = tr * (tr_final / tr) ** min(state['step'] / tr_max_steps, 1.0)
                if p.optim_info['optimizer_override'] == "fused_adam_linear_rgb_optim":
                    pass
                elif p.optim_info['optimizer_override'] in ["fused_adamtr_linear_rgb_optim", "fused_adamtr_rgb_optim"]:
                    additional_params = [p.optim_info['opacities']]
                    eps_params.append(eps_tr)
                elif p.optim_info['optimizer_override'] in ["fused_adamtr_linear_rgb_sh_optim", "fused_adamtr_rgb_sh_optim"]:
                    additional_params = [p.optim_info['features_dc'], p.optim_info['opacities']]
                    eps_params.append(eps_tr)

                _make_lazy_cuda_func(p.optim_info['optimizer_override'])(
                    p[:p.optim_info['num_splats']] if hasattr(p, 'optim_info') and 'num_splats' in p.optim_info else p,
                    p.grad,
                    state['exp_avg'],
                    state['exp_avg_sq'],
                    *additional_params,
                    lr,
                    beta1,
                    beta2,
                    *eps_params,
                    state['step']
                )

            if len(params_with_grad) == 0:
                continue
            
            # Use the step count from the first parameter (all should be the same)
            step = steps[0]
            
            # Launch fused CUDA kernel for all tensors
            _make_lazy_cuda_func("fused_adam_multi")(
                params_with_grad,
                grads,
                exp_avgs,
                exp_avg_sqs,
                lr,
                beta1,
                beta2,
                eps,
                step
            )
        
        return loss


@dataclass
class FusedAdamOptimizerConfig(OptimizerConfig):
    """Basic optimizer config with Adam"""

    _target: Type = FusedAdam

    tr: Optional[float] = 1e-6
    tr_final: Optional[float] = 1e-8
    tr_max_steps: Optional[int] = 30000



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
        betas: tuple = (0.9, 0.999),
        eps: float = 1e-8,
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
        
        # defaults = dict(lr=lr, betas=betas, eps=eps, eps_tr=eps_tr)
        defaults = dict(mode=mode, lr=lr, betas=betas, eps=eps)
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
            eps = group['eps']
            # eps_tr = group['eps_tr']
            # lr = 1e10
            
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
                    # lr,
                    1.0,
                    # min(state['step1'] / 3000.0, 1.0) ** 2,
                    beta1,
                    beta2,
                    eps,
                    # eps_tr,
                    lr,
                    state['step1'],
                    state['step2']
                )
                # zero_grad
                del p.optim_info['gradr']
                del p.optim_info['hess']

        return loss


@dataclass
class FusedNewtonOptimizerConfig(OptimizerConfig):
    """Basic optimizer config with Newton"""

    _target: Type = Fused3DGS2Tr

    mode: Optional[Literal["mean", "scale", "opacity", "quat"]] = None
