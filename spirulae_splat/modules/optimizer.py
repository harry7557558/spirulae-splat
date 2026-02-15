from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Type

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
        
        defaults = dict(lr=lr, betas=betas, eps=eps)
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
            
            # Collect all tensors for this group
            params_with_grad = []
            grads = []
            exp_avgs = []
            exp_avg_sqs = []
            
            for p in group['params']:
                if p.grad is None:
                    continue
                
                params_with_grad.append(p)
                grads.append(p.grad)
                
                state = self.state[p]
                
                # State initialization
                if len(state) == 0:
                    state['step'] = 0
                    # Exponential moving average of gradient values
                    state['exp_avg'] = torch.zeros_like(p, memory_format=torch.contiguous_format)
                    # Exponential moving average of squared gradient values
                    state['exp_avg_sq'] = torch.zeros_like(p, memory_format=torch.contiguous_format)
                
                exp_avgs.append(state['exp_avg'])
                exp_avg_sqs.append(state['exp_avg_sq'])
                
                state['step'] += 1
            
            if len(params_with_grad) == 0:
                continue
            
            # Use the step count from the first parameter (all should be the same)
            step = self.state[params_with_grad[0]]['step']
            
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



class Fused3DGS2TrMean(Optimizer):
    """
    Fully fused CUDA implementation of 3DGS^2-TR optimizer for Gaussian means.
    https://arxiv.org/abs/2602.00395
    """
    
    def __init__(
        self,
        params,
        lr: float = 1e-6,
        betas: tuple = (0.9, 0.999),
        eps: float = 1e-8,
        # eps_tr: float = 1e-6,
        **kwargs
    ):
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
            raise NotImplementedError("Fused3DGS2TrMean currently only supports weight_decay=0")
        
        # defaults = dict(lr=lr, betas=betas, eps=eps, eps_tr=eps_tr)
        defaults = dict(lr=lr, betas=betas, eps=eps)
        super(Fused3DGS2TrMean, self).__init__(params, defaults)
    
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
            # eps_tr = group['eps_tr']
            
            for p in group['params']:
                if p.grad is None:
                    continue
                
                # print(p.grad.mean().item(), p.gradr.mean().item(), p.hess.mean().item())
                
                state = self.state[p]
                
                # State initialization
                if len(state) == 0:
                    state['step1'] = 0
                    state['step2'] = 0
                    # Exponential moving average of gradient values
                    state['exp_avg'] = torch.zeros_like(p, memory_format=torch.contiguous_format)
                    # Exponential moving average of squared gradient values
                    state['exp_avg_sq'] = torch.zeros_like(p, memory_format=torch.contiguous_format)
                
                state['step1'] += 1
                state['step2'] += 1

                _make_lazy_cuda_func("fused_3dgs2tr_mean_optim")(
                    p,
                    # p.grad,
                    p.gradr,  # gradient residual product
                    p.hess,
                    p.scales,
                    p.quats,
                    p.opacities,
                    state['exp_avg'],
                    state['exp_avg_sq'],
                    # lr,
                    1.0,
                    beta1,
                    beta2,
                    eps,
                    # eps_tr,
                    lr,
                    state['step1'],
                    state['step2']
                )
                # zero_grad
                del p.gradr
                del p.hess

        return loss


@dataclass
class FusedNewtonOptimizerConfig(OptimizerConfig):
    """Basic optimizer config with Newton"""

    _target: Type = Fused3DGS2TrMean
