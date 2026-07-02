from dataclasses import dataclass, asdict
from typing import Any, Dict, List, Optional, Type, Literal


@dataclass
class OptimizerConfig:
    max_steps: Optional[int] = None
    use_scale_agnostic_mean: bool = True
    use_per_splat_bias_correction: bool = True

    # MCMC
    # means_lr: float = 1.6e-4
    # means_lr_final: Optional[float] = 1.6e-6
    # scales_lr: float = 0.005
    # scales_lr_final: Optional[float] = None
    # quats_lr: float = 0.001
    # opacities_lr: float = 0.05
    # features_dc_lr: float = 0.0025
    # features_sh_lr: float = 0.0025 / 20

    # MRNF
    means_lr: float = 1.28e-4
    # means_lr_final: Optional[float] = 1.6e-7
    means_lr_final: Optional[float] = 1.6e-6
    scales_lr: float = 0.02
    scales_lr_final: Optional[float] = 0.005
    quats_lr: float = 0.0015
    opacities_lr: float = 0.025
    features_dc_lr: float = 0.005
    features_sh_lr: float = 0.005 / 20

    background_dc_lr: float = 0.0025
    background_sh_lr: float = 0.0025 / 5
    bilagrid_lr: float = 2e-3  # 2e-3*sqrt(B) in paper
    bilagrid_lr_final: Optional[float] = 1e-4
    bilagrid_lr_warmup: int = 1000
    bilagrid_depth_lr: float = 2e-3
    bilagrid_depth_lr_final: Optional[float] = 1e-4
    bilagrid_depth_lr_warmup: int = 2000
    bilagrid_normal_lr: float = 5e-4
    bilagrid_normal_lr_final: Optional[float] = 4e-5
    bilagrid_normal_lr_warmup: int = 2000
    bilagrid_adagrad_lr: float = 4e-2
    bilagrid_adagrad_depth_lr: float = 4e-2
    bilagrid_adagrad_normal_lr: float = 1e-2
    ppisp_lr: float = 2e-3
    ppisp_lr_final: Optional[float] = 2e-5
    ppisp_lr_warmup: int = 500   # TODO: pre-warmup
    ppisp_adagrad_lr: float = 1e-1
    camera_opt_lr: float = 1e-4
    camera_opt_lr_final: Optional[float] = 5e-7
    camera_opt_lr_warmup: int = 1000

    def get_scheduled_lr(optim_config, name: str, step: int, max_steps: int):
        lr = getattr(optim_config, name+"_lr")
        lr_final = getattr(optim_config, name+"_lr_final", lr)
        warmup = getattr(optim_config, name+"_lr_warmup", None)
        pre_warmup = 0.0  # TODO
        scheduled_lr = lr
        if lr_final is not None and lr != 0.0 and lr_final != 0.0:
            scheduled_lr = lr * (lr_final / lr) ** min(step / max_steps, 1.0)
        if warmup is not None:
            scheduled_lr = min(scheduled_lr, pre_warmup + (lr - pre_warmup) * min(step / warmup, 1.0))
        return scheduled_lr
