from dataclasses import dataclass
import pytest
import torch
import torch.nn.functional as F

from spirulae_splat.modules.training_losses import SplatTrainingLosses
import spirulae_splat.splat.cuda as _C

from utils import check_close, timeit

torch.manual_seed(42)

device = torch.device("cuda:0")


@dataclass
class Config:
    mcmc_opacity_reg = 1.1
    mcmc_scale_reg = 1.2
    max_gauss_ratio = 2.0
    scale_regularization_weight = 1.3
    erank_reg = 1.4
    erank_reg_s3 = 1.5
    quat_norm_reg_weight = 0.01

    use_3dgs = False
    use_bilateral_grid = False
    use_mcmc = True
    erank_reg_warmup = 0
    
    depth_distortion_depth_degree = -1
    depth_distortion_uv_degree = -1
    depth_supervision_weight = 0.0
    normal_supervision_weight = 0.0
    alpha_loss_weight = 0.0
    alpha_loss_weight_under = 0.0
    adaptive_exposure_mode = ""


@pytest.mark.skipif(not torch.cuda.is_available(), reason="No CUDA device")
def test_compute_per_splat_losses():
    num_points = 1000

    config = Config()
    torch_impl = SplatTrainingLosses(config, 0)

    scales = torch.randn((num_points, 2+int(config.use_3dgs)), device=device)
    opacities = torch.randn((num_points, 1), device=device)
    quats = torch.randn((num_points, 4), device=device)

    _scales = torch.nn.Parameter(scales.clone())
    _opacities = torch.nn.Parameter(opacities.clone())
    _quats = torch.nn.Parameter(quats.clone())

    scales = torch.nn.Parameter(scales)
    opacities = torch.nn.Parameter(opacities)
    quats = torch.nn.Parameter(quats)

    step = 15000
    keys = "mcmc_opacity_reg mcmc_scale_reg scale_reg erank_reg quat_norm_reg".split()

    output = {}
    torch_impl.get_static_losses(
        step, quats, scales, opacities,
        output, _use_torch_impl = False
    )
    _output = {}
    torch_impl.get_static_losses(
        step, _quats, _scales, _opacities,
        _output, _use_torch_impl = True
    )
    print("test forward")
    for key in keys:
        check_close(key, output[key], _output[key])

    weights = dict(zip(keys, torch.randn(5).numpy().tolist()))
    loss = sum([weights[key]*output[key] for key in keys])
    _loss = sum([weights[key]*_output[key] for key in keys])
    loss.backward()
    _loss.backward()

    print("test backward")
    check_close('scales.grad', scales.grad, _scales.grad)
    check_close('opacities.grad', opacities.grad, _opacities.grad)
    check_close('quats.grad', quats.grad, _quats.grad)

    print()


@pytest.mark.skipif(not torch.cuda.is_available(), reason="No CUDA device")
def profile_compute_per_splat_losses():
    num_points = 2_000_000

    config = Config()
    torch_impl = SplatTrainingLosses(config, 0)

    scales = torch.randn((num_points, 2+int(config.use_3dgs)), device=device)
    opacities = torch.randn((num_points, 1), device=device)
    quats = torch.randn((num_points, 4), device=device)

    _scales = torch.nn.Parameter(scales.clone())
    _opacities = torch.nn.Parameter(opacities.clone())
    _quats = torch.nn.Parameter(quats.clone())

    scales = torch.nn.Parameter(scales)
    opacities = torch.nn.Parameter(opacities)
    quats = torch.nn.Parameter(quats)

    step = 15000
    forward = lambda: torch_impl.get_static_losses(
        step, quats, scales, opacities,
        {}, _use_torch_impl = False
    )
    _forward = lambda: torch_impl.get_static_losses(
        step, _quats, _scales, _opacities,
        {}, _use_torch_impl = True
    )
    output = forward()
    _output = _forward()

    repeat = 100

    print("profile forward")
    timeit(_forward, "torch forward", repeat=repeat)
    timeit(forward, "fused forward", repeat=repeat)
    print()

    keys = "mcmc_opacity_reg mcmc_scale_reg scale_reg erank_reg quat_norm_reg".split()
    weights = dict(zip(keys, torch.randn(5).numpy().tolist()))

    def fun(output):
        return sum([weights[key]*output[key] for key in keys])
    fun(output).backward(retain_graph=True)
    fun(_output).backward(retain_graph=True)

    print("profile backward")
    timeit(lambda: fun(_output).backward(retain_graph=True), "torch backward", repeat=repeat)
    timeit(lambda: fun(output).backward(retain_graph=True), "fused backward", repeat=repeat)
    print()


if __name__ == "__main__":
    test_compute_per_splat_losses()
    profile_compute_per_splat_losses()
