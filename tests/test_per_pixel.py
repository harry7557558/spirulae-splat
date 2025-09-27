import pytest
import traceback
import torch
import torch.nn.functional as F

import spirulae_splat.modules.per_pixel as module

from test_render_background_sh import timeit

torch.manual_seed(42)

device = torch.device("cuda:0")


def check_close(name, a, b, atol=1e-5, rtol=1e-5):
    print(name, [*a.shape], [*b.shape], a.dtype, b.dtype)
    try:
        torch.testing.assert_close(a, b, atol=atol, rtol=rtol)
    except AssertionError:
        traceback.print_exc()
        diff = torch.abs(a - b).detach()
        print(f"{diff.max()=} {diff.float().mean()=}")


def _torch_blend_background(rgb, alpha, background):
    return torch.clip(rgb + (1.0 - alpha) * background, 0.0, 1.0)


@pytest.mark.skipif(not torch.cuda.is_available(), reason="No CUDA device")
def test_blend_background():
    w, h = 1440, 1080

    rgb = torch.randn((h, 4, w), device=device).transpose(-1, -2)[..., :3]
    alpha = torch.randn((h, w+100, 1), device=device)[:, :w, :]
    background = torch.randn((h, w+100, 3), device=device)[:, -w:, :]
    # rgb, alpha, background = rgb.contiguous(), alpha.contiguous(), background.contiguous()

    _rgb = torch.nn.Parameter(rgb.clone())
    _alpha = torch.nn.Parameter(alpha.clone())
    _background = torch.nn.Parameter(background.clone())

    rgb = torch.nn.Parameter(rgb)
    alpha = torch.nn.Parameter(alpha)
    background = torch.nn.Parameter(background)

    forward = lambda: module.blend_background(rgb, alpha, background)
    _forward = lambda: _torch_blend_background(_rgb, _alpha, _background)
    output = forward()
    _output = _forward()

    print('rgb.is_contiguous():', rgb.is_contiguous())
    print('_rgb.is_contiguous():', _rgb.is_contiguous())
    print('alpha.is_contiguous():', alpha.is_contiguous())
    print('background.is_contiguous():', background.is_contiguous())
    print('output.is_contiguous():', output.is_contiguous())
    print()

    print("test forward")
    check_close("output", output, _output)
    print()

    weights = torch.randn_like(output)
    def loss_fun(output):
        return (weights*output).sum()
    loss = loss_fun(output)
    _loss = loss_fun(_output)

    loss.backward(retain_graph=True)
    _loss.backward(retain_graph=True)

    print("test backward")
    check_close('rgb.grad', rgb.grad, _rgb.grad)
    check_close('alpha.grad', alpha.grad, _alpha.grad)
    check_close('background.grad', background.grad, _background.grad)
    print()

    repeat = 100

    print("profile forward")
    timeit(_forward, "torch forward", repeat=repeat)
    timeit(forward, "fused forward", repeat=repeat)
    print()

    loss.backward(retain_graph=True)
    _loss.backward(retain_graph=True)

    print("profile backward")
    timeit(lambda: loss.backward(retain_graph=True), "torch backward", repeat=repeat)
    timeit(lambda: _loss.backward(retain_graph=True), "fused backward", repeat=repeat)
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
    test_blend_background()
    # profile_compute_per_splat_losses()
