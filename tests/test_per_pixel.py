import pytest
import torch
import torch.nn.functional as F

import spirulae_splat.splat.cuda._wrapper_per_pixel as module

from test_compute_per_splat_losses import Config
from spirulae_splat.modules.training_losses import SplatTrainingLosses

from utils import check_close, timeit

torch.manual_seed(42)

device = torch.device("cuda:0")


def _torch_blend_background(rgb, alpha, background):
    return torch.clip(rgb + (1.0 - alpha) * background, 0.0, 1.0)


@pytest.mark.skipif(not torch.cuda.is_available(), reason="No CUDA device")
def test_blend_background():
    b, w, h = 3, 1440, 1080

    rgb = torch.randn((b, h, 4, w), device=device).transpose(-1, -2)[..., :3]
    alpha = torch.randn((b, h, w+100, 1), device=device)[:, :, :w, :]
    background = torch.randn((b, h, w+100, 3), device=device)[:, :, -w:, :]
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
    timeit(lambda: _loss.backward(retain_graph=True), "torch backward", repeat=repeat)
    timeit(lambda: loss.backward(retain_graph=True), "fused backward", repeat=repeat)
    print()


@pytest.mark.skipif(not torch.cuda.is_available(), reason="No CUDA device")
def test_depth_to_normal():
    b, w, h = 4, 144, 108

    import spirulae_splat.splat.utils

    depths = torch.exp(0.1*torch.randn((b, h, w, 1), device=device))
    intrins = 0.5*(w*h)**0.5 * torch.exp(0.2*torch.randn((b, 4), device=device))

    _depths = torch.nn.Parameter(depths.clone())
    depths = torch.nn.Parameter(depths)

    forward = lambda: module.depth_to_normal(depths, "pinhole", intrins)
    _forward = lambda: torch.stack([
        spirulae_splat.splat.utils.depth_to_normal(_depths[i],
        spirulae_splat.splat.utils._Camera(h, w, "OPENCV", (intrins[i,0], intrins[i,1], intrins[i,2], intrins[i,3])),
        z_depth=False)
        for i in range(b)])
    output = forward()
    _output = _forward()

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
    check_close('depths.grad', depths.grad, _depths.grad)
    print()

    repeat = 100

    print("profile forward")
    timeit(_forward, "torch forward", repeat=repeat)
    timeit(forward, "fused forward", repeat=repeat)
    print()

    loss.backward(retain_graph=True)
    _loss.backward(retain_graph=True)

    print("profile backward")
    timeit(lambda: _loss.backward(retain_graph=True), "torch backward", repeat=repeat)
    timeit(lambda: loss.backward(retain_graph=True), "fused backward", repeat=repeat)
    print()



if __name__ == "__main__":
    # test_blend_background()
    test_depth_to_normal()
