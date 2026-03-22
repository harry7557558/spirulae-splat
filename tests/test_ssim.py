import pytest
import torch

from utils import check_close, timeit


from fused_ssim import fused_ssim
from spirulae_splat.modules.training_losses import FusedSSIM


@pytest.mark.skipif(not torch.cuda.is_available(), reason="No CUDA device")
def test_ssim():

    torch.manual_seed(42)

    num_images = 5
    h, w = 42, 67
    image = torch.rand((num_images, h, w, 3), device="cuda")
    ref = torch.rand((num_images, h, w, 3), device="cuda")

    print("Test forward")

    def forward_ssim(image1, image2):
        return fused_ssim(image1.permute(0, 3, 1, 2), image2.permute(0, 3, 1, 2))

    def my_forward_ssim(image1, image2):
        return FusedSSIM.apply(image1, image2)[0]

    image1 = torch.nn.Parameter(image.clone())
    image2 = torch.nn.Parameter(image.clone())
    out1 = forward_ssim(image1, ref)
    out2 = my_forward_ssim(image2, ref)
    check_close("forward", out1, out2, atol=1e-5, rtol=1e-5)
    print()

    print("Test backward")
    weight_out = torch.randn_like(out1)
    def backward_ssim(image, out):
        (out * weight_out).sum().backward(retain_graph=True)
        return image.grad
    def my_backward_ssim(image, out):
        (out * weight_out).sum().backward(retain_graph=True)
        return image.grad
    v_image1 = backward_ssim(image1, out1)
    v_image2 = my_backward_ssim(image2, out2)
    check_close("image", v_image1, v_image2, atol=1e-5, rtol=1e-5)
    print()


@pytest.mark.skipif(not torch.cuda.is_available(), reason="No CUDA device")
def test_memory_efficient_ssim():

    torch.manual_seed(42)

    num_images = 5
    h, w = 42, 67
    image = torch.rand((num_images, h, w, 3), device="cuda")
    ref = torch.rand((num_images, h, w, 3), device="cuda")

    print("Test forward (memory efficient)")

    def forward_ssim(image1, image2):
        return fused_ssim(image1.permute(0, 3, 1, 2), image2.permute(0, 3, 1, 2))

    def my_forward_ssim(image1, image2):
        return FusedSSIM.apply(image1, image2, False)[0]

    image1 = torch.nn.Parameter(image.clone())
    image2 = torch.nn.Parameter(image.clone())
    out1 = forward_ssim(image1, ref)
    out2 = my_forward_ssim(image2, ref)
    check_close("forward", out1, out2, atol=1e-5, rtol=1e-5)
    print()

    print("Test backward (memory efficient)")
    weight_out = torch.randn_like(out1)
    def backward_ssim(image, out):
        (out * weight_out).sum().backward(retain_graph=True)
        return image.grad
    def my_backward_ssim(image, out):
        (out * weight_out).sum().backward(retain_graph=True)
        return image.grad
    v_image1 = backward_ssim(image1, out1)
    v_image2 = my_backward_ssim(image2, out2)
    check_close("image", v_image1, v_image2, atol=1e-5, rtol=1e-5)
    print()

    # import matplotlib.pyplot as plt
    # diff = v_image1 - v_image2
    # diff = torch.log10(torch.abs(diff) + 1e-12)
    # plt.imshow(diff[0, ..., 0].cpu().numpy())
    # plt.show()


@pytest.mark.skipif(not torch.cuda.is_available(), reason="No CUDA device")
def benchmark_ssim():

    torch.manual_seed(42)

    num_images = 1
    h, w = 3840, 2160
    image = torch.rand((num_images, h, w, 3), device="cuda")
    ref = torch.rand((num_images, h, w, 3), device="cuda")

    def forward_ssim(image1, image2):
        return fused_ssim(image1.permute(0, 3, 1, 2), image2.permute(0, 3, 1, 2))

    def my_forward_ssim(image1, image2):
        return FusedSSIM.apply(image1, image2)[0]

    def my_mem_eff_forward_ssim(image1, image2):
        return FusedSSIM.apply(image1, image2, False)[0]

    def backward_ssim(image, out):
        (out * weight_out).sum().backward(retain_graph=True)
        return image.grad
    def my_backward_ssim(image, out):
        (out * weight_out).sum().backward(retain_graph=True)
        return image.grad
    def my_mem_eff_backward_ssim(image, out):
        (out * weight_out).sum().backward(retain_graph=True)
        return image.grad

    image = torch.rand((num_images, h, w, 3), device="cuda")
    ref = torch.rand((num_images, h, w, 3), device="cuda")
    image1 = torch.nn.Parameter(image.clone())
    image2 = torch.nn.Parameter(image.clone())
    image3 = torch.nn.Parameter(image.clone())
    out1 = forward_ssim(image1, ref)
    out2 = my_forward_ssim(image2, ref)
    out3 = my_mem_eff_forward_ssim(image3, ref)
    weight_out = torch.randn_like(out1)

    print("Benchmark forward")
    timeit(lambda: forward_ssim(image, ref), "fused_ssim")
    timeit(lambda: my_forward_ssim(image, ref), "spirulae")
    timeit(lambda: my_mem_eff_forward_ssim(image, ref), "spirulae_mem_eff")
    print()

    print("Benchmark backward")
    timeit(lambda: backward_ssim(image1, out1), "fused_ssim")
    timeit(lambda: my_backward_ssim(image2, out2), "spirulae")
    timeit(lambda: my_mem_eff_backward_ssim(image3, out3), "spirulae_mem_eff")
    print()



if __name__ == "__main__":
    test_ssim()
    test_memory_efficient_ssim()
    benchmark_ssim()
