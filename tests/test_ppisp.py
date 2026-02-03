import pytest
import torch

from utils import check_close, timeit


from ppisp import PPISP, PPISPConfig
from spirulae_splat.modules.training_losses import apply_ppisp, _ComputePPISPRegularization


@pytest.mark.skipif(not torch.cuda.is_available(), reason="No CUDA device")
def test_ppisp():

    torch.manual_seed(42)

    num_images = 5
    h, w = 42, 67
    image = torch.rand((num_images, h, w, 3), device="cuda")

    ppisp = PPISP(num_images, num_images).cuda()
    ppisp.config.use_controller = False
    ppisp.exposure_params.data = 0.1 * torch.randn_like(ppisp.exposure_params.data)
    ppisp.vignetting_params.data = 0.5 * torch.randn_like(ppisp.vignetting_params.data)
    ppisp.color_params.data = 2.0 * torch.randn_like(ppisp.color_params.data)
    ppisp.crf_params.data += 0.5 * torch.randn_like(ppisp.crf_params.data)
    # ppisp.exposure_params.data = 0.0 * torch.randn_like(ppisp.exposure_params.data)
    # ppisp.vignetting_params.data = 0.0 * torch.randn_like(ppisp.vignetting_params.data)
    # ppisp.color_params.data = 0.0 * torch.randn_like(ppisp.color_params.data)
    # ppisp.crf_params.data += 0.0 * torch.randn_like(ppisp.crf_params.data)

    ppisp_params = torch.nn.Parameter(torch.concatenate([
        ppisp.exposure_params.reshape((num_images, -1)),
        ppisp.vignetting_params.reshape((num_images, -1)),
        ppisp.color_params.reshape((num_images, -1)),
        ppisp.crf_params.reshape((num_images, -1))
    ], dim=-1).clone())
    # print(ppisp_params.detach().cpu().numpy().tolist())

    intrins = torch.tensor([[50.0, 50.0, w/2, h/2]], device="cuda").repeat(num_images, 1)

    print("Test forward")

    def forward_ppisp(image):
        pixel_coords = torch.stack(torch.meshgrid(
            torch.arange(h, device=image.device),
            torch.arange(w, device=image.device),
            indexing='ij'
        ), dim=-1).unsqueeze(0).repeat(len(image), 1, 1, 1)  # [B, H, W, 2]
        pixel_coords = torch.flip(pixel_coords, dims=[-1])
        camera_idx = torch.arange(len(image), device=image.device)
        frame_idx = torch.arange(len(image), device=image.device)
        return torch.stack([
            ppisp.forward(
                image[i], pixel_coords[i], (w, h), camera_idx[i], frame_idx[i]
            ) for i in range(len(image))
        ])

    def my_forward_ppisp(image):
        return apply_ppisp(image, ppisp_params[:len(image)], intrins, w, h)#.clamp(0.0, 1.0)

    image1 = torch.nn.Parameter(image.clone())
    image2 = torch.nn.Parameter(image.clone())
    out1 = forward_ppisp(image1)
    out2 = my_forward_ppisp(image2)
    check_close("forward", out1, out2, atol=1e-5, rtol=1e-5)
    print()

    # import matplotlib.pyplot as plt
    # plt.plot(image[0, h//2, :, 0].cpu().numpy(), out1[0, h//2, :, 0].detach().cpu().numpy(), 'o', label='forward_ppisp')
    # plt.plot(image[0, h//2, :, 0].cpu().numpy(), out2[0, h//2, :, 0].detach().cpu().numpy(), 'x', label='my_forward_ppisp')
    # plt.legend()
    # plt.savefig("/mnt/d/plot.png")
    # plt.close()

    print("Test backward")
    weight_out = torch.randn_like(out1)
    def backward_ppisp(image, out):
        (out * weight_out).sum().backward(retain_graph=True)
        return image.grad, torch.concatenate([
            ppisp.exposure_params.grad.reshape((num_images, -1)),
            ppisp.vignetting_params.grad.reshape((num_images, -1)),
            ppisp.color_params.grad.reshape((num_images, -1)),
            ppisp.crf_params.grad.reshape((num_images, -1))
        ], dim=-1)
    def my_backward_ppisp(image, out):
        (out * weight_out).sum().backward(retain_graph=True)
        return image.grad, ppisp_params.grad
    v_image1, v_param1 = backward_ppisp(image1, out1)
    v_image2, v_param2 = my_backward_ppisp(image2, out2)
    check_close("image", v_image1, v_image2, atol=1e-5, rtol=1e-5)
    # check_close("params", v_param1, v_param2, atol=1e-5, rtol=1e-5)
    args = {'atol':1e-4, 'rtol':1e-4}
    check_close("exposure params", ppisp.exposure_params.grad.flatten(), ppisp_params.grad[:, :1].flatten(), **args)
    check_close("vignetting params", ppisp.vignetting_params.grad.flatten(), ppisp_params.grad[:, 1:16].flatten(), **args)
    check_close("color params", ppisp.color_params.grad.flatten(), ppisp_params.grad[:, 16:24].flatten(), **args)
    check_close("crf params", ppisp.crf_params.grad.flatten(), ppisp_params.grad[:, 24:36].flatten(), **args)
    print()


@pytest.mark.skipif(not torch.cuda.is_available(), reason="No CUDA device")
def benchmark_ppisp():

    torch.manual_seed(42)

    num_images = 1
    h, w = 3840, 2160
    image = torch.rand((num_images, h, w, 3), device="cuda")

    ppisp = PPISP(num_images, num_images).cuda()
    ppisp.config.use_controller = False
    ppisp.exposure_params.data = 0.1 * torch.randn_like(ppisp.exposure_params.data)
    ppisp.vignetting_params.data = 0.5 * torch.randn_like(ppisp.vignetting_params.data)
    ppisp.color_params.data = 2.0 * torch.randn_like(ppisp.color_params.data)
    ppisp.crf_params.data += 0.5 * torch.randn_like(ppisp.crf_params.data)

    ppisp_params = torch.nn.Parameter(torch.concatenate([
        ppisp.exposure_params.reshape((num_images, -1)),
        ppisp.vignetting_params.reshape((num_images, -1)),
        ppisp.color_params.reshape((num_images, -1)),
        ppisp.crf_params.reshape((num_images, -1))
    ], dim=-1).clone())

    intrins = torch.tensor([[-1.0, -1.0, w/2, h/2]], device="cuda").repeat(num_images, 1)

    pixel_coords = torch.stack(torch.meshgrid(
        torch.arange(h, device=image.device),
        torch.arange(w, device=image.device),
        indexing='ij'
    ), dim=-1).unsqueeze(0).repeat(len(image), 1, 1, 1)  # [B, H, W, 2]
    pixel_coords = torch.flip(pixel_coords, dims=[-1])
    camera_idx = torch.arange(len(image), device=image.device)
    frame_idx = torch.arange(len(image), device=image.device)
    def forward_ppisp(image):
        return ppisp.forward(
            image[0], pixel_coords[0], (w, h), camera_idx[0], frame_idx[0]
        )[None]

    def my_forward_ppisp(image):
        return apply_ppisp(image, ppisp_params[:len(image)], intrins, w, h)#.clamp(0.0, 1.0)

    def backward_ppisp(image, out):
        (out * weight_out).sum().backward(retain_graph=True)
        return image.grad, torch.concatenate([
            ppisp.exposure_params.grad.reshape((num_images, -1)),
            ppisp.vignetting_params.grad.reshape((num_images, -1)),
            ppisp.color_params.grad.reshape((num_images, -1)),
            ppisp.crf_params.grad.reshape((num_images, -1))
        ], dim=-1)

    def my_backward_ppisp(image, out):
        (out * weight_out).sum().backward(retain_graph=True)
        return image.grad, ppisp_params.grad

    image = torch.rand((num_images, h, w, 3), device="cuda")
    image1 = torch.nn.Parameter(image.clone())
    image2 = torch.nn.Parameter(image.clone())
    out1 = forward_ppisp(image1)
    out2 = my_forward_ppisp(image2)
    weight_out = torch.randn_like(out1)

    print("Benchmark forward")
    timeit(lambda: forward_ppisp(image), "ppisp")
    timeit(lambda: my_forward_ppisp(image), "spirulae")
    print()

    print("Benchmark backward")
    timeit(lambda: backward_ppisp(image1, out1), "ppisp")
    timeit(lambda: my_backward_ppisp(image2, out2), "spirulae")
    print()


@pytest.mark.skipif(not torch.cuda.is_available(), reason="No CUDA device")
def test_ppisp_reg():

    torch.manual_seed(42)

    # num_images = 10
    num_images = 4567

    ppisp = PPISP(num_images, num_images).cuda()
    ppisp.config.use_controller = False
    ppisp.exposure_params.data = 0.1 * torch.randn_like(ppisp.exposure_params.data)
    ppisp.vignetting_params.data = 0.5 * torch.randn_like(ppisp.vignetting_params.data)
    ppisp.color_params.data = 2.0 * torch.randn_like(ppisp.color_params.data)
    ppisp.crf_params.data += 0.5 * torch.randn_like(ppisp.crf_params.data)

    # ppisp.config.exposure_mean *= 0.0
    # ppisp.config.vig_center *= 0.0
    # ppisp.config.vig_non_pos *= 0.0
    # ppisp.config.vig_channel *= 0.0
    # ppisp.config.color_mean *= 0.0
    # ppisp.config.crf_channel *= 0.0

    ppisp_params = torch.nn.Parameter(torch.concatenate([
        ppisp.exposure_params.reshape((num_images, -1)),
        ppisp.vignetting_params.reshape((num_images, -1)),
        ppisp.color_params.reshape((num_images, -1)),
        ppisp.crf_params.reshape((num_images, -1))
    ], dim=-1).clone())

    print("Test forward")

    def forward_ppisp():
        return ppisp.get_regularization_loss()

    def my_forward_ppisp():
        ppisp_loss_weights = [
            ppisp.config.exposure_mean,
            ppisp.config.vig_center,
            ppisp.config.vig_non_pos,
            ppisp.config.vig_channel,
            ppisp.config.color_mean,
            ppisp.config.crf_channel,
        ]
        return _ComputePPISPRegularization.apply(
            ppisp_params,
            ppisp_loss_weights
        ).sum()

    loss1 = forward_ppisp()
    loss2 = my_forward_ppisp()
    check_close("forward", loss1, loss2, atol=1e-5, rtol=1e-5)
    print()

    print("Test backward")
    def backward_ppisp(loss):
        loss.backward(retain_graph=True)
        return torch.concatenate([
            ppisp.exposure_params.grad.reshape((num_images, -1)),
            ppisp.vignetting_params.grad.reshape((num_images, -1)),
            ppisp.color_params.grad.reshape((num_images, -1)),
            ppisp.crf_params.grad.reshape((num_images, -1))
        ], dim=-1)
    def my_backward_ppisp(loss):
        loss.backward(retain_graph=True)
        return ppisp_params.grad
    v_param1 = backward_ppisp(loss1)
    v_param2 = my_backward_ppisp(loss2)
    # check_close("params", v_param1, v_param2, atol=1e-5, rtol=1e-5)
    args = {'atol':1e-4, 'rtol':1e-4}
    check_close("exposure params", ppisp.exposure_params.grad.flatten(), ppisp_params.grad[:, :1].flatten(), **args)
    check_close("vignetting params", ppisp.vignetting_params.grad.flatten(), ppisp_params.grad[:, 1:16].flatten(), **args)
    check_close("color params", ppisp.color_params.grad.flatten(), ppisp_params.grad[:, 16:24].flatten(), **args)
    check_close("crf params", ppisp.crf_params.grad.flatten(), ppisp_params.grad[:, 24:36].flatten(), **args)
    print()

    print("Benchmark forward")
    timeit(lambda: forward_ppisp(), "ppisp")
    timeit(lambda: my_forward_ppisp(), "spirulae")
    print()

    print("Benchmark backward")
    timeit(lambda: backward_ppisp(loss1), "ppisp")
    timeit(lambda: my_backward_ppisp(loss2), "spirulae")
    print()


if __name__ == "__main__":
    test_ppisp()
    benchmark_ppisp()
    test_ppisp_reg()
