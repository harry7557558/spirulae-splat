import torch
import torch.nn.functional as F

from spirulae_splat.splat.cuda import _make_lazy_cuda_func


@torch.no_grad()
def detect_edge(im: torch.Tensor):
    assert im.dtype == torch.float32
    return _make_lazy_cuda_func("canny_edge_filter")(im)
    # return _make_lazy_cuda_func("laplacian_edge_filter")(im)
    # return _make_lazy_cuda_func("smoothed_laplacian_edge_filter")(im)

    if im.shape[-1] in [3, 4]:
        im = 0.299 * im[..., 0] + 0.587 * im[..., 1] + 0.114 * im[..., 2]
        # TODO: mask
        im = im.unsqueeze(-1)
    im = im.permute(0, 3, 1, 2).contiguous(memory_format=torch.channels_last)

    kernel = torch.tensor([
        [0, 1, 0],
        [1, -4, 1],
        [0, 1, 0]
    ]).to(im).view(1, 1, 3, 3)

    im = F.conv2d(im, kernel, padding=1)

    im = im.permute(0, 2, 3, 1)
    im = torch.abs(im)

    return im


@torch.no_grad()
def detect_edge_ms(im: torch.Tensor):
    if im.shape[-1] in [3, 4]:
        im = 0.299 * im[..., 0] + 0.587 * im[..., 1] + 0.114 * im[..., 2]
        im = im.unsqueeze(-1)
    edge = detect_edge(im)
    sc = 1
    while im.shape[1] * im.shape[2] > 256 * 256:
        sc *= 2
        im = torch.nn.functional.avg_pool2d(im.permute(0, 3, 1, 2), 2).permute(0, 2, 3, 1)
        edge_i = detect_edge(im)
        edge_i = edge_i[:, :, None, :, None, :] \
            .repeat(1, 1, sc, 1, sc, 1).reshape(edge_i.shape[0], sc*edge_i.shape[1], sc*edge_i.shape[2], edge_i.shape[3])
        edge[:, :edge_i.shape[1], :edge_i.shape[2], :] += edge_i
    return edge


if __name__ == "__main__":
    import cv2
    import numpy as np
    path = "/mnt/d/gs/data/360_v2/bicycle_4/images/_DSC8783.JPG"
    # path = "/mnt/d/gs/data/Atrium_Godiva/images/00006.jpg"
    # path = "/mnt/d/gs/data/adr/20260221-queens_park/images/cam0/00030.jpg"
    # path = "/mnt/d/gs/data/adr/20251109-elevator-insta360-resolution-test/8k/images/cam0/00048.jpg"
    # path = "/mnt/d/gs/data/zipnerf/nyc/images/DSC02547.JPG"
    im = cv2.cvtColor(cv2.imread(path), cv2.COLOR_BGR2RGB)
    im = detect_edge(torch.from_numpy(im)[None])[0].cpu().numpy()
    cv2.imwrite("/mnt/d/temp.png", (255 * im / np.amax(im)).astype(np.uint8))

