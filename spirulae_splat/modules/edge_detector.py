import torch
import torch.nn.functional as F


@torch.no_grad()
def detect_edge_laplacian(im: torch.Tensor):
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
def detect_edge_laplacian_ms(im: torch.Tensor):
    if im.shape[-1] in [3, 4]:
        im = 0.299 * im[..., 0] + 0.587 * im[..., 1] + 0.114 * im[..., 2]
        im = im.unsqueeze(-1)
    edge = detect_edge_laplacian(im)
    sc = 1
    while im.shape[1] * im.shape[2] > 256 * 256:
        sc *= 2
        im = torch.nn.functional.avg_pool2d(im.permute(0, 3, 1, 2), 2).permute(0, 2, 3, 1)
        edge_i = detect_edge_laplacian(im)
        edge_i = edge_i[:, :, None, :, None, :] \
            .repeat(1, 1, sc, 1, sc, 1).reshape(edge_i.shape[0], sc*edge_i.shape[1], sc*edge_i.shape[2], edge_i.shape[3])
        edge[:, :edge_i.shape[1], :edge_i.shape[2], :] += edge_i
    return edge
