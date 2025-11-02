import torch
import torch.nn.functional as F

import math

from spirulae_splat.splat._camera import _Camera
from spirulae_splat.splat._torch_impl import depth_map, depth_inv_map
from spirulae_splat.splat import depth_to_normal


def normalize_image(im, degree=1, _B_cache={}):
    H, W, C = im.shape

    key = (H, W, degree)
    if key in _B_cache:
        B, invBTB = _B_cache[key]

    else:

        v, u = torch.meshgrid(
            torch.pi * (torch.arange(H, device=im.device)+0.5)/H,
            torch.pi * (torch.arange(W, device=im.device)+0.5)/W,
            indexing='ij'
        )
        v, u = v.flatten(), u.flatten()

        B = []
        for i in range(degree+1):
            for j in range(degree+1):
                ci, si = torch.cos(i*u), torch.sin(i*u)
                cj, sj = torch.cos(j*v), torch.sin(j*v)
                B.extend(
                    [ci * cj] +
                    [ci * sj] * (j != 0) +
                    [si * cj] * (i != 0) +
                    [si * sj] * (j != 0 and i != 0)
                )

        B = torch.stack(B, dim=-1)
        invBTB = torch.linalg.inv(B.H @ B)

        _B_cache[key] = (B, invBTB)

    X = im.reshape((H*W, C))
    weights = invBTB @ B.H @ X
    Y = X - B @ weights

    im = Y.reshape(im.shape)
    # return im / torch.norm(im, dim=-1, keepdim=True)
    return im / (torch.norm(im, dim=(0,1), keepdim=True) / (im.shape[0]*im.shape[1])**0.5)


class SupervisionLosses(torch.nn.Module):

    _uv_cache = {}

    def __init__(self, config):
        super().__init__()

        self.depth_degree = config.depth_distortion_depth_degree
        self.uv_degree = config.depth_distortion_uv_degree

        self.depth_weight = config.depth_supervision_weight
        self.normal_weight = config.normal_supervision_weight
        self.alpha_weight = config.alpha_loss_weight
        self.alpha_weight_under = config.alpha_loss_weight_under

    @staticmethod
    def get_alpha_loss(x, y):
        """Compute asymmetric loss for alpha, penalize only when first argument is greater than second argument

        Args:
            x: alpha output by the model
            y: reference alpha with same shape as x
        """
        # return torch.relu(x-y).mean()
        # return F.binary_cross_entropy(x, y, reduction="mean")
        return F.binary_cross_entropy(torch.fmax(x, y), y, reduction="mean")

    def forward(self, ref_depth, camera: _Camera, pred_depth, pred_normal, pred_alpha):
        # details: https://github.com/harry7557558/Graphics/blob/master/mapping/relative_depth_matching/depth_fitting_01.ipynb

        ref_depth = torch.nan_to_num(ref_depth, 0.0, 0.0, 0.0)
        ref_depth = torch.clip(ref_depth, 0.0, 1e4)

        # resize depth
        h, w, _ = pred_depth.shape
        alpha = (ref_depth > 0.0)
        ref_depth = torch.where(alpha, ref_depth, 2.0*torch.amax(ref_depth).detach())
        alpha = alpha.squeeze(-1)

        # generate UV embeddings
        if not self.uv_degree >= 0:
            uv = None

        elif (h, w) not in self._uv_cache:
            u = (torch.arange(w, dtype=torch.float32)+0.5)/w * 2.0 - 1.0
            v = (torch.arange(h, dtype=torch.float32)+0.5)/h * 2.0 - 1.0
            u = u[None, :].float().to(pred_depth.device).repeat((h, 1))
            v = v[:, None].float().to(pred_depth.device).repeat((1, w))
            u, v = u.flatten(), v.flatten()

            uv = []
            for i in range(self.uv_degree+1):
                for j in range(self.uv_degree+1):
                    uv.append(torch.cos(math.pi/2*i*u)*torch.cos(math.pi/2*j*v))
                    if j != 0:
                        uv.append(torch.cos(math.pi/2*i*u)*torch.sin(math.pi/2*j*v))
                    if i != 0:
                        uv.append(torch.sin(math.pi/2*i*u)*torch.cos(math.pi/2*j*v))
                    if i != 0 and j != 0:
                        uv.append(torch.sin(math.pi/2*i*u)*torch.sin(math.pi/2*j*v))

            uv = torch.stack(uv)
            self._uv_cache[(h, w)] = uv

        else:
            uv = self._uv_cache[(h, w)]

        # depth distortion fitting
        depth_loss = 0.0
        if self.depth_weight > 0.0:
            # normalize depths - inputs are metric linear depth
            x, y = ref_depth.flatten(), pred_depth.flatten()
            # x, y = pred_depth.flatten(), ref_depth.flatten()

            m = alpha.flatten().float()
            m_sum = alpha.sum().item()

            if (
                self.depth_degree >= 0 and
                self.uv_degree >= 0
            ):
                # normalize by mean and log
                x = x / ((x*m).sum() / m_sum)
                y = y / ((y*m).sum() / m_sum)
                x, y = depth_map(x), depth_map(y)

                # generate depth embeddings
                A = []
                for k in range(self.depth_degree+1):
                    ed = x**k / math.factorial(k) * torch.exp(-x)
                    A.append(ed * uv)
                A = torch.concatenate(A)

                # linear least squares
                mat = (alpha.reshape(1,-1) * A) @ A.T
                vec = A @ (alpha.flatten() * y)
                c = torch.linalg.solve(mat, vec)
                corr_depth = torch.relu(c @ A)
                depth_diff = (corr_depth - y).reshape((h, w))

            else:
                # normalize log by mean and std
                x, y = torch.log(torch.relu(x) + 1e-6), torch.log(torch.relu(y) + 1e-6)
                x = x - (x*m).sum() / m_sum
                y = y - (y*m).sum() / m_sum
                x = x / torch.sqrt(((x*x)*m).sum() / m_sum)
                y = y / torch.sqrt(((y*y)*m).sum() / m_sum)

                depth_diff = (x - y)**2
                # depth_diff = torch.abs(x - y)
                # depth_diff = (normalize_image(x.reshape(ref_depth.shape)) - normalize_image(y.reshape(ref_depth.shape)))**2

            # TODO: sometimes it's overfitting mean depth, need per-splat depth supervision

            depth_loss = (depth_diff.reshape(alpha.shape) * alpha).mean()

        # normal loss
        normal_loss = 0.0
        if self.normal_weight > 0.0:
            if camera.model == "OPENCV_FISHEYE":
                # disable distortion since most monocular depth models aren't trained on fisheye
                # gives much better results than accurate fisheye-distorted normal
                camera.model = ""
                pred_normal = depth_to_normal(pred_depth, camera)
                ref_normal = depth_to_normal(ref_depth, camera)
                camera.model = "OPENCV_FISHEYE"
            else:
                if pred_normal is None:
                    raise ValueError("--pipeline.model.compute_depth_normal must be enabled to use normal supervision")
                ref_normal = depth_to_normal(ref_depth, camera)
            ref_normal = torch.nan_to_num(ref_normal, 0.0, 0.0, 0.0)
            pred_normal = torch.nan_to_num(pred_normal, 0.0, 0.0, 0.0)

            # from torchvision.transforms.functional import gaussian_blur
            # s = ref_normal.shape[0]
            # ref_normal = gaussian_blur(ref_normal.permute(2,0,1), s//30+1, s/100).permute(1,2,0)
            # pred_normal = gaussian_blur(pred_normal.permute(2,0,1), s//30+1, s/100).permute(1,2,0)

            normal_loss = 1.0 - (ref_normal * pred_normal).sum(-1, True)
            # normal_loss = (ref_normal - pred_normal)**2
            # normal_loss = (normalize_image(ref_normal) - normalize_image(pred_normal))**2
            # normal_loss = torch.abs(normalize_image(ref_normal) - normalize_image(pred_normal))
            normal_loss = (normal_loss * alpha.reshape(*normal_loss.shape[:-1], 1)).mean()

        # alpha loss
        pred_alpha, alpha = pred_alpha.float().flatten(), alpha.float().flatten()
        alpha_loss = self.get_alpha_loss(pred_alpha, alpha)
        alpha_loss_over = 0.0
        if self.alpha_weight_under > 0.0:
            alpha_loss_over = self.get_alpha_loss(1.0-pred_alpha, 1.0-alpha)

        return (
            self.depth_weight * depth_loss,
            self.normal_weight * normal_loss,
            self.alpha_weight * alpha_loss +
                self.alpha_weight_under * alpha_loss_over
        )
