import torch
import torch.nn.functional as F

import math

from spirulae_splat.splat._camera import _Camera
from spirulae_splat.splat._torch_impl import depth_map, depth_inv_map
from spirulae_splat.splat import depth_to_normal


class SupervisionLosses(torch.nn.Module):

    _uv_cache = {}

    def __init__(self, config):
        super().__init__()

        self.depth_degree = config.depth_distortion_depth_degree
        self.uv_degree = config.depth_distortion_uv_degree

        self.depth_weight = config.depth_supervision_weight
        self.normal_weight = config.normal_supervision_weight
        self.alpha_weight = config.alpha_supervision_weight
        self.alpha_weight_under = config.alpha_supervision_weight_under

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

                # distortion disabled, just use it
                depth_diff = (x - y).reshape((h, w))
            # TODO: sometimes it's overfitting mean depth, need per-splat depth supervision

            # depth_loss = ((depth_diff**2)*alpha).mean()
            depth_loss = (torch.abs(depth_diff)*alpha).mean()

        # normal loss
        normal_loss = 0.0
        if self.normal_weight > 0.0:
            if pred_normal is None:
                raise ValueError("--pipeline.model.compute_depth_normal must be enabled to use normal supervision")
            ref_normal = depth_to_normal(ref_depth, camera)
            normal_loss = 1.0 - (ref_normal * pred_normal).sum(-1, True)
            # normal_loss = torch.sqrt(torch.relu(normal_loss))
            normal_loss = normal_loss.mean()

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
