import torch

from typing import Optional


class ExposureCorrection(torch.nn.Module):

    def __init__(self, config):
        super().__init__()

        self.adaptive_exposure_mode = config.adaptive_exposure_mode

    def forward(
            self,
            pred_img: torch.Tensor, gt_img: torch.Tensor,
            alpha_mask: Optional[torch.Tensor]=None
        ):
        exposure_param_reg = 0.0
        eeps = 0.5/255.0

        total_alpha = max(alpha_mask.sum().item(), 1.0) if alpha_mask is not None else 1.0
        def immean(img: torch.Tensor, dim: Optional[int]=None, keepdim=False):
            if alpha_mask is None:
                return img.mean(dim, keepdim)
            mask = alpha_mask.reshape(*img.shape[:-1], 1)
            return (img*mask).sum(dim, keepdim) / total_alpha

        # gt ~ k * pred
        if self.adaptive_exposure_mode == "linear":
            k = (immean(gt_img)+eeps) / (immean(pred_img)+eeps)
            pred_img_e = k * (pred_img+eeps) - eeps
            # print(k.item())
            exposure_param_reg = (k-1.0)**2

        # log(gt) ~ k * log(pred) + b
        elif self.adaptive_exposure_mode == "log_linear":
            x = torch.log(pred_img+eeps)
            y = torch.log(gt_img+eeps)
            x_mean, y_mean = immean(x), immean(y)
            x2_mean, xy_mean = immean(x**2), immean(x*y)
            m = (xy_mean - x_mean * y_mean) / (x2_mean - x_mean**2)
            b = y_mean - m * x_mean
            pred_img_e = torch.exp(m * x + b) - eeps
            # print(m.item(), b.item())
            exposure_param_reg = (m-1.0)**2 + b**2

        # gt ~ k * pred, per channel
        elif self.adaptive_exposure_mode == "channel":
            x = pred_img.reshape(-1, 3)
            y = gt_img.reshape(-1, 3)
            k = (immean(y, 0, True)+eeps) / (immean(x, 0, True)+eeps)
            pred_img_e = k * (pred_img+eeps) - eeps
            # print(k.detach())
            exposure_param_reg = torch.mean((k-1.0)**2)

        # log(gt) ~ k * log(pred) + b, per channel
        elif self.adaptive_exposure_mode == "log_channel":
            x = torch.log(pred_img + eeps).reshape(-1, 3)
            y = torch.log(gt_img + eeps).reshape(-1, 3)
            x_mean, y_mean = immean(x, 0, True), immean(y, 0, True)
            x2_mean, xy_mean = immean(x**2, 0, True), immean(x*y, 0, True)
            m = (xy_mean - x_mean * y_mean) / (x2_mean - x_mean**2)
            b = y_mean - m * x_mean
            pred_img_e = torch.exp(m * x + b).reshape(pred_img.shape) - eeps
            # print(m.detach(), b.detach())
            exposure_param_reg = torch.mean((m-1.0)**2 + b**2)

        # gt ~ A * pred
        elif self.adaptive_exposure_mode == "affine":
            x = (pred_img + eeps).reshape(-1, 3)
            y = (gt_img + eeps).reshape(-1, 3)
            if alpha_mask is not None:
                x = x * alpha_mask.reshape(-1, 1)
                y = y * alpha_mask.reshape(-1, 1)
            xTy, xTx = x.T @ y, x.T @ x
            xTy, xTx = xTy.cpu(), xTx.cpu()
            try:
                xTx = xTx + (1e-6*torch.trace(xTx).item()) * torch.eye(3)
                A = torch.linalg.inv(xTx) @ xTy
            except torch._C._LinAlgError:
                print("Warning: torch._C._LinAlgError")
                print(xTx)
                pred_img_e = pred_img.detach()
            else:
                A = A.to(x.device)
                pred_img_e = (x @ A).reshape(pred_img.shape) - eeps
                exposure_param_reg = torch.mean((A-torch.eye(3).to(A))**2)

        # log(gt) ~ A * log(pred) + b
        elif self.adaptive_exposure_mode == "log_affine":
            x = torch.log(pred_img + eeps).reshape(-1, 3)
            y = torch.log(gt_img + eeps).reshape(-1, 3)
            if alpha_mask is not None:
                x = x * alpha_mask.reshape(-1, 1)
                y = y * alpha_mask.reshape(-1, 1)
            mean_x = immean(x, 0, True)
            mean_y = immean(y, 0, True)
            centered_x = x - mean_x
            centered_y = y - mean_y
            xTy = centered_x.T @ centered_y
            xTx = centered_x.T @ centered_x
            xTy, xTx = xTy.cpu(), xTx.cpu()
            try:
                xTx = xTx + (1e-6*torch.trace(xTx).item()) * torch.eye(3)
                A = torch.linalg.inv(xTx) @ xTy
            except torch._C._LinAlgError:
                print("Warning: torch._C._LinAlgError")
                print(xTx)
                pred_img_e = pred_img.detach()
            else:
                A = A.to(x.device)
                B = mean_y - mean_x @ A
                pred_img_e = torch.exp(x @ A + B).reshape(pred_img.shape) - eeps
                exposure_param_reg = 0.5*(torch.mean((A-torch.eye(3).to(A))**2) + torch.mean(B**2))

        else:
            raise ValueError("Invalid adaptive_exposure_mode:", self.adaptive_exposure_mode)
        
        if alpha_mask is not None:
            pred_img_e = pred_img + (pred_img_e-pred_img) * alpha_mask
        return pred_img_e, exposure_param_reg
