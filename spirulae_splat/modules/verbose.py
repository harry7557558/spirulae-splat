import torch
from typing import Dict, Optional, List, Union
import math
from dataclasses import dataclass

from time import perf_counter


@dataclass
class Metric:
    decimals: Optional[int] = None
    sigfigs: Optional[int] = None


class TrainingVerbose:

    def __init__(self):
        self.sma_length = 100
        self.cache_skip = 10
        self.metrics = {}  # type: Dict[str, List[float]]
        self.cache = None
        self.num_new_lines = 0

    def add_metric(self, key: str, value: float, last_only: bool = False):
        if key not in self.metrics:
            self.metrics[key] = []
        self.metrics[key].append(value)
        if len(self.metrics[key]) > (1 if last_only else self.sma_length):
            del self.metrics[key][0]

    def get_metric(self, key: str, d: int = None, s: int = None) -> Union[str, float]:
        if key not in self.metrics or len(self.metrics[key]) == 0:
            if d is not None and s is not None:
                return float('nan')
            return '~'
        v = sum(self.metrics[key]) / len(self.metrics[key])
        if d is not None and s is not None:
            return v
        if v == 0.0:
            return '~'
        if not math.isfinite(v):
            return str(v)
        if d is not None:
            return f"{{:.{d}f}}".format(v)
        if s is not None:
            return f"{{:.{s}g}}".format(v)
        return str(v)

    def get_last_metric(self, key: str, d: int = None, s: int = None) -> str:
        if key not in self.metrics or len(self.metrics[key]) == 0:
            return '~'
        v = self.metrics[key][-1]
        if not math.isfinite(v):
            return str(v)
        if d is not None:
            return f"{{:.{d}f}}".format(v)
        if s is not None:
            return f"{{:.{s}g}}".format(v)
        return str(v)

    def verbose(self, step: int, max_steps: int):

        # Time
        cur_time = perf_counter()
        if not hasattr(self, 'start_time'):
            self.start_time = cur_time
            self.prev_time = cur_time
        else:
            self.add_metric('iter_time', cur_time-self.prev_time)
            self.prev_time = cur_time

        # VRAM
        vram_free, vram_total = [x / 1024**3 for x in torch.cuda.mem_get_info(None)]
        vram_used = vram_total - vram_free
        vram_used_percentage = (1 - vram_free / vram_total)*100
        vram_allocated = torch.cuda.memory_allocated() / 1024**3
        vram_reserved = torch.cuda.memory_reserved() / 1024**3
        # mem_stats = boldcyan(f"{used:.2f}") + f"\N{ZERO WIDTH SPACE}GiB " + boldcyan(f"{used_percentage:.0f}") + "%"
        # return f"{redbkg(bracket('Mem'), used_percentage/90)} {mem_stats}"

        if step % self.cache_skip != 0 and self.cache is not None:
            return

        lines = []

        def bracket(s: str):
            return f"\033[1m[\033[m{s}\033[1m]\033[m"

        def redbkg(s: str, threshold: float = 1.0):
            return f"\033[1;41m{s.replace('\033[m', '')}\033[m" if threshold >= 1.0 else s

        def orange(s: str):
            return f"\033[33m{s}\033[m"

        def boldcyan(s: str):
            if s == '~':
                return s
            return f"\033[1;36m{s}\033[m"

        iter_time = self.get_metric('iter_time', d=0, s=0)
        # lines.append('Training in progress' + '.'*(1+(step//self.cache_skip)%3))
        lines.append('Training in progress...')
        lines.append('-' * 80)
        lines.append(
            f"{bracket('Time')} step {boldcyan(step)}/{max_steps} | "
            f"{boldcyan(f"{1e3*iter_time:.1f}")} {orange('ms')}/step, "
            f"{boldcyan(f"{cur_time-self.start_time:.1f}")}{orange('s')} elapsed, "
            f"{boldcyan(f"{(max_steps-step)*iter_time:.1f}")}{orange('s')} remaining "
        )
        lines.append(
            f"{redbkg(bracket('Memory'), vram_used_percentage/90)} "
            f"{boldcyan(f"{vram_used:.2f}")}/{vram_total:.2f} {orange('GiB')} "
            f"({boldcyan(f"{vram_used_percentage:.0f}")}{orange('%')}) | "
            f"{boldcyan(f"{vram_allocated:.2f}")} {orange('GiB')} allocated, "
            f"{boldcyan(f"{vram_reserved:.2f}")} {orange('GiB')} reserved "
        )
        lines.append(
            f"    ("
            f"{boldcyan(self.get_last_metric('splat_vram', d=2))} {orange('splat')}, "
            f"{boldcyan(self.get_last_metric('image_vram', d=2))} {orange('image')}, "
            f"{boldcyan(self.get_last_metric('splat_x_image_vram', d=2))} {orange('splat x img')}, "
            f"{boldcyan(self.get_last_metric('appearance_vram', d=2))} {orange('appearance')}, "
            f"{boldcyan(self.get_last_metric('viewer_vram', d=2))} {orange('viewer')}, "
            f"{boldcyan(self.get_last_metric('other_vram', d=2))} {orange('other')}"
            f")"
        )
        lines.append('-' * 80)
        lines.append(
            f"{bracket('Splat')} {boldcyan(self.get_last_metric('num_splats'))}/{self.get_last_metric('max_num_splats')} "
            f"{orange('SH')}={boldcyan(self.get_last_metric('num_sh'))}/{self.get_last_metric('max_num_sh')} | "
            f"{orange('opac')}={boldcyan(self.get_last_metric('opac', d=3))} "
            f"{orange('scale')}={boldcyan(self.get_last_metric('scale', s=3))} "
            f"{orange('aniso')}={boldcyan(self.get_last_metric('aniso', s=3))} "
            f"{orange('erank')}={boldcyan(self.get_last_metric('erank', s=3))} "
        )
        lines.append(
            f"{bracket('Image')} {orange('loss')}={boldcyan(self.get_metric('rgb_loss', s=3))} "
            f"{orange('psnr')}={boldcyan(self.get_metric('psnr', d=2))} "
            f"{orange('ssim')}={boldcyan(self.get_metric('ssim', d=3))} | "
            f"{orange('depth')}={boldcyan(self.get_metric('depth_loss', s=3))} "
            f"{orange('normal')}={boldcyan(self.get_metric('normal_loss', s=3))} "
            f"{orange('alpha')}={boldcyan(self.get_metric('alpha_loss', s=3))} "
        )
        lines.append(
            f"{bracket('Bilagrid')} "
            f"{orange('tv')}={boldcyan(self.get_last_metric('bilagrid_tv', s=3))} "
            f"{orange('tv.depth')}={boldcyan(self.get_last_metric('bilagrid_depth_tv', s=3))} "
            f"{orange('tv.normal')}={boldcyan(self.get_last_metric('bilagrid_normal_tv', s=3))} "
        )
        lines.append(
            f"{bracket('PPISP')} "
            f"{orange('eμ')}={boldcyan(self.get_last_metric('ppisp_reg_exposure_mean', s=3))} "
            f"{orange('vc')}={boldcyan(self.get_last_metric('ppisp_reg_vig_center', s=3))} "
            f"{orange('v+')}={boldcyan(self.get_last_metric('ppisp_reg_vig_non_pos', s=3))} "
            f"{orange('vσ')}={boldcyan(self.get_last_metric('ppisp_reg_vig_channel_var', s=3))} "
            f"{orange('cμ')}={boldcyan(self.get_last_metric('ppisp_reg_color_mean', s=3))} "
            f"{orange('rσ')}={boldcyan(self.get_last_metric('ppisp_reg_crf_channel_var', s=3))} "
        )
        lines.append(
            f"{bracket('HDR')} {orange('WIP')}  "
            f"{bracket('Validation')} {orange('WIP')}  "
        )
        lines.append('-' * 80)
        lines.append("")

        if self.num_new_lines > 0:
            print("\033[F"*self.num_new_lines, end='')
        else:
            print(end='\n')
        self.num_new_lines = len(lines)

        self.cache = '    \n'.join(lines)
        print(self.cache, end='\n')



        #     f"{bracket('Train')} {orange('loss')}={fmt('image_loss', 1.0)} "
        #     f"{orange('psnr')}={fmt('psnr', 1.0, 2)} "
        #     f"{orange('ssim')}={fmt('ssim', 1.0, 3)}" + \
        #         f" lpips={fmt('lpips', 1.0, 3)}" * ('lpips' in losses),
        #     "                \n",

        # if hasattr(self, 'total_val_loss'):
        #     losses['val_total'] = self.total_val_loss
        # if hasattr(self, 'total_train_loss'):
        #     losses['train_total'] = self.total_train_loss
        # if hasattr(self, 'total_val_lpips'):
        #     losses['lpips_val'] = self.total_val_lpips
        # if hasattr(self, 'overfit_score'):
        #     losses['overfit_score'] = self.overfit_score

        # def fmt(key: str, s: float, decimals=None, sigfigs=3) -> str:
        #     if s == 0.0 or key not in losses:
        #         return '~'

        #     l = losses[key]
        #     if isinstance(l, torch.Tensor):
        #         l = l.detach().item()
        #     l = l / s
        #     if not math.isfinite(l):  # not finite
        #         return str(l)

        #     if key not in _storage or self.step % 1000 == 0:
        #         _storage[key] = abs(l)
        #     _storage[key] = max(_storage[key], abs(l))
        #     if _storage[key] == 0.0:
        #         return '~'

        #     if (abs(l) < 0.1**(sigfigs+1) or abs(l) >= 1.0-(0.1**sigfigs)) and decimals is None:
        #         s = f"{{:.{sigfigs}g}}".format(l)
        #     else:
        #         if decimals is None:  # 3 sig figs
        #             decimals = int(max(-math.log10((0.1**sigfigs)*abs(l)), 0)) if l != 0.0 else 0
        #         s = f"{{:.{decimals}f}}".format(l)
        #     if s.startswith('0.'):
        #         s = s[1:]
        #     if s.startswith('-0.'):
        #         s = '-' + s[2:]
        #     return boldcyan(s)

        # reg_mcmc = (self.config.use_mcmc and self.step < self.trainer_config.num_iterations-self.config.refine_stop_num_iter)
        # reg_2dgs = self.training_losses.get_2dgs_reg_weights()
        # opacity_floor = (
        #     f"{bracket('OpacFloor')} {self.strategy.get_opacity_floor(self.step):.3f} "
        #     f"{bracket('Hardness')} {self.strategy.get_hardness(self.step):.3f}"
        # ).replace('0.', '.') \
        #     if self.config.primitive == "opaque_triangle" else ""
        # chunks = [
        #     f"{boldcyan(self.step)} "
        #     f"{bracket('N')} {boldcyan(self.core.cur_num_splats)}",
        #     f"{redbkg(bracket('Mem'), used_percentage/90)} {mem_stats}",
        #     f"{bracket('Train')} {orange('loss')}={fmt('image_loss', 1.0)} "
        #     f"{orange('psnr')}={fmt('psnr', 1.0, 2)} "
        #     f"{orange('ssim')}={fmt('ssim', 1.0, 3)}" + \
        #         f" lpips={fmt('lpips', 1.0, 3)}" * ('lpips' in losses),
        #     "                \n",
        #     f"{bracket('RefLoss')} {orange('depth')}={fmt('depth_ref_loss', self.config.depth_supervision_weight, 3)} "
        #     f"{orange('normal')}={fmt('normal_ref_loss', self.config.normal_supervision_weight, 3)} "
        #     f"{orange('alpha')}={fmt('alpha_ref_loss', 0.5*(self.config.alpha_loss_weight+self.config.alpha_loss_weight_under), 4)}",
        #     f"{bracket('DistLoss')} {orange('depth')}={fmt('depth_dist_reg', reg_2dgs[0][0], 3)} "
        #     f"{orange('normal')}={fmt('normal_dist_reg', reg_2dgs[0][1], 3)} "
        #     f"{orange('rgb')}={fmt('rgb_dist_reg', reg_2dgs[0][2], 3)}",
        #     "                \n",
        #     f"{bracket('ImReg')} {orange('normal')}={fmt('normal_reg', reg_2dgs[1], 3)} "
        #     f"{orange('alpha')}={fmt('alpha_reg', self.training_losses.get_alpha_reg_weight(), 3)}",
        #     f"{bracket('SplatReg')} {orange('o')}={fmt('opacity_reg', self.config.opacity_reg * reg_mcmc, 3)} "
        #     f"{orange('s')}={fmt('scale_reg', self.config.scale_reg * reg_mcmc, 4)} "
        #     f"{orange('q')}={fmt('quat_norm_reg', self.config.quat_norm_reg, sigfigs=1)} "
        #     f"{orange('erank')}={fmt('erank_reg', max(self.config.erank_reg_s3, self.config.erank_reg))} "
        #     f"{orange('aniso')}={fmt('scale_reg', self.config.scale_regularization_weight)}",
        #     "                \n",
        #     f"{bracket('Bilagrid')} {orange('reg')}={fmt('bilagrid_mean_reg', self.config.bilagrid_mean_reg_weight)} "
        #     f"{orange('tv.rgb')}={fmt('tv_loss', self.config.bilagrid_tv_loss_weight)} "
        #     f"{orange('depth')}={fmt('tv_loss_depth', self.config.bilagrid_tv_loss_weight_geometry)} "
        #     f"{orange('normal')}={fmt('tv_loss_normal', self.config.bilagrid_tv_loss_weight_geometry)}"
        #     "                \n",
        #     f"{bracket('PPISP')} {orange('eμ')}={fmt('ppisp_reg_exposure_mean', self.config.ppisp_reg_exposure_mean)} "
        #     f"{orange('vc')}={fmt('ppisp_reg_vig_center', self.config.ppisp_reg_vig_center)} "
        #     f"{orange('v+')}={fmt('ppisp_reg_vig_non_pos', self.config.ppisp_reg_vig_non_pos)} "
        #     f"{orange('vσ')}={fmt('ppisp_reg_vig_channel_var', self.config.ppisp_reg_vig_channel_var)} "
        #     f"{orange('cμ')}={fmt('ppisp_reg_color_mean', self.config.ppisp_reg_color_mean)}",
        #     f"{orange('rσ')}={fmt('ppisp_reg_crf_channel_var', self.config.ppisp_reg_crf_channel_var)}",
        #     "                \n",
        #     f"{bracket('Validation')} {orange('train')}={fmt('train_total', 1.0)} "
        #     f"{orange('val')}={fmt('val_total', 1.0)} "
        #     f"{orange('lpips')}={fmt('lpips_val', 1.0)} "
        #     f"{orange('overfit')}={fmt('overfit_score', 1.0)} "
        #     f"{redbkg(f'{self.overfit_count}/{self.config.early_stop_patience}',
        #         self.overfit_count/(0.8*self.config.early_stop_patience))}"
        # ]
        # additional_chunks = []
        # if len(opacity_floor) > 0:
        #     additional_chunks.append(opacity_floor)
        # if 'bvh_time' in self.info:
        #     losses['bvh_time'] = self.info['bvh_time']
        #     additional_chunks.append(f"{bracket('BvhTime')} {fmt('bvh_time', 1.0, 1)} ms")
        # if len(additional_chunks) > 0:
        #     chunks.append("                \n")
        #     chunks.extend(additional_chunks)
        # chunks.append("                \n")
        # chunks = ' '.join(chunks).replace('\n ', '\n')
        # chunks += "\033[F"*(chunks.count('\n'))
        # print(chunks, end='')
        # # _storage['print_queue'].append(chunks)
        # _storage['chunks'] = chunks
