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

