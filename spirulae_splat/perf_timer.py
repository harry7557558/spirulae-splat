import time
import math
import torch
import numpy as np
from collections import defaultdict


class PerfTimer:
    _DISABLED = True  # Static variable to globally disable the timer
    _SYNC_INTERMEDIATE = True  # Static variable to control intermediate CUDA synchronization

    def __init__(self, name: str, ema_tau=10):
        self.name = name
        self.ema_tau = ema_tau
        self.timers = defaultdict(list)
        self.ema = None
        self.marks = []
        self.times = []

    @staticmethod
    def _get_time():
        if PerfTimer._SYNC_INTERMEDIATE:
            torch.cuda.synchronize()
        return int(1e6 * time.perf_counter() + 0.5)

    def start(self):
        if PerfTimer._DISABLED:
            return
        self.times = [self._get_time()]
        self.marks = []

    def mark(self, task):
        if PerfTimer._DISABLED:
            return
        self.times.append(self._get_time())
        self.marks.append(task[:1])

    def end(self, task):
        if PerfTimer._DISABLED:
            return
        self.mark(task)
        self.times = np.array(self.times)
        self.times = np.concatenate((
            self.times[1:]-self.times[:-1],
            self.times[-1:]-self.times[:1]
        ))
        if self.ema is None:
            self.ema = self.times
        alpha = 1 - math.exp(-1 / self.ema_tau)
        if len(self.ema) == len(self.times):
            self.ema = alpha * self.times + (1 - alpha) * self.ema
        self.summary()

    def summary(self):
        if PerfTimer._DISABLED:
            return
        durations = (self.ema+0.5).astype(np.int32)
        subtask_times = " ".join([c+str(d) for c, d in zip(self.marks, durations[:-1])])
        print(f"{self.name}: {subtask_times} -> {durations[-1]}")

    @classmethod
    def disable(cls):
        cls._DISABLED = True

    @classmethod
    def enable(cls):
        cls._DISABLED = False

    @classmethod
    def set_sync_intermediate(cls, value):
        cls._SYNC_INTERMEDIATE = value
