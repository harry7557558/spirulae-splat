import glob
import json
import os
import shutil
from subprocess import DEVNULL, call

from rich.console import Console
from torch.utils.cpp_extension import _get_build_directory, load

PATH = os.path.dirname(os.path.abspath(__file__))


from spirulae_splat import csrc as _C

__all__ = ["_C"]
