import os
import re
import glob
import sys
import platform
from pathlib import Path
from setuptools import setup, find_namespace_packages


BUILD_NO_CUDA = os.getenv("BUILD_NO_CUDA", "0") == "1"
WITH_SYMBOLS = os.getenv("WITH_SYMBOLS", "0") == "1"
# LINE_INFO = os.getenv("LINE_INFO", "0") == "1"
LINE_INFO = True


def get_ext():
    from torch.utils.cpp_extension import BuildExtension
    import multiprocessing

    class CustomBuildExtention(BuildExtension):
        def build_extensions(self):
            self.parallel = multiprocessing.cpu_count()

            # Prevent out of memory on low memory devices
            try:
                import psutil
                free_mem = psutil.virtual_memory().available
                mem_per_worker = 2 * 1024**3  # 2GB; TODO: separate libtorch like fused-bilagrid
                self.parallel = min(
                    self.parallel,
                    max(free_mem // mem_per_worker, 1)
                )
            except ImportError:
                pass

            super().build_extensions()

    return CustomBuildExtention.with_options(no_python_abi_suffix=True, use_ninja=True)


def get_extensions():
    import torch
    from torch.utils.cpp_extension import CUDAExtension

    cuda_arch_list = []
    if torch.cuda.is_available():
        for i in range(torch.cuda.device_count()):
            compute_capability = torch.cuda.get_device_capability(i)
            cuda_arch_list.append(f"{compute_capability[0]}{compute_capability[1]}")
    else:
        raise RuntimeError("CUDA is required for this extension.")

    extensions_dir = Path("spirulae_splat/splat/cuda")
    sources = (
        glob.glob(str(extensions_dir / "ins" / "*.cu")) +
        glob.glob(str(extensions_dir / "csrc" / "*.cu")) +
        glob.glob(str(extensions_dir / "csrc" / "*.cpp"))
    )
    sources = [s for s in sources if "hip" not in s]

    undef_macros = []
    define_macros = []

    if sys.platform == "win32":
        define_macros += [("spirulae_splat_EXPORTS", None)]

    extra_compile_args = {"cxx": ["-O3"]}
    if not os.name == "nt":  # Not on Windows:
        extra_compile_args["cxx"] += ["-Wno-sign-compare"]
    extra_link_args = [] if WITH_SYMBOLS else ["-s"]

    info = torch.__config__.parallel_info()
    if (
        "backend: OpenMP" in info
        and "OpenMP not found" not in info
        and sys.platform != "darwin"
    ):
        extra_compile_args["cxx"] += ["-DAT_PARALLEL_OPENMP"]
        if sys.platform == "win32":
            extra_compile_args["cxx"] += ["/openmp"]
        else:
            extra_compile_args["cxx"] += ["-fopenmp"]
    else:
        print("Compiling without OpenMP...")

    # Compile for mac arm64
    if sys.platform == "darwin" and platform.machine() == "arm64":
        extra_compile_args["cxx"] += ["-arch", "arm64"]
        extra_link_args += ["-arch", "arm64"]

    nvcc_flags = os.getenv("NVCC_FLAGS", "")
    nvcc_flags = [] if nvcc_flags == "" else nvcc_flags.split(" ")
    nvcc_flags += ["-O3", "--use_fast_math"]
    # nvcc_flags += ["-rdc=true", "-dlto"]
    # nvcc_flags += ["--extra-device-vectorization"]
    if LINE_INFO:
        nvcc_flags += ["-lineinfo", "--generate-line-info", "--source-in-ptx"]
        nvcc_flags += [
            # "-Xptxas", "-v",
            "-Xptxas", "--warn-on-double-precision-use",
            # "-Xptxas", "--warn-on-local-memory-usage",
            # "-Xptxas", "--warn-on-spills"
        ]
    if torch.version.hip:
        # USE_ROCM was added to later versions of PyTorch.
        # Define here to support older PyTorch versions as well:
        define_macros += [("USE_ROCM", None)]
        undef_macros += ["__HIP_NO_HALF_CONVERSIONS__"]
    else:
        nvcc_flags += ["--expt-relaxed-constexpr"]
    extra_compile_args["nvcc"] = nvcc_flags
    if sys.platform == "win32":
        extra_compile_args["nvcc"] += ["-DWIN32_LEAN_AND_MEAN"]

    # extra_compile_args["nvcc"] += ['-rdc=true']

    # disable compile warnings for glm
    extra_compile_args["nvcc"] += ['-Xcudafe=--diag_suppress=20012']

    # disable compile warnings for Slang generated code
    extra_compile_args["nvcc"] += ['-Xcudafe=--diag_suppress=550']

    extra_compile_args["nvcc"] += ['--threads', '0']

    # extra_compile_args["nvcc"] += ['-arch=native']
    for arch in cuda_arch_list:
        extra_compile_args["nvcc"] += [
            "-gencode", f"arch=compute_{arch},code=sm_{arch}"
        ]

    # enable host-side SIMD, etc.
    if sys.platform != "win32":
        extra_compile_args["nvcc"] += ['-Xcompiler', '-O3 -march=native']

    extension = CUDAExtension(
        "spirulae_splat.csrc",
        sources,
        include_dirs=[
            str((extensions_dir / "csrc").absolute()),
            str((extensions_dir / "csrc" / "glm").absolute()),
        ],
        define_macros=define_macros,
        undef_macros=undef_macros,
        extra_compile_args=extra_compile_args,
        extra_link_args=extra_link_args,
        # dlink=True
    )

    return [extension]


import importlib.util
if importlib.util.find_spec('torch') is None:
    raise ValueError("Please make sure you have PyTorch installed.")
no_fused_bilagrid = (importlib.util.find_spec('fused_bilagrid') is None)


setup(
    name="spirulae_splat",
    description="TODO",
    version="0.1.0",
    packages=find_namespace_packages(include=["spirulae_splat*"]),
    install_requires=[
        # "torch",
        "jaxtyping",
        "tyro",
        "opencv-python",
        "plyfile",
        "open3d",
        "matplotlib",
        "Pillow",
        "rawpy",
        "pytorch-msssim",
        "typing_extensions",
        "tabulate",
    ] + [
        "fused_bilagrid @ git+https://github.com/harry7557558/fused-bilagrid.git@dev",
    ] * no_fused_bilagrid,
    extras_require={
    },
    entry_points={
    },
    cmdclass={"build_ext": get_ext()} if not BUILD_NO_CUDA else {},
    ext_modules=get_extensions() if not BUILD_NO_CUDA else [],
)
