import glob
import os
import platform
import sys
import re

from setuptools import find_packages, setup


BUILD_NO_CUDA = os.getenv("BUILD_NO_CUDA", "0") == "1"
WITH_SYMBOLS = os.getenv("WITH_SYMBOLS", "0") == "1"
# LINE_INFO = os.getenv("LINE_INFO", "0") == "1"
LINE_INFO = True


def extract_function_declarations(code):
    # Regex to match non-inline function declarations
    function_decl_pattern = re.compile(r"""
        # Match comments before the function declaration
        (?:/\*[^*]*\*+(?:[^/*][^*]*\*+)*/|//[^\n]*?$\s*)*
        # Match return type (including template declarations)
        (?:template\s*<[^>]+>\s*)?
        # Match function attributes like __global__, __device__, etc.
        (?:inline)?\s*
        (?:__global__|__device__)?\s*
        # Match the return type
        (?:void|u?int[234]?|float[234]?|torch::Tensor|std::tuple<[\w:\s*&<>\[\],\/]+?>)
        # Match the function name
        \s+\b\w+\b\s*
        # Match the function parameters
        \(.*?\)\s
    """, re.MULTILINE | re.VERBOSE | re.DOTALL)
    
    matches = function_decl_pattern.findall(code)
    decls = []

    for m in matches:
        if re.compile(r"inline\s+(__global__|__device__)").findall(m):
            continue
        if True and re.compile(r"(__global__|__device__)").findall(m):
            continue
        decls.append(m.strip()+';')
    
    return decls


def generate_header(source_filename, header_filename):
    code = open(source_filename).read()
    decls = extract_function_declarations(code)

    splitter = "/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */\n"
    include = open(header_filename).read()
    include = include.split(splitter)[0].strip()

    header = '\n\n\n'.join([include, splitter]+decls)
    open(header_filename, "w").write(header+'\n')


def get_ext():
    from torch.utils.cpp_extension import BuildExtension
    import multiprocessing

    class CustomBuildExtention(BuildExtension):
        def build_extensions(self):
            self.parallel = multiprocessing.cpu_count()
            super().build_extensions()

    return CustomBuildExtention.with_options(no_python_abi_suffix=True, use_ninja=True)


def get_extensions():
    import torch
    from torch.__config__ import parallel_info
    from torch.utils.cpp_extension import CUDAExtension

    extensions_dir = os.path.abspath(os.path.join("spirulae_splat", "splat", "cuda", "csrc"))
    sources = glob.glob(os.path.join(extensions_dir, "*.cu")) + \
        glob.glob(os.path.join(extensions_dir, "*.cpp"))
    sources = [path for path in sources if "hip" not in path]

    undef_macros = []
    define_macros = []

    if sys.platform == "win32":
        define_macros += [("spirulae_splat_EXPORTS", None)]

    extra_compile_args = {"cxx": ["-O3"]}
    if not os.name == "nt":  # Not on Windows:
        extra_compile_args["cxx"] += ["-Wno-sign-compare"]
    extra_link_args = [] if WITH_SYMBOLS else ["-s"]

    info = parallel_info()
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
    if LINE_INFO:
        nvcc_flags += ["-lineinfo", "--generate-line-info", "--source-in-ptx"]
        # nvcc_flags += ["-Xptxas", "-v", "-Xptxas", "--warn-on-spills"]
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

    extension = CUDAExtension(
        f"spirulae_splat.csrc",
        sources,
        include_dirs=[
            extensions_dir,
            os.path.join(extensions_dir, "glm"),
        ],
        define_macros=define_macros,
        undef_macros=undef_macros,
        extra_compile_args=extra_compile_args,
        extra_link_args=extra_link_args,
    )

    return [extension]


import importlib.util
if importlib.util.find_spec('nerfstudio') is None and False:
    raise ValueError("Please make sure you have nerfstudio installed.")
if importlib.util.find_spec('torch') is None:
    raise ValueError("Please make sure you have PyTorch installed.")
no_fused_ssim = (importlib.util.find_spec('fused_ssim') is None)
no_fused_bilagrid = (importlib.util.find_spec('fused_bilagrid') is None)

path = "spirulae_splat/splat/cuda/csrc/"
generate_header(path+"SphericalHarmonics.cu", path+"SphericalHarmonics.cuh")
generate_header(path+"BackgroundSphericalHarmonics.cu", path+"BackgroundSphericalHarmonics.cuh")
generate_header(path+"PerSplatLoss.cu", path+"PerSplatLoss.cuh")
generate_header(path+"PerPixelLoss.cu", path+"PerPixelLoss.cuh")
generate_header(path+"PixelWise.cu", path+"PixelWise.cuh")
generate_header(path+"Projection.cu", path+"Projection.cuh")
generate_header(path+"ProjectionEval3D.cu", path+"ProjectionEval3D.cuh")
generate_header(path+"ProjectionEWA3DGSHetero.cu", path+"ProjectionEWA3DGSHetero.cuh")
generate_header(path+"RasterizationFwd.cu", path+"RasterizationFwd.cuh")
generate_header(path+"RasterizationBwd.cu", path+"RasterizationBwd.cuh")
generate_header(path+"RasterizationEval3DFwd.cu", path+"RasterizationEval3DFwd.cuh")
generate_header(path+"RasterizationEval3DBwd.cu", path+"RasterizationEval3DBwd.cuh")
generate_header(path+"RasterizationSortedEval3DFwd.cu", path+"RasterizationSortedEval3DFwd.cuh")
generate_header(path+"RasterizationSortedEval3DBwd.cu", path+"RasterizationSortedEval3DBwd.cuh")

setup(
    name="spirulae_splat",
    description="TODO",
    version="0.1.0",
    packages=find_packages(include=["spirulae_splat*"]),
    install_requires=[
        # "nerfstudio",
        # "torch",
        "jaxtyping",
        "rich>=12",
        "typing_extensions",
    ] + [
        # no need to consume internet bandwidth at each `pip install -e``
        "fused_ssim @ git+https://github.com/rahul-goel/fused-ssim.git",
        # "CUDA error: no kernel image is available for execution on the device" on torch==2.1.2+cu118 + A6000
        # "fused_ssim @ git+https://github.com/MrNeRF/optimized-fused-ssim.git",
    ] * no_fused_ssim + [
        "fused_bilagrid @ git+https://github.com/harry7557558/fused-bilagrid.git",
    ] * no_fused_bilagrid,
    extras_require={
        # dev dependencies. Install them by `pip install gsplat[dev]`
        "dev": [
            "black[jupyter]==22.3.0",
            "isort==5.10.1",
            "pylint==2.13.4",
            "pytest==7.1.2",
            "pytest-xdist==2.5.0",
            "typeguard>=2.13.3",
            "pyyaml==6.0",
            "build",
            "twine",
            "ninja",
        ],
    },
    entry_points={
        "nerfstudio.method_configs": [
            "spirulae = spirulae_splat.ns_config:spirulae",
            "spirulae-patched = spirulae_splat.ns_config:spirulae_patched",
            "spirulae-triangle = spirulae_splat.ns_config:spirulae_triangle",
        ]
    },
    cmdclass={"build_ext": get_ext()} if not BUILD_NO_CUDA else {},
    ext_modules=get_extensions() if not BUILD_NO_CUDA else [],
)
