import re
import os
from pathlib import Path


def strip_if_zero_blocks(code):
    lines = code.split('\n')
    result = []
    depth = 0  # 0 = active, >0 = inside a #if 0 dead block
    for line in lines:
        stripped = line.strip()
        if depth == 0:
            if re.match(r'#\s*if\s+0\b', stripped):
                depth = 1
            else:
                result.append(line)
        else:
            if re.match(r'#\s*if', stripped):
                depth += 1
            elif re.match(r'#\s*endif', stripped):
                depth -= 1
    return '\n'.join(result)


def extract_function_declarations(code):
    # Regex to match non-inline function declarations
    function_decl_pattern = re.compile(r"""
        \/\*\[AutoHeaderGeneratorExport\]\*\/\s*
        (.*?\))\s*\{
    """, re.MULTILINE | re.VERBOSE | re.DOTALL)
    
    matches = function_decl_pattern.findall(code)
    decls = []

    for m in matches:
        decls.append(m.strip()+';')
    
    return decls


def write_if_changed(path, new_text):
    old = None
    if os.path.exists(path):
        with open(path, "r") as f:
            old = f.read()
    if old != new_text:
        with open(path, "w") as f:
            f.write(new_text)
        return True
    return False


def generate_header(filename):
    path = "spirulae_splat/splat/cuda/csrc/"

    code = ""
    for source_filename in os.listdir(path):
        if source_filename.startswith(filename+"."):
            code += open(path+source_filename).read()
    decls = extract_function_declarations(strip_if_zero_blocks(code))

    splitter = "/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */\n"
    include = open(path + f"{filename}.cuh").read()
    include = include.split(splitter)[0].strip()

    header = '\n\n\n'.join([include, splitter]+decls) + "\n"
    cuh_path = path + f"{filename}.cuh"
    if write_if_changed(cuh_path, header):
        print("Generated", cuh_path)
        return True
    return False


def generate_headers():
    """Regenerate headers only if needed."""
    cuda_dir = Path("spirulae_splat/splat/cuda/csrc")

    header_names = [
        'IntersectTile',
        'BackgroundSphericalHarmonics',
        'PerSplatLoss',
        'PerPixelLoss',
        'PixelWise',
        'Projection',
        'ProjectionFwd',
        'ProjectionBwd',
        'ProjectionPackedFwd',
        'RasterizationFwd',
        'RasterizationBwd',
        'RasterizationEval3DFwd',
        'RasterizationEval3DBwd',
        'Optimizer',
        'FusedProjectionBwdOptim',
        'Densify',
        'BilagridUtils',
        'Visualizer',
        'FusedSSIM',
    ]

    num_generated = 0
    for name in header_names:
        if generate_header(name):
            num_generated += 1
    print(f"Generated {num_generated}/{len(header_names)} new headers")


generate_headers()
