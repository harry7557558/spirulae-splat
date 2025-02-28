import re


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
        (?:void|int[234]?|float[234]?|torch::Tensor|std::tuple<[\w:\s*&<>\[\],\/]+?>)
        # Match the function name
        \s+\b\w+\b\s*
        # Match the function parameters
        \([^)]*\)
    """, re.MULTILINE | re.VERBOSE | re.DOTALL)
    
    matches = function_decl_pattern.findall(code)
    decls = []

    for m in matches:
        if re.compile(r"inline\s+(__global__|__device__)").findall(m):
            continue
        decls.append(m.strip()+';')
    
    return decls


def generate_header(source_filename, header_filename):
    code = open(source_filename).read()
    decls = extract_function_declarations(code)

    splitter = "/* == AUTO HEADER GENERATOR - DO NOT CHANGE THIS LINE == */\n"
    include = open(header_filename).read()
    include = include.split(splitter)[0].strip()

    header = '\n\n\n'.join([include, splitter]+decls)
    open(header_filename, "w").write(header+'\n')


def main():
    path = "spirulae_splat/splat/cuda/csrc/"
    generate_header(path+"projection.cu", path+"projection.cuh")
    generate_header(path+"rasterization.cu", path+"rasterization.cuh")
    generate_header(path+"rasterization_sorted.cu", path+"rasterization_sorted.cuh")
    generate_header(path+"bindings.cu", path+"bindings.h")


if  __name__ == "__main__":
    main()
