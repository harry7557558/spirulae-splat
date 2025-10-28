import os
import re

cu_dir = "csrc/generated"

cu_files = [os.path.join(cu_dir, f) for f in os.listdir(cu_dir) if f.endswith(".cu")]

def process(src: str):
    src = src.split('\n')
    for i, line in enumerate(src):
        if line.startswith('__device__') and 'inline' not in line:
            line = 'inline ' + line  # make linker happy
            # line = '__forceinline__ ' + line  # why?
        if line.strip().startswith("static_assert(false,"):
            line = "//" + line
        src[i] = line
    src = '\n'.join(src)
    src = "#pragma once\n\n" + src
    return src

header_filename = os.path.join(cu_dir, 'slang.cuh')
header = open(header_filename).read()

for filename in cu_files:
    src = open(filename, 'r').read()
    if src[:len(header)] == header:
        src = src[len(header):]
        src = """#include "slang.cuh"\n\n""" + src
        src = process(src)
    open(filename, 'w').write(src)

match = re.match(r'\#include "(.*?slang-cuda-prelude.h)"', header)
if match:
    header = header.replace(match[0], open(match[1]).read())

header = process(header)
open(header_filename, 'w').write(header)
