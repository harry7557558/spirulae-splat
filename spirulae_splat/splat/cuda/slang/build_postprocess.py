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
        for w in ['double', 'longlong', 'ulonglong']:
            if f"({w})" in line or f', {w}, ' in line:
                line = "//" + line
                break
            if w in line and src[i-1] == 'template<>' and src[i+1] == '{' and src[i+3] == '};':
                for j in range(i-1, i+4):
                    src[j] = '//' + src[j]
                line = "//" + line
                break
        src[i] = line
    src = '\n'.join(src)
    # for (w0, w1) in [
    #     ('double4', 'my_double4'),
    #     ('longlong4', 'my_longlong4'),
    #     ('ulonglong4', 'my_ulonglong4'),
    # ]:
    #     src = re.sub(r'\b' + re.escape(w0) + r'\b', w1, src)
    return src

header_filename = os.path.join(cu_dir, 'slang.cuh')
header = open(header_filename).read()

replace_header = """#pragma once

// #if defined(CUDART_VERSION) && CUDART_VERSION >= 12000
// typedef double4_32a my_double4;
// typedef longlong4_32a my_longlong4;
// typedef ulonglong4_32a my_ulonglong4;
// #else
// typedef double4 my_double4;
// typedef longlong4 my_longlong4;
// typedef ulonglong4 my_ulonglong4;
// #endif

#include "slang.cuh"

"""

for filename in cu_files:
    src = open(filename, 'r').read()
    if src[:len(header)] == header:
        src = src[len(header):]
        src = process(src)
        src = replace_header + src
    open(filename, 'w').write(src)

match = re.match(r'\#include "(.*?slang-cuda-prelude.h)"', header)
if match:
    header = header.replace(match[0], open(match[1]).read())

header = process(header)
open(header_filename, 'w').write(header)
