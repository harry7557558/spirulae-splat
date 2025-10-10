import os

cu_dir = "csrc/generated"

cu_files = [os.path.join(cu_dir, f) for f in os.listdir(cu_dir) if f.endswith(".cu")]

def process(src: str):
    src = src.replace("\n__device__", "\ninline __device__")  # make linker happy
    return src

header_filename = os.path.join(cu_dir, 'slang.cuh')
header = open(header_filename).read()

for filename in cu_files:
    src = open(filename, 'r').read()
    if src[:len(header)] == header:
        src = src[len(header):]
        src = process(src)
        src = """#include "slang.cuh"\n\n""" + src
    open(filename, 'w').write(src)

header = process(header)
open(header_filename, 'w').write(header)
