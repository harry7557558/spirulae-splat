import os

cu_dir = "csrc/generated"

cu_files = [os.path.join(cu_dir, f) for f in os.listdir(cu_dir) if f.endswith(".cu")]

header = open(os.path.join(cu_dir, 'slang.cuh')).read()

for filename in cu_files:
    src = open(filename, 'r').read()
    if src[:len(header)] == header:
        src = src[len(header):]
        src = """#include "slang.cuh"\n\n""" + src
    open(filename, 'w').write(src)

