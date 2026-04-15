#!/bin/bash

# script_dir=$(dirname "$(realpath "$0")")
# build_dir=${script_dir}/.torch_extensions

# mkdir -p $build_dir
# TMPDIR=$build_dir \
#     pip install -e . --no-build-isolation -v \
#     --cache-dir $build_dir

python3 spirulae_splat/generate_headers.py
python3 spirulae_splat/generate_kernel_instantiation.py

cmake -G Ninja -B build
cmake --build build --verbose #-j8
mv build/libcsrc.so ./spirulae_splat/csrc.so
