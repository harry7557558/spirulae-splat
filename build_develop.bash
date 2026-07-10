#!/bin/bash

# script_dir=$(dirname "$(realpath "$0")")
# build_dir=${script_dir}/.torch_extensions

# mkdir -p $build_dir
# TMPDIR=$build_dir \
#     pip install -e . --no-build-isolation -v \
#     --cache-dir $build_dir

python3 spirulae_splat/generate_headers.py
python3 spirulae_splat/generate_kernel_instantiation.py
python3 spirulae_splat/generate_cli_config.py

# if [ ! -d "build" ]; then
    cmake -G Ninja -B build
# fi

echo ""

JOB_RAM_MB=750   # 750MB per job
AVAILABLE_KB=$(grep MemAvailable /proc/meminfo | awk '{print $2}')
AVAILABLE_MB=$(( AVAILABLE_KB / 1024 ))
MAX_JOBS_FROM_RAM=$(( AVAILABLE_MB / JOB_RAM_MB ))
CPU_CORES=$(nproc)
JOBS=$(( MAX_JOBS_FROM_RAM < CPU_CORES ? MAX_JOBS_FROM_RAM : CPU_CORES ))
[ "$JOBS" -lt 1 ] && JOBS=1
echo "Available RAM : ${AVAILABLE_MB} MB"
echo "CPU cores     : ${CPU_CORES}"
echo "Using jobs    : ${JOBS}"
echo ""
cmake --build build --verbose -j"${JOBS}"

mv build/libcsrc.so ./spirulae_splat/csrc.so
# keep ssplat-train's $ORIGIN lookup working after the rename
ln -sfr ./spirulae_splat/csrc.so build/libcsrc.so
