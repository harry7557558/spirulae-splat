#!/bin/bash

# Development build for Linux. Extra arguments are passed to CMake, e.g.
#   ./build_develop.bash -DSSPLAT_NO_TORCH=ON
# builds only the standalone ssplat-train (no Torch/Python needed).

# Regenerate headers/config. Skipped when python3 is unavailable -- the
# generated files are committed, so the build still works without it.
if command -v python3 >/dev/null 2>&1; then
    python3 tools/codegen/generate_headers.py
    python3 tools/codegen/generate_kernel_instantiation.py
    python3 tools/codegen/generate_cli_config.py
else
    echo "python3 not found -- skipping codegen (using committed generated files)"
fi

cmake -G Ninja -B build "$@" || exit $?

echo ""

JOB_RAM_MB=750   # 750MB per job
if [ "$(uname)" = "Darwin" ]; then
    # macOS has no MemAvailable; approximate with free+inactive pages.
    PAGE_SIZE=$(sysctl -n hw.pagesize)
    FREE_PAGES=$(vm_stat | awk -F'[:.]' '/Pages (free|inactive)/ {gsub(/ /,"",$2); sum+=$2} END {print sum}')
    AVAILABLE_MB=$(( FREE_PAGES * PAGE_SIZE / 1048576 ))
    CPU_CORES=$(sysctl -n hw.ncpu)
else
    AVAILABLE_KB=$(grep MemAvailable /proc/meminfo | awk '{print $2}')
    AVAILABLE_MB=$(( AVAILABLE_KB / 1024 ))
    CPU_CORES=$(nproc)
fi
MAX_JOBS_FROM_RAM=$(( AVAILABLE_MB / JOB_RAM_MB ))
JOBS=$(( MAX_JOBS_FROM_RAM < CPU_CORES ? MAX_JOBS_FROM_RAM : CPU_CORES ))
[ "$JOBS" -lt 1 ] && JOBS=1
echo "Available RAM : ${AVAILABLE_MB} MB"
echo "CPU cores     : ${CPU_CORES}"
echo "Using jobs    : ${JOBS}"
echo ""
# Propagate the build's exit status: without this the script always exits 0
# (the trailing libcsrc.so move succeeds regardless) and a failed build reads
# as a green one.
if ! cmake --build build --verbose -j"${JOBS}"; then
    echo "BUILD FAILED" >&2
    exit 1
fi

# Torch build only: expose the extension as spirulae_splat/csrc.so, with a
# symlink back so ssplat-train's $ORIGIN lookup keeps working. A no-torch
# build makes a static libcsrc.a and a self-contained exe -- nothing to move.
if [ -f build/libcsrc.so ] && [ ! -L build/libcsrc.so ]; then
    mv build/libcsrc.so ./spirulae_splat/csrc.so
    ln -sf ../spirulae_splat/csrc.so build/libcsrc.so
fi
