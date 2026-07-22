#!/bin/bash
# Fail if tracked source files contain machine-local absolute paths.
#
# This repo is public. Local dataset roots, private capture names and personal
# home directories should not be committed. Standard academic dataset names
# (Mip-NeRF 360 / 360_v2, ZipNeRF, ...) are fine -- it is the *paths* that
# leak a local layout.
#
# Usage:  bash tools/check_private_paths.sh

cd "$(dirname "$0")/.." || exit 1

PATTERN='/mnt/[a-z]/|/home/[a-z][a-z0-9_-]*/|[A-Za-z]:\\\\[Uu]sers|/media/[a-z][a-z0-9_-]*/'

# Exclusions:
#   viewer/js/ssv_wasm.js  - emscripten build output (single minified line)
#   docs/                  - the restructure proposal quotes the offenders it fixed
#   tools/check_private_paths.sh   - this file's own pattern
hits=$(git grep -nIE "$PATTERN" -- \
    ':!viewer/js/ssv_wasm.js' \
    ':!docs/' \
    ':!tools/check_private_paths.sh')

if [ -n "$hits" ]; then
    echo "Machine-local paths found in tracked files:"
    echo ""
    echo "$hits"
    echo ""
    echo "Replace them with a CLI argument or an environment variable."
    exit 1
fi

echo "OK: no machine-local paths in tracked files."
