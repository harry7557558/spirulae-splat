#!/usr/bin/env bash
# Build the WASM module for the standalone viewer.
# Requires the Emscripten SDK to be activated (emcc / emcmake on PATH), e.g.:
#   source ~/emsdk/emsdk_env.sh
set -euo pipefail
cd "$(dirname "$0")"

BUILD_DIR=build
emcmake cmake -S . -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release
cmake --build "$BUILD_DIR" -j

echo "Built js/ssv_wasm.js + js/ssv_wasm.wasm"
