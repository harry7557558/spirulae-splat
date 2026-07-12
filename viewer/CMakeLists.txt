# Standalone WebGL viewer — WASM module build.
#
# Configure with the Emscripten toolchain via emcmake, e.g.:
#   emcmake cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
#   cmake --build build
# Output: build/ssv_wasm.js + build/ssv_wasm.wasm, copied to js/.
#
# Two kinds of native code are linked in:
#   1. src/viewer.cpp — dependency-free, performance-critical PLY/OBJ parsing,
#      depth sorting, and statistics (kept -ffast-math / -fno-exceptions).
#   2. The trainer's dataset parsers (COLMAP / Nerfstudio / Metashape), reused
#      IN-PLACE from spirulae_splat/splat/cuda/csrc/app (referenced by relative
#      path, not copied) and driven from src/dataset_bridge.cpp over MEMFS.
# GLTF/GLB are parsed in JavaScript.

cmake_minimum_required(VERSION 3.16)
project(ssv_viewer_wasm CXX C)

if(NOT EMSCRIPTEN)
  message(FATAL_ERROR "This project must be configured with emcmake (Emscripten toolchain).")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Reused trainer sources (in-place; carried in the same git repo).
set(CSRC "${CMAKE_SOURCE_DIR}/../spirulae_splat/splat/cuda/csrc")
set(APP  "${CSRC}/app")

set(DATASET_SOURCES
  "${CMAKE_SOURCE_DIR}/src/dataset_bridge.cpp"
  "${APP}/ColmapParser.cpp"
  "${APP}/NerfstudioParser.cpp"
  "${APP}/MetashapeParser.cpp"
  "${APP}/DatasetCommon.cpp"
  "${CSRC}/external/miniz.c"
)

add_executable(ssv_wasm src/viewer.cpp ${DATASET_SOURCES})

# dataset_bridge.cpp includes DatasetParser.h / Json.h / Xml.h (the parser
# sources find their own headers via their in-tree relative includes).
target_include_directories(ssv_wasm PRIVATE "${APP}")

# Per-source flags: viewer.cpp stays fast/exception-free; the dataset parsers
# need C++ exceptions (they throw on bad input) and must NOT use -ffast-math
# (they rely on inf/NaN semantics, e.g. outlier_threshold = infinity()).
set_source_files_properties(src/viewer.cpp PROPERTIES
  COMPILE_OPTIONS "-O3;-ffast-math;-fno-exceptions;-fno-rtti")
set_source_files_properties(${DATASET_SOURCES} PROPERTIES
  COMPILE_OPTIONS "-O3;-fexceptions")

# Emscripten link options.
#  - MODULARIZE + EXPORT_NAME: load as a factory function (ES module friendly).
#  - ALLOW_MEMORY_GROWTH + MAXIMUM_MEMORY: support very large (multi-GB) models.
#  - FORCE_FILESYSTEM + FS: mount dropped dataset files into MEMFS.
#  - -fexceptions: enable exception support for the dataset parsers.
#  - EXPORTED_RUNTIME_METHODS: cwrap/ccall + heap views + FS used from JS.
#  - The exported C functions are auto-kept via EMSCRIPTEN_KEEPALIVE.
target_link_options(ssv_wasm PRIVATE
  -O3
  -fexceptions
  "-sMODULARIZE=1"
  "-sEXPORT_NAME=SSVModule"
  "-sEXPORT_ES6=1"
  "-sENVIRONMENT=web,worker"
  "-sALLOW_MEMORY_GROWTH=1"
  "-sMAXIMUM_MEMORY=4gb"
  "-sINITIAL_MEMORY=64mb"
  "-sSTACK_SIZE=1mb"
  "-sALLOW_TABLE_GROWTH=1"
  "-sFORCE_FILESYSTEM=1"
  "-sEXPORTED_RUNTIME_METHODS=['ccall','cwrap','FS','HEAPU8','HEAPU16','HEAPU32','HEAPF32','HEAP32']"
)

set_target_properties(ssv_wasm PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
  OUTPUT_NAME "ssv_wasm"
  SUFFIX ".js"
)

# Copy build artifacts next to the JS sources so index.html can load them.
add_custom_command(TARGET ssv_wasm POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
          "${CMAKE_BINARY_DIR}/ssv_wasm.js" "${CMAKE_SOURCE_DIR}/js/ssv_wasm.js"
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
          "${CMAKE_BINARY_DIR}/ssv_wasm.wasm" "${CMAKE_SOURCE_DIR}/js/ssv_wasm.wasm"
  COMMENT "Copying WASM artifacts into js/"
)
