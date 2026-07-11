# Standalone WebGL viewer — WASM module build.
#
# Configure with the Emscripten toolchain via emcmake, e.g.:
#   emcmake cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
#   cmake --build build
# Output: build/ssv_wasm.js + build/ssv_wasm.wasm, copied to js/.
#
# This builds only the performance-critical native code (parsing / sorting /
# statistics). The rest of the viewer is dependency-free static HTML/JS.

cmake_minimum_required(VERSION 3.16)
project(ssv_viewer_wasm CXX)

if(NOT EMSCRIPTEN)
  message(FATAL_ERROR "This project must be configured with emcmake (Emscripten toolchain).")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

add_executable(ssv_wasm src/viewer.cpp)

target_compile_options(ssv_wasm PRIVATE
  -O3
  -ffast-math
  -fno-exceptions
  -fno-rtti
)

# Emscripten link options.
#  - MODULARIZE + EXPORT_NAME: load as a factory function (ES module friendly).
#  - ALLOW_MEMORY_GROWTH + MAXIMUM_MEMORY: support very large (multi-GB) models.
#  - EXPORTED_RUNTIME_METHODS: cwrap/ccall + heap views used from JS.
#  - The exported C functions are auto-kept via EMSCRIPTEN_KEEPALIVE.
target_link_options(ssv_wasm PRIVATE
  -O3
  "-sMODULARIZE=1"
  "-sEXPORT_NAME=SSVModule"
  "-sEXPORT_ES6=1"
  "-sENVIRONMENT=web,worker"
  "-sALLOW_MEMORY_GROWTH=1"
  "-sMAXIMUM_MEMORY=4gb"
  "-sINITIAL_MEMORY=64mb"
  "-sSTACK_SIZE=1mb"
  "-sALLOW_TABLE_GROWTH=1"
  "-sEXPORTED_RUNTIME_METHODS=['ccall','cwrap','HEAPU8','HEAPU16','HEAPU32','HEAPF32','HEAP32']"
  "-sFILESYSTEM=0"
  "-sDISABLE_EXCEPTION_CATCHING=1"
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
