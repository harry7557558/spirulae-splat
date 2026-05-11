cmake_minimum_required(VERSION 3.18 FATAL_ERROR)
project(spirulae_splat LANGUAGES C CXX CUDA)

# This file is intended to be used with `build_develop.bash` for development purpose.


# ---------------------------------------------------------------------------
# Find Torch
# ---------------------------------------------------------------------------

execute_process(
    COMMAND python3 -c "import torch; print(torch.utils.cmake_prefix_path)"
    OUTPUT_VARIABLE PYTORCH_CMAKE_PREDIX_PATH
    OUTPUT_STRIP_TRAILING_WHITESPACE)
message(STATUS "PyTorch CMake prefix path: ${PYTORCH_PATH}")

execute_process(
    COMMAND python3 -c "import sysconfig; print(sysconfig.get_paths()['include'])"
    OUTPUT_VARIABLE PYTHON_ADDITION_INCLUDE_PATH
    OUTPUT_STRIP_TRAILING_WHITESPACE)
message(STATUS "Python additional include path: ${PYTHON_ADDITION_INCLUDE_PATH}")

set(CMAKE_PREFIX_PATH ${PYTORCH_CMAKE_PREDIX_PATH})

find_package(Torch REQUIRED)
message(STATUS "Found Torch: ${TORCH_VERSION}")

find_package(Python3 REQUIRED)

find_library(TORCH_PYTHON_LIBRARY torch_python PATHS "${TORCH_INSTALL_PREFIX}/lib")


# ---------------------------------------------------------------------------
# Find all sources
# ---------------------------------------------------------------------------
file(GLOB SPLAT_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/spirulae_splat/splat/cuda/csrc/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/spirulae_splat/splat/cuda/csrc/*.cu
    ${CMAKE_CURRENT_SOURCE_DIR}/spirulae_splat/splat/cuda/ins/*.cu
)

foreach(src ${SPLAT_SOURCES})
    if(src MATCHES "hip")
        list(REMOVE_ITEM SPLAT_SOURCES ${src})
    endif()
endforeach()

# message(STATUS "Extension sources:")
# foreach(src ${SPLAT_SOURCES})
#     message(STATUS "  ${src}")
# endforeach()

# ---------------------------------------------------------------------------
# Detect CUDA architectures
# ---------------------------------------------------------------------------
if(NOT DEFINED TORCH_CUDA_ARCH_LIST)
    execute_process(
        COMMAND python3 -c "import torch; print(' '.join(f'{a[0]}{a[1]}' for a in [torch.cuda.get_device_capability(i) for i in range(torch.cuda.device_count())]))"
        OUTPUT_VARIABLE TORCH_CUDA_ARCH_LIST
        OUTPUT_STRIP_TRAILING_WHITESPACE)

    if(TORCH_CUDA_ARCH_LIST STREQUAL "")
        message(FATAL_ERROR "CUDA is required for this extension.")
    endif()
else()
    string(REGEX REPLACE "\\." "" TORCH_CUDA_ARCH_LIST "${TORCH_CUDA_ARCH_LIST}")
endif()

message(STATUS "CUDA architecture(s): ${TORCH_CUDA_ARCH_LIST}")

set(CMAKE_CUDA_ARCHITECTURES "")
foreach(arch ${TORCH_CUDA_ARCH_LIST})
    list(APPEND CMAKE_CUDA_ARCHITECTURES "${arch}")
endforeach()

# ---------------------------------------------------------------------------
# Compiler Flags
# ---------------------------------------------------------------------------

# CXX flags
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# set(SPLAT_CXX_FLAGS "-O3")
set(SPLAT_CXX_FLAGS "-g")
if(NOT WIN32)
    list(APPEND SPLAT_CXX_FLAGS "-Wno-sign-compare")
endif()

# Detect OpenMP
find_package(OpenMP)
if(OpenMP_FOUND AND NOT APPLE)
    message(STATUS "Compiling with OpenMP")
    list(APPEND SPLAT_CXX_FLAGS "-DAT_PARALLEL_OPENMP")
    if(WIN32)
        list(APPEND SPLAT_CXX_FLAGS "/openmp")
    else()
        list(APPEND SPLAT_CXX_FLAGS "-fopenmp")
    endif()
else()
    message(STATUS "Compiling without OpenMP...")
endif()

# macOS ARM64 special case
if(APPLE AND CMAKE_SYSTEM_PROCESSOR MATCHES "arm64")
    message(STATUS "Enabling macOS ARM64 flags")
    add_compile_options("-arch" "arm64")
    add_link_options("-arch" "arm64")
endif()

# nvcc flags (mirrors setup.py)
set(SPLAT_NVCC_FLAGS
    "-O3" "--use_fast_math"
    "--expt-relaxed-constexpr"
    "-Xcudafe=--diag_suppress=20012"
    "-Xcudafe=--diag_suppress=550"
    "--threads" "0"
    "-lineinfo" "--generate-line-info" "--source-in-ptx"
    "-Xptxas" "--warn-on-double-precision-use"
)

# Add gencode flags
foreach(arch ${TORCH_CUDA_ARCH_LIST})
    list(APPEND SPLAT_NVCC_FLAGS "-gencode" "arch=compute_${arch},code=sm_${arch}")
endforeach()

# Host side optimizations
if(NOT WIN32)
    # list(APPEND SPLAT_NVCC_FLAGS "-Xcompiler=-O3,-march=native")
    list(APPEND SPLAT_NVCC_FLAGS "-Xcompiler=-g")
endif()

# ---------------------------------------------------------------------------
# Define the extension module
# ---------------------------------------------------------------------------

add_library(csrc SHARED ${SPLAT_SOURCES})

# Include directories
target_include_directories(csrc PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/spirulae_splat/splat/cuda/csrc
    ${CMAKE_CURRENT_SOURCE_DIR}/spirulae_splat/splat/cuda/csrc/glm
    ${Python3_INCLUDE_DIRS}
    ${PYTHON_ADDITION_INCLUDE_PATH}
)

# Definitions
if(WIN32)
    target_compile_definitions(csrc PRIVATE spirulae_splat_EXPORTS)
endif()

target_compile_definitions(csrc PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:TORCH_EXTENSION_NAME=csrc>
)

# Link Torch
target_link_libraries(csrc
    ${TORCH_LIBRARIES}
    ${TORCH_PYTHON_LIBRARY}
    ${Python3_LIBRARIES}
)

# Apply flags
target_compile_options(csrc PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:${SPLAT_CXX_FLAGS}>
    $<$<COMPILE_LANGUAGE:CUDA>:${SPLAT_NVCC_FLAGS}>
)

# Optional: strip symbols (setup.py uses '-s' if WITH_SYMBOLS==False)
if(DEFINED WITH_SYMBOLS AND NOT WITH_SYMBOLS)
    if(NOT WIN32)
        target_link_options(csrc PRIVATE "-s")
    endif()
endif()

# Required by PyTorch C++ extensions
set_property(TARGET csrc PROPERTY CXX_STANDARD 17)
set_property(TARGET csrc PROPERTY CUDA_STANDARD 17)

message(STATUS "spirulae_splat CUDA extension configured.")
