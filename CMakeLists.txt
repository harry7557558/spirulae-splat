cmake_minimum_required(VERSION 3.18 FATAL_ERROR)
project(spirulae_splat LANGUAGES C CXX)  # CUDA enabled below unless SSPLAT_BACKEND=vulkan

# This file is intended to be used with `build_develop.bash` for development purpose.


# ---------------------------------------------------------------------------
# Find Torch (optional)
#
# Torch + Python are only needed for the Python extension module (ext.cpp).
# When either is missing -- or SSPLAT_NO_TORCH=ON -- the extension is skipped
# and only the standalone ssplat-train CLI is built (engine is torch-free).
# ---------------------------------------------------------------------------

option(SSPLAT_BUILD_CLI "Build the standalone ssplat-train CLI" OFF)
option(SSPLAT_BUILD_GUI "Build the native GUI app ssplat-gui (fetches GLFW + Dear ImGui)" OFF)
option(SSPLAT_NO_TORCH "Skip Torch/Python even if present; build only ssplat-train" OFF)

# Debug symbols / line info are OFF by default: they bloat the binaries
# massively (nvcc host -g, CUDA cubin lineinfo/source-in-ptx, and slangc -g2
# in the embedded SPIR-V). Turn on for profiling/debugging builds. The pip
# build (setup.py) is independent and unaffected (it has its own WITH_SYMBOLS
# and LINE_INFO env gates).
option(SSPLAT_DEBUG_SYMBOLS "Emit debug symbols / line info (host -g, CUDA cubin lineinfo, SPIR-V -g2)" OFF)

# ---------------------------------------------------------------------------
# Compute backend selection
#
# cuda (default): the full build, exactly as before.
# vulkan: the portable engine layer (Engine*.cpp + host support) against the
# Vulkan compute runtime (csrc/backend/vulkan/, see its README.md), built
# WITHOUT the CUDA toolkit. Produces the backend tests, the ssplat-train CLI
# (always) and ssplat-gui (with SSPLAT_BUILD_GUI=ON); no Python extension.
# ---------------------------------------------------------------------------
set(SSPLAT_BACKEND "cuda" CACHE STRING "Compute backend: cuda | vulkan")
set_property(CACHE SSPLAT_BACKEND PROPERTY STRINGS cuda vulkan)

if(SSPLAT_BACKEND STREQUAL "vulkan")
    message(STATUS "SSPLAT_BACKEND=vulkan: portable engine layer + "
        "backend/vulkan runtime; kernel coverage is tracked in "
        "csrc/backend/vulkan/README.md.")
    # The Python extension is CUDA-only, so the CLI trainer is this build's
    # primary artifact (same forcing as the CUDA no-torch build).
    set(SSPLAT_BUILD_CLI ON)

    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)

    set(SSPLAT_CSRC ${CMAKE_CURRENT_SOURCE_DIR}/spirulae_splat/splat/cuda/csrc)

    find_package(Vulkan REQUIRED)

    # Portable engine layer (compile-checked; links against the Vulkan
    # runtime once kernel launchers exist for the symbols it needs).
    file(GLOB SSPLAT_PORTABLE_SOURCES CONFIGURE_DEPENDS ${SSPLAT_CSRC}/*.cpp)
    list(REMOVE_ITEM SSPLAT_PORTABLE_SOURCES ${SSPLAT_CSRC}/ext.cpp)

    add_library(csrc_portable OBJECT ${SSPLAT_PORTABLE_SOURCES})
    target_compile_definitions(csrc_portable PUBLIC SSPLAT_BACKEND_VULKAN)
    target_include_directories(csrc_portable PUBLIC ${SSPLAT_CSRC})

    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(csrc_portable PUBLIC OpenMP::OpenMP_CXX)
    endif()

    # ---- slangc: locate a matching compiler or fetch the pinned release ----
    # SPIR-V blobs are never committed; they are compiled at build time (one
    # slangc edge per blob, see slang/SpirvShaders.cmake) and embedded below.
    # No Python is required for the Vulkan build.
    set(SSPLAT_SLANG_VERSION "2026.12.0.1")
    set(SSPLAT_SLANGC "" CACHE FILEPATH
        "slangc to use (empty = find in PATH, fetch on miss/mismatch)")

    function(_ssplat_slangc_version exe out_var)
        set(${out_var} "" PARENT_SCOPE)
        if(NOT exe OR NOT EXISTS ${exe})
            return()
        endif()
        execute_process(COMMAND ${exe} -v
            OUTPUT_VARIABLE _v ERROR_VARIABLE _v
            OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_STRIP_TRAILING_WHITESPACE
            RESULT_VARIABLE _r)
        if(_r EQUAL 0)
            string(STRIP "${_v}" _v)
            set(${out_var} "${_v}" PARENT_SCOPE)
        endif()
    endfunction()

    set(_slangc "${SSPLAT_SLANGC}")
    if(NOT _slangc)
        find_program(_slangc_in_path slangc)
        if(_slangc_in_path)
            set(_slangc ${_slangc_in_path})
        endif()
    endif()
    _ssplat_slangc_version("${_slangc}" _slangc_ver)
    if(NOT _slangc_ver STREQUAL SSPLAT_SLANG_VERSION)
        if(_slangc)
            message(STATUS "slangc at ${_slangc} is version '${_slangc_ver}' "
                "(need ${SSPLAT_SLANG_VERSION}); fetching the pinned release")
        else()
            message(STATUS "slangc not found; fetching the pinned release "
                "${SSPLAT_SLANG_VERSION}")
        endif()
        if(WIN32)
            set(_slang_plat "windows-x86_64")
            set(_slang_ext "zip")
            set(_slang_exe "bin/slangc.exe")
        else()
            if(CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "aarch64|arm64")
                set(_slang_arch "aarch64")
            else()
                set(_slang_arch "x86_64")
            endif()
            if(CMAKE_HOST_SYSTEM_NAME STREQUAL "Darwin")
                set(_slang_plat "macos-${_slang_arch}")
            else()
                set(_slang_plat "linux-${_slang_arch}")
            endif()
            set(_slang_ext "tar.gz")
            set(_slang_exe "bin/slangc")
        endif()
        set(_slang_name "slang-${SSPLAT_SLANG_VERSION}-${_slang_plat}")
        set(_slang_dir ${CMAKE_CURRENT_BINARY_DIR}/${_slang_name})
        set(_slangc ${_slang_dir}/${_slang_exe})
        if(NOT EXISTS ${_slangc})
            set(_slang_url "https://github.com/shader-slang/slang/releases/download/v${SSPLAT_SLANG_VERSION}/${_slang_name}.${_slang_ext}")
            set(_slang_archive ${CMAKE_CURRENT_BINARY_DIR}/${_slang_name}.${_slang_ext})
            message(STATUS "Downloading ${_slang_url}")
            file(DOWNLOAD ${_slang_url} ${_slang_archive}
                SHOW_PROGRESS STATUS _dl_status)
            list(GET _dl_status 0 _dl_code)
            if(NOT _dl_code EQUAL 0)
                message(FATAL_ERROR "Failed to download slang release: "
                    "${_dl_status}. Install slang-${SSPLAT_SLANG_VERSION} "
                    "and pass -DSSPLAT_SLANGC=/path/to/slangc.")
            endif()
            file(ARCHIVE_EXTRACT INPUT ${_slang_archive}
                DESTINATION ${_slang_dir})
            file(REMOVE ${_slang_archive})
        endif()
        _ssplat_slangc_version("${_slangc}" _slangc_ver)
        if(NOT _slangc_ver STREQUAL SSPLAT_SLANG_VERSION)
            message(FATAL_ERROR "Fetched slangc reports version "
                "'${_slangc_ver}', expected ${SSPLAT_SLANG_VERSION}")
        endif()
    endif()
    message(STATUS "Using slangc ${_slangc} (${_slangc_ver})")

    # ---- compile Slang -> SPIR-V (build tree) and embed ----
    # Native (no Python): slang/spirv_tool.cpp is compiled by the module below,
    # which enumerates one slangc custom command per blob (so -j bounds RAM) and
    # embeds the result. SPIR-V debug info (slangc -g2) is opt-in via
    # SSPLAT_DEBUG_SYMBOLS; it otherwise dominates the embedded blob size.
    set(SSPLAT_CUDA_DIR ${CMAKE_CURRENT_SOURCE_DIR}/spirulae_splat/splat/cuda)
    set(SSPLAT_SPIRV_DIR ${CMAKE_CURRENT_BINARY_DIR}/spirv)
    set(SSPLAT_SPIRV_EMBED ${CMAKE_CURRENT_BINARY_DIR}/vk_shaders_embedded.cpp)
    include(${SSPLAT_CUDA_DIR}/slang/SpirvShaders.cmake)
    ssplat_setup_spirv(${SSPLAT_CUDA_DIR} ${_slangc} ${SSPLAT_SPIRV_DIR}
        ${SSPLAT_SPIRV_EMBED} ${SSPLAT_DEBUG_SYMBOLS})

    # Vulkan backend runtime (device context, memory, streams, events) and
    # kernel launchers.
    file(GLOB SSPLAT_VULKAN_SOURCES CONFIGURE_DEPENDS ${SSPLAT_CSRC}/backend/vulkan/*.cpp
        ${SSPLAT_CSRC}/backend/vulkan/kernels/*.cpp)
    list(APPEND SSPLAT_VULKAN_SOURCES ${SSPLAT_SPIRV_EMBED})
    add_library(ssplat_backend_vulkan STATIC ${SSPLAT_VULKAN_SOURCES})
    target_compile_definitions(ssplat_backend_vulkan PUBLIC SSPLAT_BACKEND_VULKAN)
    target_include_directories(ssplat_backend_vulkan PUBLIC ${SSPLAT_CSRC})
    target_link_libraries(ssplat_backend_vulkan PUBLIC Vulkan::Vulkan)

    # Backend tests (run manually; see backend/vulkan/README.md).
    file(GLOB SSPLAT_VULKAN_TESTS CONFIGURE_DEPENDS ${SSPLAT_CSRC}/backend/vulkan/tests/*.cpp)
    foreach(test_src ${SSPLAT_VULKAN_TESTS})
        get_filename_component(test_name ${test_src} NAME_WE)
        add_executable(${test_name} ${test_src})
        target_link_libraries(${test_name} PRIVATE ssplat_backend_vulkan)
    endforeach()

    # Cross-backend parity tools (same sources build in the CUDA branch).
    file(GLOB SSPLAT_PARITY_TESTS CONFIGURE_DEPENDS
        ${SSPLAT_CSRC}/backend/tests/*.cpp)
    foreach(test_src ${SSPLAT_PARITY_TESTS})
        get_filename_component(test_name ${test_src} NAME_WE)
        add_executable(${test_name} ${test_src})
        target_link_libraries(${test_name} PRIVATE ssplat_backend_vulkan)
    endforeach()

    # Engine-level parity tools: drive the real engine (forward_3dgs +
    # background + blit) through the portable layer. The backend lib
    # provides the ported kernels plus throwing stubs
    # (kernels/TrainingStubs*.cpp) for everything still CUDA-only.
    find_package(Threads REQUIRED)
    file(GLOB SSPLAT_ENGINE_TESTS CONFIGURE_DEPENDS
        ${SSPLAT_CSRC}/backend/tests/engine/*.cpp)
    foreach(test_src ${SSPLAT_ENGINE_TESTS})
        get_filename_component(test_name ${test_src} NAME_WE)
        add_executable(${test_name} ${test_src})
        target_link_libraries(${test_name} PRIVATE csrc_portable
            ssplat_backend_vulkan Threads::Threads)
    endforeach()

    # Torch/Python extension is CUDA-only; the app targets below build
    # against the portable engine + Vulkan backend instead.
    set(SSPLAT_WITH_TORCH OFF)
    set(SSPLAT_APP_LIBS csrc_portable ssplat_backend_vulkan Threads::Threads)
else()

enable_language(CUDA)

set(SSPLAT_WITH_TORCH OFF)
if(NOT SSPLAT_NO_TORCH)
    find_package(Python3 QUIET COMPONENTS Interpreter)
    if(Python3_Interpreter_FOUND)
        execute_process(
            COMMAND ${Python3_EXECUTABLE} -c "import torch; print(torch.utils.cmake_prefix_path)"
            OUTPUT_VARIABLE PYTORCH_CMAKE_PREFIX_PATH
            RESULT_VARIABLE SSPLAT_TORCH_PROBE_RESULT
            ERROR_QUIET
            OUTPUT_STRIP_TRAILING_WHITESPACE)
        message(STATUS "PyTorch CMake prefix path: ${PYTORCH_CMAKE_PREFIX_PATH}")

        execute_process(
            COMMAND ${Python3_EXECUTABLE} -c "import sysconfig; print(sysconfig.get_paths()['include'])"
            OUTPUT_VARIABLE PYTHON_ADDITION_INCLUDE_PATH
            ERROR_QUIET
            OUTPUT_STRIP_TRAILING_WHITESPACE)
        message(STATUS "Python additional include path: ${PYTHON_ADDITION_INCLUDE_PATH}")

        if(SSPLAT_TORCH_PROBE_RESULT EQUAL 0)
            set(CMAKE_PREFIX_PATH ${PYTORCH_CMAKE_PREFIX_PATH})
            find_package(Torch QUIET)
        endif()
    endif()

    if(Torch_FOUND)
        message(STATUS "Found Torch: ${TORCH_VERSION}")
        find_library(TORCH_PYTHON_LIBRARY torch_python PATHS "${TORCH_INSTALL_PREFIX}/lib")
        set(SSPLAT_WITH_TORCH ON)
    endif()
endif()

if(NOT SSPLAT_WITH_TORCH)
    message(STATUS "Torch/Python3 not available (or SSPLAT_NO_TORCH=ON) -- "
        "skipping the Python extension; building the standalone ssplat-train CLI only.")
    set(SSPLAT_BUILD_CLI ON)
endif()

find_package(CUDAToolkit REQUIRED)


# ---------------------------------------------------------------------------
# Find all sources
# ---------------------------------------------------------------------------
# CONFIGURE_DEPENDS: kernel families are routinely split into several TUs
# (Name.Part.cu -- see docs/codegen.md), and codegen adds/removes files under
# ins/. Without it a new source is silently not compiled and only surfaces as
# an undefined reference at link time.
file(GLOB SPLAT_SOURCES CONFIGURE_DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/spirulae_splat/splat/cuda/csrc/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/spirulae_splat/splat/cuda/csrc/*.cu
    ${CMAKE_CURRENT_SOURCE_DIR}/spirulae_splat/splat/cuda/ins/*.cu
)

foreach(src ${SPLAT_SOURCES})
    if(src MATCHES "hip")
        list(REMOVE_ITEM SPLAT_SOURCES ${src})
    endif()
endforeach()

if(NOT SSPLAT_WITH_TORCH)
    # ext.cpp (the pybind11 module) is the only Torch-dependent source.
    list(REMOVE_ITEM SPLAT_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/spirulae_splat/splat/cuda/csrc/ext.cpp)
endif()

# message(STATUS "Extension sources:")
# foreach(src ${SPLAT_SOURCES})
#     message(STATUS "  ${src}")
# endforeach()

# ---------------------------------------------------------------------------
# Detect CUDA architectures
# ---------------------------------------------------------------------------
if(NOT DEFINED TORCH_CUDA_ARCH_LIST)
    if(SSPLAT_WITH_TORCH)
        execute_process(
            COMMAND ${Python3_EXECUTABLE} -c "import torch; print(' '.join(f'{a[0]}{a[1]}' for a in [torch.cuda.get_device_capability(i) for i in range(torch.cuda.device_count())]))"
            OUTPUT_VARIABLE TORCH_CUDA_ARCH_LIST
            ERROR_QUIET
            OUTPUT_STRIP_TRAILING_WHITESPACE)
    endif()

    if(NOT TORCH_CUDA_ARCH_LIST)
        # No Torch (or Torch saw no device): ask the driver directly.
        execute_process(
            COMMAND nvidia-smi --query-gpu=compute_cap --format=csv,noheader
            OUTPUT_VARIABLE SSPLAT_SMI_CAPS
            RESULT_VARIABLE SSPLAT_SMI_RESULT
            ERROR_QUIET
            OUTPUT_STRIP_TRAILING_WHITESPACE)
        if(SSPLAT_SMI_RESULT EQUAL 0)
            string(REPLACE "\r" "" SSPLAT_SMI_CAPS "${SSPLAT_SMI_CAPS}")
            string(REPLACE "\n" " " TORCH_CUDA_ARCH_LIST "${SSPLAT_SMI_CAPS}")
        endif()
    endif()

    if(NOT TORCH_CUDA_ARCH_LIST)
        message(FATAL_ERROR "Could not detect a CUDA GPU. "
            "Pass -DTORCH_CUDA_ARCH_LIST=\"8.6\" (or similar) to set the target architecture.")
    endif()
endif()

# normalize: "8.6" / "8.6 7.5" / "86;75" all become a deduped list like "86;75"
string(REGEX REPLACE "\\." "" TORCH_CUDA_ARCH_LIST "${TORCH_CUDA_ARCH_LIST}")
string(REPLACE " " ";" TORCH_CUDA_ARCH_LIST "${TORCH_CUDA_ARCH_LIST}")
list(REMOVE_DUPLICATES TORCH_CUDA_ARCH_LIST)

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

if(MSVC)
    set(SPLAT_CXX_FLAGS "/O2")
else()
    set(SPLAT_CXX_FLAGS "-O3")
    # set(SPLAT_CXX_FLAGS "-g")
    list(APPEND SPLAT_CXX_FLAGS "-Wno-sign-compare")
endif()

if(WIN32)
    add_compile_definitions(_USE_MATH_DEFINES NOMINMAX _CRT_SECURE_NO_WARNINGS)
endif()

if (MSVC)
    add_compile_options(
        $<$<COMPILE_LANGUAGE:CXX>:/Zc:preprocessor>
        $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=/Zc:preprocessor>
    )
endif()

# Detect OpenMP
find_package(OpenMP)
if(OpenMP_CXX_FOUND AND NOT APPLE)
    message(STATUS "Compiling with OpenMP")
    list(APPEND SPLAT_CXX_FLAGS "-DAT_PARALLEL_OPENMP")
    if(SSPLAT_WITH_TORCH)
        # Compile flag only: the OpenMP runtime comes from torch's bundled
        # libgomp. Linking the system one too would shadow it (CMake warns
        # about the unsafe rpath). The no-torch build instead links the
        # OpenMP::OpenMP_CXX target, whose flags are language-guarded so
        # nvcc never sees them.
        list(APPEND SPLAT_CXX_FLAGS ${OpenMP_CXX_FLAGS})
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
    "-Xptxas" "--warn-on-double-precision-use"
)

# cubin/PTX line info (adds embedded source; needed for profiling but bloats
# the fatbin). Off unless SSPLAT_DEBUG_SYMBOLS.
if(SSPLAT_DEBUG_SYMBOLS)
    list(APPEND SPLAT_NVCC_FLAGS "-lineinfo" "--generate-line-info" "--source-in-ptx")
endif()

# Add gencode flags
foreach(arch ${TORCH_CUDA_ARCH_LIST})
    list(APPEND SPLAT_NVCC_FLAGS "-gencode" "arch=compute_${arch},code=sm_${arch}")
endforeach()

# Host side optimizations
if(NOT WIN32)
    # list(APPEND SPLAT_NVCC_FLAGS "-Xcompiler=-O3,-march=native")
    # Host debug symbols (nvcc -Xcompiler=-g). Off unless SSPLAT_DEBUG_SYMBOLS.
    if(SSPLAT_DEBUG_SYMBOLS)
        list(APPEND SPLAT_NVCC_FLAGS "-Xcompiler=-g")
    endif()
endif()

# ---------------------------------------------------------------------------
# Define the extension module
# ---------------------------------------------------------------------------

if(SSPLAT_WITH_TORCH)
    add_library(csrc SHARED ${SPLAT_SOURCES})
else()
    # The engine API has no dllexport annotations (Linux .so exports
    # everything by default); a static lib sidesteps that on Windows and
    # makes ssplat-train a self-contained executable.
    add_library(csrc STATIC ${SPLAT_SOURCES})
endif()

# Include directories
target_include_directories(csrc PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/spirulae_splat/splat/cuda/csrc
    # ${CMAKE_CURRENT_SOURCE_DIR}/spirulae_splat/splat/cuda/csrc/glm
    ${CUDAToolkit_INCLUDE_DIRS}
)

# Definitions
if(WIN32)
    target_compile_definitions(csrc PRIVATE spirulae_splat_EXPORTS)
endif()

if(SSPLAT_WITH_TORCH)
    target_include_directories(csrc PRIVATE
        ${Python3_INCLUDE_DIRS}
        ${PYTHON_ADDITION_INCLUDE_PATH}
    )

    target_compile_definitions(csrc PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:TORCH_EXTENSION_NAME=csrc>
    )

    # Link Torch
    target_link_libraries(csrc
        ${TORCH_LIBRARIES}
        ${TORCH_PYTHON_LIBRARY}
        ${Python3_LIBRARIES}
    )
else()
    # Torch normally drags these in for us
    find_package(Threads REQUIRED)
    target_link_libraries(csrc CUDA::cudart Threads::Threads)
    if(OpenMP_CXX_FOUND AND NOT APPLE)
        target_link_libraries(csrc OpenMP::OpenMP_CXX)
    endif()
endif()

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

set(SSPLAT_APP_LIBS csrc CUDA::cudart)
endif()  # SSPLAT_BACKEND (vulkan | cuda core; app targets below are shared)

# ---------------------------------------------------------------------------
# Standalone CLI trainer -- no Python at runtime; links the csrc shared
# library for the engine. Off by default (unless Torch is unavailable, which
# forces it on) so the normal extension build is unaffected.
# ---------------------------------------------------------------------------

if(SSPLAT_BUILD_CLI OR SSPLAT_BUILD_GUI)
    # Embed the (unchanged) web viewer client into the binary so the exe is
    # self-contained. Regenerated at configure time when viewer.html changes.
    set(SSPLAT_VIEWER_HTML ${CMAKE_CURRENT_SOURCE_DIR}/spirulae_splat/viewer/viewer.html)
    set(SSPLAT_VIEWER_HTML_HDR ${CMAKE_BINARY_DIR}/app_generated/viewer_html.h)
    file(READ ${SSPLAT_VIEWER_HTML} _ssplat_vh HEX)
    string(REGEX REPLACE "([0-9a-f][0-9a-f])" "0x\\1," _ssplat_vh_bytes ${_ssplat_vh})
    file(WRITE ${SSPLAT_VIEWER_HTML_HDR}
        "#pragma once\n"
        "// AUTO-GENERATED from spirulae_splat/viewer/viewer.html -- do not edit.\n"
        "#include <cstddef>\n"
        "inline const unsigned char kViewerHtml[] = {${_ssplat_vh_bytes}};\n"
        "inline const size_t kViewerHtmlSize = sizeof(kViewerHtml);\n")
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${SSPLAT_VIEWER_HTML})
endif()

# Cross-backend parity tools (dump reference outputs; the Vulkan build's
# twin executable compares against them). See csrc/backend/vulkan/README.md.
# The Vulkan branch builds its comparing twins unconditionally above.
option(SSPLAT_BUILD_BACKEND_TESTS "Build backend parity test tools" OFF)
if(SSPLAT_BUILD_BACKEND_TESTS AND NOT SSPLAT_BACKEND STREQUAL "vulkan")
    if(SSPLAT_WITH_TORCH)
        # Same libpython + static-libstdc++ notes as ssplat-train below.
        find_package(Python3 REQUIRED COMPONENTS Development.Embed)
    endif()
    file(GLOB SSPLAT_PARITY_TESTS CONFIGURE_DEPENDS
        spirulae_splat/splat/cuda/csrc/backend/tests/*.cpp
        spirulae_splat/splat/cuda/csrc/backend/tests/engine/*.cpp)
    foreach(test_src ${SSPLAT_PARITY_TESTS})
        get_filename_component(test_name ${test_src} NAME_WE)
        add_executable(${test_name} ${test_src})
        target_include_directories(${test_name} PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/spirulae_splat/splat/cuda/csrc
            ${CUDAToolkit_INCLUDE_DIRS})
        target_link_libraries(${test_name} PRIVATE csrc CUDA::cudart)
        if(SSPLAT_WITH_TORCH)
            target_link_libraries(${test_name} PRIVATE Python3::Python)
            if(NOT WIN32)
                target_link_options(${test_name} PRIVATE
                    -static-libgcc -static-libstdc++)
            endif()
        endif()
    endforeach()
endif()

if(SSPLAT_BUILD_CLI)
    add_executable(ssplat-train
        spirulae_splat/splat/cuda/csrc/app/main.cpp
        spirulae_splat/splat/cuda/csrc/app/TrainerCore.cpp
        spirulae_splat/splat/cuda/csrc/app/ColmapParser.cpp
        spirulae_splat/splat/cuda/csrc/app/NerfstudioParser.cpp
        spirulae_splat/splat/cuda/csrc/app/MetashapeParser.cpp
        spirulae_splat/splat/cuda/csrc/app/DatasetCommon.cpp
        spirulae_splat/splat/cuda/csrc/app/HttpServer.cpp
        spirulae_splat/splat/cuda/csrc/app/RenderWorker.cpp
        spirulae_splat/splat/cuda/csrc/app/Viewer.cpp
        spirulae_splat/splat/cuda/csrc/external/miniz.c
    )
    target_include_directories(ssplat-train PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/spirulae_splat/splat/cuda/csrc
        ${CMAKE_BINARY_DIR}
        ${CUDAToolkit_INCLUDE_DIRS}
    )
    target_link_libraries(ssplat-train PRIVATE ${SSPLAT_APP_LIBS})
    if (WIN32)
        target_link_libraries(ssplat-train PRIVATE ws2_32)
    endif()
    if(SSPLAT_WITH_TORCH)
        # libpython is only needed because csrc's pybind/libtorch_python layer
        # leaves CPython symbols undefined (a Python process normally provides
        # them at import time). Drops out with the no-torch build.
        find_package(Python3 REQUIRED COMPONENTS Development.Embed)
        target_link_libraries(ssplat-train PRIVATE Python3::Python)
    endif()
    target_compile_options(ssplat-train PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:${SPLAT_CXX_FLAGS}>
    )
    # Statically link libstdc++ into the exe so its C++ runtime calls bind at
    # link time. Without this, libtorch.so (pulled in via csrc) interposes
    # its own bundled std::filesystem symbols at load time, and e.g.
    # remove_all resolves to an ABI-incompatible copy that jumps to null.
    if(SSPLAT_WITH_TORCH AND NOT WIN32)
        target_link_options(ssplat-train PRIVATE -static-libgcc -static-libstdc++)
    endif()
    # libcsrc.so is renamed to spirulae_splat/csrc.so by build_develop.bash;
    # keep the exe's lookup pointing at the build tree copy.
    if(SSPLAT_WITH_TORCH)
        set_property(TARGET ssplat-train PROPERTY BUILD_RPATH
            "$ORIGIN;${TORCH_INSTALL_PREFIX}/lib")
    else()
        set_property(TARGET ssplat-train PROPERTY BUILD_RPATH "$ORIGIN")
    endif()
    set_property(TARGET ssplat-train PROPERTY CXX_STANDARD 17)

    # ---- standalone mesh-extraction CLI (shares the dataset parsers) ----
    # CUDA-only: the Vulkan backend stubs the meshing kernels
    # (meshing::OccupancyEvaluator), so the tool would throw at startup.
    if(NOT SSPLAT_BACKEND STREQUAL "vulkan")
    add_executable(ssplat-mesh
        spirulae_splat/splat/cuda/csrc/app/mesh_main.cpp
        spirulae_splat/splat/cuda/csrc/app/ColmapParser.cpp
        spirulae_splat/splat/cuda/csrc/app/NerfstudioParser.cpp
        spirulae_splat/splat/cuda/csrc/app/MetashapeParser.cpp
        spirulae_splat/splat/cuda/csrc/app/DatasetCommon.cpp
        spirulae_splat/splat/cuda/csrc/external/miniz.c
    )
    target_include_directories(ssplat-mesh PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/spirulae_splat/splat/cuda/csrc
        ${CUDAToolkit_INCLUDE_DIRS}
    )
    target_link_libraries(ssplat-mesh PRIVATE csrc CUDA::cudart)
    if(SSPLAT_WITH_TORCH)
        target_link_libraries(ssplat-mesh PRIVATE Python3::Python)
    endif()
    target_compile_options(ssplat-mesh PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:${SPLAT_CXX_FLAGS}>
    )
    if(SSPLAT_WITH_TORCH AND NOT WIN32)
        target_link_options(ssplat-mesh PRIVATE -static-libgcc -static-libstdc++)
    endif()
    if(SSPLAT_WITH_TORCH)
        set_property(TARGET ssplat-mesh PROPERTY BUILD_RPATH
            "$ORIGIN;${TORCH_INSTALL_PREFIX}/lib")
    else()
        set_property(TARGET ssplat-mesh PROPERTY BUILD_RPATH "$ORIGIN")
    endif()
    set_property(TARGET ssplat-mesh PROPERTY CXX_STANDARD 17)
    endif()  # NOT vulkan (ssplat-mesh)
endif()

# ---------------------------------------------------------------------------
# Native GUI (ssplat-gui) -- optional, off by default. Dear ImGui + GLFW +
# OpenGL 3 (imgui's embedded GL loader; no GLEW/GLAD). Dependencies are
# fetched at configure time via FetchContent and built statically, so the
# only system requirements are OpenGL dev headers (Linux: libgl1-mesa-dev +
# xorg-dev for GLFW; Windows/macOS: nothing extra).
# ---------------------------------------------------------------------------

if(SSPLAT_BUILD_GUI)
    include(FetchContent)

    set(GLFW_BUILD_DOCS     OFF CACHE BOOL "" FORCE)
    set(GLFW_BUILD_TESTS    OFF CACHE BOOL "" FORCE)
    set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
    set(GLFW_INSTALL        OFF CACHE BOOL "" FORCE)
    FetchContent_Declare(glfw
        GIT_REPOSITORY https://github.com/glfw/glfw.git
        GIT_TAG        3.4
        GIT_SHALLOW    TRUE)
    FetchContent_Declare(imgui
        GIT_REPOSITORY https://github.com/ocornut/imgui.git
        GIT_TAG        v1.92.8
        GIT_SHALLOW    TRUE)
    FetchContent_MakeAvailable(glfw imgui)

    find_package(OpenGL REQUIRED)

    # imgui ships no CMakeLists; build the core + glfw/gl3 backends +
    # std::string InputText helper as a small static lib.
    add_library(imgui_glfw STATIC
        ${imgui_SOURCE_DIR}/imgui.cpp
        ${imgui_SOURCE_DIR}/imgui_draw.cpp
        ${imgui_SOURCE_DIR}/imgui_tables.cpp
        ${imgui_SOURCE_DIR}/imgui_widgets.cpp
        ${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp
        ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl3.cpp
        ${imgui_SOURCE_DIR}/misc/cpp/imgui_stdlib.cpp
    )
    target_include_directories(imgui_glfw PUBLIC
        ${imgui_SOURCE_DIR}
        ${imgui_SOURCE_DIR}/backends
        ${imgui_SOURCE_DIR}/misc/cpp
    )
    target_link_libraries(imgui_glfw PUBLIC glfw)
    set_property(TARGET imgui_glfw PROPERTY CXX_STANDARD 17)

    # Embed scripts/mask.py (AI masking helper, run via external Python) so
    # the exe is self-contained. Same mechanism as the viewer.html embed.
    set(SSPLAT_MASK_PY ${CMAKE_CURRENT_SOURCE_DIR}/scripts/mask.py)
    set(SSPLAT_MASK_PY_HDR ${CMAKE_BINARY_DIR}/app_generated/mask_py.h)
    file(READ ${SSPLAT_MASK_PY} _ssplat_mp HEX)
    string(REGEX REPLACE "([0-9a-f][0-9a-f])" "0x\\1," _ssplat_mp_bytes ${_ssplat_mp})
    file(WRITE ${SSPLAT_MASK_PY_HDR}
        "#pragma once\n"
        "// AUTO-GENERATED from scripts/mask.py -- do not edit.\n"
        "#include <cstddef>\n"
        "inline const unsigned char kMaskPy[] = {${_ssplat_mp_bytes}};\n"
        "inline const size_t kMaskPySize = sizeof(kMaskPy);\n")
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${SSPLAT_MASK_PY})

    add_executable(ssplat-gui
        spirulae_splat/splat/cuda/csrc/app/gui/GuiMain.cpp
        spirulae_splat/splat/cuda/csrc/app/gui/GuiApp.cpp
        spirulae_splat/splat/cuda/csrc/app/gui/ConfigUI.cpp
        spirulae_splat/splat/cuda/csrc/app/gui/FileDialog.cpp
        spirulae_splat/splat/cuda/csrc/app/gui/Subprocess.cpp
        spirulae_splat/splat/cuda/csrc/app/gui/ColmapRunner.cpp
        spirulae_splat/splat/cuda/csrc/app/gui/FrameSelect.cpp
        spirulae_splat/splat/cuda/csrc/app/gui/TrainRunner.cpp
        spirulae_splat/splat/cuda/csrc/app/gui/GlLoader.cpp
        spirulae_splat/splat/cuda/csrc/app/gui/NavCamera.cpp
        spirulae_splat/splat/cuda/csrc/app/gui/PreviewRenderer.cpp
        spirulae_splat/splat/cuda/csrc/app/gui/ViewportPanel.cpp
        spirulae_splat/splat/cuda/csrc/app/TrainerCore.cpp
        spirulae_splat/splat/cuda/csrc/app/RenderWorker.cpp
        spirulae_splat/splat/cuda/csrc/app/HttpServer.cpp
        spirulae_splat/splat/cuda/csrc/app/Viewer.cpp
        spirulae_splat/splat/cuda/csrc/app/ColmapParser.cpp
        spirulae_splat/splat/cuda/csrc/app/NerfstudioParser.cpp
        spirulae_splat/splat/cuda/csrc/app/MetashapeParser.cpp
        spirulae_splat/splat/cuda/csrc/app/DatasetCommon.cpp
        spirulae_splat/splat/cuda/csrc/external/miniz.c
    )
    target_include_directories(ssplat-gui PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/spirulae_splat/splat/cuda/csrc
        ${CMAKE_BINARY_DIR}
        ${CUDAToolkit_INCLUDE_DIRS}
    )
    target_link_libraries(ssplat-gui PRIVATE
        ${SSPLAT_APP_LIBS} imgui_glfw OpenGL::GL)
    if(SSPLAT_WITH_TORCH)
        find_package(Python3 REQUIRED COMPONENTS Development.Embed)
        target_link_libraries(ssplat-gui PRIVATE Python3::Python)
    endif()
    if (WIN32)
        target_link_libraries(ssplat-gui PRIVATE ws2_32)
    endif()
    target_compile_options(ssplat-gui PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:${SPLAT_CXX_FLAGS}>
    )
    # Same libtorch std::filesystem interposition workaround as ssplat-train.
    if(SSPLAT_WITH_TORCH AND NOT WIN32)
        target_link_options(ssplat-gui PRIVATE -static-libgcc -static-libstdc++)
    endif()
    if(SSPLAT_WITH_TORCH)
        set_property(TARGET ssplat-gui PROPERTY BUILD_RPATH
            "$ORIGIN;${TORCH_INSTALL_PREFIX}/lib")
    else()
        set_property(TARGET ssplat-gui PROPERTY BUILD_RPATH "$ORIGIN")
    endif()
    set_property(TARGET ssplat-gui PROPERTY CXX_STANDARD 17)
endif()

message(STATUS "spirulae_splat configured (backend: ${SSPLAT_BACKEND}).")
