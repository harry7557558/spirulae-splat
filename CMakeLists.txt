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

# ---------------------------------------------------------------------------
# Compute backend selection
#
# cuda (default): the full build, exactly as before.
# vulkan: WORK IN PROGRESS -- no kernel implementations yet. Compiles the
# portable engine layer (Engine*.cpp + host support) against the backend API
# (csrc/backend/, see its README.md) with SSPLAT_BACKEND_VULKAN defined and
# WITHOUT the CUDA toolkit. Produces a compile-check object library only --
# no executables and no Python extension.
# ---------------------------------------------------------------------------
set(SSPLAT_BACKEND "cuda" CACHE STRING "Compute backend: cuda | vulkan (vulkan is WIP)")
set_property(CACHE SSPLAT_BACKEND PROPERTY STRINGS cuda vulkan)

if(SSPLAT_BACKEND STREQUAL "vulkan")
    message(STATUS "SSPLAT_BACKEND=vulkan (work in progress): compile-checking "
        "the portable engine layer; no executables are produced.")

    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)

    file(GLOB SSPLAT_PORTABLE_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/spirulae_splat/splat/cuda/csrc/*.cpp)
    list(REMOVE_ITEM SSPLAT_PORTABLE_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/spirulae_splat/splat/cuda/csrc/ext.cpp)

    add_library(csrc_portable OBJECT ${SSPLAT_PORTABLE_SOURCES})
    target_compile_definitions(csrc_portable PUBLIC SSPLAT_BACKEND_VULKAN)
    target_include_directories(csrc_portable PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/spirulae_splat/splat/cuda/csrc)

    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(csrc_portable PUBLIC OpenMP::OpenMP_CXX)
    endif()

    return()  # the rest of this file is the CUDA build
endif()

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
    target_link_libraries(ssplat-train PRIVATE csrc CUDA::cudart)
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
        csrc CUDA::cudart imgui_glfw OpenGL::GL)
    if(SSPLAT_WITH_TORCH)
        find_package(Python3 REQUIRED COMPONENTS Development.Embed)
        target_link_libraries(ssplat-gui PRIVATE Python3::Python)
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

message(STATUS "spirulae_splat CUDA extension configured.")
