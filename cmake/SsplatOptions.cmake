# Build options, backend selection, and the source-tree paths every other
# module builds on. Included first by the top-level CMakeLists.txt.

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
# SSPLAT_ROOT is the repository root (the directory holding CMakeLists.txt).
# Everything else is expressed relative to it so the modules never care where
# they were included from.
# SSPLAT_SRC is also the include root: every local #include in the native
# tree is written relative to it (e.g. `#include "core/Common.cuh"`).
set(SSPLAT_ROOT    ${CMAKE_CURRENT_SOURCE_DIR})
set(SSPLAT_SRC     ${SSPLAT_ROOT}/src)
set(SSPLAT_SHADERS ${SSPLAT_SRC}/shaders)                      # shared Slang device math
set(SSPLAT_VK_SHADERS ${SSPLAT_SRC}/backend/vulkan/shaders)   # Vulkan-only entry points

# The version the apps report with --version, read from the one place it is
# declared (pyproject.toml) so a CMake build and a pip build never disagree.
set(SSPLAT_VERSION "dev")
if(EXISTS ${SSPLAT_ROOT}/pyproject.toml)
    file(READ ${SSPLAT_ROOT}/pyproject.toml _ssplat_pyproject)
    string(REGEX MATCH "\nversion[ \t]*=[ \t]*\"([^\"]+)\"" _ssplat_ver_match
           "${_ssplat_pyproject}")
    if(CMAKE_MATCH_1)
        set(SSPLAT_VERSION "${CMAKE_MATCH_1}")
    endif()
endif()

# ---------------------------------------------------------------------------
# Options
#
# Everything the build has goes into ONE executable, `ssplat`, which dispatches
# on its first argument (`ssplat sfm auto ...`); see src/app/Tools.h. These two
# options decide what is in it: the command-line tools, the window, or both.
# ---------------------------------------------------------------------------
option(SSPLAT_BUILD_CLI "Build the command-line tools into ssplat" OFF)
option(SSPLAT_BUILD_GUI "Build the graphical application into ssplat (fetches GLFW + Dear ImGui)" OFF)
option(SSPLAT_BUILD_BACKEND_TESTS "Build backend parity test tools" OFF)

# Also build ssplat-sfm and ssplat-sam as standalone executables. Same code and
# the same dispatcher (src/app/Main.cpp reads argv[0]), but built alone neither
# links the training engine or libtorch, which is the point: ssplat-sfm is
# 24 MB against the combined binary's 61 MB. The other tools are not offered
# separately -- they would be identical to `ssplat`, and a separate GUI could
# not run reconstruction, which is this binary re-running itself.
option(SSPLAT_SEPARATE_TOOLS "Also build ssplat-sfm / ssplat-sam standalone" OFF)

# Debug symbols / line info are OFF by default: they bloat the binaries
# massively (nvcc host -g, CUDA cubin lineinfo/source-in-ptx, and slangc -g2
# in the embedded SPIR-V). Turn on for profiling/debugging builds. The pip
# build (setup.py) is independent and unaffected (it has its own WITH_SYMBOLS
# and LINE_INFO env gates).
option(SSPLAT_DEBUG_SYMBOLS "Emit debug symbols / line info (host -g, CUDA cubin lineinfo, SPIR-V -g2)" OFF)

# ---------------------------------------------------------------------------
# Compute backend selection
#
# cuda (default): the full build -- CUDA kernels,
# and the app targets.
# vulkan: the portable engine layer (Engine*.cpp + host support) against the
# Vulkan compute runtime (src/backend/vulkan/, see its README.md), built
# WITHOUT the CUDA toolkit. Produces the backend tests and `ssplat` -- with the
# command-line tools always, and the GUI with SSPLAT_BUILD_GUI=ON; no Python
# extension.
# ---------------------------------------------------------------------------
set(SSPLAT_BACKEND "cuda" CACHE STRING "Compute backend: cuda | vulkan")
set_property(CACHE SSPLAT_BACKEND PROPERTY STRINGS cuda vulkan)

if(NOT SSPLAT_BACKEND STREQUAL "cuda" AND NOT SSPLAT_BACKEND STREQUAL "vulkan")
    message(FATAL_ERROR "SSPLAT_BACKEND must be 'cuda' or 'vulkan', got '${SSPLAT_BACKEND}'")
endif()

# ---------------------------------------------------------------------------
# Structure from Motion (src/sfm/)
#
# The native replacement for the COLMAP subprocess. It is Vulkan-only and needs
# the Vulkan loader + headers, so it is on by default only for the Vulkan build;
# a CUDA build can still opt in with -DSSPLAT_BUILD_SFM=ON if the Vulkan SDK is
# present. See cmake/SsplatSfm.cmake and src/sfm/README.md.
# ---------------------------------------------------------------------------
if(SSPLAT_BACKEND STREQUAL "vulkan")
    option(SSPLAT_BUILD_SFM "Build the SfM module (`ssplat sfm`)" ON)
else()
    option(SSPLAT_BUILD_SFM "Build the SfM module (`ssplat sfm`)" OFF)
endif()

# ---------------------------------------------------------------------------
# GPU inference (src/nn/) and segmentation (src/sam/)
#
# The native replacement for the scripts/mask.py subprocess: SAM 2 / SAM 3 on
# the same Vulkan + Slang stack as the SfM module, over a reusable inference
# layer. Same rule as SfM -- Vulkan-only, on by default only for the Vulkan
# build, opt-in for CUDA if the Vulkan SDK is present. See cmake/SsplatNn.cmake
# and src/nn/README.md.
# ---------------------------------------------------------------------------
if(SSPLAT_BACKEND STREQUAL "vulkan")
    option(SSPLAT_BUILD_SAM "Build the inference layer + SAM segmentation" ON)
else()
    option(SSPLAT_BUILD_SAM "Build the inference layer + SAM segmentation" OFF)
endif()

# ---------------------------------------------------------------------------
# Patent-encumbered modules
#
# OFF by default, and deliberately so: this repository is GPLv3, and the video
# codecs are the one piece of it that carries third-party patent exposure
# (H.264 / H.265 via MPEG LA / Access Advance, AV1 via the claims asserted
# against AOMedia). With it OFF, src/video/ -- the container demuxers, the
# H.264 / H.265 / AV1 bitstream parsers and the VK_KHR_video_decode_* driver --
# is neither compiled nor linked, and everything that wanted it shells out to
# an external ffmpeg instead. That costs a subprocess and a temporary folder of
# JPEGs; no feature disappears from the GUI.
#
# Turn it ON to decode video in-process on the GPU (roughly 15x faster frame
# extraction, and no ffmpeg to install). Distributors should check what their
# jurisdiction and their users require before shipping a binary built with it.
# ---------------------------------------------------------------------------
option(SSPLAT_ENABLE_PATENTED
    "Compile patent-encumbered modules (in-process H.264/H.265/AV1 video decode)"
    OFF)
if(SSPLAT_ENABLE_PATENTED AND NOT SSPLAT_BUILD_SAM)
    message(FATAL_ERROR
        "SSPLAT_ENABLE_PATENTED=ON needs SSPLAT_BUILD_SAM=ON: the video "
        "decoder is built on the inference layer's Vulkan runtime (src/nn/vk).")
endif()

# ---------------------------------------------------------------------------
# Shared C++ settings
# ---------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ---------------------------------------------------------------------------
# Host compiler flags
#
# SPLAT_CXX_FLAGS is what every target compiled from this tree gets for C++:
# the engine (csrc), the apps, and ssplat_sfm. It lives here, not in a backend
# module, because it is backend-neutral -- when it was set in
# SsplatBackendCuda.cmake the whole Vulkan build (engine *and* the SfM
# pipeline, which is host-heavy: stb decode, RANSAC, the mapper) compiled at
# -O0, since CMAKE_BUILD_TYPE is deliberately left empty.
#
# These intentionally mirror, but do not share, setup.py's flags: the pip build
# targets a released wheel (its own WITH_SYMBOLS / LINE_INFO env gates) while
# this one targets local development. Only the *source list* is shared, via
# cmake/sources.txt. Backend modules may append (OpenMP, torch defines).
# ---------------------------------------------------------------------------
if(MSVC)
    set(SPLAT_CXX_FLAGS "/O2")
else()
    set(SPLAT_CXX_FLAGS "-O3")
    list(APPEND SPLAT_CXX_FLAGS "-Wno-sign-compare")
endif()

if(WIN32)
    add_compile_definitions(_USE_MATH_DEFINES NOMINMAX _CRT_SECURE_NO_WARNINGS)
endif()
