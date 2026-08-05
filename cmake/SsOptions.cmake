# Build options, backend selection, and the source-tree paths every other
# module builds on. Included first by the top-level CMakeLists.txt.

# ---------------------------------------------------------------------------
# Deprecated SSPLAT_* option names
#
# The prefix became SS_ when spirulae-splat became Spirula Studio. A cached
# SSPLAT_* value from an older build tree, or a scripted -DSSPLAT_*, still
# works and says so. Delete this block -- and the SSPLAT_ fallback in
# src/core/Env.h -- one release after the rename.
#
# SSPLAT_NO_TORCH / SSPLAT_WITH_TORCH deliberately get no alias: there is no
# Torch build left to select, and a warning about an unknown option is
# friendlier than silently ignoring one.
# ---------------------------------------------------------------------------
set(SS_LEGACY_OPTIONS
    BUILD_CLI BUILD_GUI BUILD_SFM BUILD_SAM BUILD_BACKEND_TESTS
    BACKEND SEPARATE_TOOLS DEBUG_SYMBOLS ENABLE_PATENTED
    SLANGC SFM_REALS SFM_LOSSES)
foreach(opt ${SS_LEGACY_OPTIONS})
    if(DEFINED SSPLAT_${opt} AND NOT DEFINED SS_${opt})
        message(DEPRECATION
            "SSPLAT_${opt} is deprecated; use SS_${opt}. "
            "The alias will be removed in the release after next.")
        set(SS_${opt} "${SSPLAT_${opt}}" CACHE STRING "" FORCE)
    endif()
endforeach()

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
# SS_ROOT is the repository root (the directory holding CMakeLists.txt).
# Everything else is expressed relative to it so the modules never care where
# they were included from.
# SS_SRC is also the include root: every local #include in the native
# tree is written relative to it (e.g. `#include "core/Common.cuh"`).
set(SS_ROOT    ${CMAKE_CURRENT_SOURCE_DIR})
set(SS_SRC     ${SS_ROOT}/src)
set(SS_SHADERS ${SS_SRC}/shaders)                      # shared Slang device math
set(SS_VK_SHADERS ${SS_SRC}/backend/vulkan/shaders)   # Vulkan-only entry points

# The version the apps report with --version. Declared here and nowhere else;
# it used to be read out of pyproject.toml, back when there was a package.
set(SS_VERSION "dev")

# ---------------------------------------------------------------------------
# Options
#
# Everything the build has goes into ONE executable, `spirula`, which dispatches
# on its first argument (`spirula sfm auto ...`); see src/app/Tools.h. These two
# options decide what is in it: the command-line tools, the window, or both.
# ---------------------------------------------------------------------------
option(SS_BUILD_CLI "Build the command-line tools into spirula" OFF)
option(SS_BUILD_GUI "Build the graphical application into spirula (fetches GLFW + Dear ImGui)" OFF)
option(SS_BUILD_BACKEND_TESTS "Build backend parity test tools" OFF)

# Also build spirula-sfm and spirula-sam as standalone executables. Same code
# and the same dispatcher (src/app/Main.cpp reads argv[0]), but built alone
# neither links the training engine, which is the point: spirula-sfm is 24 MB
# against the combined binary's 61 MB. The other tools are not offered
# separately -- they would be identical to `spirula`, and a separate GUI could
# not run reconstruction, which is this binary re-running itself.
option(SS_SEPARATE_TOOLS "Also build spirula-sfm / spirula-sam standalone" OFF)

# Debug symbols / line info are OFF by default: they bloat the binaries
# massively (nvcc host -g, CUDA cubin lineinfo/source-in-ptx, and slangc -g2
# in the embedded SPIR-V). Turn on for profiling/debugging builds.
option(SS_DEBUG_SYMBOLS "Emit debug symbols / line info (host -g, CUDA cubin lineinfo, SPIR-V -g2)" OFF)

# ---------------------------------------------------------------------------
# Compute backend selection
#
# cuda (default): the full build -- CUDA kernels,
# and the app targets.
# vulkan: the portable engine layer (Engine*.cpp + host support) against the
# Vulkan compute runtime (src/backend/vulkan/, see its README.md), built
# WITHOUT the CUDA toolkit. Produces the backend tests and `spirula` -- with
# the command-line tools always, and the GUI with SS_BUILD_GUI=ON.
# ---------------------------------------------------------------------------
set(SS_BACKEND "cuda" CACHE STRING "Compute backend: cuda | vulkan")
set_property(CACHE SS_BACKEND PROPERTY STRINGS cuda vulkan)

if(NOT SS_BACKEND STREQUAL "cuda" AND NOT SS_BACKEND STREQUAL "vulkan")
    message(FATAL_ERROR "SS_BACKEND must be 'cuda' or 'vulkan', got '${SS_BACKEND}'")
endif()

# ---------------------------------------------------------------------------
# Structure from Motion (src/sfm/)
#
# The native replacement for the COLMAP subprocess. It is Vulkan-only and needs
# the Vulkan loader + headers, so it is on by default only for the Vulkan build;
# a CUDA build can still opt in with -DSS_BUILD_SFM=ON if the Vulkan SDK is
# present. See cmake/SsSfm.cmake and src/sfm/README.md.
# ---------------------------------------------------------------------------
if(SS_BACKEND STREQUAL "vulkan")
    option(SS_BUILD_SFM "Build the SfM module (`spirula sfm`)" ON)
else()
    option(SS_BUILD_SFM "Build the SfM module (`spirula sfm`)" OFF)
endif()

# ---------------------------------------------------------------------------
# GPU inference (src/nn/) and segmentation (src/sam/)
#
# The native replacement for the scripts/mask.py subprocess: SAM 2 / SAM 3 on
# the same Vulkan + Slang stack as the SfM module, over a reusable inference
# layer. Same rule as SfM -- Vulkan-only, on by default only for the Vulkan
# build, opt-in for CUDA if the Vulkan SDK is present. See cmake/SsNn.cmake
# and src/nn/README.md.
# ---------------------------------------------------------------------------
if(SS_BACKEND STREQUAL "vulkan")
    option(SS_BUILD_SAM "Build the inference layer + SAM segmentation" ON)
else()
    option(SS_BUILD_SAM "Build the inference layer + SAM segmentation" OFF)
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
option(SS_ENABLE_PATENTED
    "Compile patent-encumbered modules (in-process H.264/H.265/AV1 video decode)"
    OFF)
if(SS_ENABLE_PATENTED AND NOT SS_BUILD_SAM)
    message(FATAL_ERROR
        "SS_ENABLE_PATENTED=ON needs SS_BUILD_SAM=ON: the video "
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
# the engine (csrc), the apps, and ss_sfm. It lives here, not in a backend
# module, because it is backend-neutral -- when it was set in
# SsBackendCuda.cmake the whole Vulkan build (engine *and* the SfM
# pipeline, which is host-heavy: stb decode, RANSAC, the mapper) compiled at
# -O0, since CMAKE_BUILD_TYPE is deliberately left empty.
#
# Backend modules may append to it (OpenMP, for one).
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
