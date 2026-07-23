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

# ---------------------------------------------------------------------------
# Options
#
# Torch + Python are only needed for the Python extension module (ext.cpp).
# When either is missing -- or SSPLAT_NO_TORCH=ON -- the extension is skipped
# and only the standalone ssplat-train CLI is built (engine is torch-free).
# ---------------------------------------------------------------------------
option(SSPLAT_BUILD_CLI "Build the standalone ssplat-train CLI" OFF)
option(SSPLAT_BUILD_GUI "Build the native GUI app ssplat-gui (fetches GLFW + Dear ImGui)" OFF)
option(SSPLAT_NO_TORCH "Skip Torch/Python even if present; build only ssplat-train" OFF)
option(SSPLAT_BUILD_BACKEND_TESTS "Build backend parity test tools" OFF)

# Debug symbols / line info are OFF by default: they bloat the binaries
# massively (nvcc host -g, CUDA cubin lineinfo/source-in-ptx, and slangc -g2
# in the embedded SPIR-V). Turn on for profiling/debugging builds. The pip
# build (setup.py) is independent and unaffected (it has its own WITH_SYMBOLS
# and LINE_INFO env gates).
option(SSPLAT_DEBUG_SYMBOLS "Emit debug symbols / line info (host -g, CUDA cubin lineinfo, SPIR-V -g2)" OFF)

# ---------------------------------------------------------------------------
# Compute backend selection
#
# cuda (default): the full build -- CUDA kernels, optional Torch extension,
# and the app targets.
# vulkan: the portable engine layer (Engine*.cpp + host support) against the
# Vulkan compute runtime (src/backend/vulkan/, see its README.md), built
# WITHOUT the CUDA toolkit. Produces the backend tests, the ssplat-train CLI
# (always) and ssplat-gui (with SSPLAT_BUILD_GUI=ON); no Python extension.
# ---------------------------------------------------------------------------
set(SSPLAT_BACKEND "cuda" CACHE STRING "Compute backend: cuda | vulkan")
set_property(CACHE SSPLAT_BACKEND PROPERTY STRINGS cuda vulkan)

if(NOT SSPLAT_BACKEND STREQUAL "cuda" AND NOT SSPLAT_BACKEND STREQUAL "vulkan")
    message(FATAL_ERROR "SSPLAT_BACKEND must be 'cuda' or 'vulkan', got '${SSPLAT_BACKEND}'")
endif()

# ---------------------------------------------------------------------------
# Shared C++ settings
# ---------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(WIN32)
    add_compile_definitions(_USE_MATH_DEFINES NOMINMAX _CRT_SECURE_NO_WARNINGS)
endif()
