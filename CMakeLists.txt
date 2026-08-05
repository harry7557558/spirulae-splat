cmake_minimum_required(VERSION 3.18 FATAL_ERROR)
project(spirulae_splat LANGUAGES C CXX)  # CUDA enabled by the cuda backend module

# Development build (see build_develop.bash / build_develop.bat and
# docs/build.md). The pip build is separate -- setup.py -- and shares only the
# source list, cmake/sources.txt.
#
# The logic lives in cmake/; this file is just the running order:
#
#   SsplatOptions.cmake        options, backend selection, tree paths
#   SsplatSources.cmake        the engine/kernel source list
#   SsplatBackendCuda.cmake    CUDA kernels + the engine library         } one
#   SsplatBackendVulkan.cmake  portable engine + Vulkan compute runtime } of
#   SsplatSfm.cmake            the SfM module (optional; Vulkan-only)
#   SsplatNn.cmake             inference layer + SAM + video (optional; Vulkan-only)
#   SsplatApps.cmake           the `ssplat` executable (every tool in one)
#
# Each backend module leaves SSPLAT_APP_LIBS behind for
# the shared app targets.

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

include(SsplatOptions)
include(SsplatSources)

# The web viewer's browser client is compiled into the engine library
# (app/webviewer/Viewer.cpp), so it must be embedded before either backend
# module defines its targets -- not alongside the app targets.
include(SsplatEmbed)
ssplat_embed_file(
    ${SSPLAT_SRC}/app/webviewer/viewer.html
    ${CMAKE_BINARY_DIR}/app_generated/viewer_html.h
    ViewerHtml)

if(SSPLAT_BACKEND STREQUAL "vulkan")
    include(SsplatBackendVulkan)
else()
    include(SsplatBackendCuda)
endif()

# Before SsplatApps: the app links the ssplat_sfm target.
if(SSPLAT_BUILD_SFM)
    include(SsplatSfm)
endif()

# Likewise for `ssplat sam` and the GUI's built-in masking.
if(SSPLAT_BUILD_SAM)
    include(SsplatNn)
endif()

include(SsplatApps)

message(STATUS "spirulae_splat configured (backend: ${SSPLAT_BACKEND}).")
