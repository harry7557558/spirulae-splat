# SSPLAT_BACKEND=vulkan: the portable engine layer built against the Vulkan
# compute runtime, without the CUDA toolkit. Kernel coverage is tracked in
# src/backend/vulkan/README.md.
#
# Defines, for the shared app targets in SsplatApps.cmake:
#   SSPLAT_WITH_TORCH = OFF   (the Python extension is CUDA-only)
#   SSPLAT_APP_LIBS           (portable engine + Vulkan runtime + threads)

message(STATUS "SSPLAT_BACKEND=vulkan: portable engine layer + "
    "backend/vulkan runtime; kernel coverage is tracked in "
    "src/backend/vulkan/README.md.")

# The Python extension is CUDA-only, so the CLI trainer is this build's
# primary artifact (same forcing as the CUDA no-torch build).
set(SSPLAT_BUILD_CLI ON)

find_package(Vulkan REQUIRED)
find_package(Threads REQUIRED)

# ---------------------------------------------------------------------------
# Portable engine layer
# ---------------------------------------------------------------------------
# Compile-checked here; links against the Vulkan runtime once kernel launchers
# exist for the symbols it needs.
ssplat_collect_portable_sources(SSPLAT_PORTABLE_SOURCES)

add_library(csrc_portable OBJECT ${SSPLAT_PORTABLE_SOURCES})
target_compile_definitions(csrc_portable PUBLIC SSPLAT_BACKEND_VULKAN)
target_include_directories(csrc_portable PUBLIC ${SSPLAT_SRC} ${CMAKE_BINARY_DIR})

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(csrc_portable PUBLIC OpenMP::OpenMP_CXX)
endif()

# ---------------------------------------------------------------------------
# Slang -> SPIR-V (build tree), then embed
# ---------------------------------------------------------------------------
# Native (no Python): backend/vulkan/shaders/spirv_tool.cpp is compiled by the
# module below,
# which enumerates one slangc custom command per blob (so -j bounds RAM) and
# embeds the result. SPIR-V debug info (slangc -g2) is opt-in via
# SSPLAT_DEBUG_SYMBOLS; it otherwise dominates the embedded blob size.
include(SsplatSlang)
ssplat_find_slangc(SSPLAT_SLANGC_EXE)

set(SSPLAT_SPIRV_DIR ${CMAKE_CURRENT_BINARY_DIR}/spirv)
set(SSPLAT_SPIRV_EMBED ${CMAKE_CURRENT_BINARY_DIR}/vk_shaders_embedded.cpp)
include(${SSPLAT_VK_SHADERS}/SpirvShaders.cmake)
ssplat_setup_spirv(${SSPLAT_SHADERS} ${SSPLAT_VK_SHADERS} ${SSPLAT_SLANGC_EXE}
    ${SSPLAT_SPIRV_DIR} ${SSPLAT_SPIRV_EMBED} ${SSPLAT_DEBUG_SYMBOLS})

# ---------------------------------------------------------------------------
# Vulkan backend runtime (device context, memory, streams, events) + kernels
# ---------------------------------------------------------------------------
file(GLOB SSPLAT_VULKAN_SOURCES CONFIGURE_DEPENDS
    ${SSPLAT_SRC}/backend/vulkan/*.cpp
    ${SSPLAT_SRC}/backend/vulkan/kernels/*.cpp)
list(APPEND SSPLAT_VULKAN_SOURCES ${SSPLAT_SPIRV_EMBED})

add_library(ssplat_backend_vulkan STATIC ${SSPLAT_VULKAN_SOURCES})
target_compile_definitions(ssplat_backend_vulkan PUBLIC SSPLAT_BACKEND_VULKAN)
target_include_directories(ssplat_backend_vulkan PUBLIC ${SSPLAT_SRC})
target_link_libraries(ssplat_backend_vulkan PUBLIC Vulkan::Vulkan)

# ---------------------------------------------------------------------------
# Tests (run manually; see backend/vulkan/README.md)
# ---------------------------------------------------------------------------
# Vulkan-only unit tests.
file(GLOB SSPLAT_VULKAN_TESTS CONFIGURE_DEPENDS ${SSPLAT_SRC}/backend/vulkan/tests/*.cpp)
foreach(test_src ${SSPLAT_VULKAN_TESTS})
    get_filename_component(test_name ${test_src} NAME_WE)
    add_executable(${test_name} ${test_src})
    target_link_libraries(${test_name} PRIVATE ssplat_backend_vulkan)
endforeach()

# Cross-backend parity tools (the same sources build in the CUDA branch, where
# they dump the reference outputs these compare against).
file(GLOB SSPLAT_PARITY_TESTS CONFIGURE_DEPENDS ${SSPLAT_SRC}/backend/tests/*.cpp)
foreach(test_src ${SSPLAT_PARITY_TESTS})
    get_filename_component(test_name ${test_src} NAME_WE)
    add_executable(${test_name} ${test_src})
    target_link_libraries(${test_name} PRIVATE ssplat_backend_vulkan)
endforeach()

# Engine-level parity tools: drive the real engine (forward_3dgs + background
# + blit) through the portable layer. The backend lib provides the ported
# kernels plus throwing stubs (kernels/TrainingStubs*.cpp) for everything still
# CUDA-only.
file(GLOB SSPLAT_ENGINE_TESTS CONFIGURE_DEPENDS ${SSPLAT_SRC}/backend/tests/engine/*.cpp)
foreach(test_src ${SSPLAT_ENGINE_TESTS})
    get_filename_component(test_name ${test_src} NAME_WE)
    add_executable(${test_name} ${test_src})
    target_link_libraries(${test_name} PRIVATE csrc_portable
        ssplat_backend_vulkan Threads::Threads)
endforeach()

# Torch/Python extension is CUDA-only; the app targets build against the
# portable engine + Vulkan backend instead.
set(SSPLAT_WITH_TORCH OFF)
set(SSPLAT_APP_LIBS csrc_portable ssplat_backend_vulkan Threads::Threads)
