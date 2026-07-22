# The native applications: ssplat-train (CLI), ssplat-mesh, ssplat-gui.
#
# Backend-agnostic: whichever backend module ran first left the engine to link
# in SSPLAT_APP_LIBS and set SSPLAT_WITH_TORCH. ssplat-mesh is the exception --
# the Vulkan backend stubs the meshing kernels, so it is CUDA-only.

include(SsplatEmbed)

# ---------------------------------------------------------------------------
# Source groups shared by several app targets
# ---------------------------------------------------------------------------
# Dataset parsing (COLMAP / Nerfstudio / Metashape) + miniz for the zipped
# Metashape .psx payloads.
set(SSPLAT_DATASET_SOURCES
    ${SSPLAT_SRC}/data/parsers/ColmapParser.cpp
    ${SSPLAT_SRC}/data/parsers/NerfstudioParser.cpp
    ${SSPLAT_SRC}/data/parsers/MetashapeParser.cpp
    ${SSPLAT_SRC}/data/parsers/DatasetCommon.cpp
    ${SSPLAT_SRC}/external/miniz.c
)

# The embedded web viewer: HTTP server + the render worker feeding it.
set(SSPLAT_VIEWER_SOURCES
    ${SSPLAT_SRC}/app/webviewer/HttpServer.cpp
    ${SSPLAT_SRC}/app/webviewer/RenderWorker.cpp
    ${SSPLAT_SRC}/app/webviewer/Viewer.cpp
)

# ---------------------------------------------------------------------------
# ssplat_configure_app(<target>)
#
# The settings every app target needs: include paths, the engine libraries,
# host flags, and the two link-time workarounds a Torch-enabled build requires.
# ---------------------------------------------------------------------------
function(ssplat_configure_app target)
    target_include_directories(${target} PRIVATE
        ${SSPLAT_SRC}
        ${CMAKE_BINARY_DIR}      # app_generated/{viewer_html,mask_py}.h
        ${CUDAToolkit_INCLUDE_DIRS}
    )
    target_link_libraries(${target} PRIVATE ${SSPLAT_APP_LIBS})
    target_compile_options(${target} PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:${SPLAT_CXX_FLAGS}>
    )
    set_property(TARGET ${target} PROPERTY CXX_STANDARD 17)

    if(SSPLAT_WITH_TORCH)
        # libpython is only needed because csrc's pybind/libtorch_python layer
        # leaves CPython symbols undefined (a Python process normally provides
        # them at import time). Drops out with the no-torch build.
        find_package(Python3 REQUIRED COMPONENTS Development.Embed)
        target_link_libraries(${target} PRIVATE Python3::Python)

        # Statically link libstdc++ into the exe so its C++ runtime calls bind
        # at link time. Without this, libtorch.so (pulled in via csrc)
        # interposes its own bundled std::filesystem symbols at load time, and
        # e.g. remove_all resolves to an ABI-incompatible copy that jumps to
        # null.
        if(NOT WIN32)
            target_link_options(${target} PRIVATE -static-libgcc -static-libstdc++)
        endif()

        # libcsrc.so is renamed to spirulae_splat/csrc.so by
        # build_develop.bash; keep the exe's lookup pointing at the build tree
        # copy (and at libtorch).
        set_property(TARGET ${target} PROPERTY BUILD_RPATH
            "$ORIGIN;${TORCH_INSTALL_PREFIX}/lib")
    else()
        set_property(TARGET ${target} PROPERTY BUILD_RPATH "$ORIGIN")
    endif()
endfunction()

# ---------------------------------------------------------------------------
# Embedded web viewer client
# ---------------------------------------------------------------------------
if(SSPLAT_BUILD_CLI OR SSPLAT_BUILD_GUI)
    # Embed the (unchanged) web viewer client into the binary so the exe is
    # self-contained. Regenerated at configure time when viewer.html changes.
    ssplat_embed_file(
        ${SSPLAT_ROOT}/spirulae_splat/viewer/viewer.html
        ${CMAKE_BINARY_DIR}/app_generated/viewer_html.h
        ViewerHtml)
endif()

# ---------------------------------------------------------------------------
# ssplat-train -- standalone CLI trainer, no Python at runtime. Off by default
# (unless Torch is unavailable, or the backend is Vulkan, which force it on).
# ---------------------------------------------------------------------------
if(SSPLAT_BUILD_CLI)
    add_executable(ssplat-train
        ${SSPLAT_SRC}/app/cli/main.cpp
        ${SSPLAT_SRC}/app/TrainerCore.cpp
        ${SSPLAT_DATASET_SOURCES}
        ${SSPLAT_VIEWER_SOURCES}
    )
    ssplat_configure_app(ssplat-train)
    if(WIN32)
        target_link_libraries(ssplat-train PRIVATE ws2_32)
    endif()

    # ---- standalone mesh-extraction CLI (shares the dataset parsers) ----
    # CUDA-only: the Vulkan backend stubs the meshing kernels
    # (meshing::OccupancyEvaluator), so the tool would throw at startup.
    if(NOT SSPLAT_BACKEND STREQUAL "vulkan")
        add_executable(ssplat-mesh
            ${SSPLAT_SRC}/app/cli/mesh_main.cpp
            ${SSPLAT_DATASET_SOURCES}
        )
        ssplat_configure_app(ssplat-mesh)
    endif()
endif()

# ---------------------------------------------------------------------------
# ssplat-gui -- optional, off by default. Dear ImGui + GLFW + OpenGL 3
# (imgui's embedded GL loader; no GLEW/GLAD). Dependencies are fetched at
# configure time via FetchContent and built statically, so the only system
# requirements are OpenGL dev headers (Linux: libgl1-mesa-dev + xorg-dev for
# GLFW; Windows/macOS: nothing extra).
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
    ssplat_embed_file(
        ${SSPLAT_ROOT}/scripts/mask.py
        ${CMAKE_BINARY_DIR}/app_generated/mask_py.h
        MaskPy)

    file(GLOB SSPLAT_GUI_SOURCES CONFIGURE_DEPENDS ${SSPLAT_SRC}/app/gui/*.cpp)
    add_executable(ssplat-gui
        ${SSPLAT_GUI_SOURCES}
        ${SSPLAT_SRC}/app/TrainerCore.cpp
        ${SSPLAT_DATASET_SOURCES}
        ${SSPLAT_VIEWER_SOURCES}
    )
    ssplat_configure_app(ssplat-gui)
    target_link_libraries(ssplat-gui PRIVATE imgui_glfw OpenGL::GL)
    if(WIN32)
        target_link_libraries(ssplat-gui PRIVATE ws2_32)
    endif()
endif()
