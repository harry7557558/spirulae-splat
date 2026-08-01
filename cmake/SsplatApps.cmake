# The native applications: ssplat-train (CLI), ssplat-mesh, ssplat-sfm,
# ssplat-gui.
#
# Backend-agnostic: whichever backend module ran first left the engine to link
# in SSPLAT_APP_LIBS and set SSPLAT_WITH_TORCH. Two exceptions: ssplat-mesh is
# CUDA-only (the Vulkan backend stubs the meshing kernels), and ssplat-sfm is
# Vulkan-only and links the SfM library instead of the engine.

# ---------------------------------------------------------------------------
# Source groups shared by several app targets
# ---------------------------------------------------------------------------
# Dataset parsing (COLMAP / Nerfstudio / Metashape) and miniz now live in the
# engine library (cmake/sources.txt), so every app target -- and the Python
# extension -- shares one build of them. Nothing to list here.

# The embedded web viewer (HTTP server + render worker) also lives in the
# engine library now -- the Python extension binds it. Nothing to list here.

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
# ssplat-train -- standalone CLI trainer, no Python at runtime. Off by default
# (unless Torch is unavailable, or the backend is Vulkan, which force it on).
# ---------------------------------------------------------------------------
if(SSPLAT_BUILD_CLI)
    add_executable(ssplat-train
        ${SSPLAT_SRC}/app/cli/main.cpp
    )
    ssplat_configure_app(ssplat-train)
    if(WIN32)
        target_link_libraries(ssplat-train PRIVATE ws2_32)
    endif()

    # ---- standalone mesh-extraction CLI ----
    # CUDA-only: the Vulkan backend stubs the meshing kernels
    # (meshing::OccupancyEvaluator), so the tool would throw at startup.
    if(NOT SSPLAT_BACKEND STREQUAL "vulkan")
        add_executable(ssplat-mesh ${SSPLAT_SRC}/app/cli/mesh_main.cpp)
        ssplat_configure_app(ssplat-mesh)
    endif()

    # ---- standalone Structure-from-Motion CLI ----
    # Links only ssplat_sfm: the SfM module carries its own Vulkan context and
    # SPIR-V, and shares nothing with the training engine. Deliberately NOT
    # ssplat_configure_app -- that pulls in the engine (and, on a Torch build,
    # libpython) for a tool that needs neither.
    if(SSPLAT_BUILD_SFM)
        add_executable(ssplat-sfm
            ${SSPLAT_SRC}/app/cli/sfm_main.cpp
            ${SSPLAT_SRC}/app/cli/sfm_ba.cpp
        )
        target_link_libraries(ssplat-sfm PRIVATE ssplat_sfm)
        target_compile_definitions(ssplat-sfm PRIVATE SSPLAT_VERSION="${SSPLAT_VERSION}")
        target_compile_options(ssplat-sfm PRIVATE
            $<$<COMPILE_LANGUAGE:CXX>:${SPLAT_CXX_FLAGS}>)
        set_property(TARGET ssplat-sfm PROPERTY CXX_STANDARD 17)
    endif()

    # ---- standalone segmentation / frame-extraction CLI ----
    # Links ssplat_sam (and ssplat_video when patented modules are enabled) for
    # the same reason ssplat-sfm links only ssplat_sfm: it needs the inference
    # layer, not the training engine.
    if(SSPLAT_BUILD_SAM)
        set(SSPLAT_SAM_CLI_SOURCES ${SSPLAT_SRC}/app/cli/sam_main.cpp)
        if(SSPLAT_ENABLE_PATENTED)
            list(APPEND SSPLAT_SAM_CLI_SOURCES
                 ${SSPLAT_SRC}/app/cli/sam_extract.cpp
                 ${SSPLAT_SRC}/app/FrameExtract.cpp)
        endif()
        add_executable(ssplat-sam ${SSPLAT_SAM_CLI_SOURCES})
        target_link_libraries(ssplat-sam PRIVATE ssplat_sam)
        if(SSPLAT_ENABLE_PATENTED)
            target_link_libraries(ssplat-sam PRIVATE ssplat_video)
        endif()
        target_compile_definitions(ssplat-sam PRIVATE SSPLAT_VERSION="${SSPLAT_VERSION}")
        target_compile_options(ssplat-sam PRIVATE
            $<$<COMPILE_LANGUAGE:CXX>:${SPLAT_CXX_FLAGS}>)
        set_property(TARGET ssplat-sam PROPERTY CXX_STANDARD 17)
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
    if(SSPLAT_ENABLE_PATENTED)
        list(APPEND SSPLAT_GUI_SOURCES ${SSPLAT_SRC}/app/FrameExtract.cpp)
    endif()
    add_executable(ssplat-gui
        ${SSPLAT_GUI_SOURCES}
    )
    ssplat_configure_app(ssplat-gui)
    target_link_libraries(ssplat-gui PRIVATE imgui_glfw OpenGL::GL)

    # In-process segmentation (interactive preview + dataset masking) and, when
    # patented modules are enabled, in-process video decoding. Both are
    # optional: DatasetPrep falls back to python + scripts/mask.py and to
    # ffmpeg, and the GUI hides what this build cannot do rather than failing
    # at run time.
    if(SSPLAT_BUILD_SAM)
        target_link_libraries(ssplat-gui PRIVATE ssplat_sam)
        target_compile_definitions(ssplat-gui PRIVATE SSPLAT_BUILD_SAM=1)
        if(SSPLAT_ENABLE_PATENTED)
            target_link_libraries(ssplat-gui PRIVATE ssplat_video)
        endif()
    endif()

    if(WIN32)
        target_link_libraries(ssplat-gui PRIVATE ws2_32)
    endif()
endif()
