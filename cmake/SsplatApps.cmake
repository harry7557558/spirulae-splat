# The native application: one `ssplat` executable holding every tool this
# build has (see src/app/Tools.h), dispatching on its first argument.
#
# One file instead of five is easier to install and to put on a PATH, and it is
# what lets the GUI run a reconstruction as a child process without a sibling
# binary having to be found next to it -- it re-runs itself. Nothing is shared
# between the tools at run time: each entry point is what used to be its own
# `main`, and a build with the GUI off produces a command-line binary that
# needs neither a display nor GL.
#
# SSPLAT_SEPARATE_TOOLS additionally builds the old one-executable-per-tool
# layout. Worth having for ssplat-sfm in particular: on its own it links
# neither the engine nor libtorch, so it stays a small, portable binary.
#
# Backend-agnostic: whichever backend module ran first left the engine to link
# in SSPLAT_APP_LIBS. Two exceptions: the mesh tool is
# CUDA-only (the Vulkan backend stubs the meshing kernels), and the SfM tool is
# Vulkan-only.

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
# host flags.
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

    set_property(TARGET ${target} PROPERTY BUILD_RPATH "$ORIGIN")
endfunction()

# ---------------------------------------------------------------------------
# Which tools this build has
#
# Each block appends the tool's sources, the macro that declares its entry
# point (src/app/Tools.h), and whatever it links. The lists are then used twice:
# once for the combined `ssplat`, and once per tool if SSPLAT_SEPARATE_TOOLS is
# on.
# ---------------------------------------------------------------------------
set(SSPLAT_TOOL_SOURCES "")
set(SSPLAT_TOOL_DEFS "")
set(SSPLAT_TOOL_LIBS "")

if(SSPLAT_BUILD_CLI)
    # ---- trainer: no Python at runtime ----
    list(APPEND SSPLAT_TOOL_SOURCES ${SSPLAT_SRC}/app/cli/main.cpp)
    list(APPEND SSPLAT_TOOL_DEFS SSPLAT_TOOL_TRAIN=1)

    # ---- mesh extraction ----
    # CUDA-only: the Vulkan backend stubs the meshing kernels
    # (meshing::OccupancyEvaluator), so the tool would throw at startup.
    if(NOT SSPLAT_BACKEND STREQUAL "vulkan")
        list(APPEND SSPLAT_TOOL_SOURCES ${SSPLAT_SRC}/app/cli/mesh_main.cpp)
        list(APPEND SSPLAT_TOOL_DEFS SSPLAT_TOOL_MESH=1)
    endif()
endif()

if((SSPLAT_BUILD_CLI OR SSPLAT_BUILD_GUI) AND SSPLAT_BUILD_SFM)
    # ---- structure from motion ----
    # The SfM module carries its own Vulkan context and SPIR-V and shares
    # nothing with the training engine.
    list(APPEND SSPLAT_TOOL_SOURCES
         ${SSPLAT_SRC}/app/cli/sfm_main.cpp
         ${SSPLAT_SRC}/app/cli/sfm_ba.cpp)
    list(APPEND SSPLAT_TOOL_DEFS SSPLAT_TOOL_SFM=1)
    list(APPEND SSPLAT_TOOL_LIBS ssplat_sfm)
endif()

if((SSPLAT_BUILD_CLI OR SSPLAT_BUILD_GUI) AND SSPLAT_BUILD_SAM)
    # ---- segmentation / frame extraction ----
    list(APPEND SSPLAT_TOOL_SOURCES ${SSPLAT_SRC}/app/cli/sam_main.cpp)
    list(APPEND SSPLAT_TOOL_DEFS SSPLAT_TOOL_SAM=1)
    list(APPEND SSPLAT_TOOL_LIBS ssplat_sam)
    if(SSPLAT_ENABLE_PATENTED)
        list(APPEND SSPLAT_TOOL_SOURCES
             ${SSPLAT_SRC}/app/cli/sam_extract.cpp
             ${SSPLAT_SRC}/app/FrameExtract.cpp)
        list(APPEND SSPLAT_TOOL_LIBS ssplat_video)
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
    list(APPEND SSPLAT_TOOL_SOURCES ${SSPLAT_GUI_SOURCES})
    list(APPEND SSPLAT_TOOL_DEFS SSPLAT_TOOL_GUI=1)
    list(APPEND SSPLAT_TOOL_LIBS imgui_glfw OpenGL::GL)

    # In-process segmentation (interactive preview + dataset masking) and, when
    # patented modules are enabled, in-process video decoding. Both are
    # optional: DatasetPrep falls back to python + scripts/mask.py and to
    # ffmpeg, and the GUI hides what this build cannot do rather than failing
    # at run time.
    if(SSPLAT_BUILD_SAM)
        list(APPEND SSPLAT_TOOL_DEFS SSPLAT_BUILD_SAM=1)
        if(SSPLAT_ENABLE_PATENTED)
            list(APPEND SSPLAT_TOOL_SOURCES ${SSPLAT_SRC}/app/FrameExtract.cpp)
        endif()
    endif()
endif()

# ---------------------------------------------------------------------------
# ssplat -- the executable
# ---------------------------------------------------------------------------
if(SSPLAT_BUILD_CLI OR SSPLAT_BUILD_GUI)
    # FrameExtract is claimed by both the segmentation tool and the GUI.
    list(REMOVE_DUPLICATES SSPLAT_TOOL_SOURCES)
    list(REMOVE_DUPLICATES SSPLAT_TOOL_LIBS)

    add_executable(ssplat ${SSPLAT_SRC}/app/Main.cpp ${SSPLAT_TOOL_SOURCES})
    ssplat_configure_app(ssplat)
    target_link_libraries(ssplat PRIVATE ${SSPLAT_TOOL_LIBS})
    target_compile_definitions(ssplat PRIVATE
        ${SSPLAT_TOOL_DEFS} SSPLAT_VERSION="${SSPLAT_VERSION}")
    if(WIN32)
        target_link_libraries(ssplat PRIVATE ws2_32)
    endif()

    # ---- the standalone CLI tools, on request ----
    # Same dispatcher, one tool compiled into each: src/app/Main.cpp reads its
    # own argv[0], so `ssplat-sfm auto ...` reaches the SfM tool with its
    # arguments untouched, --help and all.
    #
    # Only these two, and deliberately: built alone they skip
    # ssplat_configure_app(), so neither drags in the training engine, which
    # is the entire reason to want them separate -- ssplat-sfm is 24 MB
    # against the combined binary's 61 MB. A separate
    # ssplat-train or ssplat-mesh would be byte-for-byte the work `ssplat` does
    # anyway, and a separate ssplat-gui would be worse than the combined one:
    # it could not run reconstruction, since that is this binary re-running
    # itself. Symlink `ssplat` if you want those names.
    if(SSPLAT_SEPARATE_TOOLS)
        function(ssplat_tool_exe name sources defs libs)
            add_executable(${name} ${SSPLAT_SRC}/app/Main.cpp ${sources})
            target_include_directories(${name} PRIVATE ${SSPLAT_SRC} ${CMAKE_BINARY_DIR})
            target_link_libraries(${name} PRIVATE ${libs})
            target_compile_definitions(${name} PRIVATE
                ${defs} SSPLAT_VERSION="${SSPLAT_VERSION}")
            target_compile_options(${name} PRIVATE
                $<$<COMPILE_LANGUAGE:CXX>:${SPLAT_CXX_FLAGS}>)
            set_property(TARGET ${name} PROPERTY CXX_STANDARD 17)
        endfunction()

        if(SSPLAT_BUILD_SFM)
            ssplat_tool_exe(ssplat-sfm
                "${SSPLAT_SRC}/app/cli/sfm_main.cpp;${SSPLAT_SRC}/app/cli/sfm_ba.cpp"
                "SSPLAT_TOOL_SFM=1" "ssplat_sfm")
        endif()
        if(SSPLAT_BUILD_SAM)
            set(_sam_src ${SSPLAT_SRC}/app/cli/sam_main.cpp)
            set(_sam_lib ssplat_sam)
            if(SSPLAT_ENABLE_PATENTED)
                list(APPEND _sam_src ${SSPLAT_SRC}/app/cli/sam_extract.cpp
                                     ${SSPLAT_SRC}/app/FrameExtract.cpp)
                list(APPEND _sam_lib ssplat_video)
            endif()
            ssplat_tool_exe(ssplat-sam "${_sam_src}" "SSPLAT_TOOL_SAM=1" "${_sam_lib}")
        endif()
    endif()
endif()
