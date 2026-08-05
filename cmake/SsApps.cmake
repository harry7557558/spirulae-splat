# The native application: one `spirula` executable holding every tool this
# build has (see src/app/Tools.h), dispatching on its first argument.
#
# One file instead of five is easier to install and to put on a PATH, and it is
# what lets the GUI run a reconstruction as a child process without a sibling
# binary having to be found next to it -- it re-runs itself. Nothing is shared
# between the tools at run time: each entry point is what used to be its own
# `main`, and a build with the GUI off produces a command-line binary that
# needs neither a display nor GL.
#
# SS_SEPARATE_TOOLS additionally builds the old one-executable-per-tool
# layout. Worth having for spirula-sfm in particular: on its own it does not
# link the engine, so it stays a small, portable binary.
#
# Backend-agnostic: whichever backend module ran first left the engine to link
# in SS_APP_LIBS. Two exceptions: the mesh tool is CUDA-only (the Vulkan
# backend stubs the meshing kernels), and the SfM tool is Vulkan-only.
#
# Dataset parsing, miniz and the embedded web viewer live in the engine
# library (cmake/sources.txt), so every app target shares one build of them --
# there are no source groups to list here.

# ---------------------------------------------------------------------------
# ss_configure_app(<target>)
#
# The settings every app target needs: include paths, the engine libraries,
# host flags.
# ---------------------------------------------------------------------------
function(ss_configure_app target)
    target_include_directories(${target} PRIVATE
        ${SS_SRC}
        ${CMAKE_BINARY_DIR}      # app_generated/{viewer_html,mask_py}.h
        ${CUDAToolkit_INCLUDE_DIRS}
    )
    target_link_libraries(${target} PRIVATE ${SS_APP_LIBS})
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
# once for the combined `spirula`, and once per tool if SS_SEPARATE_TOOLS is
# on.
# ---------------------------------------------------------------------------
set(SS_TOOL_SOURCES "")
set(SS_TOOL_DEFS "")
set(SS_TOOL_LIBS "")

if(SS_BUILD_CLI)
    # ---- trainer: no Python at runtime ----
    list(APPEND SS_TOOL_SOURCES ${SS_SRC}/app/cli/main.cpp)
    list(APPEND SS_TOOL_DEFS SS_TOOL_TRAIN=1)

    # ---- mesh extraction ----
    # CUDA-only: the Vulkan backend stubs the meshing kernels
    # (meshing::OccupancyEvaluator), so the tool would throw at startup.
    if(NOT SS_BACKEND STREQUAL "vulkan")
        list(APPEND SS_TOOL_SOURCES ${SS_SRC}/app/cli/mesh_main.cpp)
        list(APPEND SS_TOOL_DEFS SS_TOOL_MESH=1)
    endif()
endif()

if((SS_BUILD_CLI OR SS_BUILD_GUI) AND SS_BUILD_SFM)
    # ---- structure from motion ----
    # The SfM module carries its own Vulkan context and SPIR-V and shares
    # nothing with the training engine.
    list(APPEND SS_TOOL_SOURCES
         ${SS_SRC}/app/cli/sfm_main.cpp
         ${SS_SRC}/app/cli/sfm_ba.cpp)
    list(APPEND SS_TOOL_DEFS SS_TOOL_SFM=1)
    list(APPEND SS_TOOL_LIBS ss_sfm)
endif()

if((SS_BUILD_CLI OR SS_BUILD_GUI) AND SS_BUILD_SAM)
    # ---- segmentation / frame extraction ----
    list(APPEND SS_TOOL_SOURCES ${SS_SRC}/app/cli/sam_main.cpp)
    list(APPEND SS_TOOL_DEFS SS_TOOL_SAM=1)
    list(APPEND SS_TOOL_LIBS ss_sam)
    if(SS_ENABLE_PATENTED)
        list(APPEND SS_TOOL_SOURCES
             ${SS_SRC}/app/cli/sam_extract.cpp
             ${SS_SRC}/app/FrameExtract.cpp)
        list(APPEND SS_TOOL_LIBS ss_video)
    endif()
endif()

# ---------------------------------------------------------------------------
# spirula-gui -- optional, off by default. Dear ImGui + GLFW + OpenGL 3
# (imgui's embedded GL loader; no GLEW/GLAD). Dependencies are fetched at
# configure time via FetchContent and built statically, so the only system
# requirements are OpenGL dev headers (Linux: libgl1-mesa-dev + xorg-dev for
# GLFW; Windows/macOS: nothing extra).
# ---------------------------------------------------------------------------
if(SS_BUILD_GUI)
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
    ss_embed_file(
        ${SS_ROOT}/scripts/mask.py
        ${CMAKE_BINARY_DIR}/app_generated/mask_py.h
        MaskPy)

    file(GLOB SS_GUI_SOURCES CONFIGURE_DEPENDS ${SS_SRC}/app/gui/*.cpp)
    list(APPEND SS_TOOL_SOURCES ${SS_GUI_SOURCES})
    list(APPEND SS_TOOL_DEFS SS_TOOL_GUI=1)
    list(APPEND SS_TOOL_LIBS imgui_glfw OpenGL::GL)

    # In-process segmentation (interactive preview + dataset masking) and, when
    # patented modules are enabled, in-process video decoding. Both are
    # optional: DatasetPrep falls back to python + scripts/mask.py and to
    # ffmpeg, and the GUI hides what this build cannot do rather than failing
    # at run time.
    if(SS_BUILD_SAM)
        list(APPEND SS_TOOL_DEFS SS_BUILD_SAM=1)
        if(SS_ENABLE_PATENTED)
            list(APPEND SS_TOOL_SOURCES ${SS_SRC}/app/FrameExtract.cpp)
        endif()
    endif()
endif()

# ---------------------------------------------------------------------------
# spirula -- the executable
# ---------------------------------------------------------------------------
if(SS_BUILD_CLI OR SS_BUILD_GUI)
    # FrameExtract is claimed by both the segmentation tool and the GUI.
    list(REMOVE_DUPLICATES SS_TOOL_SOURCES)
    list(REMOVE_DUPLICATES SS_TOOL_LIBS)

    add_executable(spirula ${SS_SRC}/app/Main.cpp ${SS_TOOL_SOURCES})
    ss_configure_app(spirula)
    target_link_libraries(spirula PRIVATE ${SS_TOOL_LIBS})
    target_compile_definitions(spirula PRIVATE
        ${SS_TOOL_DEFS} SS_VERSION="${SS_VERSION}")
    if(WIN32)
        target_link_libraries(spirula PRIVATE ws2_32)
    endif()

    # ---- the standalone CLI tools, on request ----
    # Same dispatcher, one tool compiled into each: src/app/Main.cpp reads its
    # own argv[0], so `spirula-sfm auto ...` reaches the SfM tool with its
    # arguments untouched, --help and all.
    #
    # Only these two, and deliberately: built alone they skip
    # ss_configure_app(), so neither drags in the training engine, which
    # is the entire reason to want them separate -- spirula-sfm is 24 MB
    # against the combined binary's 61 MB. A separate
    # spirula-train or spirula-mesh would be byte-for-byte the work `spirula` does
    # anyway, and a separate spirula-gui would be worse than the combined one:
    # it could not run reconstruction, since that is this binary re-running
    # itself. Symlink `spirula` if you want those names.
    if(SS_SEPARATE_TOOLS)
        function(ss_tool_exe name sources defs libs)
            add_executable(${name} ${SS_SRC}/app/Main.cpp ${sources})
            target_include_directories(${name} PRIVATE ${SS_SRC} ${CMAKE_BINARY_DIR})
            target_link_libraries(${name} PRIVATE ${libs})
            target_compile_definitions(${name} PRIVATE
                ${defs} SS_VERSION="${SS_VERSION}")
            target_compile_options(${name} PRIVATE
                $<$<COMPILE_LANGUAGE:CXX>:${SPLAT_CXX_FLAGS}>)
            set_property(TARGET ${name} PROPERTY CXX_STANDARD 17)
        endfunction()

        if(SS_BUILD_SFM)
            ss_tool_exe(spirula-sfm
                "${SS_SRC}/app/cli/sfm_main.cpp;${SS_SRC}/app/cli/sfm_ba.cpp"
                "SS_TOOL_SFM=1" "ss_sfm")
        endif()
        if(SS_BUILD_SAM)
            set(_sam_src ${SS_SRC}/app/cli/sam_main.cpp)
            set(_sam_lib ss_sam)
            if(SS_ENABLE_PATENTED)
                list(APPEND _sam_src ${SS_SRC}/app/cli/sam_extract.cpp
                                     ${SS_SRC}/app/FrameExtract.cpp)
                list(APPEND _sam_lib ss_video)
            endif()
            ss_tool_exe(spirula-sam "${_sam_src}" "SS_TOOL_SAM=1" "${_sam_lib}")
        endif()
    endif()
endif()
