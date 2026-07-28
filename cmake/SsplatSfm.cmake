# SSPLAT_BUILD_SFM: the Structure-from-Motion module (src/sfm/), a Vulkan-only
# subsystem that replaces the COLMAP subprocess. See src/sfm/README.md.
#
# Defines:
#   ssplat_sfm       the static library (SfM pipeline + embedded SPIR-V)
#   sfm_*_test       one executable per src/sfm/tests/*.cpp
#
# The `ssplat-sfm` CLI lives with the other app targets, in SsplatApps.cmake.
#
# This module needs Vulkan and slangc but NOT the compute backend: it carries
# its own Vulkan context (src/sfm/vk/) and its own SPIR-V blobs, so it builds
# identically under either SSPLAT_BACKEND. Converging it onto the engine's
# device is a later step (docs/notes/sfm-port-plan.md phase 6).

find_package(Vulkan REQUIRED)
find_package(Threads REQUIRED)

set(SSPLAT_SFM_SRC ${SSPLAT_SRC}/sfm)
set(SSPLAT_SFM_SHADERS ${SSPLAT_SFM_SRC}/shaders)

# ---------------------------------------------------------------------------
# Shader variant matrix
# ---------------------------------------------------------------------------
# The BA kernels are compiled as whole modules (all entry points in one blob,
# hence -fvk-use-entrypoint-name), once per (Real, Loss) configuration. SIFT and
# matching are single float-only blobs with no matrix, so a trimmed
# -DSSPLAT_SFM_REALS=df build still includes them.
#
# Trim the matrix while iterating:
#   cmake -DSSPLAT_SFM_REALS=df -DSSPLAT_SFM_LOSSES=trivial
# Both are cached, so a trimmed value sticks until you pass the full list again
# or wipe the build tree. The CLI errors out at runtime with "variant not built
# into this binary" for a combination that was trimmed away.
set(SSPLAT_SFM_REALS "float;double;df" CACHE STRING
    "SfM bundle-adjustment scalar configurations to compile")
set(SSPLAT_SFM_LOSSES "trivial;huber;cauchy" CACHE STRING
    "SfM bundle-adjustment robust losses to compile")

include(SsplatSlang)
ssplat_find_slangc(SSPLAT_SFM_SLANGC)
ssplat_build_spirv_tool(SSPLAT_SFM_SPIRV_TOOL)

set(_sfm_spirv_dir ${CMAKE_CURRENT_BINARY_DIR}/sfm_spirv)
set(_sfm_embed_cpp ${CMAKE_CURRENT_BINARY_DIR}/sfm_shaders_embedded.cpp)

# -I<shaders> lets a shader reach the shared device math by its shaders-relative
# path (`#include "common/camera.slang"`), matching the C++ convention.
set(_sfm_slang_args -target spirv -O2 -fvk-use-entrypoint-name
    -I${SSPLAT_SFM_SHADERS})

file(GLOB _sfm_ba_deps CONFIGURE_DEPENDS
    ${SSPLAT_SFM_SHADERS}/ba/*.slang ${SSPLAT_SFM_SHADERS}/common/*.slang)

set(_sfm_blobs "")
foreach(real ${SSPLAT_SFM_REALS})
    if(real STREQUAL "float")
        set(_rdef -DREAL_FLOAT)
    elseif(real STREQUAL "double")
        set(_rdef -DREAL_DOUBLE)
    elseif(real STREQUAL "df")
        set(_rdef -DREAL_DF)
    else()
        message(FATAL_ERROR "unknown SSPLAT_SFM_REALS entry '${real}' "
            "(expected float, double or df)")
    endif()
    foreach(loss ${SSPLAT_SFM_LOSSES})
        if(loss STREQUAL "trivial")
            set(_ldef "")
        elseif(loss STREQUAL "huber")
            set(_ldef -DLOSS_HUBER)
        elseif(loss STREQUAL "cauchy")
            set(_ldef -DLOSS_CAUCHY)
        else()
            message(FATAL_ERROR "unknown SSPLAT_SFM_LOSSES entry '${loss}' "
                "(expected trivial, huber or cauchy)")
        endif()

        set(_name ba_${real}_${loss})
        set(_out ${_sfm_spirv_dir}/${_name}.spv)
        if(real STREQUAL "df")
            # Compile to .raw.spv, then decorate. slangc emits no NoContraction
            # and some drivers then contract float expressions, destroying the
            # error-free transforms the emulated double-float type is built on.
            set(_raw ${_sfm_spirv_dir}/${_name}.raw.spv)
            add_custom_command(OUTPUT ${_out}
                COMMAND ${SSPLAT_SFM_SLANGC} ${SSPLAT_SFM_SHADERS}/ba/ba.slang
                        -o ${_raw} ${_sfm_slang_args} ${_rdef} ${_ldef}
                COMMAND ${SSPLAT_SFM_SPIRV_TOOL} nocontract ${_raw} ${_out}
                DEPENDS ${_sfm_ba_deps} ${SSPLAT_SFM_SPIRV_TOOL}
                COMMENT "SPIR-V ${_name} (+NoContraction)"
                VERBATIM)
        else()
            add_custom_command(OUTPUT ${_out}
                COMMAND ${SSPLAT_SFM_SLANGC} ${SSPLAT_SFM_SHADERS}/ba/ba.slang
                        -o ${_out} ${_sfm_slang_args} ${_rdef} ${_ldef}
                DEPENDS ${_sfm_ba_deps}
                COMMENT "SPIR-V ${_name}"
                VERBATIM)
        endif()
        list(APPEND _sfm_blobs ${_out})
    endforeach()
endforeach()

# Single-blob stages: <name> is the blob name the host looks up.
foreach(stage sift:sift/sift.slang match:match/bruteforce.slang)
    string(REPLACE ":" ";" _parts ${stage})
    list(GET _parts 0 _name)
    list(GET _parts 1 _rel)
    file(GLOB _deps CONFIGURE_DEPENDS
        ${SSPLAT_SFM_SHADERS}/${_name}/*.slang ${SSPLAT_SFM_SHADERS}/common/*.slang)
    set(_out ${_sfm_spirv_dir}/${_name}.spv)
    add_custom_command(OUTPUT ${_out}
        COMMAND ${SSPLAT_SFM_SLANGC} ${SSPLAT_SFM_SHADERS}/${_rel}
                -o ${_out} ${_sfm_slang_args}
        DEPENDS ${_deps}
        COMMENT "SPIR-V ${_name}"
        VERBATIM)
    list(APPEND _sfm_blobs ${_out})
endforeach()

list(LENGTH _sfm_blobs _sfm_nblobs)
message(STATUS "SfM SPIR-V: ${_sfm_nblobs} blobs "
    "(reals: ${SSPLAT_SFM_REALS}; losses: ${SSPLAT_SFM_LOSSES})")

# List file for the embed step (keeps the command line short).
set(_sfm_listfile ${_sfm_spirv_dir}/blobs.txt)
string(REPLACE ";" "\n" _sfm_listbody "${_sfm_blobs}")
file(WRITE ${_sfm_listfile} "${_sfm_listbody}\n")

add_custom_command(OUTPUT ${_sfm_embed_cpp}
    COMMAND ${SSPLAT_SFM_SPIRV_TOOL} embed --sfm ${_sfm_embed_cpp}
            --list ${_sfm_listfile}
    DEPENDS ${SSPLAT_SFM_SPIRV_TOOL} ${_sfm_blobs} ${_sfm_listfile}
    COMMENT "Embedding SfM SPIR-V (${_sfm_nblobs} blobs)"
    VERBATIM)

# ---------------------------------------------------------------------------
# Library
# ---------------------------------------------------------------------------
file(GLOB_RECURSE SSPLAT_SFM_SOURCES CONFIGURE_DEPENDS ${SSPLAT_SFM_SRC}/*.cpp)
list(FILTER SSPLAT_SFM_SOURCES EXCLUDE REGEX "/tests/")

add_library(ssplat_sfm STATIC
    ${SSPLAT_SFM_SOURCES}
    ${_sfm_embed_cpp}
    # stb_image is instantiated once for the repository. As an archive member it
    # is only pulled in when nothing else already provides it, so ssplat-gui --
    # which links the engine, and with it the same TU -- does not see a
    # duplicate definition.
    ${SSPLAT_SRC}/external/stb_image_impl.cpp
)
target_include_directories(ssplat_sfm PUBLIC ${SSPLAT_SRC})
target_link_libraries(ssplat_sfm PUBLIC Vulkan::Vulkan Threads::Threads)
target_compile_options(ssplat_sfm PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:${SPLAT_CXX_FLAGS}>)
set_property(TARGET ssplat_sfm PROPERTY CXX_STANDARD 17)

# Most of the pipeline is still header-only (the split into translation units is
# port plan phase 3). List the headers so IDEs and `ninja -t deps` see them.
file(GLOB_RECURSE SSPLAT_SFM_HEADERS CONFIGURE_DEPENDS ${SSPLAT_SFM_SRC}/*.h)
target_sources(ssplat_sfm PRIVATE ${SSPLAT_SFM_HEADERS})

# ---------------------------------------------------------------------------
# Tests -- one executable per file, as in src/backend/tests/
# ---------------------------------------------------------------------------
file(GLOB SSPLAT_SFM_TESTS CONFIGURE_DEPENDS ${SSPLAT_SFM_SRC}/tests/*.cpp)
foreach(test_src ${SSPLAT_SFM_TESTS})
    get_filename_component(test_name ${test_src} NAME_WE)
    add_executable(${test_name} ${test_src})
    target_link_libraries(${test_name} PRIVATE ssplat_sfm)
    set_property(TARGET ${test_name} PROPERTY CXX_STANDARD 17)
endforeach()
