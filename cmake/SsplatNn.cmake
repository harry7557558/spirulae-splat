# SSPLAT_BUILD_SAM: the GPU inference layer (src/nn/) and what currently sits
# on it -- SAM 2 / SAM 3 segmentation (src/sam/) and, when patented modules are
# enabled, hardware video decoding (src/video/).
#
# Defines:
#   ssplat_nn      the reusable inference layer (Vulkan runtime + tensor + ops)
#   ssplat_sam     SAM 2 / SAM 3 on top of it
#   ssplat_video   demux + VK_KHR_video_decode_*  (SSPLAT_ENABLE_PATENTED only)
#   nn_ops_test / sam_pipeline_test   one executable per src/*/tests/*.cpp
#
# The `ssplat-sam` CLI lives with the other app targets, in SsplatApps.cmake.
#
# Like cmake/SsplatSfm.cmake, this needs Vulkan and slangc but NOT the compute
# backend: src/nn/vk/ is its own Vulkan context, so the module builds
# identically under either SSPLAT_BACKEND.
#
# Shaders: one slangc edge per .slang module (all of a file's entry points land
# in one blob, which is what the pipeline cache's "<stem>.<entry>" key expects),
# then one generated TU per library that registers its blobs with
# nn/vk/EmbeddedSpirv.h's process registry. Three libraries, one pipeline cache.

find_package(Vulkan REQUIRED)
find_package(Threads REQUIRED)

include(SsplatSlang)
ssplat_find_slangc(SSPLAT_NN_SLANGC)
ssplat_build_spirv_tool(SSPLAT_NN_SPIRV_TOOL)

set(SSPLAT_NN_SHADER_ROOT ${SSPLAT_SRC}/nn/shaders)

# Serialize slangc.
#
# Several slangc processes at once crash it intermittently (0xC0000005 /
# 0xC0000409 on Windows, on a varying subset of the modules; a plain rerun
# succeeds). Whatever the shared state is, this shader set compiles in a few
# seconds, so a one-wide pool costs nothing and turns a flaky build into a
# reliable one. The SfM and engine shader sets have not hit it; revisit if a
# future slangc proves safe.
set_property(GLOBAL APPEND PROPERTY JOB_POOLS ssplat_nn_slangc=1)

# ---------------------------------------------------------------------------
# ssplat_nn_shaders(<tag> <shader_dir> <out_embed_cpp>)
#
# Compiles every .slang in <shader_dir> (files starting with '_' are
# include-only helpers) and emits the generated registration TU.
# ---------------------------------------------------------------------------
function(ssplat_nn_shaders tag shader_dir out_embed_cpp)
    set(_spv_dir ${CMAKE_CURRENT_BINARY_DIR}/nn_spirv/${tag})
    file(MAKE_DIRECTORY ${_spv_dir})

    file(GLOB _shaders CONFIGURE_DEPENDS ${shader_dir}/*.slang)
    if(NOT _shaders)
        message(FATAL_ERROR "no .slang files under ${shader_dir}")
    endif()
    # Helpers live in .slang files too and any module may include any of them,
    # so every blob depends on the whole shared set: a helper edit rebuilds all
    # of them. The set is small enough that parsing includes would not pay.
    file(GLOB _common CONFIGURE_DEPENDS ${SSPLAT_NN_SHADER_ROOT}/*.slang)

    set(_args -target spirv -O3 -fvk-use-entrypoint-name
              -matrix-layout-row-major -I${SSPLAT_NN_SHADER_ROOT})
    if(SSPLAT_DEBUG_SYMBOLS)
        list(APPEND _args -g2)
    endif()

    set(_blobs "")
    foreach(_src ${_shaders})
        get_filename_component(_stem ${_src} NAME_WE)
        if(_stem MATCHES "^_")
            continue()
        endif()
        set(_out ${_spv_dir}/${_stem}.spv)
        add_custom_command(OUTPUT ${_out}
            COMMAND ${SSPLAT_NN_SLANGC} ${_src} -o ${_out} ${_args}
            DEPENDS ${_shaders} ${_common}
            COMMENT "SPIR-V ${_stem}"
            JOB_POOL ssplat_nn_slangc
            VERBATIM)
        list(APPEND _blobs ${_out})
    endforeach()

    list(LENGTH _blobs _n)
    message(STATUS "Inference SPIR-V (${tag}): ${_n} modules")

    set(_listfile ${_spv_dir}/blobs.txt)
    string(REPLACE ";" "\n" _listbody "${_blobs}")
    file(WRITE ${_listfile} "${_listbody}\n")

    set(_embed ${CMAKE_CURRENT_BINARY_DIR}/nn_shaders_${tag}.cpp)
    add_custom_command(OUTPUT ${_embed}
        COMMAND ${SSPLAT_NN_SPIRV_TOOL} embed --nn ${tag} ${_embed}
                --list ${_listfile}
        DEPENDS ${SSPLAT_NN_SPIRV_TOOL} ${_blobs} ${_listfile}
        COMMENT "Embedding SPIR-V (${tag}: ${_n} modules)"
        VERBATIM)

    set(${out_embed_cpp} ${_embed} PARENT_SCOPE)
endfunction()

# ---------------------------------------------------------------------------
# ssplat_nn -- the reusable inference layer
# ---------------------------------------------------------------------------
ssplat_nn_shaders(nn ${SSPLAT_SRC}/nn/shaders SSPLAT_NN_EMBED)

file(GLOB_RECURSE SSPLAT_NN_SOURCES CONFIGURE_DEPENDS ${SSPLAT_SRC}/nn/*.cpp)
list(FILTER SSPLAT_NN_SOURCES EXCLUDE REGEX "/tests/")

add_library(ssplat_nn STATIC
    ${SSPLAT_NN_SOURCES}
    ${SSPLAT_NN_EMBED}
    # stb is instantiated once for the repository. As archive members these are
    # only pulled in when nothing else already provides them, so a target that
    # also links the engine (ssplat-gui) does not see a duplicate definition.
    ${SSPLAT_SRC}/external/stb_image_impl.cpp
    ${SSPLAT_SRC}/external/stb_image_write_impl.cpp
)
target_include_directories(ssplat_nn PUBLIC ${SSPLAT_SRC})
# PUBLIC and on ssplat_nn, not on ssplat_video: the flag means "this build
# contains the video decoder", and what reads it is nn/vk/Context.cpp, which
# has to request a video-decode queue at device creation -- before any decoder
# exists to ask for one.
if(SSPLAT_ENABLE_PATENTED)
    target_compile_definitions(ssplat_nn PUBLIC SSPLAT_HAVE_VIDEO=1)
endif()
target_link_libraries(ssplat_nn PUBLIC Vulkan::Vulkan Threads::Threads)
target_compile_options(ssplat_nn PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:${SPLAT_CXX_FLAGS}>)
set_property(TARGET ssplat_nn PROPERTY CXX_STANDARD 17)

# ---------------------------------------------------------------------------
# ssplat_sam -- SAM 2 / SAM 3
# ---------------------------------------------------------------------------
ssplat_nn_shaders(sam ${SSPLAT_SRC}/sam/shaders SSPLAT_SAM_EMBED)

file(GLOB_RECURSE SSPLAT_SAM_SOURCES CONFIGURE_DEPENDS ${SSPLAT_SRC}/sam/*.cpp)
list(FILTER SSPLAT_SAM_SOURCES EXCLUDE REGEX "/tests/")

add_library(ssplat_sam STATIC ${SSPLAT_SAM_SOURCES} ${SSPLAT_SAM_EMBED})
target_link_libraries(ssplat_sam PUBLIC ssplat_nn)
target_compile_options(ssplat_sam PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:${SPLAT_CXX_FLAGS}>)
set_property(TARGET ssplat_sam PROPERTY CXX_STANDARD 17)

# ---------------------------------------------------------------------------
# ssplat_aliked -- the ALIKED / LightGlue learned SfM frontend
#
# Its checkpoints are ONNX files fetched from COLMAP's releases and parsed in
# process (docs/notes/aliked-port-plan.md), so this needs no protobuf, no
# onnxruntime and no converter -- src/aliked/model/Onnx.cpp is a varint walk.
#
# The ops it needs that nn/ did not have -- deformable convolution, point-wise
# grid sample, average pooling, row L2 normalization -- are general and went
# into nn/shaders. What is here is only what could not be general: the
# detector's suppression rule and soft-argmax, and its coordinate conversions.
# ---------------------------------------------------------------------------
ssplat_nn_shaders(aliked ${SSPLAT_SRC}/aliked/shaders SSPLAT_ALIKED_EMBED)

file(GLOB_RECURSE SSPLAT_ALIKED_SOURCES CONFIGURE_DEPENDS ${SSPLAT_SRC}/aliked/*.cpp)
list(FILTER SSPLAT_ALIKED_SOURCES EXCLUDE REGEX "/tests/")

add_library(ssplat_aliked STATIC ${SSPLAT_ALIKED_SOURCES} ${SSPLAT_ALIKED_EMBED})
target_link_libraries(ssplat_aliked PUBLIC ssplat_nn)
target_compile_options(ssplat_aliked PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:${SPLAT_CXX_FLAGS}>)
set_property(TARGET ssplat_aliked PROPERTY CXX_STANDARD 17)

# ---------------------------------------------------------------------------
# ssplat_video -- container demux + VK_KHR_video_decode_*
#
# Gated on SSPLAT_ENABLE_PATENTED: H.264 / H.265 / AV1 bitstream parsing is the
# patent-encumbered code in this repository, and a default build leaves it out.
# Everything that uses it falls back to an external ffmpeg, so OFF costs a
# subprocess, not a feature.
# ---------------------------------------------------------------------------
if(SSPLAT_ENABLE_PATENTED)
    ssplat_nn_shaders(video ${SSPLAT_SRC}/video/shaders SSPLAT_VIDEO_EMBED)

    file(GLOB_RECURSE SSPLAT_VIDEO_SOURCES CONFIGURE_DEPENDS
         ${SSPLAT_SRC}/video/*.cpp)
    add_library(ssplat_video STATIC ${SSPLAT_VIDEO_SOURCES} ${SSPLAT_VIDEO_EMBED})
    target_link_libraries(ssplat_video PUBLIC ssplat_nn)
    target_compile_options(ssplat_video PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:${SPLAT_CXX_FLAGS}>)
    set_property(TARGET ssplat_video PROPERTY CXX_STANDARD 17)
endif()

# ---------------------------------------------------------------------------
# Tests -- one executable per file, as in src/backend/tests/
# ---------------------------------------------------------------------------
file(GLOB SSPLAT_NN_TESTS CONFIGURE_DEPENDS
     ${SSPLAT_SRC}/nn/tests/*.cpp ${SSPLAT_SRC}/sam/tests/*.cpp
     ${SSPLAT_SRC}/aliked/tests/*.cpp)
foreach(test_src ${SSPLAT_NN_TESTS})
    get_filename_component(test_name ${test_src} NAME_WE)
    add_executable(${test_name} ${test_src})
    # Every test links every library above it: the three are small, and one
    # rule here beats a per-directory list that drifts.
    target_link_libraries(${test_name} PRIVATE ssplat_sam ssplat_aliked)
    set_property(TARGET ${test_name} PROPERTY CXX_STANDARD 17)
    target_compile_options(${test_name} PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:${SPLAT_CXX_FLAGS}>)
endforeach()
