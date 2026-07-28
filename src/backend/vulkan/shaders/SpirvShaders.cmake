# Native Slang -> SPIR-V build, replacing the former Python scripts
# (slang/build_spirv.py + spirulae_splat/embed_spirv.py). No Python required.
#
# spirv_tool.cpp (next to this file) is compiled once (via try_compile, using the project
# toolchain, so it is portable across Windows/Linux/macOS) and used to:
#   1. discover -- scan the compute shaders and enumerate every SPIR-V blob;
#   2. embed    -- gather the compiled .spv blobs into one C++ TU.
#
# Each blob becomes its own slangc custom command (one build-graph edge), so the
# build's -j governs how many slangc run concurrently -- bounding peak RAM the
# same way source compiles are -- and Ninja prints "[n/m] SPIR-V <name>" as each
# finishes.
#
# ssplat_setup_spirv(<shared_dir> <vk_dir> <slangc> <spirv_dir> <embed_cpp> <debug>)
# defines the custom commands producing <embed_cpp>.
#   <shared_dir> = src/shaders            device math shared with the CUDA build
#   <vk_dir>     = src/backend/vulkan/shaders   Vulkan-only compute entry points
function(ssplat_setup_spirv shared_dir vk_dir slangc spirv_dir embed_cpp debug)
    # The host tool is shared with the SfM module (cmake/SsplatSfm.cmake), so
    # it is built by SsplatSlang.cmake rather than here.
    ssplat_build_spirv_tool(tool_exe)

    # -I<src> lets a Vulkan shader reach the shared device math by its
    # src-relative path (`#include "shaders/densify.slang"`), matching the C++
    # convention and disambiguating the names that exist in both directories.
    set(incdirs -I${SSPLAT_SRC} -I${shared_dir} -I${vk_dir})

    # Common slangc arguments. Debug info (-g2) is opt-in; the release path
    # strips line directives and slang's capability diagnostics (the embed
    # step's capability audit still runs).
    set(spirv_args -target spirv -stage compute -O2)
    if(debug)
        list(APPEND spirv_args -g2)
    else()
        list(APPEND spirv_args -ignore-capabilities -line-directive-mode none)
    endif()

    # Discover every blob (name / source / entry / defines / include-closure).
    file(GLOB slang_sources CONFIGURE_DEPENDS ${vk_dir}/*.slang)
    execute_process(
        COMMAND ${tool_exe} discover ${incdirs} ${slang_sources}
        OUTPUT_VARIABLE manifest
        RESULT_VARIABLE disc_rv
        ERROR_VARIABLE disc_err)
    if(NOT disc_rv EQUAL 0)
        message(FATAL_ERROR "spirv_tool discover failed (${disc_rv}):\n${disc_err}")
    endif()
    if(disc_err)
        message(STATUS "spirv_tool discover: ${disc_err}")
    endif()

    string(REPLACE "\n" ";" manifest_lines "${manifest}")
    set(blob_outputs "")
    foreach(line ${manifest_lines})
        if(line STREQUAL "")
            continue()
        endif()
        # fields: name <tab> source <tab> entry <tab> defines <tab> deps
        string(REPLACE "\t" ";" fields "${line}")
        list(GET fields 0 name)
        list(GET fields 1 source)
        list(GET fields 2 entry)
        list(GET fields 3 defines)
        list(GET fields 4 deps)

        set(define_args "")
        if(NOT defines STREQUAL "-")
            string(REPLACE " " ";" define_args "${defines}")
        endif()
        string(REPLACE " " ";" dep_list "${deps}")

        set(out ${spirv_dir}/${name}.spv)
        add_custom_command(
            OUTPUT ${out}
            COMMAND ${slangc} ${source} -entry ${entry} -o ${out}
                ${incdirs} ${spirv_args} ${define_args}
            DEPENDS ${dep_list}
            COMMENT "SPIR-V ${name}"
            VERBATIM)
        list(APPEND blob_outputs ${out})
    endforeach()

    list(LENGTH blob_outputs blob_count)
    message(STATUS "SPIR-V: ${blob_count} blobs to compile (per-blob, -j aware)")

    # List file for the embed step (avoids a very long command line).
    set(listfile ${spirv_dir}/blobs.txt)
    string(REPLACE ";" "\n" listfile_body "${blob_outputs}")
    file(WRITE ${listfile} "${listfile_body}\n")

    add_custom_command(
        OUTPUT ${embed_cpp}
        COMMAND ${tool_exe} embed ${embed_cpp} --list ${listfile}
        DEPENDS ${tool_exe} ${blob_outputs}
        COMMENT "Embedding SPIR-V (${blob_count} blobs)"
        VERBATIM)
endfunction()
