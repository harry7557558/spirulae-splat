# Collects the engine/kernel sources from cmake/sources.txt.
#
# The glob patterns live in a plain text file rather than here so that setup.py
# (which builds the same sources through torch.utils.cpp_extension, without
# CMake) reads the exact same list. See cmake/sources.txt.

# _ssplat_expand_section(<section> <out_var>)
#
# CONFIGURE_DEPENDS: kernel families are routinely split into several
# translation units, each named after what it does (see docs/codegen.md), and
# codegen adds/removes files under src/instantiations/. Without it a newly
# added source is silently not compiled and only surfaces as an undefined
# reference at link time.
function(_ssplat_expand_section section out_var)
    set(spec ${SSPLAT_ROOT}/cmake/sources.txt)
    file(STRINGS ${spec} lines)

    set(patterns "")
    set(in_section FALSE)
    foreach(line ${lines})
        string(STRIP "${line}" line)
        if(line STREQUAL "" OR line MATCHES "^#")
            continue()
        endif()
        if(line MATCHES "^\\[(.*)\\]$")
            if(CMAKE_MATCH_1 STREQUAL section)
                set(in_section TRUE)
            else()
                set(in_section FALSE)
            endif()
            continue()
        endif()
        if(in_section)
            list(APPEND patterns ${SSPLAT_ROOT}/${line})
        endif()
    endforeach()

    if(NOT patterns)
        message(FATAL_ERROR "No [${section}] patterns found in ${spec}")
    endif()

    file(GLOB sources CONFIGURE_DEPENDS ${patterns})

    # Drop ROCm sources (neither build compiles them).
    list(FILTER sources EXCLUDE REGEX "hip")

    if(NOT sources)
        message(FATAL_ERROR "cmake/sources.txt [${section}] matched no files")
    endif()

    # Reconfigure when the spec itself changes.
    set_property(DIRECTORY ${SSPLAT_ROOT} APPEND
        PROPERTY CMAKE_CONFIGURE_DEPENDS ${spec})

    set(${out_var} ${sources} PARENT_SCOPE)
endfunction()

# ssplat_collect_sources(<out_var>) -- everything the CUDA `csrc` lib compiles.
function(ssplat_collect_sources out_var)
    _ssplat_expand_section(cuda result)
    set(${out_var} ${result} PARENT_SCOPE)
endfunction()

# ssplat_collect_portable_sources(<out_var>) -- the torch-free, CUDA-free
# subset compiled as csrc_portable by the Vulkan build.
function(ssplat_collect_portable_sources out_var)
    _ssplat_expand_section(portable result)
    set(${out_var} ${result} PARENT_SCOPE)
endfunction()
