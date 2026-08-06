# Embedding data files into the executables as byte arrays, so the apps are
# self-contained (no runtime lookup of viewer.html or scripts/mask.py).

# ss_embed_file(<input> <output_header> <symbol>)
#
# Writes a header defining `k<symbol>[]` / `k<symbol>Size` holding the bytes of
# <input>. Regenerated at configure time whenever <input> changes.
function(ss_embed_file input output_header symbol)
    file(READ ${input} _hex HEX)
    string(REGEX REPLACE "([0-9a-f][0-9a-f])" "0x\\1," _bytes ${_hex})
    file(RELATIVE_PATH _rel ${SS_ROOT} ${input})
    file(WRITE ${output_header}
        "#pragma once\n"
        "// AUTO-GENERATED from ${_rel} -- do not edit.\n"
        "#include <cstddef>\n"
        "inline const unsigned char k${symbol}[] = {${_bytes}};\n"
        "inline const size_t k${symbol}Size = sizeof(k${symbol});\n")
    set_property(DIRECTORY ${SS_ROOT} APPEND
        PROPERTY CMAKE_CONFIGURE_DEPENDS ${input})
endfunction()

# ss_cjk_faces()
#
# Generates two headers from assets/fonts/cjk_faces.txt:
#
#   app_generated/cjk_faces.h    the table of downloadable FULL faces
#   app_generated/cjk_subsets.h  the four SUBSETS, embedded as byte arrays
#
# and -- for a regional build (SS_FONT_CJK=sc|tc|jp|kr|all) -- downloads the
# named full faces so the executable ships with them beside it.
#
# The subsets are embedded and the full faces are not, which is the whole
# design in one line. A subset is ~110 KB because it holds only the characters
# this program's own translations use, so all four fit in the executable and
# every language renders with nothing to download. A full face is 4-8 MB, and
# ss_embed_file() turns one byte into five characters of C source, so `all`
# would be a 130 MB array literal. A regional build therefore installs
# <exe dir>/fonts/, which is the first place Fonts.cpp looks.
function(ss_cjk_faces)
    set(spec ${SS_ROOT}/assets/fonts/cjk_faces.txt)
    set(header ${CMAKE_BINARY_DIR}/app_generated/cjk_faces.h)
    set(sub_header ${CMAKE_BINARY_DIR}/app_generated/cjk_subsets.h)
    # ENCODING UTF-8, or file(STRINGS) drops the non-ASCII bytes in the
    # comments and hands back the fragments around them as extra "lines".
    file(STRINGS ${spec} lines ENCODING UTF-8)

    set(body "")
    set(sub_arrays "")
    set(sub_table "")
    set(wanted "")
    if(SS_FONT_CJK STREQUAL "all")
        set(wanted sc tc jp kr)
    elseif(NOT SS_FONT_CJK MATCHES "^(fetch|none)$")
        set(wanted ${SS_FONT_CJK})
    endif()

    foreach(line ${lines})
        string(STRIP "${line}" line)
        if(line STREQUAL "" OR line MATCHES "^#")
            continue()
        endif()
        string(REPLACE " | " ";" f "${line}")
        list(LENGTH f n)
        if(NOT n EQUAL 6)
            message(FATAL_ERROR "${spec}: expected 6 ' | '-separated fields, got ${n}:\n  ${line}")
        endif()
        list(GET f 0 id)
        list(GET f 1 file)
        list(GET f 2 url)
        list(GET f 3 sha)
        list(GET f 4 bytes)
        list(GET f 5 label)
        string(APPEND body
            "    {\"${id}\", \"${file}\", \"${url}\",\n"
            "     \"${sha}\", ${bytes}ull, \"${label}\"},\n")

        # The committed subset that goes into the executable. Generated from
        # the catalogs by tools/make_ui_font.py; tools/check_font_coverage.py
        # fails the build when a translation outgrows it.
        string(TOUPPER ${id} ID)
        set(sub ${SS_ROOT}/assets/fonts/SpirulaCJK-${ID}.otf)
        if(NOT EXISTS ${sub})
            message(FATAL_ERROR
                "${sub} is missing.\n"
                "  It is a committed build artifact; regenerate the fonts with\n"
                "  python3 tools/make_ui_font.py")
        endif()
        file(READ ${sub} _sub_hex HEX)
        string(REGEX REPLACE "([0-9a-f][0-9a-f])" "0x\\1," _sub_bytes ${_sub_hex})
        string(APPEND sub_arrays
            "inline const unsigned char kCjkSubset${ID}[] = {${_sub_bytes}};\n")
        string(APPEND sub_table
            "    {\"${id}\", kCjkSubset${ID}, sizeof(kCjkSubset${ID})},\n")
        set_property(DIRECTORY ${SS_ROOT} APPEND
            PROPERTY CMAKE_CONFIGURE_DEPENDS ${sub})

        if(id IN_LIST wanted)
            set(dst ${CMAKE_BINARY_DIR}/fonts/${file})
            if(NOT EXISTS ${dst})
                message(STATUS "Downloading CJK font ${file} (${bytes} bytes)")
                file(DOWNLOAD ${url} ${dst}
                     EXPECTED_HASH SHA256=${sha}
                     SHOW_PROGRESS
                     STATUS dl_status)
                list(GET dl_status 0 dl_code)
                if(NOT dl_code EQUAL 0)
                    list(GET dl_status 1 dl_msg)
                    file(REMOVE ${dst})
                    message(FATAL_ERROR
                        "SS_FONT_CJK=${SS_FONT_CJK} needs ${file}, and the "
                        "download failed: ${dl_msg}\n"
                        "  Fetch it by hand into ${CMAKE_BINARY_DIR}/fonts/, "
                        "or build with SS_FONT_CJK=fetch.")
                endif()
            endif()
        endif()
    endforeach()

    if(body STREQUAL "")
        message(FATAL_ERROR "${spec} declared no faces")
    endif()

    file(RELATIVE_PATH rel ${SS_ROOT} ${spec})
    file(WRITE ${header}
        "#pragma once\n"
        "// AUTO-GENERATED from ${rel} -- do not edit.\n"
        "#include \"app/gui/Fonts.h\"\n"
        "\n"
        "namespace gui {\n"
        "inline constexpr CjkFace kCjkFaces[] = {\n"
        "${body}"
        "};\n"
        "}  // namespace gui\n")

    file(WRITE ${sub_header}
        "#pragma once\n"
        "// AUTO-GENERATED from assets/fonts/SpirulaCJK-*.otf -- do not edit.\n"
        "// Those files are themselves generated: tools/make_ui_font.py.\n"
        "#include \"app/gui/Fonts.h\"\n"
        "\n"
        "namespace gui {\n"
        "${sub_arrays}"
        "inline const CjkSubset kCjkSubsets[] = {\n"
        "${sub_table}"
        "};\n"
        "}  // namespace gui\n")

    set_property(DIRECTORY ${SS_ROOT} APPEND
        PROPERTY CMAKE_CONFIGURE_DEPENDS ${spec})
endfunction()
