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
