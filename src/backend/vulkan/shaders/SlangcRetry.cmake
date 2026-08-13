# Runs one slangc invocation, retrying it twice before giving up.
#
# slangc dies intermittently when several copies run at once -- 0xC0000005, or
# "aborted due to internal error", on a *varying* subset of the blobs, and only
# on Windows so far. cmake/SsNn.cmake serializes its own (much smaller) shader
# set for the same crash. This set is 680 blobs and serializing it costs the
# better part of an hour, so retry instead: the burst that provokes it has
# drained by the time the retry runs.
#
# The command arrives as one string with '|' between arguments, because -D
# cannot carry a list through to a script.

string(REPLACE "|" ";" _cmd "${SS_CMD}")

foreach(_attempt RANGE 2)
    execute_process(COMMAND ${_cmd} RESULT_VARIABLE _rv
                    OUTPUT_VARIABLE _out ERROR_VARIABLE _err)
    if(_rv EQUAL 0)
        if(_out)
            message("${_out}")
        endif()
        if(_err)
            message("${_err}")
        endif()
        return()
    endif()
    message(STATUS "slangc failed (${_rv}) on ${SS_BLOB}, attempt "
                   "${_attempt}; retrying")
endforeach()

message(FATAL_ERROR "slangc failed three times on ${SS_BLOB}:\n${_err}${_out}")
