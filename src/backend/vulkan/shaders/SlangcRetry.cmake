# Runs one slangc invocation, retrying it a few times before giving up.
#
# slangc dies intermittently when several copies run at once -- 0xC0000005, or
# "aborted due to internal error", on a *varying* subset of the blobs, and only
# on Windows so far. cmake/SsNn.cmake serializes its own (much smaller) shader
# set for the same crash. This set is 680 blobs and serializing it costs the
# better part of an hour, so retry instead -- after a wait that grows each
# time, since an immediate re-run sits in the same burst that killed the first.
#
# The command arrives as one string with '|' between arguments, because -D
# cannot carry a list through to a script.

string(REPLACE "|" ";" _cmd "${SS_CMD}")

# Seconds to wait before each retry; its length sets the retry count.
set(_backoff 2 5 10 20)
list(LENGTH _backoff _retries)
math(EXPR _attempts "${_retries} + 1")

foreach(_attempt RANGE ${_retries})
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
    if(_attempt LESS _retries)
        list(GET _backoff ${_attempt} _wait)
        message(STATUS "slangc failed (${_rv}) on ${SS_BLOB}, attempt "
                       "${_attempt}; retrying in ${_wait}s")
        execute_process(COMMAND ${CMAKE_COMMAND} -E sleep ${_wait})
    endif()
endforeach()

message(FATAL_ERROR
        "slangc failed ${_attempts} times on ${SS_BLOB}:\n${_err}${_out}")
