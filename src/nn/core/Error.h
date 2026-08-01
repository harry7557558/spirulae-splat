#pragma once
// One exception type for the whole library.
//
// Internally every failure throws; the public API (pipeline/Session.cpp,
// pipeline/Tracker.cpp) catches at its boundary and reports through the
// session's error string. Nothing is compiled out: a failure always produces a
// message that names what failed and, where the driver told us, why.

#include <stdexcept>
#include <string>

namespace nn {

class Error : public std::runtime_error {
public:
    explicit Error(const std::string& what) : std::runtime_error(what) {}
};

// Formats like printf and throws. Defined in core/Log.cpp.
[[noreturn]] void fail(const char* fmt, ...);

// Throws unless `cond`. Used for programmer errors and unmet device
// requirements alike -- both are unrecoverable at the call site.
#define NN_CHECK(cond, ...)                     \
    do {                                          \
        if (!(cond)) ::nn::fail(__VA_ARGS__);   \
    } while (0)

}  // namespace nn
