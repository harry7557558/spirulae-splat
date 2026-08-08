#pragma once

// The one way this program reads an environment variable.
//
// Every knob is SS_<suffix>. The old SSPLAT_<suffix> spelling still works and
// warns once per run; delete the fallback below together with the CMake option
// aliases in cmake/SsOptions.cmake.
//
//   if (const char* d = spirula::env("VK_DEVICE")) ...
//
// Header-only on purpose: the standalone tool binaries link neither the engine
// nor each other, and this is small enough not to need a home.

#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <string>

namespace spirula {

// Value of SS_<suffix>, or of the deprecated SSPLAT_<suffix>, or nullptr.
inline const char* env(const char* suffix) {
    std::string name = "SS_";
    name += suffix;
    if (const char* v = std::getenv(name.c_str())) return v;

    name = "SSPLAT_";
    name += suffix;
    const char* v = std::getenv(name.c_str());
    if (v) {
        // One shot for the whole process: enough to prompt a rename, and safe
        // to call from the worker threads that read these knobs.
        static std::atomic<bool> warned{false};
        if (!warned.exchange(true))
            // English, and one of the few things that is. This header is
            // deliberately dependency-free -- the standalone tool binaries
            // include it and link neither the engine nor the i18n library --
            // and the sentence is about the spelling of an environment
            // variable, read by whoever typed one.
            std::fprintf(stderr,
                         "warning: %s is deprecated; use SS_%s "
                         "(SSPLAT_* is reported once per run)\n",
                         name.c_str(), suffix);
    }
    return v;
}

// Set to anything other than "0" -- the shape most of the knobs want.
inline bool env_on(const char* suffix) {
    const char* v = env(suffix);
    return v && v[0] && !(v[0] == '0' && v[1] == '\0');
}

}  // namespace spirula
