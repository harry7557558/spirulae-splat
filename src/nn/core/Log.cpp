#include "nn/core/Log.h"

#include "nn/core/Error.h"

#include <chrono>
#include <cstdlib>
#include <cstring>
#include "core/Env.h"

namespace nn {

static int g_level = -1;

int log_level() {
    if (g_level < 0) {
        const char* env = spirula::env("NN_LOG");
        g_level = env ? std::atoi(env) : 2;
    }
    return g_level;
}

void set_log_level(int lvl) { g_level = lvl; }

void log_write(int level, const char* fmt, ...) {
    if (level > log_level()) return;
    std::va_list ap;
    va_start(ap, fmt);
    std::vfprintf(stderr, fmt, ap);
    va_end(ap);
}

double now_ms() {
    using clock = std::chrono::steady_clock;
    static const clock::time_point t0 = clock::now();
    return std::chrono::duration<double, std::milli>(clock::now() - t0).count();
}

void fail(const char* fmt, ...) {
    char buf[1024];
    std::va_list ap;
    va_start(ap, fmt);
    std::vsnprintf(buf, sizeof buf, fmt, ap);
    va_end(ap);
    throw Error(buf);
}

}  // namespace nn
