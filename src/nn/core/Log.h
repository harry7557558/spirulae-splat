#pragma once
// Diagnostics. Everything goes to stderr so a CLI run can pipe results to
// stdout cleanly (`ssam-cli ... 2>/dev/null > out.json`).
//
// The level is a runtime knob (SS_NN_LOG=0..3) with a compile-time ceiling
// (NN_LOG_MAX_LEVEL) so release builds can drop the calls entirely.
//
//   0 silent   1 warnings/errors   2 stage timing (default)   3 verbose

#include <cstdarg>
#include <cstdio>

#ifndef NN_LOG_MAX_LEVEL
#define NN_LOG_MAX_LEVEL 3
#endif

namespace nn {

int  log_level();            // cached read of $SS_NN_LOG, default 2
void set_log_level(int lvl);
void log_write(int level, const char* fmt, ...);

// A monotonic wall clock in milliseconds, for the stage timers.
double now_ms();

}  // namespace nn

#define NN_LOG_AT(level, ...)                              \
    do {                                                     \
        if ((level) <= NN_LOG_MAX_LEVEL)                   \
            ::nn::log_write((level), __VA_ARGS__);         \
    } while (0)

#define NN_LOG_ERROR(...) NN_LOG_AT(1, __VA_ARGS__)
#define NN_LOG_WARN(...)  NN_LOG_AT(1, __VA_ARGS__)
#define NN_LOG_INFO(...)  NN_LOG_AT(2, __VA_ARGS__)
#define NN_LOG_DEBUG(...) NN_LOG_AT(3, __VA_ARGS__)

// Scoped stage timer: logs "<name>: X.Y ms" at level 2 on destruction.
#define NN_SCOPED_TIMER(name) ::nn::ScopedTimer _nn_timer_##__LINE__(name)

namespace nn {

class ScopedTimer {
public:
    explicit ScopedTimer(const char* name) : name_(name), t0_(now_ms()) {}
    ~ScopedTimer() { NN_LOG_INFO("%s: %.1f ms\n", name_, now_ms() - t0_); }

private:
    const char* name_;
    double      t0_;
};

}  // namespace nn
