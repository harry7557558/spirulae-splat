#pragma once

// Minimal cross-platform subprocess runner for the GUI's external tools
// (colmap, ffmpeg). Streams merged stdout+stderr line-by-line to a callback
// and supports cooperative cancellation (the process is killed).

#include <atomic>
#include <functional>
#include <string>
#include <vector>

namespace gui {

// Exit codes returned by run_process in addition to the child's own.
constexpr int kSpawnFailed = -1;   // executable not found / spawn error
constexpr int kCancelled   = -2;   // cancel flag was set; process killed

// Runs argv[0] with argv[1..] as arguments, working directory `cwd` ("" =
// inherit). stderr is merged into stdout; each completed line is passed to
// on_line (without the newline). `cancel` is polled ~10x per second.
// Blocking -- call from a worker thread.
int run_process(const std::vector<std::string>& argv,
                const std::string& cwd,
                const std::function<void(const std::string&)>& on_line,
                const std::atomic<bool>& cancel);

// True when `exe` resolves to an executable (PATH search like the shell).
bool command_exists(const std::string& exe);

// Splits a free-form flag string the way a shell would for the simple cases:
// whitespace separated, with "quoted runs" kept together. Enough for
// `--max-error 2 --masks "/path/with space"`. This is what the runners' "extra
// arguments" fields go through before becoming argv entries.
std::vector<std::string> split_args(const std::string& s);

// Hands a URL to the desktop's default browser and returns immediately. False
// when nothing could be launched (a headless session, no xdg-open), which is
// the caller's cue to fall back to showing the address. Never blocks and never
// waits for the browser to exit.
bool open_url(const std::string& url);

}  // namespace gui
