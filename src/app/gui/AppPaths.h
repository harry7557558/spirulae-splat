#pragma once

// Where the app keeps things, and where it finds its own siblings.
//
// Three questions with platform-specific answers that several GUI modules ask,
// answered once: the settings directory, the cache directory (model
// checkpoints, the COLMAP vocabulary tree), and the directory this executable
// lives in -- which is how ssplat-gui locates ssplat-sfm and ssplat-sam
// without depending on PATH.

#include <string>

namespace gui {

// Created on first call. Roaming config on Windows, $XDG_CONFIG_HOME on Linux.
std::string config_dir();

// Created on first call. LOCALAPPDATA on Windows, $XDG_CACHE_HOME on Linux.
// Large downloads live here, so it is deliberately not the config directory.
std::string cache_dir();

// The directory holding the running executable, or "" if it cannot be
// determined. Resolved once.
std::string exe_dir();

// Full path to a tool shipped alongside this executable ("ssplat-sfm" ->
// "<exe_dir>/ssplat-sfm[.exe]"), or "" when it is not there. Falls back to
// nothing rather than to PATH: a stale copy of ssplat-sfm found on PATH would
// be far more confusing than a clear "not found".
std::string sibling_tool(const std::string& name);

}  // namespace gui
