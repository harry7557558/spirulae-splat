#pragma once

// Where the app keeps things, and where it finds its own siblings.
//
// Three questions with platform-specific answers that several GUI modules ask,
// answered once: the settings directory, the cache directory (model
// checkpoints, the COLMAP vocabulary tree), and where this executable is --
// which is how the GUI runs a reconstruction, since that is `spirula sfm`, i.e.
// this same binary again.

#include <string>

namespace gui {

// Created on first call. Roaming config on Windows, $XDG_CONFIG_HOME on Linux.
std::string config_dir();

// Created on first call. LOCALAPPDATA on Windows, $XDG_CACHE_HOME on Linux.
// Large downloads live here, so it is deliberately not the config directory.
std::string cache_dir();

// The running executable, and the directory holding it. Both "" if they cannot
// be determined. Resolved once.
//
// exe_path() is what the GUI spawns to run a tool -- `<exe_path> sfm auto ...`
// -- rather than looking one up by name. Nothing to find, nothing to install,
// and no chance of reaching a different build that happens to be on PATH.
std::string exe_path();
std::string exe_dir();

}  // namespace gui
