#pragma once

// The tools this program is made of, and the entry point of each.
//
// Everything ships as one executable (`ssplat`), which dispatches on its first
// argument -- `ssplat sfm auto ...`, `ssplat train ...`, `ssplat` on its own
// for the window. One file is easier to install, to put on a PATH and to copy
// onto another machine than five, and it is what lets the GUI run a
// reconstruction as a child process without depending on a sibling binary
// being present: it re-runs itself.
//
// Each entry point below is what used to be a `main`, unchanged apart from the
// name, and each still works entirely on its own -- nothing here initializes
// anything for anything else. A build with SSPLAT_BUILD_GUI off produces a
// pure command-line binary that needs no display and no GL.
//
// Which of these exist is decided at compile time by the SSPLAT_TOOL_* macros
// the build sets; see cmake/SsplatApps.cmake.

#include <string>

namespace app {

// How this tool was invoked, for its own usage and error messages: "ssplat
// sfm" as dispatched, or whatever a renamed copy or a symlink was called. Set
// once from argv[0] at the top of each entry point.
inline std::string& program_name() {
    static std::string name = "ssplat";
    return name;
}
inline void set_program_name(const char* argv0, const char* fallback) {
    program_name() = (argv0 && argv0[0]) ? argv0 : fallback;
}

// A block of help text with the program name it was written against swapped
// for the real one. The usage examples are worth keeping as readable literals
// rather than threading a %s through thirty of them.
inline std::string help_text(const char* text, const char* written_as) {
    std::string s = text;
    const std::string& prog = program_name();
    const std::string was = written_as;
    if (prog == was || was.empty()) return s;
    for (size_t p = s.find(was); p != std::string::npos;
         p = s.find(was, p + prog.size()))
        s.replace(p, was.size(), prog);
    return s;
}

// The subcommand name each tool answers to. Also matched against argv[0], so a
// symlink named ssplat-sfm behaves as `ssplat sfm` (which is how the separate
// executables of earlier releases keep working).
constexpr const char* kToolTrain = "train";
constexpr const char* kToolMesh  = "mesh";
constexpr const char* kToolSfm   = "sfm";
constexpr const char* kToolSam   = "sam";
constexpr const char* kToolGui   = "gui";

}  // namespace app

#ifdef SSPLAT_TOOL_TRAIN
int ssplat_train_main(int argc, char** argv);
#endif
#ifdef SSPLAT_TOOL_MESH
int ssplat_mesh_main(int argc, char** argv);
#endif
#ifdef SSPLAT_TOOL_SFM
int ssplat_sfm_main(int argc, char** argv);
#endif
#ifdef SSPLAT_TOOL_SAM
int ssplat_sam_main(int argc, char** argv);
#endif
#ifdef SSPLAT_TOOL_GUI
int ssplat_gui_main(int argc, char** argv);
#endif
