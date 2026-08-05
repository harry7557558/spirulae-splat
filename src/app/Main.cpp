// The one entry point -- see app/Tools.h for why there is only one.
//
//   spirula                    the window
//   spirula <file-or-folder>   the window, opening what was named
//   spirula sfm auto ...       structure from motion
//   spirula train ...          the trainer
//   spirula sam segment ...    segmentation
//   spirula mesh ...           mesh extraction
//
// A first argument that is not a subcommand goes to the GUI untouched, so
// "Open with" from a file manager and a shell alias both land on the right
// screen. On a build with no GUI it is an error naming the subcommands, which
// is the only thing such a build can do.

#include "app/Tools.h"

#include <cctype>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace {

struct Tool {
    const char* name;
    const char* summary;
    int (*run)(int, char**);
};

const std::vector<Tool>& tools() {
    static const std::vector<Tool> kTools = {
#ifdef SS_TOOL_GUI
        {app::kToolGui, "the graphical application (the default)", spirula_gui_main},
#endif
#ifdef SS_TOOL_SFM
        {app::kToolSfm, "structure from motion: photos or frames -> cameras",
         spirula_sfm_main},
#endif
#ifdef SS_TOOL_TRAIN
        {app::kToolTrain, "train a splat model on a dataset", spirula_train_main},
#endif
#ifdef SS_TOOL_SAM
        {app::kToolSam, "segmentation, tracking and frame extraction",
         spirula_sam_main},
#endif
#ifdef SS_TOOL_MESH
        {app::kToolMesh, "extract a mesh from a trained model", spirula_mesh_main},
#endif
    };
    return kTools;
}

const Tool* find_tool(const char* name) {
    if (!name) return nullptr;
    for (const Tool& t : tools())
        if (std::strcmp(t.name, name) == 0) return &t;
    return nullptr;
}

// argv[0] as a tool name: a copy or symlink called spirula-sfm runs the SfM
// tool, which is how the separately-named executables of earlier releases keep
// working. Only the basename matters, and only the part after the prefix.
// "ssplat-" is the pre-rename spelling and still resolves.
const Tool* tool_from_argv0(const char* argv0) {
    if (!argv0) return nullptr;
    std::string s = argv0;
    const size_t slash = s.find_last_of("/\\");
    if (slash != std::string::npos) s = s.substr(slash + 1);
    for (char& c : s) c = (char)std::tolower((unsigned char)c);
    if (s.size() > 4 && s.compare(s.size() - 4, 4, ".exe") == 0)
        s.resize(s.size() - 4);
    for (const char* prefix : {"spirula-", "ssplat-"})
        if (s.rfind(prefix, 0) == 0)
            return find_tool(s.c_str() + std::strlen(prefix));
    return nullptr;
}

void print_usage() {
    std::printf("Spirula Studio " SS_VERSION "\n\n");
#ifdef SS_TOOL_GUI
    std::printf("  spirula                        open the application\n");
    std::printf("  spirula <file-or-folder>       open it in the application\n");
#endif
    std::printf("  spirula <command> [options]\n\n");
    std::printf("Commands:\n");
    for (const Tool& t : tools())
        std::printf("  %-8s %s\n", t.name, t.summary);
    std::printf("\n`spirula <command> --help` describes one of them.\n");
}

}  // namespace

int main(int argc, char** argv) {
    // An explicit subcommand wins over the argv[0] hint, so a binary that was
    // renamed or symlinked still answers to every tool it holds. No subcommand
    // name collides with an argument any of them takes, so `spirula-sfm auto`
    // falls through to the name check below and reaches the SfM tool.
    if (argc > 1) {
        if (const Tool* t = find_tool(argv[1])) {
            // The tool sees "spirula sfm" as its program name, so its own usage
            // text prints a command line that can be pasted back.
            std::string prog = std::string(argv[0] ? argv[0] : "spirula") + " " + t->name;
            std::vector<char*> sub;
            sub.push_back(prog.data());
            for (int i = 2; i < argc; i++) sub.push_back(argv[i]);
            sub.push_back(nullptr);
            return t->run((int)sub.size() - 1, sub.data());
        }
    }

    // Named as a tool: hand it everything, --help and --version included, so a
    // spirula-sfm symlink behaves exactly as the separate executable did.
    if (const Tool* t = tool_from_argv0(argc > 0 ? argv[0] : nullptr))
        return t->run(argc, argv);

    if (argc > 1) {
        const std::string a = argv[1];
        if (a == "help" || a == "--help" || a == "-h") {
            print_usage();
            return 0;
        }
        if (a == "--version" || a == "-v") {
            std::printf("%s\n", SS_VERSION);
            return 0;
        }
    }

#ifdef SS_TOOL_GUI
    return spirula_gui_main(argc, argv);
#else
    if (argc > 1)
        std::fprintf(stderr, "error: unknown command '%s'\n\n", argv[1]);
    else
        std::fprintf(stderr, "error: this build has no graphical application "
                             "(-DSS_BUILD_GUI=OFF)\n\n");
    print_usage();
    return 2;
#endif
}
