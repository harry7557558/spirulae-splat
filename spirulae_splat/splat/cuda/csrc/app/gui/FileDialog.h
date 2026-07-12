#pragma once

// Minimal built-in ImGui file/folder browser (modal popup) -- no native
// dialog dependency, so it behaves identically on Linux/Windows/macOS and
// over remote X/WSLg sessions.

#include <string>
#include <vector>

namespace gui {

class FileDialog {
public:
    enum class Mode { Folder, File };

    // Arm the dialog; it opens on the next draw() call. `extensions` filters
    // File mode (lowercase, with dot, e.g. {".mp4", ".mov"}); empty = all
    // files. `start_dir` = initial directory ("" = last used / home).
    void open(const std::string& title, Mode mode,
              const std::vector<std::string>& extensions = {},
              const std::string& start_dir = "");

    // Draw the modal if open. Returns true exactly once when the user
    // confirmed a selection; the picked path is in result().
    bool draw();

    const std::string& result() const { return _result; }

private:
    void refresh();

    struct Entry { std::string name; bool is_dir = false; };

    std::string _title = "Select";
    Mode _mode = Mode::Folder;
    std::vector<std::string> _extensions;
    std::string _cwd;
    std::string _path_edit;          // editable path bar
    std::vector<Entry> _entries;
    std::string _selected;           // basename of the selected entry
    std::string _result;
    bool _want_open = false;
    bool _is_open = false;
};

}  // namespace gui
