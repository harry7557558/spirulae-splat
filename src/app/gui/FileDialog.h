#pragma once

// The one picker every screen calls. It shows the DESKTOP's file dialog when
// this session has one (src/app/gui/NativeDialog.h) and falls back to the
// minimal ImGui browser below it when it does not -- a bare X server, a
// container, a machine with neither zenity nor kdialog installed. Callers see
// the same interface either way.

#include "app/gui/NativeDialog.h"

#include <string>
#include <vector>

namespace gui {

class FileDialog {
public:
    enum class Mode { Folder, File };

    // Whether to prefer the desktop's picker. Persisted by GuiApp; the
    // built-in browser is what a user who cannot get the system one to appear
    // (or does not want it) falls back to.
    void use_native(bool on) { _use_native = on; }
    bool native_enabled() const { return _use_native; }
    bool native_available() const { return NativeDialog::available(); }

    // Arm the dialog; it opens on the next draw() call. `extensions` filters
    // File mode (lowercase, with dot, e.g. {".mp4", ".mov"}); empty = all
    // files. `start_dir` = initial directory ("" = last used / home).
    // `multi_select` (File mode only) lets several files come back at once --
    // a dataset can be built from several video clips.
    void open(const std::string& title, Mode mode,
              const std::vector<std::string>& extensions = {},
              const std::string& start_dir = "", bool multi_select = false);

    // Draw the modal if open. Returns true exactly once when the user
    // confirmed a selection; the picked paths are in results(), the first of
    // them in result().
    bool draw();

    // Armed or on screen. ImGui shows one modal at a time, so a caller that is
    // itself a modal has to step aside for this one and needs to know when it
    // is gone -- whether the user confirmed or cancelled.
    bool is_open() const { return _is_open || _want_open || _native.busy(); }

    const std::string& result() const { return _result; }
    const std::vector<std::string>& results() const { return _results; }

private:
    void refresh();
    bool is_selected(const std::string& name) const;
    void toggle(const std::string& name);

    struct Entry { std::string name; bool is_dir = false; };

    std::string _title = "Select";
    Mode _mode = Mode::Folder;
    bool _multi = false;
    std::vector<std::string> _extensions;
    std::string _cwd;
    std::string _path_edit;          // editable path bar
    std::vector<Entry> _entries;
    // Basenames of the highlighted entries, in the order they were clicked.
    // Single-select keeps at most one.
    std::vector<std::string> _selected;
    std::string _result;             // _results[0], or "" -- the common case
    std::vector<std::string> _results;
    bool _want_open = false;
    bool _is_open = false;
    NativeDialog _native;
    bool _use_native = true;
};

}  // namespace gui
