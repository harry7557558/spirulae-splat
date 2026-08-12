#pragma once

// The desktop's own file picker: Explorer on Windows, Finder on macOS, and
// whatever the session's portal helper is on Linux (kdialog / zenity, both of
// which route through xdg-desktop-portal on a modern desktop).
//
// It is worth the platform code because the built-in browser in FileDialog.h
// has no places sidebar, no typeahead, no bookmarks and no idea where the
// user keeps things. That one stays as the fallback for a session where none
// of this is reachable -- a bare X server, a container, a machine with neither
// zenity nor kdialog installed.
//
// Everything here is asynchronous: the picker is a separate process (Linux) or
// a worker thread (Windows), so the GUI keeps drawing while it is up. macOS is
// the exception -- NSOpenPanel is main-thread only and runs modal, so the
// window stops repainting until it is dismissed.

#include <atomic>
#include <memory>
#include <string>
#include <thread>
#include <vector>

namespace gui {

class NativeDialog {
public:
    enum class Mode { Folder, File };

    ~NativeDialog();

    // Whether this build and this session can show one at all.
    static bool available();

    // Starts the picker. False when it could not be launched, which is the
    // caller's cue to fall back to the built-in browser. `extensions` are
    // lowercase with the dot (".mp4"); empty means every file.
    bool open(const std::string& title, Mode mode,
              const std::vector<std::string>& extensions,
              const std::string& start_dir, bool multi);

    bool busy() const { return _job != nullptr; }

    // True exactly once, when the picker has closed. results() is empty when
    // it was cancelled.
    bool poll();

    const std::vector<std::string>& results() const { return _results; }

    // What the worker writes into. Shared with it rather than owned, because
    // quitting must not wait on a picker: a Windows IFileDialog cannot be
    // dismissed from another thread at all, so the thread is abandoned and has
    // to be able to outlive this object safely.
    struct Job {
        std::atomic<bool> done{false};
        std::atomic<bool> cancel{false};
        std::vector<std::string> picked;   // written before `done` is set
    };

private:
    std::thread _worker;
    std::shared_ptr<Job> _job;
    std::vector<std::string> _results;
};

}  // namespace gui
