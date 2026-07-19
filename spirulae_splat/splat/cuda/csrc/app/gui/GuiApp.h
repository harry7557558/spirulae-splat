#pragma once

// GuiApp -- top-level UI for Spirulae Splat: screens, layout, state wiring
// between the config editor, COLMAP runner, train runner, and the native
// viewport.

#include "../generated/cli_config.h"
#include "ColmapRunner.h"
#include "ConfigUI.h"
#include "FileDialog.h"
#include "TrainRunner.h"
#include "ViewportPanel.h"

#include <deque>
#include <string>
#include <vector>

namespace gui {

class GuiApp {
public:
    GuiApp();
    ~GuiApp();

    // Draw one frame (between ImGui::NewFrame and ImGui::Render).
    void frame();

    // Window close button pressed; may open a confirmation dialog instead
    // of quitting when training is in flight.
    void request_close();
    bool wants_exit() const { return _quit; }

    // Stop worker threads + free GL resources; call while the GL context is
    // still current.
    void shutdown();

private:
    enum class Screen { Home, Colmap, Train };
    enum class PickAction {
        None, OpenDataset, ColmapImages, ColmapVideo, ColmapWorkspace,
        OutputPrefix, VocabTree
    };
    // Session-destroying actions deferred behind the "stop training?"
    // confirmation (and executed once the stop has completed).
    enum class Pending { None, GoHome, OpenDataset, Quit };

    // ---- persistence (recents + tool paths) ----
    static std::string settings_path();
    void load_settings();
    void save_settings();
    // By value: callers pass elements of _recents, which this mutates.
    void add_recent(std::string path);

    // ---- actions ----
    // By value: callers pass elements of _recents, which open_dataset
    // mutates via add_recent (a const& here would dangle).
    void open_dataset(std::string dir, std::string image_dir = "");
    // Route for user-initiated opens: confirms first when training.
    void request_open_dataset(std::string dir);

public:
    // Drag-and-drop entry (GLFW drop callback, main thread): auto-detects
    // whether the path is an SfM dataset folder, a folder of photos, or a
    // video file, and routes it like the corresponding Home-screen action.
    void handle_drop(const std::string& path);

private:
    void request_go_home();
    void apply_preset(const std::string& preset);
    void start_training();
    bool training_busy() const;   // Preparing or Training

    // ---- screens ----
    void draw_menu_bar();
    void draw_home();
    void draw_colmap();
    void draw_colmap_options();
    void draw_train();
    void draw_train_settings();      // left panel
    void draw_basic_options();
    void draw_train_controls();
    void draw_metrics();
    void draw_status_strip();
    void draw_log_panel(float height);
    void draw_confirm_modal();
    void handle_dialog_result(const std::string& path);
    void run_pending_if_stopped();
    void append_logs();
    void log(const std::string& s);

    Screen _screen = Screen::Home;
    bool _quit = false;
    bool _open_confirm = false;      // arm the stop-training modal
    bool _confirm_shown = false;     // modal currently expected open
    bool _stop_confirmed = false;    // user chose "Stop & Save"
    Pending _pending = Pending::None;
    std::string _pending_path;       // dataset dir for Pending::OpenDataset
    bool _parse_dirty = false;       // dataparser option edited -> reload
    bool _device_locked = false;     // backend initialized -> device fixed

    // Config being edited + the preset baseline it diffs against.
    SsplatConfig _cfg;
    SsplatConfig _defaults;
    std::string _preset = "3dgs";
    ConfigUIState _cfg_ui;

    TrainRunner _runner;
    ViewportPanel _viewport;

    ColmapRunner _colmap;
    ColmapJob _colmap_job;

    FileDialog _dialog;
    PickAction _pick = PickAction::None;

    // Settings (persisted).
    std::vector<std::string> _recents;
    std::string _colmap_exe = "colmap";
    std::string _ffmpeg_exe = "ffmpeg";
#ifdef _WIN32
    std::string _python_exe = "python";
#else
    std::string _python_exe = "python3";
#endif

    // Log console.
    std::deque<std::string> _log;
    bool _log_autoscroll = true;
    bool _show_log = true;
};

}  // namespace gui
