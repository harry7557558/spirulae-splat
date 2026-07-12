#pragma once

// GuiApp -- top-level UI for Spirulae Splat Studio: screens, layout, state
// wiring between the config editor, COLMAP runner, train runner, and the
// native viewport.

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
        None, OpenDataset, ColmapImages, ColmapVideo, ColmapWorkspace
    };

    // ---- persistence (recents + tool paths) ----
    static std::string settings_path();
    void load_settings();
    void save_settings();
    void add_recent(const std::string& path);

    // ---- actions ----
    void open_dataset(const std::string& dir, const std::string& image_dir = "");
    void apply_preset(const std::string& preset);
    void start_training();

    // ---- screens ----
    void draw_menu_bar();
    void draw_home();
    void draw_colmap();
    void draw_train();
    void draw_train_settings();      // left panel
    void draw_basic_options();
    void draw_train_controls();
    void draw_metrics();
    void draw_status_strip();
    void draw_log_panel(float height);
    void draw_exit_confirm();
    void handle_dialog_result(const std::string& path);
    void append_logs();
    void log(const std::string& s);

    Screen _screen = Screen::Home;
    bool _quit = false;
    bool _confirm_exit = false;

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

    // Log console.
    std::deque<std::string> _log;
    bool _log_autoscroll = true;
    bool _show_log = true;
};

}  // namespace gui
