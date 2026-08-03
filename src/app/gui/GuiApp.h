#pragma once

// GuiApp -- top-level UI for Spirulae Splat: screens, layout, state wiring
// between the config editor, COLMAP runner, train runner, and the native
// viewport.

#include "backend/api/BackendRuntime.h"
#include "app/generated/cli_config.h"
#include "app/gui/ColmapRunner.h"
#include "app/gui/ConfigUI.h"
#include "app/gui/FileDialog.h"
#include "app/gui/ModelCache.h"
#include "app/gui/SegmentPanel.h"
#include "app/gui/SfmRunner.h"
#include "app/gui/TrainRunner.h"
#include "app/gui/ViewportPanel.h"

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
    enum class Screen { Home, NewDataset, Train };
    enum class PickAction {
        None, OpenDataset, SourceImages, SourceVideo, SourceReplace, Workspace,
        OutputPrefix, VocabTree, MaskModelFile
    };
    // Which reconstruction back end the New Dataset screen runs.
    enum class Engine { BuiltIn, Colmap };
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
    // whether each path is an SfM dataset folder, a folder of photos, or a
    // video file, and routes them like the corresponding Home-screen action.
    // Several videos dropped together become several inputs of one dataset.
    void handle_drop(const std::vector<std::string>& paths);

private:
    void request_go_home();
    void apply_preset(const std::string& preset);
    void start_training();
    bool training_busy() const;   // Preparing or Training

    // ---- dataset creation ----
    // Which engines this build and this machine can actually offer.
    bool builtin_sfm_available() const;
    bool colmap_available() const;
    Engine effective_engine() const;
    bool dataset_busy() const;
    void start_dataset_job();
    void cancel_dataset_job();
    // Copies the panel-level state into whichever job struct will run.
    void sync_dataset_jobs();
    // Path of the selected checkpoint, or "" when it is not downloaded yet.
    std::string selected_model_path() const;
    bool license_accepted(const std::string& family) const;

    // ---- screens ----
    void draw_menu_bar();
    void draw_home();
    void draw_new_dataset();
    void draw_dataset_source();       // input list / output / resume
    void draw_dataset_basics();       // the four or five knobs a beginner needs
    void draw_source_cameras();       // one lens per input, when there are several
    void draw_masking_options();
    void draw_sfm_advanced();
    void draw_colmap_options();
    void draw_tool_locations();
    void draw_license_modal();
    void draw_train();
    void draw_train_settings();      // left panel
    void draw_basic_options();
    void draw_train_controls();
    void draw_metrics();
    void draw_status_strip();
    // Right-aligned VRAM readout on the status strip; x0/avail describe the
    // strip's content region (window-local left edge and width).
    void draw_vram_readout(float x0, float avail);
    void draw_log_panel(float height);
    void draw_confirm_modal();
    void handle_dialog_result(const std::vector<std::string>& paths);
    // Take paths onto the input list, `replace` clearing what was there (a
    // fresh pick from Home) rather than adding to it (the panel's Add buttons).
    // Sets the per-input defaults and, unless the user has edited it, the
    // output folder.
    void add_sources(const std::vector<std::string>& paths, bool replace);
    // Re-derive what is a function of the list: the sub-folder each input's
    // images go into, and the default workspace.
    void refresh_sources();
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

    // Dataset creation. Both runners exist; only one runs, chosen by _engine
    // (and forced when only one is available).
    Engine _engine = Engine::BuiltIn;
    ColmapRunner _colmap;
    ColmapJob _colmap_job;
    SfmRunner _sfm;
    SfmJob _sfm_job;
    // Panel-level state, copied into whichever job runs. The inputs are kept as
    // the struct both runners take (PrepInput), so the panel edits the thing
    // that runs instead of a parallel copy of it: a video file or photo folder
    // each, plus the sub-folder and the lens that belong to it.
    std::vector<PrepInput> _sources;
    std::string _workspace;
    // The output folder this screen derived from the inputs. Kept so a folder
    // the user typed is never overwritten when the input list changes.
    std::string _workspace_auto;
    bool _resume = true;
    bool _mask_enable = false;
    MaskSettings _mask;
    SegmentPanel _segment;
    // Which input "Try the mask" runs on, and so which input the clicked
    // objects are prompts for (MaskSettings::clicks_source).
    int _mask_preview_input = 0;

    // Segmentation checkpoints.
    std::string _model_id = "sam3-q4_0";
    ModelDownload _download;
    // Families whose licence the user has accepted, persisted in the settings.
    std::vector<std::string> _accepted_licenses;
    std::string _license_prompt;      // family whose modal is open
    bool _license_tick = false;

    FileDialog _dialog;
    PickAction _pick = PickAction::None;
    int _pick_source = -1;            // which input PickAction::SourceReplace edits

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

    // VRAM readout on the status strip, polled from the backend at ~2 Hz.
    backend::MemoryUsage _vram;
    double _vram_polled_at = -1.0;
};

}  // namespace gui
