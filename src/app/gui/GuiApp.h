#pragma once

// GuiApp -- top-level UI: screens, layout, state wiring
// between the config editor, COLMAP runner, train runner, and the native
// viewport.

#include "backend/api/BackendRuntime.h"
#include "config/TrainConfig.h"
#include "app/gui/BatchTrain.h"
#include "app/gui/ColmapRunner.h"
#include "app/gui/Fonts.h"
#include "app/gui/ConfigUI.h"
#include "app/gui/FileDialog.h"
#include "app/gui/ImageCompare.h"
#include "app/gui/MeshRunner.h"
#include "app/gui/ModelCache.h"
#include "app/gui/SegmentPanel.h"
#include "app/gui/SfmRunner.h"
#include "app/gui/SplatViewer.h"
#include "app/gui/TrainPreset.h"
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

    // The font atlas. GuiMain calls ensure() between frames -- swapping a face
    // invalidates every ImFont pointer, so it cannot happen mid-frame.
    FontSet& fonts() { return _fonts; }

private:
    enum class Screen { Home, NewDataset, Train, Viewer, Batch, Mesh };
    enum class PickAction {
        None, OpenDataset, SourceImages, SourceVideo, SourceReplace, Workspace,
        OutputPrefix, VocabTree, MaskModelFile, SplatFile,
        PresetFile, PresetSaveFolder, BatchDataset, BatchOutput, BatchPresetFile,
        MeshSource, MeshPhotos, MeshOutput
    };
    // Which reconstruction back end the New Dataset screen runs.
    enum class Engine { BuiltIn, Colmap };
    // Session-destroying actions deferred behind the "stop training?"
    // confirmation (and executed once the stop has completed).
    enum class Pending { None, GoHome, OpenDataset, OpenSplat, Quit, StartBatch };

    // ---- persistence (recents + tool paths) ----
    static std::string settings_path();
    void load_settings();
    void save_settings();
    // By value: callers pass elements of _recents, which this mutates.
    void add_recent(std::string path);

    // ---- actions ----
    // By value: callers pass elements of _recents, which open_dataset
    // mutates via add_recent (a const& here would dangle).
    // Clears the log panel unless `keep_log`: what is in it belongs to
    // whatever was open before. The reconstruction handoff passes true,
    // because there the log is this dataset's own build log.
    void open_dataset(std::string dir, std::string image_dir = "",
                      std::string mask_dir = "", bool keep_log = false);
    // Route for user-initiated opens: confirms first when training.
    void request_open_dataset(std::string dir);

    // The viewer screen: a splat file (or a checkpoint / run directory) opened
    // for looking at. Takes the engine over, so it goes through the same
    // confirmation as any other session-destroying action.
    void open_splat(std::string path);
    void request_open_splat(std::string path);
    // Give the engine back and leave the screen. Called before anything that
    // needs the engine for itself.
    void close_splat();

public:
    // Drag-and-drop entry (GLFW drop callback, main thread): auto-detects
    // whether each path is an SfM dataset folder, a folder of photos, or a
    // video file, and routes them like the corresponding Home-screen action.
    // Several videos dropped together become several inputs of one dataset.
    void handle_drop(const std::vector<std::string>& paths);

private:
    void request_go_home();
    void apply_preset(const std::string& preset);
    // A preset the user saved (or a run's config.json). Replaces the whole
    // config the way a built-in preset does, but keeps the GUI-owned context
    // -- which dataset is open, where its runs go.
    void apply_user_preset(const TrainPreset& p);
    // Read `path` and apply it. Reports the outcome in the preset panel and
    // the log rather than throwing; refuses while training, because applying
    // one re-parses the dataset and would take the running session down.
    void load_preset_file(const std::string& path);
    // Re-scan preset_dir(), rate-limited: this runs while a dropdown is open.
    void refresh_presets();
    void start_training();
    // Everything that renders from the training session, released together.
    // Every path that replaces or destroys the session goes through this.
    void detach_session_views();
    // Hand a config to the runner. Everything a new session invalidates --
    // the splat viewer holding the engine, the viewport's render worker --
    // is released here, so both the Start button and the batch queue go
    // through it.
    void launch_training(const TrainConfig& cfg, const std::string& preset);
    bool training_busy() const;   // Preparing or Training

    // ---- batch ----
    // Append a row for this dataset, seeded with the preset the trainer
    // screen is on. Shared by the picker, the recents menu and the drop
    // handler.
    void add_batch_row(const std::string& dataset);
    void check_batch();                     // pre-flight every row
    void request_start_batch(bool skip_invalid);
    void start_batch(bool skip_invalid);
    // Called once per frame while a batch is live: records the row that just
    // finished and launches the next one.
    void advance_batch();
    void finish_batch();
    // Give up on the queue without waiting for the current row (the stop
    // confirmation took the session away).
    void cancel_batch();

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
    // Opens "Try the mask" on the input the combo points at, which is also the
    // input whose clicks and stencil it edits.
    void open_mask_preview();
    void draw_sfm_advanced();
    void draw_colmap_options();
    void draw_tool_locations();
    void draw_license_modal();
    void draw_language_menu();       // the picker, and the CJK font prompt
    void draw_train();
    void draw_viewer();
    void draw_mesh();
    void draw_mesh_options();
    void draw_mesh_preview(float height);
    // Set the meshing source (and derive the output name from it). Shared by
    // the picker, the drop handler and the "mesh this run" shortcut.
    void set_mesh_source(const std::string& path);
    void start_meshing();
    // Would the run about to start actually get cameras? Resolves the dataset
    // the way the child does -- the folder typed here, else the `data` entry
    // in the run's config.json -- so the screen can warn BEFORE the run that
    // the mesh is about to be the density-only, much rougher kind. Touches
    // the disk, so the answer is cached until the source or folder changes.
    bool mesh_dataset_found();
    // Load whatever the run wrote, plus the model it came from, into the two
    // preview panels. No-op (mesh side only) while training holds the engine.
    void open_mesh_preview();
    void close_mesh_preview();
    void draw_batch();
    void draw_batch_table();
    void draw_batch_preset_combo(BatchJob& job, int row);
    void draw_batch_issues();
    void draw_batch_progress();      // the running-job block on the trainer screen
    void draw_train_settings();      // left panel
    void draw_preset_picker();       // built-in + saved presets, save / load
    void draw_preset_save_modal();
    void draw_preset_delete_modal();
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
    void rescan_found_masks();
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
    bool _pending_batch_skip = false;  // Pending::StartBatch's argument
    bool _parse_dirty = false;       // dataparser option edited -> reload
    bool _device_locked = false;     // backend initialized -> device fixed

    // Config being edited + the preset baseline it diffs against.
    TrainConfig _cfg;
    TrainConfig _defaults;
    // The built-in preset the config descends from. Still set when a saved
    // preset is in use -- it is what the run's config.json records, and what
    // the options editor's "preset default" tooltips are relative to.
    std::string _preset = "3dgs";
    ConfigUIState _cfg_ui;

    // Saved presets. `_preset_file` is "" while a built-in preset is selected;
    // otherwise it is the file in use and names the row shown in the picker.
    std::vector<TrainPreset> _presets;
    double _presets_scanned_at = -1.0;
    std::string _preset_file;
    std::string _preset_display;
    std::string _preset_desc;
    // The last save / load outcome, already formatted, and whether it failed.
    std::string _preset_msg;
    bool _preset_msg_err = false;
    // The save dialog.
    bool _preset_save_open = false;    // arm the modal
    bool _preset_save_shown = false;
    // The modal stepped aside for the file dialog and wants to come back.
    bool _preset_save_reopen = false;
    std::string _preset_save_name;
    std::string _preset_save_desc;
    std::string _preset_save_path;
    // The path field follows the name until the user edits it by hand.
    bool _preset_path_edited = false;
    // The delete confirmation, which always targets _preset_file.
    bool _preset_delete_open = false;
    bool _preset_delete_shown = false;

    TrainRunner _runner;
    ViewportPanel _viewport;
    // The trainer screen's other preview: one training photograph beside the
    // render of the same camera. Reads from the session's engine, so it is
    // released together with the viewport (detach_session_views).
    ImageCompare _images;
    bool _preview_images = false;    // which of the two the trainer screen shows
    SplatViewer _splat;
    // The viewport is showing _splat rather than _runner's session.
    bool _viewing_splat = false;

    // ---- meshing ----
    // The extraction runs as a child process (MeshRunner); the preview after
    // it is two panels of the same scene -- the splats on the engine
    // renderer, the mesh on the GL one -- with their navigation linked.
    MeshRunner _mesh;
    MeshJob _mesh_job;
    SplatViewer _mesh_view;          // the produced file, as geometry
    ViewportPanel _mesh_viewport;    // the right-hand panel
    bool _mesh_preview_open = false; // the mesh side is attached
    bool _mesh_splats_open = false;  // the splat side is attached too
    bool _mesh_link_views = true;
    // The MeshRunner::run_id() whose result has already been opened. The
    // runner stays Done for the rest of the session, so without this the
    // screen would reopen the preview the frame after it is closed -- which
    // is what made "make another mesh" impossible without a restart.
    uint64_t _mesh_shown_run = 0;
    // mesh_dataset_found()'s cache: the (source, folder) it was asked about
    // and what it answered.
    std::string _mesh_data_probe_key;
    bool _mesh_data_probe_found = false;

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
    // The static stencils are kept on the inputs themselves (PrepInput::
    // stencil); this only says whether the run is given them, so that turning
    // the option off and on again does not throw away what was drawn.
    bool _border_enable = false;
    MaskSettings _mask;
    SegmentPanel _segment;
    // Which input "Try the mask" runs on, and so which input the clicked
    // objects are prompts for (MaskSettings::clicks_source) and which one's
    // stencil the panel edits.
    int _mask_preview_input = 0;

    // Segmentation checkpoints.
    std::string _model_id = "sam3-q4_0";
    ModelDownload _download;

    // Interface language and the glyphs to draw it with. The font download is
    // separate from _download so that fetching a face cannot cancel a
    // half-finished 700 MB checkpoint.
    FontSet _fonts;
    FileDownload _font_download;
    const CjkFace* _font_fetching = nullptr;
    // Families whose licence the user has accepted, persisted in the settings.
    std::vector<std::string> _accepted_licenses;
    std::string _license_prompt;      // family whose modal is open
    bool _license_tick = false;

    // Batch training. The queue is data; the driver is advance_batch(), so a
    // running batch is an ordinary training session the trainer screen shows
    // exactly as it shows a hand-started one.
    std::vector<BatchJob> _batch;
    bool _batch_dirty = false;        // edited -> persist once the widget is idle
    bool _batch_checked = false;      // a pre-flight has run since the last edit
    bool _batch_active = false;
    bool _batch_launched = false;     // a row's session is in flight
    int  _batch_current = -1;         // which row that is
    bool _batch_stop_after = false;   // finish this row, then stop
    bool _batch_stop_now = false;     // ... and record it as stopped, not done
    std::string _batch_msg;           // already formatted; "" when there is none
    bool _batch_msg_err = false;

    FileDialog _dialog;
    PickAction _pick = PickAction::None;
    int _pick_source = -1;            // which input PickAction::SourceReplace edits
    // Which batch row the pending pick edits; -1 appends a new row.
    int _pick_row = -1;

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
    // Masks sitting beside the photos are adopted automatically; this is the
    // way out for a folder whose masks/ describes something else.
    bool _use_found_masks = true;

    // VRAM readout on the status strip, polled from the backend at ~2 Hz.
    backend::MemoryUsage _vram;
    double _vram_polled_at = -1.0;
};

}  // namespace gui
