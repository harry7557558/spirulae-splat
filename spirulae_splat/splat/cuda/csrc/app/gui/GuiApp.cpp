// GuiApp.cpp -- see GuiApp.h.

#include "GuiApp.h"

#include "../Json.h"
#include "Subprocess.h"

#include "imgui.h"
#include "imgui_stdlib.h"

#include <algorithm>
#include <cfloat>
#include <cstdio>
#include <cstdlib>
#include <filesystem>

namespace fs = std::filesystem;

namespace gui {

namespace {

const ImVec4 kOk(0.35f, 0.85f, 0.45f, 1.0f);
const ImVec4 kErr(1.0f, 0.42f, 0.42f, 1.0f);
const ImVec4 kDim(0.6f, 0.6f, 0.6f, 1.0f);

std::string format_eta(double s) {
    if (s < 0) return "--:--";
    int t = (int)(s + 0.5);
    char buf[32];
    if (t >= 3600) std::snprintf(buf, sizeof buf, "%d:%02d:%02d", t/3600, (t/60)%60, t%60);
    else           std::snprintf(buf, sizeof buf, "%d:%02d", t/60, t%60);
    return buf;
}

std::string format_count(double n) {
    char buf[32];
    if (n >= 1e6) std::snprintf(buf, sizeof buf, "%.2fM", n / 1e6);
    else if (n >= 1e3) std::snprintf(buf, sizeof buf, "%.1fk", n / 1e3);
    else std::snprintf(buf, sizeof buf, "%d", (int)n);
    return buf;
}

const char* preset_help(const std::string& name) {
    for (const auto& p : kSsplatPresets)
        if (name == p.name) return p.help;
    return "";
}

}  // namespace


// ===========================================================================
// Lifecycle + persistence
// ===========================================================================

GuiApp::GuiApp() {
    load_settings();
    apply_preset("3dgs");
    _colmap_job.colmap_exe = _colmap_exe;
    _colmap_job.ffmpeg_exe = _ffmpeg_exe;
}

GuiApp::~GuiApp() = default;

void GuiApp::shutdown() {
    save_settings();
    _viewport.detach();
    _viewport.destroy_gl();
    _colmap.cancel();
    _runner.shutdown();
}

std::string GuiApp::settings_path() {
#ifdef _WIN32
    const char* base = std::getenv("APPDATA");
    fs::path dir = base ? fs::path(base) : fs::path(".");
#else
    fs::path dir;
    if (const char* x = std::getenv("XDG_CONFIG_HOME")) dir = x;
    else if (const char* h = std::getenv("HOME")) dir = fs::path(h) / ".config";
    else dir = ".";
#endif
    dir /= "spirulae-splat";
    std::error_code ec;
    fs::create_directories(dir, ec);
    return (dir / "gui.conf").string();
}

void GuiApp::load_settings() {
    FILE* f = std::fopen(settings_path().c_str(), "r");
    if (!f) return;
    char line[1024];
    while (std::fgets(line, sizeof line, f)) {
        std::string s = line;
        while (!s.empty() && (s.back() == '\n' || s.back() == '\r')) s.pop_back();
        size_t eq = s.find('=');
        if (eq == std::string::npos) continue;
        std::string k = s.substr(0, eq), v = s.substr(eq + 1);
        if (k == "recent" && !v.empty()) _recents.push_back(v);
        else if (k == "colmap_exe" && !v.empty()) _colmap_exe = v;
        else if (k == "ffmpeg_exe" && !v.empty()) _ffmpeg_exe = v;
    }
    std::fclose(f);
}

void GuiApp::save_settings() {
    FILE* f = std::fopen(settings_path().c_str(), "w");
    if (!f) return;
    for (const auto& r : _recents)
        std::fprintf(f, "recent=%s\n", r.c_str());
    std::fprintf(f, "colmap_exe=%s\n", _colmap_exe.c_str());
    std::fprintf(f, "ffmpeg_exe=%s\n", _ffmpeg_exe.c_str());
    std::fclose(f);
}

void GuiApp::add_recent(const std::string& path) {
    _recents.erase(std::remove(_recents.begin(), _recents.end(), path),
                   _recents.end());
    _recents.insert(_recents.begin(), path);
    if (_recents.size() > 10) _recents.resize(10);
}

void GuiApp::log(const std::string& s) {
    _log.push_back(s);
    while (_log.size() > 4000) _log.pop_front();
}

void GuiApp::append_logs() {
    for (auto& s : _runner.drain_log()) log(s);
    for (auto& s : _colmap.drain_log()) log(s);
}


// ===========================================================================
// Actions
// ===========================================================================

void GuiApp::apply_preset(const std::string& preset) {
    SsplatConfig fresh;
    ssplat_apply_preset(fresh, preset);
    // Keep GUI-managed context across preset switches.
    fresh.data = _cfg.data;
    fresh.image_dir = _cfg.image_dir;
    fresh.output_dir_prefix = _cfg.output_dir_prefix;
    fresh.output_dir_name = _cfg.output_dir_name;
    // The native viewport replaces the web viewer by default; it can be
    // re-enabled in Basic Options for remote monitoring.
    fresh.disable_viewer = true;
    _preset = preset;
    _cfg = fresh;
    _defaults = fresh;
    if (!_cfg.data.empty()) {
        _viewport.detach();
        _runner.load_dataset(_cfg, _preset);
    }
}

void GuiApp::open_dataset(const std::string& dir, const std::string& image_dir) {
    _cfg.data = dir;
    if (!image_dir.empty()) {
        _cfg.image_dir = image_dir;
    } else {
        // COLMAP-runner datasets record their (possibly external) image dir.
        std::error_code ec;
        fs::path marker = fs::path(dir) / "ssplat_dataset.json";
        if (fs::exists(marker, ec)) {
            try {
                JsonValue v = json_parse_file(marker.string());
                if (const JsonValue* d = v.find("image_dir"))
                    if (!d->as_string().empty()) _cfg.image_dir = d->as_string();
            } catch (...) {}
        }
    }
    _defaults.data = _cfg.data;
    _defaults.image_dir = _cfg.image_dir;
    add_recent(dir);
    save_settings();
    _viewport.detach();
    _runner.load_dataset(_cfg, _preset);
    _screen = Screen::Train;
}

void GuiApp::start_training() {
    if (_cfg.data.empty()) return;
    _viewport.detach();
    _runner.start_training(_cfg, _preset);
}

void GuiApp::request_close() {
    TrainRunner::Phase ph = _runner.phase();
    if (ph == TrainRunner::Phase::Training || ph == TrainRunner::Phase::Preparing)
        _confirm_exit = true;
    else
        _quit = true;
}

void GuiApp::handle_dialog_result(const std::string& path) {
    switch (_pick) {
        case PickAction::OpenDataset:
            open_dataset(path);
            break;
        case PickAction::ColmapImages:
            _colmap_job.input_path = path;
            _colmap_job.is_video = false;
            if (_colmap_job.workspace.empty())
                _colmap_job.workspace = path + "_dataset";
            _screen = Screen::Colmap;
            break;
        case PickAction::ColmapVideo: {
            _colmap_job.input_path = path;
            _colmap_job.is_video = true;
            fs::path p(path);
            _colmap_job.workspace = (p.parent_path() / (p.stem().string() + "_dataset")).string();
            _screen = Screen::Colmap;
            break;
        }
        case PickAction::ColmapWorkspace:
            _colmap_job.workspace = path;
            break;
        default:
            break;
    }
    _pick = PickAction::None;
}


// ===========================================================================
// Frame
// ===========================================================================

void GuiApp::frame() {
    append_logs();

    // Attach the native viewport once the engine is live.
    if (_runner.engine_ready() && !_viewport.attached() && _runner.session())
        _viewport.attach(*_runner.session());

    const ImGuiViewport* vp = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos(vp->WorkPos);
    ImGui::SetNextWindowSize(vp->WorkSize);
    ImGuiWindowFlags flags =
        ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove |
        ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoBringToFrontOnFocus |
        ImGuiWindowFlags_MenuBar;
    ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.0f);
    ImGui::Begin("##host", nullptr, flags);
    ImGui::PopStyleVar();

    draw_menu_bar();
    switch (_screen) {
        case Screen::Home:   draw_home();   break;
        case Screen::Colmap: draw_colmap(); break;
        case Screen::Train:  draw_train();  break;
    }

    if (_dialog.draw()) handle_dialog_result(_dialog.result());
    draw_exit_confirm();

    ImGui::End();
}

void GuiApp::draw_menu_bar() {
    if (!ImGui::BeginMenuBar()) return;
    if (ImGui::BeginMenu("File")) {
        if (ImGui::MenuItem("Open Dataset Folder...")) {
            _pick = PickAction::OpenDataset;
            _dialog.open("Open Dataset Folder", FileDialog::Mode::Folder);
        }
        if (ImGui::MenuItem("New Dataset from Photos...")) {
            _pick = PickAction::ColmapImages;
            _dialog.open("Select Photo Folder", FileDialog::Mode::Folder);
        }
        if (ImGui::MenuItem("New Dataset from Video...")) {
            _pick = PickAction::ColmapVideo;
            _dialog.open("Select Video File", FileDialog::Mode::File,
                         {".mp4", ".mov", ".avi", ".mkv", ".webm", ".m4v", ".insv"});
        }
        ImGui::Separator();
        if (ImGui::MenuItem("Quit")) request_close();
        ImGui::EndMenu();
    }
    if (ImGui::BeginMenu("View")) {
        ImGui::MenuItem("Show Log Panel", nullptr, &_show_log);
        ImGui::EndMenu();
    }
    if (ImGui::BeginMenu("Help")) {
        if (ImGui::BeginMenu("About")) {
            ImGui::TextUnformatted("Spirulae Splat Studio");
            ImGui::TextDisabled("Native trainer GUI for spirulae-splat.");
            ImGui::TextDisabled("github.com/harry7557558/spirulae-splat");
            ImGui::EndMenu();
        }
        ImGui::EndMenu();
    }
    ImGui::EndMenuBar();
}


// ===========================================================================
// Home
// ===========================================================================

void GuiApp::draw_home() {
    float w = 480.0f;
    ImGui::SetCursorPosX((ImGui::GetWindowWidth() - w) * 0.5f);
    ImGui::BeginChild("##home", ImVec2(w, 0));
    ImGui::Dummy(ImVec2(0, 40));

    ImGui::SetWindowFontScale(1.7f);
    ImGui::TextUnformatted("Spirulae Splat Studio");
    ImGui::SetWindowFontScale(1.0f);
    ImGui::TextDisabled("Reconstruct 3D scenes from photos with Gaussian splatting.");
    ImGui::Dummy(ImVec2(0, 24));

    if (ImGui::Button("Open a Dataset...", ImVec2(-1, 42))) {
        _pick = PickAction::OpenDataset;
        _dialog.open("Open Dataset Folder", FileDialog::Mode::Folder);
    }
    help_tooltip_on_hover(
        "A folder containing an already-processed dataset: COLMAP "
        "(sparse/0), Nerfstudio (transforms.json), or Metashape exports. "
        "The format is detected automatically.");

    if (ImGui::Button("Create Dataset from Photos...", ImVec2(-1, 42))) {
        _pick = PickAction::ColmapImages;
        _dialog.open("Select Photo Folder", FileDialog::Mode::Folder);
    }
    help_tooltip_on_hover(
        "Pick a folder of overlapping photos of a scene or object. COLMAP "
        "(installed separately) estimates the camera poses; the result "
        "opens directly in the trainer.");

    if (ImGui::Button("Create Dataset from Video...", ImVec2(-1, 42))) {
        _pick = PickAction::ColmapVideo;
        _dialog.open("Select Video File", FileDialog::Mode::File,
                     {".mp4", ".mov", ".avi", ".mkv", ".webm", ".m4v", ".insv"});
    }
    help_tooltip_on_hover(
        "Pick a video walking around a scene or object. Frames are "
        "extracted with ffmpeg, then processed with COLMAP.");

    if (!_recents.empty()) {
        ImGui::Dummy(ImVec2(0, 18));
        ImGui::SeparatorText("Recent");
        for (size_t i = 0; i < _recents.size(); i++) {
            ImGui::PushID((int)i);
            if (ImGui::Selectable(_recents[i].c_str()))
                open_dataset(_recents[i]);
            ImGui::PopID();
        }
    }

    ImGui::Dummy(ImVec2(0, 18));
    if (!command_exists(_colmap_exe))
        ImGui::TextColored(kDim, "note: 'colmap' was not found on PATH -- "
                           "dataset creation needs it (training does not).");

    ImGui::EndChild();
}


// ===========================================================================
// COLMAP screen
// ===========================================================================

void GuiApp::draw_colmap() {
    bool running = _colmap.state() == ColmapRunner::State::Running;

    if (ImGui::Button("< Home")) _screen = Screen::Home;
    ImGui::SameLine();
    ImGui::SetWindowFontScale(1.2f);
    ImGui::TextUnformatted(_colmap_job.is_video
                               ? "Create Dataset from Video"
                               : "Create Dataset from Photos");
    ImGui::SetWindowFontScale(1.0f);
    ImGui::Spacing();

    ImGui::BeginDisabled(running);

    // ---- input / output ----
    ImGui::SetNextItemWidth(-220);
    ImGui::InputText("##in", &_colmap_job.input_path);
    ImGui::SameLine();
    if (ImGui::Button("Browse...##in")) {
        if (_colmap_job.is_video) {
            _pick = PickAction::ColmapVideo;
            _dialog.open("Select Video File", FileDialog::Mode::File,
                         {".mp4", ".mov", ".avi", ".mkv", ".webm", ".m4v", ".insv"});
        } else {
            _pick = PickAction::ColmapImages;
            _dialog.open("Select Photo Folder", FileDialog::Mode::Folder);
        }
    }
    ImGui::SameLine();
    ImGui::TextUnformatted(_colmap_job.is_video ? "video file" : "photo folder");

    ImGui::SetNextItemWidth(-220);
    ImGui::InputText("##ws", &_colmap_job.workspace);
    ImGui::SameLine();
    if (ImGui::Button("Browse...##ws")) {
        _pick = PickAction::ColmapWorkspace;
        _dialog.open("Select Output Folder", FileDialog::Mode::Folder);
    }
    ImGui::SameLine();
    ImGui::TextUnformatted("output dataset folder");

    // ---- options ----
    ImGui::Spacing();
    const char* qualities[] = {"Fast", "Balanced", "High quality"};
    ImGui::SetNextItemWidth(180);
    ImGui::Combo("Quality", &_colmap_job.quality, qualities, 3);
    help_tooltip_on_hover("Feature count and image size used for matching. "
                          "Higher finds more cameras in difficult scenes but "
                          "is slower.");

    static const char* cam_models[] = {
        "OPENCV", "PINHOLE", "SIMPLE_PINHOLE", "SIMPLE_RADIAL", "OPENCV_FISHEYE"};
    int cam_idx = 0;
    for (int i = 0; i < 5; i++)
        if (_colmap_job.camera_model == cam_models[i]) cam_idx = i;
    ImGui::SetNextItemWidth(180);
    if (ImGui::Combo("Camera model", &cam_idx, cam_models, 5))
        _colmap_job.camera_model = cam_models[cam_idx];
    help_tooltip_on_hover(
        "Lens model COLMAP fits. OPENCV (default) suits most phone/camera "
        "footage; OPENCV_FISHEYE for action cams and fisheye lenses; "
        "PINHOLE only for pre-undistorted or synthetic images.");

    const char* matchers[] = {"Auto", "Exhaustive", "Sequential"};
    ImGui::SetNextItemWidth(180);
    ImGui::Combo("Matcher", &_colmap_job.matcher, matchers, 3);
    help_tooltip_on_hover(
        "How image pairs are matched. Auto: sequential for video / large "
        "sets, exhaustive otherwise. Exhaustive is thorough but O(N^2).");

    ImGui::Checkbox("Single camera", &_colmap_job.single_camera);
    help_tooltip_on_hover("All images share one set of lens parameters. "
                          "Keep enabled when everything was shot with the "
                          "same camera and zoom.");
    ImGui::SameLine();
    ImGui::Checkbox("Use GPU", &_colmap_job.use_gpu);

    if (_colmap_job.is_video) {
        ImGui::SetNextItemWidth(180);
        ImGui::InputFloat("Frames per second", &_colmap_job.video_fps, 0, 0, "%.2g");
        help_tooltip_on_hover("How many frames to extract per second of "
                              "video. 1-3 fps is typical for a slow walkthrough.");
    }

    if (ImGui::CollapsingHeader("Tool locations")) {
        bool ch = false;
        ImGui::SetNextItemWidth(300);
        ch |= ImGui::InputText("colmap executable", &_colmap_exe);
        ImGui::SetNextItemWidth(300);
        ch |= ImGui::InputText("ffmpeg executable", &_ffmpeg_exe);
        if (ch) save_settings();
    }
    _colmap_job.colmap_exe = _colmap_exe;
    _colmap_job.ffmpeg_exe = _ffmpeg_exe;

    ImGui::EndDisabled();

    // ---- run / status ----
    ImGui::Spacing();
    if (!running) {
        bool ready = !_colmap_job.input_path.empty() && !_colmap_job.workspace.empty();
        ImGui::BeginDisabled(!ready);
        if (ImGui::Button("Run COLMAP", ImVec2(180, 34)))
            _colmap.start(_colmap_job);
        ImGui::EndDisabled();
    } else {
        if (ImGui::Button("Cancel", ImVec2(180, 34))) _colmap.cancel();
        ImGui::SameLine();
        ImGui::Text("%s ...", _colmap.stage().c_str());
    }

    switch (_colmap.state()) {
        case ColmapRunner::State::Done:
            ImGui::TextColored(kOk, "Done: %s", _colmap.dataset_dir().c_str());
            ImGui::SameLine();
            if (ImGui::Button("Open in Trainer"))
                open_dataset(_colmap.dataset_dir(), _colmap.image_dir());
            break;
        case ColmapRunner::State::Failed:
            ImGui::TextColored(kErr, "Failed: %s", _colmap.error().c_str());
            break;
        case ColmapRunner::State::Cancelled:
            ImGui::TextColored(kDim, "Cancelled.");
            break;
        default:
            break;
    }

    ImGui::Spacing();
    draw_log_panel(-1.0f);
}


// ===========================================================================
// Train screen
// ===========================================================================

void GuiApp::draw_train() {
    if (ImGui::Button("< Home")) _screen = Screen::Home;
    ImGui::SameLine();
    ImGui::TextDisabled("%s", _cfg.data.c_str());

    ImGui::BeginChild("##settings", ImVec2(420, 0), ImGuiChildFlags_Borders);
    draw_train_settings();
    ImGui::EndChild();

    ImGui::SameLine();
    ImGui::BeginGroup();
    float log_h = _show_log ? 150.0f : 0.0f;
    float status_h = ImGui::GetFrameHeightWithSpacing() +
                     ImGui::GetTextLineHeightWithSpacing() + 10;
    float spacing = ImGui::GetStyle().ItemSpacing.y;
    float vp_h = -(log_h + (log_h > 0 ? spacing : 0) + status_h + spacing);
    ImGui::BeginChild("##viewport", ImVec2(0, vp_h), ImGuiChildFlags_Borders);
    _viewport.draw(_runner.phase() == TrainRunner::Phase::Training);
    ImGui::EndChild();
    draw_status_strip();
    if (_show_log) draw_log_panel(log_h);
    ImGui::EndGroup();
}

void GuiApp::draw_train_settings() {
    TrainRunner::Phase ph = _runner.phase();
    bool busy = ph == TrainRunner::Phase::Loading ||
                ph == TrainRunner::Phase::Preparing ||
                ph == TrainRunner::Phase::Training;

    // ---- dataset ----
    ImGui::SeparatorText("Dataset");
    switch (ph) {
        case TrainRunner::Phase::Loading:
            ImGui::TextColored(kDim, "Parsing dataset ...");
            break;
        case TrainRunner::Phase::LoadError:
            ImGui::PushTextWrapPos();
            ImGui::TextColored(kErr, "%s", _runner.error().c_str());
            ImGui::PopTextWrapPos();
            break;
        default:
            if (auto* s = _runner.session()) {
                if (ph == TrainRunner::Phase::Ready ||
                    ph == TrainRunner::Phase::Training ||
                    ph == TrainRunner::Phase::Done) {
                    ImGui::Text("%lld cameras (%lld views) - %s points",
                                (long long)s->ds.num_cameras,
                                (long long)s->post.n_post,
                                format_count((double)s->ds.points.num()).c_str());
                }
            } else {
                ImGui::TextColored(kDim, "no dataset loaded");
            }
            break;
    }
    ImGui::BeginDisabled(busy && ph != TrainRunner::Phase::Loading);
    if (ImGui::Button("Change...")) {
        _pick = PickAction::OpenDataset;
        _dialog.open("Open Dataset Folder", FileDialog::Mode::Folder);
    }
    ImGui::EndDisabled();

    // ---- preset + options ----
    ImGui::BeginDisabled(busy);
    ImGui::SeparatorText("Preset");
    ImGui::SetNextItemWidth(-8);
    if (ImGui::BeginCombo("##preset", _preset.c_str())) {
        for (const auto& p : kSsplatPresets) {
            bool sel = _preset == p.name;
            if (ImGui::Selectable(p.name, sel)) apply_preset(p.name);
            if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort) &&
                ImGui::BeginTooltip()) {
                ImGui::PushTextWrapPos(360);
                ImGui::TextUnformatted(p.help);
                ImGui::PopTextWrapPos();
                ImGui::EndTooltip();
            }
            if (sel) ImGui::SetItemDefaultFocus();
        }
        ImGui::EndCombo();
    }
    ImGui::PushTextWrapPos();
    ImGui::TextColored(kDim, "%s", preset_help(_preset));
    ImGui::PopTextWrapPos();

    ImGui::SeparatorText("Basic Options");
    draw_basic_options();

    ImGui::Spacing();
    if (ImGui::CollapsingHeader("All Options (Advanced)"))
        draw_config_editor(_cfg, _defaults, _cfg_ui);
    ImGui::EndDisabled();

    // ---- controls + metrics ----
    ImGui::SeparatorText("Training");
    draw_train_controls();
    draw_metrics();
}

void GuiApp::draw_basic_options() {
    const float w = 170.0f;

    ImGui::SetNextItemWidth(w);
    ImGui::InputInt("Training steps", &_cfg.num_iterations, 0, 0);
    help_tooltip_on_hover(
        "How long to optimize. 30000 is a solid default; small scenes can "
        "look good at 10000-15000, and quality saturates beyond ~30000.");

    ImGui::SetNextItemWidth(w);
    ImGui::InputInt("Max splats", &_cfg.cap_max, 0, 0);
    help_tooltip_on_hover(
        "Upper bound on the number of Gaussians. More captures more detail "
        "but uses more VRAM and renders slower. ~1M suits most scenes; "
        "large outdoor scenes may want 2-4M.");

    int ds_idx = _cfg.rescale_camera_to_fit == 2.0f ? 1
               : _cfg.rescale_camera_to_fit == 4.0f ? 2
               : _cfg.rescale_camera_to_fit == 8.0f ? 3 : 0;
    const char* ds_items[] = {"native", "1/2", "1/4", "1/8"};
    ImGui::SetNextItemWidth(w);
    if (ImGui::Combo("Image resolution", &ds_idx, ds_items, 4)) {
        const float vals[] = {0.0f, 2.0f, 4.0f, 8.0f};
        _cfg.rescale_camera_to_fit = vals[ds_idx];
    }
    help_tooltip_on_hover(
        "Train at a fraction of the input resolution. Downscaling trains "
        "much faster and saves VRAM; use it for 4K+ footage or quick "
        "previews. (Also matches pre-downscaled images_2 / images_4 "
        "folders of academic datasets.)");

    ImGui::SetNextItemWidth(w);
    ImGui::SliderInt("Color detail (SH)", &_cfg.sh_degree, 0, 4);
    help_tooltip_on_hover(
        "Spherical-harmonics degree for view-dependent color (reflections, "
        "highlights). 3 is standard; 0 gives flat colors and the smallest "
        "model.");

    ImGui::Checkbox("Split fisheye into pinhole views", &_cfg.warp_to_pinhole);
    help_tooltip_on_hover(
        "For fisheye / 360-camera footage: split each image into 5 "
        "undistorted views for training. Enabled by the 360-camera preset.");

    bool web = !_cfg.disable_viewer;
    if (ImGui::Checkbox("Also serve web viewer", &web))
        _cfg.disable_viewer = !web;
    help_tooltip_on_hover(
        "Additionally serve the browser-based viewer while training, e.g. "
        "to monitor from another device. The native viewport here works "
        "either way.");
    if (web) {
        ImGui::SameLine();
        ImGui::SetNextItemWidth(90);
        ImGui::InputInt("port", &_cfg.viewer_port, 0, 0);
    }

    ImGui::SetNextItemWidth(w);
    ImGui::InputText("Output folder", &_cfg.output_dir_prefix);
    help_tooltip_on_hover("Where run outputs (checkpoints, splat.ply, "
                          "config.json) are written.");
    ImGui::SetNextItemWidth(w);
    ImGui::InputText("Run name", &_cfg.output_dir_name);
    help_tooltip_on_hover("Subfolder name for this run. Leave empty for "
                          "<dataset>_<timestamp>.");
}

void GuiApp::draw_train_controls() {
    TrainRunner::Phase ph = _runner.phase();
    switch (ph) {
        case TrainRunner::Phase::Idle:
        case TrainRunner::Phase::Ready:
        case TrainRunner::Phase::LoadError:
        case TrainRunner::Phase::Done:
        case TrainRunner::Phase::TrainError: {
            if (ph == TrainRunner::Phase::TrainError) {
                ImGui::PushTextWrapPos();
                ImGui::TextColored(kErr, "%s", _runner.error().c_str());
                ImGui::PopTextWrapPos();
            }
            if (ph == TrainRunner::Phase::Done) {
                ImGui::TextColored(kOk, "Training complete.");
                if (auto* s = _runner.session()) {
                    std::string out = s->out_dir.string();
                    ImGui::PushTextWrapPos();
                    ImGui::TextDisabled("Saved to %s", out.c_str());
                    ImGui::PopTextWrapPos();
                }
            }
            bool can_start = ph == TrainRunner::Phase::Ready ||
                             ph == TrainRunner::Phase::Done ||
                             ph == TrainRunner::Phase::TrainError;
            ImGui::BeginDisabled(!can_start);
            if (ImGui::Button(ph == TrainRunner::Phase::Done ? "Train Again"
                                                             : "Start Training",
                              ImVec2(-8, 36)))
                start_training();
            ImGui::EndDisabled();
            break;
        }
        case TrainRunner::Phase::Loading:
            ImGui::TextColored(kDim, "Parsing dataset ...");
            break;
        case TrainRunner::Phase::Preparing:
            ImGui::TextColored(kDim, "Preparing engine (seeding splats, "
                                     "caching images) ...");
            break;
        case TrainRunner::Phase::Training: {
            bool paused = _runner.paused();
            float half = (ImGui::GetContentRegionAvail().x - 12) * 0.5f;
            if (ImGui::Button(paused ? "Resume" : "Pause", ImVec2(half, 32)))
                _runner.set_paused(!paused);
            ImGui::SameLine();
            bool stopping = _runner.session() &&
                            _runner.session()->stop_requested.load();
            ImGui::BeginDisabled(stopping);
            if (ImGui::Button(stopping ? "Stopping..." : "Stop && Save",
                              ImVec2(half, 32)))
                _runner.request_stop();
            ImGui::EndDisabled();
            help_tooltip_on_hover("Finish the current step, save a "
                                  "checkpoint, and keep the result loaded "
                                  "for viewing.");
            break;
        }
    }
}

void GuiApp::draw_metrics() {
    static std::vector<TrainRunner::MetricPoint> pts;
    _runner.get_metrics(pts);
    if (pts.empty()) return;

    ImGui::SeparatorText("Metrics");
    // Downsample to <= 240 plot points.
    int stride = std::max<size_t>(1, pts.size() / 240);
    static std::vector<float> psnr, splats;
    psnr.clear();
    splats.clear();
    for (size_t i = 0; i < pts.size(); i += stride) {
        psnr.push_back(pts[i].psnr);
        splats.push_back(pts[i].num_splats);
    }
    const auto& last = pts.back();
    char overlay[64];
    std::snprintf(overlay, sizeof overlay, "PSNR %.2f", last.psnr);
    ImGui::PlotLines("##psnr", psnr.data(), (int)psnr.size(), 0, overlay,
                     FLT_MAX, FLT_MAX, ImVec2(-8, 64));
    ImGui::Text("splats: %s   ssim: %.3f   loss: %.4f",
                format_count(last.num_splats).c_str(), last.ssim, last.rgb_loss);
}

void GuiApp::draw_status_strip() {
    TrainRunner::Phase ph = _runner.phase();
    ssplat::TrainerProgress p = _runner.latest_progress();
    if (ph == TrainRunner::Phase::Training && p.total_steps > 0) {
        float frac = (float)(p.step + 1) / (float)p.total_steps;
        char label[64];
        std::snprintf(label, sizeof label, "step %d / %d", p.step + 1, p.total_steps);
        ImGui::ProgressBar(frac, ImVec2(-8, 0), label);
        bool paused = _runner.paused();
        ImGui::Text("%s%.0f ms/step   ETA %s   %s splats",
                    paused ? "[paused]  " : "",
                    p.step_latency * 1000.0,
                    format_eta(_runner.eta_seconds()).c_str(),
                    format_count((double)p.num_splats).c_str());
    } else if (ph == TrainRunner::Phase::Done && p.total_steps > 0) {
        char label[64];
        std::snprintf(label, sizeof label, "done (%d steps)", p.step + 1);
        ImGui::ProgressBar(1.0f, ImVec2(-8, 0), label);
        ImGui::TextDisabled("explore the result in the viewport above");
    } else if (ph == TrainRunner::Phase::Preparing) {
        ImGui::ProgressBar(-(float)ImGui::GetTime(), ImVec2(-8, 0), "preparing");
        ImGui::TextDisabled(" ");
    } else {
        ImGui::ProgressBar(0.0f, ImVec2(-8, 0), "idle");
        ImGui::TextDisabled(" ");
    }
}

void GuiApp::draw_log_panel(float height) {
    ImGui::BeginChild("##log", ImVec2(0, height), ImGuiChildFlags_Borders);
    for (const auto& s : _log) ImGui::TextUnformatted(s.c_str());
    if (_log_autoscroll && ImGui::GetScrollY() >= ImGui::GetScrollMaxY() - 4)
        ImGui::SetScrollHereY(1.0f);
    ImGui::EndChild();
}

void GuiApp::draw_exit_confirm() {
    if (_confirm_exit) {
        ImGui::OpenPopup("Exit?");
        _confirm_exit = false;
    }
    if (ImGui::BeginPopupModal("Exit?", nullptr,
                               ImGuiWindowFlags_AlwaysAutoResize)) {
        ImGui::TextUnformatted("Training is in progress.");
        ImGui::TextUnformatted("Stop training (a final checkpoint is saved) "
                               "and exit?");
        ImGui::Spacing();
        if (ImGui::Button("Stop && Exit", ImVec2(140, 0))) {
            _runner.request_stop();
            _quit = true;
            ImGui::CloseCurrentPopup();
        }
        ImGui::SameLine();
        if (ImGui::Button("Keep Training", ImVec2(140, 0)))
            ImGui::CloseCurrentPopup();
        ImGui::EndPopup();
    }
}

}  // namespace gui
