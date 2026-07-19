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
#include <ctime>
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

// True when two configs parse to the same dataset: every dataparser-group
// field (generated) plus the non-dataparser fields load_dataset() consumes.
bool parse_settings_equal(const SsplatConfig& a, const SsplatConfig& b) {
    bool eq = true;
#define SSPLAT_CMP(member, cli_key, pyname, group, choices, help)              \
    if (!std::strcmp(group, "dataparser")) eq = eq && (a.member == b.member);
    SSPLAT_CONFIG_FIELDS(SSPLAT_CMP)
#undef SSPLAT_CMP
    return eq && a.data == b.data &&
           a.warp_to_pinhole == b.warp_to_pinhole &&
           a.warp_spherical_to_pinhole == b.warp_spherical_to_pinhole &&
           a.load_depths == b.load_depths &&
           a.load_normals == b.load_normals &&
           a.relative_scale == b.relative_scale &&
           a.auto_scale_poses == b.auto_scale_poses;
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
    _colmap_job.python_exe = _python_exe;
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
        if (k == "recent" && !v.empty() &&
            std::find(_recents.begin(), _recents.end(), v) == _recents.end())
            _recents.push_back(v);
        else if (k == "colmap_exe" && !v.empty()) _colmap_exe = v;
        else if (k == "ffmpeg_exe" && !v.empty()) _ffmpeg_exe = v;
        else if (k == "python_exe" && !v.empty()) _python_exe = v;
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
    std::fprintf(f, "python_exe=%s\n", _python_exe.c_str());
    std::fclose(f);
}

void GuiApp::add_recent(std::string path) {
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

bool GuiApp::training_busy() const {
    TrainRunner::Phase ph = _runner.phase();
    return ph == TrainRunner::Phase::Training ||
           ph == TrainRunner::Phase::Preparing;
}

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

void GuiApp::open_dataset(std::string dir, std::string image_dir) {
    if (dir.empty()) return;
    _cfg.data = dir;
    // image_dir: the COLMAP runner hands its (possibly external) image dir
    // over in-memory right after a run; otherwise the dataparser default
    // applies and the user can set data.image_dir under Advanced.
    _cfg.image_dir = !image_dir.empty() ? image_dir : _defaults.image_dir;
    // Default the output next to the dataset -- much easier to find than a
    // CWD-relative "outputs" for someone who launched from a desktop icon.
    // Follows the dataset unless the user customized it (i.e. it still
    // matches the previous auto default).
    if (_cfg.output_dir_prefix == "outputs" ||
        _cfg.output_dir_prefix.empty() ||
        _cfg.output_dir_prefix == _defaults.output_dir_prefix)
        _cfg.output_dir_prefix = (fs::path(dir) / "outputs").string();
    _defaults.data = _cfg.data;
    _defaults.image_dir = _cfg.image_dir;
    _defaults.output_dir_prefix = _cfg.output_dir_prefix;
    add_recent(dir);
    save_settings();
    _viewport.detach();
    _runner.load_dataset(_cfg, _preset);
    _screen = Screen::Train;
}

void GuiApp::request_open_dataset(std::string dir) {
    if (training_busy()) {
        _pending = Pending::OpenDataset;
        _pending_path = dir;
        _open_confirm = true;
        return;
    }
    open_dataset(dir);
}

void GuiApp::request_go_home() {
    if (training_busy()) {
        _pending = Pending::GoHome;
        _open_confirm = true;
        return;
    }
    _screen = Screen::Home;
}

void GuiApp::start_training() {
    if (_cfg.data.empty()) return;
    _viewport.detach();
    // Engine setup initializes the backend on the selected device; from
    // here on the device combo is display-only (one device per process).
    _device_locked = true;
    _runner.start_training(_cfg, _preset);
}

void GuiApp::request_close() {
    if (training_busy()) {
        _pending = Pending::Quit;
        _open_confirm = true;
        return;
    }
    _quit = true;
}

void GuiApp::run_pending_if_stopped() {
    if (_pending == Pending::None || training_busy()) return;
    // Only fire once the user confirmed the stop (or training was never
    // busy); a dismissed modal must not leave a delayed action armed.
    if (!_stop_confirmed) { _pending = Pending::None; return; }
    _stop_confirmed = false;
    Pending p = _pending;
    _pending = Pending::None;
    switch (p) {
        case Pending::GoHome:     _screen = Screen::Home; break;
        case Pending::OpenDataset: open_dataset(_pending_path); break;
        case Pending::Quit:       _quit = true; break;
        default: break;
    }
}

static bool is_image_ext(const fs::path& p) {
    std::string e = p.extension().string();
    for (auto& c : e) c = (char)std::tolower((unsigned char)c);
    return e == ".jpg" || e == ".jpeg" || e == ".png" || e == ".webp" ||
           e == ".tif" || e == ".tiff" || e == ".bmp";
}

void GuiApp::handle_drop(const std::string& path) {
    std::error_code ec;
    fs::path p(path);
    if (fs::is_directory(p, ec)) {
        // SfM dataset markers first (a processed dataset usually also
        // contains an images folder -- opening wins over re-running COLMAP).
        bool dataset = fs::exists(p / "transforms.json", ec) ||
                       fs::is_directory(p / "sparse", ec) ||
                       fs::is_directory(p / "colmap", ec);
        if (!dataset) {
            // Metashape export: a camera .xml next to a point-cloud .ply
            // (what MetashapeParser probes for).
            bool has_xml = false, has_ply = false;
            for (fs::directory_iterator it(p, ec), end; !ec && it != end;
                 it.increment(ec)) {
                if (!it->is_regular_file(ec)) continue;
                std::string e = it->path().extension().string();
                for (auto& c : e) c = (char)std::tolower((unsigned char)c);
                has_xml = has_xml || e == ".xml";
                has_ply = has_ply || e == ".ply";
            }
            dataset = has_xml && has_ply;
        }
        if (dataset) {
            request_open_dataset(path);
            return;
        }
        bool has_image = false;
        for (fs::recursive_directory_iterator it(
                 p, fs::directory_options::skip_permission_denied, ec), end;
             !ec && it != end; it.increment(ec))
            if (it->is_regular_file(ec) && is_image_ext(it->path())) {
                has_image = true;
                break;
            }
        if (has_image) {
            if (training_busy()) {
                log("Dropped photo folder ignored: stop training first");
                return;
            }
            _pick = PickAction::ColmapImages;
            handle_dialog_result(path);
            return;
        }
        log("Dropped folder contains no dataset or images: " + path);
        return;
    }
    if (fs::is_regular_file(p, ec)) {
        std::string ext = p.extension().string();
        for (auto& c : ext) c = (char)std::tolower((unsigned char)c);
        static const char* kVideoExt[] = {".mp4", ".mov", ".avi", ".mkv",
                                          ".insv", ".webm", ".mts", ".m2ts",
                                          ".360", ".ts", ".wmv"};
        for (const char* v : kVideoExt)
            if (ext == v) {
                if (training_busy()) {
                    log("Dropped video ignored: stop training first");
                    return;
                }
                _pick = PickAction::ColmapVideo;
                handle_dialog_result(path);
                return;
            }
        // A file from inside a dataset (transforms.json, database.db, a
        // COLMAP .bin/.txt, a Metashape camera .xml) opens the dataset it
        // belongs to.
        if (p.filename() == "transforms.json" || ext == ".db" ||
            ext == ".bin" || ext == ".txt" || ext == ".xml") {
            request_open_dataset(p.parent_path().string());
            return;
        }
        log("Unsupported dropped file: " + path);
    }
}

// Default COLMAP workspace: never point at an existing non-empty directory
// (e.g. a previous run) -- append _2, _3, ... instead of overwriting.
static std::string fresh_workspace(const std::string& base) {
    std::error_code ec;
    if (!fs::exists(base, ec) || fs::is_empty(base, ec)) return base;
    for (int i = 2; i < 1000; i++) {
        std::string cand = base + "_" + std::to_string(i);
        if (!fs::exists(cand, ec) || fs::is_empty(cand, ec)) return cand;
    }
    return base;
}

void GuiApp::handle_dialog_result(const std::string& path) {
    switch (_pick) {
        case PickAction::OpenDataset:
            request_open_dataset(path);
            break;
        case PickAction::ColmapImages:
            _colmap_job.input_path = path;
            _colmap_job.is_video = false;
            _colmap_job.workspace = fresh_workspace(path + "_dataset");
            _colmap_job.matcher = 1;   // photos default to exhaustive
            _screen = Screen::Colmap;
            break;
        case PickAction::ColmapVideo: {
            _colmap_job.input_path = path;
            _colmap_job.is_video = true;
            fs::path p(path);
            _colmap_job.workspace = fresh_workspace(
                (p.parent_path() / (p.stem().string() + "_dataset")).string());
            _colmap_job.matcher = 2;   // frames are ordered: sequential
            // 360-camera preset for Insta360 .insv files: dual-fisheye
            // tracks land in one folder per lens (one COLMAP camera each),
            // THIN_PRISM_FISHEYE fits the lenses, and the known Insta360
            // focal length (fx = fy ~ 0.269 * width on every X5 dataset
            // measured) makes fisheye mapper initialization reliable.
            // Exhaustive matching, not sequential: the two lens tracks are
            // concatenated (not temporally interleaved) so sequential
            // neighbors miss the cross-lens pairs -- on a real X5 capture
            // exhaustive + model merge registered 116/118 frames where
            // sequential + loop detection topped out at 68 (frame counts
            // at 0.5-2 fps stay well within exhaustive range).
            std::string ext = p.extension().string();
            for (auto& c : ext) c = (char)std::tolower((unsigned char)c);
            if (ext == ".insv") {
                _colmap_job.camera_model = "THIN_PRISM_FISHEYE";
                _colmap_job.camera_mode = 1;
                _colmap_job.init_focal_factor = 0.269f;
                _colmap_job.matcher = 1;
                _colmap_job.seq_loop_closure = true;   // if switched to sequential
            }
            _screen = Screen::Colmap;
            break;
        }
        case PickAction::ColmapWorkspace:
            _colmap_job.workspace = path;
            break;
        case PickAction::OutputPrefix:
            _cfg.output_dir_prefix = path;
            break;
        case PickAction::VocabTree:
            _colmap_job.vocab_tree_path = path;
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
    run_pending_if_stopped();

    // Dataset-parsing option changed: re-parse once the user finishes
    // editing (no active widget), so the preview / camera counts stay in
    // sync with the config. Never fires while training (options are locked).
    if (_parse_dirty && !training_busy() &&
        _runner.phase() != TrainRunner::Phase::Loading &&
        !ImGui::IsAnyItemActive()) {
        _parse_dirty = false;
        if (!_cfg.data.empty()) {
            log("Dataset settings changed; reloading dataset");
            _viewport.detach();
            _runner.load_dataset(_cfg, _preset);
        }
    }

    // Viewport backend transitions: GL dataset preview once parsed, engine
    // renderer once training set the engine up.
    if (_runner.engine_ready()) {
        if (!_viewport.attached() && _runner.session())
            _viewport.attach(*_runner.session());
    } else if (_runner.phase() == TrainRunner::Phase::Ready) {
        if (!_viewport.preview_active() && _runner.session())
            _viewport.attach_preview(*_runner.session());
    }

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
    draw_confirm_modal();

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
            ImGui::TextUnformatted("Spirulae Splat");
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
    ImGui::TextUnformatted("Spirulae Splat");
    ImGui::SetWindowFontScale(1.0f);
    ImGui::TextDisabled("Reconstruct 3D scenes from photos with Gaussian splatting.");
    ImGui::Dummy(ImVec2(0, 24));

    // A session (possibly still training) exists -- offer the way back.
    if (_runner.phase() != TrainRunner::Phase::Idle) {
        if (ImGui::Button(training_busy() ? "Back to Training"
                                          : "Back to Trainer", ImVec2(-1, 42)))
            _screen = Screen::Train;
        ImGui::Dummy(ImVec2(0, 10));
    }

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
        "Pick a folder of overlapping photos of a scene or object "
        "(subfolders are included). COLMAP (installed separately) estimates "
        "the camera poses; the result opens directly in the trainer.");

    if (ImGui::Button("Create Dataset from Video...", ImVec2(-1, 42))) {
        _pick = PickAction::ColmapVideo;
        _dialog.open("Select Video File", FileDialog::Mode::File,
                     {".mp4", ".mov", ".avi", ".mkv", ".webm", ".m4v", ".insv"});
    }
    help_tooltip_on_hover(
        "Pick a video walking around a scene or object. The least blurry "
        "frames are extracted with ffmpeg, then processed with COLMAP.");
    ImGui::Spacing();
    ImGui::TextDisabled("...or drop a dataset folder, photo folder, or "
                        "video file anywhere in this window");

    if (!_recents.empty()) {
        ImGui::Dummy(ImVec2(0, 18));
        ImGui::SeparatorText("Recent");
        for (size_t i = 0; i < _recents.size(); i++) {
            ImGui::PushID((int)i);
            if (ImGui::Selectable(_recents[i].c_str()))
                request_open_dataset(_recents[i]);
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

void GuiApp::draw_colmap_options() {
    const char* qualities[] = {"Fast", "Balanced", "High quality"};
    ImGui::SetNextItemWidth(180);
    ImGui::Combo("Quality", &_colmap_job.quality, qualities, 3);
    help_tooltip_on_hover("Feature count used for matching (4k / 8k / 16k). "
                          "Higher finds more cameras in difficult scenes but "
                          "matching is O(n^2) in feature count.");

    int cam_idx = 0;
    for (int i = 0; i < kNumColmapCameraModels; i++)
        if (_colmap_job.camera_model == kColmapCameraModels[i]) cam_idx = i;
    ImGui::SetNextItemWidth(180);
    if (ImGui::Combo("Camera model", &cam_idx, kColmapCameraModels,
                     kNumColmapCameraModels))
        _colmap_job.camera_model = kColmapCameraModels[cam_idx];
    help_tooltip_on_hover(
        "Lens model COLMAP fits (all models the trainer's parser supports "
        "are listed). OPENCV (default) suits most phone/camera footage; "
        "OPENCV_FISHEYE or THIN_PRISM_FISHEYE for action cams and fisheye "
        "lenses; PINHOLE only for pre-undistorted or synthetic images.");

    const char* cam_modes[] = {"one shared camera", "one camera per folder",
                               "one camera per image"};
    ImGui::SetNextItemWidth(180);
    ImGui::Combo("Camera sharing", &_colmap_job.camera_mode, cam_modes, 3);
    help_tooltip_on_hover(
        "How lens parameters are shared. \"One shared camera\" when "
        "everything was shot with the same camera and zoom. \"Per folder\" "
        "for multi-camera rigs organized as one subfolder per camera "
        "(multi-track 360 videos switch to this automatically). \"Per "
        "image\" when zoom/focus varied between shots.");

    const char* features[] = {"SIFT (classic)", "ALIKED (neural)"};
    ImGui::SetNextItemWidth(180);
    if (ImGui::Combo("Features", &_colmap_job.feature_type, features, 2))
        _colmap_job.lightglue = _colmap_job.feature_type == 1;
    help_tooltip_on_hover(
        "Keypoint detector/descriptor. SIFT is the robust default. ALIKED "
        "(learned, via ONNX; models auto-download) finds strong matches in "
        "low-texture or wide-baseline footage with far fewer keypoints, "
        "pairs with LightGlue matching, but does not support the "
        "vocabulary tree (tree indexes SIFT descriptors).");

    const char* matchers[] = {"Exhaustive", "Sequential", "Vocabulary tree"};
    int matcher_idx = _colmap_job.matcher - 1;
    if (matcher_idx < 0 || matcher_idx > 2)
        matcher_idx = _colmap_job.is_video ? 1 : 0;
    ImGui::SetNextItemWidth(180);
    ImGui::Combo("Matcher", &matcher_idx, matchers, 3);
    _colmap_job.matcher = matcher_idx + 1;
    help_tooltip_on_hover(
        "How image pairs are matched. Exhaustive tries every pair (best "
        "quality, O(n^2) -- fine up to a few hundred images). Sequential "
        "matches temporal neighbors (video). Vocabulary tree scales to "
        "thousands of unordered photos (the tree file is found next to the "
        "dataset or downloaded automatically, one-time).");
    if (_colmap_job.matcher == 2) {
        ImGui::Indent();
        ImGui::Checkbox("Loop closure detection", &_colmap_job.seq_loop_closure);
        help_tooltip_on_hover(
            "SequentialMatching.loop_detection: retrieve visually similar "
            "non-neighbor frames via the vocabulary tree and match them, "
            "closing loops when the camera revisits a spot. SIFT features "
            "only.");
        ImGui::Unindent();
    }

    if (_colmap_job.is_video) {
        ImGui::SetNextItemWidth(180);
        ImGui::InputFloat("Frames per second", &_colmap_job.video_fps, 0, 0, "%.2g");
        help_tooltip_on_hover("How many frames to keep per second of video. "
                              "1-3 fps is typical for a slow walkthrough.");
        ImGui::SetNextItemWidth(180);
        ImGui::SliderInt("Sharpness window", &_colmap_job.sharp_window, 1, 8);
        help_tooltip_on_hover(
            "Extract this many candidate frames per kept frame and keep the "
            "least motion-blurred one (variance-of-Laplacian, multithreaded). "
            "1 disables blur-aware selection.");
    }

    // ---- masking ----
    ImGui::Checkbox("Mask out objects (AI segmentation)", &_colmap_job.mask_enable);
    help_tooltip_on_hover(
        "Generate masks from text prompts with lang-segment-anything "
        "(external Python; scripts/mask.py is bundled). Masked regions are "
        "ignored by COLMAP feature extraction AND by training -- use it to "
        "remove people, cars, or the camera operator.");
    if (_colmap_job.mask_enable) {
        ImGui::Indent();
        ImGui::SetNextItemWidth(280);
        ImGui::InputTextWithHint("Mask prompt", "people; cars; shadow",
                                 &_colmap_job.mask_prompt);
        help_tooltip_on_hover("Semicolon-separated objects to mask out.");
        ImGui::SetNextItemWidth(280);
        ImGui::InputTextWithHint("Keep prompt", "person in painting",
                                 &_colmap_job.mask_negative_prompt);
        help_tooltip_on_hover("Objects matching these stay UNmasked even if "
                              "they also match the mask prompt. Optional.");
        const char* mask_models[] = {"sam2.1_hiera_large", "sam3"};
        int mi = _colmap_job.mask_model == "sam3" ? 1 : 0;
        ImGui::SetNextItemWidth(180);
        if (ImGui::Combo("Segmentation model", &mi, mask_models, 2))
            _colmap_job.mask_model = mask_models[mi];
        help_tooltip_on_hover("sam2.1 via lang-segment-anything (default); "
                              "SAM-3 is higher quality but needs the sam3 "
                              "package and model access.");
        ImGui::Unindent();
    }

    // ---- advanced ----
    if (ImGui::CollapsingHeader("Advanced")) {
        bool fisheye = _colmap_job.camera_model.find("FISHEYE") != std::string::npos;
        ImGui::SetNextItemWidth(180);
        ImGui::InputFloat("Initial focal length (x width, 0 = unknown)",
                          &_colmap_job.init_focal_factor, 0, 0, "%.4g");
        help_tooltip_on_hover(
            "Seed COLMAP with fx = fy = factor * image width (principal "
            "point centered, zero distortion) instead of its generic guess. "
            "A known focal length stabilizes mapper initialization a lot, "
            "especially for fisheye lenses. Insta360 X5: ~0.269 (set "
            "automatically for .insv input).");
        ImGui::SetNextItemWidth(280);
        ImGui::InputTextWithHint("Initial camera params",
                                 "fx,fy,cx,cy,... (overrides focal length)",
                                 &_colmap_job.camera_params);
        help_tooltip_on_hover(
            "Raw ImageReader.camera_params for the selected camera model "
            "(full calibration prior). Leave empty to use the focal-length "
            "factor above, or both empty for COLMAP's default "
            "initialization.");
        ImGui::SetNextItemWidth(180);
        ImGui::InputInt("Max features (0 = auto)", &_colmap_job.max_num_features, 0, 0);
        help_tooltip_on_hover("SiftExtraction / AlikedExtraction "
                              ".max_num_features; overrides the Quality "
                              "preset when non-zero.");
        ImGui::SetNextItemWidth(180);
        ImGui::InputInt("Max image size (0 = off)", &_colmap_job.max_image_size, 0, 0);
        help_tooltip_on_hover("FeatureExtraction.max_image_size: downscale "
                              "images beyond this for feature extraction.");
        ImGui::SetNextItemWidth(180);
        ImGui::InputInt("Sequential overlap", &_colmap_job.seq_overlap, 0, 0);
        help_tooltip_on_hover("How many neighboring frames each frame is "
                              "matched against (sequential matcher).");
        ImGui::Checkbox("Quadratic overlap", &_colmap_job.seq_quadratic_overlap);
        help_tooltip_on_hover("Additionally match frame i against frames "
                              "i +- 2^k (sequential matcher). Helps close "
                              "loops in longer captures; enabled by default.");
        if (_colmap_job.is_video) {
            ImGui::SetNextItemWidth(180);
            ImGui::InputInt("Max frames", &_colmap_job.max_frames, 0, 0);
            help_tooltip_on_hover("Upper bound on extracted frames; the "
                                  "sharpness window grows to stay under it.");
        }
        ImGui::Checkbox("LightGlue matching", &_colmap_job.lightglue);
        help_tooltip_on_hover("Neural feature matcher (FeatureMatching.type "
                              "*_LIGHTGLUE): more matches on hard pairs than "
                              "brute-force descriptor distance. Default for "
                              "ALIKED features; also works with SIFT.");
        if (_colmap_job.feature_type == 0) {
            ImGui::Checkbox("Affine SIFT + guided matching",
                            &_colmap_job.estimate_affine_shape);
            help_tooltip_on_hover("SiftExtraction.estimate_affine_shape + "
                                  "FeatureMatching.guided_matching: slower but "
                                  "more robust matching.");
        }
        const char* extra_modes[] = {"Auto", "During mapping", "Final pass only"};
        ImGui::SetNextItemWidth(180);
        ImGui::Combo("Distortion refinement", &_colmap_job.mapper_extra_params,
                     extra_modes, 3);
        help_tooltip_on_hover(
            "When distortion coefficients are optimized. \"Final pass "
            "only\" holds them fixed during mapping "
            "(Mapper.ba_refine_extra_params 0) -- more stable for "
            "low-distortion perspective lenses -- and recovers them in the "
            "final refinement pass. Auto: final-pass-only for perspective "
            "models, during mapping for fisheye.");
        ImGui::SetNextItemWidth(180);
        ImGui::InputInt("Min matches per pair (0 = default)",
                        &_colmap_job.min_num_matches, 0, 0);
        help_tooltip_on_hover("Mapper.min_num_matches (default 15): image "
                              "pairs with fewer inlier matches are ignored "
                              "by the mapper. Raise to suppress spurious "
                              "registrations, lower for sparse overlap.");
        ImGui::SeparatorText("Repetitive scenes");
        help_tooltip_on_hover(
            "Large scenes with repeating structure (several similar rooms, "
            "tiled facades) often weld physically different but "
            "similar-looking parts together. These make matching and "
            "registration stricter to suppress that; 0 = COLMAP default.");
        // Preset levels filling the five fields below (editing any field
        // afterwards shows "Custom"). Stricter = fewer wrong welds but
        // fewer registered images on genuinely weak overlap.
        struct RepLevel {
            const char* name;
            float ratio; int pair_in; int reg_in; float reg_ratio; float err;
        };
        static const RepLevel kRepLevels[] = {
            {"Off (COLMAP defaults)", 0.0f,    0,   0, 0.0f,  0.0f},
            {"Low",                   0.75f,  30,  40, 0.30f, 10.0f},
            {"Medium",                0.70f,  60,  60, 0.40f,  8.0f},
            {"High",                  0.62f, 100, 100, 0.50f,  6.0f},
        };
        int rep_idx = -1;
        for (int i = 0; i < 4; i++)
            if (_colmap_job.match_max_ratio == kRepLevels[i].ratio &&
                _colmap_job.min_inliers_per_pair == kRepLevels[i].pair_in &&
                _colmap_job.abs_pose_min_num_inliers == kRepLevels[i].reg_in &&
                _colmap_job.abs_pose_min_inlier_ratio == kRepLevels[i].reg_ratio &&
                _colmap_job.abs_pose_max_error == kRepLevels[i].err) {
                rep_idx = i;
                break;
            }
        ImGui::SetNextItemWidth(180);
        if (ImGui::BeginCombo("Repetitive level",
                              rep_idx < 0 ? "Custom" : kRepLevels[rep_idx].name)) {
            for (int i = 0; i < 4; i++)
                if (ImGui::Selectable(kRepLevels[i].name, rep_idx == i)) {
                    _colmap_job.match_max_ratio = kRepLevels[i].ratio;
                    _colmap_job.min_inliers_per_pair = kRepLevels[i].pair_in;
                    _colmap_job.abs_pose_min_num_inliers = kRepLevels[i].reg_in;
                    _colmap_job.abs_pose_min_inlier_ratio = kRepLevels[i].reg_ratio;
                    _colmap_job.abs_pose_max_error = kRepLevels[i].err;
                }
            ImGui::EndCombo();
        }
        help_tooltip_on_hover(
            "How aggressively wrong matches are suppressed; fills the "
            "fields below. Low: mild tightening, keeps registration rate. "
            "Medium: good first attempt for multi-room indoor captures. "
            "High: for heavy repetition (identical rooms/facades) -- "
            "expect fewer registered images if overlap is thin.");
        ImGui::SetNextItemWidth(180);
        ImGui::InputFloat("Match ratio test (0 = default 0.8)",
                          &_colmap_job.match_max_ratio, 0, 0, "%.3g");
        help_tooltip_on_hover(
            "SiftMatching.max_ratio, the Lowe ratio test: a feature match "
            "is kept only when its best match is this much better than the "
            "second best. LOWER is stricter -- try 0.6-0.7 when repetitive "
            "texture creates false matches. SIFT only.");
        ImGui::SetNextItemWidth(180);
        ImGui::InputInt("Min inliers per pair (0 = default 15)",
                        &_colmap_job.min_inliers_per_pair, 0, 0);
        help_tooltip_on_hover(
            "TwoViewGeometry.min_num_inliers: image pairs whose geometric "
            "verification finds fewer inliers are discarded outright. "
            "Raise to 50-100 so weakly-supported (usually false) links "
            "between similar-looking areas never enter the database.");
        ImGui::SetNextItemWidth(180);
        ImGui::InputInt("Min inliers to register (0 = default 30)",
                        &_colmap_job.abs_pose_min_num_inliers, 0, 0);
        help_tooltip_on_hover(
            "Mapper.abs_pose_min_num_inliers: minimum absolute-pose inliers "
            "to register an image into the model. Raise to 50-100 to stop "
            "images from registering onto the wrong (similar-looking) part "
            "of the scene.");
        ImGui::SetNextItemWidth(180);
        ImGui::InputFloat("Min inlier ratio to register (0 = default 0.25)",
                          &_colmap_job.abs_pose_min_inlier_ratio, 0, 0, "%.3g");
        help_tooltip_on_hover(
            "Mapper.abs_pose_min_inlier_ratio: minimum fraction of 2D-3D "
            "correspondences that must be pose inliers. Try 0.35-0.5 for "
            "stricter registration.");
        ImGui::SetNextItemWidth(180);
        ImGui::InputFloat("Max registration error px (0 = default 12)",
                          &_colmap_job.abs_pose_max_error, 0, 0, "%.3g");
        help_tooltip_on_hover(
            "Mapper.abs_pose_max_error: reprojection error threshold (px) "
            "for absolute-pose RANSAC when registering images. Lower "
            "(6-8) = stricter; combine with the inlier thresholds above.");
        ImGui::Separator();
        if (fisheye) ImGui::BeginDisabled();
        ImGui::Checkbox("GPU bundle adjustment", &_colmap_job.ba_use_gpu);
        if (fisheye) ImGui::EndDisabled();
        help_tooltip_on_hover(fisheye
            ? "Mapper.ba_use_gpu -- unavailable: COLMAP's GPU bundle "
              "adjustment does not support fisheye camera models yet."
            : "Mapper.ba_use_gpu.");
        ImGui::Checkbox("Merge partial models", &_colmap_job.merge_models);
        help_tooltip_on_hover(
            "When the mapper splits the scene into several partial models, "
            "try colmap model_merger to fuse them (kept only when the "
            "merged model registers more images). The trainer otherwise "
            "auto-picks the largest partial model.");
        ImGui::Checkbox("Final refinement pass", &_colmap_job.final_bundle_adjust);
        help_tooltip_on_hover("Run bundle_adjuster after mapping on the "
                              "largest (or merged) model, refining focal "
                              "length, principal point, and distortion "
                              "(as in scripts/run_colmap.bash).");
        ImGui::SetNextItemWidth(-160);
        ImGui::InputTextWithHint("##vocab", "vocabulary tree (auto find/download)",
                                 &_colmap_job.vocab_tree_path);
        ImGui::SameLine();
        if (ImGui::Button("Browse...##vt")) {
            _pick = PickAction::VocabTree;
            _dialog.open("Select Vocabulary Tree (.bin)", FileDialog::Mode::File,
                         {".bin"});
        }
        ImGui::SameLine();
        ImGui::TextUnformatted("vocab tree");
    }

    if (ImGui::CollapsingHeader("Tool locations")) {
        bool ch = false;
        ImGui::SetNextItemWidth(300);
        ch |= ImGui::InputText("colmap executable", &_colmap_exe);
        ImGui::SetNextItemWidth(300);
        ch |= ImGui::InputText("ffmpeg executable", &_ffmpeg_exe);
        ImGui::SetNextItemWidth(300);
        ch |= ImGui::InputText("python executable", &_python_exe);
        help_tooltip_on_hover("Python is only used for AI masking.");
        if (ch) save_settings();
    }
    _colmap_job.colmap_exe = _colmap_exe;
    _colmap_job.ffmpeg_exe = _ffmpeg_exe;
    _colmap_job.python_exe = _python_exe;
}

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

    // Resume: when the chosen workspace already holds a (partial) run,
    // completed stages can be reused instead of starting from scratch.
    {
        std::error_code ec;
        fs::path ws(_colmap_job.workspace);
        bool prior = !_colmap_job.workspace.empty() &&
                     (fs::exists(ws / "database.db", ec) ||
                      (fs::is_directory(ws / "images", ec) &&
                       !fs::is_empty(ws / "images", ec)) ||
                      fs::is_directory(ws / "sparse", ec));
        if (prior) {
            ImGui::Checkbox("Resume previous run", &_colmap_job.resume);
            help_tooltip_on_hover(
                "This folder contains a previous (possibly interrupted) "
                "run. When checked, completed work is reused: extracted "
                "video frames, masks, features and matches already in "
                "database.db, and finished sparse models -- only the "
                "missing stages run. Uncheck to require an empty folder "
                "instead (nothing is deleted automatically).");
            ImGui::SameLine();
            ImGui::TextDisabled("(previous run detected in this folder)");
        }
    }

    ImGui::Spacing();
    draw_colmap_options();

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
            if (ImGui::Button("Open in Trainer")) {
                if (training_busy()) {
                    _pending = Pending::OpenDataset;
                    _pending_path = _colmap.dataset_dir();
                    _open_confirm = true;
                } else {
                    open_dataset(_colmap.dataset_dir(), _colmap.image_dir());
                }
            }
            break;
        case ColmapRunner::State::Failed:
            ImGui::PushTextWrapPos();
            ImGui::TextColored(kErr, "Failed: %s", _colmap.error().c_str());
            ImGui::PopTextWrapPos();
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
    bool preparing = _runner.phase() == TrainRunner::Phase::Preparing;
    ImGui::BeginDisabled(preparing);   // no clean interruption point
    if (ImGui::Button("< Home")) request_go_home();
    ImGui::EndDisabled();
    if (training_busy())
        help_tooltip_on_hover("Training is in progress -- leaving stops it "
                              "(a checkpoint is saved first).");
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

    // ---- device ----
    ImGui::SeparatorText("Device");
    {
        int n_dev = backend::device_count();
        int cur = backend::device_current();
        backend::DeviceInfo curd = backend::device_info(cur);
        ImGui::BeginDisabled(_device_locked || busy || n_dev == 0);
        ImGui::SetNextItemWidth(-8);
        if (ImGui::BeginCombo("##device",
                              cur >= 0 ? curd.name : "(no device found)")) {
            for (int i = 0; i < n_dev; i++) {
                backend::DeviceInfo d = backend::device_info(i);
                char label[300];
                std::snprintf(label, sizeof(label), "%s (%s, %llu MB)%s##d%d",
                              d.name, d.type,
                              (unsigned long long)(d.vram_bytes >> 20),
                              d.usable ? "" : " [unsupported]", i);
                ImGui::BeginDisabled(!d.usable);
                bool sel = i == cur;
                if (ImGui::Selectable(label, sel)) backend::device_select(i);
                if (sel) ImGui::SetItemDefaultFocus();
                ImGui::EndDisabled();
            }
            ImGui::EndCombo();
        }
        ImGui::EndDisabled();
        if (_device_locked) {
            ImGui::PushTextWrapPos();
            ImGui::TextColored(kDim, "Device is fixed once training starts; "
                                     "restart the app to change it.");
            ImGui::PopTextWrapPos();
        }
    }

    // ---- preset + options ----
    // Snapshot: any change to a dataset-parsing option below triggers an
    // automatic reload (deferred until the edited widget loses focus).
    SsplatConfig parse_before = _cfg;
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

    if (!parse_settings_equal(parse_before, _cfg)) _parse_dirty = true;

    // ---- controls + metrics ----
    ImGui::SeparatorText("Training");
    draw_train_controls();
    draw_metrics();
}

void GuiApp::draw_basic_options() {
    const float w = 170.0f;

    // Output location first -- the thing every new user looks for.
    ImGui::SetNextItemWidth(w);
    ImGui::InputText("##outdir", &_cfg.output_dir_prefix);
    ImGui::SameLine();
    if (ImGui::Button("...##outdir")) {
        _pick = PickAction::OutputPrefix;
        _dialog.open("Select Output Folder", FileDialog::Mode::Folder,
                     {}, _cfg.output_dir_prefix);
    }
    ImGui::SameLine();
    ImGui::TextUnformatted("Output folder");
    help_tooltip_on_hover("Where run outputs (checkpoints, splat.ply, "
                          "config.json) are written. Each run gets its own "
                          "subfolder.");
    ImGui::SetNextItemWidth(w);
    ImGui::InputTextWithHint("Run name", "auto: <dataset>_<time>",
                             &_cfg.output_dir_name);
    help_tooltip_on_hover("Subfolder name for this run. Leave empty for "
                          "<dataset>_<timestamp>.");
    {
        std::string run = _cfg.output_dir_name.empty()
            ? fs::path(_cfg.data).stem().string() + "_<time>"
            : _cfg.output_dir_name;
        ImGui::PushTextWrapPos();
        ImGui::TextColored(kDim, "-> %s",
            (fs::path(_cfg.output_dir_prefix) / run).string().c_str());
        ImGui::PopTextWrapPos();
    }
    ImGui::Spacing();

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

    {
        static const char* prims[] = {"3dgs", "mip", "3dgut"};
        int pi = _cfg.primitive == "mip" ? 1 : _cfg.primitive == "3dgut" ? 2 : 0;
        ImGui::SetNextItemWidth(w);
        if (ImGui::Combo("Primitive", &pi, prims, 3))
            _cfg.primitive = prims[pi];
        help_tooltip_on_hover(
            "Splat primitive. 3dgs: standard 3D Gaussian splatting. mip: "
            "anti-aliased Mip-Splatting, reduces shimmering when zooming "
            "out. 3dgut: Unscented-Transform projection, exact for "
            "distorted (fisheye/wide) cameras; used by the meshing preset.");
    }

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
    } else if (ph == TrainRunner::Phase::Ready) {
        ImGui::ProgressBar(0.0f, ImVec2(-8, 0), "ready");
        ImGui::TextDisabled("dataset preview -- press Start Training when ready");
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

void GuiApp::draw_confirm_modal() {
    if (_open_confirm) {
        ImGui::OpenPopup("Stop training?");
        _open_confirm = false;
        _confirm_shown = true;
        _stop_confirmed = false;
    }
    if (ImGui::BeginPopupModal("Stop training?", nullptr,
                               ImGuiWindowFlags_AlwaysAutoResize)) {
        ImGui::TextUnformatted("Training is in progress.");
        const char* what = _pending == Pending::Quit ? "exit"
                         : _pending == Pending::GoHome ? "go to the home screen"
                         : "open the new dataset";
        ImGui::Text("Stop training (a final checkpoint is saved) and %s?", what);
        ImGui::Spacing();
        if (ImGui::Button("Stop && Save", ImVec2(150, 0))) {
            _runner.request_stop();
            _stop_confirmed = true;
            _confirm_shown = false;
            // run_pending_if_stopped() completes the action once the stop
            // lands (phase leaves Training).
            ImGui::CloseCurrentPopup();
        }
        ImGui::SameLine();
        if (ImGui::Button("Keep Training", ImVec2(150, 0))) {
            _pending = Pending::None;
            _pending_path.clear();
            _confirm_shown = false;
            ImGui::CloseCurrentPopup();
        }
        ImGui::EndPopup();
    } else if (_confirm_shown) {
        // Dismissed (Esc / click-away): treat as "keep training".
        _pending = Pending::None;
        _pending_path.clear();
        _confirm_shown = false;
    }
}

}  // namespace gui
