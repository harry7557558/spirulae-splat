// SegmentPanel.cpp -- see SegmentPanel.h.

#include "app/gui/SegmentPanel.h"

#include "app/gui/ConfigUI.h"   // help_tooltip_on_hover

#include "imgui.h"
#include "imgui_stdlib.h"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <filesystem>

#ifdef SSPLAT_BUILD_SAM
#include "nn/io/Image.h"
#include "sam/Masking.h"
#endif
#ifdef SSPLAT_HAVE_VIDEO
#include "video/Video.h"
#endif

namespace fs = std::filesystem;

namespace gui {

namespace {

bool is_image_file(const fs::path& p) {
    std::string e = p.extension().string();
    for (auto& c : e) c = (char)std::tolower((unsigned char)c);
    return e == ".jpg" || e == ".jpeg" || e == ".png" || e == ".webp" ||
           e == ".tif" || e == ".tiff" || e == ".bmp";
}

}  // namespace

// The GPU-side state, kept alive across runs so a prompt edit costs one
// forward pass rather than a 3-second weight upload.
struct SegmentPanel::Job {
#ifdef SSPLAT_BUILD_SAM
    sam::Masker masker;
    std::string loaded_model;       // what masker was initialized with
    std::string loaded_signature;   // and with which prompt settings
    std::string frame_key;          // the frame currently held
    nn::Image   frame;
#endif
};

SegmentPanel::SegmentPanel() = default;

SegmentPanel::~SegmentPanel() {
    _cancel = true;
    if (_worker.joinable()) _worker.join();
}

void SegmentPanel::open(const std::string& input, bool is_video,
                        const std::string& model_path) {
    _input = input;
    _is_video = is_video;
    _model_path = model_path;
    _frame_idx = 0;
    _clicks.clear();
    _frame_dirty = true;
    _needs_run = true;
    _kept_fraction = -1.0f;
    {
        std::lock_guard<std::mutex> lk(_mu);
        _status.clear();
        _error.clear();
    }
    collect_frames(input, is_video);
    _open = true;
}

void SegmentPanel::close() {
    _cancel = true;
    if (_worker.joinable()) _worker.join();
    _cancel = false;
    // Drop the model: the reconstruction that usually follows wants the VRAM.
    _job.reset();
    _open = false;
}

void SegmentPanel::destroy_gl() {
    if (_tex) {
        glDeleteTextures(1, &_tex);
        _tex = 0;
    }
}

void SegmentPanel::collect_frames(const std::string& input, bool is_video) {
    _frames.clear();
    std::error_code ec;
    if (is_video) {
        // Frames are decoded on demand; the list is just how many offers the
        // slider makes. Seeking is not supported by the decoder, so "frame N"
        // means "the Nth frame we sample while reading forward".
        for (int i = 0; i < 8; i++) _frames.push_back("");
        return;
    }
    for (fs::recursive_directory_iterator it(
             input, fs::directory_options::skip_permission_denied, ec), end;
         !ec && it != end; it.increment(ec))
        if (it->is_regular_file(ec) && is_image_file(it->path()))
            _frames.push_back(it->path().string());
    std::sort(_frames.begin(), _frames.end());
    // A dozen spread through the capture is enough to judge a prompt, and
    // keeps the slider meaningful on a 3000-photo folder.
    if (_frames.size() > 12) {
        std::vector<std::string> picked;
        for (int i = 0; i < 12; i++)
            picked.push_back(_frames[(size_t)i * (_frames.size() - 1) / 11]);
        _frames.swap(picked);
    }
}

// ---------------------------------------------------------------------------
// The worker
// ---------------------------------------------------------------------------

void SegmentPanel::start_job(const MaskSettings& s) {
#ifndef SSPLAT_BUILD_SAM
    (void)s;
    std::lock_guard<std::mutex> lk(_mu);
    _error = "this build has no built-in segmentation "
             "(-DSSPLAT_BUILD_SAM=OFF)";
#else
    if (_busy.load()) return;
    if (_worker.joinable()) _worker.join();
    if (_model_path.empty()) {
        std::lock_guard<std::mutex> lk(_mu);
        _error = "no model selected -- download one first";
        return;
    }

    const int idx = _frame_idx;
    const std::string frame_path =
        (idx >= 0 && idx < (int)_frames.size()) ? _frames[idx] : std::string();
    const std::string input = _input;
    const bool is_video = _is_video;
    const std::string model = _model_path;
    const std::vector<Click> clicks = _clicks;
    const MaskSettings settings = s;
    const bool frame_dirty = _frame_dirty;
    _frame_dirty = false;

    _busy = true;
    {
        std::lock_guard<std::mutex> lk(_mu);
        _error.clear();
        _status = "working...";
    }

    _worker = std::thread([this, settings, model, frame_path, input, is_video,
                           idx, clicks, frame_dirty] {
        auto set_error = [&](const std::string& e) {
            std::lock_guard<std::mutex> lk(_mu);
            _error = e;
            _status.clear();
        };
        try {
            if (!_job) _job = std::make_unique<Job>();
            Job& j = *_job;

            // ---- the frame ----
            const std::string key = is_video ? ("#" + std::to_string(idx))
                                             : frame_path;
            if (frame_dirty || j.frame_key != key || j.frame.empty()) {
                {
                    std::lock_guard<std::mutex> lk(_mu);
                    _status = "loading the frame...";
                }
                if (is_video) {
#ifdef SSPLAT_HAVE_VIDEO
                    // No seek: read forward and keep every (total/8)th frame.
                    video::VideoReader r;
                    if (!r.open(input)) return set_error(r.lastError());
                    const int total = std::max(r.info().frame_count, 8);
                    const int want = std::min(total - 1, idx * (total / 8));
                    nn::Image img;
                    for (int i = 0; i <= want; i++) {
                        if (_cancel.load()) return;
                        img = r.readFrame();
                        if (img.empty()) break;
                    }
                    if (img.empty())
                        return set_error("could not decode a frame from this "
                                         "video; try the ffmpeg path");
                    j.frame = std::move(img);
#else
                    return set_error("this build cannot decode video "
                                     "(-DSSPLAT_ENABLE_PATENTED=OFF); pick a "
                                     "folder of photos to preview, or run the "
                                     "job and check the masks it writes");
#endif
                } else {
                    if (frame_path.empty()) return set_error("no frame to show");
                    j.frame = nn::load_image(frame_path);
                    if (j.frame.empty())
                        return set_error("could not read " + frame_path);
                }
                j.frame_key = key;
            }
            if (_cancel.load()) return;

            // ---- the model ----
            // Every field below changes what the masker computes, so the
            // signature is what decides whether it can be reused. The weights
            // are the expensive part and only the model path moves them.
            char sig[512];
            std::snprintf(sig, sizeof sig, "%s|%s|%d|%d|%.4f|%zu",
                          settings.prompt.c_str(),
                          settings.negative_prompt.c_str(),
                          (int)settings.keep_subject, settings.max_image_size,
                          settings.threshold, clicks.size());
            if (j.loaded_model != model || j.loaded_signature != sig) {
                {
                    std::lock_guard<std::mutex> lk(_mu);
                    _status = j.loaded_model == model
                                  ? "preparing..."
                                  : "loading the model (first run is slow)...";
                }
                sam::MaskOptions mo;
                mo.model = model;
                mo.text = settings.prompt;
                mo.neg_text = settings.negative_prompt;
                mo.keep_prompted = settings.keep_subject;
                mo.max_size = settings.max_image_size;
                mo.threshold = settings.threshold;
                mo.video = false;      // one still frame, no memory bank
                for (const Click& c : clicks) {
                    if (c.positive) mo.seed.pos_points.push_back({c.x, c.y});
                    else            mo.seed.neg_points.push_back({c.x, c.y});
                }
                mo.has_seed = !mo.seed.pos_points.empty();
                std::string err;
                if (!j.masker.init(mo, err)) {
                    j.loaded_model.clear();
                    return set_error(err);
                }
                j.loaded_model = model;
                j.loaded_signature = sig;
            }
            if (_cancel.load()) return;

            // ---- run ----
            {
                std::lock_guard<std::mutex> lk(_mu);
                _status = "segmenting...";
            }
            sam::Mask mask;
            sam::Result detections;
            if (!j.masker.run(j.frame, mask, &detections))
                return set_error(j.masker.lastError());

            // ---- composite ----
            // Kept pixels stay as they are; masked-out pixels are dimmed and
            // tinted red, which reads as "this will be ignored" far better
            // than a separate black-and-white mask image next to the photo.
            const nn::Image& img = j.frame;
            std::vector<uint8_t> rgb = img.data;
            size_t kept = 0;
            const size_t n = (size_t)img.width * img.height;
            for (size_t i = 0; i < n && i * 3 + 2 < rgb.size(); i++) {
                if (i < mask.data.size() && mask.data[i] > 127) { kept++; continue; }
                rgb[i * 3 + 0] = (uint8_t)(rgb[i * 3 + 0] / 3 + 150);
                rgb[i * 3 + 1] = (uint8_t)(rgb[i * 3 + 1] / 3);
                rgb[i * 3 + 2] = (uint8_t)(rgb[i * 3 + 2] / 3);
            }
            {
                std::lock_guard<std::mutex> lk(_mu);
                _preview.swap(rgb);
                _preview_w = img.width;
                _preview_h = img.height;
                _preview_dirty = true;
                _kept_fraction = n ? (float)((double)kept / (double)n) : -1.0f;
                _status.clear();
                if (detections.detections.empty() && !settings.prompt.empty() &&
                    clicks.empty())
                    _error = "nothing matched this prompt on this frame";
            }
        } catch (const std::exception& e) {
            set_error(e.what());
        }
        _busy = false;
    });
#endif
}

// ---------------------------------------------------------------------------
// Drawing
// ---------------------------------------------------------------------------

void SegmentPanel::upload_preview() {
    std::vector<uint8_t> pixels;
    int w = 0, h = 0;
    {
        std::lock_guard<std::mutex> lk(_mu);
        if (!_preview_dirty) return;
        _preview_dirty = false;
        pixels = _preview;
        w = _preview_w;
        h = _preview_h;
    }
    if (pixels.empty() || w <= 0 || h <= 0) return;
    if (!_tex) glGenTextures(1, &_tex);
    glBindTexture(GL_TEXTURE_2D, _tex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, w, h, 0, GL_RGB, GL_UNSIGNED_BYTE,
                 pixels.data());
    _tex_w = w;
    _tex_h = h;
}

void SegmentPanel::draw_image(MaskSettings& settings) {
    const ImVec2 avail = ImGui::GetContentRegionAvail();
    const float box_h = std::max(avail.y - 8.0f, 120.0f);
    if (!_tex || _tex_w <= 0) {
        ImGui::Dummy(ImVec2(avail.x, box_h));
        return;
    }
    const float aspect = (float)_tex_w / (float)_tex_h;
    ImVec2 size(avail.x, avail.x / aspect);
    if (size.y > box_h) {
        size.y = box_h;
        size.x = box_h * aspect;
    }
    const ImVec2 origin = ImGui::GetCursorScreenPos();
    ImGui::Image((ImTextureID)(intptr_t)_tex, size);

    // Clicks land in source-image pixels, which is what the model wants.
    const bool hovered = ImGui::IsItemHovered();
    if (hovered && (ImGui::IsMouseClicked(ImGuiMouseButton_Left) ||
                    ImGui::IsMouseClicked(ImGuiMouseButton_Right))) {
        const ImVec2 m = ImGui::GetMousePos();
        Click c;
        c.x = (m.x - origin.x) / size.x * (float)_tex_w;
        c.y = (m.y - origin.y) / size.y * (float)_tex_h;
        c.positive = ImGui::IsMouseClicked(ImGuiMouseButton_Left);
        _clicks.push_back(c);
        start_job(settings);
    }
    if (hovered)
        ImGui::SetTooltip("Left-click: this is the object.  "
                          "Right-click: not this.");

    ImDrawList* dl = ImGui::GetWindowDrawList();
    for (const Click& c : _clicks) {
        const ImVec2 p(origin.x + c.x / (float)_tex_w * size.x,
                       origin.y + c.y / (float)_tex_h * size.y);
        const ImU32 col = c.positive ? IM_COL32(80, 220, 110, 255)
                                     : IM_COL32(240, 90, 90, 255);
        dl->AddCircleFilled(p, 5.0f, col);
        dl->AddCircle(p, 5.0f, IM_COL32(20, 20, 20, 200), 0, 1.5f);
    }
}

void SegmentPanel::draw(MaskSettings& settings) {
    if (!_open) return;

    upload_preview();

    const ImGuiViewport* vp = ImGui::GetMainViewport();
    ImGui::SetNextWindowSize(ImVec2(vp->WorkSize.x * 0.8f, vp->WorkSize.y * 0.8f),
                             ImGuiCond_Appearing);
    ImGui::SetNextWindowPos(vp->GetCenter(), ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
    bool open = true;
    if (!ImGui::Begin("Try the mask", &open, ImGuiWindowFlags_NoCollapse)) {
        ImGui::End();
        if (!open) close();
        return;
    }

    ImGui::TextDisabled(
        "Red = removed from the reconstruction. Adjust the prompt until only "
        "what you want gone is red.");
    ImGui::Separator();

    // ---- controls ----
    const float panel_w = 360.0f;
    ImGui::BeginChild("##segctl", ImVec2(panel_w, 0), ImGuiChildFlags_Borders);

    ImGui::TextUnformatted("What should be removed?");
    ImGui::SetNextItemWidth(-1);
    bool edited = ImGui::InputTextWithHint("##prompt", "people; cars; my shadow",
                                           &settings.prompt);
    help_tooltip_on_hover(
        "Plain words for the things to take out of the reconstruction, "
        "separated by semicolons. Anything that moved, reflected, or was not "
        "part of the scene is a good candidate.");

    ImGui::Spacing();
    ImGui::TextUnformatted("...but keep these");
    ImGui::SetNextItemWidth(-1);
    edited |= ImGui::InputTextWithHint("##negprompt", "person in a painting",
                                       &settings.negative_prompt);
    help_tooltip_on_hover(
        "Exceptions: things that match the prompt above but should stay. "
        "Optional.");

    ImGui::Spacing();
    // Stacked, not side by side: at this panel width the second label clips.
    int polarity = settings.keep_subject ? 1 : 0;
    if (ImGui::RadioButton("Remove what I named", polarity == 0)) polarity = 0;
    if (ImGui::RadioButton("Keep only what I named", polarity == 1)) polarity = 1;
    if ((polarity == 1) != settings.keep_subject) {
        settings.keep_subject = polarity == 1;
        edited = true;
    }
    help_tooltip_on_hover(
        "\"Remove\" is for distractors -- people, cars, the photographer's "
        "shadow. \"Keep only\" is for object captures, where everything but "
        "the subject should be ignored.");

    ImGui::Spacing();
    ImGui::Separator();
    if (!_clicks.empty()) {
        ImGui::Text("%d click%s on the image", (int)_clicks.size(),
                    _clicks.size() == 1 ? "" : "s");
        ImGui::SameLine();
        if (ImGui::SmallButton("Clear")) {
            _clicks.clear();
            edited = true;
        }
    } else {
        ImGui::TextDisabled("Click the image to point at an object.");
    }

    // ---- frame chooser ----
    if (_frames.size() > 1) {
        ImGui::Spacing();
        ImGui::SetNextItemWidth(-1);
        int idx = _frame_idx;
        if (ImGui::SliderInt("##frame", &idx, 0, (int)_frames.size() - 1,
                             "frame %d")) {
            _frame_idx = idx;
            _frame_dirty = true;
            _needs_run = true;
        }
        help_tooltip_on_hover("A few frames from across the capture. Check a "
                              "prompt on more than one before running.");
    }

    ImGui::Spacing();
    ImGui::BeginDisabled(_busy.load());
    if (ImGui::Button("Try it", ImVec2(-1, 30))) _needs_run = true;
    ImGui::EndDisabled();

    // Rerun once the user stops typing, so every keystroke does not queue a
    // forward pass.
    if (edited) _needs_run = true;
    if (_needs_run && !_busy.load() && !ImGui::IsAnyItemActive()) {
        _needs_run = false;
        start_job(settings);
    }

    ImGui::Spacing();
    {
        std::lock_guard<std::mutex> lk(_mu);
        if (!_error.empty()) {
            ImGui::PushTextWrapPos();
            ImGui::TextColored(ImVec4(1.0f, 0.45f, 0.45f, 1.0f), "%s",
                               _error.c_str());
            ImGui::PopTextWrapPos();
        } else if (!_status.empty()) {
            ImGui::TextDisabled("%s", _status.c_str());
        }
    }
    if (_kept_fraction >= 0.0f) {
        ImGui::Text("%.0f%% of the frame is kept", 100.0f * _kept_fraction);
        if (_kept_fraction < 0.05f)
            ImGui::TextColored(ImVec4(0.95f, 0.75f, 0.30f, 1.0f),
                               "Almost everything is masked out -- did you mean "
                               "\"Keep only what I named\"?");
    }

    ImGui::EndChild();

    ImGui::SameLine();
    ImGui::BeginChild("##segimg", ImVec2(0, 0));
    draw_image(settings);
    ImGui::EndChild();

    ImGui::End();
    if (!open) close();
}

}  // namespace gui
