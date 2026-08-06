// SegmentPanel.cpp -- see SegmentPanel.h.

#include "app/gui/SegmentPanel.h"

#include "app/gui/MaskPrompt.h"
#include "app/gui/Subprocess.h"
#include "app/gui/Ui.h"

#include "i18n/catalog/Dataset.h"

#include "imgui.h"
#include "imgui_stdlib.h"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cstdio>
#include <filesystem>

#ifdef SS_BUILD_SAM
#include "nn/io/Image.h"
#include "sam/Masking.h"
#endif
#ifdef SS_HAVE_VIDEO
#include "video/Video.h"
#endif

namespace fs = std::filesystem;

namespace dmsg = spirula::i18n::msg::dataset;

namespace gui {

namespace {

bool is_image_file(const fs::path& p) {
    std::string e = p.extension().string();
    for (auto& c : e) c = (char)std::tolower((unsigned char)c);
    return e == ".jpg" || e == ".jpeg" || e == ".png" || e == ".webp" ||
           e == ".tif" || e == ".tiff" || e == ".bmp";
}

// One colour per object, so a dot on the image and a row in the list are
// obviously the same thing. Red is reserved for negative clicks.
ImU32 object_color(int object) {
    static const ImU32 kColors[] = {
        IM_COL32(80, 220, 110, 255),  IM_COL32(90, 170, 245, 255),
        IM_COL32(245, 200, 70, 255),  IM_COL32(200, 130, 245, 255),
        IM_COL32(80, 225, 220, 255),  IM_COL32(245, 150, 90, 255),
    };
    const int n = (int)(sizeof(kColors) / sizeof(kColors[0]));
    return kColors[((object % n) + n) % n];
}

}  // namespace

// The GPU-side state, kept alive across runs so a prompt edit costs one
// forward pass rather than a 3-second weight upload.
struct SegmentPanel::Job {
#ifdef SS_BUILD_SAM
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
                        const std::string& model_path,
                        const std::string& ffmpeg_exe, bool force_ffmpeg) {
    _input = input;
    _is_video = is_video;
    _model_path = model_path;
    _ffmpeg_exe = ffmpeg_exe.empty() ? "ffmpeg" : ffmpeg_exe;
    _force_ffmpeg = force_ffmpeg;
    _frame_idx = 0;
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
    // Drop the model: the reconstruction that usually follows wants the VRAM,
    // and ~Masker -> Session::unload() hands the weights back for real rather
    // than leaving them in the inference layer's process-wide pool. The device
    // itself stays up -- reopening the panel is common, and rebuilding it
    // would cost a pipeline rebuild for the ~50 MB it holds. A dataset run
    // takes the device down at the end (DatasetPrep::run).
    _job.reset();
    _open = false;
}

void SegmentPanel::destroy_gl() {
    if (_tex) {
        glDeleteTextures(1, &_tex);
        _tex = 0;
    }
}

bool SegmentPanel::builtin_decode() const {
    return !_force_ffmpeg && backends().builtin_video;
}

void SegmentPanel::collect_frames(const std::string& input, bool is_video) {
    _frames.clear();
    _video_seconds = 0.0;
    std::error_code ec;
    if (is_video) {
        // Frames are decoded on demand; the list is just how many offers the
        // slider makes. Seeking is not supported by the decoder, so "frame N"
        // means "the Nth frame we sample while reading forward", and the
        // decoded index each offer lands on is filled in below.
        constexpr int kOffers = 8;
        long long total = 0;
#ifdef SS_HAVE_VIDEO
        if (builtin_decode()) {
            video::VideoReader r;
            if (r.open(input)) {
                total = r.info().frame_count;
                // Free here, and the only thing a later fallback would be
                // missing: ffmpeg seeks by time, not by frame.
                if (r.info().fps > 0.0)
                    _video_seconds = (double)total / r.info().fps;
            }
        }
#endif
        if (total <= 0) {
            // No in-process decoding on this machine, or it could not open the
            // file: ffmpeg is what will read this video, in the preview and in
            // the run, so ffmpeg is who to ask how long it is. The duration is
            // what the frame worker seeks by -- it has no reader to count
            // frames with.
            VideoFacts facts;
            if (ffmpeg_probe_video(_ffmpeg_exe, input, facts, _cancel)) {
                _video_seconds = facts.duration;
                total = facts.frames;
            }
        }
        total = std::max(total, (long long)kOffers);
        for (int i = 0; i < kOffers; i++) {
            Frame f;
            f.index = std::min(total - 1, (long long)i * (total / kOffers));
            f.position = (float)((double)f.index / (double)std::max(1LL, total - 1));
            _frames.push_back(f);
        }
        return;
    }
    std::vector<std::string> all;
    // follow_directory_symlink for the same reason DatasetPrep walks that way:
    // a prepared capture's images/ is often a link into the raw one, and the
    // default iterator quietly returns nothing for it.
    for (fs::recursive_directory_iterator it(
             input, fs::directory_options::skip_permission_denied |
                        fs::directory_options::follow_directory_symlink, ec), end;
         !ec && it != end; it.increment(ec))
        if (it->is_regular_file(ec) && is_image_file(it->path()))
            all.push_back(it->path().string());
    std::sort(all.begin(), all.end());
    // A dozen spread through the capture is enough to judge a prompt, and
    // keeps the slider meaningful on a 3000-photo folder. The index kept is
    // the one in the FULL list, because that is what the masking run counts.
    const size_t n = all.size();
    const size_t offers = std::min<size_t>(n, 12);
    for (size_t i = 0; i < offers; i++) {
        Frame f;
        f.index = offers > 1 ? (long long)(i * (n - 1) / (offers - 1)) : 0;
        f.path = all[(size_t)f.index];
        f.position = n > 1 ? (float)((double)f.index / (double)(n - 1)) : 0.0f;
        _frames.push_back(f);
    }
}

// ---------------------------------------------------------------------------
// The worker
// ---------------------------------------------------------------------------

void SegmentPanel::start_job(const MaskSettings& s) {
#ifndef SS_BUILD_SAM
    (void)s;
    std::lock_guard<std::mutex> lk(_mu);
    _error = "this build has no built-in segmentation "
             "(-DSS_BUILD_SAM=OFF)";
#else
    if (_busy.load()) return;
    if (_worker.joinable()) _worker.join();
    if (_model_path.empty()) {
        std::lock_guard<std::mutex> lk(_mu);
        _error = "no model selected -- download one first";
        return;
    }

    const int idx = _frame_idx;
    const Frame frame = (idx >= 0 && idx < (int)_frames.size()) ? _frames[idx]
                                                                : Frame{};
    const std::string frame_path = frame.path;
    const std::string input = _input;
    const bool is_video = _is_video;
    const std::string model = _model_path;
    const std::string ffmpeg_exe = _ffmpeg_exe;
    const bool builtin = builtin_decode();
    const double video_seconds = _video_seconds;
    // Only what was drawn on THIS frame: the preview segments one still, with
    // no memory bank, so a click made on another frame has nothing to say
    // about this one. The run is where they all come together.
    std::vector<MaskClick> clicks;
    for (const MaskClick& c : s.clicks)
        if (c.frame == frame.index) clicks.push_back(c);
    const MaskSettings settings = s;
    const bool frame_dirty = _frame_dirty;
    _frame_dirty = false;

    _busy = true;
    {
        std::lock_guard<std::mutex> lk(_mu);
        _error.clear();
        _status = "working...";
    }

    _worker = std::thread([this, settings, model, frame, frame_path, input, is_video,
                           idx, clicks, frame_dirty, ffmpeg_exe, builtin,
                           video_seconds] {
        // Every failure below leaves through `return set_error(...)`, so the
        // flag cannot be cleared at the end of the function: one early exit
        // would strand it at true, and start_job() refuses to run while it is
        // -- a panel that never shows anything again, whatever you type.
        struct BusyGuard {
            std::atomic<bool>& flag;
            ~BusyGuard() { flag = false; }
        } busy_guard{_busy};

        auto set_error = [&](const std::string& e) {
            std::lock_guard<std::mutex> lk(_mu);
            _error = e;
            _status.clear();
        };
        // The frame on its own, before any mask exists. Called as soon as it
        // is decoded so the panel shows a picture while the checkpoint uploads
        // -- and so that a user with a SAM 2 checkpoint, whose only prompt is
        // a click, has something to click on.
        auto show_frame = [&](const nn::Image& img) {
            std::lock_guard<std::mutex> lk(_mu);
            _preview = img.data;
            _preview_w = img.width;
            _preview_h = img.height;
            _preview_dirty = true;
            _kept_fraction = -1.0f;
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
                    nn::Image img;
#ifdef SS_HAVE_VIDEO
                    if (builtin) {
                        // No seek: read forward to the frame this offer names.
                        video::VideoReader r;
                        if (r.open(input)) {
                            for (long long i = 0; i <= frame.index; i++) {
                                if (_cancel.load()) return;
                                img = r.readFrame();
                                if (img.empty()) break;
                            }
                        }
                    }
#else
                    (void)builtin;
#endif
                    // The same fallback the dataset run takes when the driver
                    // cannot decode (DatasetPrep::extract_video): shell out to
                    // ffmpeg. It is a whole subprocess for one still, which is
                    // why it is not the first choice -- but it is the only
                    // thing between "no hardware decoding" and a panel that
                    // shows nothing at all, and the run will read this video
                    // the same way.
                    if (img.empty()) {
                        if (!command_exists(ffmpeg_exe))
                            return set_error(spirula::i18n::format(
                                dmsg::preview_needs_ffmpeg, {ffmpeg_exe}));
                        // ffmpeg seeks by time, and the very end of a file is
                        // past the last frame often enough to be worth backing
                        // away from.
                        const double at = std::min(
                            (double)frame.position * video_seconds,
                            std::max(0.0, video_seconds - 0.1));
                        const fs::path tmp =
                            fs::temp_directory_path() /
                            ("spirula-mask-preview-" +
                             std::to_string(
                                 std::chrono::steady_clock::now()
                                     .time_since_epoch()
                                     .count()) +
                             ".jpg");
                        const bool ok = ffmpeg_extract_frame(
                            ffmpeg_exe, input, at, tmp.string(), _cancel);
                        if (_cancel.load()) return;
                        if (ok) img = nn::load_image(tmp.string());
                        std::error_code rm;
                        fs::remove(tmp, rm);
                    }
                    if (img.empty())
                        return set_error(dmsg::preview_frame_unreadable.get());
                    j.frame = std::move(img);
                } else {
                    if (frame_path.empty()) return set_error("no frame to show");
                    j.frame = nn::load_image(frame_path);
                    if (j.frame.empty())
                        return set_error("could not read " + frame_path);
                }
                j.frame_key = key;
                show_frame(j.frame);
            }
            if (_cancel.load()) return;

            // Nothing to segment yet. The frame is up, which is the whole
            // point of getting here: with a SAM 2 checkpoint a click is the
            // only prompt there is, and it cannot be made on a blank panel.
            if (settings.prompt.empty() && clicks.empty()) {
                show_frame(j.frame);
                std::lock_guard<std::mutex> lk(_mu);
                _status = "Say what to look for above, or click the object in "
                          "the picture.";
                _error.clear();
                return;
            }

            // ---- the model ----
            // Every field below changes what the masker computes, so the
            // signature is what decides whether it can be reused. The weights
            // are the expensive part and only the model path moves them.
            std::string sig = settings.prompt + "|" + settings.negative_prompt + "|" +
                              std::to_string((int)settings.keep_subject) + "|" +
                              std::to_string(settings.max_image_size) + "|" +
                              std::to_string(settings.threshold);
            for (const MaskClick& c : clicks)
                sig += "|" + std::to_string(c.object) + ":" + std::to_string(c.x) +
                       "," + std::to_string(c.y) + (c.positive ? "+" : "-");
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
                for (const MaskClick& c : clicks) {
                    sam::SeedPrompt seed;
                    seed.object = c.object;
                    seed.frame = frame.index;
                    if (c.positive) seed.prompt.pos_points.push_back({c.x, c.y});
                    else            seed.prompt.neg_points.push_back({c.x, c.y});
                    mo.seeds.push_back(seed);
                }
                std::string err;
                // Re-initializing is cheap when only the prompt moved: the
                // session keeps the weights it already uploaded.
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
            if (!j.masker.run(j.frame, mask, &detections, frame.index))
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
                else if (detections.detections.empty() && !clicks.empty())
                    _error = "the clicks on this frame did not select anything";
            }
        } catch (const std::exception& e) {
            set_error(e.what());
        }
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

    const Frame frame = _frames.empty() ? Frame{} : _frames[(size_t)_frame_idx];

    // Clicks land in source-image pixels, which is what the model wants.
    const bool hovered = ImGui::IsItemHovered();
    if (hovered && (ImGui::IsMouseClicked(ImGuiMouseButton_Left) ||
                    ImGui::IsMouseClicked(ImGuiMouseButton_Right))) {
        const ImVec2 m = ImGui::GetMousePos();
        MaskClick c;
        c.x = (m.x - origin.x) / size.x * (float)_tex_w;
        c.y = (m.y - origin.y) / size.y * (float)_tex_h;
        c.positive = ImGui::IsMouseClicked(ImGuiMouseButton_Left);
        c.object = settings.current_object;
        c.frame = frame.index;
        c.position = frame.position;
        settings.clicks.push_back(c);
        start_job(settings);
    }
    if (hovered)
        ui::SetTooltip(dmsg::click_tooltip, {settings.current_object + 1});

    ImDrawList* dl = ImGui::GetWindowDrawList();
    for (const MaskClick& c : settings.clicks) {
        if (c.frame != frame.index) continue;
        const ImVec2 p(origin.x + c.x / (float)_tex_w * size.x,
                       origin.y + c.y / (float)_tex_h * size.y);
        const ImU32 col = c.positive ? object_color(c.object)
                                     : IM_COL32(240, 90, 90, 255);
        dl->AddCircleFilled(p, 6.0f, col);
        dl->AddCircle(p, 6.0f, IM_COL32(20, 20, 20, 200), 0, 1.5f);
        // A negative click is a cross, so the two are told apart without
        // relying on colour alone.
        if (!c.positive) {
            dl->AddLine(ImVec2(p.x - 3, p.y - 3), ImVec2(p.x + 3, p.y + 3),
                        IM_COL32(255, 255, 255, 255), 1.5f);
            dl->AddLine(ImVec2(p.x - 3, p.y + 3), ImVec2(p.x + 3, p.y - 3),
                        IM_COL32(255, 255, 255, 255), 1.5f);
        }
    }
}

// ---------------------------------------------------------------------------
// The object list
// ---------------------------------------------------------------------------

void SegmentPanel::draw_objects(MaskSettings& settings, bool& edited) {
    ui::Text(dmsg::objects_to_click);
    ui::help_on_hover(dmsg::objects_to_click_help);

    for (int o = 0; o < settings.object_count; ++o) {
        ImGui::PushID(o);
        int here = 0, elsewhere = 0;
        const long long cur = _frames.empty() ? 0 : _frames[(size_t)_frame_idx].index;
        for (const MaskClick& c : settings.clicks)
            if (c.object == o) (c.frame == cur ? here : elsewhere)++;

        const ImU32 col = object_color(o);
        ImGui::ColorButton("##col", ImGui::ColorConvertU32ToFloat4(col),
                           ImGuiColorEditFlags_NoTooltip |
                               ImGuiColorEditFlags_NoDragDrop,
                           ImVec2(12, 12));
        ImGui::SameLine();
        const std::string label =
            (here || elsewhere)
                ? spirula::i18n::format(dmsg::object_with_clicks,
                                        {o + 1, here, elsewhere})
                : spirula::i18n::format(dmsg::object_no_clicks, {o + 1});
        // PushID(o) above already separates the rows, so the label carries no
        // ID of its own.
        if (ui::RadioButtonRaw(label.c_str(), settings.current_object == o))
            settings.current_object = o;
        if (here || elsewhere) {
            ImGui::SameLine();
            if (ui::SmallButton(dmsg::object_clear)) {
                auto& v = settings.clicks;
                v.erase(std::remove_if(v.begin(), v.end(),
                                       [o](const MaskClick& c) { return c.object == o; }),
                        v.end());
                edited = true;
            }
        }
        ImGui::PopID();
    }

    if (ui::SmallButton(dmsg::object_another)) {
        settings.current_object = settings.object_count++;
    }
    ui::help_on_hover(dmsg::object_another_help);
    if (settings.object_count > 1) {
        ImGui::SameLine();
        if (ui::SmallButton(dmsg::object_clear_all)) {
            settings.clicks.clear();
            settings.object_count = 1;
            settings.current_object = 0;
            edited = true;
        }
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
    if (!ImGui::Begin(ui::detail::label(dmsg::preview_title), &open,
                      ImGuiWindowFlags_NoCollapse)) {
        ImGui::End();
        if (!open) close();
        return;
    }

    ui::TextDisabled(dmsg::preview_legend);
    ImGui::Separator();

    // ---- controls ----
    const float panel_w = 360.0f;
    ImGui::BeginChild("##segctl", ImVec2(panel_w, 0), ImGuiChildFlags_Borders);

    bool edited = false;

    // Polarity first, and the two fields below it are labelled by it: the same
    // box means "take this out" or "this is the subject" depending on the
    // radio, and a label that does not follow the switch reads as a bug.
    // Stacked, not side by side: at this panel width the second label clips.
    int polarity = settings.keep_subject ? 1 : 0;
    if (ui::RadioButton(dmsg::mask_remove_named, polarity == 0)) polarity = 0;
    if (ui::RadioButton(dmsg::mask_keep_named, polarity == 1)) polarity = 1;
    if ((polarity == 1) != settings.keep_subject) {
        settings.keep_subject = polarity == 1;
        edited = true;
    }
    ui::help_on_hover(dmsg::preview_polarity_help);
    const bool keep = settings.keep_subject;

    // English, whatever the interface language is -- see MaskPrompt.h.
    ImGui::Spacing();
    ui::Text(keep ? dmsg::preview_what_kept : dmsg::preview_what_removed);
    ImGui::SetNextItemWidth(-1);
    edited |= ui::InputTextEnglishRaw(
        "##prompt",
        keep ? "the statue; its pedestal" : "people; cars; my shadow",
        &settings.prompt);
    ui::help_on_hover(keep ? dmsg::preview_prompt_help_keep
                           : dmsg::preview_prompt_help_remove);

    ImGui::Spacing();
    ui::Text(keep ? dmsg::preview_but_remove_these
                  : dmsg::preview_but_keep_these);
    ImGui::SetNextItemWidth(-1);
    edited |= ui::InputTextEnglishRaw(
        "##negprompt",
        keep ? "the hand holding it" : "person in a painting",
        &settings.negative_prompt);
    ui::help_on_hover(keep ? dmsg::preview_negative_help_keep
                           : dmsg::preview_negative_help_remove);

    // The palette matters more here than on the dataset screen: this panel is
    // where a prompt is actually iterated on, one try at a time.
    if (spirula::i18n::current() != spirula::i18n::Lang::en)
        ui::TextDisabled(dmsg::mask_english_only);
    edited |= draw_subject_palette(settings.prompt, settings.negative_prompt,
                                   keep);

    ImGui::Spacing();
    ImGui::Separator();
    draw_objects(settings, edited);

    // ---- frame chooser ----
    if (_frames.size() > 1) {
        ImGui::Spacing();
        ImGui::SetNextItemWidth(-1);
        int idx = _frame_idx;
        // The slider's own overlay format string: the index is substituted by
        // ImGui, so this is not a place a placeholder can be used.
        const std::string fmt = spirula::i18n::format(
            dmsg::preview_frame, {(long long)_frames[(size_t)_frame_idx].index});
        if (ui::SliderIntRaw("##frame", &idx, 0, (int)_frames.size() - 1,
                             fmt.c_str())) {
            _frame_idx = idx;
            _frame_dirty = true;
            _needs_run = true;
        }
        ui::help_on_hover(dmsg::preview_frame_help);
    }

    ImGui::Spacing();
    ImGui::BeginDisabled(_busy.load());
    if (ui::Button(dmsg::preview_try_it, ImVec2(-1, 30))) _needs_run = true;
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
        // Segmentation errors and progress come from the inference layer.
        if (!_error.empty())
            ui::TextColoredWrappedRaw(ImVec4(1.0f, 0.45f, 0.45f, 1.0f), _error);
        else if (!_status.empty())
            ui::TextDisabledRaw(_status);
    }
    if (_kept_fraction >= 0.0f) {
        char pct[16];
        std::snprintf(pct, sizeof pct, "%.0f", 100.0f * _kept_fraction);
        ui::Text(dmsg::preview_kept_fraction, {pct});
        if (_kept_fraction < 0.05f)
            ui::TextColoredWrapped(ImVec4(0.95f, 0.75f, 0.30f, 1.0f),
                                   keep ? dmsg::preview_almost_nothing_kept
                                        : dmsg::preview_almost_all_masked);
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
