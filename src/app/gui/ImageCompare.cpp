// ImageCompare.cpp -- see ImageCompare.h.

#include "app/gui/ImageCompare.h"

#include "app/TrainerCore.h"
#include "app/gui/Ui.h"
#include "engine/Engine.h"
#include "external/stb_image.h"

#include "i18n/catalog/Gui.h"

#include "imgui.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>

namespace msg = spirula::i18n::msg::gui;

namespace gui {

namespace {

uint8_t to_byte(float v) {
    return (uint8_t)(std::min(std::max(v, 0.0f), 1.0f) * 255.0f + 0.5f);
}

// A split image's views are drawn as one picture, so the zoom moves over all
// of them at once instead of paging between them.
//
// Where each one goes, and how far it has to be turned to sit flush against
// its neighbours. The views are cubemap faces in the order of
// DataManager.cpp's kAxesFisheye5 / kAxesEquirect6 -- 0 front, 1 right,
// 2 down, 3 left, 4 up, and for the 6-face table 5 back -- so unfolding them
// as a cross is what makes the seams meet: five of the twelve cube edges
// join, which is as many as any flat unfolding can. The two tables differ in
// the right face's roll, which is why one turns and the other does not.
//
// This is LAYOUT only. The texture behind it packs the views tightly (see
// pack_rgb), and the cross's empty cells cost nothing: draw_pane places each
// view from the packed texture on its own quad.
const FaceCell kCross5[5] = {
    {1, 1, 0}, {2, 1, 1}, {1, 2, 2}, {0, 1, 3}, {1, 0, 0},
};
const FaceCell kCross6[6] = {
    {1, 1, 0}, {2, 1, 0}, {1, 2, 2}, {0, 1, 3}, {1, 0, 0}, {3, 1, 2},
};

void face_layout(int views, std::vector<FaceCell>& cells, int& cols, int& rows) {
    if (views == 5) { cells.assign(kCross5, kCross5 + 5); cols = 3; rows = 3; return; }
    if (views == 6) { cells.assign(kCross6, kCross6 + 6); cols = 4; rows = 3; return; }
    views = std::max(views, 1);
    cols = (int)std::ceil(std::sqrt((double)views));
    rows = (views + cols - 1) / cols;
    cells.clear();
    for (int v = 0; v < views; v++) cells.push_back({v % cols, v / cols, 0});
}

// The packed grid: square-ish, and every cell used.
void pack_shape(int views, int& cols, int& rows) {
    views = std::max(views, 1);
    cols = (int)std::ceil(std::sqrt((double)views));
    rows = (views + cols - 1) / cols;
}

// float [B, h, w, C] -> uint8 [rows*h, cols*w, C], one view per packed cell.
void pack_rgb(const float* src, int views, int h, int w, int C,
              int cols, int rows, std::vector<uint8_t>& dst) {
    const int W = cols * w;
    dst.assign((size_t)rows * h * W * C, 0);
    for (int v = 0; v < views; v++) {
        const float* s = src + (size_t)v * h * w * C;
        const int oy = (v / cols) * h, ox = (v % cols) * w;
        for (int y = 0; y < h; y++) {
            uint8_t* d = dst.data() + ((size_t)(oy + y) * W + ox) * C;
            for (int x = 0; x < w * C; x++)
                d[x] = to_byte(s[(size_t)y * w * C + x]);
        }
    }
}

// The same, for the mask -- which the engine keeps at whatever resolution the
// files are (the loss samples it rather than resizing it), so it is generally
// NOT the render's shape. Nearest, because it is boolean.
void pack_mask(const uint8_t* src, int views, int h, int w, int mh, int mw,
               int cols, int rows, std::vector<uint8_t>& dst) {
    const int W = cols * w;
    dst.assign((size_t)rows * h * W, 0);
    for (int v = 0; v < views; v++) {
        const uint8_t* s = src + (size_t)v * mh * mw;
        const int oy = (v / cols) * h, ox = (v % cols) * w;
        for (int y = 0; y < h; y++) {
            const int my = (mh == h) ? y : (int)((int64_t)y * mh / h);
            const uint8_t* sr = s + (size_t)my * mw;
            uint8_t* d = dst.data() + (size_t)(oy + y) * W + ox;
            for (int x = 0; x < w; x++)
                d[x] = sr[(mw == w) ? x : (int)((int64_t)x * mw / w)];
        }
    }
}

std::string lower(std::string s) {
    for (char& c : s) c = (char)std::tolower((unsigned char)c);
    return s;
}

// A one-line placeholder where the pair would be, laid out like the
// viewport's (ViewportPanel::draw).
void centered_note(const spirula::i18n::Msg& m) {
    ImGui::Dummy(ImVec2(0, ImGui::GetContentRegionAvail().y * 0.4f));
    const float w = ImGui::CalcTextSize(m.get()).x;
    ImGui::SetCursorPosX(std::max(0.0f, (ImGui::GetWindowWidth() - w) * 0.5f));
    ui::TextDisabled(m);
}

}  // namespace

// ---------------------------------------------------------------------------
// Attach / detach + worker
// ---------------------------------------------------------------------------

ImageCompare::~ImageCompare() { detach(); }

void ImageCompare::attach(spirula::TrainerSession& session) {
    detach();
    _session = &session;
    _index = std::clamp(_index, 0,
                        (int)std::max<int64_t>(session.ds.num_cameras - 1, 0));
    _error.clear();
    _shot = Shot{};
    _requested = Job{-1, false, false};
    _dirty = true;
    _running = true;
    _worker = std::thread([this] { worker_loop(); });
}

void ImageCompare::detach() {
    if (_running) {
        {
            std::lock_guard<std::mutex> lk(_mu);
            _running = false;
        }
        _cv.notify_all();
        if (_worker.joinable()) _worker.join();
    }
    _session = nullptr;
    _in_flight = 0;
    _has_pending = false;
    _shot = Shot{};
    _tex_dirty = true;
}

void ImageCompare::destroy_gl() {
    for (Pane* p : {&_gt, &_render})
        if (p->tex) { glDeleteTextures(1, &p->tex); p->tex = 0; }
}

void ImageCompare::worker_loop() {
    for (;;) {
        Job j;
        uint64_t id = 0;
        {
            std::unique_lock<std::mutex> lk(_mu);
            _cv.wait(lk, [this] { return _has_pending || !_running; });
            if (!_running) return;
            j = _pending;
            id = _pending_id;
            _has_pending = false;
        }
        Shot s;
        s.job = j;
        try {
            run_job(j, s);
        } catch (const std::exception& e) {
            s.error = e.what();
        }
        s.id = id;
        bool more;
        {
            std::lock_guard<std::mutex> lk(_mu);
            _result = std::move(s);
            more = _has_pending;
        }
        // Nothing here wants the engine any more; let the trainer step.
        if (!more) _session->render_pending = false;
    }
}

void ImageCompare::run_job(const Job& j, Shot& out) {
    spirula::TrainerSession& s = *_session;
    out.step = s.cur_step.load();

    const int sh_deg = out.step / std::max(s.cfg.sh_degree_warmup_every, 1);
    // Members, not locals: a 4K pair is ~200 MB of readback, and freshly
    // mapped scratch is faulted in and zeroed by the kernel on every job.
    std::vector<float>& gt = _wgt;
    std::vector<float>& render = _wrender;
    std::vector<uint8_t>& alpha = _walpha;
    alpha.clear();
    int64_t B = 0, H = 0, W = 0, C = 0;
    int64_t mh = 0, mw = 0;   // the mask's own resolution, which need not be H, W
    {
        std::lock_guard<std::mutex> lk(s.engine_mutex);
        engine_preview_forward(j.index, s.cfg.primitive, sh_deg,
                               s.cfg.packed || s.cfg.use_bvh, j.color_correct);
        auto shape = engine_get_render_rgb_shape();
        B = std::get<0>(shape); H = std::get<1>(shape);
        W = std::get<2>(shape); C = std::get<3>(shape);
        const int64_t n = B * H * W * C;
        if (n <= 0) {
            // Nothing came back and nothing threw. Say so rather than leave
            // the panel on "rendering" for the rest of the run.
            out.error = "engine_preview_forward produced no image";
            return;
        }
        render.resize((size_t)n);
        gt.resize((size_t)n);
        engine_copy_render_to_host(
            TorchTensorView{(uint64_t)(uintptr_t)render.data(), 4, {B, H, W, C}},
            TorchTensorView{0, 0, {}}, TorchTensorView{0, 0, {}},
            TorchTensorView{0, 0, {}}, TorchTensorView{0, 0, {}});
        engine_copy_gt_rgb_to_host(
            TorchTensorView{(uint64_t)(uintptr_t)gt.data(), 4, {B, H, W, C}});
        // The mask keeps its own resolution: the loss samples it rather than
        // resizing it, so an unsplit 1920x1920 fisheye trained against a
        // 1600x1600 mask is ordinary. Only the batch dimension has to agree.
        auto ashape = engine_get_gt_alpha_shape();
        if (std::get<0>(ashape) == B) {
            mh = std::get<1>(ashape);
            mw = std::get<2>(ashape);
            alpha.resize((size_t)(B * mh * mw));
            engine_copy_gt_alpha_to_host(
                TorchTensorView{(uint64_t)(uintptr_t)alpha.data(), 1,
                                {B, mh, mw, 1LL}});
        }
    }

    // Scored over the pixels the loss is computed on, which is why the mask
    // enters here and not only on screen.
    {
        double se = 0.0;
        int64_t n = 0;
        for (int64_t b = 0; b < B; b++)
            for (int64_t y = 0; y < H; y++)
                for (int64_t x = 0; x < W; x++) {
                    if (!alpha.empty()) {
                        const int64_t my = mh == H ? y : y * mh / H;
                        const int64_t mx = mw == W ? x : x * mw / W;
                        if (!alpha[(size_t)((b * mh + my) * mw + mx)]) continue;
                    }
                    const int64_t p = (b * H + y) * W + x;
                    for (int64_t c = 0; c < C; c++) {
                        const double d = (double)gt[(size_t)(p * C + c)] -
                                         (double)render[(size_t)(p * C + c)];
                        se += d * d;
                    }
                    n += C;
                }
        if (n > 0)
            out.psnr = (float)(-10.0 * std::log10(std::max(se / (double)n, 1e-12)));
    }

    out.views = (int)B;
    out.view_w = (int)W;
    out.view_h = (int)H;
    face_layout(out.views, out.cells, out.lay_cols, out.lay_rows);
    pack_shape(out.views, out.pack_cols, out.pack_rows);
    pack_rgb(gt.data(), out.views, (int)H, (int)W, (int)C,
             out.pack_cols, out.pack_rows, out.gt);
    pack_rgb(render.data(), out.views, (int)H, (int)W, (int)C,
             out.pack_cols, out.pack_rows, out.render);
    if (!alpha.empty())
        pack_mask(alpha.data(), out.views, (int)H, (int)W, (int)mh, (int)mw,
                  out.pack_cols, out.pack_rows, out.mask);

    if (j.source_gt && j.index < (int)s.ds.image_filenames.size()) {
        int w = 0, h = 0, ch = 0;
        stbi_uc* img = stbi_load(s.ds.image_filenames[(size_t)j.index].c_str(),
                                 &w, &h, &ch, 3);
        if (img) {
            out.src.assign(img, img + (size_t)w * h * 3);
            out.src_w = w;
            out.src_h = h;
            stbi_image_free(img);
        }
    }
}

void ImageCompare::submit(const Job& j) {
    {
        std::lock_guard<std::mutex> lk(_mu);
        _pending = j;
        _pending_id = _next_id++;
        _has_pending = true;
        _in_flight = _pending_id;
    }
    // Synchronously, so the training loop yields at its next boundary instead
    // of after the whole step this job would otherwise queue behind.
    _session->render_pending = true;
    _cv.notify_one();
    _submitted_at = ImGui::GetTime();
}

void ImageCompare::poll() {
    if (!_in_flight) return;
    Shot s;
    {
        std::lock_guard<std::mutex> lk(_mu);
        if (_result.id != _in_flight) return;
        s = std::move(_result);
        _result = Shot{};
    }
    _in_flight = 0;
    const double took = ImGui::GetTime() - _submitted_at;
    if (took > 0) _job_secs = _job_secs > 0 ? 0.7 * _job_secs + 0.3 * took : took;
    if (!s.error.empty()) {
        _error = s.error;
        return;
    }
    _error.clear();
    _shot = std::move(s);
    _rendered_step = _shot.step;
    _tex_dirty = true;
}

// ---------------------------------------------------------------------------
// Compositing + upload
// ---------------------------------------------------------------------------

void ImageCompare::upload(GLuint& tex, const uint8_t* rgb, int w, int h) {
    if (!tex) {
        glGenTextures(1, &tex);
        glBindTexture(GL_TEXTURE_2D, tex);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    }
    glBindTexture(GL_TEXTURE_2D, tex);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, w, h, 0, GL_RGB,
                 GL_UNSIGNED_BYTE, rgb);
}

void ImageCompare::rebuild_textures() {
    _tex_dirty = false;
    if (_shot.view_w <= 0 || _shot.view_h <= 0) return;

    const int pw = _shot.pack_cols * _shot.view_w;
    const int ph = _shot.pack_rows * _shot.view_h;

    // The excluded region is marked the same way on both sides, so a
    // difference between the panes is a difference in the model rather than
    // in what each of them chose to hide. The mask is in packed coordinates,
    // which is what the render side is in too; only a source-file photograph
    // (never split, so one packed cell) has to be scaled into it.
    auto composite = [&](const uint8_t* rgb, int w, int h,
                         std::vector<uint8_t>& out) {
        out.assign(rgb, rgb + (size_t)w * h * 3);
        if (_mask_show == MaskShow::Ignore || _shot.mask.empty()) return;
        for (int y = 0; y < h; y++) {
            const int my = (h == ph) ? y : (int)((int64_t)y * ph / h);
            for (int x = 0; x < w; x++) {
                const int mx = (w == pw) ? x : (int)((int64_t)x * pw / w);
                if (_shot.mask[(size_t)my * pw + mx]) continue;
                uint8_t* p = out.data() + ((size_t)y * w + x) * 3;
                if (_mask_show == MaskShow::Hide) p[0] = p[1] = p[2] = 0;
                else for (int c = 0; c < 3; c++) p[c] = (uint8_t)(p[c] / 4);
            }
        }
    };

    auto to_pane = [&](Pane& p, const uint8_t* rgb, int w, int h, bool packed,
                       std::vector<uint8_t>& buf) {
        composite(rgb, w, h, buf);
        upload(p.tex, buf.data(), w, h);
        p.tex_w = w;
        p.tex_h = h;
        if (packed) {
            p.view_w = _shot.view_w;
            p.view_h = _shot.view_h;
            p.pack_cols = _shot.pack_cols;
            p.cells = _shot.cells;
            p.canvas_w = _shot.lay_cols * _shot.view_w;
            p.canvas_h = _shot.lay_rows * _shot.view_h;
        } else {
            p.view_w = w;
            p.view_h = h;
            p.pack_cols = 1;
            p.cells.assign(1, FaceCell{0, 0, 0});
            p.canvas_w = w;
            p.canvas_h = h;
        }
    };

    std::vector<uint8_t> buf;
    // The source file is only offered for an unsplit image, so its canvas has
    // the same aspect as the render's and one zoom still serves both.
    const bool use_src = _shot.job.source_gt && !_shot.src.empty();
    if (use_src) to_pane(_gt, _shot.src.data(), _shot.src_w, _shot.src_h, false, buf);
    else         to_pane(_gt, _shot.gt.data(), pw, ph, true, buf);
    to_pane(_render, _shot.render.data(), pw, ph, true, buf);
}

// ---------------------------------------------------------------------------
// Selection
// ---------------------------------------------------------------------------

int ImageCompare::faces() const {
    const auto& K = _session->post.K_per_camera;
    if (K.empty() || _index >= (int)K.size()) return 1;
    return std::max(1, K[(size_t)_index]);
}

bool ImageCompare::has_masks() const {
    return !_session->ds.mask_filenames.empty() || _session->post.any_fisheye_warp;
}

void ImageCompare::select(int index) {
    const int n = (int)_session->ds.num_cameras;
    if (n <= 0) return;
    index = std::clamp(index, 0, n - 1);
    if (index == _index) return;
    _index = index;
    _dirty = true;
}

// ---------------------------------------------------------------------------
// Toolbar
// ---------------------------------------------------------------------------

void ImageCompare::draw_image_picker() {
    const auto& names = _session->ds.image_filenames;
    const int n = (int)names.size();
    auto base = [&](int i) {
        return std::filesystem::path(names[(size_t)i]).filename().string();
    };

    place(ImGui::GetFrameHeight() * 2 + ImGui::GetStyle().ItemSpacing.x);
    ImGui::BeginDisabled(n <= 1);
    if (ui::ButtonRaw("<")) select(_index - 1);
    ImGui::SameLine();
    if (ui::ButtonRaw(">")) select(_index + 1);
    ImGui::EndDisabled();
    ui::help_on_hover(msg::compare_navigate_help);

    // "12 / 240" is two numbers and a separator: the same in every language.
    char fmt[32];
    std::snprintf(fmt, sizeof fmt, "%%d / %d", std::max(n, 1));
    int one_based = _index + 1;
    place(190);
    ImGui::SetNextItemWidth(190);
    if (ui::SliderIntRaw("##imgidx", &one_based, 1, std::max(n, 1), fmt))
        select(one_based - 1);
    ui::help_on_hover(msg::compare_slider_help);

    place(260);
    ImGui::SetNextItemWidth(260);
    // File names are file names in every language.
    if (ui::BeginComboRaw("##imgname", n > 0 ? base(_index).c_str() : "")) {
        ImGui::SetNextItemWidth(-1);
        ui::InputTextWithHintRaw("##imgfilter", msg::compare_search_hint,
                                 &_filter);
        const std::string needle = lower(_filter);
        std::vector<int> hits;
        hits.reserve((size_t)n);
        for (int i = 0; i < n; i++)
            if (needle.empty() || lower(base(i)).find(needle) != std::string::npos)
                hits.push_back(i);
        ImGui::BeginChild("##imglist", ImVec2(0, 260));
        ImGuiListClipper clip;
        clip.Begin((int)hits.size());
        while (clip.Step())
            for (int r = clip.DisplayStart; r < clip.DisplayEnd; r++) {
                const int i = hits[(size_t)r];
                if (ui::SelectableRaw(base(i), i == _index)) {
                    select(i);
                    ImGui::CloseCurrentPopup();
                }
            }
        clip.End();
        ImGui::EndChild();
        ImGui::EndCombo();
    }
    ui::help_on_hover(msg::compare_pick_help);
}

// One flowing row of controls, wrapping onto as many rows as they need --
// which ones are present depends on the dataset and on the run, so a fixed
// layout would overflow on some combinations and waste a row on others.
// `used` is TRACKED rather than read back from the last item's rectangle,
// which would pick up a tooltip's width whenever one is hovered.
void ImageCompare::place(float w) {
    const ImGuiStyle& st = ImGui::GetStyle();
    if (_row_used > 0.0f && _row_used + st.ItemSpacing.x + w <= _row_w) {
        ImGui::SameLine();
        _row_used += st.ItemSpacing.x + w;
    } else {
        _row_used = w;
    }
}

void ImageCompare::draw_toolbar(bool training) {
    const ImGuiStyle& st = ImGui::GetStyle();
    _row_w = ImGui::GetContentRegionAvail().x;
    _row_used = 0.0f;
    auto check_w = [&](const spirula::i18n::Msg& m) {
        return ImGui::GetFrameHeight() + st.ItemInnerSpacing.x +
               ImGui::CalcTextSize(m.get()).x;
    };

    draw_image_picker();

    if (has_masks()) {
        const float w = ImGui::CalcTextSize(msg::compare_mask_dim.get()).x + 60;
        place(w);
        ImGui::SetNextItemWidth(w);
        int m = (int)_mask_show;
        if (ui::ComboRaw("##maskshow", &m,
                         {&msg::compare_mask_dim, &msg::compare_mask_ignore,
                          &msg::compare_mask_hide})) {
            _mask_show = (MaskShow)m;
            _tex_dirty = true;
        }
        ui::help_on_hover(msg::compare_mask_help);
    }

    if (_session->st.bilagrid_rgb_init || _session->st.ppisp_init) {
        place(check_w(msg::compare_color_match));
        if (ui::Checkbox(msg::compare_color_match, &_color_correct))
            _dirty = true;
        ui::help_on_hover(msg::compare_color_match_help);
    }

    // A split image has no single source frame its views line up with, so
    // there is nothing to offer there.
    const bool split = faces() > 1;
    place(check_w(msg::compare_source_gt));
    ImGui::BeginDisabled(split);
    if (ui::Checkbox(msg::compare_source_gt, &_source_gt)) _dirty = true;
    ImGui::EndDisabled();
    ui::help_on_hover(split ? msg::compare_source_gt_split
                            : msg::compare_source_gt_help);

    place(check_w(msg::compare_smooth));
    ui::Checkbox(msg::compare_smooth, &_smooth);
    ui::help_on_hover(msg::compare_smooth_help);

    if (training) {
        place(check_w(msg::viewport_live));
        ui::Checkbox(msg::viewport_live, &_live);
        ui::help_on_hover(msg::compare_live_help);
    }

    // A percentage is a percentage in every language.
    char z[16];
    std::snprintf(z, sizeof z, "%d%%", (int)std::lround(_zoom * 100.0f));
    place(ImGui::CalcTextSize(msg::compare_reset_zoom.get()).x +
          2.0f * st.FramePadding.x + st.ItemSpacing.x +
          ImGui::CalcTextSize(z).x);
    if (ui::Button(msg::compare_reset_zoom)) {
        _zoom = 1.0f;
        _cx = _cy = 0.5f;
    }
    ImGui::SameLine();
    ui::TextDisabledRaw(z);
    ui::help_on_hover(msg::compare_zoom_help);
}

// ---------------------------------------------------------------------------
// Panes
// ---------------------------------------------------------------------------

void ImageCompare::handle_view_input() {
    ImGuiIO& io = ImGui::GetIO();
    const float dw = std::max(_hot_w, 1.0f), dh = std::max(_hot_h, 1.0f);

    if (_hot && io.MouseWheel != 0.0f) {
        // Keep whatever the cursor is over where it is.
        const float u = _cx + (io.MousePos.x - _hot_min.x - _hot_size.x * 0.5f) / dw;
        const float v = _cy + (io.MousePos.y - _hot_min.y - _hot_size.y * 0.5f) / dh;
        const float prev = _zoom;
        _zoom = std::clamp(_zoom * std::pow(1.2f, io.MouseWheel), 1.0f, 64.0f);
        const float k = prev / _zoom;
        _cx = u + (_cx - u) * k;
        _cy = v + (_cy - v) * k;
    }

    if (_hot && ImGui::IsMouseClicked(ImGuiMouseButton_Left)) _dragging = true;
    if (_dragging) {
        if (!ImGui::IsMouseDown(ImGuiMouseButton_Left)) _dragging = false;
        else {
            _cx -= io.MouseDelta.x / dw;
            _cy -= io.MouseDelta.y / dh;
        }
    }
}

void ImageCompare::draw_pane(const Pane& p, const ImVec2& box,
                             const spirula::i18n::Msg& caption) {
    // No scrolling: the wheel is the zoom here, and a child that can scroll
    // eats it before the pane sees it.
    ImGui::BeginChild(caption.id, box, ImGuiChildFlags_Borders,
                      ImGuiWindowFlags_NoScrollbar |
                      ImGuiWindowFlags_NoScrollWithMouse);
    ui::TextDisabled(caption);
    ImVec2 avail = ImGui::GetContentRegionAvail();
    avail.x = std::max(avail.x, 16.0f);
    avail.y = std::max(avail.y, 16.0f);

    if (p.tex && p.canvas_w > 0 && p.canvas_h > 0) {
        // The canvas is aspect-fit at zoom 1 and cropped beyond it, so the
        // drawn area never overflows the pane: both panes stay the same size
        // whatever each one's resolution is, and one zoom drives both.
        const float fit = std::min(avail.x / (float)p.canvas_w,
                                   avail.y / (float)p.canvas_h);
        const float disp_w = (float)p.canvas_w * fit * _zoom;
        const float disp_h = (float)p.canvas_h * fit * _zoom;
        const float vw = std::min(avail.x, disp_w), vh = std::min(avail.y, disp_h);
        const float fu = vw / disp_w, fv = vh / disp_h;
        _cx = std::clamp(_cx, fu * 0.5f, 1.0f - fu * 0.5f);
        _cy = std::clamp(_cy, fv * 0.5f, 1.0f - fv * 0.5f);

        ImGui::SetCursorPos(
            ImVec2(ImGui::GetCursorPosX() + (avail.x - vw) * 0.5f,
                   ImGui::GetCursorPosY() + (avail.y - vh) * 0.5f));
        const ImVec2 org = ImGui::GetCursorScreenPos();
        ImGui::Dummy(ImVec2(vw, vh));

        // Canvas origin in screen space, and one cell's size there. Each view
        // is drawn straight from its slot in the packed texture onto its own
        // quad, so the cells the cross leaves empty simply are not drawn.
        const float ox = org.x - (_cx - fu * 0.5f) * disp_w;
        const float oy = org.y - (_cy - fv * 0.5f) * disp_h;
        const float cw = (float)p.view_w * fit * _zoom;
        const float ch = (float)p.view_h * fit * _zoom;
        // Half a texel in from each edge: adjacent cells come from slots that
        // are NOT adjacent in the texture, and linear filtering at the seam
        // would otherwise blend in whichever face happens to be packed next.
        // A lone view has no neighbour to bleed from and is left exact, so a
        // source-file photograph still lines up with the render texel for texel.
        const bool inset = p.cells.size() > 1;
        const float iu = inset ? 0.5f / (float)p.tex_w : 0.0f;
        const float iv = inset ? 0.5f / (float)p.tex_h : 0.0f;

        const ImGuiPlatformIO& pio = ImGui::GetPlatformIO();
        ImDrawList* dl = ImGui::GetWindowDrawList();
        // The backend binds its own LINEAR sampler object, which overrides
        // anything glTexParameteri put on the texture -- point sampling has to
        // be asked for through these two draw callbacks, in a matched pair.
        const bool point = !_smooth && pio.DrawCallback_SetSamplerNearest &&
                           pio.DrawCallback_SetSamplerLinear;
        dl->PushClipRect(org, ImVec2(org.x + vw, org.y + vh), true);
        if (point) dl->AddCallback(pio.DrawCallback_SetSamplerNearest, nullptr);
        for (size_t v = 0; v < p.cells.size(); v++) {
            const FaceCell& fc = p.cells[v];
            const float x0 = ox + fc.col * cw, x1 = x0 + cw;
            const float y0 = oy + fc.row * ch, y1 = y0 + ch;
            if (x1 < org.x || x0 > org.x + vw || y1 < org.y || y0 > org.y + vh)
                continue;
            const int pc = (int)v % p.pack_cols, pr = (int)v / p.pack_cols;
            const float u0 = (float)(pc * p.view_w) / p.tex_w + iu;
            const float u1 = (float)((pc + 1) * p.view_w) / p.tex_w - iu;
            const float v0 = (float)(pr * p.view_h) / p.tex_h + iv;
            const float v1 = (float)((pr + 1) * p.view_h) / p.tex_h - iv;
            const ImVec2 quad[4] = {{x0, y0}, {x1, y0}, {x1, y1}, {x0, y1}};
            const ImVec2 uv[4]   = {{u0, v0}, {u1, v0}, {u1, v1}, {u0, v1}};
            // A quarter turn clockwise is one step around the corner list.
            const int t = fc.turns & 3;
            dl->AddImageQuad((ImTextureID)(intptr_t)p.tex,
                             quad[0], quad[1], quad[2], quad[3],
                             uv[(4 - t) & 3], uv[(5 - t) & 3],
                             uv[(6 - t) & 3], uv[(7 - t) & 3]);
        }
        if (point) dl->AddCallback(pio.DrawCallback_SetSamplerLinear, nullptr);
        dl->PopClipRect();

        // The wheel and the drag are applied ONCE per frame, after both panes
        // (see draw_panes): handling them per pane would double every delta on
        // a drag that leaves one pane hovered.
        if (ImGui::IsItemHovered()) {
            _hot = true;
            _hot_min = org;
            _hot_size = ImVec2(vw, vh);
            _hot_w = disp_w;
            _hot_h = disp_h;
        }
    }
    ImGui::EndChild();
}

void ImageCompare::draw_panes(const ImVec2& avail) {
    const float half = std::max(64.0f,
                                (avail.x - ImGui::GetStyle().ItemSpacing.x) * 0.5f);
    const ImVec2 box(half, avail.y);
    _hot = false;
    draw_pane(_gt, box, msg::compare_pane_gt);
    ImGui::SameLine();
    draw_pane(_render, box, msg::compare_pane_render);
    handle_view_input();
}

void ImageCompare::handle_keys() {
    if (ImGui::GetIO().WantTextInput) return;
    if (!ImGui::IsWindowHovered(ImGuiHoveredFlags_ChildWindows) &&
        !ImGui::IsWindowFocused(ImGuiFocusedFlags_ChildWindows))
        return;
    if (ImGui::IsKeyPressed(ImGuiKey_LeftArrow, true))  select(_index - 1);
    if (ImGui::IsKeyPressed(ImGuiKey_RightArrow, true)) select(_index + 1);
    if (ImGui::IsKeyPressed(ImGuiKey_PageUp, true))     select(_index - 10);
    if (ImGui::IsKeyPressed(ImGuiKey_PageDown, true))   select(_index + 10);
    if (ImGui::IsKeyPressed(ImGuiKey_Home, false))      select(0);
    if (ImGui::IsKeyPressed(ImGuiKey_End, false))
        select((int)_session->ds.num_cameras - 1);
}

// ---------------------------------------------------------------------------
// Draw
// ---------------------------------------------------------------------------

void ImageCompare::draw(bool training, int step) {
    if (!_session || _session->ds.num_cameras == 0) {
        centered_note(msg::compare_waiting);
        return;
    }

    const double now = ImGui::GetTime();
    if (step >= 0 && step != _last_seen_step) {
        if (_last_seen_step >= 0 && step > _last_seen_step && _last_step_time > 0) {
            const double per = (now - _last_step_time) / (step - _last_seen_step);
            _step_secs = _step_secs > 0 ? 0.8 * _step_secs + 0.2 * per : per;
        }
        _last_seen_step = step;
        _last_step_time = now;
    }

    poll();
    draw_toolbar(training);
    handle_keys();

    // Re-request when the selection changed, and -- while training -- often
    // enough to follow the optimization without taking much of it. Paced in
    // ITERATIONS from the measured cost of a job against the cost of a step,
    // so the share it takes does not depend on how fast the dataset trains.
    const Job want{_index, _color_correct, _source_gt};
    bool follow = false;
    if (training && _live && _rendered_step >= 0 && step >= 0) {
        const int interval = (_job_secs > 0 && _step_secs > 0)
            ? (int)std::clamp(std::ceil(_job_secs / (0.08 * _step_secs)),
                              25.0, 2000.0)
            : 50;
        follow = step - _rendered_step >= interval;
    }
    if (!_in_flight && (_dirty || follow || !(want == _requested))) {
        _dirty = false;
        _requested = want;
        submit(want);
    }
    if (_tex_dirty) rebuild_textures();

    if (!_error.empty()) {
        ui::TextColored(ImVec4(1, 0.4f, 0.4f, 1), msg::viewport_render_error,
                        {_error});
        return;
    }

    if (_shot.view_w == 0) {
        centered_note(msg::viewport_rendering);
        return;
    }

    // One line of numbers: the resolution the loss sees, what the file itself
    // holds when the two differ, how many views a split image became, the
    // score, and the step the render came from.
    ui::TextDisabled(msg::compare_trained_size, {_shot.view_w, _shot.view_h});
    ui::help_on_hover(msg::compare_trained_size_help);
    if (_shot.job.source_gt && _shot.src_w > 0) {
        ImGui::SameLine();
        ui::TextDisabled(msg::compare_source_size, {_shot.src_w, _shot.src_h});
    }
    if (_shot.views > 1) {
        ImGui::SameLine();
        ui::TextDisabled(msg::compare_faces, {_shot.views});
        ui::help_on_hover(msg::compare_faces_help);
    }
    ImGui::SameLine();
    if (_shot.psnr >= 0) {
        char v[32];
        std::snprintf(v, sizeof v, "%.2f", (double)_shot.psnr);
        ui::TextDisabled(msg::compare_psnr, {v});
    } else {
        ui::TextDisabled(msg::compare_psnr_none);
    }
    ui::help_on_hover(msg::compare_psnr_help);
    if (_shot.step >= 0) {
        ImGui::SameLine();
        ui::TextDisabled(msg::compare_at_step, {(long long)_shot.step});
    }

    ImVec2 avail = ImGui::GetContentRegionAvail();
    avail.x = std::max(avail.x, 64.0f);
    avail.y = std::max(avail.y, 64.0f);
    draw_panes(avail);
}

}  // namespace gui
