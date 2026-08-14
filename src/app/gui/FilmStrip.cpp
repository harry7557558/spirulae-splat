// FilmStrip.cpp -- see FilmStrip.h.

#include "app/gui/FilmStrip.h"

#include "app/gui/Layout.h"
#include "app/gui/Ui.h"

#include <algorithm>

namespace gui {

FilmStrip::~FilmStrip() = default;

bool FilmStrip::wants(double min_interval_s) const {
    const auto now = std::chrono::steady_clock::now();
    std::lock_guard<std::mutex> lk(_mu);
    if (_seq && std::chrono::duration<double>(now - _last).count() < min_interval_s)
        return false;
    _last = now;
    return true;
}

void FilmStrip::offer(const uint8_t* rgb, int w, int h, const std::string& name,
                      const uint8_t* mask, FramePoints points) {
    if (!rgb || w <= 0 || h <= 0) return;
    const int step = std::max(1, (std::max(w, h) + kFilmMaxSide - 1) / kFilmMaxSide);
    const int dw = std::max(1, w / step), dh = std::max(1, h / step);

    // Box average over the step x step source block: a point sample of a
    // 4K frame decimated 16x aliases into noise, which reads as a bad mask.
    std::vector<uint8_t> out((size_t)dw * dh * 3);
    for (int y = 0; y < dh; y++) {
        for (int x = 0; x < dw; x++) {
            int acc[3] = {0, 0, 0};
            int n = 0, keep = 0;
            for (int sy = y * step; sy < std::min(h, (y + 1) * step); sy++) {
                for (int sx = x * step; sx < std::min(w, (x + 1) * step); sx++) {
                    const size_t si = (size_t)sy * w + sx;
                    acc[0] += rgb[si * 3 + 0];
                    acc[1] += rgb[si * 3 + 1];
                    acc[2] += rgb[si * 3 + 2];
                    keep += !mask || mask[si] > 127;
                    n++;
                }
            }
            if (!n) continue;
            const size_t di = ((size_t)y * dw + x) * 3;
            // Same tint the mask preview uses, so the two read as one answer.
            if (mask && keep * 2 < n) {
                out[di + 0] = (uint8_t)(acc[0] / n / 3 + 150);
                out[di + 1] = (uint8_t)(acc[1] / n / 3);
                out[di + 2] = (uint8_t)(acc[2] / n / 3);
            } else {
                out[di + 0] = (uint8_t)(acc[0] / n);
                out[di + 1] = (uint8_t)(acc[1] / n);
                out[di + 2] = (uint8_t)(acc[2] / n);
            }
        }
    }

    // Fractions of the picture, so the panel can draw them at whatever size it
    // ends up giving the thumbnail.
    std::vector<float> pts;
    if (points.xy && points.count) {
        const size_t thin = points.count / kMaxOverlay + 1;
        pts.reserve(points.count / thin * 2);
        for (size_t i = 0; i < points.count; i += thin) {
            pts.push_back(points.xy[i * 2 + 0] / (float)w);
            pts.push_back(points.xy[i * 2 + 1] / (float)h);
        }
    }

    std::lock_guard<std::mutex> lk(_mu);
    Slot& s = _slots[_seq % kFilmSlots];
    s.rgb.swap(out);
    s.w = dw;
    s.h = dh;
    s.name = name;
    s.pts.swap(pts);
    s.seq = ++_seq;
}

void FilmStrip::clear() {
    std::lock_guard<std::mutex> lk(_mu);
    for (Slot& s : _slots) {
        s.rgb.clear();
        s.pts.clear();
        s.w = s.h = 0;
        s.seq = 0;
        s.uploaded = 0;
    }
    _seq = 0;
}

bool FilmStrip::has_frames() const {
    std::lock_guard<std::mutex> lk(_mu);
    return _seq != 0;
}

bool FilmStrip::draw(float height) {
    std::lock_guard<std::mutex> lk(_mu);
    if (!_seq) return false;

    // Oldest first, so the newest frame is always at the right-hand end and
    // the strip does not appear to jump when the ring wraps.
    int order[kFilmSlots];
    int n = 0;
    for (int i = 0; i < kFilmSlots; i++)
        if (_slots[i].seq) order[n++] = i;
    std::sort(order, order + n,
              [this](int a, int b) { return _slots[a].seq < _slots[b].seq; });

    // A thumbnail is worth about this much height and no more; a column tall
    // enough for several rows gets several rows rather than one enormous one.
    const float row_h = std::min(height, px(150.0f));
    ImGui::BeginChild("##film", ImVec2(0, height + px(4.0f)));
    const float avail_w = ImGui::GetContentRegionAvail().x;
    const float gap = ImGui::GetStyle().ItemSpacing.x;
    float row_x = 0.0f;
    for (int i = 0; i < n; i++) {
        Slot& s = _slots[order[i]];
        if (s.uploaded != s.seq && s.w > 0) {
            if (!s.tex) glGenTextures(1, &s.tex);
            glBindTexture(GL_TEXTURE_2D, s.tex);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
            glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, s.w, s.h, 0, GL_RGB,
                         GL_UNSIGNED_BYTE, s.rgb.data());
            s.uploaded = s.seq;
        }
        if (!s.tex) continue;
        const float w = row_h * (float)s.w / (float)std::max(1, s.h);
        if (row_x > 0.0f && row_x + w <= avail_w) {
            ImGui::SameLine();
        } else if (row_x > 0.0f) {
            row_x = 0.0f;
        }
        row_x += w + gap;
        const ImVec2 at = ImGui::GetCursorScreenPos();
        ImGui::Image((ImTextureID)(intptr_t)s.tex, ImVec2(w, row_h));
        // Feature positions, drawn over the picture rather than baked into it:
        // the same thumbnail is shown before extraction reaches it.
        if (!s.pts.empty()) {
            ImDrawList* dl = ImGui::GetWindowDrawList();
            const ImU32 col = IM_COL32(120, 255, 140, 220);
            for (size_t k = 0; k + 1 < s.pts.size(); k += 2) {
                const float x = at.x + s.pts[k] * w;
                const float y = at.y + s.pts[k + 1] * row_h;
                dl->AddRectFilled(ImVec2(x, y), ImVec2(x + 1.5f, y + 1.5f), col);
            }
        }
        if (ImGui::IsItemHovered() && !s.name.empty())
            ui::SetTooltipRaw(s.name);
    }
    ImGui::EndChild();
    return true;
}

void FilmStrip::destroy_gl() {
    std::lock_guard<std::mutex> lk(_mu);
    for (Slot& s : _slots) {
        if (s.tex) glDeleteTextures(1, &s.tex);
        s.tex = 0;
        s.uploaded = 0;
    }
}

}  // namespace gui
