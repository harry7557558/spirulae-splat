// MatchMatrix.cpp -- see MatchMatrix.h.

#include "app/gui/MatchMatrix.h"

#include "app/gui/Layout.h"
#include "app/gui/Ui.h"
#include "i18n/catalog/Dataset.h"

#include <cmath>
#include <vector>

namespace dmsg = spirula::i18n::msg::dataset;

namespace gui {

MatchMatrix::~MatchMatrix() = default;

void MatchMatrix::set(const PairMatrix& m) {
    _m = m;
    _dirty = true;
}

void MatchMatrix::clear() {
    _m = PairMatrix{};
    _dirty = true;
}

void MatchMatrix::draw(float size) {
    if (_m.empty()) return;
    if (_dirty) {
        // Log scale: inlier counts are heavy-tailed -- a handful of adjacent
        // pairs carry thousands and the loop closures that matter carry tens,
        // and on a linear ramp the second kind is black.
        const float top = std::log1p((float)_m.peak);
        std::vector<uint8_t> px((size_t)_m.bins * _m.bins * 3);
        for (size_t i = 0; i < (size_t)_m.bins * _m.bins; i++) {
            const float t = top > 0.0f ? std::log1p((float)_m.counts[i]) / top : 0.0f;
            // Dark blue -> cyan -> white, which keeps the low end visible
            // against the panel and still separates the top decade.
            px[i * 3 + 0] = (uint8_t)(255.0f * std::max(0.0f, t * 1.6f - 0.6f));
            px[i * 3 + 1] = (uint8_t)(255.0f * std::min(1.0f, t * 1.5f));
            px[i * 3 + 2] = (uint8_t)(255.0f * std::min(1.0f, 0.15f + t * 1.2f));
        }
        if (!_tex) glGenTextures(1, &_tex);
        glBindTexture(GL_TEXTURE_2D, _tex);
        // NEAREST: a cell is a fact, not a sample of a continuum.
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, (GLsizei)_m.bins, (GLsizei)_m.bins,
                     0, GL_RGB, GL_UNSIGNED_BYTE, px.data());
        _dirty = false;
    }
    if (!_tex) return;

    const ImVec2 at = ImGui::GetCursorScreenPos();
    ImGui::Image((ImTextureID)(intptr_t)_tex, ImVec2(size, size));
    if (!ImGui::IsItemHovered()) return;

    const ImVec2 p = ImGui::GetIO().MousePos;
    const int c = (int)((p.x - at.x) / size * _m.bins);
    const int r = (int)((p.y - at.y) / size * _m.bins);
    if (r < 0 || c < 0 || r >= (int)_m.bins || c >= (int)_m.bins) return;
    // A cell is a block of images once the capture is larger than the matrix,
    // so the readout names the range rather than pretending to name a pair.
    const uint32_t lo_r = (uint32_t)((uint64_t)r * _m.n_images / _m.bins);
    const uint32_t hi_r = (uint32_t)(((uint64_t)r + 1) * _m.n_images / _m.bins);
    const uint32_t lo_c = (uint32_t)((uint64_t)c * _m.n_images / _m.bins);
    const uint32_t hi_c = (uint32_t)(((uint64_t)c + 1) * _m.n_images / _m.bins);
    if (hi_r - lo_r <= 1 && hi_c - lo_c <= 1)
        ui::SetTooltip(dmsg::matrix_cell_pair,
                       {(long long)lo_r, (long long)lo_c,
                        (long long)_m.at((uint32_t)r, (uint32_t)c)});
    else
        ui::SetTooltip(dmsg::matrix_cell_range,
                       {(long long)lo_r, (long long)(hi_r - 1), (long long)lo_c,
                        (long long)(hi_c - 1),
                        (long long)_m.at((uint32_t)r, (uint32_t)c)});
}

void MatchMatrix::destroy_gl() {
    if (_tex) glDeleteTextures(1, &_tex);
    _tex = 0;
    _dirty = true;
}

}  // namespace gui
