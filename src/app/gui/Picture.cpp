// Picture.cpp -- see Picture.h.

#include "app/gui/Picture.h"

#include "app/DepthColor.h"
#include "app/FrameMask.h"          // app::load_rgb, app::load_stencil
#include "external/stb_image.h"

#include <algorithm>
#include <climits>
#include <cmath>

namespace gui {

namespace {

// The masked-out region, tinted the way the mask preview tints it so that the
// two read as one answer.
void tint(uint8_t* px, int r, int g, int b) {
    px[0] = (uint8_t)(r / 3 + 150);
    px[1] = (uint8_t)(g / 3);
    px[2] = (uint8_t)(b / 3);
}

// Area average to an exact size. The row's panels arrive at three different
// resolutions and have to line up before they can be laid side by side.
void scale_to(const uint8_t* src, int sw, int sh, int dw, int dh,
              std::vector<uint8_t>& dst) {
    dst.assign((size_t)dw * dh * 3, 0);
    for (int y = 0; y < dh; y++) {
        const int sy0 = (int)((int64_t)y * sh / dh);
        const int sy1 = std::max(sy0 + 1, (int)((int64_t)(y + 1) * sh / dh));
        for (int x = 0; x < dw; x++) {
            const int sx0 = (int)((int64_t)x * sw / dw);
            const int sx1 = std::max(sx0 + 1, (int)((int64_t)(x + 1) * sw / dw));
            int acc[3] = {0, 0, 0}, n = 0;
            for (int sy = sy0; sy < sy1 && sy < sh; sy++)
                for (int sx = sx0; sx < sx1 && sx < sw; sx++) {
                    const uint8_t* p = &src[((size_t)sy * sw + sx) * 3];
                    for (int c = 0; c < 3; c++) acc[c] += p[c];
                    n++;
                }
            if (!n) continue;
            uint8_t* d = &dst[((size_t)y * dw + x) * 3];
            for (int c = 0; c < 3; c++) d[c] = (uint8_t)(acc[c] / n);
        }
    }
}

// A 16-bit depth PNG through the viewport's ramp. Read here rather than
// through app::load_rgb, which would hand back the high byte as grey.
bool load_depth_rgb(const std::string& path, int& w, int& h,
                    std::vector<uint8_t>& rgb) {
    int ch = 0;
    stbi_us* img = stbi_load_16(path.c_str(), &w, &h, &ch, 1);
    if (!img) return false;
    const size_t n = (size_t)w * h;
    std::vector<float> d(n);
    for (size_t i = 0; i < n; i++) d[i] = (float)img[i];
    stbi_image_free(img);
    rgb.resize(n * 3);
    // 0 is the trainer's "no ground truth here" and must not set the range.
    app::depth_to_rgb(d.data(), n, /*skip_zero=*/true, rgb.data());
    return true;
}

}  // namespace

void make_picture(const uint8_t* rgb, int w, int h, const uint8_t* mask,
                  int max_side, Picture& out) {
    out = Picture{};
    if (!rgb || w <= 0 || h <= 0) return;
    const int step = max_side > 0
                         ? std::max(1, (std::max(w, h) + max_side - 1) / max_side)
                         : 1;
    const int dw = std::max(1, w / step), dh = std::max(1, h / step);
    out.rgb.assign((size_t)dw * dh * 3, 0);
    out.w = dw;
    out.h = dh;
    out.src_w = w;
    out.src_h = h;
    // Full resolution answers every pane there is.
    out.made_for = max_side > 0 ? max_side : INT_MAX;

    if (step == 1) {
        // The picture itself, with the mask laid over it: no resampling at all
        // is what a pane at least as large as the file gets.
        std::copy(rgb, rgb + (size_t)w * h * 3, out.rgb.begin());
        if (mask)
            for (size_t i = 0; i < (size_t)w * h; i++)
                if (mask[i] <= 127)
                    tint(&out.rgb[i * 3], rgb[i * 3], rgb[i * 3 + 1],
                         rgb[i * 3 + 2]);
        return;
    }

    // Box average over the step x step source block: a point sample of a 4K
    // frame decimated 16x aliases into noise, which reads as a bad mask.
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
            uint8_t* px = &out.rgb[((size_t)y * dw + x) * 3];
            if (mask && keep * 2 < n) {
                tint(px, acc[0] / n, acc[1] / n, acc[2] / n);
            } else {
                px[0] = (uint8_t)(acc[0] / n);
                px[1] = (uint8_t)(acc[1] / n);
                px[2] = (uint8_t)(acc[2] / n);
            }
        }
    }
}

bool load_picture(const std::string& image_path, const std::string& mask_path,
                  int max_side, Picture& out) {
    out = Picture{};
    int w = 0, h = 0;
    std::vector<uint8_t> rgb;
    if (image_path.empty() || !app::load_rgb(image_path, w, h, rgb)) return false;

    std::vector<uint8_t> mask;
    int mw = 0, mh = 0;
    if (!mask_path.empty() && app::load_stencil(mask_path, mw, mh, mask) &&
        (mw != w || mh != h)) {
        // A mask is written at its image's size, so this is a mask that came
        // with the capture rather than one the run made. Nearest is enough for
        // a picture that is about to be averaged down anyway.
        std::vector<uint8_t> fit((size_t)w * h);
        for (int y = 0; y < h; y++) {
            const int sy = mh > 0 ? std::min(mh - 1, y * mh / h) : 0;
            for (int x = 0; x < w; x++) {
                const int sx = mw > 0 ? std::min(mw - 1, x * mw / w) : 0;
                fit[(size_t)y * w + x] = mask[(size_t)sy * mw + sx];
            }
        }
        mask.swap(fit);
    }
    make_picture(rgb.data(), w, h,
                 mask.size() == (size_t)w * h ? mask.data() : nullptr, max_side,
                 out);
    return !out.empty();
}

bool load_picture_row(const std::vector<PicturePanel>& panels, int max_side,
                      Picture& out) {
    out = Picture{};
    struct Loaded { int w = 0, h = 0; std::vector<uint8_t> rgb; };
    std::vector<Loaded> got;
    double aspect_sum = 0.0;
    int tallest = 0;
    for (const PicturePanel& p : panels) {
        Loaded l;
        const bool ok = p.depth ? load_depth_rgb(p.path, l.w, l.h, l.rgb)
                                : app::load_rgb(p.path, l.w, l.h, l.rgb);
        if (!ok || l.w <= 0 || l.h <= 0) continue;
        aspect_sum += (double)l.w / (double)l.h;
        tallest = std::max(tallest, l.h);
        got.push_back(std::move(l));
    }
    if (got.empty() || aspect_sum <= 0.0) return false;

    // The row's LONG edge is what the pane budgets, and that edge is its
    // width: at a common height h the row is h * sum(aspect) wide.
    int h0 = tallest;
    if (max_side > 0)
        h0 = std::clamp((int)std::lround(max_side / aspect_sum), 32, tallest);

    std::vector<int> widths(got.size());
    int total = 0;
    for (size_t i = 0; i < got.size(); i++) {
        widths[i] = std::max(1, (int)std::lround(
                                    (double)h0 * got[i].w / got[i].h));
        total += widths[i];
    }

    std::vector<uint8_t> row((size_t)total * h0 * 3, 0);
    std::vector<uint8_t> panel;
    int x0 = 0;
    for (size_t i = 0; i < got.size(); i++) {
        scale_to(got[i].rgb.data(), got[i].w, got[i].h, widths[i], h0, panel);
        for (int y = 0; y < h0; y++)
            std::copy(panel.begin() + (size_t)y * widths[i] * 3,
                      panel.begin() + (size_t)(y + 1) * widths[i] * 3,
                      row.begin() + ((size_t)y * total + x0) * 3);
        x0 += widths[i];
    }
    make_picture(row.data(), total, h0, nullptr, max_side, out);
    return !out.empty();
}

}  // namespace gui
