// Picture.cpp -- see Picture.h.

#include "app/gui/Picture.h"

#include "app/FrameMask.h"          // app::load_rgb, app::load_stencil

#include <algorithm>
#include <climits>

namespace gui {

namespace {

// The masked-out region, tinted the way the mask preview tints it so that the
// two read as one answer.
void tint(uint8_t* px, int r, int g, int b) {
    px[0] = (uint8_t)(r / 3 + 150);
    px[1] = (uint8_t)(g / 3);
    px[2] = (uint8_t)(b / 3);
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

}  // namespace gui
