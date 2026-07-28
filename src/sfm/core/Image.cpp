// The decode path declared in sfm/core/Image.h.
//
// stb_image is instantiated once for the whole repository, in
// src/external/stb_image_impl.cpp, which cmake/SsplatSfm.cmake adds to this
// library. Do NOT define STB_IMAGE_IMPLEMENTATION here: ssplat-gui links both
// this library and the engine, which carries that TU too.
#include "sfm/core/Image.h"

#include <algorithm>

#include "external/stb_image.h"

namespace sfm {

// Box-average an interleaved-RGB uint8 image to (dw,dh). A cheap area filter is
// enough: the color is only ever point-sampled at keypoints for the point cloud,
// not fed to SIFT. Downscaling avoids holding a full-res color buffer alongside
// the gray one, which matters for the batch decoder's memory budget.
static std::vector<uint8_t> downscaleRgb(const unsigned char* src, int w, int h, int dw, int dh) {
    std::vector<uint8_t> out((size_t)dw * dh * 3);
    for (int y = 0; y < dh; y++) {
        int sy0 = (int)((int64_t)y * h / dh), sy1 = std::max(sy0 + 1, (int)((int64_t)(y + 1) * h / dh));
        for (int x = 0; x < dw; x++) {
            int sx0 = (int)((int64_t)x * w / dw), sx1 = std::max(sx0 + 1, (int)((int64_t)(x + 1) * w / dw));
            uint32_t acc[3] = {0, 0, 0}, n = 0;
            for (int sy = sy0; sy < sy1; sy++)
                for (int sx = sx0; sx < sx1; sx++) {
                    const unsigned char* p = src + 3 * ((size_t)sy * w + sx);
                    acc[0] += p[0]; acc[1] += p[1]; acc[2] += p[2]; n++;
                }
            uint8_t* o = &out[3 * ((size_t)y * dw + x)];
            o[0] = (uint8_t)(acc[0] / n); o[1] = (uint8_t)(acc[1] / n); o[2] = (uint8_t)(acc[2] / n);
        }
    }
    return out;
}

GrayImage loadGrayImage(const std::string& path, int max_image_size, bool want_color,
                        const std::string& mask_path) {
    int w = 0, h = 0, chan = 0;
    // Force 3 channels; we do our own luma so behavior is decoder-independent.
    unsigned char* rgb = stbi_load(path.c_str(), &w, &h, &chan, 3);
    if (!rgb)
        throw std::runtime_error("cannot decode image " + path + ": " + stbi_failure_reason());

    GrayImage img;
    img.width = w;
    img.height = h;
    img.orig_width = w;
    img.orig_height = h;
    img.data.resize((size_t)w * h);
    // Rec.601 luma, matching COLMAP's FreeImage grayscale conversion.
    const float kR = 0.299f / 255.0f, kG = 0.587f / 255.0f, kB = 0.114f / 255.0f;
    for (size_t i = 0; i < img.pixels(); i++)
        img.data[i] = kR * rgb[3 * i] + kG * rgb[3 * i + 1] + kB * rgb[3 * i + 2];

    // Downscale so the long edge is at most max_image_size (COLMAP default 3200).
    int dw = w, dh = h;
    int longEdge = std::max(w, h);
    if (max_image_size > 0 && longEdge > max_image_size) {
        double scale = (double)max_image_size / longEdge;
        dw = std::max(1, (int)std::lround(w * scale));
        dh = std::max(1, (int)std::lround(h * scale));
    }
    // Keep the color companion at the *gray* (post-downscale) resolution, so a
    // keypoint's coordinates index it directly.
    if (want_color) {
        img.rgb = (dw == w && dh == h) ? std::vector<uint8_t>(rgb, rgb + (size_t)w * h * 3)
                                       : downscaleRgb(rgb, w, h, dw, dh);
    }
    stbi_image_free(rgb);
    if (dw != w || dh != h) {
        std::vector<uint8_t> keepRgb = std::move(img.rgb);
        img = resizeGray(img, dw, dh);
        img.rgb = std::move(keepRgb);
        img.orig_width = w;   // resizeGray builds a fresh image; restore the source size
        img.orig_height = h;
    }
    // Kept at the mask file's own resolution: applyMask() samples it in uv, so
    // resampling it to match `img` would only lose detail (D39).
    if (!mask_path.empty()) img.mask = loadMask(mask_path);
    img.exif = readExif(path);  // header bytes only; see sfm/core/Exif.h
    return img;
}

Mask loadMask(const std::string& path) {
    int w = 0, h = 0, chan = 0;
    // One channel: a mask is categorical, and every convention in the wild
    // (1-bit PNG, 8-bit gray, RGB white-on-black, RGBA alpha) reduces to the
    // same thing under stb's gray conversion -- except an alpha-only mask,
    // which stb would flatten to white. Masks that carry their signal in alpha
    // are handled below.
    unsigned char* px = stbi_load(path.c_str(), &w, &h, &chan, 1);
    if (!px || w <= 0 || h <= 0) {
        if (px) stbi_image_free(px);
        return Mask();
    }
    Mask m;
    m.width = w;
    m.height = h;
    m.bits.resize((size_t)w * h);
    for (size_t i = 0; i < m.bits.size(); i++) m.bits[i] = px[i] != 0 ? 1 : 0;
    stbi_image_free(px);

    // An RGBA mask that is uniformly white in RGB carries its shape in alpha
    // (what "cut out the subject" exporters produce). stb's gray conversion
    // drops alpha, so that mask decodes as all-ones; re-read the alpha channel
    // and use it instead. Only done when the gray read was fully saturated, so
    // an ordinary opaque RGBA mask is untouched.
    if (chan == 4) {
        bool all_keep = true;
        for (uint8_t b : m.bits)
            if (!b) { all_keep = false; break; }
        if (all_keep) {
            int w2 = 0, h2 = 0, c2 = 0;
            unsigned char* rgba = stbi_load(path.c_str(), &w2, &h2, &c2, 4);
            if (rgba) {
                if (w2 == w && h2 == h)
                    for (size_t i = 0; i < m.bits.size(); i++)
                        m.bits[i] = rgba[4 * i + 3] != 0 ? 1 : 0;
                stbi_image_free(rgba);
            }
        }
    }
    return m;
}

bool imageSize(const std::string& path, int& width, int& height) {
    int comp = 0;
    return stbi_info(path.c_str(), &width, &height, &comp) != 0;
}

}  // namespace sfm
