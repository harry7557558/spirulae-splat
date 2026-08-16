#pragma once
// The host-side image the inference layer takes in and hands back, plus the
// stb-backed readers and writers for it.
//
// It lives here rather than in a model's public header because everything
// above nn/ speaks it: sam/ encodes it, video/ decodes into it, and the GUI
// blits it. Pixels are 8-bit and interleaved; device tensors are nn::Tensor.

#include <cstdint>
#include <string>
#include <vector>

namespace nn {

struct Image {
    int width = 0;
    int height = 0;
    int channels = 3;              // always 3 (RGB) once loaded
    std::vector<uint8_t> data;     // row-major, `channels` interleaved
    bool empty() const { return data.empty(); }
};

// Decodes to RGB. An unreadable or unsupported file logs and returns empty().
// `gamut` / `is_linear` describe the file (core/ColorSpace.h names); pixels
// convert to sRGB, which is what every model here was trained on.
Image load_image(const std::string& path, const std::string& gamut = "",
                 bool is_linear = false);

// Writes RGB. `quality` in 0..100 selects JPEG, anything else lossless PNG --
// the same convention as reference/scripts/extract_frames.py.
bool save_image(const Image& image, const std::string& path, int quality);

// Writes a single-channel 8-bit PNG (masks, sharpness maps).
bool save_gray_png(const uint8_t* data, int width, int height,
                   const std::string& path);

}  // namespace nn
