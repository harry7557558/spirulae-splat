#pragma once

// A 16-bit grayscale PNG writer.
//
// stb_image_write only emits 8 bits, and the trainer reads depth as 16-bit
// (data/DataManager.cpp's decode_depth_into refuses anything else) because 8
// is visibly banded in a log-space loss. This is the smallest thing that
// closes that gap: one IDAT of unfiltered big-endian scanlines through miniz,
// which is already vendored and already compiled into the engine.

#include <cstdint>
#include <string>

namespace app {

// `data` is w*h host-order uint16 samples, row-major. False on any write or
// compression failure, with nothing left behind.
bool save_depth_png16(const std::string& path, const uint16_t* data, int width,
                      int height);

}  // namespace app
